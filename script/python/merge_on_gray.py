"""
Transform explicitly named frame images and composite them onto gray.png.

Usage:
    python merge_on_gray.py --gray <gray.png> --frame <frame.png> [--frame <frame.png> ...]
        [--orientation <mode>] [--swap_xy]

For each frame passed with --frame, the script will:
  1. Reuse gray.downsampled.png, or create it with width and height divided by 10
     while keeping each dimension at least 1500 pixels.
  2. Transform the image according to --orientation, then --swap_xy when requested.
  3. Resize the downsampled gray image to match the transformed image.
  4. Composite the RGBA overlay onto the gray background.
  5. Save as merged_<original_filename> in the same directory.
"""
import argparse
import os
from pathlib import Path

from PIL import Image, ImageOps
Image.MAX_IMAGE_PIXELS = None

ORIENTATION_CHOICES = ("normal", "horizontal", "vertical", "rotate")
DOWNSAMPLE_FACTOR = 10
MIN_DOWNSAMPLED_DIMENSION = 1500


def get_downsampled_gray(gray_path: str) -> Path:
    """Return the cached gray image, creating a 10x-per-axis version if needed."""
    source_path = Path(gray_path)
    downsampled_path = source_path.with_name("gray.downsampled.png")
    if downsampled_path.is_file():
        return downsampled_path

    with Image.open(source_path) as gray:
        downsampled_size = (
            max(MIN_DOWNSAMPLED_DIMENSION, gray.width // DOWNSAMPLE_FACTOR),
            max(MIN_DOWNSAMPLED_DIMENSION, gray.height // DOWNSAMPLE_FACTOR),
        )
        gray.resize(downsampled_size, resample=Image.Resampling.LANCZOS).save(
            downsampled_path
        )
    print(
        f"Downsampled gray image: {source_path} -> {downsampled_path} "
        f"({downsampled_size[0]}x{downsampled_size[1]})"
    )
    return downsampled_path

def set_opacity(img: Image.Image, opacity: float) -> Image.Image:
    """Return a copy of img with its alpha channel scaled by opacity (0-1)."""
    r, g, b, a = img.split()
    a = a.point(lambda x: round(x * opacity))
    return Image.merge("RGBA", (r, g, b, a))


def normalize_orientation(value: str) -> str:
    orientation = value.lower()
    if orientation not in ORIENTATION_CHOICES:
        raise argparse.ArgumentTypeError(
            f"orientation must be one of {', '.join(ORIENTATION_CHOICES)}, got: {value}"
        )
    return orientation


def transform_frame(frame: Image.Image, orientation: str, swap_xy: bool) -> Image.Image:
    if orientation == "horizontal":
        frame = ImageOps.mirror(frame)
    elif orientation == "vertical":
        frame = ImageOps.flip(frame)
    elif orientation == "rotate":
        frame = frame.transpose(Image.Transpose.ROTATE_180)
    if swap_xy:
        frame = frame.transpose(Image.Transpose.TRANSPOSE)
    return frame


def merge_on_gray(frame_path: str, gray_path: str, orientation: str, swap_xy: bool) -> None:
    frame = Image.open(frame_path).convert("RGBA")
    frame_transformed = transform_frame(frame, orientation, swap_xy)
    if "umap" in frame_path:
        frame_transformed = set_opacity(frame_transformed, 0.7)

    gray = Image.open(gray_path).convert("RGBA")
    gray_resized = set_opacity(gray.resize(frame_transformed.size, resample=Image.LANCZOS), 0.7)

    # Black opaque background ensures the final image is fully opaque
    black_bg = Image.new("RGBA", frame_transformed.size, (0, 0, 0, 255))
    result = Image.alpha_composite(black_bg, gray_resized)
    result = Image.alpha_composite(result, frame_transformed)

    out_dir = os.path.dirname(frame_path)
    out_name = "merged_" + os.path.basename(frame_path)
    result.save(os.path.join(out_dir, out_name))


def main() -> None:
    parser = argparse.ArgumentParser(description="Flip and merge frame plots onto gray background.")
    parser.add_argument("--gray", required=True, help="Path to gray.png background image")
    parser.add_argument(
        "--frame",
        action="append",
        required=True,
        help="Exact path to an overlay image; may be specified more than once",
    )
    parser.add_argument("--orientation", type=normalize_orientation, default="normal", help="Transform overlay before merging: normal, horizontal, vertical, or rotate")
    parser.add_argument("--swap_xy", action="store_true", help="Swap x and y axes after applying orientation before merging")
    args = parser.parse_args()

    if not os.path.isfile(args.gray):
        raise FileNotFoundError(f"gray image not found: {args.gray}")

    downsampled_gray = get_downsampled_gray(args.gray)

    for frame_path in args.frame:
        if not os.path.isfile(frame_path):
            print(f"Skipping missing frame: {frame_path}")
            continue
        merge_on_gray(
            frame_path, str(downsampled_gray), args.orientation, args.swap_xy
        )
        print(f"Merged: {frame_path} -> merged_{os.path.basename(frame_path)}")


if __name__ == "__main__":
    main()
