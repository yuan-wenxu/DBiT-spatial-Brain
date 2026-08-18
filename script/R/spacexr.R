#!/usr/bin/env Rscript

# Run spacexr/RCTD cell-type deconvolution with a 10x reference atlas.
# The spatial H5AD must store raw integer counts in CSR-encoded X, barcodes in
# the obs index, gene names in var/gene_name, and coordinates in obsm/spatial.

usage <- function(status = 0L) {
  cat(paste0(
    "usage: 01.spacexr.R [-h] --reference-dir REFERENCE_DIR\n",
    "                     --spatial-h5ad SPATIAL_H5AD --output-dir OUTPUT_DIR\n",
    "                     [--reference-barcode-column NAME]\n",
    "                     [--reference-numi-column NAME]\n",
    "                     [--cell-type-column NAME]\n",
    "                     [--dataset-column NAME]\n",
    "                     [--reference-gene-column INT]\n",
    "                     [--max-cells-per-type INT]\n",
    "                     [--reference-cache FILE] [--cores INT]\n",
    "                     [--doublet-mode MODE] [--reference-min-umi INT]\n",
    "                     [--spatial-min-umi INT] [--seed INT] [--test-mode]\n\n",
    "Run spacexr/RCTD using a 10x reference atlas and a spatial AnnData file.\n",
    "The AnnData X matrix must contain raw integer counts.\n\n",
    "options:\n",
    "  -h, --help                    show this help message and exit\n",
    "  --reference-dir REFERENCE_DIR\n",
    "                                10x reference directory containing matrix.mtx.gz,\n",
    "                                features.tsv.gz, barcodes.tsv.gz, and metadata.tsv.gz\n",
    "  --spatial-h5ad SPATIAL_H5AD   spatial AnnData file\n",
    "  --output-dir OUTPUT_DIR       new or empty output directory\n",
    "  --reference-barcode-column NAME\n",
    "                                reference metadata barcode column (default: V1)\n",
    "  --reference-numi-column NAME  reference metadata UMI column (default: nCount_RNA)\n",
    "  --cell-type-column NAME       reference metadata annotation column (default: C66_named)\n",
    "  --dataset-column NAME         reference dataset column used for balanced sampling\n",
    "                                within each cell type (default: Dataset)\n",
    "  --reference-gene-column INT   gene-name column in features.tsv (default: 2)\n",
    "  --max-cells-per-type INT      reference sampling cap per cell type (default: 200)\n",
    "  --reference-cache FILE        reuse or save a sampled spacexr Reference RDS\n",
    "  --cores INT                   number of RCTD workers (default: 8)\n",
    "  --doublet-mode MODE           RCTD mode: full, doublet, or multi (default: full)\n",
    "  --reference-min-umi INT       minimum UMI count per reference cell (default: 100)\n",
    "  --spatial-min-umi INT         minimum UMI count per spatial spot (default: 100)\n",
    "  --seed INT                    reference sampling seed (default: 1)\n",
    "  --test-mode                   run a fast RCTD smoke test\n"
  ))
  quit(save = "no", status = status)
}

parse_args <- function(args) {
  out <- list(
    reference_barcode_column = "V1",
    reference_numi_column = "nCount_RNA",
    cell_type_column = "C66_named",
    dataset_column = "Dataset",
    reference_gene_column = "2",
    max_cells_per_type = "200",
    cores = "8",
    doublet_mode = "full",
    reference_min_umi = "100",
    spatial_min_umi = "100",
    seed = "1",
    test_mode = FALSE
  )
  if (!length(args) || "--help" %in% args) usage(0L)
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) stop("Unexpected positional argument: ", token)
    key <- gsub("-", "_", substring(token, 3L), fixed = TRUE)
    if (key == "test_mode") {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
        stop("Missing value for ", token)
      }
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

find_one <- function(directory, stem) {
  candidates <- file.path(directory, c(stem, paste0(stem, ".gz")))
  found <- candidates[file.exists(candidates)]
  if (length(found) != 1L) {
    stop("Expected exactly one of: ", paste(candidates, collapse = ", "))
  }
  found[[1L]]
}

read_tsv <- function(path, header = FALSE, ...) {
  utils::read.delim(
    path,
    header = header,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    ...
  )
}

assert_integer_counts <- function(x, label) {
  if (!inherits(x, "Matrix")) stop(label, " must be a Matrix sparse matrix")
  if (anyNA(x@x) || any(x@x < 0)) stop(label, " contains NA or negative values")
  if (any(abs(x@x - round(x@x)) > .Machine$double.eps^0.5)) {
    stop(label, " is not raw integer counts")
  }
}

collapse_duplicate_genes <- function(counts) {
  genes <- rownames(counts)
  unique_genes <- unique(genes)
  if (length(unique_genes) == length(genes)) return(counts)
  message("Collapsing ", length(genes) - length(unique_genes), " duplicate gene rows")
  mapping <- Matrix::sparseMatrix(
    i = match(genes, unique_genes),
    j = seq_along(genes),
    x = 1,
    dims = c(length(unique_genes), length(genes))
  )
  out <- mapping %*% counts
  rownames(out) <- unique_genes
  colnames(out) <- colnames(counts)
  out
}

read_h5ad_string_vector <- function(file, path) {
  object <- file[[path]]
  encoding <- hdf5r::h5attr(object, "encoding-type")

  if (inherits(object, "H5D")) {
    return(as.character(object[]))
  }

  if (!inherits(object, "H5Group")) {
    stop("Unsupported H5AD string object at ", path)
  }

  if (identical(encoding, "categorical")) {
    codes <- as.integer(object[["codes"]][])
    categories <- as.character(object[["categories"]][])
    values <- rep(NA_character_, length(codes))
    present <- codes >= 0L
    values[present] <- categories[codes[present] + 1L]
    return(values)
  }

  if (identical(encoding, "nullable-string-array")) {
    values <- as.character(object[["values"]][])
    missing <- as.logical(object[["mask"]][])
    if (length(values) != length(missing)) {
      stop("Invalid nullable string encoding at ", path)
    }
    values[missing] <- NA_character_
    return(values)
  }

  stop("Unsupported H5AD string encoding at ", path, ": ", encoding)
}

read_spatial_h5ad <- function(path) {
  message("Reading spatial H5AD: ", path)
  file <- hdf5r::H5File$new(path, mode = "r")
  on.exit(file$close_all(), add = TRUE)

  required <- c("X/data", "X/indices", "X/indptr", "obs/barcode", "var/gene_name", "obsm/spatial")
  missing <- required[!vapply(required, function(key) file$exists(key), logical(1))]
  if (length(missing)) stop("H5AD is missing: ", paste(missing, collapse = ", "))
  if (!identical(hdf5r::h5attr(file[["X"]], "encoding-type"), "csr_matrix")) {
    stop("H5AD X must use CSR sparse encoding")
  }

  shape <- as.integer(hdf5r::h5attr(file[["X"]], "shape"))
  values <- as.numeric(file[["X/data"]][])
  indices <- as.integer(file[["X/indices"]][])
  indptr <- as.integer(file[["X/indptr"]][])
  if (length(shape) != 2L || length(indptr) != shape[[1L]] + 1L) {
    stop("Invalid H5AD CSR dimensions")
  }

  # A spots x genes CSR matrix has the same slots as its genes x spots CSC transpose.
  counts <- methods::new(
    "dgCMatrix",
    i = indices,
    p = indptr,
    x = values,
    Dim = as.integer(c(shape[[2L]], shape[[1L]]))
  )
  barcodes <- read_h5ad_string_vector(file, "obs/barcode")
  gene_names <- read_h5ad_string_vector(file, "var/gene_name")
  spatial <- t(file[["obsm/spatial"]][, ])
  if (length(barcodes) != ncol(counts) || length(gene_names) != nrow(counts)) {
    stop("H5AD X dimensions disagree with obs/var indices")
  }
  if (!identical(dim(spatial), c(ncol(counts), 2L))) {
    stop("obsm/spatial must have shape n_spots x 2")
  }
  if (anyDuplicated(barcodes)) stop("H5AD contains duplicate spot barcodes")
  if (anyNA(spatial) || any(!is.finite(spatial))) stop("H5AD spatial coordinates are invalid")

  rownames(counts) <- gene_names
  colnames(counts) <- barcodes
  counts <- collapse_duplicate_genes(counts)
  coords <- data.frame(x = spatial[, 1L], y = spatial[, 2L], row.names = barcodes)
  assert_integer_counts(counts, "Spatial X")
  list(counts = counts, coords = coords)
}

read_hypomap_metadata <- function(path, barcode_column, numi_column,
                                  cell_type_column, dataset_column) {
  command <- paste("gzip -dc", shQuote(path))
  suppressWarnings(data.table::fread(
    cmd = command,
    select = c(barcode_column, numi_column, cell_type_column, dataset_column),
    na.strings = c("NA", ""),
    showProgress = FALSE
  ))
}

sample_reference_cells <- function(metadata, barcode_column, numi_column,
                                   cell_type_column, dataset_column,
                                   max_per_type, min_umi, seed) {
  names(metadata)[names(metadata) == barcode_column] <- "barcode"
  names(metadata)[names(metadata) == numi_column] <- "nUMI"
  keep <- !is.na(metadata[[cell_type_column]]) &
    nzchar(metadata[[cell_type_column]]) &
    !is.na(metadata$nUMI) &
    metadata$nUMI >= min_umi
  metadata <- metadata[keep, , drop = FALSE]
  missing_dataset <- is.na(metadata[[dataset_column]]) |
    !nzchar(as.character(metadata[[dataset_column]]))
  if (any(missing_dataset)) {
    stop(
      "Reference cells passing annotation and UMI filters contain missing ",
      dataset_column, " values"
    )
  }

  set.seed(seed)
  positions_by_type <- split(
    seq_len(nrow(metadata)),
    metadata[[cell_type_column]],
    drop = TRUE
  )
  selected <- lapply(positions_by_type, function(type_positions) {
    positions_by_dataset <- split(
      type_positions,
      metadata[[dataset_column]][type_positions],
      drop = TRUE
    )
    positions_by_dataset <- lapply(
      positions_by_dataset,
      function(positions) {
        if (length(positions) <= 1L) positions else sample(positions, length(positions))
      }
    )

    selected_for_type <- integer()
    depth <- 1L
    while (length(selected_for_type) < max_per_type) {
      current <- vapply(
        positions_by_dataset,
        function(positions) {
          if (depth <= length(positions)) positions[[depth]] else NA_integer_
        },
        integer(1)
      )
      current <- unname(current[!is.na(current)])
      if (!length(current)) break
      remaining <- max_per_type - length(selected_for_type)
      selected_for_type <- c(selected_for_type, head(current, remaining))
      depth <- depth + 1L
    }
    selected_for_type
  })
  selected <- sort(unlist(selected, use.names = FALSE))
  sampled <- metadata[selected, , drop = FALSE]
  expected_counts <- pmin(
    table(metadata[[cell_type_column]]),
    max_per_type
  )
  observed_counts <- table(factor(
    sampled[[cell_type_column]],
    levels = names(expected_counts)
  ))
  if (!identical(as.integer(observed_counts), as.integer(expected_counts))) {
    stop("Dataset-balanced reference sampling failed cell-type count validation")
  }
  sampled
}

build_reference <- function(directory, barcode_column, numi_column, cell_type_column,
                            dataset_column, gene_column, max_per_type, min_umi, seed) {
  matrix_path <- find_one(directory, "matrix.mtx")
  feature_path <- find_one(directory, "features.tsv")
  barcode_path <- find_one(directory, "barcodes.tsv")
  metadata_path <- find_one(directory, "metadata.tsv")

  message("Reading HypoMap metadata and sampling reference cells across datasets")
  metadata <- read_hypomap_metadata(
    metadata_path, barcode_column, numi_column, cell_type_column, dataset_column
  )
  barcodes <- read_tsv(barcode_path, header = FALSE)[[1L]]
  if (nrow(metadata) != length(barcodes) ||
      !identical(metadata[[barcode_column]], barcodes)) {
    stop("HypoMap metadata row IDs do not match barcodes.tsv order")
  }
  sampled <- sample_reference_cells(
    metadata, barcode_column, numi_column, cell_type_column, dataset_column,
    max_per_type, min_umi, seed
  )
  selected_columns <- match(sampled$barcode, barcodes)
  message(
    "Selected ", nrow(sampled), " reference cells across ",
    data.table::uniqueN(sampled[[cell_type_column]]), " cell types and ",
    data.table::uniqueN(sampled[[dataset_column]]), " datasets"
  )
  rm(metadata, barcodes)
  gc(verbose = FALSE)

  message("Reading large HypoMap Matrix Market file; this is memory intensive on first run")
  connection <- if (grepl("\\.gz$", matrix_path)) gzfile(matrix_path, "rt") else matrix_path
  full_counts <- Matrix::readMM(connection)
  if (inherits(connection, "connection")) close(connection)
  features <- read_tsv(feature_path, header = FALSE)
  if (nrow(full_counts) != nrow(features)) stop("HypoMap matrix/features dimensions disagree")
  if (ncol(full_counts) < max(selected_columns)) stop("HypoMap selected column exceeds matrix dimensions")

  counts <- methods::as(full_counts[, selected_columns, drop = FALSE], "CsparseMatrix")
  rm(full_counts)
  gc(verbose = FALSE)
  if (gene_column < 1L || gene_column > ncol(features)) {
    stop("--reference-gene-column exceeds the number of features.tsv columns")
  }
  rownames(counts) <- as.character(features[[gene_column]])
  colnames(counts) <- sampled$barcode
  counts <- collapse_duplicate_genes(counts)
  assert_integer_counts(counts, "Reference counts")

  cell_types <- factor(sampled[[cell_type_column]])
  names(cell_types) <- sampled$barcode
  numi <- Matrix::colSums(counts)
  reference <- spacexr::Reference(
    counts,
    cell_types,
    nUMI = numi,
    min_UMI = min_umi,
    n_max_cells = max_per_type
  )
  list(
    reference = reference,
    barcode_column = barcode_column,
    numi_column = numi_column,
    cell_type_column = cell_type_column,
    dataset_column = dataset_column,
    sampling = "round_robin_within_cell_type_across_dataset",
    gene_column = gene_column,
    max_cells_per_type = max_per_type,
    seed = seed
  )
}

options(stringsAsFactors = FALSE)
args <- parse_args(commandArgs(trailingOnly = TRUE))
required_args <- c("reference_dir", "spatial_h5ad", "output_dir")
missing_args <- required_args[!vapply(required_args, function(x) !is.null(args[[x]]), logical(1))]
if (length(missing_args)) {
  stop("Missing required option(s): --", paste(gsub("_", "-", missing_args), collapse = ", --"))
}

for (package in c("Matrix", "data.table", "hdf5r", "spacexr")) {
  if (!requireNamespace(package, quietly = TRUE)) stop("Missing pixi R dependency: ", package)
}
if (!dir.exists(args$reference_dir)) stop("Reference directory not found: ", args$reference_dir)
if (!file.exists(args$spatial_h5ad)) stop("Spatial H5AD not found: ", args$spatial_h5ad)
if (dir.exists(args$output_dir) && length(list.files(args$output_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Output directory is not empty: ", args$output_dir)
}
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)

cores <- as.integer(args$cores)
max_per_type <- as.integer(args$max_cells_per_type)
reference_gene_column <- as.integer(args$reference_gene_column)
reference_min_umi <- as.integer(args$reference_min_umi)
spatial_min_umi <- as.integer(args$spatial_min_umi)
seed <- as.integer(args$seed)
if (anyNA(c(cores, max_per_type, reference_gene_column,
            reference_min_umi, spatial_min_umi, seed)) ||
    cores < 1L || reference_gene_column < 1L || max_per_type < 25L) {
  stop("Invalid numeric option; --max-cells-per-type must be at least 25")
}
if (!args$doublet_mode %in% c("full", "doublet", "multi")) {
  stop("--doublet-mode must be full, doublet, or multi")
}

spatial <- read_spatial_h5ad(args$spatial_h5ad)

reference_bundle <- NULL
if (!is.null(args$reference_cache) && file.exists(args$reference_cache)) {
  message("Loading cached spacexr reference: ", args$reference_cache)
  cached_reference_bundle <- readRDS(args$reference_cache)
  cache_settings <- c(
    barcode_column = args$reference_barcode_column,
    numi_column = args$reference_numi_column,
    cell_type_column = args$cell_type_column,
    dataset_column = args$dataset_column,
    sampling = "round_robin_within_cell_type_across_dataset",
    gene_column = reference_gene_column
  )
  cached_settings <- unlist(
    cached_reference_bundle[names(cache_settings)],
    use.names = TRUE
  )
  if (identical(cached_settings, cache_settings)) {
    reference_bundle <- cached_reference_bundle
  } else {
    message(
      "Reference cache settings differ from the requested Dataset-balanced ",
      "sampling; rebuilding and overwriting: ", args$reference_cache
    )
  }
}
if (is.null(reference_bundle)) {
  reference_bundle <- build_reference(
    args$reference_dir,
    args$reference_barcode_column,
    args$reference_numi_column,
    args$cell_type_column,
    args$dataset_column,
    reference_gene_column,
    max_per_type,
    reference_min_umi,
    seed
  )
  if (!is.null(args$reference_cache)) {
    dir.create(dirname(args$reference_cache), recursive = TRUE, showWarnings = FALSE)
    saveRDS(reference_bundle, args$reference_cache, compress = FALSE)
    message("Saved sampled reference cache: ", args$reference_cache)
  }
}
reference <- reference_bundle$reference

common_genes <- intersect(rownames(spatial$counts), rownames(reference@counts))
if (length(common_genes) < 100L) {
  stop("Only ", length(common_genes), " shared genes after Ensembl-to-symbol mapping")
}
message(
  "Input summary: ", ncol(reference@counts), " reference cells, ",
  nlevels(reference@cell_types), " cell types, ", ncol(spatial$counts),
  " spatial spots, ", length(common_genes), " shared genes"
)

puck <- spacexr::SpatialRNA(
  spatial$coords,
  spatial$counts,
  nUMI = Matrix::colSums(spatial$counts)
)
set.seed(seed)
rctd <- spacexr::create.RCTD(
  puck,
  reference,
  max_cores = cores,
  test_mode = isTRUE(args$test_mode),
  UMI_min = spatial_min_umi
)
rctd <- spacexr::run.RCTD(rctd, doublet_mode = args$doublet_mode)

saveRDS(rctd, file.path(args$output_dir, "rctd_result.rds"), compress = FALSE)
if (!is.null(rctd@results$weights)) {
  weights <- as.data.frame(
    as.matrix(rctd@results$weights),
    check.names = FALSE
  )
  weights <- cbind(
    barcode = rownames(weights),
    x = spatial$coords[rownames(weights), "x"],
    y = spatial$coords[rownames(weights), "y"],
    weights
  )
  utils::write.csv(weights, file.path(args$output_dir, "cell_type_weights.csv"), row.names = FALSE)
}

run_info <- data.frame(
  key = c(
    "reference_dir", "spatial_h5ad",
    "reference_barcode_column", "reference_numi_column", "cell_type_column",
    "dataset_column", "reference_sampling",
    "reference_gene_column",
    "reference_cells", "cell_types", "spatial_spots", "shared_genes",
    "doublet_mode", "cores", "test_mode", "seed", "spacexr_version", "R_version"
  ),
  value = c(
    normalizePath(args$reference_dir), normalizePath(args$spatial_h5ad),
    args$reference_barcode_column,
    args$reference_numi_column, args$cell_type_column, args$dataset_column,
    "round_robin_within_cell_type_across_dataset", reference_gene_column,
    ncol(reference@counts), nlevels(reference@cell_types), ncol(spatial$counts),
    length(common_genes), args$doublet_mode, cores, isTRUE(args$test_mode), seed,
    as.character(utils::packageVersion("spacexr")), R.version.string
  )
)
utils::write.csv(run_info, file.path(args$output_dir, "run_info.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(args$output_dir, "session_info.txt"))
message("Completed. Results written to: ", args$output_dir)
