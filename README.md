# DBiT-spatial-Brain

Spatial brain-domain analysis for DBiT data. The workflow runs RCTD cell-type
deconvolution, plots the RCTD weights, and performs BANKSY spatial clustering.

## Setup

Install [Pixi](https://pixi.prefix.dev/) and create the project environments:

```bash
pixi install -a
```

Copy and edit the DBiT-spatial-Brain configuration:

```bash
cp config/dbit.brain.config.example.sh config/dbit.brain.config.sh
```

Pass the sample `dbit.config.sh` from DBiT-spatial-DARLIN to `run.sh`. Runtime,
spatial-layout, and shared Slurm settings are read from that sample config;
RCTD, BANKSY, and domain-specific Slurm settings are read from
`config/dbit.brain.config.sh`.

`mrna_dir` should contain `raw/clustered.tissuefiltered.h5ad` from the mRNA
workflow. Results are written to `<mrna_output_path>/deconv`.

## Run

Set `execution_mode=local` or `execution_mode=hpc` in the sample
`dbit.config.sh`. Shared Slurm settings are also read from that file, while
domain resources are configured through the `sbatch_domain_*` fields in
`config/dbit.brain.config.sh`.

```bash
bash script/run.sh domain --config /path/to/sample/dbit.config.sh
```

In local mode, `run.sh` executes `domain.sh` directly. In HPC mode, it submits
`domain.sh` through `sbatch --wrap` using its real absolute path, together with
the configured CPU, partition, memory, time, log paths, and requeue setting.

When `scratch` is configured, temporary files are placed on that disk. A failed
run copies available RCTD outputs back before removing its temporary directory.
