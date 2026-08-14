# oxo-flow-enrichment

[![CI](https://github.com/oxo-flow-community/oxo-flow-enrichment/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-enrichment/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

Region set and gene set enrichment analysis: LOLA (region overlap),
rGREAT (genomic region enrichment of annotated terms), pycisTarget (region
TFBS motif enrichment), GSEApy over-representation analysis (ORA) and
preranked GSEA — with multiple-test correction inside each tool, per-set
enrichment plots, per-group summary plots, and reproducibility exports.

## Source

Ported from **[epigen/enrichment_analysis](https://github.com/epigen/enrichment_analysis)**,
version `v3.0.1` (MIT). This port is maintained independently and **may lag
the upstream** — check the tag above and the fidelity table below for the
exact ported state.

## Fidelity

| Upstream process/rule | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| prepare_databases | `prepare_databases_Azimuth_2023`, `prepare_databases_Reactome` | gseapy 1.1.3 | identical command; database fan-out baked as static blocks (2 default-path databases) |
| region_enrichment_analysis_LOLA | `region_enrichment_analysis_LOLA` | bioconductor-lola 1.32.0 | identical command; database fan-out baked as static block (1 default-path database) |
| region_enrichment_analysis_GREAT | `region_enrichment_analysis_GREAT_Azimuth_2023`, `region_enrichment_analysis_GREAT_Reactome` | bioconductor-rgreat 2.4.0 | identical command; upstream `great_parameters` nested dict flattened into `great_*` config keys |
| region_gene_association_GREAT | `region_gene_association_GREAT` | bioconductor-rgreat 2.4.0 | identical command; uses the first database (Azimuth_2023) as upstream |
| region_motif_enrichment_analysis_pycisTarget | `region_motif_enrichment_analysis_pycisTarget` | pycistarget 1.1 | command text verbatim (incl. upstream error-tolerance wrapper); threads=10 as upstream |
| process_results_pycisTarget | `process_results_pycisTarget` | pycistarget 1.1 | identical command |
| gene_motif_enrichment_analysis_RcisTarget | not ported | RcisTarget | zero instances on the default path: needs `.txt` gene sets in the annotation, the default annotation has none (region sets + ranked sets only) |
| gene_ORA_GSEApy | `gene_ORA_GSEApy_Azimuth_2023`, `gene_ORA_GSEApy_Reactome` | gseapy 1.1.3 | identical command; upstream genes_dict fan-out has zero default-path members, region-set fan-out kept |
| gene_preranked_GSEApy | `gene_preranked_GSEApy_Azimuth_2023`, `gene_preranked_GSEApy_Reactome` | gseapy 1.1.3 | identical command |
| plot_enrichment_result | `plot_enrichment_result_*` (8 blocks) | r-ggplot2 3.5.0, r-svglite 2.1.0 | identical command; upstream wildcard fan-out (tool × db × feature_set) baked as per-(tool,db) scatter blocks |
| aggregate | `aggregate_*` (8 blocks) | pandas 1.1.4 / 1.5.3 | identical logic; upstream wildcards group/tool/db passed as CLI args |
| visualize | `visualize_*` (8 blocks) | r-ggplot2 3.5.0, r-pheatmap 1.0.12 | identical command/logic; `cluster_summary` config key kept as upstream numeric flag |
| config_export | `config_export` | — | upstream dumps the in-memory config dict; the port copies `config/config.yaml` (effective-config mirror) |
| annot_export | `annot_export` | — | identical command |
| env_export | not ported | — | `conda env export` needs the conda CLI inside the runtime env; exact pins are already declared in `envs/*.yaml` |
| report rendering | not ported | — | oxo-flow has no report module; reproducibility exports (configs/, envs) are kept as upstream rules |

Script ports: upstream scripts run inside snakemake's `snakemake@input/...`
namespace; the port passes the same values as positional CLI arguments
(`scripts/*`), keeping every analysis step and output byte-identical.
`utils.R` is copied verbatim. Fidelity conventions: `{config.a.b}` nested
access does not exist in oxo-flow — all upstream nested config dicts
(`great_parameters`, `pycistarget_parameters`, `column_names`, `adjp_th`,
caps) are flattened into prefixed top-level keys; the pycisTarget
`annotations_to_use` list is carried as a python-list literal string so the
rendered command is byte-identical to upstream.

## Quickstart

```bash
# 1. install oxo-flow (see Requirements)
# 2. prepare data (see test/fixtures/; region beds, ranked gene lists,
#    databases, LOLA regionDB, pycisTarget context db + motif annotations)
# 3. preview the plan
oxo-flow dry-run main.oxoflow
# 4. run
oxo-flow run main.oxoflow -j 4
# 5. run a subset
oxo-flow run main.oxoflow -t aggregate_LOLA_LOLACore_ATAC
```

## Requirements

- **oxo-flow ≥ 0.11.0** — install the prebuilt binary:

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/download/v0.11.0/oxo-flow-v0.11.0-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz
sudo mv oxo-flow /usr/local/bin/
```

- Conda users may alternatively `conda install -c bioconda oxo-flow-cli`
  (note: the bioconda package currently lags the release binary at 0.10.2 —
  some 0.11.0 format features may not validate).
- Conda at runtime, per the environments declared in `main.oxoflow`
  (`envs/*.yaml`, exact pins from upstream).

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community. Upstream attribution in
[NOTICE.md](NOTICE.md).

## Community

https://oxo-flow-community.github.io/
