from __future__ import annotations

import re
from pathlib import Path


BASE = Path(r"C:\Users\zqr20\Documents\tjmeta\BIG_revision")
RMD = BASE / "permonva_revise.Rmd"
PUBLISH = BASE / "publish_Script"
TODAY = "2026-05-29"

START = "# >>> RMD_EXTRACT_ANNOTATION_START"
END = "# <<< RMD_EXTRACT_ANNOTATION_END"


TARGETS = [
    {
        "name": "permonva_revise_01_metadata_rebuild__no_figure.R",
        "category": "01_metadata_qc_baseline_Metadata_QC_baseline",
        "chunks": ["setup", "loading pkgs", "reorganize dataset with covariates in metadata"],
        "figures": "No direct submitted figure",
        "purpose": "Rebuilds the NM revision phyloseq object with updated clinical, batch, QC, sampling-time, and trajectory metadata.",
    },
    {
        "name": "permonva_revise_02_batch_effect_qc__FigureS10_FigureS11_FigureS12.R",
        "category": "01_metadata_qc_baseline_Metadata_QC_baseline",
        "chunks": ["setup", "loading pkgs", "batch effect plot"],
        "figures": "FigureS10; FigureS11; FigureS12",
        "purpose": "Plots PCoA/QC batch-effect diagnostics and QC metric distributions used for supplementary QC figures.",
    },
    {
        "name": "permonva_revise_03_generate_batch_permanova_scripts__no_figure.R",
        "category": "02_diversity_permanova_Diversity_PERMANOVA",
        "chunks": ["setup", "loading pkgs", "write batch permonva scripts"],
        "figures": "No direct submitted figure",
        "purpose": "Generates visit- and slice-specific PERMANOVA batch scripts for technical and biological covariate adjustment.",
    },
    {
        "name": "permonva_revise_04_offspring_permanova_bubbles__Figure2D_Figure2D_adjustbiological.R",
        "category": "02_diversity_permanova_Diversity_PERMANOVA",
        "chunks": ["setup", "loading pkgs", "offspring bubble"],
        "figures": "Figure2D; Figure2D_adjustbiological",
        "purpose": "Aggregates offspring PERMANOVA results and draws offspring covariate bubble plots with technical and biological adjustment variants.",
    },
    {
        "name": "permonva_revise_05_mother_permanova_bubbles__Figure2C_Figure2C_technical_and_biological.R",
        "category": "02_diversity_permanova_Diversity_PERMANOVA",
        "chunks": ["setup", "loading pkgs", "mother bubble"],
        "figures": "Figure2C; Figure2C_technical_and_biological",
        "purpose": "Aggregates maternal PERMANOVA results and draws maternal covariate bubble plots with technical/biological adjustment variants.",
    },
    {
        "name": "permonva_revise_06_maternal_exposure_permanova_summary__no_figure.R",
        "category": "02_diversity_permanova_Diversity_PERMANOVA",
        "chunks": ["setup", "loading pkgs", "other exposure MG"],
        "figures": "No direct submitted figure",
        "purpose": "Collects selected maternal exposure PERMANOVA results across MG and BM visits for exploratory summaries.",
    },
    {
        "name": "permonva_revise_07_maaslin2_all_de_table__no_figure.R",
        "category": "03_differential_abundance_Differential_abundance_trajectory",
        "chunks": ["setup", "loading pkgs", "organize maaslin2 table be one table"],
        "figures": "No direct submitted figure",
        "purpose": "Combines MaAsLin2 significant-results files into the all differential-abundance table.",
    },
    {
        "name": "permonva_revise_08_external_bpstrain_cazy_summary__no_figure.R",
        "category": "05_functional_cazy_kegg_Functional_CAZy_KEGG",
        "chunks": ["setup", "loading pkgs", "external BPstrain"],
        "figures": "No direct submitted figure",
        "purpose": "Summarizes dbCAN/CAZy GH counts for external B. pseudocatenulatum strains and draws the comparison plot.",
    },
]


def parse_chunks(text: str) -> dict[str, dict[str, object]]:
    chunks: dict[str, dict[str, object]] = {}
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        match = re.match(r"^```\{r\s*([^,}]*)", lines[i])
        if not match:
            i += 1
            continue
        name = match.group(1).strip()
        start_line = i + 1
        i += 1
        body = []
        while i < len(lines) and not lines[i].startswith("```"):
            body.append(lines[i])
            i += 1
        end_line = i + 1
        chunks[name] = {"body": "\n".join(body).strip() + "\n", "start": start_line, "end": end_line}
        i += 1
    return chunks


def main() -> None:
    text = RMD.read_text(encoding="utf-8-sig", errors="replace")
    chunks = parse_chunks(text)
    missing = sorted({chunk for target in TARGETS for chunk in target["chunks"] if chunk not in chunks})
    if missing:
        raise SystemExit(f"Missing expected chunk(s): {missing}")

    for target in TARGETS:
        folder = PUBLISH / target["category"]
        folder.mkdir(parents=True, exist_ok=True)
        parts = []
        for chunk_name in target["chunks"]:
            chunk = chunks[chunk_name]
            parts.append(
                "\n".join(
                    [
                        f"# ---- Source Rmd chunk: {chunk_name} (permonva_revise.Rmd lines {chunk['start']}-{chunk['end']}) ----",
                        str(chunk["body"]).rstrip(),
                        "",
                    ]
                )
            )
        header = "\n".join(
            [
                START,
                f"# Generated: {TODAY}",
                "# Source: C:/Users/zqr20/Documents/tjmeta/BIG_revision/permonva_revise.Rmd",
                f"# Category: {target['category']}",
                f"# Figure(s): {target['figures']}",
                f"# Purpose: {target['purpose']}",
                "# Scope note: extracted copy only; source Rmd was not modified.",
                END,
                "",
            ]
        )
        (folder / target["name"]).write_text(header + "\n".join(parts), encoding="utf-8")

    print(f"Extracted {len(TARGETS)} scripts from {RMD} into {PUBLISH}")


if __name__ == "__main__":
    main()
