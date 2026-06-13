from __future__ import annotations

import csv
import json
import os
import re
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path


BASE = Path(r"C:\Users\zqr20\Documents\tjmeta\BIG_revision")
PUBLISH = BASE / "publish_Script"
FIG_DIR = BASE / "forNM_revision"
DOC = PUBLISH / "00_documentation"
TODAY = date.today().isoformat()

START = "# >>> PUBLISH_ANNOTATION_START"
END = "# <<< PUBLISH_ANNOTATION_END"
RMD_START = "# >>> RMD_EXTRACT_ANNOTATION_START"
RMD_END = "# <<< RMD_EXTRACT_ANNOTATION_END"


INPUT_PATTERNS = [
    ("load", re.compile(r"\bload\s*\(\s*['\"]([^'\"]+)['\"]")),
    ("readRDS", re.compile(r"\breadRDS\s*\(\s*['\"]([^'\"]+)['\"]")),
    ("read.delim", re.compile(r"\bread\.delim\s*\(\s*['\"]([^'\"]+)['\"]")),
    ("read.csv", re.compile(r"\bread\.csv\s*\(\s*['\"]([^'\"]+)['\"]")),
    ("read.table", re.compile(r"\bread\.table\s*\(\s*['\"]([^'\"]+)['\"]")),
    ("read_excel", re.compile(r"\bread_excel\s*\(\s*['\"]([^'\"]+)['\"]")),
    ("read_delim", re.compile(r"\bread_delim\s*\(\s*['\"]([^'\"]+)['\"]")),
    ("read_tsv", re.compile(r"\bread_tsv\s*\(\s*['\"]([^'\"]+)['\"]")),
    ("readr::read_tsv", re.compile(r"\breadr::read_tsv\s*\(\s*['\"]([^'\"]+)['\"]")),
    ("source", re.compile(r"\bsource\s*\(\s*['\"]([^'\"]+)['\"]")),
]

OUTPUT_PATTERNS = [
    ("ggsave", re.compile(r"\bggsave\s*\(\s*(?:filename\s*=\s*)?['\"]([^'\"]+)['\"]")),
    ("pdf", re.compile(r"\bpdf\s*\(\s*(?:file\s*=\s*)?['\"]([^'\"]+)['\"]")),
    ("png", re.compile(r"\bpng\s*\(\s*(?:filename\s*=\s*)?['\"]([^'\"]+)['\"]")),
    ("jpeg", re.compile(r"\bjpeg\s*\(\s*(?:filename\s*=\s*)?['\"]([^'\"]+)['\"]")),
    ("tiff", re.compile(r"\btiff\s*\(\s*(?:filename\s*=\s*)?['\"]([^'\"]+)['\"]")),
    ("svg", re.compile(r"\bsvg\s*\(\s*(?:filename\s*=\s*)?['\"]([^'\"]+)['\"]")),
    ("save", re.compile(r"\bsave\s*\([^#\n]*?\bfile\s*=\s*['\"]([^'\"]+)['\"]")),
    ("saveRDS", re.compile(r"\bsaveRDS\s*\([^#\n]*?,\s*['\"]([^'\"]+)['\"]")),
    ("write.csv", re.compile(r"\bwrite\.csv\s*\([^#\n]*?,\s*['\"]([^'\"]+)['\"]")),
    ("write.table", re.compile(r"\bwrite\.table\s*\([^#\n]*?\bfile\s*=\s*['\"]([^'\"]+)['\"]")),
]


PURPOSE_HINTS = [
    ("QC", "Generates quality-control summaries and QC figure/table outputs."),
    ("feeding", "Organizes feeding metadata and related infant feeding summary figures."),
    ("phenotype", "Builds phenotype, baseline, and anthropometry summaries for figures/tables."),
    ("reorganizemeta", "Cleans and reorganizes cohort metadata for downstream analyses."),
    ("everydistance", "Calculates microbiome distance/ordination summaries and PCoA figure outputs."),
    ("pcoa", "Builds PCoA ordination analysis outputs."),
    ("permonva", "Runs or organizes PERMANOVA association summaries."),
    ("maslin", "Processes MaAsLin2 differential-abundance results and heatmaps."),
    ("heatmap", "Creates differential-abundance or trajectory heatmap outputs."),
    ("transmission", "Analyzes maternal-infant strain transmission patterns and related outputs."),
    ("transmisson", "Analyzes maternal-infant strain transmission patterns and related outputs."),
    ("instrain", "Processes inStrain transmission calls and related summary figures/tables."),
    ("persistance", "Calculates persistence/retention summaries and persistence figure outputs."),
    ("cazy", "Processes CAZy functional profiles and Bifidobacterium functional figures/tables."),
    ("GH_GT", "Analyzes GH/GT functional profiles and their phenotype/transmission associations."),
    ("KEGG", "Organizes KEGG pathway/function summaries and outputs."),
    ("Bifidobacterium", "Analyzes Bifidobacterium-specific profiles and figure outputs."),
    ("compareBPsize", "Benchmarks pipeline/read-size effects and exports benchmarking plots."),
    ("prevalence", "Calculates prevalence tables from microbiome profile inputs."),
]


def rel(path: Path) -> str:
    try:
        return str(path.relative_to(PUBLISH))
    except ValueError:
        return str(path)


def normalize_path(raw: str, script_dir: Path) -> tuple[str, str, bool | None]:
    expanded = raw.replace("/", "\\")
    candidates: list[Path] = []
    if re.match(r"^[A-Za-z]:\\", expanded):
        candidates.append(Path(expanded))
    elif raw.startswith("~/") or raw.startswith("~\\"):
        candidates.append(Path.home() / raw[2:])
    else:
        candidates.extend([script_dir / raw, BASE / raw])
    exists = any(p.exists() for p in candidates)
    normalized = str(candidates[0]) if candidates else raw
    scope = "BIG_revision" if str(candidates[0]).lower().startswith(str(BASE).lower()) else "external_or_personal"
    return normalized, scope, exists


def infer_content(path_text: str, role: str) -> str:
    name = Path(path_text.replace("/", "\\")).name.lower()
    if any(x in name for x in ["metadata", "meta", "sample", "phenotype", "who", "clinical"]):
        return "sample/clinical metadata"
    if any(x in name for x in ["ps_", "phylo", "otu", "taxa", "tree"]):
        return "phyloseq/taxonomy/abundance object"
    if any(x in name for x in ["dist", "pcoa", "aitchison", "axis"]):
        return "distance/ordination result"
    if any(x in name for x in ["transmission", "transmissibillity", "instrain", "strainphlan", "share"]):
        return "strain transmission result"
    if any(x in name for x in ["cazy", "gh", "gt", "kegg", "ko", "bp"]):
        return "functional annotation/profile"
    if any(x in name for x in ["persist", "persistance", "retention", "release"]):
        return "persistence/retention result"
    if any(x in name for x in ["figure", ".pdf", ".png", ".jpg", ".tif"]):
        return "figure output" if role == "output" else "figure/image input"
    if any(x in name for x in ["stat", "summary", "table", ".csv", ".tsv", ".txt"]):
        return "tabular result"
    if any(x in name for x in [".rda", ".rds"]):
        return "serialized R object"
    return "analysis input/output"


def infer_purpose(script_name: str, category: str, figures: str, outputs: list[dict]) -> str:
    hay = script_name + " " + category
    for key, purpose in PURPOSE_HINTS:
        if key.lower() in hay.lower():
            return purpose
    if figures:
        return f"Generates or supports submitted figure(s): {figures}."
    if outputs:
        return "Generates supporting analysis outputs and intermediate files."
    return "Supporting analysis script with no direct submitted figure match."


def extract_records(script: Path, role: str, patterns: list[tuple[str, re.Pattern]]) -> list[dict]:
    rows = []
    lines = script.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    for line_no, line in enumerate(lines, start=1):
        if line.lstrip().startswith("#"):
            continue
        for call, pat in patterns:
            for m in pat.finditer(line):
                raw_path = m.group(1)
                norm, scope, exists = normalize_path(raw_path, script.parent)
                rows.append(
                    {
                        "script": script.name,
                        "relative_script": rel(script),
                        "category": script.parent.name,
                        "line": line_no,
                        "role": role,
                        "function": call,
                        "path_raw": raw_path,
                        "path_normalized": norm,
                        "scope": scope,
                        "exists_now": exists,
                        "extension": Path(raw_path).suffix.lower(),
                        "content_inferred": infer_content(raw_path, role),
                    }
                )
    return rows


def figure_ids_from_name(name: str) -> str:
    stem = Path(name).stem
    label_part = stem.split("__", 1)[1] if "__" in stem else stem
    if "no_figure" in label_part.lower():
        return "no direct figure"
    ids = re.findall(r"FigureS?\d+[A-Za-z]?", label_part)
    return "; ".join(dict.fromkeys(ids))


def table_ids_from_name(name: str) -> str:
    stem = Path(name).stem
    label_part = stem.split("__", 1)[1] if "__" in stem else stem
    labels = re.findall(r"Supplementary Table \d+", label_part)
    if "table1baseline" in label_part:
        labels.append("table1baseline")
    if "highquality_GHGTstat_summary" in label_part:
        labels.append("highquality_GHGTstat_summary")
    return "; ".join(dict.fromkeys(labels))


def compact(items: list[str], limit: int = 4) -> str:
    clean = [x for x in dict.fromkeys(items) if x]
    if not clean:
        return "none detected"
    if len(clean) <= limit:
        return "; ".join(clean)
    return "; ".join(clean[:limit]) + f"; plus {len(clean) - limit} more"


def main() -> None:
    DOC.mkdir(parents=True, exist_ok=True)
    scripts = sorted(p for p in PUBLISH.rglob("*.R") if "00_documentation" not in p.parts)
    submitted_figures = {p.name.lower(): p.name for p in FIG_DIR.glob("*") if p.suffix.lower() in {".pdf", ".png", ".jpg", ".jpeg", ".tif", ".tiff"}}

    manifest = []
    io_rows = []
    issue_rows = []

    for script in scripts:
        text = script.read_text(encoding="utf-8-sig", errors="replace")
        clean_text = re.sub(
            rf"\A{re.escape(START)}.*?{re.escape(END)}\r?\n?",
            "",
            text,
            flags=re.S,
        )
        clean_text = re.sub(
            rf"\A{re.escape(RMD_START)}.*?{re.escape(RMD_END)}\r?\n?",
            "",
            clean_text,
            flags=re.S,
        )
        figures = figure_ids_from_name(script.name)
        tables = table_ids_from_name(script.name)
        inputs = extract_records(script, "input", INPUT_PATTERNS)
        outputs = extract_records(script, "output", OUTPUT_PATTERNS)
        io_rows.extend(inputs + outputs)

        libraries = sorted(set(re.findall(r"^\s*(?:library|require)\s*\(\s*([A-Za-z0-9_.]+)", clean_text, flags=re.M)))
        purpose = infer_purpose(script.name, script.parent.name, figures, outputs)
        output_names = [Path(r["path_raw"].replace("/", "\\")).name for r in outputs]
        input_names = [Path(r["path_raw"].replace("/", "\\")).name for r in inputs]

        issues_for_script = []
        for rec in inputs:
            if rec["exists_now"] is False and not any(ch in rec["path_raw"] for ch in ["*", "$", "{"]):
                issues_for_script.append(("Missing input", rec["line"], f"Referenced input not found now: {rec['path_raw']}"))
            if "OneDrive - BGI Hong Kong Tech" in rec["path_raw"]:
                issues_for_script.append(("Personal path", rec["line"], "Input uses a personal OneDrive path."))
            if rec["path_raw"].startswith("~/") or rec["path_raw"].startswith("~\\"):
                issues_for_script.append(("Home-relative path", rec["line"], "Input uses a home-relative path."))
        for rec in outputs:
            if "OneDrive - BGI Hong Kong Tech" in rec["path_raw"]:
                issues_for_script.append(("Personal output path", rec["line"], "Output uses a personal OneDrive path."))
            if re.search(r"C:/Users/zqr20/Documents/tjmeta/(?!BIG_revision)", rec["path_raw"]):
                issues_for_script.append(("Output outside BIG_revision", rec["line"], f"Output may be outside revision folder: {rec['path_raw']}"))
        for name, count in Counter(output_names).items():
            if name and count > 1:
                issues_for_script.append(("Duplicate output", "", f"Output filename appears {count} times and may be overwritten: {name}"))

        for out in output_names:
            if re.match(r"Figure", out, flags=re.I) and out.lower() not in submitted_figures:
                issues_for_script.append(("Figure mismatch", "", f"Code output is not present in forNM_revision: {out}"))
            if out == "Figure 1C.pdf":
                issues_for_script.append(("Figure name spacing", "", "Code writes Figure 1C.pdf, while submitted folder uses Figure1C.pdf."))

        for kind, line, detail in issues_for_script:
            issue_rows.append(
                {
                    "script": script.name,
                    "relative_script": rel(script),
                    "category": script.parent.name,
                    "line": line,
                    "issue_type": kind,
                    "detail": detail,
                }
            )

        manifest.append(
            {
                "script": script.name,
                "relative_script": rel(script),
                "category": script.parent.name,
                "figure_ids": figures or "no direct figure",
                "table_ids": tables or "no direct table",
                "purpose": purpose,
                "libraries": "; ".join(libraries),
                "input_count": len(inputs),
                "output_count": len(outputs),
                "issue_count": len(issues_for_script),
                "main_inputs": compact(input_names),
                "main_outputs": compact(output_names),
            }
        )

        header = [
            START,
                f"# Generated: {TODAY}",
                f"# Category: {script.parent.name}",
                f"# Figure(s): {figures or 'No direct submitted figure matched from filename'}",
                f"# Table(s): {tables or 'No direct submitted table matched from filename'}",
                f"# Purpose: {purpose}",
            f"# Main input(s): {compact(input_names)}",
            f"# Main output(s): {compact(output_names)}",
            f"# Static warning(s): {len(issues_for_script)} potential issue(s); see publish_Script/ERROR_REPORT.md and static_issue_report.csv",
            "# Scope note: annotation added to this published copy only; original script body below is unchanged.",
            END,
            "",
        ]
        script.write_text("\n".join(header) + clean_text.lstrip("\ufeff"), encoding="utf-8")

    with (PUBLISH / "script_annotations_manifest.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(manifest[0].keys()))
        writer.writeheader()
        writer.writerows(manifest)

    with (PUBLISH / "input_output_inventory.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(io_rows[0].keys()))
        writer.writeheader()
        writer.writerows(io_rows)

    issue_fields = ["script", "relative_script", "category", "line", "issue_type", "detail"]
    with (PUBLISH / "static_issue_report.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=issue_fields)
        writer.writeheader()
        writer.writerows(issue_rows)

    by_cat = Counter(r["category"] for r in manifest)
    by_issue = Counter(r["issue_type"] for r in issue_rows)
    by_role = Counter(r["role"] for r in io_rows)

    annotations = [
        "# publish_Script annotations",
        "",
        f"Updated on {TODAY}. These annotations reflect the current manually sorted `publish_Script` folder. Original scripts in `new_script` were not changed.",
        "",
        "## Category summary",
        "",
    ]
    for cat, count in sorted(by_cat.items()):
        annotations.append(f"- {cat}: {count} script(s)")
    annotations.extend(["", "## Script annotations", ""])
    for row in manifest:
        annotations.append(f"- **{row['relative_script']}**")
        annotations.append(f"  - Figures: {row['figure_ids']}")
        annotations.append(f"  - Tables: {row['table_ids']}")
        annotations.append(f"  - Purpose: {row['purpose']}")
        annotations.append(f"  - Inputs: {row['main_inputs']}")
        annotations.append(f"  - Outputs: {row['main_outputs']}")
        annotations.append(f"  - Static issues: {row['issue_count']}")
    (PUBLISH / "ANNOTATIONS.md").write_text("\n".join(annotations) + "\n", encoding="utf-8")

    error_md = [
        "# Static issue report",
        "",
        f"Updated on {TODAY}. Static scan of the current published R copies only.",
        "",
        "## Summary",
        "",
        f"- Scripts scanned: {len(manifest)}",
        f"- Input references detected: {by_role.get('input', 0)}",
        f"- Output references detected: {by_role.get('output', 0)}",
        f"- Potential issues detected: {len(issue_rows)}",
        "",
        "## Issue types",
        "",
    ]
    for kind, count in sorted(by_issue.items()):
        error_md.append(f"- {kind}: {count}")
    error_md.extend(["", "## High-priority issue list", ""])
    high_priority = [r for r in issue_rows if r["issue_type"] not in {"Personal path", "Personal output path", "Home-relative path"}]
    for r in high_priority[:120]:
        loc = r["relative_script"] if not r["line"] else f"{r['relative_script']}:{r['line']}"
        error_md.append(f"- `{loc}` - **{r['issue_type']}**: {r['detail']}")
    if len(high_priority) > 120:
        error_md.append(f"- Plus {len(high_priority) - 120} additional high-priority rows in `static_issue_report.csv`.")
    error_md.extend(["", "## Portability note", ""])
    error_md.append("Personal OneDrive and home-relative paths are listed in `static_issue_report.csv`; these may work on your machine but are the largest reproducibility risk for a published script bundle.")
    (PUBLISH / "ERROR_REPORT.md").write_text("\n".join(error_md) + "\n", encoding="utf-8")

    summary = {
        "scripts_scanned": len(manifest),
        "io_references": len(io_rows),
        "issues": len(issue_rows),
        "categories": dict(sorted(by_cat.items())),
        "issue_types": dict(sorted(by_issue.items())),
    }
    (DOC / "publish_scan_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
