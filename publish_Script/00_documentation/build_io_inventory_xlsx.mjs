import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const publishDir = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/publish_Script";
const outPath = path.join(publishDir, "script_input_output_inventory.xlsx");

function parseCsv(text) {
  const rows = [];
  let row = [];
  let cell = "";
  let quoted = false;
  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    const next = text[i + 1];
    if (quoted) {
      if (ch === '"' && next === '"') {
        cell += '"';
        i++;
      } else if (ch === '"') {
        quoted = false;
      } else {
        cell += ch;
      }
    } else if (ch === '"') {
      quoted = true;
    } else if (ch === ",") {
      row.push(cell);
      cell = "";
    } else if (ch === "\n") {
      row.push(cell.replace(/\r$/, ""));
      rows.push(row);
      row = [];
      cell = "";
    } else {
      cell += ch;
    }
  }
  if (cell.length || row.length) {
    row.push(cell.replace(/\r$/, ""));
    rows.push(row);
  }
  const header = rows.shift() ?? [];
  return rows
    .filter((r) => r.some((v) => v !== ""))
    .map((r) => Object.fromEntries(header.map((h, idx) => [h.replace(/^\ufeff/, ""), r[idx] ?? ""])));
}

function matrixFromObjects(rows, columns) {
  return [columns, ...rows.map((row) => columns.map((col) => row[col] ?? ""))];
}

function writeSheet(workbook, name, rows, columns, widths = []) {
  const sheet = workbook.worksheets.add(name);
  const matrix = matrixFromObjects(rows, columns);
  sheet.getRangeByIndexes(0, 0, matrix.length, columns.length).values = matrix;
  sheet.freezePanes.freezeRows(1);
  sheet.showGridLines = false;
  const header = sheet.getRangeByIndexes(0, 0, 1, columns.length);
  header.format.fill = { color: "#1F4E78" };
  header.format.font = { color: "#FFFFFF", bold: true };
  header.format.wrapText = true;
  const body = sheet.getRangeByIndexes(0, 0, matrix.length, columns.length);
  body.format.wrapText = true;
  body.format.borders = {
    insideHorizontal: { style: "continuous", color: "#D9E2F3" },
    insideVertical: { style: "continuous", color: "#D9E2F3" },
    edgeBottom: { style: "continuous", color: "#8EA9DB" },
  };
  widths.forEach((w, idx) => {
    sheet.getRangeByIndexes(0, idx, Math.max(matrix.length, 1), 1).format.columnWidthPx = w;
  });
  return sheet;
}

function countBy(rows, key) {
  const map = new Map();
  for (const row of rows) {
    const value = row[key] || "blank";
    map.set(value, (map.get(value) ?? 0) + 1);
  }
  return [...map.entries()].sort((a, b) => String(a[0]).localeCompare(String(b[0]))).map(([name, count]) => ({ name, count }));
}

const manifest = parseCsv(await fs.readFile(path.join(publishDir, "script_annotations_manifest.csv"), "utf8"));
const inventory = parseCsv(await fs.readFile(path.join(publishDir, "input_output_inventory.csv"), "utf8"));
const issues = parseCsv(await fs.readFile(path.join(publishDir, "static_issue_report.csv"), "utf8"));

const workbook = Workbook.create();

const categoryCounts = countBy(manifest, "category");
const roleCounts = countBy(inventory, "role");
const issueCounts = countBy(issues, "issue_type");
const existenceCounts = countBy(inventory, "exists_now");

const summaryRows = [
  { metric: "Scripts scanned", value: manifest.length, note: "Current R files in publish_Script only" },
  { metric: "Input/output references", value: inventory.length, note: "Literal read/load/write/save/figure paths detected from code" },
  { metric: "Input references", value: roleCounts.find((r) => r.name === "input")?.count ?? 0, note: "Detected inputs" },
  { metric: "Output references", value: roleCounts.find((r) => r.name === "output")?.count ?? 0, note: "Detected outputs" },
  { metric: "Potential issues", value: issues.length, note: "Static findings; see Static_Issues tab" },
  { metric: "Missing input references", value: issues.filter((r) => r.issue_type === "Missing input").length, note: "Literal inputs not found at scan time" },
  { metric: "Personal/home path warnings", value: issues.filter((r) => /Personal|Home-relative/.test(r.issue_type)).length, note: "Portability risks, not always runtime errors on your machine" },
];

writeSheet(workbook, "Summary", summaryRows, ["metric", "value", "note"], [220, 100, 520]);
const summary = workbook.worksheets.getItem("Summary");
summary.getRange("A1:C1").format.fill = { color: "#17365D" };
summary.getRange("A1:C20").format.wrapText = true;

writeSheet(
  workbook,
  "Script_Manifest",
  manifest,
  ["category", "script", "relative_script", "figure_ids", "table_ids", "purpose", "main_inputs", "main_outputs", "input_count", "output_count", "issue_count", "libraries"],
  [260, 260, 420, 150, 190, 420, 320, 320, 90, 90, 90, 360],
);

writeSheet(
  workbook,
  "Input_Output",
  inventory,
  ["category", "script", "line", "role", "function", "path_raw", "path_normalized", "exists_now", "scope", "extension", "content_inferred"],
  [260, 260, 70, 75, 100, 430, 500, 90, 150, 80, 220],
);

writeSheet(
  workbook,
  "Static_Issues",
  issues,
  ["category", "script", "relative_script", "line", "issue_type", "detail"],
  [260, 260, 420, 70, 170, 560],
);

writeSheet(workbook, "Category_Counts", categoryCounts, ["name", "count"], [420, 100]);
writeSheet(workbook, "Issue_Counts", issueCounts, ["name", "count"], [240, 100]);
writeSheet(workbook, "Path_Existence", existenceCounts, ["name", "count"], [160, 100]);

for (const sheetName of ["Input_Output", "Static_Issues"]) {
  const sheet = workbook.worksheets.getItem(sheetName);
  sheet.freezePanes.freezeRows(1);
}

const summaryInspect = await workbook.inspect({
  kind: "table",
  range: "Summary!A1:C8",
  include: "values",
  tableMaxRows: 10,
  tableMaxCols: 4,
});
console.log(summaryInspect.ndjson);

const issueInspect = await workbook.inspect({
  kind: "table",
  range: "Issue_Counts!A1:B20",
  include: "values",
  tableMaxRows: 20,
  tableMaxCols: 3,
});
console.log(issueInspect.ndjson);

await workbook.render({ sheetName: "Summary", range: "A1:C8", scale: 1, format: "png" });
await workbook.render({ sheetName: "Input_Output", range: "A1:K18", scale: 1, format: "png" });
await workbook.render({ sheetName: "Static_Issues", range: "A1:F18", scale: 1, format: "png" });

const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(outPath);
console.log(`Saved ${outPath}`);
