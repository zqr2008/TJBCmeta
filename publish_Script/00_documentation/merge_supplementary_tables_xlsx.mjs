import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const forNmDir = "C:/Users/zqr20/Documents/tjmeta/BIG_revision/forNM_revision";
const outPath = path.join(forNmDir, "All_Supplementary_Tables_merged_with_academic_titles.xlsx");

const titles = {
  1: "Supplementary Table 1. Model selection for group-based trajectory modeling of infant weight-for-length/height z-score trajectories",
  2: "Supplementary Table 2. Multivariable associations between infant anthropometric growth trajectories and perinatal or early-life covariates",
  3: "Supplementary Table 3. Pairwise comparisons of global PCoA axes across maternal and infant microbiome visits and growth trajectories",
  4: "Supplementary Table 4. Visit-stratified pairwise comparisons of PCoA scores across infant growth trajectories",
  5: "Supplementary Table 5. Correlations between visit-stratified PCoA axes and maternal clinical or metabolic variables",
  6: "Supplementary Table 6. Maternal gut microbiome PERMANOVA associations with infant growth trajectories and maternal covariates across pregnancy visits",
  7: "Supplementary Table 7. Infant gut microbiome PERMANOVA associations with infant growth trajectories and early-life covariates across postnatal visits",
  8: "Supplementary Table 8. Genus-level MaAsLin2 associations in maternal gut microbiomes across pregnancy visits and growth trajectories",
  9: "Supplementary Table 9. Sensitivity analyses of differential-abundance associations across growth trajectories in maternal and infant microbiomes",
  10: "Supplementary Table 10. Species-level MaAsLin2 associations in maternal gut microbiomes across pregnancy visits and growth trajectories",
  11: "Supplementary Table 11. Species-level maternal-offspring strain transmission differences across infant growth trajectories",
  12: "Supplementary Table 12. One-time offspring strain persistence rates across infant growth trajectories",
  13: "Supplementary Table 13. High-persistence offspring strain differences across infant growth trajectories",
  14: "Supplementary Table 14. Pairwise comparisons of glycoside hydrolase gene counts in high-quality Bifidobacterium pseudocatenulatum genomes across infant growth trajectories",
  15: "Supplementary Table 15. Pairwise comparisons of glycoside hydrolase gene counts after inclusion of medium-quality Bifidobacterium pseudocatenulatum genomes",
  16: "Supplementary Table 16. Growth of Bifidobacterium pseudocatenulatum isolates on human milk oligosaccharides and control media with glycoside hydrolase gene counts",
  17: "Supplementary Table 17. Trajectory-associated differences in glycoside hydrolase family abundance among high-quality Bifidobacterium pseudocatenulatum genomes",
};

function parseDelimited(text, delimiter) {
  const rows = [];
  let row = [];
  let cell = "";
  let quoted = false;
  const clean = text.replace(/^\ufeff/, "");

  for (let i = 0; i < clean.length; i += 1) {
    const ch = clean[i];
    const next = clean[i + 1];
    if (quoted) {
      if (ch === '"' && next === '"') {
        cell += '"';
        i += 1;
      } else if (ch === '"') {
        quoted = false;
      } else {
        cell += ch;
      }
    } else if (ch === '"') {
      quoted = true;
    } else if (ch === delimiter) {
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

  while (rows.length && rows.at(-1).every((value) => value === "")) {
    rows.pop();
  }

  const width = Math.max(1, ...rows.map((r) => r.length));
  return rows.map((r) => [...r, ...Array(width - r.length).fill("")]);
}

function sheetNumber(name) {
  const match = name.match(/^Supplementary Table (\d+)/);
  return match ? Number(match[1]) : Number.POSITIVE_INFINITY;
}

function columnWidth(rows, colIndex) {
  const sample = rows.slice(0, 1200);
  const maxLen = Math.max(
    8,
    ...sample.map((row) => String(row[colIndex] ?? "").length),
  );
  return Math.max(70, Math.min(260, Math.round(maxLen * 7 + 18)));
}

const files = (await fs.readdir(forNmDir, { withFileTypes: true }))
  .filter((entry) => entry.isFile())
  .map((entry) => entry.name)
  .filter((name) => /^Supplementary Table \d+.*\.(csv|tsv)$/i.test(name))
  .sort((a, b) => sheetNumber(a) - sheetNumber(b) || a.localeCompare(b));

const seen = new Set();
const numberedFiles = [];
for (const file of files) {
  const number = sheetNumber(file);
  if (!seen.has(number)) {
    seen.add(number);
    numberedFiles.push({ number, file });
  }
}

const expected = Array.from({ length: 17 }, (_, idx) => idx + 1);
const missing = expected.filter((n) => !seen.has(n));
if (missing.length) {
  throw new Error(`Missing Supplementary Table file(s): ${missing.join(", ")}`);
}

const workbook = Workbook.create();
const report = [];

for (const { number, file } of numberedFiles) {
  const fullPath = path.join(forNmDir, file);
  const text = await fs.readFile(fullPath, "utf8");
  const delimiter = file.toLowerCase().endsWith(".tsv") ? "\t" : ",";
  const rows = parseDelimited(text, delimiter);
  const title = titles[number];
  if (!title) {
    throw new Error(`Missing academic title for Supplementary Table ${number}`);
  }
  const sheetName = `Supplementary Table ${number}`;
  const sheet = workbook.worksheets.add(sheetName);
  const titledRows = [[title, ...Array(rows[0].length - 1).fill("")], ...rows];
  sheet.showGridLines = true;
  sheet.getRangeByIndexes(0, 0, titledRows.length, titledRows[0].length).values = titledRows;
  sheet.freezePanes.freezeRows(2);

  const used = sheet.getRangeByIndexes(0, 0, titledRows.length, titledRows[0].length);
  used.format.wrapText = false;
  used.format.font = { name: "Arial", size: 10 };
  const titleRow = sheet.getRangeByIndexes(0, 0, 1, rows[0].length);
  titleRow.format.font = { name: "Arial", size: 11, bold: true };
  titleRow.format.fill = { color: "#E2F0D9" };
  titleRow.format.wrapText = true;
  const header = sheet.getRangeByIndexes(1, 0, 1, rows[0].length);
  header.format.font = { name: "Arial", size: 10, bold: true };
  header.format.fill = { color: "#D9EAF7" };

  for (let col = 0; col < rows[0].length; col += 1) {
    sheet.getRangeByIndexes(0, col, titledRows.length, 1).format.columnWidthPx = col === 0 ? Math.max(260, columnWidth(rows, col)) : columnWidth(rows, col);
  }

  report.push({
    sheet: sheetName,
    title,
    source: file,
    rows: rows.length,
    cols: rows[0].length,
  });
}

for (const { sheet } of report) {
  await workbook.render({ sheetName: sheet, range: "A1:H21", scale: 1, format: "png" });
}

const firstInspect = await workbook.inspect({
  kind: "table",
  range: "Supplementary Table 1!A1:H21",
  include: "values",
  tableMaxRows: 21,
  tableMaxCols: 8,
});
console.log(firstInspect.ndjson);
console.log(JSON.stringify({ output: outPath, sheets: report }, null, 2));

const xlsx = await SpreadsheetFile.exportXlsx(workbook);
await xlsx.save(outPath);
