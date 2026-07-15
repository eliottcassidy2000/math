#!/usr/bin/env node
"use strict";

// Execute the exact public complement-line API and atlas-indexing source from
// the browser explorer in a DOM-free VM.  This catches drift between the HTML
// implementation and the committed address atlas without duplicating the API.

const fs = require("fs");
const vm = require("vm");

const htmlPath = "03-artifacts/visualizations/tournament-tiling-explorer.html";
const atlasPath = "05-knowledge/results/tournament_tiling_metagraph_address_codex_S4.json";
const html = fs.readFileSync(htmlPath, "utf8");
const data = JSON.parse(fs.readFileSync(atlasPath, "utf8"));

function segment(start, end) {
  const a = html.indexOf(start);
  const b = html.indexOf(end, a);
  if (a < 0 || b < 0) throw new Error(`source anchor missing: ${start} / ${end}`);
  return html.slice(a, b);
}

const apiSource = segment("function addressSizeDataFor", "function addressPrefix");
const indexSource = segment("    addressAtlas=data;", "    buildOverview();onBitsChanged();");
const context = {data, window: {}};
vm.createContext(context);
vm.runInContext(
  `let N=5,addressAtlas=null,flowAtlas=null;\n${apiSource}\n${indexSource}\n` +
  "globalThis.audits=mergedMetagraphSupportedSizes().map(auditComplementLineAPI);",
  context
);

const expectedLines = {3: 1, 4: 4, 5: 32, 6: 512};
let failures = 0;
for (const audit of context.audits) {
  const mismatch = audit.lineCount !== expectedLines[audit.n];
  failures += audit.failures.length + Number(mismatch);
  process.stdout.write(
    `n=${audit.n} lines=${audit.lineCount} blue=${audit.blueLines} ` +
    `black=${audit.blackLines} failures=${audit.failures.length + Number(mismatch)}\n`
  );
  if (audit.failures.length) process.stdout.write(`  ${audit.failures.join(",")}\n`);
}
process.stdout.write(`TOTAL FAILURES: ${failures}\n`);
if (failures) process.exitCode = 1;
