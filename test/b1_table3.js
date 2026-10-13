// Rebuild each sample row of Zhu et al. 2023 (PLoS NTD) Table 3 on the Type X B1 backbone and check the tool's call
// against the paper's "MLST Type" column. Uncertain cells ("S/G") are tested with each alternative.
const M = require('../src/core.js'), fs = require('fs'), path = require('path');
const root = path.join(__dirname, '..');
const P = JSON.parse(fs.readFileSync(path.join(root, 'src/primers.json')));
const panel = M.buildPanel(JSON.parse(fs.readFileSync(path.join(root, 'src/panel.json'))), P);
const KM = panel.markers.B1.refs.find(r => r.name === 'Type X'); const backbone = KM.seq, OFF = 171;
const pos = [256, 330, 360, 366, 368, 428, 486, 495, 504, 513, 533, 545, 609];
const rows = [
  ['76,111,118,…', 'A T C T C T T G G T A C A', 'I'], ['77', 'A T C T C T T G G T A C G', 'I variant A'],
  ['87, 129', 'A T C T C T T G S T A C A', 'II/III'], ['92, 113, 354', 'A T C T C T T G C T A C A', 'X'],
  ['128', 'A T C T T T T G C T A C A', 'X var B'], ['151', 'A T C T C T T G S/G T/Y A C R/G', 'II/III variant A'],
  ['155', 'A T C T C T T G S T A C R', 'II/III variant A'], ['159', 'A T C T C T T G G T A C R', 'I variant B'],
  ['161', 'A T C T C T T G G T A C R', 'I variant B'], ['163', 'A T C T/Y C T T G G T A/M C R/G', 'I variant B'],
  ['188', 'R Y C T C Y Y R G T A C G', 'I variant C'], ['204', 'A T C Y C T T G S T M C A', 'II/III variant B'],
  ['279', 'A T C C C T T G G T C C/Y A', 'Atypical 1 variant B'], ['296', 'A T C C C T T G C T A C A', 'Atypical 2'],
];
const expand = cells => cells.reduce((acc, c) => acc.flatMap(a => c.split('/').map(x => [...a, x])), [[]]);
let n = 0, agree = 0; const byRow = {};
for (const [id, cellStr, paperType] of rows) {
  for (const alt of expand(cellStr.split(' '))) {
    const s = backbone.split(''); pos.forEach((p, k) => s[p - OFF - 1] = alt[k]);
    const seq = s.join('').slice(25, 470); // roughly the sequenced inner product
    const read = { fileName: 'S' + id, recName: 'S' + id, format: 'text', seq, qual: null, traces: null, ploc: null };
    const r = M.analyze([read], panel, {}, { ['S' + id]: { sample: id } }).results[0];
    const call = r.status === 'novel' ? `No exact match (closest ${r.best.ref.name}, ${r.best.events} diff)` : r.call;
    const norm = t => t.replace(/^Type /, '').replace(/var(iant)? /, 'variant ');
    const ok = norm(r.call) === norm(paperType) || (r.status === 'novel' && /variant/.test(paperType) && norm(paperType).startsWith(norm(r.best.ref.name)));
    n++; if (ok) agree++; (byRow[id] = byRow[id] || { uncertain: cellStr.includes('/'), ok: 0, n: 0 }); byRow[id].n++; if (ok) byRow[id].ok++;
    console.log((id + (cellStr.includes('/') ? ' [' + alt.join('') + ']' : '')).padEnd(32), call.padEnd(48), ('paper: ' + paperType).padEnd(30), ok ? 'agree' : 'DIFFERENT');
  }
}
const sure = Object.values(byRow).filter(r => !r.uncertain), unsure = Object.entries(byRow).filter(([, r]) => r.uncertain);
console.log(`\nRows with definite bases: ${sure.filter(r => r.ok === r.n).length}/${sure.length} agree with the paper.`);
console.log('Rows with uncertain bases (each reading tested): ' + unsure.map(([id, r]) => `${id}: ${r.ok}/${r.n} readings agree`).join('; '));
