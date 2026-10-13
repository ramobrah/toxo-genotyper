// Re-analyse Battisti et al. 2018 (Vet Parasitol 253:43-47) GenBank sequences MG587956-MG588017 and compare with their Table 2.
const M = require('../src/core.js'), fs = require('fs'), path = require('path');
const root = path.join(__dirname, '..');
const mode = process.argv[2] || 'panel';   // 'panel' = built-in primers; 'nested' = the study's inner primers for SAG1 and alt.SAG2
const P = JSON.parse(fs.readFileSync(path.join(root, 'src/primers.json')));
if (mode === 'nested') { P['alt.SAG2'] = P['alt.SAG2'].nested; P.SAG1 = { fwd: 'CAATGTGCACCTGTAGGAAGC', rev: 'GTGGTTCTCCGTCGGTGTGAG' }; }
const panel = M.buildPanel(JSON.parse(fs.readFileSync(path.join(root, 'src/panel.json'))), P);
const recs = M.parseFasta(fs.readFileSync(path.join(root, 'data/study/seqs.fasta'), 'utf8'));
const reads = [], overrides = {};
for (const r of recs) {
  const acc = r.header.split(' ')[0];
  const m = r.header.match(/isolate (\S+?)_/); const iso = m ? m[1] : '817/11 (unlabelled BTUB record)';
  reads.push({ fileName: acc, recName: acc, format: 'text', seq: r.seq, qual: null, traces: null, ploc: null, header: r.header });
  overrides[acc] = { sample: iso + ' ' + acc }; // one row per GenBank record
}
const out = M.analyze(reads, panel, {}, overrides);
// RFLP-type view of each call: which canonical strain shares the called allele?
const typeOf = (res) => {
  if (!res.consensus) return 'fail';
  const mk = panel.markers[res.marker];
  const names = [res.best.ref.name, ...(res.sameAs || [])];
  const t = []; if (names.includes('GT1')) t.push('I'); if (names.includes('ME49')) t.push('II'); if (names.includes('VEG')) t.push('III');
  return t.length ? t.join('/') : 'other (' + res.call + ')';
};
// Table 2 of the paper
const T = `12/11 I I na na na na|116/11 I I na na na na|324/11 na I na na na na|335/11 na III na na na II|339/11 I I na na na na|697/11 III na na na na na|698/11 I III na na na na|699/11 I III na na na na|700/11 I I na na na na|723/11 na I na na na na|730/11 I I na na na na|766/11 III na na na na na|780/11 I I na na na na|783/11 na I na na na na|806/11 I I na na na na|807/11 na I na na na na|817/11 II II II/III II II II|846/11 III I na na na na|847/11 I III na na na na|869/11 na I na na na na|876/11 II I II/III na III II|890/11 I na na na na na|A4 na I na na na na|A7 na I na na na na|B8 na I na na na na|C5 na I na na na na|P51 na I na na na na|P54 na I na na na na|P55 I I na na na na|P60 I I na na na na|Q64 III II na na na na|Q67 I I na na na na|Q68 na I na na na na|Q69 na I na na na na|Q70 na I na na na na|Q72 I I na na na na|Q73 I I na na na na|6 I I na na na na|41 I I na na na na|60 I I na na na na|63 III na na na III na|64 na II na na na na|91 na I na na na na`;
const cols = ['alt.SAG2', 'GRA6', 'SAG1', 'BTUB', '5SAG2', 'c22-8'];
const paper = {}; for (const row of T.split('|')) { const [id, ...v] = row.split(' '); paper[id] = Object.fromEntries(cols.map((c, i) => [c, v[i]])); }
const rows = [];
for (const res of out.results) {
  const acc = res.files.join(','); const hdr = recs.find(r => r.header.startsWith(res.files[0])).header;
  const declared = /5'SAG2/.test(hdr) ? '5SAG2' : /altSAG2/i.test(hdr) ? 'alt.SAG2' : /gra6/i.test(hdr) ? 'GRA6' : /SAG1/.test(hdr) ? 'SAG1' : /C22-8/i.test(hdr) ? 'c22-8' : /tubulin/.test(hdr) ? 'BTUB' : '?';
  const id = res.sample.split(' ')[0];
  const pub = paper[id] ? paper[id][declared] : '(not in Table 2)';
  const ours = typeOf(res);
  const agree = ours === pub ? 'yes' : 'NO';
  rows.push({ diffs: res.best ? res.best.events : null, acc, id, declared, assigned: res.marker, call: res.call, status: res.status, ours, pub, agree, flags: res.flags.filter(f => f.level !== 'info').map(f => f.text).join(' | ') });
}
rows.sort((a, b) => a.declared.localeCompare(b.declared) || a.id.localeCompare(b.id, undefined, { numeric: true }));
if (process.argv[3] === 'json') { console.log(JSON.stringify(rows)); return; }
for (const r of rows) console.log([r.declared.padEnd(8), r.id.padEnd(8), r.acc.padEnd(11), ('→' + r.assigned).padEnd(10), r.call.slice(0, 22).padEnd(23), r.status.padEnd(7), (r.ours + (r.status === 'novel' ? ' (nearest, ' + r.diffs + ' diff)' : '')).padEnd(26), ('paper:' + r.pub).padEnd(14), r.agree, r.flags.slice(0, 110)].join(' '));
const tally = {}; for (const r of rows) { const k = r.declared; tally[k] = tally[k] || { n: 0, agree: 0, ok: 0 }; tally[k].n++; if (r.agree === 'yes') tally[k].agree++; if (r.status === 'ok') tally[k].ok++; if (r.status === 'novel') tally[k].novel = (tally[k].novel || 0) + 1; }
console.log('\nSummary by marker (agree with paper / called cleanly / total):'); for (const [k, v] of Object.entries(tally)) console.log(' ', k.padEnd(8), v.agree + '/' + v.n, 'same type as paper (nearest reference),', v.ok + '/' + v.n, 'exact matches,', (v.novel || 0), 'with 1+ differences');
