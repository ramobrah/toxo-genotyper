/* Electronic PCR-RFLP (no DOM). Cuts the nested PCR product of each RFLP locus in silico with the enzymes of
 * Su et al. 2006/2010 and names the banding pattern by comparison with the reference panel (Type I = GT1,
 * II = ME49, III = VEG). Needs the MLST core (core.js). Runs in the browser and in Node. */
(function (root) {
'use strict';
const M = typeof module !== 'undefined' && module.exports ? require('./core.js') : root.MLST;
const SET = { A: 'A', C: 'C', G: 'G', T: 'T', R: 'AG', Y: 'CT', S: 'CG', W: 'AT', K: 'GT', M: 'AC', B: 'CGT', D: 'AGT', H: 'ACT', V: 'ACG', N: 'ACGT' };

/** Recognition sites in seq. A site is 'definite' when every base certainly fits, 'possible' when an ambiguity code
 *  or N could fit (partial digest, mixed infection or an unread base). Returns top-strand cut positions. */
function findSites(seq, enzName, enz) {
  const out = [];
  const scan = (site, cut) => {
    const L = site.length;
    for (let i = 0; i + L <= seq.length; i++) {
      let ok = true, sure = true;
      for (let k = 0; k < L; k++) {
        const q = SET[seq[i + k]] || 'ACGT', s = SET[site[k]];
        let any = false, all = true;
        for (const c of q) { if (s.includes(c)) any = true; else all = false; }
        if (!any) { ok = false; break; }
        if (!all) sure = false;
      }
      if (ok) { const at = i + cut; if (at > 0 && at < seq.length) out.push({ enzyme: enzName, start: i, end: i + L, at, sure }); }
    }
  };
  scan(enz.site, enz.cut);
  const rc = M.revcomp(enz.site);
  if (rc !== enz.site) scan(rc, enz.revCut !== undefined ? enz.revCut : enz.site.length - enz.cut); // non-palindromic (type IIS) sites on the other strand
  return out;
}
function fragments(len, cuts) {
  const at = [...new Set(cuts.map(c => c.at))].sort((a, b) => a - b);
  const f = []; let prev = 0;
  for (const a of at) { f.push(a - prev); prev = a; }
  f.push(len - prev);
  return f.filter(x => x > 0).sort((a, b) => b - a);
}
function digest(seq, locus, cfg) {
  const sites = locus.enzymes.flatMap(e => findSites(seq, e, cfg.enzymes[e]));
  const sure = sites.filter(s => s.sure), maybe = sites.filter(s => !s.sure);
  return { len: seq.length, sites, sure, maybe, bands: fragments(seq.length, sure), bandsMax: maybe.length ? fragments(seq.length, sites) : null };
}
/** Same pattern on a gel: the same visible bands, each within ~5% (what 2.5–3% agarose can separate). */
function sameGel(a, b, minBand) {
  const va = a.filter(x => x >= minBand), vb = b.filter(x => x >= minBand);
  if (va.length !== vb.length) return false;
  return va.every((x, i) => Math.abs(x - vb[i]) <= Math.max(4, 0.05 * Math.max(x, vb[i])));
}
/** The nested PCR product cut from a sequence that contains both inner primer sites, forward primer on the left. */
function productOf(seq, locus) {
  const h = M.locatePrimers(seq, locus);
  if (!h.fwd || !h.rev || !h.product) return null;
  return { lo: Math.min(h.fwd.pos, h.rev.pos), hi: Math.max(h.fwd.end, h.rev.end), flip: h.fwd.strand === '-' };
}

const withPrimers = (seq, L) => L.fwd + seq.slice(L.fwd.length, seq.length - L.rev.length) + M.revcomp(L.rev);
const strainOf = n => n.split(' (')[0];
/** Reference patterns for each locus: group panel references by gel pattern and name the groups. */
function buildRflp(panel, cfg) {
  const loci = {};
  for (const L of cfg.loci) {
    const from = [L.marker, ...(L.alsoFrom || [])].find(m => panel.markers[m]);
    const mk = from && panel.markers[from];
    if (!mk) continue;
    const refs = [];
    for (const r of mk.refs) {
      const p = productOf(r.seq, L);
      if (p) refs.push({ name: r.name, ref: r, amp: r.seq.slice(p.lo, p.hi), flip: p.flip });
    }
    // references trimmed a few bases into the primers: complete the product ends from a reference that has them
    const tpl = refs.find(x => x.ref === mk.refs[0]) || refs[0];
    const minP = Math.min(L.fwd.length, L.rev.length);
    if (tpl) for (const r of mk.refs) {
      if (refs.some(x => x.ref === r)) continue;
      const al = M.align(tpl.amp, r.seq);
      const both = al.cols.filter(c => c[0] >= 0 && c[1] >= 0);
      if (!both.length) continue;
      const [a0, r0] = both[0], [a1, r1] = both[both.length - 1];
      if (a0 > minP || tpl.amp.length - 1 - a1 > minP) continue; // missing more than the primer: leave it out
      refs.push({ name: r.name, ref: r, amp: tpl.amp.slice(0, a0) + r.seq.slice(r0, r1 + 1) + tpl.amp.slice(a1 + 1), flip: tpl.flip, completed: a0 + tpl.amp.length - 1 - a1 });
    }
    for (const x of refs) {
      x.seq = withPrimers(x.flip ? M.revcomp(x.amp) : x.amp, L);
      x.dig = digest(x.seq, L, cfg);
    }
    refs.sort((a, b) => a.ref.order - b.ref.order);
    const groups = [];
    for (const x of refs) {
      let g = groups.find(g => sameGel(g.bands, x.dig.bands, cfg.minBand));
      if (!g) { g = { bands: x.dig.bands, members: [] }; groups.push(g); }
      g.members.push(x.name); x.group = g;
    }
    for (const g of groups) {
      const arch = Object.entries(cfg.archetypes).filter(([s]) => g.members.includes(s)).map(([, t]) => t);
      g.archetypes = arch;
      if (arch.length) g.label = arch.join(' or ');
      else {
        // non-archetypal: use a published name (u-1, u-2 ...) when a panel strain in this group has one
        const pub = g.members.map(m => cfg.published[strainOf(m)] && cfg.published[strainOf(m)][L.name]).filter(v => v && /^u-/.test(v));
        g.label = pub.length ? pub[0] : 'new';
        g.labelNote = pub.length ? `named ${pub[0]} in Su et al.` : 'not I, II or III';
      }
    }
    loci[L.name] = { ...L, refs, groups };
  }
  return { cfg, loci };
}

/** Rebuild the sample's nested product: the sample's own bases where it was sequenced, a reference outside. */
function sampleProduct(res, L, entry) {
  if (!entry) { for (const c of res.comps) { entry = L.refs.find(x => x.name === c.ref.name); if (entry) break; } }
  if (!entry) return null;
  const refAmp = entry.amp; // panel orientation of the marker the locus was built from
  let cons = res.consensus;
  let al = M.align(refAmp, cons.seq);
  // the consensus may come from another marker whose panel faces the other way (5′-SAG2 typed from alt.SAG2 reads)
  const alR = M.align(refAmp, M.revcomp(cons.seq));
  if (alR.score > al.score) { al = alR; cons = { seq: M.revcomp(cons.seq), cols: [...cons.cols].reverse() }; }
  const both = al.cols.map((c, k) => c[0] >= 0 && c[1] >= 0 ? k : -1).filter(k => k >= 0);
  if (!both.length) return null;
  const first = both[0], last = both[both.length - 1];
  let s = '', fromRef = 0, own = 0; const filled = [];
  al.cols.forEach(([ri, ci], k) => {
    const inside = k >= first && k <= last;
    if (!inside) { if (ri >= 0) { s += refAmp[ri]; fromRef++; filled.push(s.length - 1); } return; }
    if (ri >= 0 && ci >= 0) {
      const col = cons.cols[ci];
      if (col && col.status === 'indel-conflict') { s += refAmp[ri]; return; } // base in one read only: excluded, keep reference
      if (cons.seq[ci] === 'N') { s += refAmp[ri]; fromRef++; filled.push(s.length - 1); return; } // unread base: treat like an unsequenced one
      s += cons.seq[ci]; own++;
    } else if (ci >= 0) { if (!(cons.cols[ci] && cons.cols[ci].status === 'indel-conflict')) { s += cons.seq[ci]; own++; } }
    // ri >= 0 && ci < 0: deletion in the sample
  });
  let seq = s, fill = filled;
  if (entry.flip) { seq = M.revcomp(s); fill = filled.map(i => s.length - 1 - i); }
  // the product ends are the primers themselves; bases under the primers do not count as unsequenced
  const lo = L.fwd.length, hi = seq.length - L.rev.length;
  fill = fill.filter(i => i >= lo && i < hi);
  return { seq: withPrimers(seq, L), own, fromRef: fill.length, filled: new Set(fill), filler: entry.name };
}

function nameBands(R, L, bands) {
  const g = R.loci[L].groups.find(g => sameGel(g.bands, bands, R.cfg.minBand));
  return g ? g : null;
}
/** e-RFLP result for one allele-calling result and one RFLP locus. */
function rflpForResult(R, res, locusName) {
  const L = R.loci[locusName];
  if (!L || !res.consensus || !res.best) return null;
  const sp = sampleProduct(res, L);
  if (!sp) return { locus: locusName, label: '—', status: 'none', notes: ['The sequence does not overlap this PCR product.'] };
  const dig = digest(sp.seq, L, R.cfg);
  const notes = [];
  const g = nameBands(R, locusName, dig.bands);
  let label = g ? g.label : 'new', status = g ? 'ok' : 'novel';
  let groupMax = null;
  if (dig.bandsMax && !sameGel(dig.bands, dig.bandsMax, R.cfg.minBand)) {
    groupMax = nameBands(R, locusName, dig.bandsMax);
    const l2 = groupMax ? groupMax.label : 'new';
    label = l2 === label ? label : `${label} / ${l2}?`;
    status = 'review';
    notes.push(`${dig.maybe.length} site(s) depend on an ambiguous or unread base (${dig.maybe.map(s => s.enzyme + ' at ' + (s.start + 1)).join(', ')}): uncut ${g ? g.label : 'new pattern'}, cut ${l2}. A double peak here can mean a mixed infection.`);
  }
  // did the unsequenced part matter? refill it from every reference allele and see if the name changes
  if (sp.fromRef > 0) {
    const inFill = dig.sites.filter(s => { for (let i = s.start; i < s.end; i++) if (sp.filled.has(i)) return true; return false; });
    const labels = new Set();
    for (const grp of L.groups) {
      const rep = L.refs.find(x => x.name === grp.members[0]);
      const alt = sampleProduct(res, L, rep);
      if (!alt) continue;
      const gg = nameBands(R, locusName, digest(alt.seq, L, R.cfg).bands);
      labels.add(gg ? gg.label : 'new');
    }
    if (labels.size > 1) { status = 'incomplete'; notes.push(`${sp.fromRef} bp of this ${L.size} bp product were not sequenced, and the pattern there differs between references (${[...labels].join(', ')}). Not determined; sequence the whole ${L.name} PCR product or check the unread bases.`); label = 'ND'; }
    else notes.push(`${sp.fromRef} bp of the ${L.size} bp product were not sequenced or not readable and were taken from ${sp.filler}; every reference gives the same pattern there.`);
  }
  if (status === 'novel') notes.push('This banding pattern matches none of the reference strains.');
  const lit = g && g.labelNote ? g.labelNote : null;
  const typesOf = gr => gr ? (gr.archetypes.length ? gr.archetypes : [gr.label]) : ['new'];
  const types = status === 'incomplete' ? [] : [...new Set([...typesOf(g), ...(groupMax !== null || (dig.bandsMax && !sameGel(dig.bands, dig.bandsMax, R.cfg.minBand)) ? typesOf(groupMax) : [])])];
  return { locus: locusName, label, status, types, bands: dig.bands, bandsMax: dig.bandsMax, sites: dig.sites, product: sp, group: g, notes, lit };
}

/** Combine loci typed together, e.g. SAG2 (5′+3′): the type is the archetype(s) both patterns allow. */
function combine(R, parts, names) {
  if (parts.some(p => p && p.status === 'incomplete')) return { label: 'ND', status: 'incomplete', types: [] };
  if (parts.some(p => p && p.status !== 'none' && !p.group)) return { label: 'new', status: 'novel', types: ['new'] };
  const have = parts.map((p, i) => p && p.group ? i : -1).filter(i => i >= 0);
  if (!have.length) return { label: '—', status: 'none', types: [] };
  if (have.length < parts.length) { const p = parts[have[0]]; return { label: `${p.group.archetypes.join(' or ') || p.label} (${names[have[0]]} only)`, status: 'partial', types: p.group.archetypes.length ? p.group.archetypes : [p.label] }; }
  let set = null;
  for (const p of parts) { const a = new Set(p.group.archetypes); set = set ? new Set([...set].filter(x => a.has(x))) : a; }
  const review = parts.some(p => p.status === 'review');
  if (!set.size) return { label: 'new', status: 'novel', types: ['new'] };
  return { label: [...set].join(' or '), status: review ? 'review' : 'ok', types: [...set] };
}

/** All RFLP loci for one sample: results = allele-calling results of that sample keyed by marker. */
function rflpForSample(R, byMarker) {
  const out = {};
  for (const L of Object.values(R.loci)) {
    const from = [L.marker, ...(L.alsoFrom || [])].find(m => byMarker[m]);
    const x = from ? rflpForResult(R, byMarker[from], L.name) : null;
    if (x) x.fromMarker = from;
    out[L.name] = x;
  }
  for (const c of R.cfg.combined) out[c.name] = c.of.every(n => out[n] === null) ? null : combine(R, c.of.map(n => out[n]), c.of);
  return out;
}
/** Column order as in Su et al. 2010 Table 3. */
function columns(R) {
  const names = Object.keys(R.loci);
  const cols = [];
  for (const n of names) {
    const c = R.cfg.combined.find(c => c.of[0] === n);
    if (c) { cols.push(c.name); continue; }
    if (R.cfg.combined.some(c => c.of.includes(n))) continue;
    cols.push(n);
  }
  return cols;
}

/** Match a multilocus e-RFLP profile to the ToxoDB PCR-RFLP genotype table.
 *  loci: output of rflpForSample. Returns matching genotype numbers, or the closest ones when none match. */
function matchGenotype(R, loci) {
  const G = R.cfg.genotypes; if (!G) return null;
  const typed = G.columns.filter(c => loci[c] && loci[c].types && loci[c].types.length);
  if (!typed.length) return null;
  const fits = (col, v, T) => {
    if (!v || v === 'nd' || v === '-') return true;
    const want = col === 'SAG1' && v === 'II' ? ['II', 'III'] : [v]; // the table writes II for SAG1 types II and III
    return T.some(t => want.includes(t) || (t === 'new' && /^u-/.test(v)));
  };
  const scored = G.rows.map(r => {
    const mm = typed.filter(c => !fits(c, r[1 + G.columns.indexOf(c)], loci[c].types));
    return { id: r[0], mm };
  });
  const matches = scored.filter(x => !x.mm.length).map(x => x.id);
  const viaNew = typed.filter(c => loci[c].types.includes('new'));
  let closest = [];
  if (!matches.length) { const m = Math.min(...scored.map(x => x.mm.length)); closest = scored.filter(x => x.mm.length === m).slice(0, 5); }
  const n = matches.length;
  const label = n === 1 ? '#' + matches[0] : n > 1 && n <= 6 ? matches.map(i => '#' + i).join(', ') : n > 6 ? `${n} possible` :
    `none (closest ${closest.map(c => '#' + c.id).join(', ')}: ${closest[0].mm.length} locus${closest[0].mm.length > 1 ? 'i' : ''} differ)`;
  const notes = [];
  if (typed.length < G.columns.length) { const miss = G.columns.filter(c => !typed.includes(c)); notes.push(typed.length <= miss.length ? `Only ${typed.join(', ')} typed (${typed.length} of ${G.columns.length} loci).` : `${typed.length} of ${G.columns.length} loci typed (missing: ${miss.join(', ')}).`); }
  if (viaNew.length && n) notes.push(`${viaNew.join(', ')}: a pattern not seen in the reference panel, taken as matching any non-archetypal (u-) allele in the table.`);
  if (!n) notes.push(`No ToxoDB genotype has this profile. Closest: ${closest.map(c => `#${c.id} (differs at ${c.mm.join(', ')})`).join('; ')}.`);
  return { typed, matches, closest, label, notes, status: n === 1 && typed.length === G.columns.length && !viaNew.length ? 'ok' : n ? 'review' : 'novel' };
}

const API = { matchGenotype, buildRflp, rflpForResult, rflpForSample, digest, findSites, sameGel, productOf, columns };
if (typeof module !== 'undefined' && module.exports) module.exports = API; else root.RFLP = API;
})(typeof window !== 'undefined' ? window : globalThis);
