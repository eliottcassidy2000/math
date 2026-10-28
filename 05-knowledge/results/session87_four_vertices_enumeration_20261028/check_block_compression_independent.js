#!/usr/bin/env node
// Independent bit-vector/minimax implementation; no Python source imported.
// node check_block_compression_independent.js > independent-javascript-receipt.json
'use strict';
const need = (c, m) => { if (!c) throw new Error(m); };
const span = cols => cols.reduce((a, c) => [...new Set(a.concat(a.map(x => x ^ c)))], [0]);
const wt = (x, r, k) => {
  let w = 0;
  for (let i = 0; i < r; i++) if ((x >>> (i*k)) & ((1 << k)-1)) w++;
  return w;
};
const im = (cols, x) => cols.reduce((a, c, i) => a ^ (((x >>> i) & 1) ? c : 0), 0);
const spaces = [new Map([['0', [0]]])];
for (let d = 0; d < 3; d++) {
  const next = new Map();
  for (const sp of spaces[d].values()) for (let x = 1; x < 64; x++) {
    if (!sp.includes(x)) {
      const a = [...new Set(sp.concat(sp.map(t => t ^ x)))].sort((a,b) => a-b);
      next.set(a.join(','), a);
    }
  }
  spaces.push(next);
}
const good = spaces.map(level => [...level.values()].filter(sp =>
  sp.every(x => x === 0 || wt(x, 3, 2) >= 3)));
need(JSON.stringify(spaces.map(s => s.size)) === '[1,63,651,1395]', 'space counts');
need(JSON.stringify(good.map(s => s.length)) === '[1,27,36,0]', 'good counts');
const gl = [];
for (let a = 1; a < 4; a++) for (let b = 1; b < 4; b++) if (a !== b) gl.push([a,b]);
const diag = new Set(gl.flatMap(a => gl.map(b => [0,1,2,3]
  .map(x => x | (im(a,x) << 2) | (im(b,x) << 4)).sort((a,b) => a-b).join(','))));
need(good[2].every(x => diag.has(x.join(','))) && diag.size === 36, 'diagonal');
const T = m => {
  const sizes = (m%2 ? [3] : []).concat(Array(Math.floor((m-(m%2 ? 3 : 0))/2)).fill(2));
  let offset = 0;
  const cols = [];
  for (const n of sizes) {
    cols.push(...(n === 2 ? [2,3] : [2,4,3]).map(x => x << offset));
    offset += n;
  }
  need(span(cols).length === (1 << m), 'T invertible');
  need(span(cols.map((x,i) => x ^ (1 << i))).length === (1 << m), 'T+I invertible');
  return cols;
};
const spread = s => {
  if (s <= 3) return [[1,2]];
  const m = s-2, t = T(m);
  return Array.from({length: 1 << m}, (_,a) =>
    [a | (1 << m), im(t,a) | (1 << (m+1))]).concat(spread(m));
};
const max = s => ((1 << s) - (s%2 ? 5 : 1))/3;
const construction = [];
for (let s = 2; s <= 10; s++) {
  const lines = spread(s), used = lines.flatMap(b => span(b).filter(Boolean));
  need(new Set(used).size === used.length && used.length === 3*max(s), 'packing');
  construction.push({s, lines: lines.length,
    holes: Array.from({length: (1 << s)-1}, (_,i) => i+1).filter(x => !used.includes(x))});
}
const pop = x => { let n = 0; while (x) { x &= x-1; n++; } return n; };
const adaptive = (r,k) => {
  const states = [];
  for (let i = 0; i < r; i++) for (let a = 1; a < (1 << k); a++) states.push(a << (k*i));
  need(states.length < 31, '32-bit subset mask range');
  const queries = states.map(w => states.reduce((q,x,i) => q | (pop(x&w)%2 ? (1 << i) : 0), 0));
  const memo = new Map();
  const solve = mask => {
    if (pop(mask) <= 1) return 0;
    if (memo.has(mask)) return memo.get(mask);
    let best = states.length;
    for (const q of queries) {
      const a = mask&q, b = mask&~q;
      if (a && b) best = Math.min(best, 1 + Math.max(solve(a), solve(b)));
    }
    memo.set(mask,best);
    return best;
  };
  const ans = solve((1 << states.length)-1);
  need(ans === (r === 1 ? k : r*k-1), 'adaptive');
  return {r,k,depth: ans,nontrivial_states: memo.size};
};
const ad = [];
for (let r = 1; r <= 6; r++) ad.push(adaptive(r,2));
for (let r = 1; r <= 3; r++) ad.push(adaptive(r,3));
for (let r = 1; r <= 2; r++) ad.push(adaptive(r,4));
const report = {
  engine: 'independent JavaScript in functions runtime; no Python execution',
  subspace_counts: spaces.map(s => s.size),
  valid_kernel_counts: good.map(s => s.length),
  twisted_diagonal_count: diag.size,
  constructions: construction,
  nine_lines: spread(5).map(b => span(b).filter(Boolean).sort((a,b) => a-b)),
  minima: Array.from({length: 10}, (_,i) => {
    const r = i+1;
    for (let s = 2; s <= 6; s++) if (max(s) >= r) return {r,s};
  }),
  adaptive_minimax: ad
};
console.log(JSON.stringify(report, null, 2));
