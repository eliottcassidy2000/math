// Finite checks of the digit bound, two-run examples, and endpoint distinction.
"use strict";
const assert = require("node:assert/strict");

function norm(xs, b) {
  const c = new Map();
  for (const [e, s] of xs) c.set(e, (c.get(e) || 0) + s);
  for (;;) {
    const hit = [...c].find(([, d]) => Math.abs(d) >= b);
    if (!hit) return [...c].filter(([,d]) => d);
    const [e, d] = hit, q = Math.trunc(d/b);
    c.set(e, d % b);
    c.set(e+1, (c.get(e+1) || 0) + q);
  }
}
let checked = 0;
for (let b = 2; b <= 7; b++) {
  const signed = Array.from({length:10},(_,i) => [i>>1, i&1 ? -1:1]);
  function rec(a, start) {
    const x = a.reduce((v,[e,s])=>v+s*b**e,0);
    if (x > 0) {
      const v = norm(a,b);
      assert.equal(v.reduce((z,[e,s])=>z+s*b**e,0),x);
      const e = Math.max(...v.map(([e])=>e));
      assert(x >= b**(e-a.length+1));
      checked++;
    }
    if (a.length < 4)
      for (let i = start; i < signed.length; ++i) rec([...a,signed[i]],i);
  }
  rec([],0);
}
function pathCount(a,b) {
  const n = a+b+1, apex = n-1;
  const out = Array(n).fill(0);
  for (let i=0;i<apex;i++) {
    for (let j=i+1;j<apex;j++) out[i] |= (1<<j);
    if (i<a) out[apex] |= (1<<i);
    else out[i] |= (1<<apex);
  }
  const dp=Array.from({length:1<<n},()=>Array(n).fill(0));
  for(let i=0;i<n;i++) dp[1<<i][i]=1;
  for(let s=1;s<1<<n;s++)
    for(let i=0;i<n;i++) if(dp[s][i])
      for(let j=0;j<n;j++) if(!(s&(1<<j)) && (out[i]&(1<<j)))
        dp[s|(1<<j)][j]+=dp[s][i];
  return dp.at(-1).reduce((a,b)=>a+b,0);
}
for(let k=2;k<=9;k++)
  for(let a=1;a<k;a++) assert.equal(pathCount(a,k-a), 3 - 2**a - 2**(k-a) + 2**k);
for(let k=2;k<=8;k++) assert.equal(pathCount(1,k),1+2**k);
for(let e of [2n,3n,5n]) {
  const h=(2n**e)-1n;
  for(let c of [2n,3n,5n]){
    const next=(2n**(e*c))-1n;
    assert(h**c<=next && next<(h+1n)**c);
    assert.equal(next+1n,(h+1n)**c);
  }
}
// Base 4 has two limiting phases: 2 is a repeated unit digit and sqrt(4).
for(let k=1n;k<8n;k++)
  assert.equal(2n**(2n*k+1n),4n**k+4n**k);
console.log({signedSums:checked,twoRunChecks:43,endpointChecks:9});
