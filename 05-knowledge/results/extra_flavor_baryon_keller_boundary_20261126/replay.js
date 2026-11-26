// node replay.js
"use strict";
const assert = {
 equal(a,b){if(a!==b) throw Error(String(a)+" != "+String(b));},
 deepEqual(a,b){if(a.size!==b.size||[...a].some(([k,v])=>b.get(k)!==v)) throw Error("Unequal polynomial maps");}
};
const fact = n => n < 2 ? 1n : BigInt(n)*fact(n-1);
function* tuples(m,n) {
  if(m===1){yield [n];return;}
  for(let i=0;i<=n;i++) for(const a of tuples(m-1,n-i)) yield [i,...a];
}
function* perms(a) {
  if(!a.length){yield [];return;}
  for(let i=0;i<a.length;i++) for(const b of perms(a.filter((_,j)=>j!==i))) yield [a[i],...b];
}
const key = a => a.join(",");
const unpack = k => k.split(",").map(Number);
function add(p,k,c) {
  const v=(p.get(k)||0n)+c;
  if(v) p.set(k,v); else p.delete(k);
}
function mul(p,q) {
  const r=new Map();
  for(const [s,a] of p) for(const [t,b] of q) {
    const x=unpack(s), y=unpack(t);
    add(r,key(x.map((v,i)=>v+y[i])),a*b);
  }
  return r;
}
function diff(p,q) {
  const r=new Map();
  for(const [s,a] of p) for(const [t,b] of q) {
    const x=unpack(s), y=unpack(t);
    if(y.some((v,i)=>v>x[i])) continue;
    const z=x.map((v,i)=>v-y[i]);
    const f=x.reduce((acc,v,i)=>acc*fact(v)/fact(z[i]),1n);
    add(r,key(z),a*b*f);
  }
  return r;
}
const scale = (p,c) => new Map([...p].map(([k,v])=>[k,v*c]).filter(([,v])=>v));
function minor(d,c) {
  const m=d+1,r=new Map(), cols=Array.from({length:m},(_,i)=>i).filter(i=>i!==c);
  for(const p of perms(cols)) {
    const x=Array(d*m).fill(0);
    let inv=c;
    for(let i=0;i<d;i++) {
      x[i*m+p[i]]=1;
      for(let j=0;j<i;j++) inv+=Number(p[j]>p[i]);
    }
    add(r,key(x),inv%2?-1n:1n);
  }
  return r;
}
for(const d of [2,3]){
  const m=d+1, M=Array.from({length:m},(_,i)=>minor(d,i)), cache=new Map();
  cache.set(key(Array(m).fill(0)),new Map([[key(Array(d*m).fill(0)),1n]]));
  let norms=0,lowerings=0;
  for(let n=0;n<=4;n++) for(const a of tuples(m,n)){
    const idx=a.findIndex(v=>v);
    const b=a.slice(); if(idx>=0) b[idx]--;
    const p=idx<0 ? cache.get(key(a)) : mul(cache.get(key(b)),M[idx]);
    cache.set(key(a),p);
    const norm=[...p].reduce((s,[k,c])=>s+c*c*unpack(k).reduce((x,v)=>x*fact(v),1n),0n);
    let H=1n,g=1n;
    for(let j=1;j<d;j++){ H*=fact(n+j)/fact(j);g*=BigInt(n+j); }
    assert.equal(norm,H*a.reduce((s,v)=>s*fact(v),1n));norms++;
    for(let j=0;j<m;j++){
      const expected=new Map(), b=a.slice(); b[j]--;
      assert.deepEqual(diff(p,M[j]),a[j]?scale(cache.get(key(b)),BigInt(a[j])*g):expected);
      lowerings++;
    }
  }
  console.log(`SU(${d}): ${norms} ambient Wick norms; ${lowerings} minor lowering identities`);
}
// Verify the count of spanning monomials in the central-fiber ideal
// and which of the trace powers actually survive. The argument for equality
// of spaces is the bigrading proof in the note, not this finite check.
for(const d of [2,3,4,5]) {
  const m=d+1;
  for(let a=0;a<=d+1;a++) for(let b=0;b<=2;b++) {
    const ideal=new Set();
    if(a>=d) for(const u of tuples(m,a-d)) for(const v of tuples(m,a+b-d+1))
      for(let i=0;i<m;i++) {
        const x=u.slice();x[i]++;ideal.add(key(x)+"|"+key(v));
      }
    const expected=a<d?0:Number(fact(m+a-d)/(fact(m-1)*fact(a-d+1)))*Number(fact(m+a+b-d)/(fact(m-1)*fact(a+b-d+1)));
    assert.equal(ideal.size,expected);
  }
  console.log(`d=${d}: central-fiber bigraded ideal and trace-nilpotency thresholds checked`);
}
