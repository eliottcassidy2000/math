"use strict";
// A monomial with bases [b_1,...,b_d] is product b_j^w_j; 0^0 = 1.
const patterns = [
  [3,[[0,1],[0,2]]], [3,[[0,1],[1,2]]], [3,[[0,2],[1,2]]],
  [4,[[0,1],[2,3]]], [4,[[0,2],[1,3]]], [4,[[0,3],[1,2]]],
];
const formula = [
  [[1,1,1],[-1,0,1],[-1,1,0],[1,2,1],[2,2,2],[-1,1,2]],
  [[1,1,1],[-1,0,1],[-1,1,0],[1,0,0],[1,2,1],[1,1,2],[1,2,2]],
  [[1,1,1],[-1,0,1],[-1,1,0],[2,2,2],[-1,2,1],[1,1,2]],
  [[1,1,1,1],[-1,0,1,1],[1,2,1,1],[-1,1,1,0],[1,0,1,0],
   [-1,2,1,0],[1,1,1,2],[-1,0,1,2],[1,2,1,2]],
  [[1,1,1,1],[2,2,2,1],[-1,2,1,0],[2,1,2,2],[-1,0,1,2],[2,2,3,2]],
  [[1,1,1,1],[-1,1,0,1],[4,2,2,2],[-2,2,1,2],[1,1,2,1],[2,2,3,2]],
];
function canonical(raw) {
  const out = new Map();
  for (const [c,...bs] of raw) {
    const s = bs.join();
    out.set(s,(out.get(s)||0)+c);
  }
  return [...out].filter(x=>x[1]).sort();
}
function compile(n, es) {
  const raw = [];
  // Enumerate ordered lists of actually traversed reversed edges.
  const orders = [[],[es[0]],[es[1]],es,es.slice().reverse()];
  for (const ls of orders) {
    const R = ls.length+1;
    const low = [-1,...ls.map(e=>e[0])], high = [...ls.map(e=>e[1]),n];
    for (let mask=0; mask<R**n; mask++) {
      const a = Array.from({length:n},(_,i)=>Math.floor(mask/R**i)%R);
      if (ls.some(([u,v],i)=>a[v]!==i || a[u]!==i+1)) continue;
      if (a.some((r,i)=>i<low[r] || i>high[r])) continue;
      const demand = es.filter(([u,v])=>a[u]===a[v] &&
        !a.slice(u+1,v).some(r=>r===a[u]));
      // Inclusion-exclusion for removed forward adjacencies.
      for (let m=0; m<(1<<demand.length); m++) {
        let sign=1;
        const ok=Array.from({length:n-1},(_,j)=>Array.from({length:R},(_,r)=>low[r]<=j && high[r]>j));
        for (let i=0; i<demand.length; i++) if (m>>i&1) {
          sign=-sign;
          const [u,v]=demand[i];
          for (let j=u; j<v; j++) ok[j][a[u]]=false;
        }
        raw.push([sign,...ok.map(x=>x.filter(Boolean).length)]);
      }
    }
  }
  return canonical(raw);
}
for (let i=0; i<patterns.length; i++) {
  const got=compile(...patterns[i]), want=canonical(formula[i]);
  if (JSON.stringify(got)!==JSON.stringify(want)) throw Error(JSON.stringify({i,got,want}));
}
console.log("Six symbolic formulas verified");

const E = (1n<<60n)*((1n<<61n)-1n);
const cap = 157;
const mods = [65521,65519];
const expected = [[1,1],[0,0],[1,1],[0,0],[50,46],[0,0]];
const target = mods.map(m=>[-1n,1n].map(d=>Number((E+d)%BigInt(m))));
const pow = mods.map(m=>Array.from({length:4},(_,b)=>{
  const a=[1];
  for(let w=1;w<=cap;w++) a.push(a[w-1]*b%m);
  return a;
}));
function modValue(f,w,k) {
  const m=mods[k]; let val=0;
  for (const [cf,...bases] of f) {
    let z=cf;
    for (let j=0; j<bases.length; j++) z=z*pow[k][bases[j]][w[j]]%m;
    val+=z;
  }
  return (val%m+m)%m;
}
function exactValue(f,w) {
  let val=0n;
  for (const [cf,...bases] of f) {
    let z=BigInt(cf);
    for (let j=0;j<bases.length;j++) z*=BigInt(bases[j])**BigInt(w[j]);
    val+=z;
  }
  return val;
}
function smallDP([n,es],w) {
  const marks=[0];
  for (let j=0; j<w.length; j++) marks.push(marks[j]+w[j]+1);
  const N=1+marks[n-1], full=(1<<N)-1;
  const row=Array.from({length:N},(_,j)=>full^((1<<(j+1))-1));
  for (const [i,j] of es) {
    const u=marks[i],v=marks[j];
    row[u]^=1<<v; row[v]^=1<<u;
  }
  const d=new Uint32Array((full+1)*N);
  for(let i=0;i<N;i++)d[(1<<i)*N+i]=1;
  for(let s=1;s<=full;s++) for(let bits=s;bits;bits&=bits-1) {
    const i=31-Math.clz32(bits&-bits), q=d[s*N+i];
    if(q) for(let b=row[i]&(full^s);b;b&=b-1) {
      const t=b&-b,j=31-Math.clz32(t);
      d[(s|t)*N+j]+=q;
    }
  }
  let val=0;
  for(let i=0;i<N;i++)val+=d[full*N+i];
  return BigInt(val);
}
for(let i=0;i<6;i++) {
  const [n]=patterns[i], f=formula[i];
  for(let x=0;x<=2;x++) for(let y=0;y<=2;y++)
    for(let z=0;z<=(n===3?0:2);z++) {
      const w=n===3?[x,y]:[x,y,z];
      if(exactValue(f,w)!==smallDP(patterns[i],w)) throw Error("DP mismatch");
    }
  const first=[0,0],second=[0,0];
  const w=Array(n-1).fill(0);
  function sweep(j) {
    if(j===w.length) {
      const got=modValue(f,w,0);
      for(let k=0;k<2;k++) if(got===target[0][k]) {
        first[k]++;
        if(modValue(f,w,1)===target[1][k]) second[k]++;
      }
      return;
    }
    const L=(i===3&&j===1)?0:cap;
    for(let z=0;z<=L;z++){w[j]=z;sweep(j+1);}
  }
  sweep(0);
  if(first.join()!==expected[i].join()||second.some(k=>k)) throw Error("Sieve mismatch");
  console.log({pattern:patterns[i][1], first, both:second});
}
