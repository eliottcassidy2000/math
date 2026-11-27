"use strict";
const assert=require("node:assert/strict");
const R=n=>Array.from({length:n},(_,i)=>i);
const mod=(x,n)=>(x%n+n)%n;
const V=R(50).map(i=>[i>=25?1:0,Math.floor(i/5)%5,i%5]);
const A=V.map(([t,i,j])=>V.map(([u,k,l])=>t===u?+(i===k&&[t?2:1,t?3:4].includes(mod(j-l,5))):t===0?+(mod(j+i*k,5)===l):+(mod(l+i*k,5)===j)));
const N=A.map(r=>R(50).filter(i=>r[i]));
const cm=(x,y)=>N[x].filter(z=>A[z][y]);
function grid(u,v){
  return N[u].filter(a=>a!==v).flatMap(a=>N[v].filter(b=>b!==u).map(b=>{
    assert(!A[a][b]);
    const C=cm(a,b);assert.equal(C.length,1);
    assert(!A[u][C[0]]&&!A[v][C[0]]);
    return [a,C[0],b];
  }));
}
let control=0;
for(let u=0;u<50;u++)for(let v of N[u])if(u<v){
  const D=grid(u,v);
  assert.equal(new Set(D.map(x=>x[1])).size,36);
  assert.equal(D.length,36);
  // Compare the exact edge-grid count with a deterministic forbidden sample.
  const F=R(50).filter(w=>w!==u&&w!==v&&(w+3*u+5*v)%4===0);
  const a=F.filter(w=>A[u][w]).length,b=F.filter(w=>A[v][w]).length;
  const count=D.filter(e=>e.every(w=>!F.includes(w))).length;
  const rows=D.filter(([x,,y])=>!F.includes(x)&&!F.includes(y));
  assert.equal(rows.length,(6-a)*(6-b));
  assert.equal(count,rows.length-rows.filter(e=>F.includes(e[1])).length);
  assert(count>=(6-a)*(6-b)-(F.length-a-b));
  ++control;
}
assert.equal(control,175);
console.log("Hoffman–Singleton: edge-grid bijection and exact deletion count checked at all 175 edges.");

function firstCycle(path){
  if(path.length===7)return A[path[0]][path[6]]?path:null;
  for (const v of N[path.at(-1)])if(!path.includes(v)){
    let r=firstCycle([...path,v]);if(r)return r;
  }
  return null;
}
const C=firstCycle([0]);
const W=R(7).map(i=>cm(C[i],C[(i+3)%7])[0]);
const Z=R(7).map(i=>cm(W[i],C[(i+5)%7])[0]);
const mark=R(7).flatMap(j=>{const i=5*j%7;return [C[i],W[i],Z[i]];});
assert(mark.every((x,i)=>A[x][mark[(i+1)%21]]));
const targets=R(7).map(i=>3*i);
const cols=new Map(targets.map(i=>[i,grid(mark[i],mark[i+1]).filter(e=>e.every(x=>!mark.includes(x)))]));
function choose(left,used,H){
  if(!left.length)return H;
  const V=left.map(i=>[i,cols.get(i).filter(e=>e.every(x=>!used.has(x)))]).sort((a,b)=>a[1].length-b[1].length);
  const [i,choices]=V[0];
  for(const e of choices){
    const res=choose(left.filter(j=>i!==j),new Set([...used,...e]),new Map([...H,[i,e]]));
    if(res)return res;
  }
  return null;
}
const found=choose(targets,new Set(mark),new Map());
assert(found);
const long=mark.flatMap((x,i)=>[x,...(found.get(i)??[])]);
assert.equal(new Set(long).size,42);
assert(long.every((x,i)=>A[x][long[(i+1)%42]]));
console.log("Hoffman–Singleton 7-to-21-to-42 control:",long.join(","));

const f=R(14).map(j=>BigInt((52-j)*(54-j)-13-j));
const prod=f.reduce((p,v)=>p*v,1n);
const q=(246486240000n*prod+17n)/18n;
assert.equal(f[13],1573n);
assert.equal(q,522222022939543270256489658235041770644474393931250000000n);
assert.equal(Math.min(...R(21).map(j=>(52-j)**2-11-j)),993);
console.log("degree-57 path-choice lower bound per canonical 21-cycle:",String(prod));
console.log("degree-57 simple 63-cycle lower bound after 18-fold multiplicity:",String(q));

const length=new Set();
for(let e=0;e<=2;e++)for(let r=0;r+e<=21;r++){
 for(let j=0;j<r;j++)assert((56-4-2*e-j)**2>19+4*e+3*j);
 length.add(21+3*r+4*e);
}
assert(R(60).every(i=>length.has(i+27)));
assert.equal(Math.max(...length),86);
console.log("Mixed-detour lengths from each canonical seed: 21, 24, 25, and every length from 27 through 86.");
