"use strict";
const assert = require("node:assert/strict");
const R=n=>Array.from({length:n},(_,i)=>i);
const choose2=n=>n*(n-1)/2;
function occupancy(edges,slots) {
  const q=Math.floor(edges/slots),r=edges%slots;
  return slots*choose2(q)+r*q;
}
function data(k,h) {
  return [choose2(k)*h*(h-k+2),(57-k)*occupancy(k*h,56)];
}
assert.deepEqual(data(14,13),[1183,9030]);
assert.deepEqual(data(14,15),[4095,12642]);
assert.equal(data(9,28)[1]-data(9,28)[0],336);
const k8=[7,9,11,13,15,17,19,21].map(h=>[h,data(8,h)[0]-data(8,h)[1]]);
assert.deepEqual(k8,[[7,196],[9,-28],[11,-28],[13,196],[15,252],[17,140],[19,252],[21,588]]);
const lower=k=>Math.ceil((data(k,28)[1]-data(k,28)[0])/(2*k-5));
assert.deepEqual([9,14,15].map(lower),[26,427,541]);
assert.deepEqual([14,15].map(k=>Math.ceil(lower(k)/28)),[16,20]);
console.log("14 x 13: missing and forced outside pairs",data(14,13).join(", "));
console.log("8 branches, half-size and slack:",JSON.stringify(k8));
console.log("Sign-preserving direct edge minima, k=9/14/15:",[9,14,15].map(lower).join(", "));

function sizes(A) {
  const seen=new Set(), out=[];
  for(let i=0;i<A.length;i++) if(!seen.has(i)){
    const visit=[i]; seen.add(i);
    for(const a of visit) for(const b of A[a]) if(!seen.has(b)){
      seen.add(b); visit.push(b);
    }
    out.push(visit.length);
  }
  return out.sort((a,b)=>a-b);
}
function plane(q) {
  const V=[];
  for(let x=0;x<q;x++)for(let y=0;y<q;y++)V.push([x,y,1]);
  for(let x=0;x<q;x++)V.push([x,1,0]);
  V.push([1,0,0]);
  const has=(l,p)=>l.reduce((s,v,i)=>s+v*p[i],0)%q===0;
  const O=[0,0,1], atInfinity=V.filter(p=>has(O,p));
  const rays=V.filter(l=>has(l,O));
  const corners=rays.map(l=>atInfinity.find(p=>has(l,p)));
  const P=V.filter(p=>p.some((v,i)=>v!==O[i])&&!has(O,p) && rays.findIndex(l=>has(l,p))!==q);
  const L=V.filter(l=>l.some((v,i)=>v!==O[i]) && !has(l,O) && corners.findIndex(p=>has(l,p))!==q);
  const pg=P.map(p=>rays.findIndex(l=>has(l,p)));
  const lg=L.map(l=>corners.findIndex(p=>has(l,p)));
  assert.equal(P.length,q*(q-1)); assert.equal(L.length,P.length);
  const N=P.map(p=>L.map(l=>+has(l,p)));
  const A=R(P.length).map(()=>[]), B=R(L.length).map(()=>[]);
  for(let i=0;i<N.length;i++)for(let j=0;j<N.length;j++){
    const dot=N[i].reduce((s,v,k)=>s+v*N[j][k],0);
    const dual=N.reduce((s,r)=>s+r[i]*r[j],0);
    assert(dot===(i===j?q-1:0)||pg[i]!==pg[j]&&dot===1);
    assert(dual===(i===j?q-1:0)||lg[i]!==lg[j]&&dual===1);
    if(pg[i]!==pg[j]&&!dot) A[i].push(j);
    if(lg[i]!==lg[j]&&!dual) B[i].push(j);
  }
  for(const [X,G] of [[A,pg],[B,lg]]) {
    assert(X.every((a,i)=>a.length===q-1&&new Set(a.map(j=>G[j])).size===q-1&&!a.some(j=>G[j]===G[i])));
    assert.deepEqual(sizes(X),Array(q-1).fill(q));
  }
  console.log(`Deleted PG(2,${q}): both leaves are ${q-1} copies of K${q}.`);
}
[3,5,7].forEach(plane);

// Only tests the abstract graph/PSD conditions, not 0/1 factorability.
const fake=R(182).map(x=>R(14).flatMap(j=>{
  const i=Math.floor(x/13),t=x%13;
  return i===j?[]:[13*j+(t+(i===0&&j===1?1:i===1&&j===0?12:0))%13];
}));
assert.deepEqual(sizes(fake),[182]);
console.log("An abstract connected matching leave passes automatic Gram PSD but is not an incidence certificate.");
