"use strict";
const assert=require("node:assert/strict");
const R=n=>Array.from({length:n},(_,i)=>i);
function *perms(a){if(!a.length){yield [];return;}for(let i=0;i<a.length;i++)for(const b of perms(a.filter((_,j)=>j!==i)))yield [a[i],...b];}
function compile(a,F,U){
 if(F.length>=31)throw Error("implementation uses safe integer trits");
 for(let j=0;j<U.length;j++)for(let h=j+1;h<U.length;h++)assert(a[U[j]][U[h]]);
 const P=R(F.length+1).map(i=>3**i),rows=[];
 function build(p){
  const L=U.map(v=>P.map((_,s)=>s===0||a[p[s-1]][v]));
  const D=U.map(v=>P.map((_,s)=>(s===p.length||a[v][p[s]])?1:2));
  const canEmpty=P.map((_,s)=>s===0||s===p.length||a[p[s-1]][p[s]]);
  const cache=R(U.length+1).map(()=>new Map());
  function put(k,z,i){
   const c=Math.floor(z/P[i])%3;
   return c===0&&!L[k][i]?-1:z+(D[k][i]-c)*P[i];
  }
  function count(k,z){
   if(cache[k].has(z))return cache[k].get(z);
   let w=0n;
   if(k===U.length){
    w=1n;
    for(let i=0;i<P.length;i++){
     const c=Math.floor(z/P[i])%3;
     if(c===2||(c===0&&!canEmpty[i])){w=0n;break;}
    }
   }else for(let i=0;i<P.length;i++){const t=put(k,z,i);if(t>=0)w+=count(k+1,t);}
   cache[k].set(z,w);return w;
  }
  return {p,put,count,H:count(0,0)};
 }
 for(const p of perms(F))rows.push(build(p));
 const H=rows.reduce((a,b)=>a+b.H,0n);
 function unrank(r){
  if(r<0n||r>=H)throw Error("rank out of range");
  let row;
  for(const e of rows){if(r<e.H){row=e;break;}r-=e.H;}
  const list=P.map(()=>[]);let z=0;
  for(let k=0;k<U.length;k++){
   let done=false;
   for(let i=0;i<P.length;i++){
    const t=row.put(k,z,i);if(t<0)continue;
    const n=row.count(k+1,t);
    if(r<n){list[i].push(U[k]);z=t;done=true;break;}r-=n;
   }
   assert(done);
  }
  return list.flatMap((v,i)=>i<row.p.length?[...v,row.p[i]]:v);
 }
 return {H,unrank};
}
function tour(n,B){const A=R(n).map(i=>R(n).map(j=>i<j));for(const [i,j] of B){A[j][i]=true;A[i][j]=false;}return A;}
function check(A,F){
 const V=R(A.length), t=compile(A,F,V.filter(i=>!F.includes(i)));
 const paths=[...perms(V)].filter(p=>p.slice(1).every((x,j)=>A[p[j]][x]));
 assert.equal(t.H,BigInt(paths.length));
 assert.deepEqual(paths.map((_,j)=>t.unrank(BigInt(j)).join(",")).sort(),paths.map(x=>x.join(",")).sort());
}
let cases=0;
for(let m=0;m<1024;m+=3){const B=R(5).flatMap(j=>R(j).map(i=>[i,j])).filter((_,i)=>m>>i&1);check(tour(5,B),[...new Set(B.flat())]);cases++;}
for(let m=0;m<4096;m+=37){const B=R(4).flatMap(j=>[0,1,2].map(i=>[i,j+3])).filter((_,i)=>m>>i&1);check(tour(7,B),[2,0,1]);cases++;}
const T=tour(80,[[0,79]]), table=compile(T,[0,79],R(78).map(i=>i+1));
assert.equal(table.H,1n+2n**78n);
for(const i of [0n,1n,table.H/2n,table.H-1n]){const p=table.unrank(i);assert.equal(new Set(p).size,80);assert(p.slice(1).every((v,j)=>T[p[j]][v]));}
assert.equal(compile(tour(8,[[0,3],[0,4],[0,7]]),[0],R(7).map(i=>i+1)).H,63n);
const cross=(p,r,q)=>1n+2n**BigInt(p+r+1)+2n**BigInt(q+r+1)+2n**BigInt(p+q+1)*3n**BigInt(r);
function v(x,p){let c=0;while(x%p===0n){x/=p;c++;}return c;}
function decode(H){if(H<=1n)return null;let u=H-1n,a=v(u,2n);u-=2n**BigInt(a);if(!u)return null;let b=v(u,2n);u-=2n**BigInt(b);if(!u)return null;let r=v(u,3n),p=a-r-1,q=b-r-1;return 1<=r&&r<p&&p<q&&cross(p,r,q)===H?[p,r,q]:null;}
let triple=0;
for(let q=3;q<=24;q++)for(let p=2;p<q;p++)for(let r=1;r<p;r++){assert.deepEqual(decode(cross(p,r,q)),[p,r,q]);triple++;}
for(let q=3;q<=6;q++)for(let p=2;p<q;p++)for(let r=1;r<p;r++){
 const E=[0,p+1,p+r+2,p+r+q+3],U=R(E[3]+1).filter(i=>!E.includes(i));
 assert.equal(compile(tour(E[3]+1,[[E[0],E[2]],[E[1],E[3]]]),E,U).H,cross(p,r,q));
}
console.log("Exact path unranking checked on "+cases+" small labeled inputs, plus an 80-vertex near-transitive control.");
console.log("Three-arc witness has H=63. Two-arc crossing decoder replayed at "+triple+" triples.");
module.exports={compile,tour,decode};
