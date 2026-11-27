"use strict";
const assert = require("node:assert/strict");
const range = n => Array.from({length:n}, (_,i) => i);
const mod = (x,p) => (x%p+p)%p;
function matrix(n, adj) { return range(n).map(i=>range(n).map(j=>+adj(i,j))); }
const A5 = matrix(5,(i,j)=>[1,4].includes(mod(i-j,5)));
const P = range(5).flatMap(i=>range(i).map(j=>[i,j]));
const A10 = matrix(10,(i,j)=>P[i].every(x=>!P[j].includes(x)));
const V=range(50).map(i=>[+(i>=25),Math.floor(i/5)%5,i%5]);
const A50 = matrix(50,(a,b)=>{
 const [t,i,j]=V[a],[u,k,l]=V[b];
 return t===u ? i===k&&[t?2:1,t?3:4].includes(mod(j-l,5)) :
 t===0 ? mod(j+i*k,5)===l : mod(l+i*k,5)===j;
});
const C42=[0,30,9,29,4,34,8,7,39,10,41,43,3,33,31,14,38,36,1,46,48,16,26,11,35,22,23,24,37,13,28,18,19,15,25,27,2,42,20,21,47,45];
const dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
const add=(a,b)=>a.map((x,i)=>x+b[i]);
function closure(A,C) {
 let S=new Set(C),counts=[S.size];
 for (;;) {
  const next=range(A.length).filter(v=>!S.has(v)&&[...S].filter(w=>A[v][w]).length>=2);
  if(!next.length)return counts;
  for(const v of next)S.add(v);
  counts.push(S.size);
 }
}
function check(A,C){
 const n=A.length,d=A[0].reduce((x,y)=>x+y,0),I=range(n);
 for(const v of I)for(const w of I)
  assert.equal(dot(A[v],A[w])+A[v][w],v===w?d:1);
 const D=x=>I.map(i=>2*dot(A[i],x)+x[i]);
 const E=I.flatMap(i=>I.filter(j=>i<j&&A[i][j]).map(j=>[i,j]));
 const F=C.map((x,i)=>[x,C[(i+1)%C.length]]);
 for (const [x,y] of F)assert(A[x][y]);
 const source=d===7?F:E;
 const g=source.flatMap(([a,b])=>source.filter(([c,e])=>new Set([a,b,c,e]).size===4&&
  !A[a][c]&&!A[a][e]&&!A[b][c]&&!A[b][e]).map(f=>[[a,b],f]))[0];
 if(d>2)assert(g);
 const vs=(g??[F[0]]).flatMap(([a,b])=>{
  const u=I.map(i=>+(i===a)-(i===b)),w=add(D(u),u);
  assert.deepEqual(D(w),w.map((x,i)=>x+(4*d-4)*u[i]));
  return [u,w];
 });
 for(const i of range(vs.length))for(const j of range(vs.length))
  assert.equal(dot(vs[i],vs[j]),i!==j?0:i%2===0?2:2*(4*d-4));
 const rot=A.map((r,i)=>{
  const L=I.filter(j=>r[j]),k=C.indexOf(i);
  if(k<0)return L;
  const pre=C[mod(k-1,C.length)],post=C[(k+1)%C.length];
  return [pre,post,...L.filter(x=>x!==pre&&x!==post)];
 });
 const seen=new Set(),faces=[];
 for(const a of I)for(const b of rot[a]){
  let u=a,v=b,len=0; if(seen.has(u*n+v))continue;
  do{
   assert(!seen.has(u*n+v));seen.add(u*n+v);len++;
   [u,v]=[v,rot[v][(rot[v].indexOf(u)+1)%d]];
  }while(u!==a||v!==b);
  faces.push(len);
 }
 assert(faces.includes(C.length)&&faces.every(x=>x>=5));
 const genus=(2-n+d*n/2-faces.length)/2;
 const floor=Math.ceil((2-n+d*n/2-1-Math.floor((d*n-C.length)/5))/2);
 assert(Number.isInteger(genus)&&genus>=floor);
 console.log(JSON.stringify({n,d,cycle:C.length,closure:closure(A,C),orthogonalPairs:g??null,genus,floor}));
}
check(A5,[0,1,2,3,4]);
function find5(A, path){
 if(path.length===5)return A[path[0]][path[4]]?path:null;
 for(const j of range(A.length))if(A[path.at(-1)][j]&&!path.includes(j)){
  const out=find5(A,[...path,j]);if(out)return out;
 }
 return null;
}
check(A10,find5(A10,[0]));
check(A50,C42);
const n=3250,s2=37*52**2+26*51**2;
assert.equal(37*52+26*51,n);
assert.equal(Math.floor((5*s2-n+272*n)/66),26047);
assert.equal(n*57/2-26047,66578);
assert.equal(1+Math.floor((n*57-63)/5),37038);
assert.equal(Math.ceil((2-n+n*57/2-37038)/2),26170);
for(const k of [21,63]){
 assert(!range(58).some(a=>(k*a-(50+7*k))%15===0));
}
console.log("Degree-57 conditional arithmetic checked.");
