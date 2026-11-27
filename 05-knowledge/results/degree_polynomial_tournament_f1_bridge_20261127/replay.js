"use strict";
const assert=require("node:assert/strict");
const R=n=>Array.from({length:n},(_,i)=>i);
const tr=A=>A[0].map((_,j)=>A.map(r=>r[j]));
const mul=(A,B)=>{const V=tr(B);return A.map(r=>V.map(c=>r.reduce((s,x,i)=>s+x*c[i],0)));};
const sum=a=>a.reduce((s,x)=>s+x,0);
const eq=(A,B)=>A.every((r,i)=>r.every((x,j)=>x===B[i][j]));
const I=(n,d=1)=>R(n).map(i=>R(n).map(j=>i===j?d:0));
const off2=G=>sum(G.map((r,i)=>sum(r.map((x,j)=>i===j?0:x*x))));
const sub=(A,B)=>A.map((r,i)=>r.map((x,j)=>x-B[i][j]));
function verifyParity(C){
  const n=C.length,G=mul(C,tr(C)),P=C.map((r,i)=>r.reduce((s,x,j)=>i===j?s:s*x,1));
  assert(n%4===2);
  assert(G.every((r,i)=>r.every((x,j)=>i===j?x===n-1:(x-(P[i]===P[j]?2:0))%4===0)));
  assert(off2(G)>=2*n*(n-2));
  return off2(G);
}
function rankMod(A,p){
  const n=A.length,M=A.map(r=>r.map(x=>(x%p+p)%p));let r=0;
  for(let j=0;j<n;j++){
    let i=r;while(i<n&&!M[i][j])i++;
    if(i===n)continue;
    [M[i],M[r]]=[M[r],M[i]];
    const inv=R(p).find(v=>v*M[r][j]%p===1);
    for(let k=r+1;k<n;k++){
      const b=M[k][j]*inv%p;
      for(let h=j;h<n;h++) M[k][h]=(M[k][h]-b*M[r][h]%p+p)%p;
    }
    ++r;
  }
  return r;
}
function paley(p){
  const sq=new Set(R(p-1).map(k=>(k+1)**2%p));
  const S=R(p).map(i=>R(p).map(j=>i===j?0:sq.has((j-i+p)%p)?1:-1));
  return [[0,...R(p).map(()=>1)],...S.map(r=>[-1,...r])];
}
for(const p of [3,7,11,19,23,59]){
  const C=paley(p),n=p-1;
  assert(eq(mul(C,tr(C)),I(p+1,p)));
  assert.equal(rankMod(C,p),(p+1)/2);
  const X=C.slice(0,n).map(r=>r.slice(0,n)),G=mul(X,tr(X));
  const B=C.slice(n).map(r=>r.slice(0,n));
  assert(eq(mul(B,tr(B)),I(2,p-1)));
  assert(eq(G,sub(I(n,p),mul(tr(B),B))));
  assert.equal(verifyParity(X),2*n*(n-2));
  assert.equal(rankMod(X,p),(p+1)/2);
  console.log("Paley",p,"deleted Gram energy",off2(G),"mod-p rank",rankMod(X,p));
}
function oneBack(n,a,b){
  return R(n).map(i=>R(n).map(j=>i===j?0:(i<j?1:-1)*((i===a&&j===b||i===b&&j===a)?-1:1)));
}
function reroot(A,v){
  const n=A.length,L=R(n-1).map(k=>(v+1+k)%n);
  return [[0,...L.map(()=>1)],...L.map(i=>[-1,...L.map(j=>A[v][i]*A[v][j]*A[i][j])])];
}
function hp(S){
  const n=S.length,e=(1<<n)-1; 
  function f(u,U){if(U===e)return 1;let x=0;for(let j=0;j<n;j++)if(!(U>>j&1)&&S[u][j]>0)x+=f(j,U|1<<j);return x;}
  return sum(R(n).map(u=>f(u,1<<u)));
}
for(const n of [6,8,10,12,14,16]){
  const X=oneBack(n,1,3),Y=reroot(X,2);
  assert(eq(Y,oneBack(n,1,n-1)));
  const del=A=>A.slice(1).map(r=>r.slice(1));
  assert.equal(hp(del(X)),3);
  assert.equal(hp(del(Y)),1+2**(n-3));
  console.log("Same degree polynomial, tournament order",n-1,"counts",hp(del(X)),hp(del(Y)));
}
const P=R(6).map(i=>[0,1,2,3,4,5][i]);
for(let m=0;m<1<<15;m++){
 const C=I(6,0);let x=m;
 for(const i of P)for(let j=i+1;j<6;j++){
   C[i][j]=(x&1)?-1:1;C[j][i]=-C[i][j];x>>=1;
 }
 verifyParity(C);
}
console.log("Checked the parity bound in all 32,768 labelled 6 by 6 examples.");

function fromBack(n,pairs){
 const A=oneBack(n,-1,-1);
 for(const [i,j] of pairs){A[i][j]=-1;A[j][i]=1;}
 return A;
}
function evenPolynomial(A){
 const G=mul(A,tr(A)),n=Math.floor(A.length/2);
 const e=[1],p=[0];let X=I(A.length);
 for(let k=1;k<=n;k++){
   X=mul(X,G);p.push(sum(X.map((r,j)=>r[j]))/2);
   e.push(sum(R(k).map(j=>(j%2?-1:1)*e[k-1-j]*p[j+1]))/k);
 }
 return e;
}
for(const [arcs,h] of [
  [[[0,2],[0,4],[0,6]],31],
  [[[0,2],[0,6],[1,3]],41]
]){
 const A=fromBack(7,arcs),C=[[0,...R(7).map(()=>1)],...A.map(r=>[-1,...r])];
 assert(A[6][0]>0);
 assert.equal(hp(A),h);
 assert.equal(evenPolynomial(A).join(","),"1,21,99,71");
 assert.equal(evenPolynomial(C).join(","),"1,28,182,252,49");
}
console.log("Strong pointed pair: both spectral polynomials match, Hamiltonian counts 31 and 41.");
