"use strict";
const assert=require("node:assert/strict");
const {cores,dp,flip,F}=require("../overnight_chanpak_moore_20261126/commuting_unit_forks");

function weights(m){return Array.from({length:m},(_,j)=>j===0||j===m-1?1:2**(j-1));}
function start(m,incoming){
 const A=Array.from({length:m},(_,j)=>((1<<m)-1)^((1<<(j+1))-1));
 A[0]^=1<<(m-1);A[m-1]|=1;
 for(let j=0;j<m;j++)if(incoming>>j&1)A[j]|=1<<m;
 A.push(((1<<m)-1)^incoming);
 return A;
}
function direct(A){
 const n=A.length, u=n-1, [P,D]=dp(A),[Z]=dp(A,u), M=P.length-1;
 let R=0,V=0,Q=0;
 for(let s=0;s<=M;s++){
  Q+=P[s]*P[M^s];
  if(s>>u&1){R+=D[s*n+u]*P[M^s];V+=Z[s]*P[M^s];}
 }
 return {h:P[M],E:Z[M],e:D[M*n+u],R,V,Q};
}
function* realizations(g, extra){
 const m0=Math.max(3,3+Math.floor(Math.log2(g)));
 for(let m=m0;m<=m0+extra;m++){
  const w=weights(m);
  for(let j=1;j<1<<m;j++){
   let s=0;for(let k=j;k;k&=k-1)s+=w[31-Math.clz32(k&-k)];
   if(s===g)yield {A:start(m,j),m,j};
  }
 }
}
function index(a,b){
 // If either sign is reached, starting at an appropriate corner orients it as required.
 let z=a,t=b;while(t){const r=z%t;z=t;t=r;}if(z!==1)return null;
 let u=b%a;
 for(let i=1;i<=a;i++){
  u=(2*u)%a;if(u===2||u===a-2)return [i,u===2?2:-2];
  if(u===b%a)break;
 }
 return null;
}
function unitDeep(max,extra){
 const no=[], samples=[], tally=[];
 for(let g=1;g<=max;g++){
  let found=null;
  for(const {A,m,j} of realizations(g,([260,426].includes(g)&&extra>=2)?extra+1:extra)){
   const f=direct(A);if(f.e!==g)throw Error("terminal formula");
   for(let k=0;k<cores.length;k++){
    const [,,,H,R,E]=cores[k], a=H*f.h,b=H*f.R+R*f.Q;
    const q=index(a,b);
    if(q){found={g,m,j,k,q,h:f.h,a,b,d:H*f.E+R*f.V+E*g,A};break;}
   }
   if(found)break;
  }
  if(!found)no.push(g);else {if(samples.length<12||g===max)samples.push(found);tally.push(found);}
 }
 return {hits:tally.length,no,samples};
}

const extraCores=[
 [[2,72,3,21,39,15,61],5,[4,6],59,87,11],
 [[2,68,17,7,43,15,61],6,[3,4],71,72,14],
 [[34,20,9,67,13,30,55],3,[4,5],103,129,11],
 [[36,65,42,19,7,90,29],6,[0,4],143,173,23]
];
function localProfile(A,r){
 const n=A.length, [P,D]=dp(A),[Z]=dp(A,r),M=P.length-1;
 let R=0;for(let s=0;s<=M;s++)if(s>>r&1)R+=D[s*n+r]*P[M^s];
 return [P[M],Z[M],R];
}
function orbit(a,b){
 let t=1,i=0;const q=new Set();
 do{q.add(b*t%a);q.add((a-b*t%a)%a);t=2*t%a;i++;}while(t!==1);
 return [i,q.size];
}
function density(){
 const cs=[...cores,...extraCores], profiles=[];
 for(const [A,r,[u,v],H,R,E] of cs){
  assert.deepEqual(localProfile(A,r),[H,E,R]);
  assert.deepEqual(localProfile(flip(A,u,v),r),[H,E+2,R]);
  const seen=orbit(H,R);profiles.push([H,R,...seen]);
 }
 assert.deepEqual(profiles.map(q=>q[3]),[22,20,28,58,70,102,120]);
 return profiles;
}
function proofExample(){
 const t=unitDeep(7,2).samples[6];
 const N=BigInt(t.a)*F(1680n)+BigInt(t.b)*1681n-BigInt(t.d+14);
 const W=BigInt(t.a)*F(1688n)+BigInt(t.b)*1689n-BigInt(t.d+14);
 assert.equal(W,N-2n);
 const counts=[N,N+14n,W,W+14n].map(String);
 return {c:1681,counts};
}

const out={orbit:density(),example:proofExample(),census:unitDeep(512,2)};
assert.equal(out.census.hits,512);
assert.deepEqual(out.census.no,[]);
console.log(JSON.stringify(out,null,2));
