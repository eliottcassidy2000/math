"use strict";
const assert=require("node:assert/strict");
const {cores,flip,dp,seam}=require("../overnight_chanpak_moore_20261126/commuting_unit_forks");

function inflate(P,pos,n){
 const blocks=[];let m=0;
 for(let i=0;i<P.length;i++){let a=[];for(let j=0;j<(i===0?n:1);j++)a.push(m++);blocks.push(a);}
 const T=new Array(m).fill(0);
 for(let i=0;i<P.length;i++)for(const u of blocks[i]){
  for(let j=0;j<P.length;j++)if(P[i]>>j&1)for(const v of blocks[j])T[u]|=1<<v;
  for(const v of blocks[i])if(u<v)T[u]|=1<<v;
 }
 return [T,blocks[pos][0]];
}
function pr(A,u){
 const [P,D]=dp(A),[Z]=dp(A,u),n=A.length,M=P.length-1;
 let R=0;for(let s=0;s<=M;s++)if(s>>u&1)R+=D[s*n+u]*P[M^s];
 return [P[M],Z[M],R];
}
function compile(P,r,i=-1,j=-1){
 const p=i<0?P:flip(P,i,j);
 const e=[1,2,3].map(n=>pr(...inflate(p,r,n)));
 const coef=e[0].map((x,i)=>[x,(e[1][i]-x)/2]);
 if(e[2].some((x,i)=>x!==coef[i][0]+6*coef[i][1]))throw Error("not linear");
 return coef;
}
const S=[46,12,56,48,35,2],T=[110,8,122,112,99,66,2],r=4;
function probe(){
 const out=[];
 for(const P of [S,T]){
  const J=compile(P,r), ed=[];
  for(let j=2;j<P.length;j++)for(let i=1;i<j;i++){
   const d=compile(P,r,i,j).map((row,k)=>row.map((x,l)=>x-J[k][l]));
   ed.push({i,j,d});
  }
  out.push({J,ed});
 }
 return out;
}
function join(A,u,B,v){
 const n=A.length,r=n+v,T=A.concat(B.map(x=>x<<n));
 for(let y=n;y<T.length;y++)for(let x=0;x<n;x++){
  if(x===u&&y===r)T[x]|=1<<y;else T[y]|=1<<x;
 }return T;
}
function distinctStrong(A){
 const M=(1<<A.length)-1;
 for(let i=0;i<A.length;i++){let s=1<<i,o=0;while(s!==o){o=s;for(let j=0;j<A.length;j++)if(o>>j&1)s|=A[j];}
 if(s!==M)return false;}
 return true;
}
function hp(A){const [P]=dp(A);return P[P.length-1];}
function forward(P){return P.map((x,i)=>P.reduce((y,z,j)=>y|((z>>i&1)<<j),0));}
const out=probe();const j0=out[0].ed.find(x=>x.i===3&&x.j===5);
const j1=out[1].ed.find(x=>x.i===3&&x.j===5);
assert.deepEqual(j0.d,[[0,0],[-2,0],[-2,0]]);
assert.deepEqual(j1.d,[[-2,0],[-2,0],[-2,0]]);
function hullProbe(){
 const [a,u,[x,y]]=cores[0], rows=[];
 const left=[[6,4,0],[4,5,0],[0,5,1]];
 // Root 0 takes the three consecutive positions in an ordered three-vertex host.
 for(const C of [forward(a),forward(flip(a,x,y))]){
  const [h,e,R]=pr(forward(C),u);
  const N=left.map(A=>hp(join(A,0,C,u)));
  assert.deepEqual(N,[h-e+R,h+2*R,h+4*R]);
  rows.push([N,[h,e,R]]);
 }
 assert.deepEqual(rows,[[[39,71,119],[23,8,24]],[[37,71,119],[23,10,24]]]);
 return rows;
}
function historicalControls(){
 const x=[];
 for(const n of [1,2,3,4]){
  const p=[];
  for(const K of [S,T]){
   const A=inflate(K,r,n),B=inflate(flip(K,3,5),r,n);
   const v=[A,B].map(([M,u])=>join(M,u,[0],0));
   assert(v.every(distinctStrong));
   p.push(...v.map(hp));
  }
  const N=30*2**n-1;assert.deepEqual(p,[N,N,N,N-2]);
  x.push([n,p]);
 }
 return x;
}
function graphControls(){
 const bank=[...cores];
 const examples=[
  {H:289,T:[0],j:2,c:3,pad:0,Eflip:1},
  {H:889,T:[0],j:0,c:3,pad:2,Eflip:0},
  {H:9011,T:[14,28,8,16,5],j:0,c:1,pad:2,Eflip:0}
 ];
 return examples.map(({H,T,j,c,pad,Eflip})=>{
  const [X,u,[a,b]]=bank[j],n=T.length+pad,mask=(1<<n)-1;
  const Z=Array.from({length:pad},(_,k)=>mask^((1<<(k+1))-1)).concat(T.map(y=>y<<pad));
  const B=seam(Z,n-1,c,0);
  const G=[X,flip(X,a,b)].map(L=>join(L,u,B,B.length-1));
  const P=G.map(hp);assert.equal(P[Eflip],H);assert.equal(P[0]-P[1],2);assert(G.every(distinctStrong));
  return {H,count:P,order:G[0].length};
 });
}
console.log(JSON.stringify({profiles:[j0.d,j1.d],probes:hullProbe(),repairs:historicalControls(),exactTargets:graphControls()},null,2));
