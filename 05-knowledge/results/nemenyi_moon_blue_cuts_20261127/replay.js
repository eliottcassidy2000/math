"use strict";
const a=require("node:assert/strict");
const bin=(n,k)=>{if(k<0||k>n)return 0;let c=1;for(let j=1;j<=k;j++)c=c*(n+1-j)/j;return c;};
const range=n=>Array.from({length:n},(_,i)=>i);
const mod=(x,n)=>(x%n+n)%n;
function datum(n,m){
 const I=range(n),bit=i=>(m>>i)&1;
 const A=I.map(i=>I.map(j=>i!==j&&((i<j)===(bit(i)===bit(j)))));
 const d=A.map(r=>r.filter(Boolean).length);
 if(d.includes(0)||d.includes(n-1))return null;
 const C=I.filter(i=>!bit(i)).concat(I.filter(bit));
 if(!C.every((x,i)=>A[x][C[(i+1)%n]]))return null;
 const D=C.map(i=>d[i]);
 const rots=I.filter(s=>I.every(i=>D[(i+s)%n]===D[i]));
 const refs=I.filter(r=>I.every(i=>D[mod(r-i,n)]===n-1-D[i]));
 const sc=!!refs.length,profile=range(n+1).map(k=>k<3?0:bin(n,k)-d.reduce((s,v)=>s+bin(v,k-1),0));
 const blue=(p)=>I.every(i=>I.every(j=>A[p[i]][p[j]]===A[p[n-1-j]][p[n-1-i]]));
 const cuts=I.map(i=>C.slice(i).concat(C.slice(0,i)));
 a.equal(cuts.filter(blue).length,rots.length*(sc?(n%2?1:2):0));
 const subsets=range(1<<n).map(m=>C.filter(v=>(m>>v)&1)).filter(S=>S.length>=3);
 const profile2=range(n+1).fill(0),rp=I.map(()=>range(n+1).fill(0));
 for(const S of subsets){
  const cyc=S.every((v,j)=>A[v][S[(j+1)%S.length]]);
  if(cyc){profile2[S.length]++; for(const v of S)rp[v][S.length]++;}
  else a(S.some(v=>S.every(w=>v===w||A[v][w])));
 }
 a.deepEqual(profile2,profile);
 for(const v of I)for(let k=3;k<=n;k++)
  a.equal(rp[v][k],bin(n-1,k-1)-bin(d[v],k-1)-I.filter(i=>A[i][v]).reduce((s,i)=>s+bin(d[i]-1,k-2),0));
 return {n,A,d,C,D,sc,profile,blue,rots};
}
function summary(T){
 const {n,A,d,D,sc,profile,blue,rots}=T,I=range(n);
 let H=0,b=0,close=0;
 function rec(p,m){
  if(!m){H++;b+=+blue(p);close+=+A[p.at(-1)][p[0]];return;}
  for(const i of I)if((m>>i)&1&&(!p.length||A[p.at(-1)][i]))rec([...p,i],m^(1<<i));
 }
 rec([], (1<<n)-1);
 return {H,blueTilings:b/rots.length,closingPaths:close,sc,cyclicDegrees:D,score:d.slice().sort((a,b)=>a-b),profile:profile.slice(3)};
}
let count=0,first=null;
for(let n=3;n<=8;n++){
 const seen=new Map();
 for(let m=0;m<1<<(n-1);m++){
  const T=datum(n,m);if(!T)continue;count++;
  const key=T.d.slice().sort((a,b)=>a-b).join(",");
  if(seen.has(key)&&seen.get(key)[0]!==T.sc&&!first)first=[n,seen.get(key)[1],m];
  if(!seen.has(key)) seen.set(key,[T.sc,m]);
 }
}
a.deepEqual(first,[8,12,18]);
const n=8;
const X=range(n).map(i=>[2,3].includes(i)?2:[1,4].includes(i)?5:8);
for(const [j,mask,r4] of [[2,12,1],[5,18,4]]){
 const lam=X.map((x,i)=>4*(10000*(x-j)**2+(x*i)**2)-10000*r4);
 a.equal(lam.reduce((m,l,i)=>m+(l<0?1<<i:0),0),mask);
 // z=(100x+i x sqrt(-1))/100; scale h=z^2 by 10000.
 const Q=X.map((x,i)=>[lam[i]*x*x*(10000-i*i),lam[i]*x*x*200*i].map(BigInt));
 const T=datum(n,mask);
 for(const u of range(n))for(const v of range(n))
  a.equal(Q[u][0]*Q[v][1]-Q[u][1]*Q[v][0]>0n,T.A[u][v]);
 const S=summary(T);a.equal(S.H,mask===12?393:347);
 a.equal(S.blueTilings,mask===12?23:0);
 console.log(JSON.stringify({mask,...S}));
}
console.log("Round score, rooted profiles, and blue canonical cuts checked on "+count+" strong sign words through order 8.");
