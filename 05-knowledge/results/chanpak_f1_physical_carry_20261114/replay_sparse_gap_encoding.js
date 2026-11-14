// node replay_sparse_gap_encoding.js
"use strict";
const bit=(n,i)=>(n>>>i)&1;
function check(x,m){if(!x)throw Error(m)}
function each(n,f){
 for(let g=0;g<2**(n*(n-1)/2);g++){
  let i=0, r=Array(n).fill(0);
  for(let a=0;a<n;a++)for(let b=a+1;b<n;b++,i++){
   if(bit(g,i))r[a]|=1<<b;else r[b]|=1<<a;
  }
  f(r);
 }
}
function triangle(r,a,b,c){
 return bit(r[a],b)===bit(r[b],c)&&bit(r[b],c)===bit(r[c],a);
}
function order(r,L){
 if(!L.length)return [];
 const [v,...R]=L;
 return [
  ...order(r,R.filter(i=>bit(r[i],v))),
  v,
  ...order(r,R.filter(i=>bit(r[v],i)))
 ];
}
function backward(r,O){
 let t=0;
 for(let i=0;i<O.length;i++)for(let j=i+1;j<O.length;j++)t+=bit(r[O[j]],O[i]);
 return t;
}
function endings(r){
 const n=r.length, F=(1<<n)-1;
 const E=Array.from({length:F+1},()=>Array(n).fill(0));
 for(let i=0;i<n;i++)E[1<<i][i]=1;
 for(let s=1;s<=F;s++)for(let i=0;i<n;i++)if(bit(s,i))
  for(let j=0;j<n;j++)if(bit(r[j],i))E[s][i]+=E[s^(1<<i)][j];
 return E;
}
function H(r){
 if(!r.length)return 1;
 return endings(r).at(-1).reduce((a,b)=>a+b,0);
}
function induced(r,L){return L.map(i=>L.reduce((a,j,k)=>a|(bit(r[i],j)<<k),0))}
let tcount=0, contacts=0;
for(let n=0;n<=6;n++)each(n,r=>{
 const L=Array.from({length:n},(_,i)=>i);
 let t=0;
 for(let a=0;a<n;a++)for(let b=a+1;b<n;b++)for(let c=b+1;c<n;c++)t+=triangle(r,a,b,c);
 check(backward(r,order(r,L))<=t,"pivot encoding");tcount++;
 if(n<2||n>5)return;
 const F=(1<<n)-1,end=endings(r),start=endings(r.map((v,i)=>F^(1<<i)^v));
 const ep=Array.from({length:n},()=>Array(n).fill(0));
 const bwd=Array.from({length:n},()=>Array(n).fill(0));
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){
  for(let m=1;m<F;m++)ep[i][j]+=end[m][i]*start[F^m][j]+end[m][j]*start[F^m][i];
  const l=L.filter(k=>k!==i&&k!==j), W=induced(r,l), v=W.length;
  function marked(or){
   const flags=l.map(k=>or?bit(r[k],i)|bit(r[k],j):bit(r[k],i)&bit(r[k],j));
   return W.map((row,k)=>row|(flags[k]<<v)).concat(flags.reduce((s,b,k)=>s|((1-b)<<k),0));
  }
  check(2*Math.min(H(marked(true)),H(marked(false)))<=ep[i][j],"pair average");
  bwd[i][j]=backward(W,order(W,Array.from({length:v},(_,k)=>k)));
 }
 for(let u=1;u<=F;u++){
  const A=L.filter(i=>bit(u,i));if(A.length<2)continue;
  let E=0,a=-1,b=-1;
  for(const i of A)for(const j of A)if(i<j){
   E+=ep[i][j];if(a<0||ep[i][j]<ep[a][b]){a=i;b=j;}
  }
  const p=A.length*(A.length-1)/2;
  check(2*p*(2*bwd[a][b]+1)<=E,"refined bound after deleting contact pair");
  contacts++;
 }
});
check(tcount===33868&&contacts>20000,"full controls");
console.log("PASS:",tcount,"labelled pivot orders and",contacts,"multi-contact marks");
