// node replay_gap_two_63.js
"use strict";
const check=(c,m)=>{if(!c)throw Error(m)};
const bits=(n,i)=>(n>>>i)&1;
const pop=n=>{let s=0;while(n){n&=n-1;s++}return s};
function paths(rows){
 const n=rows.length,F=(1<<n)-1;
 const E=Array.from({length:F+1},()=>Array(n).fill(0)),H=Array(F+1).fill(0);
 H[0]=1;for(let i=0;i<n;i++)E[1<<i][i]=1;
 for(let s=1;s<=F;s++){
  for(let j=0;j<n;j++)if(bits(s,j)){
   for(let i=0;i<n;i++)if(bits(rows[i],j))E[s][j]+=E[s^(1<<j)][i];
  }
  H[s]=E[s].reduce((a,b)=>a+b,0);
 }
 return {H,E,F};
}
function decorate(rows,U){
 const n=rows.length;
 return rows.map((r,i)=>r|((1-bits(U,i))<<n)).concat(U);
}
function strong(rows){
 for(let root=0;root<rows.length;root++){
  let s=1<<root,t=0;
  while(t!==s){t=s;rows.forEach((r,i)=>{if(bits(t,i))s|=r})}
  if(s!==((1<<rows.length)-1))return false;
 }
 return true;
}
function canonical(rows){
 const N=rows.length,deg=rows.map(pop),d=deg.slice().sort((a,b)=>a-b);
 let best=Infinity;
 function visit(p,M){
  if(p.length<N){
   for(let j=0;j<N;j++)if(!bits(M,j)&&deg[j]===d[p.length])visit(p.concat(j),M|(1<<j));
   return;
  }
  let k=0;
  for(let i=0;i<N;i++)for(let j=i+1;j<N;j++)k=2*k+bits(rows[p[i]],p[j]);
  best=Math.min(best,k);
 }
 visit([],0);return best;
}
function profile(rows,U){
 const {H,E,F}=paths(rows),A=paths(rows.map((r,i)=>F^(1<<i)^r)).E;
 let P=0,Q=0,f=0,T=0;
 for(let s=0;s<=F;s++){
  f+=H[s]*H[F^s];
  for(let i=0;i<rows.length;i++){
   if(bits(U,i))P+=A[s][i]*H[F^s];
   else Q+=E[s][i]*H[F^s];
  }
 }
 for(let i=0;i<rows.length;i++)if(!bits(U,i))T+=A[F][i];
 return {h:H[F],f,P,Q,T,K:paths(decorate(rows,U)).H.at(-1),
  miss:E[F].reduce((acc,b,i)=>acc+(bits(U,i)?0:b),0)};
}
function g(x){
 let h=1;
 while(x){const p=x%2;x=Math.floor(x/2);if(!p)h+=x+x%2}
 return h;
}
function wing(x){
 return x===0?[]:x.toString(2).split("").map(Number);
}
function literal(B,U,x){
 const w=wing(x),a=w.length,b=B.length,n=a+b;
 const D=w.map((_,i)=>((1<<n)-1)^((2<<i)-1)).concat(B.map(r=>r<<a));
 const V=w.reduce((m,t,i)=>m|(t<<i),U<<a);
 return decorate(decorate(D,V),1<<n);
}
const expected=[1,1,2,4,12,56,456];
let layer=[[0]],tested=0,frames=0;
const result=new Set();
for(let b=1;b<=7;b++){
 check(layer.length===expected[b-1],"number of iso classes at "+b);
 for(const B of layer.filter(strong)){
  for(let U=0;U<(1<<b);U++){
   if(b-pop(U)>2||!U)continue;
   const {h,f,P,Q,T,K,miss}=profile(B,U);
   if(miss!==2)continue;
   frames++;
   for(let x=0;(f+Q)*x+K+P+h-2<=63;x++){
    const H=h*g(x)+(f+Q)*x+K+P-2+T*(x%2);
    // Avoid explosive full DP outside the proof-led finite target window.
    if(H>63)continue;
    const rows=literal(B,U,x);
    check(paths(rows).H.at(-1)===H&&strong(rows),"literal inverse mismatch");
    result.add(H);tested++;
   }
  }
 }
 if(b<7){
  const next=new Map();
  for(const B of layer)for(let U=0;U<(1<<b);U++){
   const rows=decorate(B,U),k=canonical(rows);next.set(k,rows);
  }
  layer=[...next.values()];
 }
}
check([...result].sort((a,b)=>a-b).join(",")==="9,23,29,37,51,53,59","one-contact list");
const rows=literal([0],0,19);
check(rows.join(",")==="190,252,248,176,160,192,25,64"&&strong(rows),"eight-vertex witness");
const old=rows.slice(0,-1).map(r=>r&127);
const {H,E,F}=paths(old), A=paths(old.map((r,i)=>F^(1<<i)^r)).E;
let P=0;for(let s=1;s<=F;s++)P+=A[s][6]*H[F^s];
check([H[F],P,E[F][6],paths(rows).H.at(-1)].join(",")==="25,39,1,63","gap-one profile");
console.log("PASS: complete marked terminal search through order seven,",frames,"eligible frames,",tested,"successful literal realizations; (63,64) verified");
