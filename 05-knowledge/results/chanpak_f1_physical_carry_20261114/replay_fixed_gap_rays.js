// node replay_fixed_gap_rays.js
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

function graphs(n, use) {
 const N=n*(n-1)/2;
 for(let k=0;k<2**N;k++){
  const r=Array(n).fill(0);let b=0;
  for(let i=0;i<n;i++)for(let j=i+1;j<n;j++,b++){
   if(bits(k,b))r[i]|=1<<j;else r[j]|=1<<i;
  }
  use(r);
 }
}
function witness(rows){
 const {H,E,F}=paths(rows),A=paths(rows.map((r,i)=>F^(1<<i)^r)).E;
 const n=rows.length;
 const p=Array(n).fill(0), e=Array.from({length:n},()=>Array(n).fill(0));
 let f=0;
 for(let s=0;s<=F;s++){
  f+=H[s]*H[F^s];
  for(let i=0;i<n;i++)if(bits(s,i)){
   p[i]+=A[s][i]*H[F^s];
   for(let j=0;j<n;j++)if(!bits(s,j))e[i][j]+=A[s][i]*E[F^s][j];
  }
 }
 return U=>{
  let P=0,S=0,L=0;
  for(let i=0;i<n;i++)if(bits(U,i)){
   P+=p[i];S+=E[F][i];
   for(let j=0;j<n;j++)if(bits(U,j))L+=e[i][j];
  }
  return {h:H[F],f,P,S,E:L,K:H[F]+P-S-L};
 };
}
function induced(rows,Y){return Y.map(i=>Y.reduce((r,j,k)=>r|(bits(rows[i],j)<<k),0));}
function comps(rows,Y){
 let remaining=Y.reduce((a,i)=>a|(1<<i),0);
 const reach={};
 for(const i of Y){
  let s=1<<i,t=0;
  while(s!==t){t=s;for(const j of Y)if(bits(t,j))s|=rows[j]&remaining}
  reach[i]=s;
 }
 const out=[];
 while(remaining){
  const i=Y.find(j=>bits(remaining,j)&&(reach[j]&remaining)===remaining);
  const C=Y.filter(j=>bits(remaining,j)&&bits(reach[j],i));
  for(const j of C)remaining^=1<<j;
  out.push(C);
 }
 return out;
}
function peel(rows,U){
 const Z=rows.map((_,i)=>i), W=Z.filter(i=>!bits(U,i)), C=comps(rows,W);
 const G=[],D=[];
 while(C.length&&C[0].every(i=>(rows[i]&U)===U))G.push(...C.shift());
 while(C.length&&C.at(-1).every(i=>(rows[i]&U)===0))D.push(...C.pop());
 const core=Z.filter(i=>!G.includes(i)&&!D.includes(i));
 return [G,core,D];
}
let pairs=0, multi=0;
const rays=new Map(), singles=new Map();
for(let n=2;n<=6;n++)graphs(n,rows=>{
 const opt=witness(rows);
 for(let U=3;U<(n<=5?1<<n:4);U++){
  const k=pop(U);if(k<2)continue;
  const [G,C,D]=peel(rows,U);
  const {h,f,P,S,E,K}=opt(U);
  if(k===2&&U===3)pairs++;
  if(k>2)multi++;
  check(C.length<=E+k,"fixed contact core bound");
  const uC=C.reduce((r,i,j)=>r|(bits(U,i)<<j),0);
  const m=witness(induced(rows,G))(0).h;
  const v=witness(induced(rows,C))(uC);
  const B=witness(induced(rows,D))(0);
  check(h+P===m*(B.h*v.h+B.f*v.P) &&
     E+S===m*(D.length?B.h*v.E:v.E+v.S),"wing factorization");
  if(n<=4)check(K===paths(decorate(rows,U)).H.at(-1),"literal K");
  if(k===2&&U===3&&!G.length&&!D.length){
   if(E+S===4)singles.set(K,rows);
   if(E===4)rays.set((v.h-4)+","+v.P,rows);
  }
 }
});
check(pairs===33866,"marked pair count");
const lists=[...singles.keys()].sort((a,b)=>a-b), rr=[...rays.keys()].sort();
check(multi===16712&&lists.join(",")==="5,15","gap-four terminal list");
check(rr.join(";")==="-1,15;-1,9;-3,9;1,11;1,12;1,16;5,22;5,25;5,30;5,31;5,34;7,30;7,33;9,33","gap-four rays");
console.log("PASS: 33,866 marked pairs, 16,712 multi-contact marks, and complete gap-four two-contact menu");
