// node replay_one_contact_upper.js
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

const key=(H,R)=>H+","+R;
let layer=[[0]], out=new Map();
for(let a=1;a<=7;a++){
 for (const A of layer.filter(strong)){
  const {H,E,F}=paths(A);
  for(let u=0;u<a;u++){
   const S=E[F][u], v=profile(A,1<<u), H1=v.K, R1=v.h+v.P, e=2*H1-R1;
   check(R1-H1===S,"unique-gap identity");
   if(a>1){const t=pop(A[u]);check(e>=2*(a-1)+2**t-1-t,"upper-deficit bound");}
   if(e>0&&e<=13)out.set(key(H1,R1),A.concat(u));
  }
 }
 if(a<7){
 const next=new Map();
 for(const A of layer)for(let U=0;U<(1<<a);U++){const G=decorate(A,U);next.set(canonical(G),G);}
 layer=[...next.values()];
 }
}
const primitive=[...out].map(([p,rows])=>[...p.split(",").map(Number),rows]);
const full=primitive.map(([H,R])=>[H,R,2*H-R]);
for(const [H,R] of primitive)if(3*(2*H-R)<=13)full.push([3*H,3*R,3*(2*H-R)]);
full.sort((a,b)=>a[2]-b[2]||a[0]-b[0]);const actual=[...new Set(full.map(row=>row.join(",")))];
check(actual.join(";")==="5,6,4;9,11,7;9,10,8;11,12,10;15,19,11;13,14,12;15,18,12","upper corridor");
console.log("PASS: exact strip 1 <= 2H-R <= 13; seven pairs");
