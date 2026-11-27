
"use strict";
const assert=q=>{if(!q)throw Error("check failed")};
const car=n=>Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?0:(j-i+n)%n<=n>>1?1:-1));
const mul=(A,B)=>A.map(r=>r.map((_,j)=>r.reduce((s,a,k)=>s+a*B[k][j],0)));
const eq=(a,b)=>JSON.stringify(a)===JSON.stringify(b);
const sw=(A,e)=>A.map((r,i)=>r.map((x,j)=>e[i]*e[j]*x));
function subsets(n,k){
 const out=[];
 function go(a,j){if(a.length===k){out.push(a);return}for(let i=j;i<=n-k+a.length;i++)go(a.concat(i),i+1)}
 go([],0);return out
}
function trans(C,[a,b,c]){
 return +(C[a][b]!==C[b][c]||C[b][c]!==C[c][a])
}
function calc(C,G,[a,b,c]){
 const H=C[a].reduce((s,x,j)=>s+x*C[b][j]*C[c][j],0);
 const J=G[a][b]*(C[c][a]+C[c][b])+G[a][c]*(C[b][a]+C[b][c])+G[b][c]*(C[a][b]+C[a][c]);
 return (3*H-J)/2
}
const key=t=>t.join(",");
const parity=x=>x&1;
function test(C,checkCarousel){
 const n=C.length,m=n>>1,G=mul(C,C),triples=subsets(n,3),four=subsets(n,4);
 assert(C.every(r=>r.reduce((a,b)=>a+b,0)===0));
 const flags=new Map(), degrees=new Map(),hist=new Map();
 for(const T of triples){
  const A=calc(C,G,T);assert(A%2===0);
  flags.set(key(T),parity(A/2));degrees.set(key(T),0);
  hist.set(Math.abs(A),(hist.get(Math.abs(A))||0)+1);
 }
 for(const S of four){
  const Ts=S.map((_,i)=>S.filter((_,j)=>i!==j));
  const d=parity(Ts.reduce((v,T)=>v+trans(C,T),0));
  assert(d===parity(Ts.reduce((v,T)=>v+flags.get(key(T)),0)));
  if(d)for(const T of Ts)degrees.set(key(T),degrees.get(key(T))+1);
 }
 for(const [T,flag] of flags)assert(flag===parity(degrees.get(T)));
 if(checkCarousel){
  let sum=0;
  for(let k=1;k<=m-2;k++){
   const f=2*n*Math.floor((m-k)/2);
   assert(hist.get(4*k)===f);
   sum+=f;
  }
  assert(sum===triples.length-hist.get(0));
 }
}
function bank(Q,blocks){
 const I=blocks.flatMap((B,k)=>B.map((_,i)=>[k,i]));
 return I.map(([p,i])=>I.map(([q,j])=>p===q?blocks[p][i][j]:Q[p][q]))
}
function flipBlock(A,S){
 const set=new Set(S);
 return A.map((r,i)=>r.map((a,j)=>set.has(i)&&set.has(j)?-a:a))
}
function paths(A,k){
 const n=A.length;
 const out=A.map(r=>r.reduce((s,a,i)=>a>0?s|1<<i:s,0));
 function walk(i,used,left){
  let dest=out[i]&~used;
  if(!left)return 1;
  let result=0;
  while(dest){const b=dest&-dest;dest^=b;result+=walk(31-Math.clz32(b),used|b,left-1)}
  return result
 }
 return out.reduce((s,_,i)=>s+walk(i,1<<i,k-1),0)
}
function coeff(A,S,k){
 let ans=0;
 for(let m=0;m<1<<S.length;m++){
  const e=Array(A.length).fill(1);let sign=1;
  S.forEach((i,j)=>{if(m>>j&1){e[i]=-1;sign=-sign}});
  ans+=sign*paths(sw(A,e),k)
 }
 return ans/2**S.length
}
for(let p=3;p<=35;p+=2){
 const C=car(p);test(C,true);
 const triples=subsets(p,3);
 let st=173+p;
 for(let step=0;step<12*p;step++){
  st=(Math.imul(st,1664525)+1013904223)>>>0;
  const T=triples[st%triples.length];
  if(!trans(C,T))for(let a of T)for(let b of T)C[a][b]*=-1;
 }
 test(C,false)
}
console.log("Both diamond incidence identities and carousel census: every odd size 3..35.");
const Q2=[[0,1],[-1,0]];
for(const kind of ["carousel","Paley"]){
 const C=kind==="carousel"?car(7):Array.from({length:7},(_,i)=>Array.from({length:7},(_,j)=>i===j?0:[1,2,4].includes((j-i+7)%7)?1:-1));
 const A=bank(Q2,[C,[[0]]]),B=flipBlock(A,[0,1,2,3,4,5,6]),G=mul(C,C);
 assert(eq(mul(A,A),mul(B,B)));
 const f=[];
 for(let mask=0;mask<128;mask++){
  const e=Array.from({length:8},(_,i)=>mask>>i&1?-1:1);
  f.push(paths(sw(B,e),6)-paths(sw(A,e),6))
 }
 let detectors=0,flags=0;
 for(const T of subsets(7,3)){
  const A=calc(C,G,T);
  const four=f.reduce((s,n,mask)=>s+(T.reduce((z,i)=>mask>>i&1?-z:z,1))*n,0)/128;
  assert(four===-A);detectors+=+(A!==0);flags+=parity(A/2);
 }
 console.log(kind,"order-seven host: actual path detector and flag counts",detectors,flags);
}
const C=car(3);C[0][1]*=-1;C[1][0]*=-1;
const A=bank(Q2,[C,[[0]]]),B=flipBlock(A,[0,1,2]);
assert(coeff(B,[0,3],3)-coeff(A,[0,3],3)===2);
const Q4=[[0,1,1,1],[-1,0,1,-1],[-1,-1,0,1],[-1,1,-1,0]];
const X=bank(Q4,[car(7),car(3),car(3),car(3)]);
const Y=X.map((r,i)=>r.map((x,j)=>i>=7&&j>=7&&((i-7)/3|0)!==((j-7)/3|0)?-x:x));
assert(eq(mul(X,X),mul(Y,Y)));
const T=[7,8,10,0];
const delta=coeff(Y,T,6)-coeff(X,T,6);
assert(delta!==0);
console.log("Equality of squares under reversing the 3-small-block quotient, followed by an exact six-path coefficient change:",delta);
