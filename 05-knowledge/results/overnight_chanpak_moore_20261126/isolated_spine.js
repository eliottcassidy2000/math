"use strict";
const assert=require("node:assert/strict");

function poly(n,i,j) {
 const N=1<<n,U=N-1, mask=(1<<i)|(1<<j), inner=((1<<j)-1)^((1<<(i+1))-1);
 const low=x=>x&-x, high=x=>x?1<<(31-Math.clz32(x)):0;
 const F=new Int32Array(N),E=new Int32Array(N),R=new Int32Array(N);
 for(let T=U^mask;;T=(T-1)&(U^mask)){
  const S=U^T,I=S&inner,plo=S&((1<<i)-1),phi=S&~((1<<(j+1))-1);
  const u=low(T),v=high(T);
  function charge(a,b,z){
   F[b]+=z;if(u)F[b|u]-=z;
   F[v]+=z;F[v|a]-=z;
   R[b]+=z;R[v]+=z;
   if(T===0){E[0]+=z;E[a]-=z;}
  }
  if(!I)charge(low(S),high(S),-1);
  for(let A=I;;A=(A-1)&I){
   charge(low(plo|A|(1<<j)),high(phi|(I^A)|(1<<i)),1);
   if(!A)break;
  }
  if(!T)break;
 }
 for(let s=1;s<N;s*=2)for(let j=0;j<N;j+=2*s)for(let k=j;k<j+s;k++){
  F[k+s]+=F[k];E[k+s]+=E[k];R[k+s]+=R[k];
 }
 return [F,E,R];
}
function scan(n){
 const results=[], prof=new Map();
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){
  const [F,E,R]=poly(n,i,j);
  for(let v=0;v<1<<n;v++)if(F[v]===0&&R[v]===0&&E[v]!==0){
   const d=[v+1,i,j,E[v]];results.push(d);
   prof.set(E[v],(prof.get(E[v])||0)+1);
  }
 }
 return {n, hits:results.length, sizes:[...prof], unit:results.filter(z=>Math.abs(z[3])===2),preview:results.slice(0,60)};
}
const n=Number(process.argv[2]||12);
assert(n>=4&&n<=16);
const got=scan(n);
const expected=[];
for(let i=1;2*i+1<n;i++)expected.push([2**(2*i+1)+1,i,2*i+1,2**i]);
assert.deepEqual(got.preview,expected);assert.equal(got.hits,expected.length);
console.log(JSON.stringify(got,null,2));
