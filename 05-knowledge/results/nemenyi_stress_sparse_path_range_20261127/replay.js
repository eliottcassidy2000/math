"use strict";
const assert=require("node:assert/strict"),limit=4097;
const R=n=>Array.from({length:n},(_,i)=>i);
const values=new Set([1]),samples=new Map();
for(let n=3;n<=13;n++){
 const N=1<<n,I=R(n),seen=new Set(),dp=new Uint32Array(N*n);let pass=0;
 for(let m=0;m<1<<(n-1);m++){
  const bit=i=>m>>i&1,c=I.filter(i=>!bit(i)).concat(I.filter(bit));
  const adj=I.map(i=>I.reduce((b,j)=>b|(i!==j&&((i<j)===(bit(i)===bit(j)))?1<<j:0),0));
  const d=c.map(i=>I.reduce((b,j)=>b+(adj[i]>>j&1),0));
  if(d.includes(0))continue;
  const key=I.map(j=>d.slice(j).concat(d.slice(0,j)).join(",")).sort()[0];
  if(seen.has(key))continue;
  seen.add(key);dp.fill(0);
  for(const i of I)dp[(1<<i)*n+i]=1;
  for(let u=1;u<N-1;u++)for(let b=u;b;b&=b-1){
   const i=31-Math.clz32(b&-b),v=dp[u*n+i];if(!v)continue;
   for(let c=adj[i]&~u;c;c&=c-1){
    const e=c&-c,j=31-Math.clz32(e),o=(u|e)*n+j;dp[o]=Math.min(limit+1,v+dp[o]);
   }
  }
  const h=I.reduce((b,i)=>b+dp[(N-1)*n+i],0);
  assert(h>=1+2**(n-2));
  if(h<=limit){values.add(h);pass++;if(!samples.has(h))samples.set(h,[n,m]);}
 }
 console.log(JSON.stringify({n,classes:seen.size,pass}));
}
assert.equal(values.size,75);
// At order 14 the floor is exactly 4097, attained by reversing the extreme edge.
// The floor excludes all orders >14 and all other heights at order 14.
assert(!values.has(4097));
values.add(4097);
const ordered=[...values].sort((a,b)=>a-b);
const products=new Uint8Array(limit+1);products[1]=1;
for(let i=2;i<=limit;i++)products[i]=+ordered.some(h=>h>1&&i%h===0&&products[i/h]);
const holes=R(6).map(i=>7*3**i);
for(const h of holes)assert(!values.has(h)&&!products[h]);
for(const h of [13,49,343])assert(!values.has(h));
console.log(JSON.stringify({single:ordered.length,joins:products.reduce((a,b)=>a+b,0),holes,witness23:samples.get(23),witness65:samples.get(65)}));
console.log(ordered.join(","));
for(let k=25;k<1000;k+=2){
 const m=Math.floor(k/12),t=[0,Math.floor(k/3),Math.floor(2*k/3)];
 for(let i=0;i<3;i++){
  const gap=(t[(i+1)%3]-t[i]+k)%k;
  assert(gap-m+1>=1&&gap+m-1<=(k-1)/2);
 }
}
