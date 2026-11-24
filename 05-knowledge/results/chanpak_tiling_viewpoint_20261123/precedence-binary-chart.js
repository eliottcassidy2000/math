'use strict';
const assert=require('node:assert/strict');
const idx=w=>31-Math.clz32(w&-w);
function profile(A){
 const n=A.length,N=1<<n,h=new Float64Array(N),D=new Float64Array(N*n);h[0]=1;
 for(let s=1;s<N;s++)for(let w=s;w;w&=w-1){
  let i=idx(w),t=s^(1<<i),v=+!t;
  for(let z=t&A[i];z;z&=z-1)v+=D[t*n+idx(z)];
  D[s*n+i]=v;h[s]+=v;
 }
 let Q=0;for(let s=0;s<N;s++)Q+=h[s]*h[N-1-s];return[h[N-1],Q];
}
function scan(A){
 const n=A.length,N=1<<n,W=new Int32Array(n),p=[];
 let H=0;
 function go(s,v){
  if(s===N-1){
   H++;for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)W[p[i]]^=1<<p[j];return;
  }
  for(let z=A[v]&~s;z;z&=z-1){const i=idx(z);p.push(i);go(s|1<<i,i);p.pop();}
 }
 for(let i=0;i<n;i++){p.push(i);go(1<<i,i);p.pop();}
 const scores=Array.from(W,w=>w.toString(2).replace(/0/g,'').length);
 return [H,new Set(scores).size===n,scores.every((s,i)=>s===n-1-i)];
}
function wrap(M,L,B){
 const n=L+M.length+1,A=Array(n).fill(0),a=n-1;
 for(let i=0;i<L;i++){
  for(let j=i+1;j<a;j++)A[i]|=1<<j;
  if(B>>(L-1-i)&1)A[a]|=1<<i;else A[i]|=1<<a;
 }
 for(let i=0;i<M.length;i++)A[L+i]=(M[i]<<L)|(1<<a);
 return A;
}
function defect(L,B){let E=0;for(let i=0;i<L;i++)if(B>>i&1)for(let j=0;j<i;j++)if(B>>j&1)
 E+=i-j===1?2:2**(i-j-1);return E;}
const J=[6,8,10,1], seed=[246,248,250,241,96,128,160,16];
assert.deepEqual(profile(J),[5,32]);
assert.deepEqual(profile(seed),[25,1024]);
assert.deepEqual(scan(seed),[25,true,true]);
for(let L=1;L<=7;L++)for(let B=0;B<1<<L;B++){
 const [H,p,ord]=scan(wrap([6,4,0],L,B));
 assert.equal(H,1+8*B-defect(L,B));
 const gap=(k)=>!(B&(B>>k));const ok=gap(1)&&gap(2);
 assert.equal(p,ok);assert.equal(ord,ok&&gap(3));
}
assert.deepEqual(scan(wrap(J,8,136)),[4317,true,true]);
assert.deepEqual(scan(wrap(seed,1,1)),[1049,true,true]);
assert.deepEqual(scan(wrap(seed,8,136)),[139089,true,true]);
assert.deepEqual(scan(wrap(seed,4,9)),[9141,true,false]);
for(let h of [1049,1048601,4194329])for(let p=2;p*p<=h;p++)assert(h%p);
console.log('checked exact binary chart, gap-three criterion, and 25 seed');
