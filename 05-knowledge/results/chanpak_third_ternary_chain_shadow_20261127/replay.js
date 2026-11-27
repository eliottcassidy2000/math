"use strict";
const assert=require("node:assert/strict");
const R=n=>Array.from({length:n},(_,i)=>i);
function word(w){
 let x=0n,c=0n;
 for(const b of w){if(b==="0")c+=(x+1n)/2n;x=2n*x+BigInt(b);}
 return 2n*c+1n;
}
function paths(adj){
 const n=adj.length,S=1<<n,P=new Float64Array(S*n);
 for(let j=0;j<n;j++)P[(1<<j)*n+j]=1;
 for(let m=1;m<S;m++)for(let j=0;j<n;j++)if(m>>j&1){
  const val=P[m*n+j]; if(!val)continue;
  let next=adj[j]&(~m);
  while(next){const bit=next&-next,k=31-Math.clz32(bit);next^=bit;P[(m|bit)*n+k]+=val;}
 }
 return R(n).reduce((s,j)=>s+P[(S-1)*n+j],0);
}
function apex(w){
 const n=w.length;
 return R(n+1).map(i=>R(n+1).reduce((m,j)=>{
  if(i===j)return m;
  const a=i===n?w[j]==="1":j===n?w[i]==="0":i<j;
  return m+(a?(1<<j):0);
 },0));
}
for(let b=1;b<=8;b++){
 const w="110"+"1".repeat(b)+"0", H=5n+7n*2n**BigInt(b);
 assert.equal(word(w),H);
 assert.equal(paths(apex(w)),Number(H));
}
const core=[238,220,216,16,225,206,8,72];
assert.equal(paths(core),63);
const pa=core.map(r=>r+(7<<8)).concat([1<<9,1<<10,1<<8]);
assert.equal(paths(pa),189);
function mexp(b,m){let x=2n,r=1n;while(b){if(b&1n)r=r*x%m;x=x*x%m;b>>=1n;}return r;}
const replay=[];
for(let k=1,b=1n,q=3n,period=2n;k<=12;k++,q*=3n,period*=3n){
 const prev=k===1?2n:period/3n;
 const choose=[b,b+prev,b+2n*prev].find(u=>(5n+7n*mexp(u,q))%q===7n%q);
 assert(choose!==undefined);b=choose;
 assert.equal((5n+7n*mexp(b+3n,q))%q,21n%q);
 assert.equal((5n+7n*mexp(b+period,q))%q,7n%q);
 replay.push(String(b));
}
console.log("3-adic shadow exponents at k=1,...,12:",replay.join(", "));
console.log("Hamiltonian path DP independently checked the first eight sparse shadow words, the 63 core, and its ordinal 189 child.");
for(let h=1;h<=4000;h+=2){
 let a=h,k=0;while(a%3===0){a/=3;++k;}
 assert.equal(a*3**k,h);assert(a%3!==0);
 const t=((h-1)/2).toString(3);
 assert(t.endsWith("1".repeat(k)));
 if(k<t.length)assert(t[t.length-1-k]!=="1");
}
console.log("All odd targets below 4000 satisfy the ternary trailing-ones and root rules.");
