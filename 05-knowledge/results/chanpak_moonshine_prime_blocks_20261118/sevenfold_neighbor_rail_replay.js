const assert = require("node:assert/strict");
const T = [
 [80,121,227,5,12,25,56,123],
 [96,121,203,33,13,20,56,123],
 [240,249,203,161,140,20,184,32]
];
function expand(A,m,k){
 const L=[...Array(m).fill(0),...Array(k).fill(1),2,3,4,5,6,7];
 return L.map((a,i)=>L.reduce((mask,b,j)=>mask|((a===b?+(i<j):(A[a]>>b)&1)<<j),0));
}
function H(A){
 const n=A.length,L=1<<n;
 const D=Array.from({length:L},()=>new Float64Array(n));
 for(let i=0;i<n;i++)D[1<<i][i]=1;
 for(let s=1;s<L;s++)for(let i=0;i<n;i++){
  if(!D[s][i])continue;
  for(let w=A[i]&~s;w;w&=w-1){
   const bit=w&-w,j=31-Math.clz32(bit);
   D[s|bit][j]+=D[s][i];
  }
 }
 return D[L-1].reduce((a,b)=>a+b,0);
}
function pair(m,k){
 m=BigInt(m);k=BigInt(k);
 return (16n*2n**m+24n*3n**m-16n)*2n**k+3n*(2n**m-1n);
}
function prime(n){
 if(n<=1n)return false;
 for(let d=2n;d*d<=n;d++)if(n%d===0n)return false;
 return true;
}
function f(x){
 let z=1n;
 while(x){
  let y=x/2n;
  if(x%2n===0n)z+=y+y%2n;
  x=y;
 }
 return z;
}
for(let m=1;m<=6;m++)for(let k=1;k<=3;k++){
 assert.equal(pair(m,k)%7n===0n,[1,2].includes(m%6)&&k%3===0);
}
for(let m=1;m<=19;m+=6){
 let exc=0;
 for(let j=0;j<7;j++)exc+=+(pair(m,3+3*j)%49n===0n);
 assert.equal(exc,1);
}
assert.deepEqual(T.map(A=>H(expand(A,1,3))),[705,707,709]);
assert.equal(H([8,13,49,84,99,11,39]),101);
assert.ok(prime(101n)&&prime(709n)&&prime(411941n));
assert.equal(pair(1,15),2883587n);
const b=3n,e=84871n;
assert.equal(b*e+f(b-1n)*f(e-1n)-(b%2n)*(e%2n),411941n);
console.log("PASS: periodic eligibility, distinct v7 filter, literal 705/707/709 and 101, and exact F1 411941");
