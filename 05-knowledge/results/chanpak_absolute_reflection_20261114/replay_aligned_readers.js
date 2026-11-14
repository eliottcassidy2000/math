// node replay_aligned_readers.js
"use strict";
function check(q,m){if(!q)throw Error(m);}
function Hsub(rows){
 const n=rows.length,F=(1<<n)-1,pre=rows.map((_,i)=>rows.reduce((s,r,j)=>s|(((r>>i)&1)<<j),0));
 const H=new Float64Array(F+1),A=new Float64Array(n*(F+1));H[0]=1;
 for(let i=0;i<n;i++)A[(1<<i)*n+i]=1;
 for(let s=1;s<=F;s++){
  let z=s;
  while(z){
   const b=z&-z,i=31-Math.clz32(b),m=s^b;
   let p=pre[i]&m,v=A[s*n+i];
   while(p){const d=p&-p;p^=d;v+=A[m*n+31-Math.clz32(d)];}
   A[s*n+i]=v;H[s]+=v;z^=b;
  }
 }
 return H;
}
function TT(n){return Array.from({length:n},(_,i)=>((1<<n)-1)^((2<<i)-1));}
function join(A,B){return A.map(r=>r|(((1<<B.length)-1)<<A.length)).concat(B.map(r=>r<<A.length));}
function tri(A,B){
 const a=A.length,n=a+B.length;
 return join(A,B).map((r,i)=>r|((i>=a?1:0)<<n)).concat((1<<a)-1);
}
const cores=[
 [91,2176,[252,41,10,16,38,12,62,126]],
 [59,940,[12,17,2,102,13,23,55]],
 [413,6392,[204,81,34,38,13,19,188,62]],
 [409,6352,[204,49,50,70,73,25,38,126]],
 [251,4520,[252,105,234,144,70,152,168,18]],
 [149,3176,[188,73,10,16,38,14,189,62]],
 [53,904,[92,41,10,16,6,93,30]]
];
for(const [h,f,R] of cores){
 const H=Hsub(R),F=H.length-1;
 check(H[F]===h&&H.reduce((s,v,i)=>s+v*H[F^i],0)===f,"core profile");
}
const r=i=>cores[i][2];
const A=tri(join(r(0),r(1)),TT(2)), B=tri(join(r(2),TT(6)),TT(4));
check(A.length===18&&B.length===19&&Hsub(A).at(-1)===6130951&&Hsub(B).at(-1)===6130951,"literal neutral pair");
function prime(n){if(n<2)return false;for(let i=2;i*i<=n;i++)if(n%i===0)return false;return true;}
check(prime(6130951),"neutral base prime");
const a=3*6352*4520,b=15*3176*904*2,c=409*251,d=13*149*53;
check(a===b&&d-c===2&&prime(a*64-c)&&prime(a*64-d),"aligned twin-prime repair");
console.log("PASS: seven literal source profiles, two full neutral tournaments, and the prime checks");
