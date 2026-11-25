'use strict';
const assert=require('node:assert/strict');
const gcd=(a,b)=>b?gcd(b,a%b):(a<0n?-a:a);
const lcm=(a,b)=>a/gcd(a,b)*b;
const A=[1000,1000,1001,1008,1044,1146,1303,1517,1760,2065,2412,2820,3299,3858,4518,5288,6189,7233,8477,9915,11584,13508,15907,18467,20748].map(BigInt);
const Z=[788,862,953,1053,1152,1234,1287,1291,1309,1300,1302,1304,1303,1306,1305,1304,1301,1301,1296,1289,1280,1271,1213,1134,1032].map(BigInt);
const B=1000n, m=130n, D=B*A.reduce(lcm), N=3251n;
function ldlsigns(X){
 const W=X.map(r=>r.slice());let p=1n;
 const signs=[];
 for(let k=0;k<X.length;k++){
  const v=W[k][k];assert.notEqual(v,0n);
  signs.push(v*p>0n);
  for(let i=k+1;i<X.length;i++)for(let j=k+1;j<X.length;j++){
   const x=v*W[i][j]-W[i][k]*W[k][j];assert.equal(x%p,0n);W[i][j]=x/p;
  }
  p=v;
 }
 return signs;
}
const G=A.map((a,i)=>A.map((b,j)=>i===j?0n:i>j?
 Z[i]*a*(D/(B*b))-Z[j]*b*(D/(B*a)):
 Z[j]*b*(D/(B*a))-Z[i]*a*(D/(B*b))));
const X=G.slice(0,24).map((r,i)=>r.slice(0,24).map((g,j)=>g-G[i][24]-G[24][j]));
assert.equal(ldlsigns(X).filter(Boolean).length,1);
const shift=(q,d)=>X.map((r,i)=>r.map((x,j)=>q*D*(i===j?2n:1n)-d*x));
let cap=0n;
for(let q=40n;q<=50n;q++){
 if(ldlsigns(shift(q,10n)).every(Boolean)){cap=q;break;}
}
assert.equal(cap,46n);
const S=G.flat().reduce((p,v)=>p+v);
let f=1n, sum=0n;
for(let i=1n;i<=13n;i++)f*=i;
for(let i=0n,g=1n;i<=13n;i++){if(i>0n)g*=i;sum+=f/g}
const E=271828183n, dE=100000000n, Elow=27182818n, dElow=10000000n;
assert(sum*dElow>f*Elow);
assert((13n*sum+1n)*dE<13n*f*E); // tail <= 1/(13*13!)
let fact=1n;
for(let k=1n;k<=N;k++){
 fact*=k;
 assert(fact*E**k<((k+5n)*dE)**k);
 if(k<=59n)assert(fact*(10n*E)**k<((10n*k+31n)*dE)**k);
}
// Scale every weighted arc/correction by 100*D, retaining all decimals.
const L=100n*D;
let base=L*20n;
for(let i=0;i<25;i++){
 const a=A[i],z=Z[i];
 base+=m*(4n*B*(L/a) + z*(m-1n)*(L/(2n*B)) + z*a*(L/(4n*B*B)) 
    +31n*z*a*(L/(10n*B*B))+19n*z*(L/(10n*B)));
 for(let j=i+1;j<25;j++)base+=m*m*z*a*(L/(B*A[j]));
}
// Cost <= (m/2)*(65*S/(25D)+15*cap/10).
assert.equal(L%(25n*D),0n);
base+=(m/2n)*(65n*S*(L/(25n*D))+15n*cap*(L/10n));
let spine=0n;
for(let i=0;i<24;i++)spine+=100n*G[i][i+1];
const prod=Z.reduce((p,v)=>p*v,1n)**m;
const den=B**3250n;
const qs=[];
for(let add of [0n,spine]){
 const a=(base+add)*dElow*1000n,b=L*N*Elow;
 const q=(a+b-1n)/b;
 qs.push(q);
 for(let bits=26894n;bits<=26894n;bits++){
  if(5n*q**N*den < (4n*B**N*prod<<bits)){
   console.log(String(cap),String(q),String(bits));break;
  }
 }
 assert(5n*q**N*den < (4n*B**N*prod<<26894n));
}
