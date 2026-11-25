'use strict';
const A=require('node:assert/strict');
const w=[1000,1000,1000,1000,1000,1000,1000,1016,1078,1143,1213,1287,1366,1449,1537,1631,1730,1836,1948,2067,2193,2326,2468,2619,2779,2948,3128,3319,3521,3736,3964,4206,4463,4735,5024,5331,5656,6002,6368,6757,7169,7607,8072,8564,9087,9642,10231,10856,11519,12223,12970,13762,14603,15495,16442,17446,18513,19644,20845,22118,23470,24905,26427,28043,29757].map(BigInt),z=[2105,2198,2299,2411,2533,2669,2820,2942,2942,2942,2943,2943,2943,2944,2944,2944,2945,2945,2945,2946,2946,2946,2946,2947,2947,2947,2948,2948,2948,2949,2949,2949,2950,2950,2950,2951,2951,2951,2951,2952,2952,2952,2953,2953,2953,2954,2954,2954,2955,2955,2955,2956,2956,2956,2956,2957,2957,2957,2958,2958,2958,2959,2959,2959,2960].map(BigInt);
const B=1000n,M=w.length,m=50n,N=3251n,p=7941n;
const gcd=(a,b)=>b?gcd(b,a%b):a;
const lcm=(a,b)=>a/gcd(a,b)*b;
const D=50n*B*w.reduce(lcm);
function g(i,j){
 if(i===j)return 0n;
 if(i>j)return g(j,i);
 return 51n*z[j]*w[j]*(D/(50n*B*w[i]))-z[j]*(D/(50n*B))-z[i]*w[i]*(D/(B*w[j]));
}
const X=w.slice(1).map((_,i)=>w.slice(1).map((_,j)=>g(i,j)-g(i,64)-g(64,j)));
function ldlsigns(X){
 const T=X.map(x=>x.slice()),s=[];let p=1n;
 for(let k=0;k<T.length;k++){
  const x=T[k][k];A(x);
  s.push(x*p>0n);
  for(let i=k+1;i<T.length;i++)for(let j=i;j<T.length;j++){
   const a=x*T[i][j]-T[i][k]*T[k][j];
   A.equal(a%p,0n);T[i][j]=T[j][i]=a/p;
  }
  p=x;
 }return s;
}
A.deepEqual(ldlsigns(X).flatMap((v,i)=>v?[i]:[]),[13]);
A(ldlsigns(X.map((r,i)=>r.map((v,j)=>1047n*D*(i===j?2n:1n)-10n*v))).every(Boolean));
const EL=27182818n,ELD=10n**7n,EU=271828183n,EUD=10n**8n;
let f=1n,t=1n;
for(let i=1n;i<=13n;i++){f*=i;t=i*t+1n;}
A(t*ELD>f*EL);A((13n*t+1n)*EUD<13n*f*EU);
f=1n;
for(let i=1n;i<=N;i++){
 f*=i;
 A(f*EU**i<((i+5n)*EUD)**i);
 if(i<=59n)A(f*(500n*EU)**i<((510n*i+1087n)*EUD)**i);
}
const U=lcm(lcm(lcm(D*65n,2n*B),50n*B*p),500n*B*B);
const div=d=>{A.equal(U%d,0n);return U/d;};
let sum=5n*p*div(B),spine=0n,ss=0n;
for(let i=0;i<M;i++){
 const a=w[i], b=z[i];
 sum+=m*(p*div(a)+b*(49n*div(2n*B)+a*div(B*p)+1087n*a*div(500n*B*B)+1413n*div(500n*B)
                    +(a>p?(a-p)*div(50n*B*p):0n)));
 for(let j=i+1;j<M;j++){
  sum+=m*m*b*a*div(B*w[j]);ss+=2n*g(i,j);
  if(j===i+1)spine+=g(i,j)*div(D);
 }
}
sum+=(m/2n)*(65n*ss*div(65n*D)+15n*1047n*div(10n));
const zz=z.reduce((a,b)=>a*b,1n)**m;
const means=[];
for(const cost of [sum,sum+spine]){
 const den=N*EL*U,num=1000n*ELD*cost;
 const q=(num+den-1n)/den;
 means.push(Number(q));
 A(8n*q**N*B**3250n<(7n*B**N*zz<<26854n));
}
A.deepEqual(means,[885325,885327]);
console.log(means);
