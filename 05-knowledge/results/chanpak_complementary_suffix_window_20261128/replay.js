"use strict";
const assert=require("node:assert/strict");

function f(n){
  let p=0n,z=0n;
  for(const c of n.toString(2)){
    const d=BigInt(c);
    if(!d)z+=p+p%2n;
    p=2n*p+d;
  }
  return z;
}
function ceilroot(n,k){
  if(n<2n)return n;
  let lo=0n,hi=1n<<BigInt(Math.ceil(n.toString(2).length/Number(k)));
  while(hi-lo>1n){const mid=(lo+hi)/2n;if(mid**k>=n)hi=mid;else lo=mid;}
  return hi;
}
const L=(x,y)=>111n+129n*x+162n*y+320n*x*y+29n*(f(x)+f(y));
function decode(N){
  assert(N>=10n**12n);
  const D=1n<<BigInt(Math.floor((N.toString(2).length-1)/6)-4);
  const T=N+80n*ceilroot(N,3n);
  const A=2560n,B=603n,C=735n,M=A*(2n*T-280n)+B*C;
  const m0=(ceilroot(4n*M,2n)-B-C+A-1n)/A;
  const m=m0+((D-1n-m0)%D+D)%D;
  const S=A*m+B+C;
  const a=(S-2n*C+A-ceilroot(S*S-4n*M,2n))/(2n*A),b=m-a;
  assert(a>0n&&b>0n&&(a+b)%D===D-1n);
  const x=2n*a,y=2n*b,h=L(x,y);
  assert(h>N&&(h+2n-N)**3n<160n**3n*N);
  assert((h>T?h-T:T-h)**3n<71n**3n*N);
  return {x,y,h,D};
}
function suffix(k,u,v){
  const D=1n<<BigInt(k), rev=BigInt("0b"+v.toString(2).padStart(k,"0").split("").reverse().join(""));
  const Y=D-1n-rev;
  const e=(v<D/2n?1n:0n)*(u%2n);
  assert.equal(f(D*u+v),f(u)+u*Y+f(v)+e);
}
for(let k=1;k<=8;k++){
  const D=1n<<BigInt(k);
  for(let v=0n;v<D;v++)for(let u=0n;u<4n;u++)suffix(k,u,v);
}
for(let n=10n**12n;n<10n**12n+10000n;n++) decode(n);
for(let k=40n;k<=1540n;k++){
 const n=2n**k;
 for(let d of [-2n,-1n,0n,1n,2n]) decode(n+d);
}
function firstShell(h){
 const N=ceilroot(h**5n,3n);
 assert(N**3n>=h**5n && (N-1n)**3n<h**5n);
 const o=decode(N);
 assert(o.h**3n>h**5n && (o.h+3n)**3n<(h+1n)**5n);
 return {parent:h,x:o.x,y:o.y,counts:[o.h,o.h+2n],n:o.x.toString(2).length+o.y.toString(2).length+7};
}
console.log(decode(10n**12n));
console.log(firstShell(2n**63n+1n));
console.log("2,040 suffix identities and 17,505 targets passed");

function decode29(N,R){
  const D=1n<<BigInt(Math.floor((N.toString(2).length-1)/6)-4);
  const T=N+12000n*ceilroot(N,3n);
  const A=2560n,B=603n,C=735n,M=A*(2n*T-280n)+B*C;
  const m0=(ceilroot(4n*M,2n)-B-C+A-1n)/A;
  let m=m0+((D-1n-m0)%D+D)%D;
  while((111n+353n*m-R)%29n!==0n)m+=D;
  const S=A*m+B+C;
  const a=29n*((S-2n*C+29n*A-ceilroot(S*S-4n*M,2n))/(58n*A)),b=m-a;
  assert(a>0n&&b>0n&&(a+b)%D===D-1n);
  const h=L(2n*a,2n*b);
  assert(h>N&&(h+2n-N)**3n<25000n**3n*N&&(h-R)%29n===0n);
}
for(const N of [10n**12n,10n**12n+17n,2n**61n,10n**100n])
  for(let R=0n;R<29n;R++)decode29(N,R);
console.log("116 exact residue-and-window checks passed");
