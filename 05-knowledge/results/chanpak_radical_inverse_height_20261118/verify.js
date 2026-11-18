"use strict";
function check(v,msg){if(!v)throw Error(msg);}
function F(u){
 let h=1n;
 for(;u>0n;u>>=1n)if(!(u&1n)){let q=u>>1n;h+=q+(q&1n);}
 return h;
}
for(let u=1n;u<20000n;u+=2n){
 let w=[...u.toString(2)].reverse().map(BigInt),r=w.length;
 let s=0n,B=0n;
 for(let i=0;i<r;i++)if(w[i]){
  B+=1n<<BigInt(r-1-i);
  for(let j=0;j<i;j++)if(w[j])s+=1n<<BigInt(i-j-1);
  if(i&&w[i-1])s++;
 }
 const z=(u*B)-(s<<BigInt(r));
 check(s===u-F(u)&&z>0n&&z<=(BigInt(r)<<BigInt(r)),"scaled binary identity");
}
const X=10**7,Ft=new Uint32Array(1+Math.floor(X/24));
Ft[0]=1;
for(let u=1;u<Ft.length;u++){let b=u>>1;Ft[u]=Ft[b]+((u&1)?0:b+(b&1));}
let sum=0;
for(let k=0;k<=30;k++){
 const h=12*2**k-1,R=36*3**k-h;
 const c=Math.log1p(h/(2*R))/h; sum+=c;
 if(k<4){
  let actual=0;
  for(let u=1;24*u*3**k<=X;u+=2)if(R*u+h*Ft[u]+2*h-5<=X)actual++;
  check(Math.abs(actual-X*c)<9,"finite recurrence control");
  console.log(k,actual,X*c);
 }
}
check(Math.abs(2*sum-0.052301466288545)<1e-13,"series");
const H=u=>25n*u+11n*F(u)+17n;
for(let L=1;L<=9;L++){
 let seen=new Set();
 for(let a=0;a<1<<(L-1);a++){
  let u=1n,r=1;
  for(let b=1;b<L;b++)if(a>>(b-1)&1){u+=1n<<BigInt(r+b);r+=b+1;}
  check(u.toString(2).length===r,"width law");
  seen.add(H(u)% (1n<<BigInt(L)));
 }
 check(seen.size===1<<(L-1),"decision bijection");
}
console.log("PASS: binary height identity, address tallies, and greedy width distribution.");
