// node replay_absolute_reflection.js
"use strict";
const check=(c,m)=>{if(!c)throw Error(m)};
const bits=(n,i)=>(n>>>i)&1;
const pop=n=>{let s=0;while(n){n&=n-1;s++}return s};
function paths(rows){
 const n=rows.length,F=(1<<n)-1;
 const E=Array.from({length:F+1},()=>Array(n).fill(0)),H=Array(F+1).fill(0);
 H[0]=1;for(let i=0;i<n;i++)E[1<<i][i]=1;
 for(let s=1;s<=F;s++){
  for(let j=0;j<n;j++)if(bits(s,j)){
   for(let i=0;i<n;i++)if(bits(rows[i],j))E[s][j]+=E[s^(1<<j)][i];
  }
  H[s]=E[s].reduce((a,b)=>a+b,0);
 }
 return {H,E,F};
}
function decorate(rows,U){
 const n=rows.length;
 return rows.map((r,i)=>r|((1-bits(U,i))<<n)).concat(U);
}

function color(rows,q) {
 const {H,F}=paths(rows);
 let v=Array(F+1).fill(0);v[0]=1;
 while(q--){
   const next=Array(F+1).fill(0);
   for(let s=0;s<=F;s++){
      for(let j=s;;j=(j-1)&s){
        next[s]+=H[j]*v[s^j];if(j===0)break;
      }
   }
   v=next;
 }
 return v[F];
}
function TT(n){return Array.from({length:n},(_,i)=>((1<<n)-1)^((2<<i)-1));}
function join(A,B){return A.map(r=>r|(((1<<B.length)-1)<<A.length)).concat(B.map(r=>r<<A.length));}
function tri(A,B){
 const a=A.length,n=a+B.length;
 return join(A,B).map((r,i)=>r|((i>=a?1:0)<<n)).concat((1<<a)-1);
}
const bank=[[[0],1,2], [TT(3),1,8],[[4,1,2],3,12],[[12,1,2,6],5,32],
 [[44,49,2,38,45,4],19,256],[[236,113,18,70,73,220,4,94],217,4096],
 [[236,241,82,6,201,220,8,76],223,4096],[[44,241,50,38,73,208,141,29],599,8192]];
for(const [r,h,f] of bank)check(color(r,1)===h&&color(r,2)===f,"dyadic palette");
for(const [A,h,f] of bank)for(const B of bank.slice(0,4).map(t=>t[0])) {
 const T=tri(A,B),g=color(B,1),z=color(B,2),N=color(T,1);
 check(N===(f-h)*(z-g)+2*h*g&&3*N-2*f*z===(f-3*h)*(z-3*g),"literal reader and factorization");
 if(T.length<=10)for(let q=1;q<=3;q++){
   const a=color(A,q),b=color(B,q);
   const p=q*(a*b+(color(A,q+1)-a)*(color(B,q+1)-b)+(color(A,q-1)-a)*(color(B,q-1)-b));
   check(color(T,q)===p,"full color transfer");
 }
}
const A=join(bank[3][0],TT(6));
check(color(tri(A,[0]),1)===2053&&color(tri(A,TT(2)),1)===6139,"same source prime and 7p");
function prime(n){if(n<2)return false;for(let j=2;j*j<=n;j++)if(n%j===0)return false;return true;}
check(prime(2053)&&prime(877)&&7*877===6139,"prime pair");
const N1=[12,17,2,38,13,23], N2=[124,9,42,80,6,90,22];
check(color(N1,1)===29&&color(N1,2)===320&&color(tri(bank[3][0],N1),1)===8147&&prime(8147),"nondyadic prime repair");
check(color(N2,1)===109&&color(N2,2)===1364&&color(tri(TT(2),N2),1)===3983&&prime(569)&&7*569===3983,"nondyadic 7p repair");
console.log("PASS: literal shifted products and q=1,2,3, five dyadic cores, and a physical prime/7p pair");
