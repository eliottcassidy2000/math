"use strict";
const assert=require("node:assert/strict");
const mod=(a,b)=>((a%b)+b)%b;
function gcd(a,b){while(b){[a,b]=[b,a%b];}return a<0n?-a:a;}
function inv(a,b){
 let x=1n,y=0n,s=mod(a,b),t=b;
 while(t){const k=s/t;[s,t]=[t,s-k*t];[x,y]=[y,x-k*y];}
 assert(s===1n);return mod(x,b);
}
function rootCeil(n,k){
 if(n<2n)return n;
 let x=1n<<BigInt(Math.ceil(n.toString(2).length/Number(k)));
 for(;;){const y=((k-1n)*x+n/x**(k-1n))/k;if(y>=x)break;x=y;}
 while(x**k<n)x++;return x;
}
function f(n){let p=0n,z=0n;for(const c of n.toString(2)){const d=BigInt(c);if(!d)z+=p+p%2n;p=2n*p+d;}return z;}
function L(a,b){return 111n+1280n*a*b+287n*a+353n*b+29n*(f(a)+f(b)+a%2n+b%2n);}

function prep(M){
 let e=0,m=M;while(m%2n===0n){m/=2n;e++;}
 const g=gcd(m,87n), z=m/g;
 const n=z===1n?0:(z-1n).toString(2).length;
 let period=2;
 if(m>1n){let v=4n%m;while(v!==1n){v=4n*v%m;period+=2;}}
 const d=period*Math.ceil((2*e+n+6)/period);
 const width=d+3+2*(e+2),ops=[], groups=[], bits=[];
 const tw=l=>1n<<BigInt(l);
 function block(l,reverse=false){
   const bit=i=>tw(width-1-i);
   const fixed=bit(e+3+l);
   let left=fixed|bit(e+3),right=fixed|bit(e+3+d);
   if(reverse)[left,right]=[right,left];
   const k=ops.length;ops.push([left,right]);return k;
 }
 const c3=m%3n===0n?[block(1),block(1)]:[];
 for(let l=2;l<e;l++)bits.push([l,block(l)]);
 for(let i=0;i<n;i++){const l=e+2+i;groups.push([block(l+1),block(l+1),block(l,true)]);}
 const K=7+width*ops.length,P=tw(K);
 const word=(on,tail)=>{let v=0n;for(let i=0;i<ops.length;i++)v=(v<<BigInt(width))+ops[i][on.has(i)?1:0];return (v<<5n)+tail;};
 const I=l=>29n*(tw(l)-tw(d-l));
 return {P,K,e,m,g,z,c3,bits,groups,word,I,adj:29n*(4n-tw(d-1))};
}
function decode(N,M,R,C){
 const {P,e,m:odd,g,z:oddTarget,c3,bits,groups,word,I,adj}=C;
 assert(N>=P**6n && (M%2n===1n||R%2n!==0n));
 const D=1n<<BigInt(Math.floor((N.toString(2).length-1)/6)),T=N+100000n*P*rootCeil(N,3n);
 const A=2560n,B=603n,Cc=735n,F=A*(2n*T-280n)+B*Cc;
 const m0=(rootCeil(4n*F,2n)-B-Cc+A-1n)/A;
 let m=m0+mod(D-1n-m0,D);
 const use29=M%29n===0n;
 if(use29)while(mod(111n+353n*m-R,29n))m+=D;
 const S=A*m+B+Cc,delta=rootCeil(S*S-4n*F,2n);
 const on=new Set(),v0=word(on,0n);
 const lattice=use29?29n:1n;
 const u0=use29?mod(-v0*inv(P,29n),29n):0n;
 const cell=P*u0+v0,step=lattice*P;
 const z=(S-2n*Cc-2n*A*cell+A*step-delta)/(2n*A*step);
 const u=u0+lattice*z;
 let a=P*u+v0, b=m-a, tail=0n;
 if(M%4n===0n && mod(L(a,b)-R,4n))tail=29n;
 a=P*u+word(on,tail);b=m-a;
 let base=L(a,b);
 if(c3.length){
   assert(mod(adj,3n)===1n);const k=Number(mod(R-base,3n));
   for(let i=0;i<k;i++){on.add(c3[i]);base+=adj;}
 }
 for(const [l,k] of bits){
   const divisor=1n<<BigInt(l+1);
   if(mod(R-base,divisor)){on.add(k);base+=I(l);}
 }
 assert(mod(R-base,1n<<BigInt(e))===0n && (R-base)%g===0n);
 if(oddTarget>1n){
   const slope=(87n*(1n<<BigInt(e+2)))/g;
   const target=mod((R-base)/g*inv(slope,oddTarget),oddTarget);
   for(let i=0;i<groups.length;i++)if(target>>BigInt(i)&1n)
     for(const k of groups[i])on.add(k);
 }
 a=P*u+word(on,tail);b=m-a;const h=L(a,b);
 assert(a>0n&&b>0n&&mod(h-R,M)===0n);
 assert(h>N&&(h+2n-N)**3n<(200000n*P)**3n*N);
 return {a,b,h};
}

let checks=0;const samples={};
for(const M of [...Array.from({length:90},(_,i)=>BigInt(i+1)),116n,232n,792n,841n,2310n]){
 const C=prep(M),N=C.P**6n;
 for(const R of [1n,37n,-1n]){
    const r=decode(N+abs(R),M,R,C);
    if(R===-1n){assert(gcd(r.h*(r.h+2n),M)===1n);}
    checks++;
 }
 if(M>=792n)samples[String(M)]={K:C.K,thresholdBits:6*C.K+1};
}
function abs(x){return x<0n?-x:x;}
assert(checks===285 && samples["792"].K===910 && samples["841"].K===12292);
console.log({checks,samples});
