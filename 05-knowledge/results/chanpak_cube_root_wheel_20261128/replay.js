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
 let d=2;
 if(m>1n){let z=4n%m;while(z!==1n){z=4n*z%m;d+=2;}}
 const G=gcd(M,116n),J=Number(M/G-1n),w=d+4+2*(e+2),K=7+J*w,P=1n<<BigInt(K);
 const bit=i=>1n<<BigInt(w-1-i);
 const v0=bit(e+3)|bit(e+4),v1=bit(e+3)|bit(e+4+d);
 const word=(j,tail)=>{let v=0n;for(let i=0;i<J;i++)v=(v<<BigInt(w))+(i<j?v1:v0);return (v<<5n)+tail;};
 return {G,J,K,P,word,I:29n*(4n-(1n<<BigInt(d+1)))};
}
function decode(N,M,R,C){
 const {G,J,K,P,word,I}=C;
 assert(N>=P**6n && (M%2n===1n||R%2n!==0n));
 const D=1n<<BigInt(Math.floor((N.toString(2).length-1)/6)),T=N+100000n*P*rootCeil(N,3n);
 const A=2560n,B=603n,Cc=735n,F=A*(2n*T-280n)+B*Cc;
 const m0=(rootCeil(4n*F,2n)-B-Cc+A-1n)/A;
 let m=m0+mod(D-1n-m0,D);
 const use29=M%29n===0n;
 if(use29)while(mod(111n+353n*m-R,29n))m+=D;
 const S=A*m+B+Cc,delta=rootCeil(S*S-4n*F,2n);
 const v0=word(0,0n);
 const lattice=use29?29n:1n;
 const u0=use29?mod(-v0*inv(P,29n),29n):0n;
 const cell=P*u0+v0,step=lattice*P;
 const z=(S-2n*Cc-2n*A*cell+A*step-delta)/(2n*A*step);
 const u=u0+lattice*z;
 let a=P*u+v0, b=m-a, tail=0n;
 if(M%4n===0n && mod(L(a,b)-R,4n))tail=29n;
 a=P*u+word(0,tail);b=m-a;
 const diff=R-L(a,b);assert(diff%G===0n);
 const j=J===0?0:Number(mod(diff/G*inv(I/G,M/G),M/G));
 a=P*u+word(j,tail);b=m-a;const h=L(a,b);
 assert(a>0n&&b>0n&&mod(h-R,M)===0n);
 assert(h>N&&(h+2n-N)**3n<(200000n*P)**3n*N);
 return {j,h};
}
for(let M=1n;M<=100n;M++){
 const C=prep(M),{J,P,word,I}=C;
 assert(gcd(M,I)===C.G);
 for(const [u,w] of [[0n,0n],[1n,2n],[9n,7n]])for(const tail of [0n,29n]){
   let v=word(0,tail),base=L(P*u+v,P*w+P-1n-v);
   for(let j=1;j<=Math.min(J,3);j++){
     v=word(j,tail);
     assert(mod(L(P*u+v,P*w+P-1n-v)-base-BigInt(j)*I,M)===0n);
   }
 }
}
const targets=[1n,2n,4n,7n,8n,12n,29n,58n,116n,792n,841n];
let largest=0,checks=0;
for(const M of targets){
 const C=prep(M),N=C.P**6n;
 const choices=M<=116n?[1n,3n,5n,17n]:[1n,5n,17n,337n];
 for(const R of choices){
   decode(N+R,M,R,C);
   checks++;
 }
 largest=Math.max(largest,N.toString(2).length);
}
console.log({checks,largest,switchRange:100});
