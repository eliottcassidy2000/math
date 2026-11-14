// node replay_fixed_triple_height.js
"use strict";
const check=(x,m)=>{if(!x)throw Error(m)};
const pw=k=>1n<<BigInt(k);
const mod=(x,m)=>(x%m+m)%m;
function g(x){let h=1n;while(x){const b=x&1n;x>>=1n;if(!b)h+=x+(x&1n);}return h;}
function rev(x,n){let y=0n;while(n--){y=(y<<1n)|(x&1n);x>>=1n;}return y;}
function join(A,B,b){return A*pw(b)+B;}
function portrait(r,k,Z,t){
 const x=r-1n,l=x.toString(2).length, A=pw(k)-1n-rev(x,k), B=x;
 const P=join(join(A,Z,t),B,k), n=2*k+t;
 const H=g(P), R=1n+P, Rc=pw(n)-rev(P,n);
 check(mod(R-r,pw(k))===0n&&mod(Rc-r,pw(k))===0n,"frozen endpoints");
 const Y=pw(t)-1n-Z,L=pw(k-l),rho=rev(x,l);
 check(mod(H-2n*g(x)-g(Y),L)===0n,"collar reduction");
 if(t>0){
  const a=rho*L,Q=pw(t)-1n-rev(Z,t);
  const val=2n*g(x)+g(Y)+a*(2n*pw(t)-Z-Q)+a*a*pw(t);
  check(mod(H-val,pw(k))===0n,"unreduced collar identity");
 }
 return {n,H,R,Rc,P,Y};
}
let cases=0;
for(const r of [3n,5n,7n,9n,13n,17n,33n,63n,65n,129n]){
 const l=(r-1n).toString(2).length;
 for(let k=l+2;k<=l+8;k++)for(let t=0;t<=8;t++)for(let z=0n;z<pw(t);z++){portrait(r,k,z,t);cases++;}
}
check(g(109n)===33n&&g(43n)===9n,"middle intercept");
for(let k=15;k<=96;k++){
 const L=pw(k-7),t=k+6;
 const Y=join(join(109n,0n,k-7),43n,6), Z=pw(t)-1n-Y;
 const {H,n}=portrait(65n,k,Z,t);
 check(n===3*k+6&&g(Y)===1199n*L-67n,"inserted middle");
 check(mod(H-63n,pw(k))===0n,"exact residue triple");
 check(H>=pw(3*k-9)+1n&&H<=pw(3*k+6),"cubic window");
}
console.log("PASS: "+cases+" general collar words and 82 explicit fixed-triple witnesses");
