"use strict";
function ck(x,msg){if(!x)throw Error(msg)}


function h(a){let n=a.length,N=1<<n,dp=new Int32Array(N*n);for(let i=0;i<n;i++)dp[(1<<i)*n+i]=1;for(let s=1;s<N;s++)for(let i=0;i<n;i++){let w=dp[s*n+i];if(!w)continue;let b=a[i]&(N-1-s);while(b){let x=b&-b;b-=x;dp[(s|x)*n+31-Math.clz32(x)]+=w}}let z=0;for(let i=0;i<n;i++)z+=dp[(N-1)*n+i];return z}
function graph(n,p){let a=Array(n).fill(0);for(let j=1;j<n;j++)for(let i=0;i<j;i++){if(p&1)a[i]|=1<<j;else a[j]|=1<<i;p>>>=1;}return a}
function blue(n){let orb=[];for(let i=0;i<n;i++)for(let j=i+2;j<n;j++){let k=n-1-j,l=n-1-i;if(i*n+j<=k*n+l)orb.push(i==k?[[i,j]]:[[i,j],[k,l]])}let frames=[];for(let m=0;m<(1<<orb.length);m++){let a=Array.from({length:n},(_,i)=>i+1<n?1<<(i+1):0);orb.forEach((o,u)=>o.forEach(([i,j])=>{if((m>>u)&1)a[i]|=1<<j;else a[j]|=1<<i}));frames.push(a)}return frames}
function mate(a){let n=a.length;return a.map((x,i)=>((1<<n)-1)^x^(1<<i)^(i>0?1<<(i-1):0)^(i+1<n?1<<(i+1):0));}

function gcd(a,b){while(b){let t=a%b;a=b;b=t}return a<0n?-a:a;}
const F=(a,b=1n)=>{let g=gcd(a,b);return[a/g,b/g]},A=([a,b],[c,d])=>F(a*d+b*c,b*d),D=([a,b],[c,d])=>A([a,b],[-c,d]),M=([a,b],c)=>F(a*c,b),L=([a,b],[c,d])=>a*d<c*b;
function at(a,m){let x=[0n,1n];for(let i=0;i<m;i++)x=A(x,F(i%2?-1n:1n,BigInt(2*i+1)*a**BigInt(2*i+1)));let y=A(x,F(m%2?-1n:1n,BigInt(2*m+1)*a**BigInt(2*m+1)));return L(x,y)?[x,y]:[y,x]}
let [u,v]=at(5n,12),[w,z]=at(239n,12);let lo=D(M(u,16n),M(z,4n)),hi=D(M(v,16n),M(w,4n));let q0=F(355n,113n),dl=D(q0,hi),dh=D(q0,lo);

const expected={"41":[42,40,35,68,111,72,7],"99":[12,105,50,68,75,89,5],"113":[44,37,88,2,75,28,43],"223":[84,197,48,55,34,3,60,125],"311":[38,228,248,227,203,80,129,33],"349":[156,213,216,18,64,31,169,56],"355":[248,137,203,64,46,14,50,120],"541":[244,41,202,193,78,156,162,18],"599":[150,152,234,193,76,27,35,112],"977":[446,416,10,482,238,132,39,324,116]};
for(const [n,a] of Object.entries(expected)){
  ck(a.every((r,i)=>!(r&(1<<i))&&r<(1<<a.length)),"row");
  for(let i=0;i<a.length;i++)for(let j=0;j<i;j++)
    ck(((a[i]>>j)&1)+((a[j]>>i)&1)==1,"tournament");
  ck(h(a)===+n,"literal witness "+n);
}
ck(349*599==209051&&3*41*541==66543&&h([2,4,1])==3,"product");
const f=blue(7), seed=[46,68,8,82,103,78,5],b=mate(seed);
ck(f.length===512&&f.some(a=>String(a)===String(seed))&&f.some(a=>String(a)===String(b)),"blue symmetry");
ck(h(seed)===69&&h(b)===71&&String(b)==="82,60,121,49,32,65,26","native repair");
ck(String(f.filter(a=>h(a)===71).map(a=>h(mate(a))).sort((a,b)=>a-b))==="13,45,69,93,133,137,145","mate list");
let mx=0;for(let i=0;i<32768;i++)mx=Math.max(mx,h(graph(6,i)));
ck(mx===45,"birth bound");
ck(L(F(311n,99n),lo)&&L(hi,q0)&&L(F(977n,311n),lo),"initial bracket and next lower");
ck(L(F(208341n,66317n),lo)&&L(hi,F(209051n,66543n))&&L(M(dh,66543n*113n),F(4n)),"next upper");
ck(355n*99n-311n*113n===2n&&333n*99n-311n*106n===1n&&355n*106n-333n*113n===1n,"Farey determinants");
console.log("PASS: 32768 order-six counts, 512 blue frames, literal ratio factors and exact Machin intervals.");
