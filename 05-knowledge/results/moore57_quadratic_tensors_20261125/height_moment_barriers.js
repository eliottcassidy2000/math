
const model=[[0,0,4],[0,2,4],[0,3,4],[1,0,14],[1,2,2],[2,0,3],[2,2,1]];
const gcd=(a,b)=>b?gcd(b,a%b):a;
const totals=[2,3,4,5].map(k=>model.reduce((s,[r,p,m])=>s+m*k**r*gcd(k,p||1)**2,0));
if(totals.join()!=="96,128,224,192")throw Error("Selmer totals");
if(model.reduce((s,[r,p,m])=>s+m*r,0)!==24)throw Error("rank total");
if(model.reduce((s,[r,p,m])=>s+m*(r>=2),0)!==4)throw Error("tail");
if(model.reduce((s,[r,p,m])=>s+m*(r%2),0)!==16)throw Error("odd");
const z=(q)=>{let s=0;for(let n=1;n<=1e6;n++)s+=n**(-q);return s};
const c=5*z(10)/(4*(Math.PI**2/6)*z(5));
console.log("Moment totals:",totals,"; small-base coefficient:",c);

function iroot(X,p){
 let a=0n,b=1n;
 const k=BigInt(p);
 while(b**k<=X)b*=2n;
 while(b-a>1n){let m=(a+b)/2n;if(m**k<=X)a=m;else b=m;}
 return a;
}
function sieve(n){
 const t=Array(n+1).fill(1),p=Array(n+1).fill(true);
 for(let k=2;k<=n;k++)if(p[k]){
  for(let j=k;j<=n;j+=k){t[j]=-t[j];if(j>k)p[j]=false;}
  for(let j=k*k;j<=n;j+=k*k)t[j]=0;
 }
 return t;
}
function base(x){
 const a=iroot(x,3),b=iroot(x,2),d=Number(iroot(x,6)), mu=sieve(d);
 let Y=0n;
 for(let i=1;i<=d;i++){const k=BigInt(i);Y+=BigInt(mu[i])*(a/(k*k))*(b/(k*k*k));}
 return 2n*Y-(x>=27n?1n:0n);
}
for(const X of [10n**6n,10n**12n,10n**18n]){
 const D=Number(iroot(X,6)),mu=sieve(D),B=iroot(X,2);
 let all=0n,small=0n;
 for(let d=1;d<=D;d++)if(mu[d]){
  const T=X/BigInt(d)**6n;
  all+=2n*base(T);small+=2n*base(T<B?T:B);
 }
 const scaled=Number(small)/Number(all)/((Number(B)/Number(X))**(2/3));
 console.log(String(X),String(small),String(all),scaled);
}