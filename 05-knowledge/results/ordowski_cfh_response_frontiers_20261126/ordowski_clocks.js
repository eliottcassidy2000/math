"use strict";

const H0=[[23n,24n],[25n,28n],[29n,31n],[59n,87n],[71n,72n],[103n,129n],[143n,173n]];
const mod=(a,b)=>(a%b+b)%b;
function gcd(a,b){while(b){const r=a%b;a=b;b=r;}return a;}
function eg(a,b){if(b===0n)return [a,1n,0n];const [g,x,y]=eg(b,a%b);return[g,y,x-(a/b)*y];}
function crt(a,m,b,n){
 const [g,x]=eg(m,n), k=b-a;
 if(k%g)return null;
 const M=m*(n/g);return [mod(a+m*mod(k/g*x,n/g),M),M];
}
function order(a,m){let t=1n,k=0;do{t=t*a%m;k++;}while(t!==1n);return k;}
function pow(a,b,m){let t=1n;for(;b;b>>=1n,a=a*a%m)if(b&1n)t=t*a%m;return t;}
function local(H,R){
 const d=order(2n,H), opt=[];
 for(let i=0;i<d;i++){
  const z=R*pow(2n,BigInt(i),H)%H;
  if(z===2n||z===H-2n)opt.push([i,z===2n?2:-2]);
 }
 return {H,R,d,opt};
}
function solve(rows){
 let seen=[[0n,1n]];
 for(const x of rows){
  seen=seen.flatMap(([i,d])=>x.opt.flatMap(([j])=>{
   const r=crt(i,d,BigInt(j),BigInt(x.d));return r?[r]:[];
  }));
 }
 return seen;
}
function valuations(p,a){
 const d=order(2n,p), t=2n**BigInt(d)-1n;
 let u=t,v=0;while(!(u%p)){u/=p;v++;}
 const pm=p**BigInt(a), o=d*Number(p**BigInt(Math.max(0,a-v)));
 if(pow(2n,BigInt(o),pm)!==1n)throw Error("order value");
 for(let k=2;k<=o;k++)if(o%k===0&&pow(2n,BigInt(o/k),pm)===1n)throw Error("order minimal");
 return [Number(p),a,d,v,o];
}
function run(){
 const ds=H0.map(([a,b])=>local(a,b)), max=[];
 let best=0;
 for(let t=1;t<(1<<ds.length);t++){
  const rows=ds.filter((_,i)=>t>>i&1), m=rows.length,z=solve(rows);
  if(!z.length||m<best)continue;
  if(m>best){best=m;max.length=0;}
  max.push({moduli:rows.map(x=>Number(x.H)),clock:z.map(([u,d])=>[String(u===0n?d:u),String(d)]).sort()});
 }
 let magnitudes=[[0n,1n]];
 for(const [h,r] of H0)magnitudes=magnitudes.flatMap(([a,d])=>[r,-r].map(b=>crt(a,d,b,h)));
 const smallest=magnitudes.map(x=>x[0]).sort((a,b)=>a<b?-1:1)[0], modM=magnitudes[0][1];
 const phase=ds.reduce((a,b)=>a*BigInt(b.d)/gcd(a,BigInt(b.d)),1n);
 let lifted=[];
 for(const [p,a] of [[5n,4],[11n,3],[13n,3],[71n,2],[1093n,2]])lifted.push(valuations(p,a));
 let neg=[0n,1n];
 for(const [h,r] of H0)neg=crt(neg[0],neg[1],-r,h);
 function F(x){let z=1n;while(x){const d=x>>1n;if(!(x&1n))z+=d+(d&1n);x=d;}return z;}
 const physical=H0.map(([h,R])=>{
   const beta=(2n*(R+neg[0]))/h+2n; if(beta&1n)throw Error("parity");
   const c=4n*beta+1n,s=2n*R+h*(F(c+1n)-F(c-1n));
   if(s!==-2n*neg[0])throw Error("not the shared shift");
   return [String(h),String(c),c.toString(2).length];
 });
 if(mod(pow(2n,364n,1093n**3n)-1n,1093n**3n)===0n)throw Error("valuation");
 if(pow(1095n,364n,1093n**2n)===1n)throw Error("bad witness");
 return {negativeCommonM:String(neg[0]),physical,table:ds.map(x=>[Number(x.H),Number(x.R),x.d,x.opt]),best,max,
   smallest:String(smallest),amplitude:String(2n*smallest),mModulus:String(modM),commonPeriod:String(phase),lifted};
}

const out=run();
if(out.best!==5||out.max.length!==1||out.max[0].clock[0][0]!=="1091454")throw Error("CRT");
console.log(JSON.stringify(out,null,2));
