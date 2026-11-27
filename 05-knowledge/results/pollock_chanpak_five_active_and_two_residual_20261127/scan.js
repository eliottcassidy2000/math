function assert(x){if(!x)throw Error("failed replay")}
const F=x=>x*(x-1)*(x+1)/6;
function H(v){let p=0,q=0;for(const b of v.toString(2)){if(b==="1")q+=Math.ceil(p/2);p=2*p+Number(b)}return 1662*v-58*q+427;}
function payloads(hi,lo=0){const r=new Map();for(let v=Math.max(1,Math.ceil((lo-427)/1662));1604*v+427<=hi;v++)if(v.toString(2).startsWith("1010")){const h=H(v);if(lo<=h&&h<=hi){if(!r.has(h))r.set(h,[]);r.get(h).push(v)}}return r}
function solver(max){
 const f=[0,0];for(let i=2;F(i)<=max;i++)f.push(F(i));
 const bank=new Map();
 for(let a=2;a<f.length&&f[a]*2<=max;a++)for(let b=a;b<f.length&&f[a]+f[b]<=max;b++){
   const t=f[a]+f[b];if(!bank.has(t))bank.set(t,[]);bank.get(t).push([a,b]);
 }
 return n=>{
  const out=[];
  for(let a=2;5*f[a]<=n;a++)for(let b=a;f[a]+4*f[b]<=n;b++)for(let c=b;f[a]+f[b]+3*f[c]<=n;c++){
   const row=bank.get(n-f[a]-f[b]-f[c]);if(row)for(const [d,e] of row)if(d>=c)out.push([a,b,c,d,e]);
  }
  return out
 }
}
const key=a=>a.reduce((s,i)=>s+i,0)+":"+a.reduce((s,i)=>s+i*i,0);
function matching(A,B){
 const m=new Map(),out=[];
 for(const a of A){const k=key(a);if(!m.has(k))m.set(k,[]);m.get(k).push(a)}
 for(const b of B)for(const a of (m.get(key(b))||[]))if(a.every(x=>!b.includes(x)))out.push([a,b]);
 return out
}
const chart=payloads(500001), small=solver(250001), hits=[];
for(const [h,vs] of [...chart].sort((a,b)=>a[0]-b[0])){
 const n=(h-1)/2;for(const t of matching(small(n),small(n+1)))hits.push([h,vs,...t]);
}
assert(chart.size===31&&hits.map(x=>x[0]).join(",")==="273689,273689,279295");
console.log({chart:chart.size,hits});
for(let n=0;n<=124;n++){const h=matching(small(n),small(n+1));assert(Boolean(h.length)===(n===124));if(h.length)console.log({firstUnrestricted:n,h});}
const big=273299016,bigH=2*big+1, bigV=payloads(bigH,bigH).get(bigH), solveBig=solver(big+1), A=solveBig(big), B=solveBig(big+1), bigHits=matching(A,B);
assert(A.length===59326&&B.length===27104&&bigHits.length===1&&bigV.join(",")==="329528,333037");
console.log({bigH,bigV,lowerFibre:A.length,upperFibre:B.length,bigHits});
function cubic(t){
 const z=t*t;
 if(t%2===0)return [t**3,2*(z-t+1),z/2,t*(t+1)];
 if(t%4===1)return [t**3,2*z,t*(t+1)/2,z-t+1];
 return [t**3-1,2*z,t*(t-1)/2,z+t+1];
}
for(let n=11;n<=10000;n++){
 const t=BigInt(n),[k,S,a,D]=cubic(n).map(BigInt), f=x=>x*(x*x-1n)/6n;
 const x=(S+D-a)/2n,y=(S-D+a)/2n,z=t*t/2n;
 assert(x+y===S&&x-y+a===D&&x>=y&&y-a>=2n&&a*S*D===k*(k+1n));
 assert(f(x)+f(y)+2n*f(z)<(k+1n)*(k+2n)/2n);
 assert(f(k+1n)+f(x)+f(y)===f(k)+f(x+a)+f(y-a));
}
console.log("Cubic shell identities: all t=11..10000");
