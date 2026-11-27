const assert=require('node:assert/strict');
const bin=(n,k)=>{let y=1n;for(let i=1;i<=k;i++)y=y*BigInt(n-i+1)/BigInt(i);return y;};
const atom=(d,n)=>bin(d+n,d)-1n;
function apery(d){
 const D=BigInt(d), w=Array(d).fill(null), used=[];w[0]=0n;
 for(let n=2;;n++){
  const a=atom(d,n);
  if(w.every(x=>x!==null&&x<a))return {w,used,stop:n};
  const r=Number(a%D);
  if(w[r]!==null&&w[r]<=a)continue;
  used.push(a);
  for(let i=0;i<d;i++)for(let j=0;j<d;j++)if(w[j]!==null){
   const l=(j+r)%d,z=w[j]+a;
   if(w[l]===null||w[l]>z)w[l]=z;
  }
 }
}
function thresholds(d){
 const a=[];
 for(let p=2;p<=d;p++){
  let prime=true;for(let q=2;q*q<=p;q++)if(p%q===0)prime=false;
  if(prime)for(let t=p;d%t===0;t*=p)a.push([t,p]);
 }
 return a.sort((a,b)=>a[0]-b[0]);
}
let exceptions=[],families=0;
for(let d=2;d<=150;d++){
 const {w,used,stop}=apery(d),D=BigInt(d);
 const wmax=w.reduce((u,v)=>u>v?u:v),F=wmax-D,G=(w.reduce((u,v)=>u+v)-D*(D-1n)/2n)/D;
 {
  families++;
  const gs=thresholds(d).map(([t,p])=>[atom(d,t),p]);
  assert.equal(JSON.stringify(used.map(String)),JSON.stringify(gs.map(([g])=>String(g))));
  const w2=Array(d).fill(null);
  for(let r=0;r<d;r++){let q=r,z=0n;
   for(const [g,p] of gs){z+=BigInt(q%p)*g;q=Math.floor(q/p)}
   const s=Number(z%D);assert.equal(w2[s],null);w2[s]=z;
  }
  for(let r=0;r<d;r++)assert.equal(w2[r],w[r]);
  assert.equal(F,gs.reduce((s,[g,p])=>s+BigInt(p-1)*g,0n)-D);
  assert.equal(2n*G,F+1n);
 }
 if(2n*G!==F+1n)exceptions.push({d,F:String(F),G:String(G),stop,used:used.map(String),w:d<=30?w.map(String):[]});
 if([2,3,4,5,8,9,12,16,25,30].includes(d))console.log(d,used.map(String).join(','),String(F),String(G));
}
console.log('all degrees verified',families);
console.log('nonsymmetric cases',exceptions.length,'first',JSON.stringify(exceptions.slice(0,4)));
