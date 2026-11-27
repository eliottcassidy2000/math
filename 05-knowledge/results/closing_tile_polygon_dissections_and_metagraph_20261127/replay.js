const assert=require('node:assert/strict');
function mk(n,m,close){
 const a=Array(n).fill(0);let j=0;
 for(let x=0;x<n;x++)for(let y=x+1;y<n;y++){
  const b=y===x+1?0:x===0&&y===n-1?+close:m>>j++&1;
  if(b)a[y]|=1<<x;else a[x]|=1<<y;
 }return a;
}
function read(a,start=null){
 const n=a.length,U=(1<<n)-1,dp=new Int32Array((U+1)*n);
 for(let S=1;S<=U;S++)for(let J=S;J;J&=J-1){
  const j=31-Math.clz32(J&-J),R=S^(1<<j);
  let h=R?0:(start===null||start===j?1:0);
  for(let T=R;T;T&=T-1){const i=31-Math.clz32(T&-T);if(a[i]>>j&1)h+=dp[R*n+i];}
  dp[S*n+j]=h;
 }
 return [...dp.slice(U*n)];
}
const sum=(a)=>a.reduce((x,y)=>x+y,0);
function profile(a){
 const n=a.length,A=read(a,0);return [sum(read(a)),A.reduce((v,h,i)=>v+(a[i]&1)*h,0)];
}
function tri(lo,hi){
 if(hi-lo<2)return [[]];
 const out=[];
 for(let k=lo+1;k<hi;k++)for(const a of tri(lo,k))for(const b of tri(k,hi))out.push([[lo,k,hi],...a,...b]);
 return out;
}
const sign=(a,i,j)=>(a[i]>>j&1)?1:-1;
function charge(a,[i,j,k]){return (sign(a,i,j)+sign(a,j,k)+sign(a,k,i)-1)/2;}
let checks=0;
for(let n=3;n<=7;n++){
 let seen=new Map,types=new Map,examples={};
 const all=tri(0,n-1),fan=Array.from({length:n-2},(_,i)=>[0,i+1,i+2]);
 for(let m=0;m<1<<((n-1)*(n-2)/2-1);m++){
  const a=mk(n,m,0),b=mk(n,m,1),x=profile(a),y=profile(b),dh=y[0]-x[0],dc=y[1]-x[1];
  assert.equal(dc,read(a,0)[n-1]-read(a,n-1)[0]);
  const t=[Math.sign(dh),Math.sign(dc)].join();types.set(t,(types.get(t)||0)+1);
  if(dh===0&&dc!==0&&!examples.neutralPath)examples.neutralPath={m,x,y};
  if(dh!==0&&dc===0&&!examples.neutralCycle)examples.neutralCycle={m,x,y};
  if(dh*dc<0&&!examples.opposed)examples.opposed={m,x,y};
  const sig=fan.map((z)=>[sign(b,z[0],z[1]),sign(b,z[1],z[2]),sign(b,z[2],z[0])].join()).join(';');
  if(seen.has(sig)&&!examples.sameFan){const v=seen.get(sig);if(v.dh!==dh||v.dc!==dc)examples.sameFan=[v,{m,x,y,dh,dc}];}
  if(!seen.has(sig))seen.set(sig,{m,x,y,dh,dc});
  if(n<=6||m<128)for(const ds of all){assert.equal(ds.reduce((s,t)=>s+charge(a,t),0),0);assert.equal(ds.reduce((s,t)=>s+charge(b,t),0),1);checks+=2;}
 }
 console.log(JSON.stringify({n,types:[...types].sort(),examples}));
}
console.log("triangulation charge checks",checks);

const choose=(n,k)=>{let s=1;for(let i=1;i<=k;i++)s=s*(n-i+1)/i;return s;};
for(const n of [3,6,9,12]){
 let count=0;
 for(const ts of tri(0,n-1)){
  const mp=new Map;const neigh=ts.map(()=>[]);let root;
  for(let i=0;i<ts.length;i++){const [a,b,c]=ts[i];for(const [x,y] of [[a,b],[a,c],[b,c]]){
   if(x===0&&y===n-1)root=i;
   if(y===x+1||x===0&&y===n-1)continue;
   const k=n*x+y;if(mp.has(k)){const j=mp.get(k);neigh[i].push(j);neigh[j].push(i);}else mp.set(k,i);
  }}
  const visit=(i,parent,d)=>!(d%2&&neigh[i].length!==3)&&neigh[i].filter(j=>j!==parent).every(j=>visit(j,i,d+1));
  count+=+visit(root,-1,0);
 }
 const t=n/3-1;assert.equal(count,2*choose(4*t+2,t)/(4*t+2));console.log("all-cyclic compatible dissections",n,count);
}
