const assert=require('node:assert/strict');
function tournament(n, back=0){
 const o=Array(n).fill(0);let p=0;
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++,p++)
  if(back>>p&1)o[j]|=1<<i;else o[i]|=1<<j;
 return o;
}
function point(A,O){
 const n=A.length;return A.map((e,i)=>e|(((O>>i&1)^1)<<n)).concat(O);
}
function chainapex(w){
 return point(tournament(w.length),[...w].reduce((m,b,i)=>m|(Number(b)<<i),0));
}
function cover(A){
 const n=A.length,M=1<<n,df=new Int32Array(M*n),db=new Int32Array(M*n),all=new Int32Array(M);
 all[0]=1;
 for(let i=0;i<n;i++)df[(1<<i)*n+i]=db[(1<<i)*n+i]=1;
 for(let a=1;a<M;a++)for(let v=a;v;v&=v-1){
   const i=31-Math.clz32(v&-v),b=a^(1<<i);if(b)
   for(let z=b;z;z&=z-1){
     const j=31-Math.clz32(z&-z);
     if(A[j]>>i&1)df[a*n+i]+=df[b*n+j];
     if(A[i]>>j&1)db[a*n+i]+=db[b*n+j];
   }
   all[a]+=df[a*n+i];
 }
 return {n,M,df,db,all};
}
function height(A){
 // Float64 integers are exact for all small instances checked below.
 const n=A.length,M=1<<n,f=new Float64Array(M*n);
 for(let i=0;i<n;i++)f[(1<<i)*n+i]=1;
 for(let a=1;a<M-1;a++)for(let v=a;v;v&=v-1){
   const i=31-Math.clz32(v&-v),v0=f[a*n+i];if(!v0)continue;
   for(let z=A[i]&~a;z;z&=z-1){const j=31-Math.clz32(z&-z);f[(a|(1<<j))*n+j]+=v0}
 }
 return f.slice((M-1)*n).reduce((a,b)=>a+b,0)
}
function powersums(A){
 const {M,n,db,all}=cover(A),h=all[M-1],out=[];
 for(let O=0;O<M;O++){
  let val=0;for(let i=0;i<n;i++)if(O>>i&1)val+=db[(M-1)*n+i];
  if(out[val]===undefined)out[val]=O;
 }
 assert.equal(out.length,h+1);
 for(let i=0;i<=h;i++)assert.notEqual(out[i],undefined)
 return out;
}
function sumEnds(A){
 const {M,n,df}=cover(A);let p=0;
 return Array.from({length:n},(_,i)=>df[(M-1)*n+i]);
}
function join(...parts){
 const n=parts.flat().length;let pos=0,out=[];
 for(const p of parts){
  const tail=((1<<n)-1)^((1<<(pos+p.length))-1);
  for(const r of p)out.push((r<<pos)|tail);
  pos+=p.length;
 }
 return out;
}
const Y=[tournament(6,177),tournament(6,181)], OY=19;
const U=tournament(3,2), V=chainapex('100'), OV=powersums(V)[2];
const EU=sumEnds(U),OU=((1<<U.length)-1)^3;
assert.equal(EU[0]+EU[1],2)
const ex=Y.map(y=>height(point(join(U,y,V),OU|(OY<<3)|(OV<<9))));
assert.deepEqual(ex,[33575,33577]);console.log('literal 14-site heights',ex);
for(let d=1;d<=24;d++){
 const L=Math.max(1,Math.ceil(Math.log2(Math.max(1,d-1)))), V=chainapex('1'+'0'.repeat(L)), O=powersums(V)[d];
 const [g,h]=Y.map(y=>height(point(join(y,V),OY|(O<<6))));
 assert.equal(h-g,2*d);
}
for(let L=2;L<=8;L++)for(let bits=0;bits<(1<<(L-2));bits++){
 const w='1'+(L==2?'':bits.toString(2).padStart(L-2,'0'))+'0';
 powersums(chainapex(w));
}
function capRatio(x,r){
 let n=0,d=1;
 for(let a=1;a<=x;a++){
  const b=Math.min(x,Math.floor(a*(a-1)/(r*(r-1))));
  if(b*d>n*a){n=b;d=a}
 }
 return [n,d]
}
function exclude(M,start){
 let checked=0;
 for(let r=start;r<=56;r++)
  for(let x=Math.max(2*r+1,r*(115-3*r));x<=Math.floor(M/r);x++){
   const [p,q]=capRatio(x,r),e=Math.floor((50*x*(3250-x)+3249)/3250);
   assert((e-r*x)*q>p*(M-r*x),JSON.stringify({M,r,x}));
   checked++;
  }
 return checked;
}
assert.equal(exclude(1408,2),639);
assert.equal(exclude(2824,3),984);
console.log('signed-support distance tests: 639 and 984');
