const assert=require('node:assert/strict');
function word(v) {return v.toString(2);}
function h(v) {if(!v)return 1;let u=v>>1;return h(u)+(v&1?0:u+(u&1));}
const charts=[{m:[177,181],out:19,R:158,J:162},{m:[11524,11588],out:7,R:191,J:129}];
function formula(B,i){
 const g=charts[i];
 return (640+g.J)*B+29*(h(B)-1)+2*g.R+111;
}
function graph(B,i,endpoint,K=1){
 const z=word(B),g=charts[i],n=z.length+8+K,m=n-1, v=Array(n).fill(0);
 for(let a=0;a<m;a++)for(let b=a+1;b<m;b++)v[a]|=1<<b;
 const off=z.length;
 for(let a=0,t=0;a<6;a++)for(let b=a+1;b<6;b++,t++)if(g.m[endpoint]>>t&1){
   v[a+off]^=1<<(b+off);v[b+off]^=1<<(a+off);
 }
 for(let a=0;a<m;a++){
  let bit=a<off?+z[a]:a<off+6?g.out>>(a-off)&1:a<m-1?1:0;
  if(bit)v[m]|=1<<a;else v[a]|=1<<m;
 }
 return v;
}
function dp(A){
 const n=A.length,M=1<<n,f=new Float64Array(M*n);
 for(let i=0;i<n;i++)f[(1<<i)*n+i]=1;
 for(let a=1;a<M-1;a++)for(let v=a;v;v&=v-1){
  const i=31-Math.clz32(v&-v),p=f[a*n+i];if(!p)continue;
  for(let w=A[i]&~a;w;w&=w-1){
   const j=31-Math.clz32(w&-w),k=(a|(1<<j))*n+j;
   f[k]+=p;
   if(!Number.isSafeInteger(f[k]))throw Error('overflow')
  }
 }
 return f.slice((M-1)*n).reduce((a,b)=>a+b,0)
}
function prime(v){if(v<2)return false;for(let d=2;d*d<=v;d++)if(v%d===0)return false;return true}

function graph8(B,O,endpoint,tail){
 const w=word(B),n=w.length+9+tail.length,m=n-1,A=Array(n).fill(0),off=w.length;
 for(let i=0;i<m;i++)for(let j=i+1;j<m;j++)A[i]|=1<<j;
 const flip=(i,j)=>{A[i+off]^=1<<(j+off);A[j+off]^=1<<(i+off);}
 let k=0;
 for(let i=0;i<6;i++)for(let j=i+1;j<6;j++,k++)if((endpoint?181:177)>>k&1)flip(i,j);
 for(let i=0;i<6;i++)for(let j=6;j<8;j++)if(189>>(2*i+j-6)&1)flip(i,j);
 for(let i=0;i<m;i++){
   const bit=i<off?+w[i]:i<off+8?O>>(i-off)&1:+tail[i-off-8];
   if(bit)A[m]|=1<<i;else A[i]|=1<<m;
 }
 return A;
}
const q=v=>{let a=0n;for(let u=v;u;u>>=1n)if(u&1n)a+=((u>>1n)+1n)/2n;return a;}
const P=(v,u)=>35232n*v*u+13536n*v+5382n*u-1302n*(q(v)+q(u))+2313n;
assert.deepEqual([0,1].map(e=>dp(graph8(10,36,e,'10'))),[250233,250235]);
assert.deepEqual([0,1].map(e=>dp(graph8(10,219,e,'10'))),[217619,217617]);
for(let L=1;L<=4;L++)for(let t=1<<(L-1);t<(1<<L);t++){
 const w=word(t),u=BigInt('0b'+[...w].reverse().map(x=>1-Number(x)).join(''))/2n;
 assert.deepEqual([0,1].map(e=>BigInt(dp(graph8(10,36,e,w)))),[P(5n,u),P(5n,u)+2n]);
}
const rd=(v)=>[798n*v-29n*q(v),1471n*v-29n*q(v),24384n*v-651n*q(v)];
function read([U,V,N]){
 const v=33457n*(V-U)-120n*(29n*N-651n*U);
 return [v,-202n*(798n*v-U)+9n*(24384n*v-N)];
}
const mod=(x,m)=>(x%m+m)%m;
for(let v=1n;v<=2048n;v++){
 assert.deepEqual(read(rd(v)),[v,q(v)]);
 for(let p of [2n,3n,5n,11n,29n,651n,673n,1000n])
  assert.deepEqual(read(rd(v).map(x=>mod(x,p))).map(x=>mod(x,p)),[v%p,q(v)%p]);
}
console.log('eight-site exact arc pairs and three fixed linear channels replayed');
