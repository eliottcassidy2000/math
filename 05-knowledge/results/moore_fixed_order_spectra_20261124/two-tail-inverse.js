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

function q(v){let o=0n;while(v){const u=v>>1n;if(v&1n)o+=(u+1n)>>1n;v=u}return o}
function H(v,u){return 1280n*v*u+382n*v+316n*u-58n*(q(v)+q(u))+111n}
function nilgraph(B,i,e,tail){
 const z=word(B),g=charts[i],n=z.length+7+tail.length,m=n-1,v=Array(n).fill(0),off=z.length;
 for(let a=0;a<m;a++)for(let b=a+1;b<m;b++)v[a]|=1<<b;
 for(let a=0,t=0;a<6;a++)for(let b=a+1;b<6;b++,t++)if(g.m[e]>>t&1){
  v[a+off]^=1<<(b+off);v[b+off]^=1<<(a+off);
 }
 for(let a=0;a<m;a++){
  const bit=a<off?+z[a]:a<off+6?g.out>>(a-off)&1:+tail[a-off-6];
  if(bit)v[m]|=1<<a;else v[a]|=1<<m;
 }
 return v;
}
for(let L=1;L<=5;L++)for(let t=1<<(L-1);t<(1<<L);t++){
 const w=word(t),j=BigInt('0b'+[...w].reverse().map(x=>1-Number(x)).join('')),v=5n,u=j/2n;
 for(let i=0;i<2;i++){
  const a=i?H(u,v):H(v,u);
  assert.deepEqual([0,1].map(e=>BigInt(dp(nilgraph(10,i,e,w)))),[a,a+2n])
 }
}
assert.deepEqual([0,1].map(e=>dp(nilgraph(10,0,e,'100'))),[22053,22055]);
assert.deepEqual([0,1].map(e=>dp(nilgraph(10,1,e,'100'))),[21921,21923]);
function reverse(u,k){let v=0n;for(let i=0;i<k;i++){v=(v<<1n)|(u&1n);u>>=1n}return v}
function isqrt(a) {
 if(a<2n)return a;let x=1n<<BigInt((a.toString(2).length+1)>>1);
 while(true){const y=(x+a/x)>>1n;if(y>=x)return x;x=y}
}
function inverse(target){
 if(target<111n||target%2n===0n)return [];
 const out=new Set();
 for(let s=0n,max=isqrt((target-111n)/1280n);s<=max;s++){
  for(let j=0;j<2;j++){
   const L=1280n*s+(j?316n:382n),D=111n+(j?382n:316n)*s-58n*q(s),T=target-D;
   if(T<0n)continue;
   let lo=(T+L-1n)/L,hi=T/(L-29n);
   if(lo<s)lo=s;
   if(hi<lo)continue;
   const store=t=>{if(D+L*t-58n*q(t)===target)out.add([j?s:t,j?t:s].join(','))};
   if(lo===0n){store(0n);lo=1n}
   for(let k=lo.toString(2).length;k<=hi.toString(2).length;k++){
    const start=1n<<BigInt(k-1),end=(start<<1n)-1n;
    const a=lo>start?lo:start,b=hi<end?hi:end;
    if(a>b)continue;
    let z=0,B=1n;while(start>=(L+29n)*4n*B*B){z++;B*=2n}
    for(let p=a/B;p<=b/B;p++){
     let l=0n,r=B;
     while(l<r){
      const m=(l+r)>>1n,t=p*B+reverse(m,z);
      if(D+L*t-58n*q(t)>target)l=m+1n;else r=m;
     }
     if(l<B){const t=p*B+reverse(l,z);if(t>=a&&t<=b)store(t)}
    }
   }
  }
 }
 return [...out].sort()
}
for(let v=0n;v<30n;v++)for(let u=0n;u<30n;u++)
 assert(inverse(H(v,u)).includes([v,u].join(',')));
console.log('62 full subset DP pairs, 900 exact inverses');
const exps=(process.argv[2]||'7,13,17').split(',').map(x=>BigInt(x));
for(let p of exps){
 const E=(1n<<(p-1n))*((1n<<p)-1n);
 for(let d of [-1n,1n]){
  const t=E+d;console.log('target',String(t),'first endpoint',inverse(t),'second endpoint',inverse(t-2n))
 }
}

function mod(z,m){return (z%m+m)%m}
function inverseMod(a,m){
  let b=m,x=1n,y=0n;
  while(b){const t=a/b;[a,b]=[b,a-t*b];[x,y]=[y,x-t*y]}
  assert.equal(a,1n);return mod(x,m)
}
for(let M=1n;M<=128n;M++){
 let m=M;while(m%2n===0n)m/=2n;while(m%5n===0n)m/=5n;
 const a=m===1n?0n:mod(-32n*inverseMod(640n,m),m),u=1n+a;
 const slope=33n+640n*a;
 const invs=inverseMod(slope,2n*M);
 for(let v of [5n,41n,250n,1234n]){
  const A=mod((H(v,u)-316n*u+58n*q(u)-111n)/2n,2n*M);
  const B=mod((H(1n,v)-493n)/2n,2n*M);
  const c=mod(4n*v*v+9n*v/2n-2n*q(v),2n*M);
  const v0=mod((A-B)*invs,2n*M);
  assert.equal(v0,mod(v,2n*M));
  assert.equal(mod((798n*v0-B)-14n*(4n*v0*v0+9n*v0/2n-c),M),mod(q(v),M));
 }
}
console.log('adaptive three-channel inverse for M=1..128')
