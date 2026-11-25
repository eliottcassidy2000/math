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
function poly(v){return 4n*v*v+(9n*v-(v&1n))/2n}
function channels(v,q){
 return [831n*v-29n*q,798n*v-29n*q,1471n*v-29n*q,poly(v)-2n*q];
}
function decode([X,Y,Z,W]){
 let v=28n*(Z-X)-543n*(X-Y);
 return [v,831n*v-X-14n*(poly(v)-W)]
}
function mod(x,m){return (x%m+m)%m}
for(let v=1n;v<=4096n;v++){
 const c=channels(v,q(v));
 assert.deepEqual(decode(c),[v,q(v)]);
 const v0=Number(v),Hc=2*(2*v0)*(2*v0+1)+h(2*v0)+h(2*v0+1)-1+2*v0;
 assert.equal(2n*c[3]+1n,BigInt(Hc));
 for(let m of [2n,3n,5n,11n,29n,673n,100n]){
  const [u,k]=decode(c.map(a=>mod(a,2n*m)));
  assert.equal(mod(u,2n*m),v%(2n*m));assert.equal(mod(k,m),q(v)%m)
 }
}
function pointword(w){
 const n=w.length,A=Array(n+1).fill(0);
 for(let a=0;a<n;a++)for(let b=a+1;b<n;b++)A[a]|=1<<b;
 for(let a=0;a<n;a++)if(w[a]=='1')A[n]|=1<<a;else A[a]|=1<<n;
 return A;
}
const first=[[0,1],[1,1],[0,2]].map(([i,K])=>[0,1].map(e=>dp(graph(10,i,e,K))));
assert.deepEqual(first,[[8679,8681],[8415,8417],[15395,15397]]);
assert.deepEqual(['101000010','101010010'].map(w=>dp(pointword(w))),[241,243]);
function inv(a,p){for(let i=1n;i<p;i++)if(i*a%p===1n)return i;throw Error('nonunit')}
for(let omit=0;omit<4;omit++){
 let p=[673n,5n,3n,29n][omit],m=2n*p,left,right;
 if(omit===3){
   left=[4n,0n];right=[4n,2n];
 }else{
   let a=[798n,831n,831n][omit],seen=new Map();
   for(let r=0n;r<p;r++){
     const v=4n*r; // fixes v modulo 4 as well as surfing all residues mod p
     const q=2n*mod(inv(58n%p,p)*a*v,p); // even q, and a*v=29*q mod p
     const k=channels(v,q).filter((x,i)=>i!==omit).map(x=>String(mod(x,m))).join(',');
     if(seen.has(k)){left=seen.get(k);right=[v,q];break}
     seen.set(k,[v,q]);
   }
 }
 assert(left&&right);
 assert(left.some((x,i)=>mod(x,i===0?m:p)!==mod(right[i],i===0?m:p)));
 assert.deepEqual(channels(...left).filter((x,i)=>i!==omit).map(x=>mod(x,m)),
                  channels(...right).filter((x,i)=>i!==omit).map(x=>mod(x,m)));
 console.log('omit',omit,'modulus',String(p),'collision states',left.map(String),right.map(String))
}
console.log('four full subset-DP replays and universal decoder identity passed')
