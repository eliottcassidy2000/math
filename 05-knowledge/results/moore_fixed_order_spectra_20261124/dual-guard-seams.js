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
for(const [i,B,H] of [[0,18,15269],[0,82,67931],[1,74,58907],[1,664,528509]]){
 const X=[graph(B,i,0),graph(B,i,1)];
 assert.equal(X[0].filter((x,j)=>x!==X[1][j]).length,2);
 assert.deepEqual(X.map(dp),[H,H+2]);assert.equal(formula(B,i),H);
 assert(prime(H)&&prime(H+2));
 // Removing the only guard collapses the seam count.
 assert.equal(dp(graph(B,i,0,0)),dp(graph(B,i,1,0)));
 const v=B/2, p=B===82||B===664?'1010':'', vmax=Math.floor((H-2*charts[i].R-111)/(2*(640+charts[i].J)));
 for(let u=0;u<=Math.max(v-1,vmax);u++){
   if(u===v||!word(u).startsWith(p))continue;
   const F=formula(2*u,i);assert(!((u<v||F<=H)&&prime(F)&&prime(F+2)));
 }
 console.log('full DP seam',i,B,word(B).length+9,H,H+2)
}
for(const [v,w,H] of [[1274,1339,2104127],[5180,5463,8574323]]){
 assert.equal(formula(v*2,0),H);assert.equal(formula(w*2,1),H+2);
}
console.log('scalar three-level overlaps checked')
