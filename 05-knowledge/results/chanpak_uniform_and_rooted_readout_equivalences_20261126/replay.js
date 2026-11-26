const assert = require('node:assert/strict');
function graph(n,m,cycle=false){
 let a=Array(n).fill(0),k=0;
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){
  let b=cycle&&i===0&&j===n-1;
  if(!(cycle&&(b||j===i+1)))b=!!(m>>k++&1);
  if(b)a[j]|=1<<i;else a[i]|=1<<j;
 }return a;
}
function inflate(a,sizes){
 const I=sizes.flatMap((m,i)=>Array.from({length:m},(_,j)=>[i,j]));
 return I.map(([i,j])=>I.reduce((h,[k,l],v)=>h|((i===k?+(j<l):(a[i]>>k&1))<<v),0));
}
function paths(a){
 const n=a.length,u=(1<<n)-1,h=new Float64Array(u+1),d=new Float64Array((u+1)*n);h[0]=1;
 for(let s=1;s<=u;s++)for(let b=s;b;b&=b-1){
  const v=31-Math.clz32(b&-b),t=s^(b&-b);let x=t?0:1;
  for(let c=t;c;c&=c-1){let j=31-Math.clz32(c&-c);if(a[j]>>v&1)x+=d[t*n+j];}
  h[s]+=d[s*n+v]=x;
 }return h;
}
function triangleHost(a,k){
 const n=a.length, b=[...a,...Array(2*k).fill(0)];
 for(let i=0;i<b.length;i++)for(let j=0;j<b.length;j++){
  const t=Math.min(2,Math.max(0,1+Math.floor((i-n)/k))),u=Math.min(2,Math.max(0,1+Math.floor((j-n)/k)));
  if((t+1)%3===u||t===u&&t>0&&i<j)b[i]|=1<<j;
 }return b;
}
function chainCovers(k,j){
 if(j<0||j>k)return 0n;
 let dp=Array(k+1).fill(0n);dp[0]=1n;
 for(let r=1;r<=k;r++)for(let i=r;i>=0;i--)dp[i]=BigInt(i)*(dp[i]+(i?dp[i-1]:0n));
 return dp[j];
}
function coeff(k,j){
 const p=[[0,0],[0,1],[1,1]],s=[[0,0],[1,0],[1,1]];let c=0n;
 for(const [a,b] of p)for(const [x,y] of s)c+=chainCovers(k,j-1+a+x)*chainCovers(k,j-1+b+y);
 return c;
}
function fact(n){let x=1n;for(let i=1;i<=n;i++)x*=BigInt(i);return x;}
function radix(n, digits){
 const N=BigInt(n), M=fact(n)*(2n**BigInt(n-1));let k=1,a= N+1n,b=N;
 while(a<2n*(N+1n)*b||a*a<=36n*N*M*b*b){k++;a*=N+1n;b*=N;}
 const c=Array.from({length:n},(_,i)=>coeff(k,i+1));
 let h=c.reduce((s,x,i)=>s+x*digits[i],0n),v=h;
 for(let j=n-1;j>=0;j--){assert.equal(v/c[j],digits[j]);v%=c[j];}
 return [k,h.toString().length];
}
for(const [n,m] of [[5,0],[5,1],[6,0],[6,11]]){
 const a=graph(n,m,true),d=Array(n+1).fill(0n),h=paths(a);
 let table=[new Float64Array(h.length)];table[0][0]=1;
 for(let j=1;j<=n;j++){const t=new Float64Array(h.length);
  for(let s=1;s<h.length;s++)for(let r=s;r;r=(r-1)&s)t[s]+=h[r]*table[j-1][s^r];
  table.push(t);d[j]=BigInt(t.at(-1));}
 const r=[Number(d[1])];
 for(let k=1;k<=3;k++){const x=paths(triangleHost(a,k)).at(-1);
  assert.equal(BigInt(x),d.slice(1).reduce((s,y,i)=>s+y*coeff(k,i+1),0n));r.push(x);}
 const tri=paths(inflate(a,[2,2,2,...Array(n-3).fill(1)])).at(-1);
 console.log(JSON.stringify({n,m,r,tri,radix:radix(n,d.slice(1))}));
}
for(let n=3;n<=5;n++){
 const triples=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)for(let l=j+1;l<n;l++)triples.push([i,j,l]);
 const hand=new Map(),rev=new Map();let lo=Infinity,hi=0;
 for(let m=0;m<(1<<(n*(n-1)/2));m++){
  const a=graph(n,m),t=[],p=[];
  for(const s of triples){const cic=(s.every((i)=>s.filter(j=>a[i]>>j&1).length===1));t.push(+cic);
   const h=paths(inflate(a,a.map((_,i)=>1+Number(s.includes(i))))).at(-1);p.push(h);
   if(cic)lo=Math.min(lo,h);else hi=Math.max(hi,h);
  }const x=t.join(''),y=p.join(',');
  if(hand.has(x))assert.equal(hand.get(x),y);
  hand.set(x,y);if(rev.has(y))assert.equal(rev.get(y),x);rev.set(y,x);
 }
 console.log(JSON.stringify({n,hypergraphs:hand.size,minCyclic:lo,maxTransitive:hi}));
}
