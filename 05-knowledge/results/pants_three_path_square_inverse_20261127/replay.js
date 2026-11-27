
const eq=(a,b)=>JSON.stringify(a)===JSON.stringify(b);
const diag=(n)=>Array.from({length:n},()=>Array(n).fill(0));
function make(n,mask,root=true){
 const A=diag(n);let b=0;
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){A[i][j]=(root&&i===0)||!(mask>>b&1)?1:-1;A[j][i]=-A[i][j];if(!root||i>0)b++}
 return A
}
function mul(A,B){return A.map(r=>r.map((_,j)=>r.reduce((s,a,k)=>s+a*B[k][j],0)))}
function sw(A,d){return A.map((r,i)=>r.map((a,j)=>a*d[i]*d[j]))}
function canon(A){
 const G=mul(A,A),d=Array(A.length).fill(0),roots=[];
 for(let i=0;i<A.length;i++)if(!d[i]){roots.push(i);d[i]=1;const q=[i];for(const u of q)for(let j=0;j<A.length;j++)if(j!==u&&!d[j]&&G[u][j]){d[j]=Math.sign(G[u][j])*d[u];q.push(j)}}
 const M=sw(G,d);
 return [M.map((r,i)=>r.slice(i+1)).flat().join(','),d,roots.length];
}
function H(A){
 const n=A.length,S=1<<n,dp=Array.from({length:S},()=>Array(n).fill(0)),end=S-1;
 for(let i=0;i<n;i++)dp[1<<i][i]=1;
 for(let s=1;s<S;s++)for(let i=0;i<n;i++)if(dp[s][i])for(let j=0;j<n;j++)if(!(s>>j&1)&&A[i][j]===1)dp[s|1<<j][j]+=dp[s][i];
 return dp[end].reduce((a,b)=>a+b)
}
function h(G){
 const n=G.length, m=n>>1;
 const tr=M=>M.reduce((a,r,i)=>a+r[i],0);
 const a=-tr(G)/2;
 if(m===1)return [a];
 const G2=mul(G,G), b=(a*a-tr(G2)/2)/2;
 if(m===2)return [a,b];
 const c=-(tr(mul(G2,G))/2+a*tr(G2)/2-b*a)/3;
 return [a,b,c]
}
function sf(f){const [a,b,c]=f;return f.length===1|| (f.length===2?a*a-4*b!==0:a*a*b*b-4*b*b*b-4*a*a*a*c-27*c*c+18*a*b*c!==0);}

function P(A,k){
 const n=A.length;
 function go(last,left,used){
  if(!left)return 1;
  let total=0;
  for(let j=0;j<n;j++)if(!(used>>j&1)&&(last<0||A[last][j]>0))total+=go(j,left-1,used|1<<j);
  return total
 }
 return go(-1,k,0)
}
function shortFormula(A,k){
 const n=A.length,G=mul(A,A),G2=mul(G,G);
 let total=Array.from({length:k},(_,i)=>n-i).reduce((a,b)=>a*b,1)/2**(k-1);
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)
  total+=k===3?G[i][j]/2:k===4?(n-3)*G[i][j]/2:(G2[i][j]+3*(n*n-6*n+10)*G[i][j])/8;
 return total
}
function census(n){
 const N=1<<((n-1)*(n-2)/2), M=new Map(), agg={};
 for(let k=0;k<N;k++){
  const [key,_,c]=canon(make(n,k));
  if(!M.has(key))M.set(key,[0,c]);
  M.get(key)[0]++;
 }
 for(const [a,c] of M.values()){
  const m=a*2**(c-1);
  agg[m]=(agg[m]||0)+a;
 }
 return agg
}
function raw(n){
 const counts=new Map(),ans={};
 for(let i=0;i<1<<(n*(n-1)/2);i++){
  const k=JSON.stringify(mul(make(n,i,false),make(n,i,false)));
  counts.set(k,(counts.get(k)||0)+1)
 }
 for(const a of counts.values())ans[a]=(ans[a]||0)+a;
 return ans
}
function assert(q){if(!q)throw Error("replay check failed")}
function poly(G){
 const n=G.length,m=n>>1,X=G.map(r=>r.map(BigInt)),tr=[0n],h=[1n];
 for(let j=1,Y=X;j<=m;j++,Y=mulBig(Y,X)){
  const sum=Y.reduce((s,r,i)=>s+r[i],0n);
  assert(sum%2n===0n);tr.push(sum/2n);
 }
 for(let j=1;j<=m;j++){
  let v=0n;for(let k=1;k<=j;k++)v-=h[j-k]*tr[k];
  assert(v%BigInt(j)===0n);h.push(v/BigInt(j));
 }
 return h.map(Number)
}
function mulBig(A,B){return A.map(r=>r.map((_,j)=>r.reduce((s,a,k)=>s+a*B[k][j],0n)))}
function polyMul(a,b){const r=Array(a.length+b.length-1).fill(0);a.forEach((v,i)=>b.forEach((w,j)=>r[i+j]+=v*w));return r}
function pairWitness(n,maskA,dA,maskB,dB){
 const A=sw(make(n,maskA),dA),B=sw(make(n,maskB),dB),G=mul(A,A);
 assert(eq(G,mul(B,B)));
 let hvals=[H(A),H(B)];
 assert(eq(hvals,n===6?[17,13]:[111,109]));
 for(let b=0;b<1<<(n-1);b++){
  const d=Array.from({length:n},(_,i)=>i&&(b>>(i-1)&1)?-1:1);
  for(let k=3;k<=5;k++){
   const X=sw(A,d),Y=sw(B,d);
   assert(P(X,k)===P(Y,k)&&P(X,k)===shortFormula(X,k));
  }
 }
 console.log("same-square witness",n,"path counts (A,B)",Array.from({length:n-2},(_,i)=>[P(A,i+3),P(B,i+3)]),"h",poly(G));
 return G
}
function bank(mask){
 const sz=[3,5,7],index=[];
 sz.forEach((l,i)=>{for(let x=0;x<l;x++)index.push([i,x])});
 const A=diag(15);
 for(let a=0;a<15;a++)for(let b=a+1;b<15;b++){
  const [i,x]=index[a],[j,y]=index[b];
  A[a][b]=(1-2*(mask>>(i===j?i:3)&1))*(i===j?(-1)**(x+y):1);
  A[b][a]=-A[a][b]
 }
 return A
}
const table=[
 {"2":2},{"2":6,"16":2},{"2":24,"4":40},{"2":480,"4":384,"8":40,"24":120},
 {"2":17520,"4":14448,"8":560,"240":240}
];
for(let n=3;n<=7;n++){
 const c=census(n);assert(eq(c,table[n-3]));
 if(n<=6){
  const d=raw(n);for(const m in c)assert(d[m]===c[m]*2**(n-1));
 }
 console.log("square census",n,c)
}
pairWitness(6,41,[1,-1,1,-1,1,1],262,[1,-1,1,1,-1,1]);
const S=pairWitness(7,73,[1,1,-1,-1,-1,1,1],4614,[1,1,-1,-1,-1,1,1]);
assert(eq(poly(S),polyMul([1,7],[1,14,17])));
const A=bank(0),G=mul(A,A),cuts=[],keys=new Set();
for(let b=0;b<24;b++){
 const d=Array.from({length:15},(_,i)=>((b*5371+273*(b%4))>>i&1)?-1:1);
 cuts.push([d,P(sw(A,d),3)]);
}
for(let k=0;k<16;k++){
 const B=bank(k);assert(eq(mul(B,B),G));keys.add(JSON.stringify(B));
 for(const [d,p] of cuts)assert(p===P(sw(B,d),3));
}
assert(keys.size===16);
const expected=[[1,3],[1,10,5],[1,21,35,7],[1,71]].reduce(polyMul);
assert(eq(poly(G),expected));
assert(eq(poly(mul(sw(make(6,0),[1,1,1,1,1,1]),make(6,0))),[1,15,15,1]));
const E=make(6,0,false);E[1][4]=-1;E[4][1]=1;
assert(eq(poly(mul(E,E)),[1,15,47,1]));
const pf=E[1][3]*E[4][5]-E[1][4]*E[3][5]+E[1][5]*E[3][4];
assert(Math.abs(pf)===3);
console.log("15-vertex example:",keys.size,"same-square orientations; reduced polynomial",expected);
console.log("irreducible nontransitive-switching example verified");
