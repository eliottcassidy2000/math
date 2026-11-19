const a=require("node:assert/strict");
function f(v){
 let z=1n;
 while(v){
  let q=v/2n;
  if(v%2n===0n)z+=q+q%2n;
  v=q;
 }
 return z;
}
function count(b,e){
 b=BigInt(b);e=BigInt(e);
 return b*e+f(b-1n)*f(e-1n)-(b%2n)*(e%2n);
}
function C(v){
 const n=v?v.toString(2).length:0;
 return [...Array.from({length:n},(_,i)=>((1<<n)-(1<<(i+1)))+(((v>>i)&1)<<n)),((1<<n)-1)^v];
}
function join(b,e){
 const A=C(b-1),B=C(e-1),l=A.length,r=B.length;
 const R=[...A,...B.map((_,i)=>(B.reduce((s,x,j)=>s|(((x>>i)&1)<<j),0)<<l)|((1<<l)-1))];
 R[l+r-1]^=1<<(l-1);R[l-1]|=1<<(l+r-1);
 return R;
}
function H(A){
 const n=A.length,L=1<<n;
 const D=Array.from({length:L},()=>new Float64Array(n));
 for(let i=0;i<n;i++)D[1<<i][i]=1;
 for(let s=1;s<L;s++)for(let i=0;i<n;i++)if(D[s][i]){
  for(let w=A[i]&~s;w;w&=w-1){
   const d=w&-w;
   D[s|d][31-Math.clz32(d)]+=D[s][i];
  }
 }
 return D[L-1].reduce((x,y)=>x+y,0);
}
function prime(n){
 for(let d=2n;d*d<=n;d++)if(n%d===0n)return false;
 return n>1n;
}
const F=[
 [[3,17],[3,206],[3,167]],
 [[35,123],[3,11826],[37,1118]],
 [[3,84871],[1,2283890],[3,666715]],
 [[3,5713007],[1,105909313],[1,108439665]]
];
let p=101n,step=6336n;
for(let j=0;j<F.length;j++){
 a.equal(p,(704n*64n**BigInt(j)+3n)/7n);
 let row=[p,8n*p-1n,8n*p+1n];
 for(let i=0;i<3;i++){
  const [b,e]=F[j][i];
  a.equal(count(b,e),row[i]);
  if(j===0)a.equal(BigInt(H(join(b,e))),row[i]);
 }
 a.equal(prime(p),j===0||j===2);
 p+=step;step*=64n;
}
const A=[96,121,203,33,13,20,56,123], L=[0,1,1,1,2,3,4,5,6,7];
a.equal(H(L.map((u,i)=>L.reduce((s,v,j)=>s|((u===v?+(i<j):(A[u]>>v)&1)<<j),0))),707);
console.log("PASS four exact rows, first literal p / 7p / odd sums, prime tests and recurrence");

const Q=[[3,17],null,[35,123],[1,36458],[3,84871],[53,48470],[3,5713007],[1,156879846]];
let q=101n,d=704n;
const [u,v,w]=[2,3,1],L35=[120,113,115,54,32,64,24],
 L23=[...Array(u).fill(0),...Array(v).fill(1),...Array(w).fill(2)];
const R23=L23.map((x,i)=>L23.reduce((s,y,j)=>s|((x===y?+(i<j):+((y-x+3)%3===1))<<j),0));
a.equal(H(L35),35);a.equal(H(R23),23);
for(let i=0;i<Q.length;i++){
 a.equal(q,(704n*8n**BigInt(i)+3n)/7n);
 a.equal(7n*q,88n*2n**BigInt(3+3*i)+3n);
 a.equal(Q[i]?count(...Q[i]):BigInt(H(L35)*H(R23)),q);
 a.equal(q+d,8n*q-3n);
 q+=d;d*=8n;
}
console.log("PASS first eight consecutive short-return quotient witnesses");
