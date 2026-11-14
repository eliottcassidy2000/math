// node select_family.js TARGET
// Exact acceptance on the listed charts, never a scalar nonexistence test.
"use strict";
const Q5=[28n,25n,10n,16n,4n], Q6=[60n,9n,42n,16n,6n,26n], Q3=[2n,4n,1n];
const ones=n=>(1n<<BigInt(n))-1n;
const power=(a,k)=>a**BigInt(k);
const log=n=>n.toString(2).length-1;
const pow2=n=>n>=2n&&(n&(n-1n))===0n ? log(n) : 0;
function TT(n){return Array.from({length:n},(_,i)=>ones(n)^ones(i+1))}
function blow(Q, blocks){
 const starts=[0];for(const b of blocks)starts.push(starts.at(-1)+b.length);
 return blocks.flatMap((B,i)=>{
  let common=0n;
  for(let j=0;j<Q.length;j++)if(Q[i]>>BigInt(j)&1n)common|=ones(blocks[j].length)<<BigInt(starts[j]);
  return B.map(r=>common|(r<<BigInt(starts[i])));
 });
}
function select(N){
 if(N===1n)return {chart:"singleton",widths:[1],rows:[0n]};
 if(N<3n||!(N&1n))return null;
 let widths;
 const emit=(chart,Q,V)=>({chart,widths:V,rows:blow(Q,V.map(TT))});
 for(let c=1,M=2n;M<N;M*=2n,c++){
  const D=N+M-2n;
  for(let a=1,U=3n;U<D;U=2n*U-1n,a++){
   if(D%U)continue;
   const V=D/U-1n;
   if(V%(M-1n))continue;
   const b=pow2(V/(M-1n));
   if(b)return emit("prime quotient, three widths",Q5,[a,1,1,b,c]);
  }
 }
 const D=4n*N-11n;
 for(let a=1,U=5n;U<D;U=2n*U+3n,a++){
  if(D%U)continue;
  const b=pow2(D/U+3n)-2;
  if(b>=1)return emit("prime quotient, two widths",Q5,[1,a,1,1,b]);
 }
 if((N+15n)%56n===0n){
  const n=pow2((N+15n)/56n);
  if(n)return emit("six-vertex quotient",Q6,[1,2,1,n,1,1]);
 }
 for(let n=1;n<=log(N);n++){
  if(4n*power(3n,n)+5n*power(2n,n)-11n===N)return emit("three-slot cyclic",Q3,[2,2,n]);
 }
 if(N===159n)return {chart:"three-slot cyclic with triangle",widths:[2,2,3],rows:blow(Q3,[TT(2),TT(2),Q3])};
 for(let a=0;a<log(N);a++){
  const six=6n*power(2n,a),three=2n*power(3n,a),p=power(2n,a);
  const offset=p+(a===0?1n:0n)-1n;
  let errors=[[0,0,six+three-2n+offset-2n],[1,0,six+three-5n]];
  const A=six+three-2n;
  if(N>offset&&(N-offset)%A===0n){
   const b=pow2((N-offset)/A);if(b)errors.push([0,b,N]);
  }
  const B=six+three-3n;
  if((N+1n)%B===0n){
   const b=pow2((N+1n)/B);if(b)errors.push([1,b,N]);
  }
  for(const [i,b,total] of errors)if(total===N){
   widths=[1,2,a+b+1];const rows=blow(Q3,widths.map(TT)),u=1+i,v=3+a;
   rows[u]^=1n<<BigInt(v);rows[v]^=1n<<BigInt(u);
   return {chart:"single reversed cyclic seam",widths,edge:[u,v],rows};
  }
 }
 return null;
}
function H(rows){
 const n=rows.length,F=(1<<n)-1;
 const x=rows.map(Number),dp=Array.from({length:F+1},()=>Array(n).fill(0));
 for(let i=0;i<n;i++)dp[1<<i][i]=1;
 for(let s=1;s<=F;s++)for(let j=0;j<n;j++)if((s>>j)&1)
  for(let i=0;i<n;i++)if((x[i]>>j)&1)dp[s][j]+=dp[s^(1<<j)][i];
 return dp[F].reduce((a,b)=>a+b,0);
}
for(const N of [1,5,11,13,45,59,63,83,97,119,137,159,881]){
 const S=select(BigInt(N));if(!S||H(S.rows)!==N)throw Error("target "+N);
}
if(typeof process!=="undefined"&&process.argv[2]){
 const out=select(BigInt(process.argv[2]));
 console.log(JSON.stringify(out,(_,v)=>typeof v==="bigint"?v.toString():v,2));
 if(!out)console.log("No witness in the tested charts; this is not a scalar exclusion.");
}else console.log("PASS: portfolio returns literal witnesses for thirteen exact targets");
