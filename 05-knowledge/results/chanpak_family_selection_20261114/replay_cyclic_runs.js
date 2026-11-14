// node replay_cyclic_runs.js
"use strict";
const check=(c,m)=>{if(!c)throw new Error(m)};
const bit=(x,i)=>(x>>i)&1;
function paths(R){
 const n=R.length,z=1<<n,E=Array.from({length:z},()=>Array(n).fill(0)),H=Array(z).fill(0);H[0]=1;
 for(let i=0;i<n;i++)E[1<<i][i]=1;
 for(let S=1;S<z;S++)for(let j=0;j<n;j++)if(bit(S,j)){
  for(let i=0;i<n;i++)if(bit(R[i],j))E[S][j]+=E[S^(1<<j)][i];
  H[S]+=E[S][j];
 }
 return H;
}
function covers(R){
 const H=paths(R),n=R.length,F=H.length-1,out=[0];
 let P=Array(H.length).fill(0);P[0]=1;
 for(let k=1;k<=n;k++){
  const Q=Array(H.length).fill(0);
  for(let S=1;S<=F;S++)for(let W=S;W;W=(W-1)&S)Q[S]+=H[W]*P[S^W];
  out.push(Q[F]);P=Q;
 }
 return out;
}
function canon(R){
 const n=R.length,Ds=R.map(x=>x.toString(2).split("1").length-1);
 const d=Ds.slice().sort((a,b)=>a-b);let best=Infinity;
 function walk(p,u){
  if(p.length===n){let w=0;
   for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)w=2*w+bit(R[p[i]],p[j]);
   best=Math.min(best,w);return;
  }
  for(let i=0;i<n;i++)if(!bit(u,i)&&Ds[i]===d[p.length])walk([...p,i],u|(1<<i));
 }
 walk([],0);return best;
}
function extend(R,U){
 const n=R.length;return R.map((r,i)=>r|((1-bit(U,i))<<n)).concat(U);
}
function cyclic(A,B,C){
 const blocks=[A,B,C],s=[0,A.length,A.length+B.length];
 return blocks.flatMap((T,i)=>{
  const j=(i+1)%3;
  return T.map(r=>(r<<s[i])|(((1<<blocks[j].length)-1)<<s[j]));
 });
}
function read(A,B,C){
 const t=(r,i)=>r[i]||0;let ans=0;
 for(let k=1;k<Math.max(A.length,B.length,C.length);k++){
  ans+=3*t(A,k)*t(B,k)*t(C,k);
  for(let u=1;u<7;u++)ans+=t(A,k+bit(u,0))*t(B,k+bit(u,1))*t(C,k+bit(u,2));
 }
 return ans;
}
const TT=n=>Array.from({length:n},(_,i)=>((1<<n)-1)^((2<<i)-1));
let L=[[0]],M=[];
for(let n=1;n<=7;n++){
 check(L.length===[0,1,1,2,4,12,56,456][n],"isomorphism count");
 M.push(...L.map(R=>[R,covers(R)]));
 if(n<7){
  const next=new Map();
  for(const R of L)for(let U=0;U<(1<<n);U++){const B=extend(R,U);next.set(canon(B),B);}
  L=[...next.values()];
 }
}
let controls=0;
for(const [A,a] of M)for(const [B,b] of M){
 if(A.length+B.length>=7)continue;
 for(const [C,c] of M)if(A.length+B.length+C.length<=7){
  check(paths(cyclic(A,B,C)).at(-1)===read(a,b,c),"run decomposition");
  controls++;
 }
}
const exclude=[45,63,97,137,159];
for(const N of exclude)for(const [A,a] of M)for(const [B,b] of M){
 const x=a[1],y=b[1],p=a[2]||0,q=b[2]||0;
 check(3*x*y+p*y+x*q+p*q!==N,"singleton exclusion");
}
for(const [B,N] of [[TT(2),45],[TT(3),137],[[4,1,2],159]]){
 check(paths(cyclic(TT(2),B,TT(2))).at(-1)===N,"literal fallback");
}
for(let n=1;n<=9;n++){
 check(read(covers(TT(2)),covers(TT(n)),covers(TT(2)))===4*3**n+5*2**n-11,"unbounded formula");
}

const split=cyclic(TT(2),TT(2),TT(1));
split[0]^=1<<2;split[2]^=1<<0;
function templates(j) {
 const ans=[0];
 function visit(previous,used,left){
  if(used===15&&left===0){ans[0]++;return;}
  for(let v=0;v<=4;v++){
   if(v===4?left===0:(used>>v)&1)continue;
   if(previous>=0&&!bit(split[previous],v))continue;
   if(v===4)visit(v,used,left-1);else visit(v,used|(1<<v),left);
  }
 }
 visit(-1,0,j);return ans[0];
}
check([1,2,3,4].map(templates).join()==="13,12,2,0","split skeleton coefficient");
for(let n=1;n<=7;n++){
 const a=cyclic(TT(2),TT(2),TT(n));a[0]^=1<<2;a[2]^=1<<0;
 check(paths(a).at(-1)===2*3**n+6*2**n-5,"split seam");
}
const z=cyclic(TT(1),TT(1),TT(6));
for(const j of [4,5]){z[0]^=1<<j;z[j]^=1;}
check(z.join()==="50,252,249,241,224,192,129,1"&&paths(z).at(-1)===63,"two-contact repair at 63");
console.log("PASS: "+controls+" literal uniform controls; targets 45,63,97,137,159; and all-order singleton exclusions for those five");

