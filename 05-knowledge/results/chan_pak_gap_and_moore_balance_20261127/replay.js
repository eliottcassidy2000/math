
function mat(n,code){
 const A=Array.from({length:n},()=>Array(n).fill(0));let k=0;
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){
  A[i][j]=code>>k&1?-1:1;A[j][i]=-A[i][j];k++;
 }
 return A
}
function square(A){
 return A.map((r,i)=>r.map((v,j)=>r.reduce((s,w,k)=>s+w*A[k][j],0)))
}
function H(A){
 const n=A.length,all=(1<<n)-1,dp=new Int32Array((all+1)*n),out=A.map(r=>r.reduce((s,x,i)=>x>0?s|1<<i:s,0));
 for(let i=0;i<n;i++)dp[(1<<i)*n+i]=1;
 for(let m=1;m<=all;m++)for(let i=0;i<n;i++)if(dp[m*n+i]){
  let d=out[i]&~m;
  while(d){let x=d&-d;d^=x;dp[(m|x)*n+(31-Math.clz32(x))]+=dp[m*n+i]}
 }
 return dp.slice(-n).reduce((a,b)=>a+b)
}
function strong(A){
 const n=A.length,s=A.map(r=>r.reduce((a,b)=>a+(b===1),0)).sort((a,b)=>a-b);
 let sum=0;
 for(let i=0;i<n-1;i++){sum+=s[i];if(sum===(i+1)*i/2)return false;}
 return true
}
function C3(A){let t=0;for(let i=0;i<A.length;i++)for(let j=i+1;j<A.length;j++)for(let k=j+1;k<A.length;k++)if(A[i][j]===A[j][k]&&A[j][k]===A[k][i])t++;return t}
function scan(n){
 const K=1<<(n*(n-1)/2),M=new Map(),stat={},st={};let haveCross=null,triangleCross=null,scoreCross=null;
 for(let mask=0;mask<K;mask++){
  const A=mat(n,mask),cnt=H(A),key=JSON.stringify(square(A)),isStrong=strong(A),t=C3(A);
  stat[cnt]=(stat[cnt]||0)+1;if(isStrong)st[cnt]=(st[cnt]||0)+1;
  if(!M.has(key))M.set(key,{});
  M.get(key)[cnt]=mask;
 }
 const patterns={};
 for(const x of M.values()){
  const k=Object.keys(x).join(',');
  patterns[k]=(patterns[k]||0)+1;
  if(x[19]!==undefined&&x[23]!==undefined)haveCross=[x[19],x[23]];
 }
 return {stat,st,haveCross,patterns}
}

const assert=(x)=>{if(!x)throw Error("check failed")};
const list=o=>Object.keys(o).map(Number);
const rep=[3,4,5,6].map(scan);
assert(rep.every(x=>x.stat[7]===undefined&&x.stat[21]===undefined));
for(let n=3;n<=6;n++){
 const r=rep[n-3];console.log("order",n,"strong spectrum",list(r.st));
}
assert(Object.keys(rep[3].patterns).filter(s=>s.includes(',')).join('')==="13,17");
console.log("six-vertex exact-square fibers: the only nonconstant Hamiltonian spectrum is {13,17}");
function cycles(A,k){
 const n=A.length;
 function go(first,last,used,left){
  if(!left)return +(A[last][first]>0);
  let x=0;
  for(let j=first+1;j<n;j++)if(!(used>>j&1)&&A[last][j]>0)x+=go(first,j,used|1<<j,left-1);
  return x
 }
 return A.reduce((sum,_,i)=>sum+go(i,i,1<<i,k-1),0)
}
function d33(A){
 let r=0;
 if(A.length===6)for(let m=1;m<64;m+=2)
  if(m.toString(2).split("1").length===4){
   const X=A.filter((_,i)=>m>>i&1).map(row=>row.filter((_,i)=>m>>i&1)),
         Y=A.filter((_,i)=>!(m>>i&1)).map(row=>row.filter((_,i)=>!(m>>i&1)));
   r+=cycles(X,3)*cycles(Y,3);
  }
 return r
}
const c5=A=>cycles(A,5);
const pack=A=>[C3(A),c5(A),d33(A)];
const deg=A=>A.map(r=>r.filter(x=>x>0).length);
const eq=(a,b)=>JSON.stringify(a)===JSON.stringify(b);
for(const [n,m,t,cn,p] of [[5,4,260,[5,9],[[2,0,0],[3,1,0]]],[6,20,2054,[19,23],[[5,4,0],[5,4,1]]]]){
 const A=mat(n,m),B=mat(n,t);
 assert(eq([H(A),H(B)],cn)&&eq([pack(A),pack(B)],p));
 if(n===6)assert(eq(deg(A),deg(B)));
 console.log("gap crossing",cn,"(cycles 3,5; disjoint 3+3)",p,"core strong",strong(A),strong(B));
}
function orderedJoin(A,B){
 const a=A.length,b=B.length;
 return A.map(r=>r.concat(Array(b).fill(1))).concat(B.map(r=>Array(a).fill(-1).concat(r)));
}
const tri=mat(3,2),X=orderedJoin(tri,mat(4,0)),Y=orderedJoin(tri,tri);
assert(H(X)===3&&H(Y)===9);
function blow(Q,M,r){
 const ix=Q.flatMap((_,i)=>i===r?M.map((_,j)=>[i,j]):[[i,0]]);
 return ix.map(([a,i])=>ix.map(([b,j])=>a===b?(a===r?M[i][j]:0):Q[a][b]))
}
const b4=new Set();
for(let mask=0;mask<64;mask++)if(strong(mat(4,mask)))for(let r=0;r<4;r++)b4.add(H(blow(mat(4,mask),tri,r)));
assert(H(blow(tri,tri,0))===15&&eq([...b4].sort(),[27,33]));
console.log("small strong tournament with proper regular module: order 5 -> 15; order 6 -> 27 or 33");
const pairs5=[];
for(let i=0;i<5;i++)for(let j=i+1;j<5;j++)pairs5.push([i,j]);
const pet=pairs5.map(x=>pairs5.map(y=>+x.every(v=>!y.includes(v))));
const hs=Array.from({length:50},(_,a)=>Array.from({length:50},(_,b)=>{
 const ar=(a%25)/5|0,br=(b%25)/5|0,x=a%5,y=b%5;
 if(a>=25&&b<25)return +((x-y-ar*br+25)%5===0);
 if(a<25&&b>=25)return +((y-x-ar*br+25)%5===0);
 if(ar!==br)return 0;
 const s=(x-y+5)%5,k=a<25?1:2;
 return +(s===k||s===5-k);
}));
function paley(p){
 const R=new Set(Array.from({length:p-1},(_,i)=>(i+1)**2%p));
 return Array.from({length:p},(_,i)=>Array.from({length:p},(_,j)=>i===j?0:R.has((j-i+p)%p)?1:-1));
}
function avg(G,A,i,j){
 let s=0;
 for(let a of A[i])for(let b of A[j])s+=G[a][b];
 return s;
}
function moore(B,d,Q){
 const n=B.length,m=d-1,X=square(B);
 assert(n===d*d+1);
 for(let i=0;i<n;i++)for(let j=0;j<n;j++)assert(X[i][j]+B[i][j]===(i===j?d:1));
 const lead=B[0].map((x,i)=>x?i:-1).filter(i=>i>=0);
 const parts=lead.map(i=>B[i].map((x,j)=>x&&j!==0?j:-1).filter(j=>j>=0));
 assert(parts.every(p=>p.length===m));
 const N=parts.flat(),A=N.map((a,i)=>N.map((b,j)=>{
  const x=i/m|0,y=j/m|0;
  return x===y?+(i<j)-+(i>j):Q[x][y]*(1-2*B[a][b])
 }));
 for(let i=0;i<d;i++)for(let j=0;j<d;j++)if(i!==j)
  for(let a of parts[i])for(let b of parts[j]){
   let count=0;
   for(let k=0;k<d;k++)if(k!==i&&k!==j)for(let x of parts[k])count+=B[a][x]*B[x][b];
   assert(count+B[a][b]===1)
  }
 const G=square(A),R=square(Q),indices=parts.map((_,i)=>Array.from({length:m},(_,j)=>i*m+j));
 let max=0;
 for(let i=0;i<d;i++)for(let j=i+1;j<d;j++){
  const s=avg(G,indices,i,j);
  assert(s===m*(m-2)**2*R[i][j]);max=Math.max(max,Math.abs(s)/m**2)
 }
 console.log("degree",d,"root transport verified; max absolute cross-block mean",max);
}
moore(pet,3,paley(3));moore(hs,7,paley(7));
const S=paley(59).slice(2).map(r=>r.slice(2)),G=square(S),h={};
for(let i=0;i<57;i++)for(let j=i+1;j<57;j++)h[G[i][j]]=(h[G[i][j]]||0)+1;
assert(h[3]===378&&Math.max(...Object.keys(h).map(n=>Math.abs(+n)))===3);
console.log("57-vertex quotient from Paley 59: unordered Gram entries",h,
 "optimal max cross-block mean",3*54**2/56);

function diamonds(A){
 const n=A.length;let d=0;
 for(let a=0;a<n;a++)for(let b=a+1;b<n;b++)for(let c=b+1;c<n;c++)for(let e=c+1;e<n;e++){
  const p=A[a][b]*A[c][e]-A[a][c]*A[b][e]+A[a][e]*A[b][c];
  if(Math.abs(p)===3)d++;
 }
 return d
}
for(const [n,x,y,a,b] of [[5,4,260,2,2],[6,20,2054,3,5]]){
 const A=mat(n,x),B=mat(n,y);
 assert(diamonds(A)===a&&diamonds(B)===b);
 console.log("diamond counts across gap",n===5?7:21,a,b);
}
