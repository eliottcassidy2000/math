// Node: exact check of the prefix-amplification inequality on two shared-cycle factors.
function need(x, k){if(!x) throw Error(JSON.stringify(k));}
function make(p,q,mask){
 const n=p+q+1, r=Array(n).fill(0);
 const A=[0,...Array.from({length:p},(_,i)=>1+i)],
       B=[0,...Array.from({length:q},(_,i)=>1+p+i)];
 for(const L of [A,B]) for(let j=0;j<L.length;j++) for(let k=j+1;k<L.length;k++){
  const [x,y]=j===0&&k===L.length-1?[L[k],L[j]]:[L[j],L[k]];r[x]|=1<<y;
 }
 for(let i=1;i<=p;i++)for(let j=1;j<=q;j++){
  if(mask>>((i-1)*q+j-1)&1)r[A[i]]|=1<<B[j]; else r[B[j]]|=1<<A[i];
 }
 return [r,A,B];
}
function h(r){
 const n=r.length,N=(1<<n)-1,dp=Array.from({length:N+1},()=>Array(n).fill(0));
 for(let i=0;i<n;i++)dp[1<<i][i]=1;
 for(let S=1;S<N;S++)for(let i=0;i<n;i++)if(dp[S][i]){
  let m=r[i]&~S; while(m){const b=m&-m;m^=b;dp[S|b][31-Math.clz32(b)]+=dp[S][i];}
 }
 return dp[N].reduce((a,b)=>a+b,0);
}
let mat=0,local=0,tight=0;
for(let p=2;p<=3;p++)for(let q=2;q<=4;q++)for(let M=0;M<2**(p*q);M++){
 const [r,A,B]=make(p,q,M), e=(x,y)=>r[x]>>y&1;
 const S=h(r)-(2**(p-1)+1)*(2**(q-1)+1);
 for(let i=1;i<=p;i++)for(let j=2;j<=q;j++){
  if(!e(A[i],B[1]) || !e(B[j],A[i]) || i>1&&!e(A[i-1],B[1]))continue;
  const L=1+(j<q&&e(A[i],B[j+1])?q-j:0)+(i<p&&e(B[j],A[i+1])?p-i:0);
  need(S>=L,{p,q,M,i,j,S,L});
  local++;tight+=S===L;
 }
 mat++;
}
console.log({mat,local,tight});
