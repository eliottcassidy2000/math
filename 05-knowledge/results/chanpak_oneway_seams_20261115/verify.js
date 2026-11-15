
function chain(ms){let t=0,revs=[];for(const m of ms){revs.push([t,t+m-1]);t+=m-1;}const R=Array(t+1).fill(0);for(let i=0;i<=t;i++)for(let j=i+1;j<=t;j++)if(revs.some(z=>z[0]===i&&z[1]===j))R[j]|=1<<i;else R[i]|=1<<j;return R;}
function data(R){const n=R.length,N=(1<<n)-1,S=Array.from({length:N+1},()=>Array(n).fill(0)),E=S.map(z=>z.slice());for(let i=0;i<n;i++)S[1<<i][i]=E[1<<i][i]=1;for(let m=1;m<=N;m++)for(let i=0;i<n;i++)if(m>>i&1&&m!==(1<<i)){const o=m^(1<<i);for(let j=0;j<n;j++)if(o>>j&1){if(R[i]>>j&1)S[m][i]+=S[o][j];else E[m][i]+=E[o][j];}}
function add(z,r){return z.reduce((s,k,i)=>s+(r>>i&1?k:0),0);}
return R.map((out,v)=>{const m=N^(1<<v), a0=add(E[m],out);let a1=0;for(let s=(m-1)&m;s;s=(s-1)&m)a1+=add(E[s],out)*add(S[m^s],out);return [a0,a1,S[N][v],E[N][v]]});}

function join(a,b,i,j) {
 const A=Array.from({length:a.length},(_,u)=>u), B=Array.from({length:b.length},(_,u)=>u===j?i:a.length+u-(u>j?1:0));
 const R=Array(a.length+b.length-1).fill(0);
 for(const [rr,L] of [[a,A],[b,B]])for(let x=0;x<L.length;x++)for(let y=0;y<L.length;y++)if(rr[x]>>y&1)R[L[x]]|=1<<L[y];
 for(const u of A)for(const v of B)if(u!==i&&v!==i)R[u]|=1<<v;
 return R;
}
function h(R){const n=R.length,N=(1<<n)-1,dp=Array.from({length:N+1},()=>Array(n).fill(0));
 for(let i=0;i<n;i++)dp[1<<i][i]=1;
 for(let S=1;S<N;S++)for(let i=0;i<n;i++)if(dp[S][i]){
  let m=R[i]&~S;while(m){const b=m&-m;m^=b;dp[S|b][31-Math.clz32(b)]+=dp[S][i];}
 }
 return dp[N].reduce((x,y)=>x+y,0);}
function need(t,z){if(!t)throw Error(JSON.stringify(z));}
function conf(ms){const e=ms.reduce((a,b)=>a+(b===3?1:b===5?2:0),0);return e<=2;}
function words(ms){
 let ans=ms.length&&conf(ms)?[ms]:[];
 if(ms.length<4)for(const b of [3,4,5])if(ms.reduce((a,b)=>a+b-1,1)+b-1<=13&&conf([...ms,b]))ans=ans.concat(words([...ms,b]));
 return ans;
}
function op(R){return R.map((_,i)=>R.reduce((v,a,j)=>v+((a>>i&1)<<j),0));}
const W=words([]);let profiles=0,joins=0,small=0;
for(const ms of W){const r=chain(ms),d=data(r),h0=h(r);
 for(let i=0;i<r.length;i++){const [a,b,s,e]=d[i],k=a+b;
  let target=null;
  if(ms.length===1)target=({3:[0,0,0],4:[2,4,1,0],5:[6,16,7,1,0]})[ms[0]][i];
  else if(ms.length===2&&ms[0]===3&&ms[1]===3)target=[12,14,2,2,0][i];
  else if(i===r.length-1)target=0;
  need(target===null?k>=3:k===target,[ms,i,k,target]);
  need(9*Math.max(s,e)<=4*h0,[ms,i,s,e,h0]);
  profiles++;small+=k<=2;
 }
}
for(const x of W.filter(m=>m.length<=2))for(const y of W.filter(m=>m.length<=2)){
 const A=chain(x),B=chain(y);if(A.length+B.length>12)continue;
 const da=data(A),db=data(op(B)),ha=h(A),hb=h(B);
 for(let i=0;i<A.length;i++)for(let j=0;j<B.length;j++){
  const [a0,ap,,ea]=da[i], [b0,bp,,sb]=db[j];
  const L=(a0+ap)*(hb-sb)+(b0+bp)*(ha-ea)+ap*bp;
  const real=h(join(A,B,i,j))-ha*hb;
  need(L===real,[x,y,i,j,L,real]);joins++;
 }
}
console.log({chains:W.length,profiles,small,joins});
