const d=(i,j)=>Math.min(Math.abs(i-j),10-Math.abs(i-j));
const m=x=>(x%5+5)%5;
let Q=Array.from({length:10},(_,i)=>Array.from({length:10},(_,j)=>d(i,j)===1?2:0)),a=Array(10).fill(3),sol=[],nn=0;
function dfs(i,j){
 if(i===10){nn++;if(Q.every((r,i)=>r.every((v,j)=>v+r.reduce((s,x,k)=>s+x*Q[k][j],0)===(i===j?11:5))))sol.push(Q.map(r=>r.slice()));return}
 if(j===10){if(!a[i])dfs(i+1,i+2);return}
 dfs(i,j+1);if(a[i]&&a[j]&&[2,4,5].includes(d(i,j))){Q[i][j]=Q[j][i]=1;a[i]--;a[j]--;dfs(i,j+1);Q[i][j]=Q[j][i]=0;a[i]++;a[j]++}
}
dfs(0,1);
if(nn!==814||sol.length!==2)throw Error('quotient tally');
const edges=Q=>Q.flatMap((r,i)=>r.flatMap((x,j)=>x===1&&i<j?[[i,j]]:[]));
let shift=new Set(edges(sol[0]).map(([i,j])=>[ (i+1)%10,(j+1)%10].sort((a,b)=>a-b).join()));
if(!edges(sol[1]).every(e=>shift.has(e.join())))throw Error('rotation');
console.log('814 projected degree completions, 2 projected Moore matrices related by rotation');
for(const Q of sol){
 const E=edges(Q), n=E.length+1;
 const vec=()=>Array(n).fill(0);
 const C=Q.map(()=>Q.map(vec));
 E.forEach(([i,j],k)=>{C[i][j][k]=1; C[j][i][k]=-1});
 // Last variable is center voltage on the closing ring edge 9 -> 0.
 C[9][0][n-1]=1;C[0][9][n-1]=-1;
 const L=[];
 for(let i=0;i<10;i++)for(let j=i+1;j<10;j++){
  const row=C[i][j].map(x=>x*Q[i][j]);
  for(let k=0;k<10;k++) for(let s=0;s<n;s++)row[s]+=Q[i][k]*Q[k][j]*(C[i][k][s]+C[k][j][s]);
  L.push(row.map(m));
 }
 function space(L,n){
  L=L.map(x=>x.slice());
  let rank=0,P=[];
  for(let k=0;k<n;k++){
   let p=L.findIndex((r,i)=>i>=rank&&r[k]); if(p<0)continue;
   [L[p],L[rank]]=[L[rank],L[p]];p=rank++;
   const s=[0,1,3,2,4][L[p][k]];L[p]=L[p].map(x=>m(x*s));P.push(k);
   for(let i=0;i<L.length;i++)if(i!==p){let t=L[i][k];L[i]=L[i].map((x,j)=>m(x-t*L[p][j]));}
  }
  const out=[],free=Array.from({length:n},(_,i)=>i).filter(i=>!P.includes(i));
  for(let mask=0;mask<5**free.length;mask++){
   let h=mask,b=Array(n).fill(0);
   for(const k of free){b[k]=h%5;h=Math.floor(h/5)}
   P.forEach((k,i)=>b[k]=m(-L[i].reduce((z,x,j)=>z+x*b[j],0)));
   out.push(b);
  }
  return [rank,out];
 }
 for(let parity=0;parity<2;parity++){
  const [rank,V]=space(L,n),good=[];
  for(const v of V){
   const U=Q.map((r,i)=>r.map((x,j)=>{
    if(!x)return [];
    const c=m(C[i][j].reduce((z,y,k)=>z+y*v[k],0));
    if(x===1)return [c];
    const r=1+(( (j===(i+1)%10?i:j) + parity)%2);
    return [m(c-r),m(c+r)];
   }));
   let ok=true;
   for(let i=0;i<10;i++)for(let j=0;j<10;j++){
    let b=Array.from({length:5},(_,k)=>Number(U[i][j].includes(k)));
    for(let k=0;k<10;k++)for(let x of U[i][k])for(let y of U[k][j])b[m(x+y)]++;
    if(b.some((x,k)=>x!==1+6*(i===j&&k===0)))ok=false;
   }
   if(ok)good.push(v);
  }
  console.log('first simple neighbors',E.filter(e=>e[0]===0).map(e=>e[1]).join(), 'color parity',parity,
   'linear rank',rank,'candidates',V.length,'exact',JSON.stringify(good));
  if(good.length!==1||good[0][n-1]!==0)throw Error('lift tally');
 }
}

const perms=[];
function perm(a,k=1){if(k===5){perms.push(a.slice());return}for(let j=k;j<5;j++){[a[k],a[j]]=[a[j],a[k]];perm(a,k+1);[a[k],a[j]]=[a[j],a[k]]}}
perm([0,1,2,3,4]);
let D=[[0,0,0,0,0],[0,1,2,3,4]],A=[];
function matrix(){
 if(D.length===5){A.push(D.slice());return}
 for(const p of perms)if(D.every(q=>new Set(p.map((v,i)=>m(v-q[i]))).size===5)){D.push(p);matrix();D.pop()}
}
matrix();
if(A.length!==6||A.some(M=>M.some((r,i)=>r.some((x,j)=>x!==m(r[1]*j)))))throw Error('cyclic block');
console.log('cyclic block: 6 normalized matrices, all row reorderings of ij');
const moments=[];let x=2,y=49;
for(let j=0;j<=6;j++){
 moments.push(x+28*(-3)**j+21*2**j+1679*7**j+1470*(-8)**j);
 [x,y]=[y,49*y+406*x];
}
if(moments.join()!=='3200,0,179900,0,20048700,541731400,32163036500')throw Error('moments');
console.log('3200-site moments:',moments.join(', '));
