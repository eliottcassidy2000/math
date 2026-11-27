
"use strict";
const assert = p=>{if(!p)throw Error("check failed")};
const sq=A=>A.map((r,i)=>r.map((_,j)=>r.reduce((s,v,k)=>s+v*A[k][j],0)));
const eq=(a,b)=>JSON.stringify(a)===JSON.stringify(b);
const range=n=>Array.from({length:n},(_,i)=>i);
const pairs=[];
for(let i=0;i<5;i++)for(let j=i+1;j<5;j++)pairs.push([i,j]);
const pet=pairs.map(a=>pairs.map(b=>+a.every(v=>!b.includes(v))));
const hs=range(50).map(a=>range(50).map(b=>{
 const i=((a%25)/5)|0,j=((b%25)/5)|0,x=a%5,y=b%5;
 if(a<25&&b>=25)return +((y-x-i*j+25)%5===0);
 if(a>=25&&b<25)return +((x-y-i*j+25)%5===0);
 const k=a<25?1:2,t=(x-y+5)%5;
 return +(i===j&&(t===k||t===5-k));
}));
function paley(n){
 const r=new Set(range(n-1).map(i=>(i+1)**2%n));
 return range(n).map(i=>range(n).map(j=>i===j?0:r.has((j-i+n)%n)?1:-1))
}
function transport(B,root){
 const parents=range(B.length).filter(i=>B[root][i]);
 const parts=parents.map(i=>range(B.length).filter(j=>j!==root&&B[i][j]));
 const P=parts.map((part,i)=>parts.map((next,j)=>i===j?null:part.map(a=>next.findIndex(b=>B[a][b]))));
 for(let i=0;i<P.length;i++)for(let j=0;j<P.length;j++)if(i!==j)assert(new Set(P[i][j]).size===parts[0].length&&P[i][j].every(x=>x>=0));
 return P
}
const comp=(f,g)=>g.map(x=>f[x]), inv=f=>f.map((_,i)=>f.indexOf(i));
function group(P,b){
 const n=P.length,m=n-1,one=range(m),F=P[b].map(p=>p||one);
 const gen=[];
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){
  const g=comp(inv(F[j]),comp(P[i][j],F[i]));
  gen.push(g)
 }
 const seen=new Set([one.join(",")]),queue=[one];
 for(const f of queue)for(const g of gen){
  const h=comp(f,g),k=h.join(",");
  if(!seen.has(k)){seen.add(k);queue.push(h)}
 }
 return seen.size
}
for(const [B,n,order] of [[pet,3,2],[hs,7,720]]){
 for(let r=0;r<B.length;r++){
  const P=transport(B,r);
  for(let b=0;b<n;b++)assert(group(P,b)===order);
 }
 console.log("degree",n,"root/reference choices",n*B.length,"transport group order",order)
}
function checkExpansion(B,d){
 const P=transport(B,0),m=d-1,S=paley(d),g=sq(S);
 const C=range(d).map(i=>range(m).map(a=>range(m).map(b=>
   a===b?0:(a<b?1:-1)*((i+a+b)%3===0?-1:1))));
 const r=C.map(B=>B.map(row=>row.reduce((s,x)=>s+x,0)));
 const A=range(d*m).map(a=>range(d*m).map(b=>{
  const i=a/m|0,j=b/m|0;
  return i===j?C[i][a%m][b%m]:S[i][j]*(P[i][j][a%m]===b%m?-1:1)
 }));
 const G=sq(A);let k=0;
 for(let i=0;i<d;i++)for(let j=i+1;j<d;j++)for(let a=0;a<m;a++)for(let b=0;b<m;b++){
  let q=0,count=+(P[i][j][a]===b);
  for(let h=0;h<d;h++)if(h!==i&&h!==j&&P[h][j][P[i][h][a]]===b){q+=S[i][h]*S[h][j];count++}
  assert(count===1);
  const pred=(m-4)*g[i][j]+4*q+S[i][j]*(r[i][a]-r[j][b]-2*(C[i][a][inv(P[i][j])[b]]+C[j][P[i][j][a]][b]));
  assert(pred===G[i*m+a][j*m+b]);k++;
 }
 return k
}
console.log("Hoffman–Singleton exact expanded cross-block entries",checkExpansion(hs,7));
const D=paley(59),S=D.slice(0,57).map(r=>r.slice(0,57)),G=sq(S);
const F=G.map((r,i)=>r.map((v,j)=>i!==j&&Math.abs(v)===3&&S[i][j]*v>0));
const u=F.map(row=>row.filter(Boolean).length);
const v=F.map((_,i)=>F.filter(row=>row[i]).length);
assert(Math.min(...u)>0&&Math.min(...v)>0);
const walk=[16,8,6,11,0,4];
assert(walk.slice(1).every((i,k)=>F[walk[k]][i]));
const grp={};
for(let i=0;i<57;i++){
 const key=D[i].slice(57).join(",");
 grp[key]=(grp[key]||0)+1;
 for(let j=i+1;j<57;j++)
  assert(G[i][j]===1+D[i][57]*D[j][57]+D[i][58]*D[j][58])
}
console.log("Paley-59 deletion pattern classes",grp,"min outgoing/incoming defect degrees",Math.min(...u),Math.min(...v),"walk",walk);
const round=x=>2*Math.ceil(x/2);
assert(round((3*54**2+52*5)/56)===162&&round(5*54**2/56)===262);
assert(148+55-(-55)===258);
console.log("score-sensitive bounds checked: R>=5 -> 162; two transitive interiors -> 258.");
