// Independent marked fixture reconstruction, cycle census, and lift checks.
const fs=require("node:fs");

const N=56, I=Array.from({length:N},(_,i)=>i), S=I.filter(i=>i!==0&&i!==28), centers=[56,0,28];
const mod=(x,n=56)=>((x%n)+n)%n, inv=f=>{const g=[];f.forEach((y,x)=>g[y]=x);return g};
function board(text) {
 const L=text.trim().split(/\n/).map(x=>x.split(/\s+/)), m=+L[0][1], trades=L.slice(1,m+1).map(r=>r.map(Number)), selected=[L[m+1].slice(1).map(Number),...L.slice(m+2).map(r=>r.map(Number))];
 const shift=t=>I.map(x=>mod(x+t)), s=[];
 for(const [u,v,w] of [[9,10,19],[8,12,20],[7,11,18],[6,15,21],[5,22,27],[4,13,17],[3,23,26],[2,14,16],[1,24,25]]){
  const R=[u,v,-w];for(let i=0;i<3;i++)for(const a of [-1,1])s[mod(a*R[i])]=mod(-a*R[(i+1)%3]);
 }
 function chart(steps) {
  const p=[],q=[],C=Array.from({length:57},()=>[]);
  for(const h of S){p[h]=shift(-h);q[h]=shift(-s[h]);}
  for(const [h,k,start] of steps){
   const ip=inv(p[k]),iq=inv(q[k]),o=new Set([start]);
   for(const x of o){o.add(ip[p[h][x]]);o.add(iq[q[h][x]]);}
   if(o.size===56)throw Error("bad trade");
   for(const x of o){[p[h][x],p[k][x]]=[p[k][x],p[h][x]];[q[h][x],q[k][x]]=[q[k][x],q[h][x]];}
  }
  function put(a,b,f){C[a][b]=f;C[b][a]=inv(f);}
  for(const h of I)put(56,h,I);
  put(0,28,shift(28));
  for(const h of S){put(h,0,p[h]);put(h,28,q[h]);}
  return C;
 }
 const E=chart(trades), C=chart([...trades,[10,45,2],[27,49,0]]);
 const triple=(m,h,x)=>centers.map(c=>m[h][c][x]);
 const loc=new Map();
 for(const h of S)for(const x of I)loc.set(triple(C,h,x).join(),[h,x]);
 const a=Array.from({length:3250},()=>new Set()), id=(h,x)=>58+56*h+x;
 const join=(A,u,v)=>{A[u].add(v);A[v].add(u);};
 for(let h=0;h<57;h++){
  join(a,0,1+h);
  for(const x of I)join(a,1+h,id(h,x));
  for(let j=h+1;j<57;j++)if(C[h][j])for(const x of I)join(a,id(h,x),id(j,C[h][j][x]));
 }
 for(let d=0;d<56;d+=4){
  const [src,...targets]=selected.map(([h,x])=>loc.get(triple(E,h,x).map(v=>mod(v+d)).join()));
  for(const t of targets)join(a,id(...src),id(...t));
 }
 const baseEdges=a.reduce((n,x)=>n+x.size,0)/2;
 if(baseEdges!==13231)throw Error("base");
 const v=id(1,1), mark=id(2,7);
 function legal(A,v,w){
  const Q=new Set([v,...A[v]]);
  for(const u of [...Q])for(const z of A[u])Q.add(z);
  return !Q.has(w)&&![...A[w]].some(j=>Q.has(j));
 }
 if(!legal(a,v,mark))throw Error("marker");
 join(a,v,mark);
 const R=[];
 for(const h of S.filter(h=>h!==1&&h!==2))for(const x of I){
  if(!legal(a,v,id(h,x)))continue;
  R.push({h,x,t:triple(C,h,x)});
 }
 if(R.length!==1741)throw Error("full board");
 return {a,R,v,id,join,legal};
}

function replay(result,wit,C){
 const n=C.length, ok=(x,msg)=>{if(!x)throw Error(msg);};
 const A=Array.from({length:n},()=>new Int8Array(n)), sign=(a,b,k)=>1^(a.t[k]>=28)^(b.t[k]>=28);
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){
  let g=-1, hits=0;
  for(let k=0;k<3;k++)if(C[i].t[k]%28===C[j].t[k]%28){g=k;hits++;}
  A[i][j]=A[j][i]=C[i].h===C[j].h||hits>1?-1:hits===0?0:sign(C[i],C[j],g)+1;
 }
 let hist=[0,0];
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)if(A[i][j]===0){
  const X=[];for(let k=i+1;k<n;k++)if(A[i][k]>0&&A[j][k]>0)X.push(k);
  for(let x=0;x<X.length;x++)for(let y=x+1;y<X.length;y++)if(A[X[x]][X[y]]===0){
   const k=X[x],l=X[y];
   hist[(A[i][k]^A[j][k]^A[j][l]^A[i][l])&1]++;
  }
 }
 ok(hist.join()==="487750,488133","4-cycle counts"); // least cycle vertex is i
 const W=wit.map(([h,y])=>C.findIndex(e=>e.h===h&&e.x===y));
 const L=W.length,P=W.map(i=>C[i]);let xor=0;
 for(let i=0;i<L;i++)for(let j=i+1;j<L;j++){
  const a=A[W[i]][W[j]], near=j===i+1||(i===0&&j===L-1);
  ok(near?a>0:a===0,"induced");if(near)xor^=a-1;
 }
 ok(xor===1,"frustration");
 for(let d=0;d<L;d++){
  let b=0;
  const adj=result.a.map(s=>new Set(s));
  for(let j=1;j<L;j++){
   const i=(d+j)%L, e=P[i];
   if(j>1)b^=A[W[(i+L-1)%L]][W[i]]-1;
   const v=result.id(e.h,e.x+28*b);
   ok(result.legal(adj,result.v,v),"illegal lift");result.join(adj,result.v,v);
  }
  ok(adj[result.v].size===L+4,"physical degree");
 }
 const S=P.map(p=>p.t.map(z=>z%28));
 const sizes=[0,1,2].map(g=>new Set(S.map(t=>t[g])).size);
 return {hist,sizes,P};
}

function extend(result,W){
 const rest=result.R, capacity= [0,1,2].map(g=>Array.from({length:28},(_,r)=>new Set(rest.filter(t=>t.t[g]%28===r).map(t=>t.t[g])).size));
 const pri=W.map(([h,y])=>rest.find(e=>e.h===h&&e.x===y));
 const range=[...new Set(rest.map(x=>x.h))];
 let tot=0;
 function run(P,fold){
  const sz=fold?28:56, V=[0,1,2].map(g=>fold?[...capacity[g]]:Array.from({length:56},(_,r)=>+rest.some(e=>e.t[g]===r)));
  const H=new Map(range.map(h=>[h,rest.filter(e=>e.h===h).map(e=>({h:e.h,x:e.x,s:e.t.map(z=>z%sz)}))]));
  const W=P.map(p=>({h:p.h,x:p.x,s:p.t.map(z=>z%sz)}));
  const dec=(e,c)=>e.s.forEach((v,g)=>V[g][v]+=c);
  for(const p of W){dec(p,-1);H.delete(p.h);}
  let count=0,dead=[];
   for(const [h,D] of H)if(!D.some(e=>e.s.every((s,g)=>V[g][s]>0)))dead.push(h);
  function f(H){
   count++;
   if(!H.size)return true;
   let Q,b=-1;
   for(const [h,D] of H){
    const A=D.filter(e=>e.s.every((s,g)=>V[g][s]>0));
    if(!A.length)return false;
    if(!Q||A.length<Q.length){Q=A;b=h;}
   }
   const J=new Map(H);J.delete(b);
   for(const e of Q){dec(e,-1);if(f(J))return true;dec(e,1);}
   return false;
  }
  const good=f(H);tot+=count;return [good,count,dead];
 }
 const F=run(pri,true);
 const phys=[];
 for(let l=0;l<pri.length;l++)for(let b0=0;b0<2;b0++){
  let b=b0, A=[];
  for(let j=1;j<pri.length;j++){
   let i=(l+j)%pri.length,e=pri[i],p=pri[(i+pri.length-1)%pri.length];
   if(j>1){const g=e.t.findIndex((v,k)=>v%28===p.t[k]%28);b^=1^(e.t[g]>=28)^(p.t[g]>=28);}
   A.push(rest.find(p=>p.h===e.h&&p.x===e.x+28*b));
  }
  phys.push(run(A,false));
 }
 return [F,phys,tot];
}
const old=[[15,19],[4,8],[12,8],[50,18],[8,18],[33,4],[37,15],[9,15],[23,22],[29,22],[42,11],[6,17],[18,21],[30,5],[7,9],[32,9],[39,12],[46,13],[5,13],[34,20],[26,20],[10,26],[41,26],[19,14],[17,6],[31,6],[20,3],[35,3],[44,23],[52,23],[11,16],[38,16],[24,2],[48,24],[43,19]];
const best=[[38,16],[24,2],[48,24],[43,19],[15,19],[4,8],[12,8],[50,18],[8,18],[33,4],[37,15],[9,15],[23,22],[29,22],[42,11],[6,17],[18,21],[30,5],[7,9],[32,9],[39,12],[46,13],[21,13],[44,23],[52,23],[3,3],[20,3],[31,6],[17,6],[19,14],[55,27],[41,27],[34,20],[26,20],[10,26],[13,16]];

const result=board(fs.readFileSync(process.argv[2]??"period_fourteen_site.txt","utf8"));
const C=result.R.filter(e=>e.x<28&&result.R.some(f=>e.h===f.h&&f.x===e.x+28&&f.t.every((v,j)=>v===(e.t[j]+28)%56)));
if(C.length!==677)throw Error("two-phase options");
const cap=[0,1,2].map(g=>Array.from({length:28},(_,r)=>new Set(result.R.filter(t=>t.t[g]%28===r).map(t=>t.t[g])).size));
if(cap.flat().filter(n=>n===2).length!==73)throw Error("double capacities");
for(const W of [old,best]){
 const {hist,sizes}=replay(result,W,C);
 const [check,phys,nodes]=extend(result,W);
 if(check[0]||phys.some(x=>x[0]||!x[2].length))throw Error("unexpected extension");
 if(W.length===36&&sizes.join()!=="23,24,25")throw Error("cycle residues");
 console.log("Length",W.length,"four-cycle counts",hist,"used residues",sizes,
  "residue completions",Number(check[0]),"search nodes",check[1],
  "physical deletions embedded",W.length,"and",phys.length,"oriented deletions rejected by a dead branch");
}
