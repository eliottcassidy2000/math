function example(k){
 const n=k*k+1, g=Array.from({length:n},()=>Array(n).fill(0));
 function add(i,j){g[i][j]=g[j][i]=1}
 if(k==7) {
  for(let i=0;i<5;i++)for(let x=0;x<5;x++){
   add(i*5+x,i*5+(x+1)%5);add(25+i*5+x,25+i*5+(x+2)%5);
   for(let j=0;j<5;j++)add(i*5+x,25+j*5+(x+i*j)%5)
  }
 }else for(let x=0;x<5;x++){
  add(x,(x+1)%5);add(5+x,5+(x+2)%5);add(x,5+x)
 }
 const head=g[0].flatMap((e,i)=>e?[i]:[]);
 const br=head.map(h=>g[h].flatMap((e,i)=>(e&&i!=0)?[i]:[]));
 for(let i=1;i<k;i++)br[i]=br[0].map(u=>br[i].find(v=>g[u][v]));
 const D=br.map((b,i)=>br.map((c,j)=>i==j?null:b.map(u=>c.findIndex(v=>g[u][v]))));
 return D
}
function perms(a){
 if(a.length===0)return [[]];
 return a.flatMap((v,i)=>perms(a.filter((_,j)=>j!==i)).map(w=>[v,...w]));
}
function exact(D){
 const k=D.length,w=k-1, labels=perms(Array.from({length:w},(_,i)=>i)), branches=labels.map(a=>[0,...a.map(x=>x+1)]);
 const hist=Array(k*(k-1)/2+1).fill(0);let min=Infinity,witness;
 for(const T of labels){
  const Z=D.map(row=>row.map(col=>col&&T.map(x=>col[x])));
  const inv=[];for(let x=0;x<w;x++)inv[T[x]]=x;
  for(const I of branches){
   let edges=0;
   for(let i=1;i<k;i++)for(let j=i+1;j<k;j++){
    let a=Z[I[i]][I[j]],b=D[i][j];
    for(let x=0;x<w;x++)if(a[x]!=T[b[x]]){edges++;break;}
   }
   hist[edges]++;
   if(edges&&edges<min){min=edges;witness={I,T}};
  }
 }
 return {k,hist,min,witness};
}
const assert=require('node:assert/strict');
const p=exact(example(3)),h=exact(example(7));
assert.equal(p.hist[0],4);
assert.equal(h.min,8);
assert.deepEqual(h.hist.slice(0,16),[720,0,0,0,0,0,0,0,10800,0,0,0,72000,0,216000,218880]);
console.log(p,h);
