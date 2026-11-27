"use strict";
const assert=require("node:assert/strict");
const R=n=>Array.from({length:n},(_,i)=>i);
const mod=(x,n)=>(x%n+n)%n;
function hoffman(){
 const V=R(50).map(i=>[i>=25?1:0,Math.floor(i/5)%5,i%5]);
 return V.map(([t,i,j])=>V.map(([u,k,l])=>t===u?+(i===k&&[t?2:1,t?3:4].includes(mod(j-l,5))):t===0?+(mod(j+i*k,5)===l):+(mod(l+i*k,5)===j)));
}
function petersen(){
 const V=R(10).map(i=>[i>=5?1:0,i%5]);
 return V.map(([t,i])=>V.map(([u,j])=>t===u?+[t?2:1,t?3:4].includes(mod(i-j,5)):+(i===j)));
}

function counts(A){
 const n=A.length,N=A.map(r=>R(n).filter(i=>r[i])),d=N[0].length;
 let local=0, allFixed=0, mob=0, six=0,seven=0, bSum=0;
 function closure(S){
  const Z=new Set(S),seq=[Z.size];
  for(;;){
   const add=R(n).filter(x=>!Z.has(x)&&N[x].filter(y=>Z.has(y)).length>=2);
   if(!add.length)return seq;
   for(const x of add) Z.add(x);
   seq.push(Z.size);
  }
 }
 const cm=(a,b)=>N[a].filter(x=>A[b][x]);
 for(let o=0;o<n;o++){
  const B=N[o].map(i=>N[i].filter(j=>j!==o));
  const m=(i,j,x)=>B[j].find(y=>A[x][y]);
  for(let i=0;i<d;i++)for(let j=i+1;j<d;j++)for(let k=j+1;k<d;k++){
   const h=x=>m(k,i,m(j,k,m(i,j,x)));
   assert(B[i].every(x=>h(x)!==x));
   const not2=B[i].filter(x=>h(h(x))!==x).length;
   local+=(d-1-not2)/2;
   mob+=3*not2;
   allFixed+=!not2;
  }
 }
 function walk(path, first){
  let len=path.length,x=path.at(-1);
  if(len>=6&&path[1]<x&&A[x][first]){
   if(len===6)six++;
   if(len===7){
    seven++;
    const W=R(7).map(i=>cm(path[i],path[(i+3)%7])[0]);
    const b=R(n).filter(y=>!path.includes(y)&&!N[y].some(x=>path.includes(x))&&W.filter(w=>A[y][w]).length===3).length;
    bSum+=b;
    const seq=closure(path);assert.equal(seq.slice(0,3).join(),[7,14,35-2*b].join());
    if(d===7)assert.equal(seq.join(),"7,14,35,50");
   }
  }
  if(len===7)return;
  for(const y of N[x])if(y>first&&!path.includes(y))walk([...path,y],first);
 }
 for(let x=0;x<n;x++)walk([x],x);
 assert.equal(six,n*d*(d-1)**2*(d-2)/12);
 assert.equal(bSum,mob);
 assert.equal(six,local);
 assert.equal(local%10,0);
 console.log("degree",d,"hexagons",six,"Petersens",local/10,"involutional root triangles",allFixed,"heptagons",seven,"triple witnesses",bSum);
}
counts(petersen());counts(hoffman());
const n=3250,d=57,localMax=(d-1)*d*(d-1)*(d-2)/12,globalMax=n*localMax/10;
assert.equal(localMax,819280);
assert.equal(globalMax,266266000);
console.log("degree-57 maxima: rooted",localMax,"total",globalMax);
