"use strict";
const assert=require("node:assert/strict");
const R=n=>Array.from({length:n},(_,i)=>i);
const c=(k,s)=>Math.max(0,Math.ceil(k*(k-8)*(56-s)/112));
const t=(k,r)=>Math.max(0,Math.ceil((r-1)*k*(k-8)/(2*r)));
function bounds(f){
 const B=[0,0],S=[0,0];
 for(let k=2;k<=15;k++){
  S[k]=Math.min(...R(k-1).map(i=>{const a=i+1,b=k-a;return a*b+B[a]+B[b];}));
  B[k]=Math.min(f(k),S[k]);
 }
 return [B,S];
}
for(const [label,f,want] of [
 ["arbitrary",k=>c(k,27),[22,28,30,36]],
 ["scores <= 3",k=>c(k,26),[23,29,31,37]],
 ["scores <= 3 and positive holonomy",k=>t(k,4),[32,40,38,46]]
]){
 const [B,S]=bounds(f);
 const got=[B[14],B[15],S[14],S[15]];
 assert.equal(got.join(),want.join());
 console.log(label,"cliques",B[14],B[15],"split",S[14],S[15],"total",3*B[14]+B[15]);
}
const V=R(50).map(i=>[Math.floor(i/25),Math.floor(i/5)%5,i%5]);
const mod=x=>(x%5+5)%5;
const A=V.map(([t,i,j])=>V.map(([u,k,l])=>{
 if(t===u)return +(i===k&&(mod(j-l)===(t?2:1)||mod(l-j)===(t?2:1)));
 if(t===0)return +(mod(j+i*k)===l);
 return +(mod(l+i*k)===j);
}));
for(let i=0;i<50;i++)for(let j=0;j<50;j++){
 const s=A[i].reduce((v,x,k)=>v+x*A[k][j],0);
 assert.equal(s,i===j?7:1-A[i][j]);
}
const branch=R(50).filter(i=>A[0][i]).map(i=>R(50).filter(j=>j!==0&&A[i][j]));
const m=(i,j,v)=>branch[j].find(w=>A[v][w]);
let checked=0;
for(let i=0;i<7;i++)for(let j=i+1;j<7;j++)for(let k=j+1;k<7;k++){
 const h=v=>m(k,i,m(j,k,m(i,j,v)));
 assert(branch[i].every(v=>h(h(v))===v&&h(v)!==v));
 const O=[branch[i][0],h(branch[i][0])], Y=O.map(v=>m(i,j,v)), Z=O.map(v=>m(i,k,v));
 assert(Y.every(v=>Z.includes(m(j,k,v))));
 const x=R(50).map(v=>[i,j,k].some(u=>branch[u].includes(v))?(O.includes(v)||Y.includes(v)||Z.includes(v)?4:-2):0);
 assert(x.every((v,i)=>A[i].reduce((s,w,j)=>s+w*x[j],0)===2*v));
 ++checked;
}
assert.equal(checked,35);
console.log("Hoffman–Singleton: all 35 branch triangles have coherent two-slot orbits and attain eigenvalue 2.");
