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
const calc=(n,d)=>n*d*(d-1)**2*(d-2)*(d-3)/14;
function run(A){
 const n=A.length, d=A[0].reduce((a,b)=>a+b), N=A.map(r=>R(n).filter(i=>r[i]));
 function common(i,j){return N[i].filter(k=>A[k][j]);}
 for(let i=0;i<n;i++){assert.equal(N[i].length,d);for(let j=i+1;j<n;j++)assert.equal(common(i,j).length,1-A[i][j]);}
 let counts=new Map(), C21=new Set(), cycle=[], use=new Set(), tot=0;
 function add(alt){
  let s=alt.indexOf(Math.min(...alt)), a=alt.slice(s).concat(alt.slice(0,s));
  if(a[1]>a.at(-1))a=[a[0],...a.slice(1).reverse()];
  assert(a.every((x,i)=>A[x][a[(i+1)%a.length]])&&new Set(a).size===21);
  const k=a.join(",");assert(!C21.has(k));C21.add(k);
 }
 function record(){
  ++tot;
  let f=i=>cycle[mod(i,7)];
  let W=R(7).map(i=>common(f(i),f(i+3))[0]), Z=R(7).map(i=>common(W[i],f(i+5))[0]);
  let all=[...cycle,...W,...Z];
  assert.equal(new Set(all).size,21);
  const count=(v,S)=>S.filter(w=>A[v][w]).length;
  assert(W.every((v,i)=>count(v,cycle)===2&&count(v,W)===0&&count(v,Z)===1&&A[v][Z[i]]));
  assert(Z.every((v,i)=>count(v,cycle)===1&&A[v][f(i+5)]&&Z.every((w,j)=>!A[v][w]||[3,4].includes(mod(j-i,7)))));
  for(const a of [0,3]){
   const alt=[];
   let at=0;
   do{const i=mod(at-a,7);alt.push(f(at),W[i],Z[i]);at=mod(i+5,7);}while(at!==0);
   add(alt);
  }
  const Y=R(n).filter(v=>!cycle.includes(v)&&count(v,cycle)===0);
  const x=R(4).map(i=>Y.filter(v=>count(v,W)===i).length);
  const b=Z.reduce((sum,v)=>sum+count(v,Z),0)/2;
  assert.equal(Y.length,d*d-7*d+15);
  assert.equal(x.reduce((a,b)=>a+b),Y.length);
  assert.equal(x[2]+3*x[3],14);assert.equal(x[1]+2*x[2]+3*x[3],7*(d-3));
  assert(x[3]<=3);
  for(const v of Y){
   const k=count(v,W);
   assert.equal(count(v,Y),d-7+k);
  }
  for(let i=0;i<7;i++){
   const U=Y.filter(v=>A[v][W[i]]), L=Y.filter(v=>A[v][Z[i]]);
   assert.equal(U.length,d-3);assert.equal(L.length,d-4);
   for(const v of Y)assert.equal(+U.includes(v)+count(v,U)+(+L.includes(v)),1);
   for(let j=0;j<7;j++){
    const sep=mod(j-i,7);
    const cross=i===j||[2,5].includes(sep)||A[Z[i]][Z[j]]?0:1;
    assert.equal(U.filter(v=>A[v][Z[j]]).length,cross);
   }
  }
  if(d===7){
   const seen=new Set(), parts=[];
   for(const y of Y)if(!seen.has(y)){
    const Q=[y];seen.add(y);for(const z of Q)for(const u of Y)if(A[z][u]&&!seen.has(u)){seen.add(u);Q.push(u);}
    parts.push(Q.length);
   }
   const k=x[3]+","+b+"; "+parts.sort((a,b)=>a-b).join("+");
   counts.set(k,(counts.get(k)||0)+1);
  }
 }
 function dfs(first,last){
  if(cycle.length===7){if(A[last][first]&&cycle[1]<last) record();return;}
  for(const z of N[last])if(z>first&&!use.has(z)){
   use.add(z);cycle.push(z);dfs(first,z);cycle.pop();use.delete(z);
  }
 }
 for(let i=0;i<n;i++){cycle=[i];use=new Set([i]);dfs(i,i);}
 assert.equal(tot,calc(n,d));assert.equal(C21.size,2*tot);
 console.log("degree",d,"7-cycles",tot,"distinct canonical 21-cycles",C21.size,"(triple witnesses, special edges; residual component sizes):",JSON.stringify([...counts].sort()));
}
run(petersen());run(hoffman());
assert.equal(calc(3250,57),123243120000);
console.log("conditional degree 57: 123243120000 induced 7-cycles; at least 246486240000 ordinary 21-cycles.");

for(const [r,target] of [[1,50],[2,98],[3,146],[4,193],[5,239]]){
 const c=r===5?7:r+1;
 assert(1849*c>target*(c+42));
}
assert(50-367*53/2498>42.2);
assert(422n**4n-10000n*422n*1000n>2730000n*10000n);
assert.equal(5*50-8-3,239);
console.log("Degree-57 cut arithmetic: 50, 98, 146, 193, 239; the induced-path count clears the three bad-pair bounds.");
