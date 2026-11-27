
"use strict";
const assert=x=>{if(!x)throw Error("replay failure")};
const range=n=>Array.from({length:n},(_,i)=>i);
const car=n=>range(n).map(i=>range(n).map(j=>i===j?0:(j-i+n)%n<=n/2?1:-1));
const m=56;
const C=car(57).slice(0,m).map(r=>r.slice(0,m)),opp=C.map(r=>r.map(x=>-x));
const score=A=>A.map(r=>r.reduce((s,v)=>s+v,0));
assert(score(C).every(x=>Math.abs(x)===1));
const q=range(m).map(a=>range(m).map(b=>a===b?0:(b-a+m)%m<=29?1:-1));
for(let a=0;a<m;a++){
 assert(q[a].reduce((s,x)=>s+x,0)===3);
 assert(q.reduce((s,r)=>s+r[a],0)===3)
}
function block(X,Y,q) {
 const x=score(X),y=score(Y);
 return X.map((r,a)=>r.map((_,b)=>156+4*q[a][b]+x[a]-y[b]-2*(X[a][b]+Y[a][b])));
}
const mx=V=>Math.max(...V.flat().map(Math.abs));
assert(mx(block(C,opp,q))===162);
assert(mx(block(C,C,q))===166);
console.log("Synthetic single routed block: anti-aligned max",mx(block(C,opp,q)),"equal aligned max",mx(block(C,C,q)));
const phases=new Map();
function type(K){
 if(K<=-2)return 0;
 if(K<=0)return 1;
 if(K<=2)return 2;
 return 3;
}
for(const B of [160,162,164,166]){
 let dec=[];
 for(let q of [-1,0,1])for(let d=-6;d<=6;d+=2){
  const K=(156+4*q+d-B)/2,t=type(K);
  for(const a of [-1,1])for(const b of [-1,1]){
   const pred=t===0||t===1&&(a===1||b===1)||t===2&&a===1&&b===1;
   assert(pred===(156+4*q+d-2*(a+b)<=B));
  }
  if(q===1&&Math.abs(d)<=2)dec.push(t)
 }
 phases.set(B,dec);
}
assert(JSON.stringify([...phases.values()])==="[[1,1,2],[0,1,1],[0,0,1],[0,0,0]]");
console.log("Semiregular thresholds at positive routes, score differences -2,0,2:",JSON.stringify([...phases]));
const tail=[1,1,-1,-1],h=[2,3,0,1];
const inv=h.map((_,i)=>h.indexOf(i));
assert(h.every((x,i)=>h[x]===i&&tail[x]===-tail[i]));
// Explicit signed-orbit compatibility test on all ordered pairs.
function signedOrbit(n,gens){
 let seen=new Map();
 for(let a=0;a<n;a++)for(let b=a+1;b<n;b++){
  const id=a*n+b;
  if(seen.has(id))continue;
  seen.set(id,1);
  const Q=[[a,b,1]];
  for(const [x,y,sgn] of Q)for(const [p,s] of gens){
   let z=p[x],w=p[y],t=sgn*s;
   if(z>w){[z,w]=[w,z];t=-t}
   const m=z*n+w;
   if(seen.has(m)){if(seen.get(m)!==t)return false}
   else {seen.set(m,t);Q.push([z,w,t])}
  }
 }
 return true
}
assert(signedOrbit(4,[[h,-1]]));
assert(!signedOrbit(4,[[h,1]]));
assert(!signedOrbit(4,[[range(4),-1]]));
console.log("Signed transport orbit controls passed.");

function paley(p){
 const sq=new Set(range(p-1).map(x=>(x+1)**2%p));
 return range(p).map(i=>range(p).map(j=>i===j?0:sq.has((j-i+p)%p)?1:-1));
}
const P59=paley(59),classes=new Map();
for(let i=0;i<57;i++){
 const k=P59[i].slice(57).join(",");
 if(!classes.has(k))classes.set(k,[]);
 classes.get(k).push(i);
}
assert(JSON.stringify([...classes.values()].map(a=>a.length).sort())==="[14,14,14,15]");
for(const group of classes.values())for(let i of group)for(let j of group)if(i!==j){
 let s=0;for(let k=0;k<57;k++)s+=P59[i][k]*P59[k][j];
 assert(s===3)
}
const parts=[];
for(let a=13;a<=28;a+=2)for(let b=a;b<=28;b+=2)if(a+b===28)parts.push([2*a,2*b]);
assert(JSON.stringify(parts)==="[[26,30]]");
assert(!range(4).some(x=>range(4).some(y=>x*x+y*y===14)));
console.log("Paley defect classes 14,14,14,15; conditional orbit split 26+30; order 14 is nonsquare-sum.");
function validatePlane(q){
 const V=[];
 for(let x=0;x<q;x++)for(let y=0;y<q;y++)V.push([x,y,1]);
 for(let x=0;x<q;x++)V.push([x,1,0]);
 V.push([1,0,0]);
 const has=(l,p)=>l.reduce((s,x,i)=>s+x*p[i],0)%q===0;
 const inf=V.findIndex(v=>v.join(",")==="0,0,1"), lInf=inf;
 const oldPoints=range(V.length).filter(i=>i!==inf&&!has(V[lInf],V[i]));
 const groupLines=range(V.length).filter(i=>has(V[i],V[inf]));
 const oldLines=range(V.length).filter(i=>i!==lInf&&!has(V[i],V[inf]));
 const corners=groupLines.map(l=>range(V.length).find(i=>has(V[l],V[i])&&has(V[lInf],V[i])));
 const pointGroup=oldPoints.map(i=>groupLines.findIndex(l=>has(V[l],V[i])));
 const lineGroup=oldLines.map(l=>corners.findIndex(i=>has(V[l],V[i])));
 assert(groupLines.length===q+1 && oldPoints.length===(q-1)*(q+1) && oldLines.length===oldPoints.length);
 const bare=oldLines.map((l,j)=>oldPoints.flatMap((p,i)=>has(V[l],V[p])?[i]:[]));
 for(let j=0;j<bare.length;j++)assert(!bare[j].some(i=>pointGroup[i]===lineGroup[j])&&new Set(bare[j].map(i=>pointGroup[i])).size===q);
 for(let i=0;i<oldPoints.length;i++)for(let j=i+1;j<oldPoints.length;j++)
  assert(bare.filter(l=>l.includes(i)&&l.includes(j)).length===+(pointGroup[i]!==pointGroup[j]));
 const K=oldPoints.length;
 const aug=bare.map((l,j)=>l.concat(K+lineGroup[j]));
 for(let j=0;j<=q;j++)aug.push(range(K).filter(i=>pointGroup[i]===j).concat(K+j,K+q+1));
 aug.push(range(q+1).map(i=>K+i));
 assert(aug.length===q*q+q+1&&aug.every(l=>l.length===q+1));
 for(let i=0;i<aug.length;i++)for(let j=i+1;j<aug.length;j++)
  assert(aug.filter(l=>l.includes(i)&&l.includes(j)).length===1);
 console.log("Projective-plane completion cross-check, order",q,"points",aug.length)
}
validatePlane(3);validatePlane(5);
