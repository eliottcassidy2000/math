// Run: node replay.js
const assert = require("node:assert/strict");
const V = [0,1,2,3], free = ["02","03","13"];
function perms(A) {
  return A.length ? A.flatMap((x,i) =>
    perms(A.filter((_,j)=>i!==j)).map(p=>[x,...p])) : [[]];
}
const P = perms(V);
const phi = x => x.split("").map(d=>d==="1"?"0":"1").join("");
const tau = x => x.split("").reverse().join("");
function arc(x,u,v) {
  if (u>v) return !arc(x,v,u);
  return x[free.indexOf(""+u+v)]!=="1";
}
function cert(x) {
  const H = P.filter(p=>p.slice(1).every((v,j)=>arc(x,p[j],v)));
  const C = H.filter(p=>arc(x,p[3],p[0])).length/4;
  const score=V.map(u=>V.filter(v=>v!==u&&arc(x,u,v)).length).sort();
  const code=P.map(p=>p.flatMap((u,i)=>
    p.slice(i+1).map(v=>Number(arc(x,u,v)))).join("")).sort()[0];
  let t=0;
  for(let i=0;i<4;i++) for(let j=i+1;j<4;j++) for(let k=j+1;k<4;k++)
    if(arc(x,i,j)===arc(x,j,k)&&arc(x,j,k)===arc(x,k,i))t++;
  return {H:H.length,C,score,code,t};
}
const rows=Array.from({length:8},(_,i)=>i.toString(2).padStart(3,"0"));
const classes = new Map(), orbits = new Set(), flip5 = new Set();
for (const x of rows) {
  const y = cert(x), r = cert(tau(x)), f = cert(phi(x));
  const [a,b,c]=[...x].map(Number);
  assert.equal(y.t,a+c+b*(2-a-c));
  assert.equal(y.H,1+2*y.t);
  assert.equal(y.H,r.H);
  assert.equal(y.C,r.C);
  assert.equal(phi(tau(x)),tau(phi(x)));
  assert.equal(phi(phi(x)),x);
  assert.equal(tau(tau(x)),x);
  assert.notEqual(x,phi(x));
  assert.notEqual(x,tau(phi(x)));
  if(y.H===5)flip5.add(f.H);
  classes.set(y.code,(classes.get(y.code)??0)+1);
  orbits.add([...new Set([x,tau(x),phi(x),tau(phi(x))])].sort().join("/"));
  console.log(x, y.H, y.C, y.score.join(","));
}
assert.deepEqual([...classes.values()].sort(),[1,1,1,5]);
assert.deepEqual([...flip5].sort(),[1,3,5]);
assert.deepEqual([...orbits].sort(),
                 ["000/111","001/011/100/110","010/101"]);
assert.equal(cert("101").C,1);
assert.equal(cert("111").C,1);
