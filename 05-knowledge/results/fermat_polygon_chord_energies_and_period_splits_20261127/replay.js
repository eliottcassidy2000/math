const assert=require('node:assert/strict');
const p=17,g=10,sq=Math.sqrt(17);
const units=[1];for(let j=1;j<16;j++)units.push(g*units[j-1]%17);
assert.equal(new Set(units).size,16);assert.equal(g*g%17,15);
const C=(j,k)=>units.filter((_,i)=>i%2**j===k);
function count(X,Y){
 let f=Array(17).fill(0);for(const x of X)for(const y of Y)f[(x+y)%17]++;
 return f;
}
function row(j,k){
 const f=count(C(j+1,k),C(j+1,k+2**j));
 const r=[f[0]];for(let a=0;a<2**j;a++){
  const c=C(j,a).map(x=>f[x]);assert.ok(c.every(x=>x===c[0]));r.push(c[0]);
 }return r;
}
assert.deepEqual(row(0,0),[0,4]);
assert.deepEqual(row(1,0),[0,1,1]);assert.deepEqual(row(1,1),[0,1,1]);
assert.deepEqual(row(2,0),[0,0,0,0,1]);
for(let j=0;j<=2;j++)for(let k=0;k<2**j;k++)console.log(j,k,JSON.stringify(row(j,k)));
const S=(sq-1)/2,N=(-sq-1)/2,D=Math.sqrt(S*S+4),E=2*sq/D;
assert.ok(Math.abs(E*E-(N*N+4))<1e-12);
const A=(S+D)/2,B=(N+E)/2,W=(A+Math.sqrt(A*A-4*B))/2;
assert.ok(Math.abs(W-2*Math.cos(2*Math.PI/17))<1e-12);
const energies=[0,1].map(j=>C(1,j).reduce((s,x)=>s+2-2*Math.cos(2*Math.PI*x/17),0));
assert.ok(energies.every((e,j)=>Math.abs(e-(17+(-1)**(j+1)*sq))<1e-12));
console.log('chord-energy vector',JSON.stringify(energies),'cos',W/2,'side',Math.sqrt(2-W));
let M=Array.from({length:17},(_,i)=>Array.from({length:17},(_,j)=>+C(1,0).includes((j-i+17)%17)));
for(let i=0;i<17;i++)for(let j=0;j<17;j++){
 const m2=M[i].reduce((s,x,k)=>s+x*M[k][j],0);
 assert.equal(m2+M[i][j]-4*(i===j),4);
}
console.log('Paley adjacency identity checked in all 289 entries');
