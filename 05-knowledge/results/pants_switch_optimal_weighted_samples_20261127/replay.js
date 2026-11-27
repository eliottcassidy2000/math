// Run with Node.js. No floating-point computations are needed.
const assert=require("node:assert/strict");
const pop = x => x===0?0:1+pop(x&(x-1));
const chi = (S,x) => pop(S&x)&1?-1:1;
function nonlinear(n, h, expected) {
  // Representatives with last sign +1.
  const N=1<<(n-1), full=(1<<n)-1;
  const g=Array.from({length:N},(_,m)=>chi(full,m)*h(m));
  assert(g.every(v=>v===1||v===-1));
  const size=g.filter(v=>v===1).length;
  assert.equal(size,N/2);
  const spectrum=[];
  for(let S=0;S<=full;S++){
    const w=pop(S); if(w&1)continue;
    const W=g.reduce((s,v,m)=>s+v*chi(S,m),0);
    if(!W)continue;
    assert(w>=2*Math.floor(n/3)+2);
    assert.equal(2*Math.abs(W),N);
    spectrum.push(w);
  }
  spectrum.sort((a,b)=>a-b);
  assert.deepEqual(spectrum,expected);
  console.log("order",n,"nonlinear sample",size,"Fourier weights",spectrum.join(","));
}
nonlinear(8,m=>(chi(3,m)+chi(5,m)+chi(10,m)-chi(12,m))/2,[6,6,6,6]);
nonlinear(11,m=>(chi(1,m)+chi(2,m)+chi(4,m)-chi(7,m))/2,[8,10,10,10]);
for(let n=2;n<=128;n++){
  const d=2*Math.floor(n/3)+2;
  assert(2*(n-d)<=d-2);
  if(n!==3)assert(d<=n-n%2);
  if(n>=3)assert(4*Math.floor(n/3)>=n-n%2);
}
for(let n=3;n<=64;n++){
  const m=2**Math.ceil(Math.log2(n));
  // Column index 0 is constant; columns 1, 2, 4, ... generate distinct rows.
  const rows=Array.from({length:m},(_,x)=>Array.from({length:n},(_,j)=>chi(j,x)));
  assert.equal(new Set(rows.map(x=>x.join(","))).size,m);
  for(let j=0;j<n;j++)for(let k=0;k<n;k++)
    assert.equal(rows.reduce((s,v)=>s+v[j]*v[k],0),j===k?m:0);
}
console.log("Sylvester design: verified every order 3 through 64");
