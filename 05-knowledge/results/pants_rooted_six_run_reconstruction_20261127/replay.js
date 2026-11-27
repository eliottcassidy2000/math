// Exact recovery on deterministically switched bordered Paley specimens.
const assert=require("node:assert/strict");
const range=n=>Array.from({length:n},(_,i)=>i);
const pairs=a=>a.flatMap((x,i)=>a.slice(i+1).map(y=>[x,y]));
const triples=a=>a.flatMap((x,i)=>pairs(a.slice(i+1)).map(([y,z])=>[x,y,z]));
const mod8=x=>(x+800)%8;
function paley(q) {
  const sq=new Set(range(q-1).map(i=>(i+1)**2%q));
  return range(q+1).map(i=>range(q+1).map(j=>i===j?0:i===0?1:j===0?-1:
    sq.has((j-i+q)%q)?1:-1));
}
const cyclic=(B,a,b,c)=>B[a][b]===B[b][c]&&B[b][c]===B[c][a];
function decode(A) {
  const n=A.length,core=range(n).slice(1), b=1, rest=core.slice(1);
  const cache=new Map();
  function observe(i,j,k) {
    const key=[i,j,k].sort((a,b)=>a-b).join(",");
    if(!cache.has(key))cache.set(key,range(n).reduce((s,x)=>s+A[0][x]*A[i][x]*A[j][x]*A[k][x],0));
    return cache.get(key);
  }
  const flag=(i,j,k)=>mod8(observe(i,j,k))===mod8(n-4);
  const C=range(n).map(()=>range(n).map(()=>false));
  for(const [u,v] of pairs(rest)) C[u][v]=C[v][u]=flag(b,u,v);
  const colors=new Map([[rest[0],0]]), stack=[rest[0]];
  for(let p=0;p<stack.length;p++){
    const u=stack[p];
    for(const v of rest)if(C[u][v]){
      if(!colors.has(v)){colors.set(v,1-colors.get(u));stack.push(v);}
      else assert.notEqual(colors.get(u),colors.get(v));
    }
  }
  assert.equal(colors.size,rest.length);
  const O=rest.filter(v=>colors.get(v)===0), I=rest.filter(v=>colors.get(v)===1);
  const B=range(n).map(()=>range(n).map(()=>0));
  function set(u,v,x){B[u][v]=x;B[v][u]=-x;}
  for(const v of core)set(0,v,1);
  for(const v of O)set(b,v,1);
  for(const v of I)set(b,v,-1);
  for(const u of O)for(const v of I)set(u,v,C[u][v]?1:-1);
  for(const [same,other] of [[O,I],[I,O]])for(const [u,v] of pairs(same)){
    const w=other.find(w=>B[u][w]!==B[v][w]);
    assert(w!==undefined);
    const cycle=flag(u,v,w);
    // If u->w->v, a cycle occurs precisely when v->u.
    set(u,v,cycle?-B[u][w]:B[u][w]);
  }
  const stage1=cache.size, pivots=[];
  let rank=0;
  const must=(u,v,w)=>cyclic(B,u,v,w)===(n%8===0);
  function add(T) {
    const original=T.reduce((m,i)=>m|(1n<<BigInt(i)),0n);
    let m=original,p=n-1;
    for(;p>=1;p--)if(m>>BigInt(p)&1n){
      if(pivots[p])m^=pivots[p][0];else break;
    }
    if(p<1)return;
    const q=observe(...T), q0=range(n).reduce((s,x)=>s+B[0][x]*T.reduce((z,i)=>z*B[i][x],1),0);
    assert(q!==0&&q0!==0&&Math.abs(q)===Math.abs(q0));
    let tmp=original, rhs=+(q!==q0);
    for(let j=n-1;j>p;j--)if(tmp>>BigInt(j)&1n){tmp^=pivots[j][0];rhs^=pivots[j][1];}
    assert.equal(tmp,m);
    pivots[p]=[m,rhs];rank++;
  }
  for(const [u,v] of pairs(rest))if(must(b,u,v))add([b,u,v]);
  const starRank=rank;
  assert.equal(starRank,n-1-(n%8===0?2:1));
  const Ts=triples(core).filter(T=>must(...T));
  for(const T of Ts)if(rank<n-1)add(T);
  assert.equal(rank,n-1);
  const e=[1];
  let x=0n;
  for(let p=1;p<n;p++){
    const [m,rhs]=pivots[p];
    let z=rhs;
    for(let j=1;j<p;j++)if((m&x)>>BigInt(j)&1n)z^=1;
    if(z)x|=1n<<BigInt(p);
    e.push(z?-1:1);
  }
  const R=B.map((row,i)=>row.map((v,j)=>v*e[i]*e[j]));
  const polarity=R[0][1]*A[0][1];
  assert(R.every((row,i)=>row.every((v,j)=>v===polarity*A[i][j])));
  assert(cache.size<=stage1+(n%8===0?2:1));
  console.log("order",n,"mandatory rows",Ts.length,"star rank",starRank,
    "total distinct queries",cache.size,"reconstructed");
}
for(const n of [8,12,20,24,32,44]){
  const A=paley(n-1),e=range(n).map(i=>i%5===1||i%7===2?-1:1);
  for(let i=0;i<n;i++)for(let j=0;j<n;j++) {
    const g=A[i].reduce((s,v,x)=>s+v*A[j][x],0);
    assert.equal(g,i===j?n-1:0);
  }
  decode(A.map((row,i)=>row.map((v,j)=>v*e[i]*e[j])));
}
