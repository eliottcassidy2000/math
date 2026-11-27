// Exact controls; run with Node.js.
const assert=require("node:assert/strict");
const range=n=>Array.from({length:n},(_,i)=>i), bit=i=>1n<<BigInt(i);
const full=n=>bit(n)-1n;
function pop(x){let w=0;while(x){x&=x-1n;w++;}return w;}
const mask=(n,f)=>range(n).reduce((x,i)=>f(i)?x|bit(i):x,0n);
function paley(q){
  const Q=new Set(range(q-1).map(i=>(i+1)**2%q));
  return range(q+1).map(i=>range(q+1).map(j=>i===j?0:i===q?1:j===q?-1:
    Q.has((j-i+q)%q)?1:-1));
}
function double(A){
  const n=A.length;
  return range(2*n).map(i=>range(2*n).map(j=>{
    const v=A[i%n][j%n],e=+(i%n===j%n);
    return i<n?(j<n?v:v+e):(j<n?v-e:-v);
  }));
}
function compare(A,r){
  const n=A.length, ones=full(n);
  const normal=x=>x&bit(n-1)?x^ones:x;
  const plus=range(n).map(i=>mask(n,j=>A[i][j]<0));
  const minus=range(n).map(i=>mask(n,j=>A[i][j]>0));
  const C=plus.map(x=>normal(x^plus[r])),D=minus.map(x=>normal(x^minus[r]));
  const basis=[];
  function add(x,c){
    for(let i=n-2;i>=0;i--)if(x&bit(i)){
      if(!basis[i]){basis[i]=[x,c];return;}
      x^=basis[i][0];c^=basis[i][1];
    }
    assert.equal(c,0n);
  }
  for(const x of C)add(x,x);
  assert.equal(basis.filter(Boolean).length,n/2-1);
  for(const x of D)add(x,0n);
  assert.equal(basis.filter(Boolean).length,n-2);
  function proj(x){
    x=normal(x);let c=0n;
    for(let i=n-2;i>=0;i--)if(x&bit(i)){
      assert(basis[i]);x^=basis[i][0];c^=basis[i][1];
    }
    assert.equal(x,0n);return c;
  }
  const e=range(n).map(i=>i===r?1:A[r][i]);
  for(let i=0;i<n;i++){
    let c=proj(bit(i)^bit(r));if(c&bit(i))c^=ones;
    assert.equal(c,mask(n,j=>A[i][j]*e[i]*e[j]<0));
    for(let j=i+1;j<n;j++)assert.equal(pop(proj(bit(i)^bit(j))),n/2);
  }
  function words(L){
    const out=new Set([0n,ones]);
    for(const a of L)for(const b of [...out])out.add(a^b);
    return [...out];
  }
  return n===24?[words(C),words(D)]:null;
}
let pair;
for(const A of [paley(7),double(paley(7)),paley(23),paley(31),paley(47)]){
  for(const r of [0,A.length-1])pair=compare(A,r)??pair;
  console.log("order",A.length,"all projected weights and both rooted reconstructions passed");
}
const [C,D]=pair;
for(const L of pair){
  const w={};for(const x of L)w[pop(x)]=(w[pop(x)]??0)+1;
  assert.deepEqual(w,{0:1,8:759,12:2576,16:759,24:1});
}
const dSet=new Set(D);
assert.deepEqual(C.filter(x=>dSet.has(x)),[0n,full(24)]);
const Cs=C.filter(x=>pop(x)===8), Ds=D.filter(x=>pop(x)===8);
const target=[13,96,228,256,130,32,4,0,0];
for(const c of Cs){const h=Array(9).fill(0);
  for(const d of Ds)h[pop(c&d)]++;
  assert.deepEqual(h,target);
}
console.log("opposite octads intersection vector",target.join(","));
const mod=x=>(x%23+23)%23, inv=x=>range(23).find(y=>mod(x*y)===1);
const trans=range(24).map(x=>x===23?23:mod(x+1));
const invert=range(24).map(x=>x===23?0:x===0?23:mod(-inv(x)));
const neg=range(24).map(x=>x===23?23:mod(-x));
const perm=(x,p)=>mask(24,i=>!!(x&bit(p.indexOf(i))));
const cSet=new Set(C);
for(const p of [trans,invert])assert(C.every(x=>cSet.has(perm(x,p)))&&D.every(x=>dSet.has(perm(x,p))));
assert(C.every(x=>dSet.has(perm(x,neg)))&&D.every(x=>cSet.has(perm(x,neg))));
let seen=new Set(),queue=[range(24)];
seen.add(queue[0].join());
for(let j=0;j<queue.length;j++)for(const g of [trans,invert]){
  const p=queue[j].map(x=>g[x]),k=p.join();
  if(!seen.has(k)){seen.add(k);queue.push(p);}
}
assert.equal(seen.size,6072);
console.log("generated ordered group",seen.size,"unordered",2*seen.size,
            "orbit",244823040/seen.size);
const short=[];
for(let i=0;i<24;i++)for(let j=i+1;j<24;j++)for(const a of [-4,4])for(const b of [-4,4])
  short.push(range(24).map(k=>k===i?a:k===j?b:0));
for(let i=0;i<24;i++)for(const a of [-1,1])short.push(range(24).map(k=>k===i?-3*a:a));
assert.equal(short.length,1152);
const counts={};
for(const v of short){
  let d=0;for(const w of short){
    let p=0;for(let i=0;i<24;i++)p+=v[i]*w[i];
    if(Math.abs(p)===16)d++;
  }
  counts[d]=(counts[d]??0)+1;
}
assert.deepEqual(counts,{92:48,176:552,180:552});
console.log("common shortest vectors, absolute-inner-product-two degrees",JSON.stringify(counts));
