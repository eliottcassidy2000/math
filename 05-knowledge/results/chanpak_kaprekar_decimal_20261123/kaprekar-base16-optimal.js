"use strict";
const idx = b => 31-Math.clz32(b);
const assert = (ok,m) => {if(!ok)throw Error(m);};
function height(r) {
  const n=r.length,F=(1<<n)-1,d=new Uint32Array((F+1)*n);
  for(let i=0;i<n;i++)d[(1<<i)*n+i]=1;
  for(let u=1;u<=F;u++)for(let a=u;a;a&=a-1) {
    const i=idx(a&-a), p=d[u*n+i];
    if(p)for(let b=r[i]&(F^u);b;b&=b-1){
      const j=b&-b; d[(u|j)*n+idx(j)]+=p;
    }
  }
  let H=0;
  for(let i=0;i<n;i++)H+=d[F*n+i];
  return H;
}
function mark(n){
  const es=[],bits=Array.from({length:n},()=>Array(n).fill(-1));
  for(let j=2;j<n;j++)for(let i=0;i<j-1;i++){
    bits[i][j]=es.length; es.push([i,j]);
  }
  return {es,bits};
}
function rows(n,m,es=mark(n).es){
  const r=Array.from({length:n},(_,i)=>i<n-1?1<<(i+1):0);
  for(let k=0;k<es.length;k++){
    const [a,b]=es[k];
    if(m>>k&1)r[a]|=1<<b; else r[b]|=1<<a;
  }
  return r;
}
function strong(r){
  let s=1,old=0;
  while(old!==s){
    old=s;
    for(let i=1;i<r.length;i++)if(r[i]&s)s|=1<<i;
  }
  return s===(1<<r.length)-1;
}
function flip(r,i,j){r[i]^=1<<j;r[j]^=1<<i;}
function iso(a,b){
  const n=a.length,pop=x=>x?1+pop(x&(x-1)):0,da=a.map(pop),db=b.map(pop),p=[];
  function rec(i,used){
    if(i===n)return true;
    for(let j=0;j<n;j++)if(!(used>>j&1)&&da[i]===db[j]&&
      p.every((k,l)=>(a[i]>>l&1)===(b[j]>>k&1))){
      p.push(j);if(rec(i+1,used|(1<<j)))return true;p.pop();
    }
    return false;
  }
  return rec(0,0);
}
// Prefix/deletion pruning from monotonicity under insertion.
function deck(K,maxN,wanted) {
  const reps=new Map(), framed9=[];
  let prev=new Set([0]);
  for(let n=2;n<=maxN;n++){
    const {es,bits}=mark(n), old=mark(n-1).es, next=new Set();
    const maps=Array.from({length:n-1},(_,v)=>{
      const rem=Array.from({length:n},(_,i)=>i).filter(i=>i!==v);
      return {
        sep:v>0?bits[v-1][v+1]:-1,
        select:old.map(([i,j])=>bits[rem[i]][rem[j]])
      };
    });
    const late=1<<old.length, newArcs=1<<Math.max(0,n-2);
    let sc=0;
    for(const p of prev) outer:for(let x=0;x<newArcs;x++){
      const m=p+x*late;
      for(const {sep,select} of maps){
        if(sep!==-1&&!(m>>sep&1))continue; // deletion not path-framed
        let code=0;
        for(let j=0;j<select.length;j++)code|=((m>>select[j])&1)<<j;
        if(!prev.has(code))continue outer;
      }
      const r=rows(n,m,es), s=strong(r);
      if(s)sc++;
      if(n===9&&!s)continue;
      const h=height(r);
      if(h>K)continue;
      next.add(m);
      if(s&&wanted.has(h)){
        if(n===9)framed9.push(r);
        else {
          if(!reps.has(h))reps.set(h,[]);
          const out=reps.get(h);
          if(!out.some(a=>a.length===n&&iso(a,r)))out.push(r);
        }
      }
    }
    console.log("bank",K,n,next.size,sc);
    if(K===75&&n===7)assert(next.size===8033,"order-seven prefix");
    if(K===75&&n===8)assert(next.size===18497,"order-eight prefix");
    if(K===75&&n===9)assert(sc===13100,"order-nine strong candidates");
    prev=next;
  }
  return {reps,framed9};
}

// Only the height-15 start deck is needed. The complete strong bound limits it to n<=6.
const {reps}=deck(15,6,new Set([3,5,15]));
const c15=reps.get(15);
assert(c15.filter(x=>x.length===5).length===2&&c15.filter(x=>x.length===6).length===2,"height-15 core deck");
const tri=reps.get(3)[0], five=reps.get(5)[0];
const starts=[...c15.map(c=>[c]),[tri,five],[five,tri]];
function padded(n,h=15){
  const out=[];
  function go(cs,i,k,path){
    if(i===cs.length){
      const blocks=path.concat(Array(k).fill([0])), r=[];
      let l=0;
      for(const b of blocks){
        const tail=((1<<n)-1)^((1<<(l+b.length))-1);
        for(const w of b)r.push((w<<l)|tail);
        l+=b.length;
      }
      out.push(r);return;
    }
    for(let a=0;a<=k;a++)go(cs,i+1,k-a,path.concat(Array(a).fill([0]),[cs[i]]));
  }
  for(const cs of (h===15?starts:h===5?[[five]]:h===3?[[tri]]:[[]])){
    const k=n-cs.reduce((s,c)=>s+c.length,0);
    if(k>=0)go(cs,0,k,[]);
  }
  return out;
}
const exactMin=[1,1,Infinity,3,5,9,15,25,45,75,125,225,375];
function hits(r,target){
  const n=r.length,S=r.map((row,i)=>row|(1<<i)), C=[];
  for(let j=0;j<n;j++)for(let i=0;i<n;i++)if(S[i]>>j&1)S[i]|=S[j];
  let todo=(1<<n)-1, min=1;
  while(todo){
    const i=idx(todo&-todo), c=[];
    for(let a=todo&S[i];a;a&=a-1){
      const b=a&-a,j=idx(b);if(S[j]>>i&1){c.push(j);todo^=b;}
    }
    min*=exactMin[c.length];if(min>target)return false;
    if(c.length>1)C.push(c);
  }
  for(const c of C){
    const sub=c.map(i=>c.reduce((m,j,k)=>m|((r[i]>>j&1)<<k),0));
    const n=height(sub);
    if(target%n!==0)return false;
    target/=n;
  }
  return target===1;
}
function stages(n, lengths, targets,h=15){
  const rows=padded(n,h), es=[];
  for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)es.push([i,j]);
  const arrivals=targets.map(()=>0), tests=targets.map(()=>0);
  function visit(r,k){
    if(k===targets.length)return;
    for(let i=0;i<es.length;i++){
      flip(r,...es[i]);
      if(lengths[k]===1)check();
      else for(let j=i+1;j<es.length;j++){
        flip(r,...es[j]);check();flip(r,...es[j]);
      }
      flip(r,...es[i]);
    }
    function check(){
      tests[k]++;
      if(hits(r,targets[k])){arrivals[k]++;visit(r,k+1);}
    }
  }
  for(const r of rows){
    assert(hits(r,h),"source");
    visit(r,0);
  }
  console.log({n,lengths,targets,sources:rows.length,tests,arrivals});
  return [rows.length,...tests,...arrivals].join();
}
assert(stages(12,[1,1],[135,195])==="72,4752,11484,174,0","all-order suffix");
assert(stages(9,[2,1,2,1],[225,195,135,15])==="30,18900,1584,1260,144,44,2,4,0","first schedule");
assert(stages(9,[2,1,1,2],[225,195,135,15])==="30,18900,1584,72,0,44,2,0,0","second schedule");

for(const [n,h,u,v,check] of [
  [10,15,135,195,"42,1890,1530,34,0"],
  [8,5,45,65,"5,140,56,2,0"],
  [7,3,27,39,"5,105,210,10,0"],
  [5,1,9,13,"1,10,10,1,0"],
]){
  assert(stages(n,[1,1],[u,v],h)===check,"focused suffix");
}
