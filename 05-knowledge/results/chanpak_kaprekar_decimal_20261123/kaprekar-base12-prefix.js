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
const {reps}=deck(55,8,new Set([3,5,9,11,33,55]));
const eq=(a,b)=>String(a)===String(b);
function withProfile(r){
  const n=r.length,F=(1<<n)-1;
  function dp(a){
    const d=Array.from({length:F+1},()=>Array(n).fill(0));
    for(let i=0;i<n;i++)d[1<<i][i]=1;
    for(let s=1;s<=F;s++)for(let i=0;i<n;i++)if(d[s][i])
      for(let m=a[i]&(F^s);m;m&=m-1){
        const b=m&-m;d[s|b][idx(b)]+=d[s][i];
      }
    return d;
  }
  const e=dp(r),s=dp(r.map((row,i)=>F^row^(1<<i)));
  const H=e.map(a=>a.reduce((u,v)=>u+v,0));H[0]=1;
  const U=Array(n).fill(0),V=Array(n).fill(0);
  let q=0;
  for(let m=0;m<=F;m++){
    const a=H[F^m];q+=a*H[m];
    for(let i=0;i<n;i++){ U[i]+=a*s[m][i]; V[i]+=a*e[m][i]; }
  }
  const internal=[];
  for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){
    const a=r.slice();flip(a,i,j);internal.push(height(a));
  }
  return {r,n,h:H[F],q,U,V,internal,start:s[F],end:e[F]};
}
const one=withProfile([0]);
const cached=Object.fromEntries([...reps].map(([h,rr])=>[h,rr.map(withProfile)]));
function sourceDeck(lists,n){
  let out=[];
  for(const list of lists)for(const template of select(list,[])){
    const left=n-template.reduce((v,c)=>v+c.n,0);
    if(left>=0)pad(template,0,left,[]);
  }
  return out;
  function* select(list,base){
    if(!list.length){yield base;return;}
    for(const c of cached[list[0]])yield* select(list.slice(1),base.concat(c));
  }
  function pad(t,i,k,b){
    if(i===t.length){out.push(b.concat(Array(k).fill(one)));return;}
    for(let l=0;l<=k;l++)pad(t,i+1,k-l,b.concat(Array(l).fill(one),[t[i]]));
  }
}
function fastNeighbors(bs){
  const H=bs.reduce((p,b)=>p*b.h,1);
  const ans=[];
  for(let a=0;a<bs.length;a++){
    const A=bs[a];
    for(const h of A.internal)ans.push(h*(H/A.h));
    let q=1,inner=1;
    for(let b=a+1;b<bs.length;b++){
      const B=bs[b],R=H/(A.h*B.h*inner);
      for(let u=0;u<A.n;u++)for(let v=0;v<B.n;v++){
        const gained=A.U[u]*q*B.V[v],
         lost=b===a+1?A.end[u]*B.start[v]:0;
        ans.push(H+R*(gained-lost));
      }
      q*=B.q;inner*=B.h;
    }
  }
  return ans;
}
function run(n,lists,goals){
  const shapes=sourceDeck(lists,n), counts=goals.map(()=>0);
  let both=0;
  for(const s of shapes){
    const h=fastNeighbors(s);
    const c=goals.map((v,i)=>h.filter(x=>x===v).length);
    c.forEach((v,i)=>counts[i]+=v);
    both+=+c.every(v=>v>0);
  }
  const r=[shapes.length,...counts,both];console.log(n,goals,r);
  return r;
}
function join(bs){
  const N=bs.reduce((s,b)=>s+b.n,0),F=(1<<N)-1,r=[];
  let l=0;
  for(const b of bs){
    const tail=F^((1<<(l+b.n))-1);
    for(const w of b.r)r.push((w<<l)|tail);
    l+=b.n;
  }
  return r;
}
function edges(n){
  const E=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)E.push([i,j]);
  return E;
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

const t11=[[11]], t33=[[33],[3,11],[11,3]], t55=[[55],[5,11],[11,5]];
assert(reps.get(11).length===1&&reps.get(11)[0].length===5,"11 core");
assert(eq(reps.get(33).map(r=>r.length),[6,6,6,6,6,7,7,7]),"33 cores");
assert(eq(reps.get(55).map(r=>r.length),[7,7,7,7]),"55 cores");
assert(eq(run(7,t11,[77]),[3,0,0]),"11 to 77 arc");
assert(eq(run(13,t33,[77,55]),[103,12,194,0]),"33 fork");
function pairCheck(n,src,goal){
  let visits=0;
  for (const s of sourceDeck([[src]],n)){
    const r=join(s), es=edges(n);
    for(let i=0;i<es.length;i++){
      flip(r,...es[i]);
      for(let j=i+1;j<es.length;j++){flip(r,...es[j]);visits+=hits(r,goal);flip(r,...es[j]);}
      flip(r,...es[i]);
    }
  }
  return visits;
}
assert(pairCheck(6,11,77)===0,"order-6 11 to 77");
assert(eq(run(9,[[33]],[77]),[29,4,3]),"strong 33");
for(const b of sourceDeck([[33]],9))if(fastNeighbors(b).includes(77))assert(b.some(s=>s.h===33&&s.n===7),"33 strong source size");
assert(eq(run(9,[[3,11],[11,3]],[77]),[6,0,0]),"reducible 33");
assert(eq(run(9,[[3]],[77]),[7,0,0]),"3 to 77");
{const a=run(9,[[5]],[77]);assert(a[0]===6&&a[1]===4,"5 to 77");}
assert(eq(run(9,[[9],[3,3]],[77]),[20,0,0]),"9 to 77");
assert(eq(run(9,t11,[77]),[5,0,0]),"11 to 77 at 9");
function sccRows(r) {
  const n=r.length,S=r.map((row,i)=>row|(1<<i)), C=[];
  for(let j=0;j<n;j++)for(let i=0;i<n;i++)if(S[i]>>j&1)S[i]|=S[j];
  let todo=(1<<n)-1;
  while(todo){
    const i=idx(todo&-todo), c=[];
    for(let a=todo&S[i];a;a&=a-1){
      const b=a&-a,j=idx(b);if(S[j]>>i&1){c.push(j);todo^=b;}
    }
    C.push(c);
  }
  return C;
}
const peaks10=sourceDeck(t55,10);
let n99=0,n121=0;
for(const bs of peaks10){
 const r=join(bs), es=edges(10);
 for(const f of es){
   flip(r,...f);
   if(hits(r,99)){
     n99++;
     for(const g of es){
       flip(r,...g);
       if(hits(r,121)){
         n121++;
         const non=sccRows(r).filter(s=>s.length>1);
         assert(non.length===2&&non.every(s=>s.length===5 &&
           height(s.map(i=>s.reduce((bits,j,k)=>bits|((r[i]>>j&1)<<k),0)))===11),"peak type");
       }
       flip(r,...g);
     }
   }
   flip(r,...f);
 }
}
assert(eq([peaks10.length,n99,n121],[22,8,12]),"peak census");
console.log({padded55:peaks10.length,n99,n121});

const direct77=[102,108,104,113,103,64,0];
assert(height(direct77)===11,"11 to 77 witness start");
flip(direct77,0,5);assert(height(direct77)===29,"11 to 77 witness middle");
flip(direct77,1,6);assert(height(direct77)===77,"11 to 77 witness end");
const H11=withProfile([18,4,9,19,6]);
assert(eq(H11.U,[20,12,13,27,17])&&eq(H11.V,[17,27,13,12,20])&&Math.min(...H11.internal)===5,"rooted 11 guards");

assert(cached[5][0].q===32 && eq(cached[3][0].U,[3,3,3]) &&
  eq(cached[3][0].V,[3,3,3]),"triangle and four-core cuts");
