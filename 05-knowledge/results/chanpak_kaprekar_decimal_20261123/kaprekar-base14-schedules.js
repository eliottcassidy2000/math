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

const {reps,framed9}=deck(13,5,new Set([11,13]));
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
const exactMin=[1,1,Infinity,3,5,9,15,25,45,75,125,225,375,625,1125,1875,3125,5625,9375];
const cacheH=new Map();
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
    const k=sub.join();let n=cacheH.get(k);if(n===undefined){n=height(sub);if(cacheH.size<300000)cacheH.set(k,n);}
    if(target%n!==0)return false;
    target/=n;
  }
  return target===1;
}



const core143=[
 [2,4,201,147,39,207,155,19],[2,196,9,19,231,207,141,13],
 [322,4,393,499,487,455,390,259,2],[450,4,393,499,487,391,422,258,2]];
cached[143]=core143.map(withProfile); // cache object holds only n<=9 cores
const d13=[[13]],d143=[[143]], split143=[[11,13],[13,11]];
function pairLoop(r,es,accept,interior=null){
 let count=0;
 for(let i=0;i<es.length;i++){
  flip(r,...es[i]);
  for(let j=i+1;j<es.length;j++){
   if(interior&&!interior[i]&&!interior[j])continue;
   flip(r,...es[j]);count++;accept(r);flip(r,...es[j]);
  }
  flip(r,...es[i]);
 }
 return count;
}
function prefixA(){
 const source=sourceDeck(d13,15),es=edges(15),p=[0,0,0];let trials=0,next=0;
 for(const s of source)trials+=pairLoop(join(s),es,r=>{
  if(!hits(r,169))return;
  const comp=sccRows(r).filter(a=>a.length>1);
  if(comp.length===2)p[2]++;
  else p[comp[0].length-8]++;
  for(const a of es){
   flip(r,...a);
   if(hits(r,143)){
    const f=sccRows(r).filter(c=>c.length>1).map(c=>
      height(c.map(i=>c.reduce((z,j,k)=>z|((r[i]>>j&1)<<k),0))));
    assert(eq(f.sort((x,y)=>x-y),[11,13]),"schedule A exit type");next++;
   }
   flip(r,...a);
  }
 });
 assert(eq([source.length,trials,...p,next],[11,60060,16,28,42,168]),"schedule A prefix census");
 console.log("A prefix",source.length,trials,p,next);
}
function suffixA(n,expected){
 const es=edges(n),source=sourceDeck(split143,n);let trials=0,good=0;
 for(const s of source)trials+=pairLoop(join(s),es,r=>{good+=hits(r,91);});
 assert(eq([source.length,trials,good],expected),"schedule A suffix");
 console.log("A suffix",n,source.length,trials,good);
}
function scheduleB(n,expected){
 const es=edges(n), source=sourceDeck(d143,n);
 const peaks=new Set();
 let first=0,second=0,good=0;
 for(const s of source)first+=pairLoop(join(s),es,r=>{if(hits(r,169))peaks.add(r.join());});
 for(const s of peaks){
  const r=s.split(",").map(Number), id=Array(n);
  sccRows(r).forEach((c,k)=>c.forEach(v=>id[v]=k));
  second+=pairLoop(r,es,T=>{good+=hits(T,13);},es.map(([i,j])=>id[i]===id[j]));
 }
 console.log("B",n,source.length,first,peaks.size,second,good);
 if(expected)assert(eq([source.length,first,peaks.size,second,good],expected),"schedule B census");
}
prefixA();suffixA(10,[2,1980,0]);suffixA(11,[6,8910,0]);
scheduleB(14,[26,106470,80,180540,0]);
scheduleB(18,[42,488376,128,531880,0]);

function subgraph(r,c){return c.map(i=>c.reduce((a,j,k)=>a|((r[i]>>j&1)<<k),0));}
for(const r of core143){
 let hits6=0;
 for(const a of edges(r.length)){
  flip(r,...a);
  const core=sccRows(r).filter(x=>x.length>1);
  if(core.reduce((z,x)=>z+x.length,0)===6){
   assert(core.length===1&&height(subgraph(r,core[0]))===23,"mass-six local height");
   hits6++;
  }else assert(core.reduce((z,x)=>z+x.length,0)>6,"mass-six minimum");
  flip(r,...a);
 }
 assert(hits6===1,"unique mass-six edit");
}
for(const [i,edge,UU,VV,sinks] of [
 [2,[1,2],[29,57,119,169,108,35,34],[143,45,29,34,54,89,157],true],
 [3,[2,3],[89,29,45,157,143,54,34],[35,119,57,34,29,108,169],false]
]){
 const r=core143[i];
 let count=0; 
 for(const a of edges(9)){
  flip(r,...a);const comp=sccRows(r).filter(x=>x.length>1);
  if(comp.reduce((z,x)=>z+x.length,0)===7)count++;
  flip(r,...a);
 }
 assert(count===1,"unique order-nine mass-seven edit");
 flip(r,...edge);
 const comp=sccRows(r),c=comp.find(x=>x.length===7);
 assert(c&&comp.length===3,"order-seven block");
 const s=withProfile(subgraph(r,c));
 assert(s.h===29&&eq(s.U,UU)&&eq(s.V,VV),"rooted 29");
 for(const t of comp.filter(x=>x.length===1))assert(c.every(v=>Boolean(r[v]>>t[0]&1)===sinks),"released pads same side");
 flip(r,...edge);
}
scheduleB(11,[14,20790,44,51120,0]);
