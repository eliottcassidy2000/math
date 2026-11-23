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

const {reps,framed9}=deck(91,9,new Set([3,11,13,39,91]));
assert(framed9.length===182,"order-nine 91 framed bank");
for(const a of framed9)if(!reps.get(91).some(b=>b.length===9&&iso(a,b)))reps.get(91).push(a);
assert(String(reps.get(91).map(r=>r.length).reduce((m,x)=>(m[x]=(m[x]||0)+1,m),[]))===",,,,,,,15,12,2","91 core classes");
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


const t13=[[13]],t39=[[39],[3,13],[13,3]],t91=[[91]];
assert(reps.get(13).length===1&&reps.get(13)[0].length===5,"13 core bank");
assert(eq(reps.get(39).map(r=>r.length),[7,7]),"39 core bank");
assert(eq(run(8,t13,[169]),[4,0,0]),"13 to 169");
assert(eq(run(7,t13,[91]),[3,0,0]),"13 to 91");
assert(eq(run(8,t39,[117]),[6,0,0]),"short orbit at order eight");
assert(eq(run(9,t39,[117]),[12,2,2]),"short orbit at order nine");
const core143W=[
 [2,4,201,147,39,207,155,19],[2,196,9,19,231,207,141,13],
 [322,4,393,499,487,455,390,259,2],[450,4,393,499,487,391,422,258,2]];
const t = sourceDeck(t91,10), arrived=[], arrivals=[0,0,0];
let no169=0;
for(const a of t){
 const r=join(a);
 for(const [i,j] of edges(10)){
  flip(r,i,j);
  if(hits(r,143)){
    const c=sccRows(r).filter(c=>c.length>1);assert(c.length===1,"direct 143 core");
    const v=c[0].map(i=>c[0].reduce((z,j,k)=>z|((r[i]>>j&1)<<k),0));
    assert(v.length===a.find(b=>b.h===91).n,"direct predecessor same order");
    arrivals[v.length-7]++;
    if(!arrived.some(x=>x.length===v.length&&iso(x,v)))arrived.push(v);
    for(const [k,l] of edges(10)){
     flip(r,k,l);no169+=hits(r,169);flip(r,k,l);
    }
  }
  flip(r,i,j);
 }
}
assert(t.length===100&&eq(arrivals,[0,6,4])&&no169===0&&arrived.length===4,"middle chain census");
assert(core143W.every(a=>arrived.some(b=>a.length===b.length&&iso(a,b))),"complete 143 precursor bank");
console.log({t:t.length,arrivals,endpoint143:arrived.length,no169});
function mass(r){return sccRows(r).reduce((s,c)=>s+(c.length>1?c.length:0),0);}
function minMass(r){
 let m=r.length;for(const [i,j] of edges(r.length)){
  flip(r,i,j);m=Math.min(m,mass(r));flip(r,i,j);
 }
 return m;
}
assert(core143W.every(r=>minMass(r)===6)&&reps.get(13).every(r=>minMass(r)===4)&&
 reps.get(11).every(r=>minMass(r)===4),"one-edit surviving support");
const r65=[18,20,9,19,36,79,31];
let count65=0;for(let m=0;m<1<<10;m++)count65+=height(rows(6,m))===65;
assert(height(r65)===65&&count65===0,"fixed orbit minimum");
const r39=[430,508,504,496,481,448,385,256,0];
assert(height(r39)===39,"short witness source");
flip(r39,6,8);assert(height(r39)===117,"short witness target");
const long=[234,244,249,242,225,192,128,0],orig=long.slice(), trace=[height(long)];
for(const [i,j] of [[3,7],[4,6],[1,3],[1,4],[1,3],[4,6],[1,4],[3,7]]){
 flip(long,i,j);trace.push(height(long));
}
assert(eq(trace,[13,85,169,99,143,203,91,85,13])&&eq(long,orig),"long literal witness");
console.log({trace});

const e7=edges(7);
let checks2=0,checks3=0;
for(const s of sourceDeck(t13,7)){
 const r=join(s);
 for(let a=0;a<e7.length;a++){
  flip(r,...e7[a]);
  for(let b=a+1;b<e7.length;b++){
   flip(r,...e7[b]);checks2++;assert(!hits(r,169),"order7 pair");
   for(let c=b+1;c<e7.length;c++){
    flip(r,...e7[c]);checks3++;assert(!hits(r,169),"order7 triple");flip(r,...e7[c]);
   }
   flip(r,...e7[b]);
  }
  flip(r,...e7[a]);
 }
}
assert(checks2===630&&checks3===3990,"order7 census totals");
for(let k=6;k<=8;k++){
 const x=2**k-2;
 for(const [a,tag,w,t] of [
  [[6,4,8,19,39,15],37,k-3,13*x+1],
  [[34,28,9,17,37,78,31],56,k-5,15*x+1]
 ]){
  const n=a.length+w+1,F=(1<<n)-1,l=w+1;
  const r=Array.from({length:n},(_,i)=>F^((1<<(i+1))-1));
  for(let i=0;i<a.length;i++){
   r[l+i]=a[i]<<l;
   if(tag>>i&1)flip(r,0,l+i);
  }
  assert(height(r)===t,"literal lift factor");
 }
 console.log("lift",k,195*x+14,13*x+1,15*x+1);
}

const l13=withProfile([6,4,8,19,39,15]), l15=withProfile([34,28,9,17,37,78,31]);
function signedStats(r,mask){
 const n=r.length,F=(1<<n)-1;
 function en(a){
  const d=Array.from({length:F+1},()=>Array(n).fill(0));
  for(let i=0;i<n;i++)d[1<<i][i]=1;
  for(let s=1;s<=F;s++)for(let i=0;i<n;i++)if(d[s][i])
    for(let m=a[i]&(F^s);m;m&=m-1){const t=m&-m;d[s|t][idx(t)]+=d[s][i];}
  return d;
 }
 const e=en(r),st=en(r.map((x,i)=>F^x^(1<<i)));
 const pick=a=>a.reduce((z,x,i)=>z+(mask>>i&1)*x,0);
 let G=0,L=0;
 for(let s=1;s<=F;s++){
  const t=F^s,u=pick(e[s]);G+=u*(t?st[t].reduce((a,b)=>a+b,0):1);L+=u*pick(st[t]);
 }
 return [height(r),G,L];
}
assert(eq(signedStats(l13.r,37),[15,104,40])&&eq(signedStats(l15.r,56),[89,480,118]),"affine lift coefficients");
