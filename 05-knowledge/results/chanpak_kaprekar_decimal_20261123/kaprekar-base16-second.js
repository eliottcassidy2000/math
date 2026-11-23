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
const {reps,framed9}=deck(75,9,new Set([3,5,9,15,25,45,75]));
assert(framed9.length===225,"order-nine 75 frames");
for(const [h,rs] of reps)console.log("cores",h,rs.map(r=>r.length));

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
const l45=[[45],[3,15],[15,3],[5,9],[9,5],[3,3,5],[3,5,3],[5,3,3]];
const l75=[[75],[3,25],[25,3],[5,15],[15,5],[3,5,5],[5,3,5],[5,5,3]];
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
assert(eq(run(17,l45,[165,105]),[1069,160,859,0]),"incompatible fork");
assert(eq(run(13,l75,[105]),[271,0,0]),"75 cross and small internal");
assert(eq(run(9,l45,[165]),[53,0,0]),"order-nine first jump");
let hits9=0,checks9=0;
for(const r of framed9)for(let i=0;i<9;i++)for(let j=i+1;j<9;j++){
  flip(r,i,j);checks9++;hits9+=height(r)===105;flip(r,i,j);
}
assert(checks9===8100&&hits9===0,"order-nine internal");
console.log({checks9,hits9});
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
const s7=sourceDeck(l45,7);let pairs=0,matches=0;
for(const s of s7) {
  const r=join(s), es=edges(7);
  assert(eq(fastNeighbors(s).sort(),es.map(([i,j])=>{
    flip(r,i,j);const h=height(r);flip(r,i,j);return h;
  }).sort()),"fast gain formula");
  for(let i=0;i<es.length;i++) {
    flip(r,...es[i]);
    assert(height(r)!==165,"order-seven first jump");
    for(let j=i+1;j<es.length;j++){flip(r,...es[j]);pairs++;matches+=height(r)===165;flip(r,...es[j]);}
    flip(r,...es[i]);
  }
}
assert(s7.length===9&&pairs===1890&&matches===0,"order-seven two-flip");
const row0=[250,252,249,112,160,192,144,8], row=row0.slice(),
  trace=[height(row)];
for(const [u,v] of [[0,3],[1,7],[1,2],[0,3],[1,2],[1,7]]){
  flip(row,u,v);trace.push(height(row));
}
assert(eq(trace,[45,87,165,75,33,105,45])&&eq(row,row0),"literal loop");
const solitary=[1006,1020,1016,1008,993,960,384,768,512,64];
assert(height(solitary)===45,"solitary source");
flip(solitary,1,5);assert(height(solitary)===165,"solitary target");
console.log({pairs,matches,trace});

function shellStats(r,tag) {
  const n=r.length, F=(1<<n)-1;
  function dp(rows) {
    const d=Array.from({length:F+1},()=>Array(n).fill(0));
    for(let i=0;i<n;i++)d[1<<i][i]=1;
    for(let s=1;s<=F;s++) for(let i=0;i<n;i++) if(d[s][i])
      for(let j=0;j<n;j++)if((rows[i]>>j&1)&&!(s>>j&1))d[s|(1<<j)][j]+=d[s][i];
    return d;
  }
  const end=dp(r), start=dp(r.map((x,i)=>F^(1<<i)^x));
  let D=end[F].reduce((a,b)=>a+b,0), G=0,L=0;
  for(let s=1;s<=F;s++) {
    const a=end[s].reduce((z,v,i)=>z+((tag>>i&1)?v:0),0);
    const b=F^s,c=b?start[b].reduce((z,v)=>z+v,0):1;
    const e=start[b].reduce((z,v,i)=>z+((tag>>i&1)?v:0),0);
    G+=a*c; L+=a*e;
  }
  return {D,G,L};
}

for(const [a,tag,expect] of [
  [[70,20,8,19,37,79,30],6,[81,480,50]],
  [[6,12,8,17,39,79,31],2,[35,136,0]]
]){
  const {D,G,L}=shellStats(a,tag);
  assert(eq([D,G,L],expect),"signed shell coefficient");
  for(let w=1;w<=4;w++){
    const n=a.length+w+1,F=(1<<n)-1,start=w+1;
    const r=Array.from({length:n},(_,i)=>F^((1<<(i+1))-1));
    for(let i=0;i<a.length;i++){
      r[start+i]=a[i]<<start;
      if(tag>>i&1)flip(r,0,start+i);
    }
    assert(height(r)===G*2**w+D-L,"literal signed shell");
  }
}
console.log("lifted shell coefficients verified");
