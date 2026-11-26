
const fs=require("node:fs"),assert=require("node:assert/strict");
const {boards,catalogue}=require("./minimum_trades");
const fixture=fs.readFileSync(process.argv[2],"utf8");
const Q0=[], names=new Map();
const labels=Array.from({length:56},(_,i)=>i).filter(i=>i!==0&&i!==28);
const bit=Array.from({length:56},(_,i)=>1n<<BigInt(i));
const pop=x=>{let n=0;for(;x;x&=x-1)n++;return n;};
const S=(a)=>a.reduce((b,i)=>b|bit[i],0n);
catalogue(boards(fixture).full,(q,old,other)=>{
 const k=q.join();let r=names.get(k);
 if(!r){r={q, set:S(q), map:new Map(q.map(a=>[a,new Map()]))};names.set(k,r);Q0.push(r);}
 for(const a of q){
  const x=old.find(z=>z[0]===a)[1],m=r.map.get(a);
  assert(!m.has(x%28));m.set(x%28,x<28?old:other);
 }
});
const Q=Q0.filter(r=>r.map.get(r.q[0]).size>=3);
assert.equal(Q.length,148);
const triples=new Map();
for(let t=0;t<Q.length;t++)for(const j of Q[t].q){
 const key=S(Q[t].q.filter(i=>i!==j));
 if(!triples.has(key))triples.set(key,[]);
 triples.get(key).push(j);
}
function enumK4(){
 let count=0,ok=0;
 for(let a=0;a<54;a++)for(let b=a+1;b<54;b++)for(let c=b+1;c<54;c++)for(let d=c+1;d<54;d++){
  count++;const L=[labels[a],labels[b],labels[c],labels[d]],s=S(L);
  if(L.every(i=>(triples.get(s^bit[i])||[]).some(j=>!(s&bit[j]))))ok++;
 }
 return {count,ok};
}
function candidates(l,d){
 const out=new Set();
 function pad(s,p,n){
  if(n===l){out.add(s);return;}
  for(let j=p;j<labels.length;j++)if(!(s&bit[labels[j]]))pad(s|bit[labels[j]],j+1,n+1);
 }
 if(d===4){
  function get(s,n,st,k){
   if(!k){pad(s,0,n);return;}
   for(let i=st;i<Q.length;i++){
    let b=0,ns=s|Q[i].set;
    for(const v of Q[i].q)if(!(s&bit[v]))b++;
    if(n+b<=l)get(ns,n+b,i+1,k-1);
   }
  }
  get(0n,0,0,l-4);
 }else {
  assert(l===5&&d===8);
  for(const s of triples.keys())pad(s,0,3);
 }
 return out;
}
function enumGraphs(L,d,collect){
 const l=L.length, s=S(L);
 const rem=[];
 // Map each possible projection to its eligible whole quartets.
 for(let i=0;i<l;i++){
  const m=new Map();
  for(let k=0;k<Q.length;k++)if(!(Q[k].set&bit[L[i]])){
   const zz=Q[k].set&s;let z=0;
   for(let j=0;j<l;j++)if(zz&bit[L[j]])z|=1<<j;
   const need=4-pop(z);if(need>d)continue;
   if(!m.has(z))m.set(z,[]);
   m.get(z).push(k);
  }
  rem.push([...m].map(([z,q])=>({z,q,need:4-pop(z)})).sort((a,b)=>a.need-b.need));
 }
 const chosen=[];
 const minima=rem.map(r=>r[0]?.need??100), suf=[0];
 for(let i=l-1;i>=0;i--)suf.unshift(suf[0]+minima[i]);
 if(suf[0]>d)return;
 function go(i,loss,expect){
  if(i===l){collect(chosen);return;}
  const below=(1<<i)-1;
  for(const obj of rem[i]){
   if(loss+obj.need+suf[i+1]>d)break;
   if((obj.z&below)!==expect[i])continue;
   const n=expect.slice();
   for(let k=i+1;k<l;k++)if(obj.z>>k&1)n[k]|=1<<i;
   chosen.push(obj);go(i+1,loss+obj.need,n);chosen.pop();
  }
 }
 go(0,0,new Array(l).fill(0));
}
function hall(c){
 const l=c.length;
 for(let card=1;card<=l;card++)for(let m=1;m<(1<<l);m++)if(pop(m)===card){
  const q=new Set();let n=0;
  for(let j=0;j<l;j++)if((m>>j&1)&&c[j].z){n++;for(const k of c[j].q)q.add(k);}
  if(q.size<n)return card;
 }
 return 0;
}
const K4=enumK4();assert.deepEqual(K4,{count:316251,ok:0});
const first=[];
for(const l of [5,6,7]){
 const sets=candidates(l,4), found={l,sets:sets.size,withGraph:0,graphs:0,pass:0,witness:[0,0,0,0,0,0,0,0]};
 for(const s of sets){
  const L=labels.filter(j=>s&bit[j]);let hit=false;
  enumGraphs(L,4,c=>{hit=true;found.graphs++;const k=hall(c);found.witness[k]++;if(!k)found.pass++;});
  if(hit)found.withGraph++;
 }
 first.push(found);
}
assert.deepEqual(first.map(r=>[r.sets,r.withGraph,r.graphs,r.pass]),[[7362,15,15,0],[2169,19,20,0],[786,82,91,0]]);
assert.equal(first.reduce((n,r)=>n+r.witness[2],0),126);
assert.equal(first.reduce((n,r)=>n+r.witness[3],0),0);
console.log(JSON.stringify({K4,first},null,2));

function localAssignments(L,c,p){
 const n=L.length, selected=[], k=p/2, every=(1<<28)-1;
 let anyLocal=false, anyProp=false, anyExact=false, packets=0;
 const use=c.map(o=>o.q.filter(q=>Q[q].map.get(Q[q].q[0]).size>=k));
 const pairLab=(row,a)=>row.find(x=>x[0]===a)[1]%28;
 function shared(q,r,a,x){
  const A=Q[q].map.get(a).get(x), B=Q[r].map.get(a).get(x);
  return A&&B&&!A.some(t=>t[0]!==a&&B.some(u=>t[0]===u[0]&&t[1]===u[1]));
 }
 const masks=new Map(), co=new Map();
 function available(q,a){
  const m=64*q+a;if(!masks.has(m)){
   let z=0;for(const x of Q[q].map.get(a).keys())z|=1<<x;
   masks.set(m,z);
  }return masks.get(m);
 }
 function compat(q,r,a){
  const l=(q*Q.length+r)*64+a;
  if(!co.has(l)){
   let z=0;for(const x of Q[q].map.get(a).keys())if(shared(q,r,a,x))z|=1<<x;
   co.set(l,z);
  }return co.get(l);
 }
 function test(active){
  anyLocal=true;
  // A low's full active set maps bijectively onto the active set of every closed target.
  const routes=selected.flatMap((q,i)=>{
   const t=[];for(let a=0;a<n;a++)if(c[i].z>>a&1)t.push(a);
   return t.slice(1).map(b=>[q,t[0],b]);
  });
  let changed=true;
  while(changed){
   changed=false;
   for(const [q,a,b] of routes){
    let x1=0,y1=0;
    for(const [x,row] of Q[q].map.get(L[a])){
     const y=pairLab(row,L[b]);
     if((active[a]>>x&1)&&(active[b]>>y&1)){x1|=1<<x;y1|=1<<y;}
    }
    if(pop(x1)<k||pop(y1)<k)return;
    if(x1!==active[a]||y1!==active[b]){changed=true;active[a]=x1;active[b]=y1;}
   }
  }
  anyProp=true; packets++;
  // Exact mask synchronization: take whole connected components in the lift
  // of the active-label transport relations. Each layer must take exactly k.
  const parent=Array.from({length:28*n},(_,i)=>i);
  const rt=i=>{while(parent[i]!==i)i=parent[i];return i;};
  for(const [q,a,b] of routes)for(const [x,row] of Q[q].map.get(L[a]))if(active[a]>>x&1){
   const y=pairLab(row,L[b]);assert(active[b]>>y&1);
   parent[rt(28*a+x)]=rt(28*b+y);
  }
  const comp=new Map();
  for(let a=0;a<n;a++)for(let x=0;x<28;x++)if(active[a]>>x&1){
   const u=rt(28*a+x);
   if(!comp.has(u))comp.set(u,new Array(n).fill(0));
   comp.get(u)[a]++;
  }
  const K=k+1, pow=Array.from({length:n},(_,i)=>K**i);
  const target=k*pow.reduce((a,b)=>a+b,0);
  let dp=new Set([0]);
  for(const row of comp.values()){
   if(row.some(x=>x>k))continue;
   const add=row.reduce((t,x,j)=>t+x*pow[j],0), next=new Set(dp);
   for(const x of dp)if(row.every((t,j)=>Math.floor(x/pow[j])%K+t<=k))next.add(x+add);
   dp=next;
  }
  if(dp.has(target))anyExact=true;
 }
 function go(i,old){
  if(i===n){test(old);return;}
  for(const q of use[i]){
   if(selected.includes(q))continue;
   const active=old.slice();
   let ok=true;
   for(let a=0;a<n;a++)if(c[i].z>>a&1){
    active[a]&=available(q,L[a]);
    for(let j=0;j<i;j++)if(c[j].z>>a&1)active[a]&=compat(q,selected[j],L[a]);
    if(pop(active[a])<k){ok=false;break;}
   }
   if(ok){selected.push(q);go(i+1,active);selected.pop();}
  }
 }
 go(0,new Array(n).fill(every));
 return [anyLocal,anyProp,anyExact,packets];
}
if(process.argv.includes("--eight")){
 const c=candidates(5,8);assert.equal(c.size,605011);
 let ng=0,pass=0,labelsPassed=0;
 const detail={8:[0,0,0,0],10:[0,0,0,0]}, survivors={8:[],10:[]};
 for(const s of c){
  const L=labels.filter(i=>s&bit[i]);let here=false;
  enumGraphs(L,8,cs=>{
   ng++;if(hall(cs))return;
   assert(cs.every(x=>x.z));
   here=true;pass++;
   for(const p of [8,10]){
    const row=localAssignments(L,cs,p);
    for(let i=0;i<4;i++)detail[p][i]+=row[i];
    if(row[2])survivors[p].push([L,cs.map(x=>x.z)]);
   }
  });labelsPassed+=here;
 }
 assert.deepEqual([ng,pass,labelsPassed],[257570,19755,8673]);
 assert.deepEqual(detail[8].slice(0,2),[48,45]);assert.equal(detail[10][1],45);
 assert.deepEqual([detail[8][2],detail[10][2]],[0,0]);
 console.log(JSON.stringify({ng,pass,labelsPassed,detail,survivors}));
}
