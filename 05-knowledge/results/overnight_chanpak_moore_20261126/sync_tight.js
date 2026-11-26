"use strict";
const fs=require("node:fs"),assert=require("node:assert/strict");
const {boards,catalogue}=require("./minimum_trades");
const seed=fs.readFileSync(process.argv[2],"utf8"), Q=[], byName=new Map();
const pop=x=>{let n=0;for(;x;x&=x-1)n++;return n;};
function rec(q,old,newer){
 const name=q.join(",");
 let o=byName.get(name);
 if(!o){o={num:Q.length,q, map:new Map(q.map(a=>[a,new Map()]))};Q.push(o);byName.set(name,o);}
 for(const a of q){
   const a0=old.find(r=>r[0]===a)[1], map=o.map.get(a);
   const ori=a0<28?old:newer, x=a0%28;
   assert(!map.has(x));map.set(x,ori);
 }
}
const all=catalogue(boards(seed).full,rec);
assert.equal(all.total,1057);assert.equal(Q.length,152);
function mask(q,a){
 let m=0;for(const x of q.map.get(a).keys())m|=1<<x;return m;
}
function compat(q,r,a){
 let m=0;
 for(const [x,L] of q.map.get(a)){
  const R=r.map.get(a).get(x);if(!R)continue;
  if(!L.some(t=>t[0]!==a&&R.some(u=>t[0]===u[0]&&t[1]===u[1])))m|=1<<x;
 }
 return m;
}
const pattern = Array.from({length:57},()=>[]), eligible6=[];
for(let a=0;a<57;a++){
 const use=Q.filter(q=>q.q.includes(a)&&q.map.get(a).size>=6), n=use.length;
 const ms=use.map(q=>mask(q,a)), cp=use.map(q=>use.map(r=>compat(q,r,a)));
 let valid6=false;
 function walk(first,selected,active){
  if(selected.length===5){
   valid6=true;const ct=pop(active);
   if(ct>=7){assert.equal(ct,7);pattern[a].push({ids:selected.map(i=>use[i].num),active});}
   return;
  }
  for(let i=first;i<=n-(5-selected.length);i++){
   let m=active&ms[i];for(const x of selected)m&=cp[i][x];
   if(pop(m)>=6)walk(i+1,[...selected,i],m);
  }
 }walk(0,[],(1<<28)-1);
 if(valid6)eligible6.push(a);
}
assert.equal(eligible6.length,24);
assert.deepEqual(eligible6,pattern.map((v,i)=>[v,i]).filter(([v])=>v.length).map(([,i])=>i));
const sumPat=pattern.reduce((n,p)=>n+p.length,0);assert.equal(sumPat,573);
function incidence(S,t,J){
 const z=[];
 for(const q of Q)if(q.map.get(q.q[0]).size>=t){
  let s=S.filter(a=>q.q.includes(a)).length;if(s)z.push(s);
 }
 z.sort((a,b)=>b-a);return z.slice(0,J).reduce((s,a)=>s+a,0);
}
function image(q,a,b,active){
 let z=0;for(const [x,row] of q.map.get(a))
   if(active>>x&1)z|=1<<(row.find(r=>r[0]===b)[1]%28);
 return z;
}
function agree(P,T,a,b){
 const p=P.ids.filter(x=>Q[x].q.includes(b)), t=T.ids.filter(x=>Q[x].q.includes(a));
 return p.join()===t.join()&&p.every(x=>image(Q[x],a,b,P.active)===T.active);
}
let traversed=0,best=0,bestAt=0,pass=0,packed=0;
function combine(S,start){
 if(S.length===4){
  traversed++;
  const n=incidence(S,6,7);if(n>best){best=n;bestAt=1;}else if(n===best)bestAt++;
  if(incidence(S,7,9)<20)return;
  pass++;
  function fill(i,chosen,union){
   if(i===4){packed++;return;}
   for(const p of pattern[S[i]])if(chosen.every((q,j)=>agree(p,q,S[i],S[j]))){
    const u=new Set([...union,...p.ids]);if(u.size<=9)fill(i+1,[...chosen,p],u);
   }
  }fill(0,[],new Set());
  return;
 }
 for(let j=start;j<=eligible6.length-(4-S.length);j++)combine([...S,eligible6[j]],j+1);
}
combine([],0);
function gap(p,J){
 const m=4*p;
 const q=a=>{const k=Math.floor(m/a);return a*k*k+(m%a)*(2*k+1);};
 const b=Math.max(0,7-J);
 return 56*(J*q(p)+b*q(14)+(14-J-b)*q(28)+39*q(56))-53*m*m-44*m*(56-p);
}
const cap=[12,14].map(p=>Array.from({length:15},(_,j)=>gap(p,j)));
assert.deepEqual([cap[0][8],cap[1][10]],[5376,3136]);
assert(cap.every(a=>a.every((x,i)=>i===0||x>=a[i-1])));
assert.deepEqual([traversed,best,bestAt,pass,packed],[10626,18,1,115,0]);
console.log(JSON.stringify({trades:all.total,quartets:Q.length,eligible:eligible6.length,sumPat,traversed,best,bestAt,pass,packed,cap},null,2));
