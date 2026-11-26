"use strict";
const fs=require("node:fs"),assert=require("node:assert/strict");
const {boards,catalogue}=require("./minimum_trades");

function catalogueCheck(cat,rows){
 const qs=new Map();
 const seen=new Set();
 const r=cat(rows,(q,old,newer)=>{
  const k=q.join(), v=qs.get(k)||{m:0,mods:q.map(()=>Array.from({length:3},()=>new Set()))};
  for(const arr of [old,newer]){
   const sorted=[...arr].sort((a,b)=>a[0]-b[0]);
   const key=sorted.map(a=>a.join(",")).join(";");
   if(seen.has(key))throw Error("duplicate full profile");
   seen.add(key);
   sorted.forEach((row,i)=>{for(let j=0;j<3;j++)v.mods[i][j].add(row[j+1]%4);});
  }
  v.m++;qs.set(k,v);
 });
 const hist={};
 for(const x of qs.values()){
  const a=x.mods.flat().map(s=>s.size);
  if(a.some(z=>z!==(x.m===14?2:1)))throw Error("residue inventory");
  hist[x.m]=(hist[x.m]||0)+1;
 }
 if(r.total!==1057||seen.size!==2114||JSON.stringify(hist)!=='{"1":4,"6":4,"7":141,"14":3}')throw Error("bad catalogue");
 return {profiles:seen.size,hist};
}
const Q=(m,a)=>{const k=Math.floor(m/a);return a*k*k+(m%a)*(2*k+1);};
function maxCaps(f,n){return [...f,...Array(n).fill(14)].sort((a,b)=>b-a).slice(0,n);}
function pooled(h,s,l,p,t,improved=true){
 const m=s*p,n=57-s,z=4*(h+s-54)+2*t,L=54-h-l;
 const d=(p>=6&&improved?[0,4,6,6,6,6,6,6]:[0,4,6,6,4,4,2,2])[l];
 const J=Math.max(0,Math.ceil((5*(s-l)-z+Math.max(0,d-l*(s-l)))/4));
 if(d>z||J>L)return null;
 const m2=J*Q(m,p)+(n-L)*Q(m,56);
 const R=n*m*m+(8-s)*(s+7)*m*(56-p);
 const gap=56*(m2+maxCaps([28,28,28],L-J).reduce((x,a)=>x+Q(m,a),0))-R;
 return {gap,m,n,J,L,z,m2,R};
}
const haveEq=(v,s,l,p)=>v.gap===0&&v.J===0&&p===8&&s*p<=56&&l>=5&&v.L>=3;
function screen(t,mode){
 const out=[];
 for(let h=0;h<=2*t;h++)for(let s=Math.max(2,54-h-Math.floor(t/2));s<=7;s++)
 for(let l=0;l<=s;l++)for(let p=2;p<=14;p+=2) {
  if(p<s-1)continue;
  const v=pooled(h,s,l,p,t);if(!v)continue;
  if(mode===0&&v.gap<=0)out.push([h,s,l,p,v.gap,haveEq(v,s,l,p)]);
  if(mode===1&&v.gap<=0&&!haveEq(v,s,l,p))out.push([h,s,l,p,v.gap]);
 }
 return out;
}
function rowCount(xs){return new Set(xs.map(x=>x.slice(0,2).join())).size;}

const state=catalogueCheck(catalogue,boards(fs.readFileSync(process.argv[2],"utf8")).full);
const rows128=[[39,7,6,8],[39,7,7,8],[40,6,6,6],[40,6,6,8],[40,6,6,10],[40,6,6,12],[40,7,6,8],[40,7,7,8]];
const gaps128=rows128.map(row=>[row,pooled(...row,20).gap]);
assert.deepEqual(gaps128.map(r=>r[1]),[28224,18816,9408,0,3360,40096,18816,9408]);
const A=screen(20,0),B=screen(21,0),C=screen(21,1);
assert.deepEqual(A,[[40,6,6,8,0,true]]);
assert.equal(screen(20,1).length,0);
assert.equal(B.length,114);
const eq129=B.filter(x=>x[5]).map(x=>x.slice(0,4));
assert.deepEqual(eq129,[[40,6,6,8],[41,6,5,8],[41,7,7,8],[42,7,6,8]]);
assert.equal(C.length,110);assert.equal(rowCount(C),11);
const summary={};for(const row of C){const k=row.slice(0,2).join();summary[k]=(summary[k]||0)+1;}
for(const m of [48,56]){
 assert(Q(m,27)>Q(m,28));assert(Q(m,14)>Q(m,28));
 assert.equal(Q(m,56),m);
}
console.log(JSON.stringify({state,gaps128,eq129,summary,tuples129:C},null,2));
