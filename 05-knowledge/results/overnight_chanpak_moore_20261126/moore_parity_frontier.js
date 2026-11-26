"use strict";
const fs=require("node:fs"),assert=require("node:assert/strict");
const {boards,catalogue}=require("./minimum_trades");

function replay(cat,rows){
 const m=new Map(), allE=new Set();
 const total=cat(rows,(q,old,newer)=>{
  const id=q.join(), v=m.get(id)||{q, rows:[], M:0};
  v.M++;v.rows.push(old,newer);m.set(id,v);
  for(const z of [old,newer]){
   const u=[...z].sort((a,b)=>a[0]-b[0]);
   for(let i=0;i<4;i++)for(let j=i+1;j<4;j++){
    const e=[u[i][0],u[i][1],u[j][0],u[j][1]].join();
    if(allE.has(e))throw Error("pair repeated");allE.add(e);
   }
  }
 });
 if(total.total!==1057||allE.size!==12684)throw Error("source inventory");
 const Q=[...m.values()].filter(x=>x.M>=2);
 const labels=Array.from({length:56},(_,i)=>i).filter(i=>i!==0&&i!==28);
 const pairs=Array.from({length:56},()=>Array.from({length:56},()=>[]));
 for(let i=0;i<Q.length;i++){
  const {q,rows}=Q[i];
  for(const a of q)for(const b of q)if(a!==b)pairs[a][b].push(i);
  Q[i].maps=new Map();Q[i].mods=new Map();
  for(const a of q)for(const b of q)if(a!==b){
   let z=new Int8Array(56).fill(-1);
   for(const r of rows){
    const x=r.find(y=>y[0]===a)[1], y=r.find(y=>y[0]===b)[1];
    if(z[x]>=0)throw Error("not a partial bijection");z[x]=y;
    Q[i].mods.set(a,(Q[i].mods.get(a)||0)|(1<<(x%4)));
   }
   Q[i].maps.set(56*a+b,z);
  }
 }
 const keys=new Set(), labelSets=new Set(), emptyHist=[0,0,0,0,0];
 for(let i=0;i<Q.length;i++)for(let j=i+1;j<Q.length;j++){
  const x=Q[i].q,y=Q[j].q, shared=x.filter(t=>y.includes(t));
  for(let ai=0;ai<shared.length;ai++)for(let bi=ai+1;bi<shared.length;bi++){
   const a=shared[ai],b=shared[bi];
   for(const d of x.filter(z=>!y.includes(z)))for(const c of y.filter(z=>!x.includes(z))){
    const opts=pairs[c][d].filter(k=>k!==i&&k!==j&&!Q[k].q.includes(a)&&!Q[k].q.includes(b));
    for(const oa of opts)for(const ob of opts)if(oa!==ob){
     const mp=[[a,oa],[b,ob],[c,i],[d,j]].sort((a,b)=>a[0]-b[0]);
     const key=mp.map(t=>t.join(",")).join(";");
     if(keys.has(key))continue;keys.add(key);labelSets.add(mp.map(t=>t[0]).join());
     let e=0;
     for(const [target] of mp){
      let bits=15;for(const [,src] of mp)if(Q[src].q.includes(target))bits&=Q[src].mods.get(target);
      if(!bits)e++;
     }
     emptyHist[e]++;
    }
   }
  }
 }
 const allTriangles={};let tested=0;
 for(let ai=0;ai<labels.length;ai++)for(let bi=ai+1;bi<labels.length;bi++)for(let ci=bi+1;ci<labels.length;ci++){
  const a=labels[ai],b=labels[bi],c=labels[ci];
  for(const qa of pairs[b][c])if(!Q[qa].q.includes(a))
  for(const qb of pairs[c][a])if(!Q[qb].q.includes(b)&&qa!==qb)
  for(const qc of pairs[a][b])if(!Q[qc].q.includes(c)&&qa!==qc&&qb!==qc){
   tested++;
   const ab=Q[qc].maps.get(56*a+b),bc=Q[qa].maps.get(56*b+c),ca=Q[qb].maps.get(56*c+a);
   const D=Array.from({length:56},(_,x)=>ab[x]<0||bc[ab[x]]<0?-1:ca[bc[ab[x]]]);
   function cycles(k,step){
    const all=[];
    for(let x=0;x<k;x++){
     let u=step(x),len=1;
     while(u>x){u=step(u);len++;}
     if(u===x)all.push(len);
    }
    return all.sort((a,b)=>a-b).join(",");
   }
   const p=cycles(28,x=>D[x]<0?-1:D[x]%28), q=cycles(56,x=>D[x]);
   const key=p+";"+q;allTriangles[key]=(allTriangles[key]||0)+1;
  }
 }
 return {quartets:Q.length,profiles:total.total*2,pairs:allE.size,maps:keys.size,labelSets:labelSets.size,emptyHist,tested,allTriangles};
}

const r=replay(catalogue,boards(fs.readFileSync(process.argv[2],"utf8")).full);
assert.deepEqual([r.maps,r.labelSets,r.tested],[918,219,8914]);
assert.deepEqual(r.emptyHist,[0,0,96,262,560]);
assert.deepEqual(r.allTriangles,{";":8617,"1,1,1,1,1,1,1;2,2,2,2,2,2,2":122,"7;7,7":104,"7;14":68,"1,1,1,1,1,1;2,2,2,2,2,2":3});
const s=fs.readFileSync(require("node:path").join(__dirname,"pooled_capacity.out"),"utf8");
const xs=s.split("# (h,s,l,p,gap)\n")[1].trim().split("\n").map(x=>x.split(" ").map(Number));
const next=xs.filter(([h,s,l,p])=>!((h===40&&s===5)||(h===41&&s===4&&l===4)
 ||(l===3&&(p===6||p===10)&&((h===41&&s===4)||(h===42&&s===3)))));
assert.equal(next.length,82);
assert.equal(new Set(next.map(x=>x.slice(0,2).join())).size,10);
console.log(JSON.stringify(r,null,2)+"\n\n# (h,s,l,p,gap)\n"+next.map(x=>x.join(" ")).join("\n"));
