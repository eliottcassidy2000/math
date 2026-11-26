const A=require("node:assert/strict");
function pairs(n){return Array.from({length:n},(_,i)=>Array.from({length:n-i-1},(_,j)=>[i,i+j+1])).flat()}
function graph(n,free,mask){
 const adj=Array(n).fill(0), flip=new Set(free.filter((_,i)=>(mask>>i)&1).map(([i,j])=>i*n+j));
 for(const [i,j] of pairs(n))if((i===0&&j===n-1)!==flip.has(i*n+j))adj[j]|=1<<i;else adj[i]|=1<<j;
 return adj;
}
function paths(a){
 const n=a.length,U=(1<<n)-1,h=new Int32Array(U+1),dp=new Int32Array((U+1)*n);h[0]=1;
 for(let S=1;S<=U;S++)for(let b=S;b;b&=b-1){let v=31-Math.clz32(b&-b),R=S^(1<<v),x=R?0:1;
 for(let J=R;J;J&=J-1){let u=31-Math.clz32(J&-J);if((a[u]>>v)&1)x+=dp[n*R+u];}
 h[S]+=dp[n*S+v]=x;}
 return h;
}
function canon(a){
 const n=a.length,E=pairs(n),p=[],d=a.map(x=>x.toString(2).split("1").length-1),s=[...d].sort((a,b)=>a-b);
 let best=Infinity,aut=0,orig=0;for(let k=0;k<E.length;k++){const [i,j]=E[k];orig|=((a[i]>>j)&1)<<k;}
 function f(used){
 if(p.length===n){let v=0;for(let k=0;k<E.length;k++){const[i,j]=E[k];v|=((a[p[i]]>>p[j])&1)<<k;}
 if(v<best){best=v;aut=0;}if(v===best)aut++;
 return;}
 for(let i=0;i<n;i++)if(!(used>>i&1)&&d[i]===s[p.length]){p.push(i);f(used|1<<i);p.pop();}
 }
 f(0);return [best,aut];
}
function packs(a){
 const n=a.length,U=(1<<n)-1,cyc=new Int32Array(U+1);
 function walk(root,v,S){
 if((a[v]>>root&1)&&S.toString(2).split("1").length%2===0)cyc[S]++;
 for(let R=a[v]&U&~S&~((1<<root)-1);R;R&=R-1){let j=31-Math.clz32(R&-R);walk(root,j,S|1<<j);}
 }
 for(let i=0;i<n;i++)walk(i,i,1<<i);
 const P=Array.from({length:U+1},()=>new Map);P[0].set("",1);
 for(let S=1;S<=U;S++){
 const b=S&-S;
 for(const [k,v] of P[S^b])P[S].set(k,v);
 for(let R=S;R;R=(R-1)&S)if((R&b)&&cyc[R]){
 let L=R.toString(2).split("1").length-1;
 for(const [k,v] of P[S^R]){
 const t=[...k,L].sort().join("");P[S].set(t,(P[S].get(t)||0)+v*cyc[R]);
 }}
 }
 return [...P[U]].sort((a,b)=>a[0].localeCompare(b[0]));
}
const bank={};
for(let n=3;n<=7;n++){
 const free=pairs(n).filter(([i,j])=>j!==i+1&&!(i===0&&j===n-1));
 const tab=new Map;
 for(let m=0;m<1<<free.length;m++){
  const a=graph(n,free,m),h=paths(a).at(-1);
  if(n===7&&h>25)continue;
  if(!tab.has(h))tab.set(h,new Map);
  const [k,aut]=canon(a),row=tab.get(h);
  if(!row.has(k)){const p=packs(a),u=(1<<n)-1,op=a.map((x,i)=>u^x^(1<<i)),self=canon(op)[0]===k;A.equal(1+p.filter(([k])=>k).reduce((v,[k,c])=>v+c*2**k.length,0),h);row.set(k,{m,aut,self,p});}
 }
 console.log(n,JSON.stringify([...tab].sort((a,b)=>a[0]-b[0]).map(([h,r])=>[h,r.size])));
 bank[n]=Object.fromEntries([...tab].filter(([h,r])=>r.size).map(([h,r])=>[h,[...r.values()]]));
}
for(let n=3;n<=7;n++)for(const [h,rs] of Object.entries(bank[n]))if(+h<=31)console.log('bank',n,h,JSON.stringify(rs));
for(let n=3;n<=6;n++)for(const [h,rs] of Object.entries(bank[n])){
 for(const x of rs)for(const y of rs)if(x.m<y.m&&x.aut!==y.aut&&JSON.stringify(x.p)===JSON.stringify(y.p))console.log('same profile different aut',n,h,x.m,y.m,x.aut,y.aut,JSON.stringify(x.p));
}

const v=new Map;
for(const rs of Object.values(bank[6]))for(const r of rs){
 const p=Object.fromEntries(r.p),k=[p["3"],p["33"]||0].join(",");
 if(!v.has(k))v.set(k,new Set);v.get(k).add(p["5"]||0);
}
console.log('six-vertex cycle certificate');
for(const [k,vv] of [...v].sort((a,b)=>a[0].localeCompare(b[0])))console.log(k,[...vv].sort((a,b)=>a-b));

const expected=[{3:1},{5:1},{9:2,11:1,13:1,15:2},{15:2,17:2,19:2,23:4,25:2},{25:1}];
for(let n=3;n<=7;n++)A.deepEqual(Object.fromEntries(Object.entries(bank[n]).filter(([h])=>+h<=25).map(([h,r])=>[h,r.length])),expected[n-3]);
const specialize=p=>{
 const d={};for(const [k,c] of p){const j=[...k].reduce((sum,l)=>sum+(+l-1)/2,0);d[j]=(d[j]||0)+c*2**k.length;}return d;
};
const x=bank[6][17], b=bank[5][15].find(r=>r.m===4),c=bank[6][15].find(r=>r.m===1);
A.deepEqual(specialize(x[0].p),specialize(x[1].p)); A.notDeepEqual(x[0].p,x[1].p);
A.deepEqual(b.p,c.p);A.equal(b.aut,3);A.equal(c.aut,1);
console.log("nested equivalence witnesses verified");
