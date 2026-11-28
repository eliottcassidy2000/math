"use strict";
const eq=(a,b)=>{if(a!==b)throw new Error(`${a} != ${b}`);};
const lengths=w=>Array.from(w.matchAll(/(0+|1+)/g),x=>x[0].length);
function intervals(L){
  let h=BigInt(L.length+1);
  for(let i=0;i<L.length;i++){
    let n=0;
    for(let j=i;j<L.length;j++){
      n+=L[j];
      h+=(j-i)%2===0?-(1n<<BigInt(n)):(1n<<BigInt(n));
    }
  }
  return h;
}
function prefix(w){
  let x=1n,b=1n,h=1n;
  for(const c of w.slice(1)){
    const e=BigInt(c);
    h+=(1n-e)*(x+b);
    x=2n*x+e;b=e;
  }
  return h;
}
function decode(h,r){
  let z=h-BigInt(r+1), l=[];
  for(let i=0;i<r;i++){
    let e=0, v=z;
    if(v===0n) return null;
    while(v%2n===0n){v/=2n;e++;}
    l.push(e); z+=1n<<BigInt(e);
  }
  if(l[0]<1||l.some((v,i)=>i>0&&v<=l[i-1])||l[r-1]>=l[0]+l[1])return null;
  return intervals(l)===h?l:null;
}
function subsetDP(w) {
  const n=w.length+1, p=n-1, full=(1<<n)-1;
  const edge=(i,j)=>i===p?w[j]==="1":j===p?w[i]==="0":i<j;
  const table=new Float64Array((full+1)*n);
  for(let i=0;i<n;i++) table[(1<<i)*n+i]=1;
  for(let s=1;s<full;s++)for(let i=0;i<n;i++){
    const v=table[s*n+i];
    if(v) for(let j=0;j<n;j++)if(!(s&(1<<j))&&edge(i,j))table[(s|1<<j)*n+j]+=v;
  }
  return BigInt(table.slice(full*n).reduce((a,b)=>a+b,0));
}
const target=new Map(Array.from({length:14},(_,t)=>[7*3**t,t]));
const hits=Array.from({length:14},()=>[]);
const upper=7*3**13, cap=Math.floor(Math.log2(upper-1))+2;
let nodes=0;
function walk(x,b,h,n){
  if(h>upper)return;
  if(!b){
    nodes++;
    if(target.has(h))hits[target.get(h)].push(x.toString(2));
  }
  if(n===cap)return;
  walk(2*x,0,h+x+b,n+1);
  walk(2*x+1,1,h,n+1);
}
walk(1,1,1,1);
eq(nodes,7735713);
eq(hits.map((w,t)=>w.length?t:null).filter(t=>t!==null).join(","),"2,8,9,12,13");
for(let t=0;t<hits.length;t++)if(hits[t].length){
  eq(hits[t].length,2);
  eq(hits[t][0],Array.from(hits[t][1]).reverse().map(e=>e==="0"?"1":"0").join(""));
  for(const w of hits[t]){
    eq(intervals(lengths(w)),BigInt(7*3**t)); eq(prefix(w),BigInt(7*3**t));
    if(w.length<=16)eq(subsetDP(w),BigInt(7*3**t));
  }
  console.log(t,7*3**t,...hits[t]);
}
let controls=0;
for(let n=2;n<=16;n++)for(let x=2**(n-1);x<2**n;x+=2){
  const w=x.toString(2), h=prefix(w);
  eq(intervals(lengths(w)),h);
  if(n<=8)eq(subsetDP(w),h);
  controls++;
}
let trials=0;
for(let r=2;r<=20;r+=2)for(let a=0;a<150;a++){
  const m=8*r+a, d=Math.floor(m/(2*r));
  const l=Array.from({length:r},(_,i)=>m+i*d+((i*13+a*7)%d));
  eq(decode(intervals(l),r).join(","),l.join(","));
  trials++;
}
console.log({nodes,controls,trials});
