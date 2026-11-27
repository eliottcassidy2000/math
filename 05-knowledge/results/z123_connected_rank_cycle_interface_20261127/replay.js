// Run: node replay.js
const assert = require("node:assert/strict");
function C(n,k) {
  if(k<0 || k>n) return 0n;
  let v=1n;
  for(let i=1;i<=k;i++) v=v*BigInt(n-i+1)/BigInt(i);
  return v;
}
const S1=[[1,[-1]]], S2=[[1,[-2]],[1,[-3,-1]],[-1,[-3,-2]]];
const S3=[[1,[-3]],[1,[-4,-1]],[-1,[-4,-3]],[1,[-5,-2]],[-1,[-5,-3]],
          [1,[-5,-3,-1]],[-1,[-5,-3,-2]],[-1,[-5,-4,-1]],[1,[-5,-4,-3]]];
function splits(d,es) {
  const arr=[];
  function go(a=[],h=[]) {
    if(a.length===es.length) {
      const key=[...new Set(a.filter(i=>i>0 && i<4))].join("");
      const w=h.slice(1).reduce((p,k,i)=>p*C(h[i],k),1n);
      arr.push([key,h[0],w]); return;
    }
    for(let i=a.at(-1)??0;i<=4;i++) {
      const k=d+i-4-es[a.length];
      if(k>=0&&(h.length===0||k<=h.at(-1))) go([...a,i],[...h,k]);
    }
  }
  go(); return arr;
}
function polynomial(d,terms) {
  const z=new Map();
  for(const [sign,es] of terms)
    for(const [key,k,w] of splits(d,es)) {
      const id=key+":"+k;
      z.set(id,(z.get(id)??0n)+BigInt(sign)*w);
    }
  return z;
}
function packing(q) {
  const c2=m=>m*(m-1n)/2n,c3=m=>m*(m-1n)*(m-2n)/6n;
  const m=q*q+q+1n,P=q**4n,L=q**3n*(q**3n+q*q+q+1n),B=q*q*(q*q+1n)*m;
  const D=P*c2(m)+c3(m);
  return {"":1n,1:P+m,2:L+P*m+c2(m),3:B+D,
          12:q*L+2n*(P*m+c2(m)),13:q*q*B+m*q*q*(q*q+1n)+3n*D,
          23:(q*q+q)*B+m*q**4n*(q*q+1n)+3n*D,
          123:q*(q*q+q)*B+2n*m*q**4n*(q*q+1n)+6n*D};
}
function z(q,s,d) {
  const f=packing(BigInt(q));
  return [S1,S2,S3].map(T=>[...polynomial(d,T)].reduce((sum,[id,w])=>{
    const [k,m]=id.split(":"); return sum+w*f[k]*C(s,Number(m));
  },0n));
}
function gap(q,s,d=2) {
  const [a,b,c]=z(q,s,d);return b*b-a*c;
}
for(let d=1;d<=20;d++){
  const [p1,p2,p3]=[S1,S2,S3].map(T=>polynomial(d,T));
  assert.equal(p1.get("3:"+d),1n);
  assert.equal(p2.get("3:"+(d+2)),C(d+2,2)-1n);
  assert.equal(p3.get("3:"+(d+4))??0n,BigInt((d+4)*(d+5)*d*(d+1)/12));
  assert.equal(p3.get("123:"+(d+2)),C(d+1,2));
  const maxDegree=(p,key)=>Math.max(-1,...[...p].filter(([id,w])=>w!==0n&&id.startsWith(key+":")).map(([id])=>+id.split(":")[1]));
  for (const [p,j] of [[p1,1],[p2,2],[p3,3]]) {
    assert(maxDegree(p,"")<=d+2*j-1);
    for(const k of ["1","2","12","13","23"]) assert(maxDegree(p,k)<=d+2*j-3);
  }
}
const root=q=>{
  let [a,b,c]=z(q,4,2);
  b-=3n;c-=2n*a;
  return a*a*b*b-4n*b**3n-4n*a**3n*c-27n*c*c+18n*a*b*c;
};
const primepower=q=>{
  for(let p=2;p<=q;p++){
    if(q%p!==0)continue;
    while(q%p===0)q/=p;
    return q===1;
  }
  return false;
};
for(const f of [root,q=>gap(q,4)]){
  const q=Array.from({length:500},(_,i)=>i+2).find(q=>q%2&&primepower(q)&&f(q)<0n);
  console.log("first odd prime-power with negative test at rank 6:",q, "n =",q**4+q*q+q+5);
}
for(const s of [500,550,570,580,632,664])console.log("gap sign, q=100003, s="+s+":", gap(100003,s)>0n?"+":"-");
