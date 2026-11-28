"use strict";
const assert=require("node:assert/strict");
function F(x){let h=1,b=0,z=0;for(const d of x.toString(2)){let v=+d;if(!v)h+=z+b;z=2*z+v;b=v;}return h;}
function H(w){let h=1;for(let i=0;i<w.length;i++)for(let j=i+1;j<w.length;j++)if(w[i]==="1"&&w[j]==="0")h+=2**(j-i-1)+(j===i+1?1:0);return h;}
const rc=s=>[...s].reverse().map(x=>x==="1"?"0":"1").join("");
function L(p,a,b){let x=2*p-a,y=2*p+2-b;return 2*x*y+2*p+F(x)+F(y)-1;}
function bounds(p){return [L(p,1,1),L(p,0,0)];}
assert.deepEqual([L(177,0,0),L(178,1,1)],[252965,254033]);
assert(252965<503**2&&504**2<254033);
for(let p=1;p<=5000;p++){
  const q=[L(p,1,1),L(p,0,1),L(p,1,0),L(p,0,0)];
  assert.equal(bounds(p+1)[0]-q[3],6*p+6);
  assert(q[0]<Math.min(q[1],q[2])&&q[3]>Math.max(q[1],q[2]));
  if(p<=250)for(let a=0;a<2;a++)for(let b=0;b<2;b++){
    const u=(2*p-a).toString(2),v=rc((2*p+2-b).toString(2));
    const x=L(p,a,b);
    assert.equal(H(u+"0"+v),x);assert.equal(H(u+"1"+v),x+2);
  }
}
function bruteMiss(X){
 let ls=[], p=1;
 while(L(p,1,1)<(X+1)**2){
   ls.push(...[L(p,1,1),L(p,1,0),L(p,0,1),L(p,0,0)].sort((a,b)=>a-b));p++;
 }
 let i=0,miss=0;
 for(let h=5;h<=X;h+=2){
   while(ls[i]<=h*h)i++;
   if(!(ls[i]<(h+1)**2-3))miss++;
 }
 return miss;
}
function interMiss(X){
 let hits=0;
 for(let p=1;;p++){
   const z=L(p,0,0);
   if(z>X*X)break;
   let h=Math.ceil(Math.sqrt(z));if(h%2===0)h++;
   if(h<=X&&(h+1)**2<=L(p+1,1,1)+3)hits++;
 } return hits;
}
assert.equal(bruteMiss(100000),1077);
assert.equal(interMiss(1000000),10752);
function small(){
 let common=0;
 for(let n=3;n<=6;n++){
   const edges=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)edges.push([i,j]);
   const ngraphs=1<<edges.length, ap=new Uint8Array(ngraphs);
   for(let m=0;m<ngraphs;m++){
      const a=Array(n).fill(0);
      for(let k=0;k<edges.length;k++){const [i,j]=edges[k];if(m>>k&1)a[i]|=1<<j;else a[j]|=1<<i;}
      let apices=(1<<n)-1;
      for(let i=0;i<n;i++)for(let j=0;j<n;j++)if(a[i]>>j&1)
        for(let k=0;k<n;k++)if((a[j]>>k&1)&&(a[k]>>i&1))apices&=(1<<i)|(1<<j)|(1<<k);
      ap[m]=apices;
   }
   for(let m=0;m<ngraphs;m++)if(ap[m])for(let k=0;k<edges.length;k++)if(!(m>>k&1)){
      const z=m|1<<k;
      if(ap[z]){assert(ap[m]&ap[z]);common++;}
   }
 }
 return common;
}
const h=503,x=h-1,y=h+1;
assert(H(x.toString(2)+rc(y.toString(2)))>h*h);
console.log({fullMiss100k:bruteMiss(100000),interGapMiss1m:interMiss(1000000),commonApexEdgesThrough6:small()});
