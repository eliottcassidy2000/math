'use strict';
const assert=require('node:assert/strict');
function chainA(w){
 const A=Array(w.length+1).fill(0),t=w.length;
 for(let i=0;i<t;i++){
  for(let j=i+1;j<t;j++)A[i]|=1<<j;
  if(w[i]==='1')A[t]|=1<<i;else A[i]|=1<<t;
 }
 return A;
}
function forkA(mask){
 const A=Array(6).fill(0);
 for(let i=0,k=0;i<6;i++)for(let j=i+1;j<6;j++,k++){
  if(mask>>k&1)A[j]|=1<<i;else A[i]|=1<<j;
 }
 return A;
}
function lanes(A,out){
 const n=A.length,state=Array(9).fill(0);
 const last=i=>(out>>i&1)?2:1,first=i=>(out>>i&1)?1:2;
 function go(used,active,l,r){
  if(active>=0&&l<0)go(used,-1,last(active),0);
  if(used===(1<<n)-1){
   if(l>=0)state[l*3+r]++;
   return;
  }
  for(let i=0;i<n;i++)if(!(used>>i&1)&&(active<0||(A[active]>>i&1)))
   go(used|(1<<i),i,l,l>=0?(r||first(i)):0);
 }
 go(0,-1,-1,0);go(0,-1,0,0);
 return state;
}
function join(X,Y){
 const s=Array(9).fill(0);
 for(let l=0;l<3;l++)for(let r=0;r<3;r++)for(let m=0;m<3;m++)for(let t=0;t<3;t++)
  s[3*(m||l)+(r||t)]+=X[3*l+r]*Y[3*m+t];
 return s;
}
function count(X){return X[0]+X[1]+X[3]+X[4]}
const U=lanes(chainA('1100'),15),V=lanes(chainA('1001110'),47);
const hl=X=>X[3]+X[6],hr=X=>X[1]+X[2];
assert.deepEqual([hl(U),U[3],hr(V),V[1]],[11,1,45,41]);
const pair=[11524,11588].map(x=>count(join(join(U,lanes(forkA(x),7)),V)));
assert.deepEqual(pair,[16029731,16029733]);
for(const H of pair)for(let d=2;d*d<=H;d++)assert(H%d);
console.log(pair);
function q(n){
 let o=0;
 while(n){
  const b=n%2;n=(n-b)/2;
  if(b)o+=Math.floor((n+1)/2);
 }
 return o;
}
function search(D,L,K,T){
 const n=T-D;if(n<0)return false;
 for(let u=Math.ceil(n/L);u<=Math.floor(n/(L-K));u++)
  if(L*u-2*K*q(u)===n)return true;
 return false;
}

function chain(w){
 const a=Array(w.length).fill(0);let out=0;
 for(let i=0;i<w.length;i++){
  if(w[i]==='1')out|=1<<i;
  for(let j=i+1;j<w.length;j++)a[i]|=1<<j;
 }
 return lanes(a,out);
}
const ext=[];
for(const mask of [11524,11588]){
 const body=join(join(U,lanes(forkA(mask),7)),V);
 assert.deepEqual([hl(body),body.reduce((x,y)=>x+y),body.filter((v,i)=>i%3<2).reduce((x,y)=>x+y),
                  body[3]+body[4]+body[5],body[3],body[1]],
 [14355,45312000,37347060,19413465,5104,13050]);
 ext.push(count(join(join(chain('1011'),body),chain('100'))));
}
assert.deepEqual(ext,[3444253511,3444253513]);
for(const x of ext)for(let d=2;d*d<=x;d++)assert(x%d);
console.log(ext);
