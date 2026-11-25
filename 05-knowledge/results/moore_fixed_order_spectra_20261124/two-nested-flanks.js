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
const Ts=[pair[0]-2,...pair];
for(const [h,a,b,c,k] of [[29,640,191,158,111],[651,17616,6768,2691,2313]])
 for(const T of Ts)for(let v=0;v<=Math.floor((T-k)/(2*b-h));v++)
  assert(!search(k+2*b*v-2*h*q(v),2*(a*v+c),h,T));
function push(old,w){
 const out=Array(18).fill(0),ix=(a,b,c)=>(3*a+b)*2+c;
 for(let a=0;a<3;a++)for(let b=0;b<3;b++)for(let c=0;c<2;c++){
  const x=old[ix(a,b,c)];
  out[ix(w?2:1,b,c)]+=x;
  out[ix(a,b||(w?1:2),c)]+=x;
  out[ix(a,b,1)]+=x;
 }
 return out;
}
let num=0;
function loop(s,n){
 const r=push(s,0);let C=0,R=0,h=0;num++;
 for(let a=0;a<=1;a++)for(let b=0;b<=1;b++)for(let c=0;c<=1;c++){
  const x=r[(3*a+b)*2+c];
  C+=((1-b)+(1-c)*a)*x;
  R+=(1+a)*x;
  h+=(1-c)*x;
 }
 for(const [J,S,F] of [[162,11,158],[129,16,191]]){
  const D=h*(111-S)+29*(C-1)+(R-h)*J+S,L=2*(320*(R-h)+h*F);
  for(const T of Ts)assert(!search(D,L,29*h,T));
 }
 if(n<16){loop(r,n+1);loop(push(s,1),n+1);}
}
assert(129*(2**17-1)>Math.max(...Ts));
loop(push([1,...Array(17).fill(0)],1),1);
assert.equal(num,65535);
console.log(num,'strong words');
