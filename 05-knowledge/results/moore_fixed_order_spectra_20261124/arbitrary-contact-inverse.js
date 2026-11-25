'use strict';
const A=require('node:assert/strict');
const C=16029731,Q=45312000,h=14355,J=19413465,R=37347060,s=13050,e=5104;
function q(x){
 let z=0;
 while(x){
  const t=x%2;x=(x-t)/2;
  if(t)z+=Math.floor((x+1)/2);
 }
 return z;
}
function H(b,j){
 return C+Q*b*j+(J+h)*b+R*j-2*h*(q(b)+q(j))-s*(b%2)-e*(j%2);
}
A.equal(H(11,6),3444253511);
function rev(x,b){
 let r=0;
 while(b--){const t=x%2;x=(x-t)/2;r=2*r+t;}
 return r;
}
function inverse(target){
 if(target<C)return [];
 const out=new Set();
 for(let a=0;Q*a*a<=target-C;a++)for(let chart=0;chart<2;chart++){
  const D=chart?H(a,0):H(0,a);
  const L=Q*a+(chart?R:J+h),P=chart?e:s;
  for(let r=0;r<2;r++){
   const T=target-D-L*r+P*r;
   if(T<0)continue;
   const slope=2*L-r*h;
   const lo=Math.max(0,Math.ceil((a-r)/2),Math.ceil(T/slope));
   const hi=Math.floor((T+r*h)/(slope-h));
   if(lo>hi)continue;
   let B=1,bits=0;
   while(lo>=2*B&&h*(Math.floor(lo/(2*B))-1)>(slope+h)*(2*B-1)+r*h){B*=2;bits++;}
   const f=x=>slope*x-2*h*q(x)-r*h*(x%2);
   for(let prefix=Math.floor(lo/B);prefix<=Math.floor(hi/B);prefix++){
    let l=0,u=B;
    while(l<u){
     const m=Math.floor((l+u)/2);
     if(f(prefix*B+rev(m,bits))>T)l=m+1;else u=m;
    }
    if(l===B)continue;
    const x=prefix*B+rev(l,bits);
    if(x<lo||x>hi||f(x)!==T)continue;
    const y=2*x+r,pt=chart?[a,y]:[y,a];
    A.equal(H(...pt),target);
    out.add(pt.join(','));
   }
  }
 }
 return [...out].sort();
}
for(let b=0;b<30;b++)for(let j=0;j<30;j++)A(inverse(H(b,j)).includes([b,j].join(',')));
for(let p of [13,17,19,23]){
 const E=2**(p-1)*(2**p-1);
 console.log(p,[-3,-1,1].map(d=>inverse(E+d)));
}