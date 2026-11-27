const assert=require('node:assert/strict');
const N=200000;
let primes=[],div=new Int32Array(N+1),rad=new Int32Array(N+1);rad[1]=1;
for(let i=2;i<=N;i++){
 if(!div[i]){div[i]=i;primes.push(i);for(let j=i;j<=N;j+=i)if(!div[j])div[j]=i;}
 const p=div[i],q=i/p;rad[i]=rad[q]*(q%p? p:1);
}
let lp=0,c=1;
for(const p of primes){lp+=3*Math.log(p)/(p*p-1);c*=1-3/(p*(p+1));}
console.log('prime cutoff',N,'mean-log-D',lp.toFixed(8),'P(D=1)',c.toFixed(8));
const density=d=>{let p=1;for(let a=d;a>1;){const l=div[a];p*=3*(l-1)/(l*l+l-3);while(a%l===0)a/=l;}return c*p/d;};
const g=(a,b)=>b?g(b,a%b):a;
for(const X of [100,400,1500]){
 let n=0,logD=0,logR=0,count=new Map;
 for(let a=1;a<=X;a++)for(let b=1;b<=X;b++)if(g(a,b)===1){
  const D=(a/rad[a])*(b/rad[b])*((a+b)/rad[a+b]);
  const R=rad[a]*rad[b]*rad[a+b];
  assert.equal(R*D,a*b*(a+b));++n;logD+=Math.log(D);logR+=Math.log(R);
  if(D<=10)count.set(D,1+(count.get(D)||0));
 }
 console.log('X',X,'primitive',n,'avg log D',(logD/n).toFixed(8),'avg log R -3logX',(logR/n-3*Math.log(X)).toFixed(8));
 for(const d of [1,2,3,4,5,6,8,9,10]){
  const emp=(count.get(d)||0)/n;
  if(X===1500)assert.ok(Math.abs(emp-density(d))<0.002);
  console.log(d,emp.toFixed(7),density(d).toFixed(7));
 }
}
for(const p of [2,3,5,7,11])for(const k of [1,2,3]){
 const q=p**k;
 let n=0,e=0,f=0;
 for(let a=0;a<q;a++)for(let b=0;b<q;b++)if(a%p||b%p){
  n++; if(a===0||b===0||(a+b)%q===0)e++;
  if(a===0)f++;
 }
 assert.equal(e*p**(k-1)*(p+1),3*n);
 assert.equal(f*3,e);
}
console.log('15 local congruence counts checked');
