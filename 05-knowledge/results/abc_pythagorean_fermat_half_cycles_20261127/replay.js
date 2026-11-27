const assert=require('node:assert/strict');
const mod=(x,p)=>((x%p)+p)%p;
function pow(a,k,p){let x=1;a=mod(a,p);while(k){if(k&1)x=x*a%p;a=a*a%p;k>>=1;}return x;}
for(const p of [17,257,65537]){
 const L=(p-1)/2,g=3*pow(2,p-2,p)%p,l=g*g%p,i=Array.from({length:p},(_,j)=>j).find(j=>j*j%p===p-1);
 let all=new Set;
 for(const a of [1,g]){
  let x=a,V=new Set,s=0,t=0,hs=0,ht=0;
  for(let j=0;j<L;j++){
   assert.ok(!V.has(x)&&!all.has(x));V.add(x);all.add(x);
   let X=mod((x+pow(x,p-2,p))*pow(2,p-2,p),p), Y=mod((x-pow(x,p-2,p))*pow(2*i,p-2,p),p);
   assert.equal(mod(X*X+Y*Y,p),1); assert.equal(mod(X+i*Y,p),x);
   s+=Math.cos(2*Math.PI*x/p);t+=Math.sin(2*Math.PI*x/p);
   hs+=Math.cos(2*Math.PI*j/L);ht+=Math.sin(2*Math.PI*j/L);
   x=x*l%p;
  }
  assert.equal(x,a);assert.ok(Math.hypot(hs,ht)<1e-6&&Math.abs(t)<1e-6);
  assert.ok(Math.abs(2*s+1-(a===1?1:-1)*Math.sqrt(p))<1e-6);
  console.log('phase',p,a,'length',L,'additive sum',s.toFixed(8));
 }
 assert.equal(all.size,p-1);
}
for(const p of [3,5,7,13])for(const j of [1,2]){
 const q=p**j;let n=0,s=0;
 for(let a=0;a<q;a++)for(let b=0;b<q;b++){
  if(a%p===0&&b%p===0)continue;n++;
  if([a,b,a-b,a+b,a*a+b*b].some(t=>t%q===0))s++;
 }
 assert.equal(s*(p+1)*p**(j-1),(p%4===1?6:4)*n);
}
console.log('eight odd-prime local counts checked');
const N=2000000,lp=new Int32Array(N+1),rad=new Int32Array(N+1);rad[1]=1;
let c=.5,mean=2*Math.log(2);
for(let i=2;i<=N;i++){
 if(!lp[i]){lp[i]=i;for(let j=i+i;j<=N;j+=i)if(!lp[j])lp[j]=i;
  if(i>2){let r=i%4===1?6:4;c*=1-r/(i*(i+1));mean+=r*Math.log(i)/(i*i-1);}
 }
 const q=i/lp[i];rad[i]=rad[q]*(q%lp[i]?lp[i]:1);
}
const gcd=(a,b)=>b?gcd(b,a%b):a;
console.log('predicted P(D=2)',c.toFixed(8),'mean log D',mean.toFixed(8));
for(const M of [100,400,1000]){
 let count=0,both=0,t=0;
 for(let m=2;m<=M;m++)for(let n=1;n<m;n++){
  if((m+n)%2===0||gcd(m,n)!==1)continue;
  count++;const a=m*m-n*n,b=2*m*n,z=m*m+n*n;
  t+=Math.log(a/rad[a])+Math.log(b/rad[b])+Math.log(z/rad[z]);
  both+=+(a===rad[a]&&b===2*rad[b]&&z===rad[z]);
 }
 console.log('primitive triangles',M,count,'P(D=2)',(both/count).toFixed(8),'avg log D',(t/count).toFixed(8));
}
