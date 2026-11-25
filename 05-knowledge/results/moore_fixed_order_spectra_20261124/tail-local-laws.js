const assert=require('node:assert/strict');
const charts=[[29,640,191,158,111],[651,17616,6768,2691,2313]];
function mod(x,p){return (x%p+p)%p}
for(let [h,a,b,c,k] of charts){
 for(let p of [3,7,29,31])if(h%p===0){
  const w=Array(p).fill(0),W=Array(p*p).fill(0);
  for(let v=0;v<p;v++)for(let u=0;u<p;u++){
   const K=2*(a*v*u+b*v+c*u)+k;
   w[mod(K,p)]++;
   for(let q=0;q<p;q++)W[mod(K-2*h*q,p*p)]++;
  }
  for(let j=0;j<p*p;j++)assert.equal(W[j],w[j%p]);
  const biased=w.map((n,i)=>[i,n]).filter(([i,n])=>n!==p-1);
  assert.deepEqual(biased,(p===3?[[0,9],[1,0],[2,0]]:[[({7:1,29:6,31:3})[p],2*p-1]]));
  console.log(h,p,biased);
 }
}
assert.equal(29**2-2*28,785);
assert.equal((29-2)*29,783);
for(let p=3n;p<=301n;p+=2n){
 const E=(1n<<(p-1n))*((1n<<p)-1n);
 assert.equal(E%3n,1n);
}
console.log('all base supports and p² fibers exact');
