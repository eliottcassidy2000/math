const assert=require('node:assert/strict');
const p=17,m=8,U=(1<<p)-1;
const rep=x=>Math.min(x,p-x);
function mask(S,k){let t=0;for(let x of S){x=x*k%p;if(x>m)t|=1<<(p-x-1);}return t;}
const ways=[];
function length(S){
 let a=0,b=1,at=[0,1];
 for(let q=1;q<=2*Math.max(...S);q++)for(let r=0;r<q;r++){
  let x=Math.min(...S.map(v=>{let z=v*r%q;return Math.min(z,q-z);}));
  if(x*b>a*q){a=x;b=q;at=[r,q];}
 }
 const gcd=(a,b)=>b?gcd(b,a%b):a;const d=gcd(a,b);return [a/d,b/d,at];
}
function paths(S){
 let P=Array.from({length:p},(_,i)=>S.reduce((s,t)=>s|(1<<((i+p-t)%p)),0));
 let dp=new Float64Array((U+1)*p);
 let H=0;for(let s=1;s<=U;s++)for(let b=s;b;b&=b-1){
  let i=31-Math.clz32(b&-b),r=s^(b&-b),x=r?0:1;
  for(let q=r&P[i];q;q&=q-1)x+=dp[r*p+31-Math.clz32(q&-q)];
  dp[s*p+i]=x;if(s===U)H+=x;
 }
 const h=new Float64Array((1<<m)*p);h[0]=1;
 for(let a=0;a<(1<<m);a++)for(let v=0;v<p;v++)if(h[a*p+v]){
  for(const s of S){const w=(v+s)%p;if(!w)continue;
    const b=1<<(rep(w)-1);if(!(a&b))h[(a|b)*p+w]+=h[a*p+v];}
 }
 let F=0;for(let v=1;v<p;v++)F+=h[((1<<m)-1)*p+v];
 assert.equal(H%p,0);assert.equal((H/p-F)%2,0);
 return [H,F,(H/p-F)/2];
}
for(let s=0;s<1<<m;s++){
 const S=Array.from({length:m},(_,i)=>s>>i&1?p-i-1:i+1);
 const orbit=new Set(Array.from({length:p-1},(_,i)=>mask(S,i+1)));assert.equal(orbit.size,16);
 if(s===Math.min(...orbit))ways.push({s,S,paths:paths(S),lonely:length(S)});
}
assert.equal(ways.length,16);
for(const w of ways)console.log(JSON.stringify(w));
const N=p-1;
for(const v of [[1,2,3,4,5,6,7,8],[1,1,1,1,1,1,1,1],[1,2,3,4,5,11,7,9]])for(let b=1;b<=8;b++){
 const B=Array.from({length:p-1},(_,i)=>i+1).filter(x=>rep(x)<=b);
 const sq=new Set(Array.from({length:m},(_,i)=>(i+1)**2%p)),chi=x=>sq.has(x%p)?1:-1;
 for(let e of [1,-1]){
  const A=Array.from({length:p-1},(_,i)=>i+1).filter(x=>chi(x)===e);
  const L=A.reduce((u,t)=>u+v.filter(w=>B.includes(w*t%p)).length,0);
  assert.equal(2*L,v.length*B.length+e*v.reduce((a,w)=>a+chi(w),0)*B.reduce((a,w)=>a+chi(w),0));
 }
}
console.log('square-class cover identities checked');
console.log('same tournament via reflection',JSON.stringify(length([9,10,11,12,13,14,15,16])));
