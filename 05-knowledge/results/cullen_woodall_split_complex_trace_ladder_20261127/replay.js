const mod=(a,m)=>(a%m+m)%m;
function traceSeed(t,k,m){
  // alpha^2=t*alpha+1, store a+b*alpha.
  const mul=([a,b],[c,d])=>[mod(a*c+b*d,m),mod(a*d+b*c+t*b*d,m)];
  let r=[1n,0n], z=[0n,1n];
  for(;k;k>>=1n,z=mul(z,z))if(k&1n)r=mul(r,z);
  return mod(2n*r[0]+t*r[1],m);
}
function ladder(n,t){
  let k=n,e=n;while(k%2n===0n){k/=2n;e++;}
  const X=n<<n,C=X+1n,W=X-1n,M=C*W;
  let r=traceSeed(t,k,M);
  r=mod(r*r+2n,M);
  for(let j=1n;j<e-2n;j++)r=mod(r*r-2n,M);
  const next=mod(r*r-2n,M),g=mod((X-1n)*r+(X+1n)*next,M);
  return [r%C===0n,next%W===0n,g===0n];
}
const assert=(b)=>{if(!b)throw Error("replay assertion failed")};
const gcd=(a,b)=>b?gcd(b,a%b):a;
function psmall(n){if(n<2n)return false;for(let i=2n;i*i<=n;i++)if(n%i===0n)return false;return true;}
for(let n=2n;n<=18n;n++){
  let X=n<<n,C=X+1n,W=X-1n;
  const c=psmall(C),w=psmall(W);
  let A=0n,B=0n;
  for(let t=0n;t<(n<=6n?C:128n);t++){
    const [a,b,g]=ladder(n,t);assert(g===(a&&b));
    if(n>=3n)assert(!a||c);
    assert(!b||w);
    if(a)A++;if(b&&t<W)B++;
  }
  if(n<=6n){
    if(c)assert(A===(C-1n)/4n);
    if(w)assert(B===(W+1n)/2n);
  }
}
assert(ladder(2n,4n)[0]);assert(!psmall(9n));
for(const n of [141n,4713n]){
  const X=n<<n,C=X+1n;
  assert(gcd(5n,C)===1n);
  function invmod(a,m){let [r,s,x,y]=[a,m,1n,0n];while(s){let q=r/s;[r,s,x,y]=[s,r%s,y,x-q*y];}assert(r===1n);return mod(x,m);}
  const t=mod(5n-invmod(5n,C),C);
  assert(ladder(n,t)[0]);console.log('Cullen certificate, n='+n);
}
for(const n of [2n,3n,6n,30n,75n,81n,115n,123n,249n,362n,384n]){
  let t=0n;while(t<100n&&!ladder(n,t)[1])t++;
  assert(t<100n);console.log('Woodall certificate, n='+n+', smallest seed t='+t+', golden seed='+ladder(n,1n)[1]);
}
console.log('All checks passed, including the Cullen false positive at 9.');
