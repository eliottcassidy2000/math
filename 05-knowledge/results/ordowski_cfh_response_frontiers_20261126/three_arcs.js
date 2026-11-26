
const assert=require("node:assert/strict");
const {cores,flip,dp,F,seam}=require("../overnight_chanpak_moore_20261126/commuting_unit_forks");
const [L,u,[i,j]]=cores[0], C=[14,28,24,16,33,15], r=5;
function profile(A,u){
 const [P,D]=dp(A),[Z]=dp(A,u),n=A.length,M=P.length-1;
 let R=0,V=0,Q=0;
 for(let s=0;s<=M;s++){
  let z=P[M^s];Q+=P[s]*z;
  if(s>>u&1){R+=D[s*n+u]*z;V+=Z[s]*z;}
 }
 return [P[M],D[M*n+u],Z[M],R,V,Q];
}
function join(A,u,B,v){
 const n=A.length,r=n+v,T=A.concat(B.map(x=>x<<n));
 for(let y=n;y<T.length;y++)for(let x=0;x<n;x++){
  if(x===u&&y===r)T[x]|=1<<y;else T[y]|=1<<x;
 }return T;
}
function strong(A){
 const m=(1<<A.length)-1;
 for(let r=0;r<A.length;r++){
  let s=1<<r,t=0;
  while(s!==t){t=s;for(let j=0;j<A.length;j++)if(t>>j&1)s|=A[j];}
  if(s!==m)return false;
 }return true;
}
const base=profile(C,r), edits=profile(flip(C,1,4),r);
assert.deepEqual(base,[23,1,8,24,75,296]);assert.deepEqual(edits,[23,1,10,24,87,296]);
assert.deepEqual(profile(C.concat(63),r),[23,1,0,48,75,592]);
assert.deepEqual(profile(flip(C,1,4).concat(63),r),[23,1,0,48,87,592]);
const e2=1n<<190n, q2=((529n+15312n)*e2+2n)/529n, c2=2n*e2*q2+1n;
assert.equal((15312n*e2+2n)%529n,0n);
assert.equal(15312n*e2+529n*(F(c2+e2-1n)-F(c2-1n)),-2n);
const a=529n,b=7656n,d=1992n;
const exact=c=>a*F(c-1n)+b*c-d;
const literal=[];
for(let v=0;v<8;v++){
 const X=v&1?flip(L,i,j):L,Y=v&2?flip(C,1,4):C, c=v&4?7:5;
 const B=seam(Y,r,c),A=join(X,u,B,B.length-1),[P]=dp(A);
 assert.equal(A.length,16);assert(strong(A));
 assert.equal(BigInt(P[P.length-1]), exact(BigInt(c))-BigInt((v&1?2:0)+(v&2?334:0)));
 literal.push(P[P.length-1]);
}
assert.deepEqual(literal,[38933,38931,38599,38597,54245,54243,53911,53909]);
const signed=new Set();let z=1n,period=0;
do{signed.add(String(b*z%a));signed.add(String((a-b*z%a)%a));period++;z=2n*z%a;}while(z!==1n);
assert.equal(period,253);assert.equal(signed.size,506);
const test=[];
for(let i=191;i<1000;i+=253){
 const ell=1n<<BigInt(i),beta=((a+b)*ell+2n)/a,c=2n*ell*beta+1n;
 assert.equal((b*ell+2n)%a,0n);assert.equal(beta%2n,0n);
 assert.equal(exact(c+ell),exact(c)-2n);
 test.push({i,order:13+(c-1n).toString(2).length,center:String(exact(c+ell))});
}
assert.equal(test[0].order,400);
console.log(JSON.stringify({base,edits,period,signed:signed.size,literal,test},null,2));
