
const assert=require("node:assert/strict");
const {cores,F,reverse}=require("../overnight_chanpak_moore_20261126/commuting_unit_forks");
const {roots,inside}=require("../overnight_chanpak_moore_20261126/nested_fork_inverse");
const bank=[...cores,
 [[2,72,3,21,39,15,61],5,[4,6],59,87,11],
 [[2,68,17,7,43,15,61],6,[3,4],71,72,14],
 [[34,20,9,67,13,30,55],3,[4,5],103,129,11],
 [[36,65,42,19,7,90,29],6,[0,4],143,173,23]
];
const mod=(a,b)=>(a%b+b)%b;
function eg(a,b){if(!b)return [a,1n,0n];const [g,x,y]=eg(b,a%b);return[g,y,x-a/b*y];}
function sqrt(x){
 if(x<2n)return x;
 let t=1n<<BigInt(Math.ceil(x.toString(2).length/2));
 for(;;){const s=(t+x/t)>>1n;if(s>=t)return t;t=s;}
}
const ceil=(a,b)=>a<=0n?a/b:(a+b-1n)/b;
// v must be even, so the original right address c=v+1 is odd.
function invert(K,A,B){
 if(K<A+B)return {v:[],suffixes:0,reads:0};
 const gcd=eg(A,B)[0];
 if(K%gcd)return {v:[],suffixes:0,reads:0};
 const want=sqrt(gcd*K)/B;
 let k=1,ell=2n;while(ell<want){ell*=2n;k++;}
 const [g,j]=eg(B*ell,A);assert.equal(g,gcd);
 const M=A/g,t=j;
 const found=[];let reads=0;
 for(let r=0n;r<ell;r+=2n){
  const Z=(ell-1n)^reverse(r,k), eps=Z&1n;
  const w=B*ell+A*Z, C=B*(r+1n)+A*(F(r)-1n), D=K-C;
  if(D<A)continue;
  const high=(D-A)/w, low=ceil(D-2n*A,w+A);
  const start=low<0n?0n:low, phase=mod(t*(D/g),M);
  for(let q=start+mod(phase-start,M);q<=high;q+=M){
   reads++;
   const cnt=w*q+C+A*(F(q)+eps*(q&1n));
   if(cnt===K){const v=ell*q+r;assert.equal(B*(v+1n)+A*F(v),K);found.push(v);}
  }
 }
 return {v:found,suffixes:Number(ell/2n),reads};
}
function controls(){
 let checks=0;
 for(let A=1n;A<=17n;A+=2n)for(const B of [A,A+1n,2n*A]){
  for(let v=0n;v<1024n;v+=2n){
   const y=A*F(v)+B*(v+1n), r=invert(y,A,B);
   assert(r.v.includes(v));checks++;
  }
 }
 return checks;
}
const targets=[13,27,289,889,1489,9011,47184047,53084231,34379738367];
function search(max){
 const data=Object.fromEntries(targets.map(t=>[t,{hits:0,reads:0,suffixes:0,min:null,examples:[]}]));
 const seen=new Set();
 const g=(obj,v,k,e,n,T,j,A,B,D)=>{
   obj.hits++;
   const order=bank[j][0].length+n+k+1+Math.max(1,v===0n?0:v.toString(2).length);
   const row={n,T,j,c:String(v+1n),pad:k,Eflip:e,A:String(A),B:String(B),D:String(D),order};
   if(!obj.min||order<obj.min.order)obj.min=row;
   if(obj.examples.length<5)obj.examples.push(row);
 };
 const read=(N,n,T,j,A,B,D,k,e)=>{
  const obj=data[N], x=invert(BigInt(N)+D,A,B);
  obj.reads+=x.reads;obj.suffixes+=x.suffixes;
  for(const v of x.v)g(obj,v,k,e,n,T,j,A,B,D);
 };
 for(const [n,T] of roots(max)){
  const [h,e,r,V,Q]=inside(T);
  for(let j=0;j<bank.length;j++){
   const [,,,a,b,c]=bank[j],H=BigInt(a),R=BigInt(b),E=BigInt(c);
   const A=H*h,B=H*r+R*Q,D=H*e+R*V+E,D0=R*V+E;
   const key=[j,A,B,D,D0].join();if(seen.has(key))continue;seen.add(key);
   for(const N of targets)for(const sw of [0,1]){
    const z=BigInt(sw*2);
    read(N,n,T,j,A,B,D+z,0,sw);
    for(let k=1,W=2n*B;W<=BigInt(N)+D0+z;W*=2n,k++)read(N,n,T,j,A,W,D0+z,k,sw);
   }
  }
 }
 return {max,distinct:seen.size,data};
}
console.log(JSON.stringify({controls:controls(),...search(Number(process.argv[2]||8))},null,2));
