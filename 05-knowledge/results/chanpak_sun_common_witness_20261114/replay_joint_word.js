// node replay_joint_word.js [modulus H R R_converse parity]
"use strict";
const check=(c,m)=>{if(!c)throw Error(m)};
const mod=(a,b)=>(a%b+b)%b;
const pow=(n)=>1n<<BigInt(n);
function reverse(x,n){let y=0n;for(let j=0;j<n;j++){y=2n*y+(x&1n);x>>=1n}return y;}
function g(x){let h=1n;while(x){const b=x&1n;x>>=1n;if(!b)h+=x+(x&1n)}return h;}
function params(M,par){
 let m=M,k=0;while(!(m&1n)){m>>=1n;k++;}
 let ord=1;for(let a=2n%m;m!==1n&&a!==1n;a=2n*a%m)ord++;
 const q=ord*Math.ceil(Math.max(2,k)/ord), b=Math.max(1,k);
 const ell=b+Math.max(1,k)+q+k+2;
 const ds=[...Array(Number(m-1n)).fill(1),...Array.from({length:Math.max(k-1,0)},(_,i)=>i+2)];
 const pairs=ds.map(d=>{
  const low=pow(b)+pow(b+d),high=pow(b)+pow(b+d+q);
  return {d,low,high,delta:g(high)-g(low)};
 });
 const anchors=Math.max(1,k), W=q*Number(m);
 const len=2*anchors+W+ell*pairs.length;
 const pad=(par-len%2+2)%2,n=len+pad;
 return {m,k,q,ell,pairs,anchors,W,pad,n};
}
function build(M,ht,r,rv,par){
 if(M<1n||!(par===0||par===1))throw Error("require positive modulus and parity 0/1");
 const {m,k,q,ell,pairs,anchors,W,pad,n}=params(M,par);
 if(k&&(!(ht&1n)||!(r&1n)||!(rv&1n)))return null;
 if(m%3n===0n&&mod(rv+(par===1?1n:-1n)*r,3n)!==0n)return null;
 let P=k?mod(r-1n,pow(k)):0n;
 const leadingRev=k?mod(pow(n)-rv,pow(k)):1n;
 P|=reverse(leadingRev,anchors)<<BigInt(n-anchors);
 const offset=anchors+W+pad;
 for(let j=0;j<pairs.length;j++)P|=pairs[j].low<<BigInt(offset+j*ell);
 if(m!==1n){
  const c=pow(n-1)%m;
  const A=mod(r-1n-P,m), B=mod(pow(n)-rv-reverse(P,n),m);
  let b2=0n;while(mod(c*b2-B,m)!==0n)b2++;
  let y=0n;while(y<m&&mod(3n*y-2n*(A-b2),m)!==0n)y++;
  check(y<m,"compatibility of odd controls");
  let x=mod(A-2n*y,m);
  for(let j=anchors;j<anchors+W;j++){
   if(j%q===0&&x){P|=pow(j);x--;}
   if(j%q===1&&y){P|=pow(j);y--;}
  }
  check(x===0n&&y===0n,"enough phase controls");
 }
 let H=g(P);
 function select(j){
  const {low,high,delta}=pairs[j];
  P^=(low^high)<<BigInt(offset+j*ell);
  H+=delta;
 }
 const rOdd=Number(mod(ht-H,m));
 for(let j=0;j<rOdd;j++)select(j);
 if(k){
  const remaining=mod(ht-H,pow(k));
  for(let s=1;s<k;s++)if(remaining&pow(s))select(Number(m)-1+s-1);
 }
 H=g(P);
 const R=P+1n,Rc=pow(n)-reverse(P,n);
 check(P>=pow(n-1)&&P<pow(n)&&!(P&1n),"strong endpoints");
 check(mod(H-ht,M)===0n&&mod(R-r,M)===0n&&mod(Rc-rv,M)===0n,"simultaneous output");
 return {n,H,R,Rc,P,word:P.toString(2)};
}
function literal(P,n){
 const rows=Array.from({length:n},(_,i)=>(pow(n)-pow(i+1))|((1n-((P>>BigInt(n-1-i))&1n))<<BigInt(n)));
 rows.push(reverse(P,n));
 const N=n+1,F=(1<<N)-1;
 const ends=Array.from({length:F+1},()=>Array(N).fill(0));
 for(let i=0;i<N;i++)ends[1<<i][i]=1;
 for(let s=1;s<=F;s++)for(let j=0;j<N;j++)if((s>>j)&1){
  for(let i=0;i<N;i++)if((rows[i]>>BigInt(j))&1n)ends[s][j]+=ends[s^(1<<j)][i];
 }
 return BigInt(ends[F].reduce((a,b)=>a+b,0));
}
let pass=0;
for(let v=1n;v<=7n;v++)for(let p=0;p<2;p++){
 for(let h=0n;h<v;h++)for(let r=0n;r<v;r++)for(let b=0n;b<v;b++){
  const x=build(v,h,r,b,p);
  const allowed=(v%2n!==0n||!!(h&r&b&1n))&&(v%3n!==0n||mod(b+(p?1n:-1n)*r,3n)===0n);
  check(!!x===allowed,"classification");pass+=!!x;
  if(x&&x.n<=9)check(literal(x.P,x.n)===x.H,"literal augmented count");
 }
}
for(const v of [2n,3n,4n,5n,6n,8n,9n,12n,16n,25n,30n,35n,42n]){
 for(const p of [0,1]){
  const x=params(v,p);
  for(const z of x.pairs){
   check(mod(z.high-z.low,v)===0n&&mod(reverse(z.high,x.ell)-reverse(z.low,x.ell),v)===0n,"both block addresses match");
   check(mod(z.delta,x.m)===(z.d===1?1n%x.m:0n),"odd part selector");
   if(z.d>1)check(mod(z.delta,pow(x.k))===pow(z.d-1),"binary selector");
   if(x.ell<=12)check(literal(z.low,x.ell)===g(z.low)&&literal(z.high,x.ell)===g(z.high),"literal switch pair");
  }
  for(let a=1n;a<=15n;a++)build(v,2n*a+1n,6n*a+1n,6n*a+(p?5n:1n),p);
 }
}
console.log("PASS: "+pass+" compatible simultaneous triples through modulus 7, plus composite-modulus and literal controls");
if(typeof process!=="undefined"&&process.argv[2]){
 const [M,H,R,Rc]=process.argv.slice(2,6).map(BigInt);
 const res=build(M,H,R,Rc,Number(process.argv[6]));
 console.log(JSON.stringify(res,(_,v)=>typeof v==="bigint"?v.toString():v,2));
}
