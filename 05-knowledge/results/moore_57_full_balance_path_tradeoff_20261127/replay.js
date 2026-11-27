
"use strict";
const assert=p=>{if(!p)throw Error("check failed")};
const range=n=>Array.from({length:n},(_,i)=>i);
const sq=A=>A.map((r,i)=>r.map((_,j)=>r.reduce((s,x,k)=>s+x*A[k][j],0)));
const mod=(a,p)=>(a%p+p)%p;
const f3=a=>mod(a,3);
const digits=x=>[x%3,((x/3)|0)%3,(x/9)|0];
function diff27(a,b){
 const x=digits(a),y=digits(b);
 return x.reduce((s,v,i)=>s+3**i*f3(v-y[i]),0)
}
function mul27(a,b){
 const x=digits(a),y=digits(b),z=[0,0,0,0,0];
 for(let i=0;i<3;i++)for(let j=0;j<3;j++)z[i+j]+=x[i]*y[j];
 for(let i=4;i>=3;i--){z[i-3]+=z[i];z[i-2]+=z[i]}
 return z.slice(0,3).reduce((s,v,i)=>s+3**i*f3(v),0)
}
const R=new Set(range(26).map(i=>mul27(i+1,i+1)));
assert(R.size===13);
const P=range(28).map(i=>range(28).map(j=>i===j?0:i===27?1:j===27?-1:
 R.has(diff27(j,i))?1:-1));
const C=range(56).map(i=>range(56).map(j=>{
 const a=i%28,b=j%28,u=P[a][b],d=+(a===b);
 if(i<28&&j<28)return u;
 if(i>=28&&j>=28)return -u;
 return u+(i<28?d:-d)
}));
const check=(A,x)=>assert(sq(A).every((r,i)=>r.every((v,j)=>v===(i===j?x:0))));
check(P,-27);check(C,-55);
function scores(A){return A.map(r=>r.reduce((s,x)=>s+x,0))}
function maskSwitch(A,mask){
 return A.map((r,i)=>r.map((x,j)=>x*(mask>>BigInt(i)&1n?-1:1)*(mask>>BigInt(j)&1n?-1:1)))
}

const scaledC=maskSwitch(C,0x1cd175b1704278n);
check(scaledC,-55);
const norm=C=>Math.max(...scores(C).map(Math.abs));
assert(norm(scaledC)===11);
console.log("Conference 28 and doubled 56 verified; explicit switched 56 has max row sum",norm(scaledC));
function paley(p){
 const r=new Set(range(p-1).map(i=>(i+1)**2%p));
 return range(p).map(i=>range(p).map(j=>i===j?0:r.has(mod(j-i,p))?1:-1));
}
const P59=paley(59), Q=P59.slice(0,57).map(r=>r.slice(0,57)), D=P59.slice(0,56).map(r=>r.slice(0,56));
const h={};
const G=sq(D);
for(let i=0;i<56;i++)for(let j=i+1;j<56;j++)h[G[i][j]]=(h[G[i][j]]||0)+1;
assert(norm(D)===3 && h[-2]===196 && h[0]===588 && h[2]===588 && h[4]===168);
console.log("Paley deletion to 56: row-sum bound",norm(D),"square census",h);
assert(norm(Q)===2&&sq(Q).every((r,i)=>r.every((x,j)=>i===j||Math.abs(x)<=3)));
const ceil=(x,y)=>(x+y-1n)/y;
function pack(n,R){
 let fac=1n,prod=1n,pow=1n,H=1n;
 for(let t=1;n-3*(t-1)>=3;t++){
  const s=n-3*(t-1),b=Math.min(s-1,R+n-s);
  const num=s*(s*s-1)-3*(s-s%2)*b*b;
  if(num<=0)break;
  prod*=ceil(BigInt(num),24n);fac*=BigInt(t);pow*=2n;
  H+=pow*ceil(prod,fac);
 }
 return H
}
function hp(A) {
 const n=A.length,N=1<<n,table=new Int32Array(N*n),s=new Int32Array(n);
 for(let i=0;i<n;i++)table[(1<<i)*n+i]=1;
 for(let m=1;m<N;m++)for(let i=0;i<n;i++)if(table[m*n+i])
  for(let j=0;j<n;j++)if(A[i][j]===1&&!(m>>j&1))table[(m|1<<j)*n+j]+=table[m*n+i];
 return BigInt(table.slice(-n).reduce((a,b)=>a+b));
}
function mk(n,b){
 let A=range(n).map(()=>range(n).map(()=>0)),k=0;
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++){A[i][j]=b>>k++&1?-1:1;A[j][i]=-A[i][j]}
 return A
}
for(let n=3;n<=6;n++)for(let b=0;b<1<<n*(n-1)/2;b++){
 const A=mk(n,b);assert(pack(n,norm(A))<=hp(A));
}
const L=pack(56,3),Lq=pack(57,2),J=Lq*L*ceil(L,2n)**56n;
assert(L===10089842422501174206273n&&Lq===60601093599556238992453n&&J>10n**1260n);
console.log("Branch and quotient triangle-packing bounds",String(L),String(Lq));
console.log("Contiguous path lower bound digits and prefix",String(J).length,String(J).slice(0,25));
function cap(g,r){return 52*g+4+2*r+4}
assert(cap(55,11)===2890&&cap(3,11)===186&&cap(3,3)===170);
console.log("Conference witness cross caps: any quotient 2890; Paley benchmark 186.");
