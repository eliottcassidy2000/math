// node replay_prime_quotient.js [positive target]
"use strict";
const Q5=[28,25,10,16,4], Q6=[60,9,42,16,6,26];
const bit=(n,i)=>(n>>i)&1;
const check=(q,t)=>{if(!q)throw Error(t)};
function primeQ(Q){
 let full=(1<<Q.length)-1;
 for(let S=1;S<full;S++){
  if((S&(S-1))&&Q.every((r,i)=>bit(S,i)||(r&S)===0||(r&S)===S))return false;
 }
 return true;
}
function runTable(Q,limits){
 const n=Q.length, a=Array(n).fill(0), table={};
 function search(last){
  if(a.every(x=>x>0)){const id=a.join("");table[id]=(table[id]??0)+1;}
  for(let i=0;i<n;i++)if(a[i]<limits[i]&&(last===-1||bit(Q[last],i))){
    a[i]++;search(i);a[i]--;
  }
 }
 search(-1);return table;
}
function TT(n){return Array.from({length:n},(_,i)=>((1<<n)-1)^((2<<i)-1));}
function blow(Q, sizes){
 const starts=[0];for(const n of sizes)starts.push(starts.at(-1)+n);
 const out=[];
 for(let i=0;i<Q.length;i++){
  let cross=0;
  for(let j=0;j<Q.length;j++)if(bit(Q[i],j))cross|=((1<<sizes[j])-1)<<starts[j];
  for(const r of TT(sizes[i]))out.push((r<<starts[i])|cross);
 }
 return out;
}
function count(rows){
 const n=rows.length,F=(1<<n)-1,pre=rows.map((_,i)=>rows.reduce((s,r,j)=>s|(bit(r,i)<<j),0));
 const dp=Array.from({length:F+1},()=>Array(n).fill(0));
 for(let i=0;i<n;i++)dp[1<<i][i]=1;
 for(let S=1;S<=F;S++)for(let i=0;i<n;i++)if(bit(S,i)){
  let s=pre[i]&S;
  while(s){const b=s&-s;s^=b;dp[S][i]+=dp[S^(1<<i)][31-Math.clz32(b)];}
 }
 return dp[F].reduce((s,v)=>s+v,0);
}
const table3={"11111":9,"11112":5,"11121":3,"11122":3,"21111":3,"21112":2,"21121":1,"21122":1};
const table2={"11111":9,"11112":5,"12111":5,"12112":4};
const table6={"111111":31,"111211":14,"121111":33,"121211":21,"131111":8,"131211":8};
function same(a,b){return Object.keys(a).length===Object.keys(b).length&&Object.keys(b).every(k=>a[k]===b[k]);}
check(primeQ(Q5)&&primeQ(Q6),"quotients must have no nontrivial proper module");
check(same(runTable(Q5,[3,1,1,3,3]),table3),"three variable slots");
check(same(runTable(Q5,[1,4,1,1,4]),table2),"two variable slots");
check(same(runTable(Q6,[1,5,1,5,1,1]),table6),"six-quotient return");
const cases=[[Q5,[2,1,1,2,2],63],[Q5,[1,2,1,1,3],97],[Q6,[1,2,1,4,1,1],881]];
for(const [Q,V,N] of cases)check(count(blow(Q,V))===N,"literal target "+N);
for(let a=1;a<=3;a++)for(let b=1;b<=3;b++)for(let c=1;c<=3;c++){
 const V=[a,1,1,b,c], m=2**c;
 check(count(blow(Q5,V))+m-2===(2**a+1)*(2**b*(m-1)+1),"three-width identity");
}
function logpow(x){
 if(x<2n||(x&(x-1n)))return 0;
 return x.toString(2).length-1;
}
function decode(N){
 if(N<3n)return [];
 const hits=[];
 for(let c=1,M=2n;M<N;M*=2n,c++){
   const D=N+M-2n;
   // Every candidate first factor is a Fermat-type number 1+2^a.
   // Enumerate its exponent rather than factor D at large heights.
   for(let a=1,U=3n;U<D;U=(U-1n)*2n+1n,a++){
      if(D%U!==0n)continue;
      const V=D/U-1n;
      if(V%(M-1n)!==0n)continue;
      const b=logpow(V/(M-1n));if(b)hits.push([a,1,1,b,c]);
   }
 }
 return hits;
}
check(decode(63n).some(x=>x.join()==="2,1,1,2,2"),"target-first decode for 63");
console.log("PASS: prime quotients and bounded run tables; 63, 97, 881 literal; 27 exact three-width tests");
if(typeof process!=="undefined"&&process.argv[2]){
 const N=BigInt(process.argv[2]);console.log(JSON.stringify(decode(N)));
}
