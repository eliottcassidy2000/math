"use strict";

const cores=[
 [[6,60,56,49,33,1],1,[0,3],23,24,8,[[1,11,2]]],
 [[12,53,56,50,33,1],1,[2,4],25,28,8,[[4,20,-2],[14,20,2]]],
 [[12,45,56,48,35,1],1,[3,5],29,31,9,[[14,28,-2],[0,28,2]]]
];
function flip(A,i,j){let B=A.slice();B[i]^=1<<j;B[j]^=1<<i;return B;}
function dp(A,root=-1){
 const n=A.length,M=(1<<n)-1,D=new Float64Array((M+1)*n);
 for(let i=0;i<n;i++)if(i===root||root<0)D[(1<<i)*n+i]=1;
 const P=new Float64Array(M+1);P[0]=1;
 for(let s=1;s<=M;s++)for(let j=0;j<n;j++){
  const z=D[s*n+j];if(!z)continue;P[s]+=z;
  for(let m=A[j]&~s;m;m&=m-1){const k=m&-m;D[(s|k)*n+31-Math.clz32(k)]+=z;}
 }
 return [P,D];
}
function profile(A,r){
 const [P,D]=dp(A),[E]=dp(A,r), n=A.length,M=P.length-1;
 let R=0;for(let s=0;s<=M;s++)if(s>>r&1)R+=D[s*n+r]*P[M^s];
 return [P[M],E[M],R];
}
function F(x){
 let f=1n;
 while(x){const q=x>>1n;if(!(x&1n))f+=q+(q&1n);x=q;}
 return f;
}
function reverse(x,i){
 let y=0n;
 for(let k=0;k<i;k++,x>>=1n)y=(y<<1n)|(x&1n);
 return y;
}
function seam(A,r,c,force=1){
 const n=A.length, v=c-1, m=Math.max(force,32-Math.clz32(v)), s=n+m;
 const B=[...A,...new Array(m+1).fill(0)];
 for(let i=n;i<=s;i++)for(let j=0;j<n;j++){
  if(i===s&&j===r)B[j]|=1<<i;else B[i]|=1<<j;
 }
 for(let i=n;i<s;i++){
  for(let j=n;j<i;j++)B[i]|=1<<j;
  if(v>>(i-n)&1)B[s]|=1<<i;else B[i]|=1<<s;
 }
 return B;
}
function strongly(A){
 const M=(1<<A.length)-1;
 for(let i=0;i<A.length;i++){
  let s=1<<i,o=0;
  while(s!==o){o=s;for(let j=0;j<A.length;j++)if(s>>j&1)s|=A[j];}
  if(s!==M)return false;
 }
 return true;
}
function count(c,h,R,E){return R*c+h*F(c-1n)-E*(c&1n);}
function expect(ok){if(!ok)throw Error("replay");}
function run(){
 const sample=[], big=[];
 for(const [A,u,[a,b],H,R0,Eu,classes] of cores){
  const K=flip(A,a,b), h=BigInt(H),R=BigInt(R0),E=BigInt(Eu);
  expect(profile(A,u).join()===[H,Eu,R0].join());
  expect(profile(K,u).join()===[H,Eu+2,R0].join());
  expect(strongly(A)&&strongly(K));
  for(const c of [1,3,5,9,17])for(const [G,e] of [[A,E],[K,E+2n]]){
   const B=seam(G,u,c),[P]=dp(B);expect(BigInt(P[P.length-1])===count(BigInt(c),h,R,e));
  }
  for(const [start,period,sig] of classes){
   const first=start||period,s=BigInt(sig);
   for(let i=first;i<=90;i+=period){
    const t=1n<<BigInt(i);
    expect((R*t-s)%h===0n);
    const echelon=[0n,2n*(t/8n),(t-2n)&~1n];
    for(const r of new Set(echelon)){
     const z=(t-1n)^reverse(r,i), beta=(R*t-s)/h+z+(z&1n);
     expect((beta&1n)===0n);
     for(const q of [beta,beta-1n]){
      const c=2n*t*q+r+1n;
      const e=s===2n?E:E+2n, N=count(c,h,R,e);
      const total=[N,count(c,h,R,e+s),count(c+t,h,R,e),count(c+t,h,R,e+s)];
      expect(total.join()===[N,N-s,N+s,N].join());
      if(H===23&&i===1&&r===0n){
       const B=seam(A,u,Number(c)),T=flip(B,a,b), W=flip(B,A.length+i,B.length-1);
       const ls=[B,T,W,flip(W,a,b)].map(g=>{expect(strongly(g));const [P]=dp(g);return BigInt(P[P.length-1]);});
       expect(ls.join()===total.join());
       sample.push({c:Number(c),size:B.length,counts:ls.map(Number)});
      }
      if(H===23&&i===23&&r===0n&&q===beta){
       expect((N-1n)%9n===1n);
       big.push({i,c:String(c),N:String(N),lowerMidMod9:String((N-1n)%9n)});
      }
     }
    }
   }
  }
 }
 return {sample,big};
}function samples(){
 const [h,R,E]=[23n,24n,8n],m=[9n,27n,81n,7n,11n,63n,99n,231n,135n];
 const rows=[];
 for(const i of [1,12,23]){
  const l=1n<<BigInt(i),set=m.map(()=>new Set()),heights=new Set();
  const limit=i===23?4096n:l;
  for(let r=0n;r<limit;r+=2n){
   const z=(l-1n)^reverse(r,i),beta=(R*l-2n)/h+z+(z&1n);
   for(const q of [beta,beta-1n]){
    const c=2n*l*q+r+1n,N=count(c,h,R,E);
    heights.add(N);
    m.forEach((a,k)=>set[k].add(Number(N%a)));
   }
  }
  rows.push([i,heights.size,m.map((a,k)=>[Number(a),set[k].size])]);
 }
 return rows;
}

function coverTest(){
 let tested=0;
 for(const [T,u] of cores){
  const n=T.length, [P,D]=dp(T), [A]=dp(T,u), full=P.length-1;
  let V=0,Q=0;
  for(let s=0;s<=full;s++){
   Q+=P[s]*P[full^s];if(s>>u&1)V+=A[s]*P[full^s];
  }
  for(let c=1;c<=64;c++){
   const B=seam(T,u,c),r=B.length-1,M=(1<<B.length)-1,[S,U]=dp(B),[Z]=dp(B,r);
   let G=0;for(let s=1<<r;s<=M;s++)G+=Z[s]*S[M^s];
   expect(G===c*Q-V&&U[M*B.length+r]===D[full*n+u]);
   tested++;
  }
 }
 return tested;
}

if(require.main!==module){module.exports={cores,flip,dp,F,reverse,seam,expect};}else{
const result={seams:run(),residues:samples(),grammarChecks:coverTest()};
expect(result.grammarChecks===192);
expect(result.residues[1][1]===4096);
expect(result.residues[1][2].every(([m,n])=>m===n));
console.log(JSON.stringify(result,null,2));

}
