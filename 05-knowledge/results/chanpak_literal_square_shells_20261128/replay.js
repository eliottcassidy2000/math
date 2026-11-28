"use strict";
const assert=require("node:assert/strict");
const core=[16,5,9,19,38,15], tags=35;
const MAX=2**21;
const cache=new Uint32Array(4*MAX);
for(let n=2;n<cache.length;n++){
  const m=Math.floor(n/2);
  cache[n]=cache[m]+(n%2?0:m+m%2);
}
function f(n){if(n<cache.length)return cache[n];const m=Math.floor(n/2);return f(m)+(n%2?0:m+m%2);}
function L(x,y){return 111+129*x+162*y+320*x*y+29*(f(x)+f(y));}
function flipLower(p,a,b){
  const x=2*p-a,y=2*p+2-b;return 2*x*y+2*p+f(x)+f(y)+1;
}
const cert=(h,L)=>BigInt(L)>h*h && BigInt(L)+3n<(h+1n)**2n;
const good=(h,l)=>l>h*h&&l+3<(h+1)**2;
function masks(n,m){
 const a=Array(n).fill(0);let k=0;
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++,k++)
   if(m>>k&1)a[i]|=1<<j;else a[j]|=1<<i;
 return a;
}
function build(x,y,flipped){
 const left=x?x.toString(2):"";
 const right=y?[...y.toString(2)].reverse().map(c=>c==="0"?"1":"0").join(""):"";
 const a=left.length+6+right.length, rows=Array(a+1).fill(0), o=left.length;
 const isTag=i=>i<o?left[i]==="1":i<o+6?Boolean(tags>>(i-o)&1):right[i-o-6]==="1";
 for(let i=0;i<a;i++){
   if(isTag(i)) rows[a]|=1<<i; else rows[i]|=1<<a;
   for(let j=i+1;j<a;j++){
     let forward=i<o||j>=o+6||Boolean(core[i-o]>>(j-o)&1);
     if(flipped&&i===o+1&&j===o+4) forward=!forward;
     if(forward) rows[i]|=1<<j; else rows[j]|=1<<i;
   }
 }
 return rows;
}
function strong(rows){
 const n=rows.length;
 for(let i=0;i<n;i++){
   let s=1<<i;
   for(let k=0;k<n;k++)for(let j=0;j<n;j++)if(s>>j&1)s|=rows[j];
   assert.equal(s,(1<<n)-1);
 }
}
function pathCounts(rows){
 const n=rows.length,S=1<<n,dp=new Float64Array(n*S);
 for(let i=0;i<n;i++)dp[(1<<i)*n+i]=1;
 for(let m=1;m<S;m++){
   let z=m;
   while(z){
     const bit=z&-z,i=31-Math.clz32(bit),u=dp[m*n+i];z^=bit;
     if(!u)continue;
     let s=rows[i]&~m;
     while(s){const b=s&-s;s^=b;dp[(m|b)*n+31-Math.clz32(b)]+=u;}
   }
 }
 const H=new Float64Array(S);H[0]=1;
 for(let m=1;m<S;m++)for(let i=0;i<n;i++)H[m]+=dp[m*n+i];
 return {H,dp};
}
function converse(rows){
 const a=Array(rows.length).fill(0);
 for(let i=0;i<rows.length;i++)for(let j=0;j<rows.length;j++)
   if(rows[i]>>j&1)a[j]|=1<<i;
 return a;
}
function profile(rows){
 const G=63,tag=35,non=G^tag;
 const H=pathCounts(rows),R=pathCounts(converse(rows));
 let q=0,c=0,t=0,e=0,r=0;
 for(let m=0;m<=G;m++){
   q+=H.H[m]*H.H[G^m];
   for(let v=0;v<6;v++){
     if(tag>>v&1)c+=R.dp[6*m+v]*H.H[G^m];
     if(non>>v&1)t+=H.dp[6*m+v]*H.H[G^m];
   }
 }
 for(let v=0;v<6;v++){
   if(tag>>v&1)e+=H.dp[6*G+v];else r+=R.dp[6*G+v];
 }
 return [H.H[G],q,c,t,e,r];
}
strong(core);
const changed=core.slice();changed[1]^=1<<4;changed[4]^=1<<1;
strong(changed);
assert.deepEqual(profile(core),[29,320,162,129,16,11]);
assert.deepEqual(profile(changed),[29,320,162,129,18,13]);
function checkLiteral(a,b,ha,hb){
 assert.equal(a.length,b.length);
 assert.equal(a.reduce((s,r,i)=>s+((r^b[i])?1:0),0),2);
 assert.equal(a.reduce((s,r,i)=>s+((r^b[i]).toString(2).match(/1/g)?.length||0),0),2);
 for(const [r,h] of [[a,ha],[b,hb]]){
   strong(r);assert.equal(pathCounts(r).H.at(-1),h);
 }
}
for(const [x,y] of [[0,0],[2,0],[0,2],[2,4],[20,38]]){
 const l=L(x,y); checkLiteral(build(x,y,true),build(x,y,false),l,l+2);
}
assert.equal(L(20,38),253265);
for(const [h,n,m,z,l] of [[5,6,80,336,27],[25,8,98750095,98749967,657]]){
 checkLiteral(masks(n,m),masks(n,z),l,l+2);assert(good(h,l));
}
function trial(h){
 for(const Y of [4,2,8,16]){
   const A=320*Y*Y,B=158*Y-29,C=191*Y-29,p=C/A,q=B/A,E=111-B*C/A;
   for(const c of [0.6,0.3,0.9,1.2,1.5,0.1]){
     const R=(h*h+c*h-E)/A,m0=Math.ceil(2*Math.sqrt(R)-p-q);
     for(let m=m0;m<m0+3;m++){
       const s=m+p+q;
       for(const sign of [-1,1]){
         const a0=Math.round(s/2-p+sign*Math.sqrt(Math.max(0,s*s/4-R)));
         for(const i of [0,-1,1,-2,2]){
           const a=a0+i,b=m-a;
           if(a>0&&b>0&&good(h,L(Y*a,Y*b)))return [Y*a,Y*b];
         }
       }
     }
   }
 }
 return null;
}
function finite(){
 const P=741461, values=new Float64Array(4*P);
 for(let p=1;p<=P;p++)for(let a=0;a<=1;a++)for(let b=0;b<=1;b++)
   values[4*(p-1)+2*a+b]=flipLower(p,a,b);
 values.sort();assert(flipLower(P+1,1,1)>MAX**2);
 let idx=0;
 const tally=[0,0,0,0], left=[];
 for(let h=5;h<MAX;h+=2){
   while(values[idx]<=h*h)idx++;
   if(good(h,values[idx])){tally[0]++;continue;}
   if(trial(h)){tally[1]++;continue;}
   left.push(h);
   if(h===5||h===25){tally[3]++;continue;}
   let y=2;while(111+162*y<(h+1)**2 && !good(h,L(0,y)))y+=2;
   assert(good(h,L(0,y)),`failed parent ${h}`);tally[2]++;
 }
 assert.deepEqual(tally,[1026067,22484,21,2]);
 assert.deepEqual(left,[5,25,271,469,1425,1657,2019,3651,4327,6075,6239,7195,7427,8103,9223,13525,14877,15519,20927,21555,22839,23467,24751]);
 return tally;
}
function isqrt(n){
 assert(n>=0n);if(n<2n)return n;
 let x=1n<<BigInt(Math.ceil(n.toString(2).length/2));
 for(;;){const y=(x+n/x)/2n;if(y>=x)return x;x=y;}
}
function ceilrt(n,d){const x=isqrt(n/d);return x*x*d===n?x:x+1n;}
function bf(n){
 let u=0n,z=0n;
 for(const c of n.toString(2)){const d=BigInt(c);if(!d)z+=u+u%2n;u=2n*u+d;}
 return z;
}
function bL(x,y){return 111n+129n*x+162n*y+320n*x*y+29n*(bf(x)+bf(y));}
function decode(t,d,Y){
 const A=320n*Y*Y,B=158n*Y-29n,C=191n*Y-29n;
 const M=A*(t-111n*d)+B*C*d;
 const m=(ceilrt(4n*M,d)-B-C+A-1n)/A;
 const S=A*m+B+C;
 const a=(S-2n*C+A-ceilrt(S*S*d-4n*M,d))/(2n*A);
 const b=m-a;assert(a>0n&&b>0n);
 return [Y*a,Y*b];
}
function powerWindow(N){
 const Y=1n<<BigInt(Math.floor((N.toString(2).length-1)/10));
 const goal=N*N;
 let a=0n,b=1n<<BigInt(Math.ceil(goal.toString(2).length/5));
 while(b-a>1n){const m=(a+b)/2n;if(m**5n>=goal)b=m;else a=m;}
 const xy=decode(N+256n*b,1n,Y); const L=bL(...xy);
 assert(L>N&&(L+2n-N)**5n<512n**5n*N**2n);
 return xy;
}
for(const h of [BigInt(MAX)+1n,10n**10n+1n,10n**100n+1n]){
 assert(cert(h,bL(...decode(5n*h*h+3n*h,5n,4n))));
}
for(const N of [10n**10n,10n**11n+9n,10n**100n])powerWindow(N);
console.log({finite:finite(),core503:[L(20,38),L(20,38)+2],exactLargeDecoders:"passed"});
