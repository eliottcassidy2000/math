"use strict";
function ck(x,msg){if(!x)throw Error(msg);}

function h(a){let n=a.length,N=1<<n,dp=new Int32Array(N*n);for(let i=0;i<n;i++)dp[(1<<i)*n+i]=1;for(let s=1;s<N;s++)for(let i=0;i<n;i++){let w=dp[s*n+i];if(!w)continue;let b=a[i]&(N-1-s);while(b){let x=b&-b;b-=x;dp[(s|x)*n+31-Math.clz32(x)]+=w}}let z=0;for(let i=0;i<n;i++)z+=dp[(N-1)*n+i];return z}
function graph(n,p){let a=Array(n).fill(0);for(let j=1;j<n;j++)for(let i=0;i<j;i++){if(p&1)a[i]|=1<<j;else a[j]|=1<<i;p>>>=1;}return a}
function blue(n){let orb=[];for(let i=0;i<n;i++)for(let j=i+2;j<n;j++){let k=n-1-j,l=n-1-i;if(i*n+j<=k*n+l)orb.push(i==k?[[i,j]]:[[i,j],[k,l]])}let frames=[];for(let m=0;m<(1<<orb.length);m++){let a=Array.from({length:n},(_,i)=>i+1<n?1<<(i+1):0);orb.forEach((o,u)=>o.forEach(([i,j])=>{if((m>>u)&1)a[i]|=1<<j;else a[j]|=1<<i}));frames.push(a)}return frames}
function mate(a){let n=a.length;return a.map((x,i)=>((1<<n)-1)^x^(1<<i)^(i>0?1<<(i-1):0)^(i+1<n?1<<(i+1):0));}

const wit={"5":[4,9,10,1],"11":[2,12,25,17,3],"17":[58,52,1,54,4,20],"19":[54,52,32,7,44,8],"29":[48,53,17,7,8,28],"41":[24,37,17,38,10,21],"47":[72,53,1,54,37,5,62],"49":[86,52,80,87,64,29,34],"51":[50,16,59,19,64,26,47],"57":[124,13,24,80,66,94,6],"59":[106,92,81,20,1,94,24],"67":[22,96,82,55,2,21,57],"71":[122,104,99,52,38,64,24],"97":[34,96,27,3,11,28,61],"179":[190,8,98,164,14,18,59,118],"239":[218,120,123,112,192,17,160,46],"397":[16,141,41,97,142,83,23,109],"617":[164,237,360,433,423,128,441,260,35]};
const edge=[[110,100,32,6,79,88,12],[110,100,96,6,79,88,8],[6,2]];
function strong(a){
 for(let i=0;i<a.length;i++){let b=1<<i;for(let _ of a)for(let j=0;j<a.length;j++)if(b>>j&1)b|=a[j];
  if(b!==(1<<a.length)-1)return false;}
 return true;
}
for(const [n,a] of Object.entries(wit)){
 for(let i=0;i<a.length;i++){ck(!(a[i]>>i&1)&&a[i]<(1<<a.length),"matrix");
  for(let j=0;j<i;j++)ck(((a[i]>>j)&1)+((a[j]>>i)&1)===1,"arc");}
 ck(h(a)===Number(n),"count "+n);
}
ck(179*617===110443&&5*397===1985,"ordinal products");
ck(h(edge[0])===49&&h(edge[1])===51&&edge[0].every((r,i)=>r===edge[1][i]||(edge[2].includes(i)&&r===(edge[1][i]^(1<<edge[2].find(j=>j!==i))))),"49/51 arc");
const expected={ "0,1,2,3":[1,1],"1,1,1,3":[3,3],"0,2,2,2":[3,3],"1,1,2,2":[5,1] };
let seen={};
for(let k=0;k<64;k++){
 let a=graph(4,k);const p=(i,j)=>a[i]>>j&1?1:-1;
 let pf=p(0,1)*p(2,3)-p(0,2)*p(1,3)+p(0,3)*p(1,2);
 let sc=String(a.map(v=>v.toString(2).replaceAll("0","").length).sort());
 ck(String([h(a),Math.abs(pf)])===String(expected[sc]),"four classes");
 seen[sc]=(seen[sc]||0)+1;
}
ck(Object.keys(seen).length===4,"all classes");
function TT(n){return Array.from({length:n},(_,i)=>(1<<n)-(1<<(i+1)));}
function host(D,q,a,t){
 let tag=a.toString(2)+"0".repeat(t),n=tag.length,d=D.length,R=TT(n).map((r,i)=>r|(((1<<d)-1)<<n)|((1-+tag[i])<<(n+d)));
 return R.concat(D.map((v,i)=>(v<<n)|((1-(q>>i&1))<<(n+d))),
                  [(q<<n)|tag.split("").reduce((r,b,i)=>r+(+b<<i),0)]);
}
for(let a of [3,5]) for(let t=1;t<=6;t++){
 let x=host([6,4,0],2,a,t), y=host([2,4,1],6,a,t);
 ck(strong(x)&&strong(y)&&h(x)===6*a*2**t+1&&h(y)===6*a*2**t-1,"two-arc counts");
 const diff=x.map((r,i)=>(r^y[i]).toString(2).replaceAll("0","").length).reduce((x,y)=>x+y,0);
 ck(diff===4,"two arcs");
}

const mul=([a,b],[c,d])=>[a*c-b*d,a*d+b*c];
const bar=([a,b])=>[a,-b];
function pw(z,n){let out=[1n,0n];for(let i=0;i<n;i++)out=mul(out,z);return out;}
function div([a,b],[c,d]){
 let e=c*c+d*d,x=a*c+b*d,y=b*c-a*d;
 if(x%e||y%e)return null;
 return [x/e,y/e];
}
const supports=[[1n,2n],[2n,3n],[24n,25n]];
for(const [n,w] of [[49,[0,0,1]],[57,[-3,1,0]],[239,[0,4,0]],[110443,[8,-1,-1]]]){
 let g=div([BigInt(n),1n],[1n,1n]);ck(!!g,"ramified prime");
 for(let i=0;i<w.length;i++){let v=w[i];g=div(g,pw(v<0?bar(supports[i]):supports[i],Math.abs(v)));ck(!!g,"factor");}
 ck(g[0]*g[0]+g[1]*g[1]===1n,"no other support");
}
function gcd(a,b){while(b){let t=a%b;a=b;b=t;}return a<0n?-a:a;}
const F=(a,b=1n)=>{let g=gcd(a,b);return [a/g,b/g];},add=([a,b],[c,d])=>F(a*d+b*c,b*d),lt=([a,b],[c,d])=>a*d<c*b, times=([a,b],c)=>F(a*c,b);
function identity(pairs){
 let A=[1n,0n],B=[1n,0n],lo=F(0n),hi=F(0n);
 for(const [n0,k0] of pairs){
  const n=BigInt(n0), k=BigInt(k0), s=F(1n,n),t=add(s,F(-1n,3n*n**3n));
  lo=add(lo,times(k>0?t:s,k));hi=add(hi,times(k>0?s:t,k));
  if(k>0)A=mul(A,pw([n,1n],k0));else B=mul(B,pw([n,1n],-k0));
 }
 const Z=mul(A,bar(B));ck(Z[0]>0n&&Z[0]===Z[1],"diagonal");
 ck(lt(F(1n,2n),lo)||String(lo)===String(F(1n,2n)),"phase below");
 ck(lt(hi,F(1n)),"phase above");
}
for(const f of [[[2,1],[3,1]],[[2,2],[7,-1]],[[3,2],[7,1]],[[5,4],[239,-1]],[[4,3],[20,1],[1985,1]],[[49,12],[57,32],[239,-5],[110443,12]]])identity(f);
const P=new Int32Array(8193), non=new Uint8Array(8193);
for(let i=2;i<P.length;i++){
 if(!non[i])for(let j=i*i;j<P.length;j+=i)non[j]=1;
 P[i]=P[i-1]+(non[i]===0?1:0);
}
const D=x=>P[x]-P[x>>1];
const R=n=>{for(let x=P.length-1;x>=0;x--)if(D(x)<n)return x+1;};
ck(String(Array.from({length:10},(_,i)=>R(i+1)))==="2,11,17,29,41,47,59,67,71,97","Ramanujan list");
ck(R(17)===179&&R(22)===239&&R(23)===241&&R(49)===607&&R(50)===641&&D(617)===50&&D(626)===49,"persistent cutoffs");
ck(239**2+1===57122&&2*57122===114244&&239%4===3&&57122%4===2,"Paley sizes");
console.log("PASS: all six angle identities, Gaussian coefficients, four classes, literal tournament controls and Ramanujan cutoffs.");
