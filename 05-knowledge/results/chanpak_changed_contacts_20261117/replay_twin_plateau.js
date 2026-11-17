function need(x,msg){if(!x)throw Error(msg);}
function pathTable(rows,root=-1){
  const n=rows.length,F=(1<<n)-1;
  for(let i=0;i<n;i++){need(!(rows[i]>>i&1),"loop");for(let j=0;j<i;j++)need((rows[i]>>j&1)+(rows[j]>>i&1)===1,"not tournament");}
  const d=new Float64Array((F+1)*n), h=new Float64Array(F+1);h[0]=1;
  for(let i=0;i<n;i++)if(root<0||root===i)d[(1<<i)*n+i]=1;
  for(let s=1;s<=F;s++)for(let i=0;i<n;i++){
    const v=d[s*n+i];if(!v)continue;h[s]+=v;
    for(let u=rows[i]&~s;u;u&=u-1){const b=u&-u,j=31-Math.clz32(b);d[(s|b)*n+j]+=v;}
  }
  return {d,h,F,n};
}
function inflate(D,q,B){
  const n=D.length,k=B.length;
  return [...D.map((a,i)=>a|(i===q?0:((1<<k)-1)<<n)),...B.map(a=>(a<<n)|(1<<q))];
}
function trans(n){return Array.from({length:n},(_,j)=>((1<<n)-1)^((1<<(j+1))-1));}
function ordinal(X,Y){return [...X.map(a=>a|(((1<<Y.length)-1)<<X.length)),...Y.map(b=>b<<X.length)];}
function full(rows){const H=pathTable(rows);return H.h[H.F];}
function ports(D,q){
  const s=pathTable(D),t=pathTable(D,q);
  const G=s.h[s.F],J=s.h[s.F^(1<<q)],e=s.d[s.F*s.n+q],L=G-e;
  let W=-J;
  for(let m=1;m<=s.F;m++)W+=t.h[m]*s.h[s.F^m];
  return [G,J,e,L,W];
}



const A=[286,276,8,402,420,463,31,327,68],q=8;
function flip(T,u,v){const c=[...T];c[u]^=1<<v;c[v]^=1<<u;return c;}
function diff(i,a,b){
 const B=trans(a+b+1),T=inflate(A,q,B);
 return full(flip(T,i,9+a))-(903*2**B.length-7);
}
function vec(f,s){if(!s)return [f(0)];const a=f(1),b=(f(2)-a)/2,c=(f(3)-a-6*b)/6;return [a,b,c];}
const cert=Array.from({length:9},(_,i)=>{
 const tab=Array.from({length:4},(_,a)=>Array.from({length:4},(_,b)=>diff(i,a,b)));
 const t=(a,b)=>tab[a][b];
 const C=[vec(b=>t(0,b),false),vec(b=>t(0,b),true),vec(a=>t(a,0),true),
  [0,1,2].map(j=>vec(a=>vec(b=>t(a,b),true)[j],true))];
 return C;
});

const lower=inflate(A,q,trans(3));
const T=flip(lower,1,6),N=full(lower);
need(full(T)===N,"plateau");
let equal=[],bad=[];
for(let u=0;u<T.length;u++)for(let v=u+1;v<T.length;v++){
 const step=full(flip(T,u,v))-N;
 if(step===0)equal.push([u,v]);if(step>0&&step<=14)bad.push([u,v,step]);
}

need(bad.length===0 && JSON.stringify(equal)==="[[1,6],[9,10],[10,11]]","second plateau neighbors");
function F(s){return [1,2**s-2,3**s-3*2**s+3];}
function value(i,a,b){
 const [c,u,v,M]=cert[i];if(!a&&!b)return c[0];
 const dot=(x,y)=>x.reduce((acc,z,j)=>acc+z*y[j],0);
 if(!a)return dot(u,F(b));if(!b)return dot(v,F(a));
 return dot(F(b),M.map(row=>dot(row,F(a))));
}
for(let i=0;i<q;i++){
 const [c,u,v,M]=cert[i],parts=[c,u,v,M.flat()];
 need(parts.flat().every(x=>Number.isInteger(x)&&x>=0),"nonnegative coefficient certificate");
 need([c[0],u[0],v[0],M[0][0]].every(x=>x>=168),"constant bound");
}
need(JSON.stringify(cert[q])==="[[-984],[-586,199,0],[-1806,-903,0],[[-1408,-903,0],[199,0,0],[0,0,0]]]","root coefficient certificate");
for(let i=0;i<=q;i++)for(const [a,b] of [[4,1],[1,4],[0,5],[5,0]])need(diff(i,a,b)===value(i,a,b),"out-of-grid verification");

const F4=[[0,4],[0,5],[4,7],[5,7]];
const S=Array.from({length:16},(_,m)=>{
 let t=[...A];F4.forEach(([u,v],j)=>{if(m>>j&1)t=flip(t,u,v);});
 const [G,J,e,L,W]=ports(t,q);return [W,L+J-W];
});
const states=[[903,-7],[821,-11],[1074,11],[1188,17],[853,-13],[823,-21],[905,-3],[1094,-1],[1033,-19],[991,-23],[1025,15],[1132,17],[781,-13],[809,-19],[703,7],[903,7]];
need(JSON.stringify(S)===JSON.stringify(states),"four-edit cube");
for(let i=0;i<8;i++)for(let j=i+1;j<9;j++){
 const [G,J,e,L,W]=ports(flip(A,i,j),q),K=L+J-W,d=W-903;
 need(d!==0,"exterior slope");
 const lo=d>0?(-7-K)/d:(K-7)/(-d),hi=d>0?(7-K)/d:(K+7)/(-d);
 for(const r of [2,4,8,16,32]){
  const inBand=r>=lo&&r<=hi;
  need(!inBand||(i===1&&j===6&&r===8&&d*r+K===-7),"unexpected dyadic exterior move");
 }
 need(hi<=8,"large exterior move");
}
globalThis.console?.log(JSON.stringify(cert));
globalThis.console?.log("All-order contact coefficients, plateau neighbors and four-state excursion bound: PASS");
