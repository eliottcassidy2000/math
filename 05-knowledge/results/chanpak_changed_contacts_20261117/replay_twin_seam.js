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


const minus=[286,276,8,402,420,463,31,327,68],plus=[466,388,393,403,294,463,158,16,192],q=8;
const start=[270,268,264,256,239,463,271,335,16], finish=[270,140,264,272,167,463,31,333,82];
const E=[[1,7],[1,8],[3,4],[6,8],[4,6]];
function change(base,flips,mask){
  const a=[...base];
  flips.forEach(([u,v],i)=>{if(mask>>i&1){a[u]^=1<<v;a[v]^=1<<u;}});
  return a;
}
function coe(a){const [G,J,e,L,W]=ports(a,8);return [W,L+J-W];}
function perm(a,fx,i=0){if(i===a.length)return fx(a);for(let j=i;j<a.length;j++){[a[i],a[j]]=[a[j],a[i]];perm(a,fx,i+1);[a[i],a[j]]=[a[j],a[i]];}}
function relabel(a,p){return p.map(i=>p.reduce((r,j,k)=>r|(((a[i]>>j)&1)<<k),0));}
function flips(a,b){const f=[];for(let i=0;i<a.length;i++)for(let j=i+1;j<a.length;j++)if((a[i]^b[i])>>j&1)f.push([i,j]);return f;}
function prime(p){for(let i=2;i*i<=p;i++)if(p%i===0)return false;return p>1;}
function check(x,y,msg){need(JSON.stringify(x)===JSON.stringify(y),msg);}
check(ports(minus,q),[815,199,118,697,903],"negative profile");
check(ports(plus,q),[817,221,128,689,903],"positive profile");
for(const [B,h,P,x,y] of [[trans(3),1,8,1031,1033],
                           [trans(5),1,32,4127,4129],
                           [[24,1,3,6,14],11,100,12889,12911]]){
  const s=pathTable(B);
  const p=s.h.reduce((r,v,i)=>r+v*s.h[s.F^i],0);
  need(h===s.h[s.F]&&P===p&&x===129*P-h&&y===129*P+h,"module profile");
  need(prime(x)&&prime(y),"prime quotients");
  need(full(inflate(minus,q,B))===7*x&&full(inflate(plus,q,B))===7*y,"literal counts");
}
let smallest=Infinity,groups=[];
perm([0,1,2,3,4,5,6,7],p=>{
 const m=[...p,q],b=relabel(plus,m),f=flips(minus,b);
 if(f.length<smallest){smallest=f.length;groups=[];}
 if(f.length===smallest)groups.push([m,f]);
});
check(groups,[[[3,1,7,4,2,0,6,5,8],[[0,4],[0,5],[4,7],[5,7]]]],"pointed relabel minimum");
const [p,F]=groups[0],b=relabel(plus,p);
const frame=[8,2,3,7,0,1,4,5,6];
need(frame.slice(1).every((v,i)=>(minus[frame[i]]&b[frame[i]])>>v&1),"shared frame");
check(change(minus,F,15),b,"four exterior reversals");

const neighbors=[];
for(let i=0;i<9;i++)for(let j=i+1;j<9;j++){
  const [W,K]=coe(change(minus,[[i,j]],1)),d=W-903;
  need(d!==0,"neutral slope");
  const lo=d>0?[-7-K,d]:[K-7,-d],hi=d>0?[7-K,d]:[K+7,-d];
  if(hi[0]>=2*hi[1] && Math.max(2,lo[0]/lo[1])<=hi[0]/hi[1])
    neighbors.push([[i,j],lo,hi]);
}
check(neighbors,[[[1,6],[74,11],[88,11]],[[3,6],[96,46],[110,46]]],"allowed bracket entries");

const S=Array.from({length:32},(_,i)=>coe(change(start,E,i)));
check(change(start,E,31),finish,"five-arc endpoints");
let count=0,best=-Infinity,bestIt=[];
perm([0,1,2,3,4],xs=>{
 let m=0,min=Infinity;
 for(let i of xs){
  const n=m|(1<<i),d=S[n][0]-S[m][0],c=S[n][1]-S[m][1];
  if(d<0||2*d+c<=0)return;
  min=Math.min(min,2*d+c);m=n;
 }
 count++;
 if(min>best){best=min;bestIt=[];}
 if(min===best)bestIt.push(xs.join(""));
});
need(count===29&&best===60&&bestIt.includes("03214"),"optimal order");
const tri=[[1,6,7],[3,4,7],[4,5,6]], edge=new Set();
for(const [i,j,k] of tri){
 need((finish[i]>>j&1)===(finish[j]>>k&1)&&(finish[i]>>j&1)!==(finish[i]>>k&1),"not triangle");
 for(const [u,v] of [[i,j],[i,k],[j,k]])edge.add(9*u+v);
}
need(edge.size===9,"triangles share an edge");
globalThis.console?.log("Two twin-prime pairs, nontransitive prime pair, marked four-edit path, and five-arc optimality: PASS");
