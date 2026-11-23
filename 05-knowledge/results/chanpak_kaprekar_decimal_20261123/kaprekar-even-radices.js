function H(r){let n=r.length,N=1<<n,d=new Uint32Array(n*N);for(let i=0;i<n;i++)d[(1<<i)*n+i]=1;for(let s=1;s<N;s++)for(let i=0;i<n;i++){let h=d[s*n+i];if(h)for(let w=r[i]&~s;w;w&=w-1){let b=w&-w;d[(s|b)*n+31-Math.clz32(b)]+=h;}}let h=0;for(let i=0;i<n;i++)h+=d[(N-1)*n+i];return h;}
function anchor(n,z){let r=Array.from({length:n},(_,i)=>i<n-1?1<<(i+1):0);for(let i=0,k=0;i<n;i++)for(let j=i+2;j<n;j++,k++)if(z>>k&1)r[i]|=1<<j;else r[j]|=1<<i;return r;}
function strong(r){let S=1;for(let j=0;j<r.length;j++)for(let i=0;i<r.length;i++)if(r[i]&S)S|=1<<i;return S===(1<<r.length)-1;}
function iso(r,s){const n=r.length;let pop=(x)=>x===0?0:1+pop(x&(x-1));const R=r.map(pop),S=s.map(pop);const p=[];function rec(i,vis){if(i===n)return true;for(let j=0;j<n;j++)if(!(vis>>j&1)&&R[i]===S[j]&&p.every((q,k)=>((r[i]>>k&1)===(s[j]>>q&1)))){p.push(j);if(rec(i+1,vis|1<<j))return true;p.pop();}return false;}return rec(0,0);}
function join(B){let R=[];for(let r of B){let n=R.length;R=R.map(x=>x|(((1<<r.length)-1)<<n)).concat(r.map(x=>x<<n));}return R;}
function reps(n,h){let V=[];for(let z=0;z<(1<<((n-1)*(n-2)/2));z++){let R=anchor(n,z);if(strong(R)&&H(R)===h&&!V.some(S=>iso(R,S)))V.push(R);}return V;}
function shapes(n,parts){let V=[];for(let B of parts){let spare=n-B.reduce((x,y)=>x+y.length,0);for(let a=0;a<=spare;a++)if(B.length===1)V.push(join([...Array(a).fill([0]),B[0],...Array(spare-a).fill([0])]));else for(let b=0;b<=spare-a;b++)V.push(join([...Array(a).fill([0]),B[0],...Array(b).fill([0]),B[1],...Array(spare-a-b).fill([0])]));}return V;}
function flip(r,i,j){r[i]^=1<<j;r[j]^=1<<i;}
function tests(V,t,two){let count=0,hits=0;for(let r of V){let E=[];for(let i=0;i<r.length;i++)for(let j=i+1;j<r.length;j++)E.push([i,j]);for(let i=0;i<E.length;i++){flip(r,...E[i]);if(two){for(let j=i+1;j<E.length;j++){flip(r,...E[j]);count++;hits+=H(r)===t;flip(r,...E[j]);}}else {count++;hits+=H(r)===t;}flip(r,...E[i]);}}return [V.length,count,hits];}


function assert(x){if(!x)throw Error('audit assertion failed');}
function exact(R,h){assert(R.every((r,i)=>!(r>>i&1)&&r>>R.length===0));for(let i=0;i<R.length;i++)for(let j=0;j<i;j++)assert((R[i]>>j&1)+(R[j]>>i&1)===1);assert(H(R)===h);}
function literal(A,edits,expected,p){const R=A.slice(),trace=[H(R)];for(const e of edits){flip(R,...e);trace.push(H(R));}assert(String(trace)===String(expected));const Q=p||A.map((_,i)=>i);for(let i=0;i<A.length;i++)for(let j=0;j<A.length;j++)assert((A[i]>>j&1)===(R[Q[i]]>>Q[j]&1));return trace;}
function endpoint(R){const n=R.length,N=1<<n,d=Array.from({length:N},()=>Array(n).fill(0));for(let i=0;i<n;i++)d[1<<i][i]=1;for(let s=1;s<N;s++)for(let i=0;i<n;i++)if(d[s][i])for(let w=R[i]&~s;w;w&=w-1){const b=w&-w;d[s|b][31-Math.clz32(b)]+=d[s][i];}return d;}
function run(log){
 const W=[131070,131068,262136,196592,262112,65472,196480,130816,196096,64512,129024,126976,122880,180224,98304,196608,9064,89763];
 const U=[16382,32764,16376,16368,49120,65472,65408,48896,65024,48128,14336,12288,24576,16384,36509,15375];
 const V=[4094,4092,8184,6128,4064,8128,6016,7936,3584,7168,4096,5192,275];
 for(const [a,h] of [[W,5360341],[U,765763],[V,21879]]){
  exact(a,h);log('common factor',a.length,h);
  const n=a.length-2;for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)assert(a[i]>>j&1);
 }
 const A=[214,76,72,49,102,199,136,286,127],I=460,n=A.length,F=(1<<n)-1,D=endpoint(A),DR=endpoint(A.map((r,i)=>r^F^(1<<i)));
 let G=0,L=0;
 for(let s=1;s<=F;s++){
  const t=F^s,h=D[s].reduce((acc,x,i)=>acc+x*(I>>i&1),0);
  G+=h*(t?DR[t].reduce((x,y)=>x+y,0):1);
  L+=h*DR[t].reduce((acc,x,i)=>acc+x*(I>>i&1),0);
 }
 assert(String([H(A),G,L])==='881,10478,5276');
 const f=((1<<9)-1)<<(n+1),RR=A.map((r,i)=>r|((I>>i&1)<<n));RR.push((F^I)|f);
 for(let i=n+1;i<n+10;i++)RR.push(F|(f&~((1<<(i+1))-1)));
 exact(RR,5360341);log('shell',H(A),G,L,RR.length);
 const s6=literal([54,60,56,49,32,0],[[0,4],[3,5],[0,4],[3,5]],[5,11,25,15,5]);log('base6 labelled',s6);
 const R=[118,124,120,113,96,64,0];
 log('base6 unlabelled',literal(R,[[3,6],[5,6],[3,6]],[5,25,15,5],[0,1,2,3,4,6,5]));
 const B5=reps(4,5);assert(B5.length===1&&String(tests(shapes(6,B5.map(r=>[r])),25,false))==='3,45,0');
 for(let i=0;i<64;i++)assert(H(anchor(5,i))<=15);
 for(let i=0;i<1024;i++)assert(H(anchor(6,i))!==21);
 const s8=literal(R,[[0,4],[4,6],[3,5],[4,6],[0,4],[3,5]],[5,11,35,55,25,15,5]);log('base8 proxy',s8);
 const at26=BigInt(25*3)/5n, r26=s8.map(t=>t*Number(at26));assert(String(r26)==='75,165,525,825,375,225,75');log('base26',r26);
 const start=[444,401,138,322,76,414,423,24,148],slot=start.slice(),read=[H(slot)];
 for(const e of [[1,2],[1,4],[0,2],[5,7],[5,8]]){flip(slot,...e);read.push(H(slot));}
 assert(String(slot)==='440,389,137,322,78,30,423,56,180'&&String(read)==='1105,765,955,1517,2265,2431');log('selector',read);
 for(const y of [1n,3n,5n,7n,9n,15n,27n]){
  const m=109395n*y,h=m-2n,F=y*h,d8=m/9n,K8=d8*h/5n;
  assert(K8===2431n*F);
  const f=g=>2n*g>m?2n*g-m:m-2n*g;
  assert(f(d8)===7n*d8&&f(7n*d8)===5n*d8&&f(5n*d8)===d8&&f(3n*d8)===3n*d8);
  log('base',String(m-1n),'denominator-nine trace',s8.map(a=>String(BigInt(a)*K8)));
  if(y===7n){const K6=m*h/35n;assert(K6===21879n*765763n);assert(f(m/7n)===5n*m/7n&&f(5n*m/7n)===3n*m/7n&&f(3n*m/7n)===m/7n);log('base',String(m-1n),'denominator-seven trace',s6.map(a=>String(BigInt(a)*K6)));}
 }
}

run(console.log);

function spineRows(n,a,b,ab){
 const full=(1<<n)-1,r=[];
 for(let i=0;i<n;i++)r.push((full^((1<<(i+1))-1))|((1-(a>>i&1))<<n)|((1-(b>>i&1))<<(n+1)));
 r.push(a|(ab?1<<(n+1):0),b|(ab?0:1<<n));
 return r;
}
function spineCount(r) {
 const n=r.length-2, edge=(u,v)=>Boolean(r[u]>>v&1), ix=(a,b,c)=>(a*3+b)*2+c;
 let ans=0n;
 for(const [x,y] of [[n,n+1],[n+1,n]]){
  let q=Array(18).fill(0n);q[0]=1n;
  for(let v=0;v<n;v++){
   const z=Array(18).fill(0n);
   for(let a=0;a<3;a++)for(let b=0;b<3;b++)for(let c=0;c<2;c++){
    const w=q[ix(a,b,c)];if(!w)continue;
    z[ix(edge(v,x)?1:2,b,c)]+=w;
    if(b||edge(x,v))z[ix(a,edge(v,y)?1:2,c)]+=w;
    if(c||edge(y,v))z[ix(a,b,1)]+=w;
   }
   q=z;
  }
  for(let a=0;a<2;a++)for(let b=0;b<2;b++)if(b||edge(x,y))
    for(let c=0;c<2;c++)ans+=q[ix(a,b,c)];
 }
 return ans;
}
for(const [n,a,b,ab,h] of [
 [16,9064,24227,false,5360341],[13,829,533,false,765763],[11,1096,275,true,21879]
]){
 const r=spineRows(n,a,b,ab);exact(r,h);assert(spineCount(r)===BigInt(h));
}
for(const [h,r] of [
 [7735,[670,200,322,468,166,847,17,356,595,222]],
 [5355,[584,509,913,580,649,733,276,841,569,66]],
 [3927,[160,245,625,199,681,328,17,100,735,235]],
 [3003,[602,564,977,246,128,661,562,67,763,152]],
 [17017,[1494,2028,48,997,682,1921,564,580,724,5,988]],
 [21879,[1024,1017,539,1377,265,789,821,381,1541,1177,246]]
])exact(r,h);
const seed=[118,124,120,113,96,64,0];
const rel=literal(seed,[[0,4],[4,6],[0,2],[4,5],[0,4]],[5,11,35,29,25,5],[1,2,3,0,5,6,4]);
assert(String(rel.map(a=>15*a))==="75,165,525,435,375,75");
for(let i=0;i<1024;i++)assert(H(anchor(6,i))!==35);
const b5=reps(4,5).map(r=>[r]);
assert(String(tests(shapes(6,b5),35,false))==="3,45,0");
const bank=shapes(11,b5),es=[];
for(let i=0;i<11;i++)for(let j=i+1;j<11;j++)es.push([i,j]);
let pairs=0,peaks=0,next=0,hits=0;
for(const r of bank)for(let a=0;a<es.length;a++){
 flip(r,...es[a]);
 for(let b=a+1;b<es.length;b++){
  flip(r,...es[b]);pairs++;
  if(H(r)===35){
   peaks++;
   for(const e of es){flip(r,...e);next++;hits+=H(r)===25;flip(r,...e);}
  }
  flip(r,...es[b]);
 }
 flip(r,...es[a]);
}
assert(String([bank.length,pairs,peaks,next,hits])==="8,11880,50,2750,0");
console.log("normalized minimum census",bank.length,pairs,peaks,next,hits);
