"use strict";
const assert=require("node:assert/strict");
const {boards,catalogue}=require("./minimum_trades");

function F(v){let k=1;while(v){let u=v>>1;if(!(v&1))k+=u+(u&1);v=u;}return k;}
const ff=(b,c)=>b*c+F(b-1)*F(c-1)-(b&1)*(c&1);
function exactChange(b,n,i,j) {
 const U=(1<<n)-1, kw=b-1, mask=(1<<i)|(1<<j), inner=((1<<j)-1)^((1<<(i+1))-1);
 const low=x=>x&-x, high=x=>1<<(31-Math.clz32(x));
 const a=x=>+(x===0||!(low(x)&kw)), d=x=>+(x===0||!!(high(x)&kw));
 let df=0,de=0,dr=0;
 for(let S=U;S;S=(S-1)&U)if((S&mask)===mask){
  const T=U^S,I=S&inner,plo=S&((1<<i)-1),phi=S&~((1<<(j+1))-1);
  let dk=I?0:-1, da=I?0:-a(S), dd=I?0:-d(S);
  for(let A=I;;A=(A-1)&I){
   const L=plo|A|(1<<j),R=phi|(I^A)|(1<<i);
   dk++;da+=a(L);dd+=d(R);if(!A)break;
  }
  if(S===U)de=da;
  df+=dd*a(T)+d(T)*da;
  dr+=dd+dk*d(T);
 }
 return [df,de,dr];
}
function literal(b,c,i,j,spoke=-1){
 const l=Math.max(1,32-Math.clz32(b-1),j+1,spoke+1),r=Math.max(1,32-Math.clz32(c-1));
 const x=l,y=l+r+1, n=y+1, A=new Uint32Array(n);
 const add=(i,j)=>A[i]|=1<<j;
 for(let u=0;u<l;u++){
   for(let v=u+1;v<l;v++)if(u===i&&v===j)add(v,u);else add(u,v);
   if((((b-1)>>u)&1)^(u===spoke))add(u,x);else add(x,u);
 }
 for(let i=x+1;i<y;i++){
   for(let j=i+1;j<y;j++)add(j,i);
   if(((c-1)>>(i-x-1))&1)add(y,i);else add(i,y);
 }
 for(let i=0;i<=x;i++)for(let j=x+1;j<=y;j++)if(i===x&&j===y)add(i,j);else add(j,i);
 const M=(1<<n)-1,D=new Float64Array((M+1)*n);
 for(let i=0;i<n;i++)D[(1<<i)*n+i]=1;
 for(let s=1;s<M;s++)for(let i=0;i<n;i++){
  const k=D[s*n+i];if(!k)continue;
  for(let m=A[i]&~s;m;m&=m-1){const j=m&-m;D[(s|j)*n+31-Math.clz32(j)]+=k;}
 }
 let k=0;for(let i=0;i<n;i++)k+=D[M*n+i];return k;
}
function forks(){
 const ws=[[7,1,0,2,0,11,13,9],[7,1,0,2,1,11,13,9],[9,1,0,3,0,17,19,15],
 [9,1,1,3,1,17,15,19],[5,3,0,2,0,29,31,27],[25,1,2,4,1,47,49,45],
 [73,1,4,6,2,141,143,139],[289,1,6,8,3,573,575,571]];
 for(const [b,c,i,j,z,a,q,w] of ws){
  const [df,de,dr]=exactChange(b,j+1,i,j);
  const got=[ff(b,c),ff(b,c)+c*dr+F(c-1)*df-(c&1)*de,literal(b,c,-1,-1,z)];
  if(got.join()!==[a,q,w].join()||literal(b,c,i,j)!==q)throw Error("fork");
 }
 let n=0;
 for(let b=1;b<=64;b++)for(let i=0;i<5;i++)for(let j=i+2;j<6;j++) {
  const [df,de,dr]=exactChange(b,Math.max(j+1,32-Math.clz32(b-1)),i,j);
  for(const c of [1,2,3,5]) {
   if(literal(b,c,i,j)!==ff(b,c)+c*dr+F(c-1)*df-(c&1)*de)throw Error("mixed interval");
   n++;
  }
 }
 return {profile:exactChange(9,4,1,3),witnesses:ws.length,mixedChecks:n};
}

function gap(s,p,J,B,L){
 function Q(m,a){const k=Math.floor(m/a);return a*k*k+(m%a)*(2*k+1);}
 const m=s*p,n=57-s,b=Math.max(0,B-J);
 if(J>L)throw Error("degree capacity");
 return 56*(J*Q(m,p)+b*Q(m,14)+(L-J-b)*Q(m,28)+(n-L)*Q(m,56))
  -n*m*m-(8-s)*(s+7)*m*(56-p);
}
function tight(qs){
 const M=qs.map(q=>q.reduce((m,v)=>m|(1n<<BigInt(v)),0n));
 const vs=Array.from({length:57},(_,i)=>i).filter(i=>![0,28,56].includes(i));
 const bit=x=>1n<<BigInt(x);
 const pop=x=>{let n=0;for(;x;x&=x-1n)n++;return n;};
 const sets6=new Set(),sets7=new Set();
 function pad(m, d, j){if(d===0){(pop(m)===6?sets6:sets7).add(m);return;}
 for(let i=j;i<vs.length;i++)if(!(m&bit(vs[i])))pad(m|bit(vs[i]),d-1,i+1);}
 for(let i=0;i<M.length;i++)for(let j=i+1;j<M.length;j++){
  const m=M[i]|M[j],k=pop(m);
  if(k<=6)pad(m,6-k,0);if(k<=7)pad(m,7-k,0);
 }
 const filter=sets=>[...sets].filter(s=>{
  const good=M.filter(m=>(m&s)===m);if(good.length<2)return false;
  let common=s;for(const g of good)common&=g;return common===0n;
 });
 const [f6,f7]=[filter(sets6),filter(sets7)], key=new Set(M);
 let best=0;
 for(const s of f7){
  const a=vs.filter(v=>s&bit(v));
  function perm(j){
   if(j===7){
    for(const len of [7,3]){
     let ok=0;
     for(let i=0;i<7;i++){
      const lo=i<len?0:len, d=i<len?len:7-len;
      const ne=s^bit(a[i])^bit(a[lo+(i-lo+d-1)%d])^bit(a[lo+(i-lo+1)%d]);
      if(key.has(ne))ok++;
     }
     best=Math.max(best,ok);
    }
   }else for(let i=j;i<7;i++){[a[i],a[j]]=[a[j],a[i]];perm(j+1);[a[i],a[j]]=[a[j],a[i]];}
  }perm(0);
 }
 const g=[[37,7,[8]],[38,6,[8,10]],[39,5,[6,8,10,12,14]],[40,4,[4,6,8,10,12,14]]]
  .map(([h,s,ps])=>[h,s,ps.map(p=>gap(s,p,Math.ceil(5*s/4),Math.max(0,54-h-Math.floor((55-h)/2)),54-h))]);
 return {sizes:[sets6.size,sets7.size,f6.length,f7.length,best],g};
}
function slack(){
 const o=[];
 for(let h=37;h<=40;h++)for(let s=44-h;s<=7;s++){
  const z=4*(h+s-44),E=Math.floor((55-h)/2);
  for(let l=0;l<=s;l++){
   if(!z&&l)continue;
   let min=l*Math.max(0,5-l);if(l===5)min=4;if(l===6||l===7)min=2;
   if(min>z)continue;
   const L=54-h-l,B=Math.max(0,L-E),J=Math.max(0,Math.ceil((5*(s-l)-z+Math.max(0,min-l*(s-l)))/4));
   const ps=[];for(let p=2;p<=(B?14:28);p+=2)if(p>=s-1&&gap(s,p,J,B,L)<=0)ps.push(p);
   if(ps.length)o.push({h,s,l,J,ps});
  }
 }
 return o;
}
function bank(n){
 const E=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)E.push([i,j]);
 const m=1<<E.length, full=(1<<n)-1, root=n-1, r=1<<root;
 const H=new Int32Array(m), Z=new Int32Array(m), V=new Int32Array(m), D=new Int32Array((full+1)*n),C=new Int32Array(D.length),A=new Int32Array(n),P=new Int32Array(full+1);
 for(let w=0;w<m;w++){
  A.fill(0);D.fill(0);C.fill(0);
  for(let i=0;i<E.length;i++){let [a,b]=E[i];if(w&(1<<i))A[b]|=1<<a;else A[a]|=1<<b;}
  for(let i=0;i<n;i++)D[(1<<i)*n+i]=1;
  C[r*n+root]=1;P[0]=1;
  for(let s=1;s<=full;s++){
   let u=0;
   for(let i=0;i<n;i++){
    const a=D[s*n+i],b=C[s*n+i];u+=a;if(!a&&!b)continue;
    for(let j=A[i]&~s;j;j&=j-1){let k=j&-j,idx=(s|k)*n+31-Math.clz32(k);D[idx]+=a;C[idx]+=b;}
   }
   P[s]=u;
  }
  let f=0;
  for(let s=r;s<=full;s++){
    let k=0;for(let i=0;i<n;i++)k+=C[s*n+i];
    f+=k*P[full^s];
  }
  H[w]=P[full];Z[w]=D[full*n+root];V[w]=f;
 }
 let cnt=0,from=0,byType=[0,0];const uniq=new Map(),ePairs={};
 for(let w=0;w<m;w++)if(Z[w]===1){
  from++;
  for(let i=0;i<E.length;i++){
   const u=w^(1<<i),d=9*(H[u]+V[u]-H[w]-V[w])-Z[u]+1;
   ePairs[Z[u]]=(ePairs[Z[u]]||0)+1;
   if(d===2){cnt++;byType[+(E[i][1]===root)]++;const k=[9*(H[w]+V[w])-1,...[H[w],Z[w],V[w],H[u],Z[u],V[u]]];uniq.set(k.join(','),[k,w,E[i]]);}
  }
 }
 return {n,from,cnt,byType,ePairs,types:[...uniq.values()]};
}

const seed=require("node:fs").readFileSync(process.argv[2],"utf8");
const qs=[...catalogue(boards(seed).full).table.keys()].map(q=>q.split(",").map(Number));
assert.equal(qs.length,152);
const t=tight(qs);
assert.deepEqual(t.sizes,[2373,67272,0,43,5]);
assert.deepEqual(t.g,[[37,7,[156800]],[38,6,[85120,98336]],[39,5,[35392,37632,30912,35392,61152]],[40,4,[11200,10304,8064,4480,-3136,-15680]]]);
const rows=slack();
assert.deepEqual(rows.map(x=>[x.h,x.s,x.l,x.ps]),[[38,7,7,[8]],[39,6,5,[8,10]],[39,6,6,[8,10]],[39,7,6,[8]],[39,7,7,[8]],[40,4,0,[12,14]],[40,5,4,[6,8,10,12,14]],[40,5,5,[6,8,10,12,14]],[40,6,5,[8,10]],[40,6,6,[6,8,10,12]],[40,7,6,[8]],[40,7,7,[8]]]);
const f=forks();
assert.deepEqual(f,{profile:[0,2,0],witnesses:8,mixedChecks:2560});
const r=[];
for(let n=1;n<=(process.argv.includes("--order7")?7:6);n++){
 const v=bank(n);assert.equal(v.cnt,0);
 assert.equal(v.from,[0,1,1,4,30,312,4080,64080][n]);
 r.push([v.n,v.from,v.cnt]);
}
console.log(JSON.stringify({t,rows,f,r},null,2));
