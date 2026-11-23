"use strict";
const same=(a,b,msg)=>{if(JSON.stringify(a)!==JSON.stringify(b))throw Error(msg+": "+JSON.stringify(a));};
const prime=n=>{for(let d=2;d*d<=n;d++)if(n%d===0)return false;return n>1;};
function profile(w,n){
  let D=1,Q=2**(n+1);
  for(let i=0;i<n;i++)if(w>>i&1)
    for(let j=i+1;j<n;j++)if(!(w>>j&1)){
      const d=j-i-1;
      D+=2**d+(d===0?1:0);
      Q+=2**(n-1-d)*(3**d+1);
    }
  return [D,Q];
}
function words(qmax,fn){
  for(let n=2;2*3**(n-2)<=qmax;n++)
    for(let w=1;w<1<<(n-1);w+=2){
      const [D,Q]=profile(w,n);
      if(Q<=qmax)fn(D,Q);
    }
}
function ordClosure(qmax){
  let desc=1,g=new Map([["1,2",[1,2]]]);
  words(qmax,(D,Q)=>{desc++;g.set([D,Q].join(),[D,Q]);});
  let b=[[1,1]],s=new Set(["1,1"]);
  for(let k=0;k<b.length;k++){
    const [D,Q]=b[k];
    for(const [d,q] of g.values())if(Q*q<=qmax){
      const p=[D*d,Q*q], z=p.join();
      if(!s.has(z)){s.add(z);b.push(p);}
    }
  }
  return {desc,g:g.size,b};
}
const A=(D,Q)=>6*Q+5*D, B=(D,Q)=>15*Q+D;
const b=ordClosure(825);
same([b.desc,b.g,b.b.length],[28,16,50],"closure exact");
for(const p of [1747,4951]) {
  same(prime(p),true,"prime");
  same(b.b.some(([D,Q])=>A(D,Q)===p||B(D,Q)===p),false,"ordinal exclusion");
}
const X=1000000,s=[new Set(),new Set()];
function add(D,Q){
  for(let q=Q;A(D,q)<=X||B(D,q)<=X;q*=2){
    if(A(D,q)<=X)s[0].add(A(D,q));
    if(B(D,q)<=X)s[1].add(B(D,q));
  }
}
add(1,2);
words(Math.floor(X/6),(D,Q)=>add(D,Q));
same(s.map(x=>[x.size,[...x].filter(prime).length]),[[1159,193],[655,84]],"single-module census");
same(new Set([...s[0],...s[1]]).size,1813,"union");
function table(rows) {
  const n=rows.length,N=1<<n,d=new Float64Array(N*n),h=new Float64Array(N);
  h[0]=1;
  for(let i=0;i<n;i++){
    same(rows[i]>>i&1,0,"loop");
    for(let j=0;j<i;j++)same((rows[i]>>j&1)+(rows[j]>>i&1),1,"orientation");
    d[(1<<i)*n+i]=1;
  }
  for(let s=1;s<N;s++)for(let i=0;i<n;i++)if(d[s*n+i]){
    h[s]+=d[s*n+i];
    for(let z=rows[i]&~s;z;z&=z-1){
      const b=z&-z,j=31-Math.clz32(b);
      d[(s|b)*n+j]+=d[s*n+i];
    }
  }
  return h;
}
function literal(tags,b,M){
  const m=M.length,s=m+6, pos=[...Array(7).keys()].flatMap(i=>i===b?M.map((_,j)=>[i,j]):[[i,0]]);
  const rows=Array(s+2).fill(0);
  for(let i=0;i<s+2;i++)for(let j=i+1;j<s+2;j++){
    let fw;
    if(j<s)fw=pos[i][0]!==pos[j][0]||Boolean(M[pos[i][1]]>>pos[j][1]&1);
    else if(i<s)fw=!(tags[pos[i][0]]>>(j-s)&1);
    else fw=false;
    if(fw)rows[i]|=1<<j;else rows[j]|=1<<i;
  }
  return rows;
}
for(const [M,t,k,marks,want,dq] of [
  [[30,60,40,48,4,17],[2,0,3,0,3,3,3],1,[[0,5],[0,6]],[1747,5831,9915,12229],[23,272]],
  [[30,44,24,48,2,21],[1,0,1,3,2,3,3],2,[[0,3],[1,6]],[4951,6795,25313,34657],[31,328]]
]){
  const h=table(M),u=h.length-1;
  same([h[u],h.reduce((a,b,i)=>a+b*h[u^i],0)],dq,"module profile");
  for(let j=0;j<M.length;j++)same(h[u^(1<<j)]>1,true,"deletion cyclic");
  for(let c=0;c<4;c++){
    const tags=t.slice();
    marks.forEach(([a,j],i)=>{if(c>>i&1)tags[j]^=1<<a;});
    same(table(literal(tags,k,M)).at(-1),want[c],"full literal height");
  }
}

for(const [p,q] of [[23,296],[31,644]]){
  const hits=new Set();
  for(let n=2;2**(n-2)<p;n++)for(let w=1;w<1<<(n-1);w+=2){
    const r=profile(w,n);if(r[0]===p)hits.add(r[1]);
  }
  same([...hits],[q],"sharp odd profile ray");
}

const glued=new Set();
for(const [D,Q] of ordClosure(114).b)if(Q>1){glued.add(A(D,Q));glued.add(B(D,Q));}
const escape=new Set();
for(let n=1;n<=6;n++){
  const other=[];
  for(let i=0;i<n;i++)for(let j=i+2;j<n;j++)other.push([i,j]);
  for(let w=0;w<1<<other.length;w++){
    const rows=Array(n).fill(0);
    for(let i=0;i+1<n;i++)rows[i]|=1<<(i+1);
    other.forEach(([i,j],k)=>{if(w>>k&1)rows[i]|=1<<j;else rows[j]|=1<<i;});
    const h=table(rows),u=h.length-1,D=h[u],Q=h.reduce((a,b,i)=>a+b*h[u^i],0);
    for(const val of [A(D,Q),B(D,Q)])if(val<=689&&!glued.has(val))escape.add(val);
  }
}
same([...escape],[689],"first scalar escape");

const primeGlued=new Set();
for(const [D,Q] of ordClosure(291).b)if(Q>1){primeGlued.add(A(D,Q));primeGlued.add(B(D,Q));}
const primeEscapes=new Set();
for(let n=1;n<=7;n++){
  const other=[];
  for(let i=0;i<n;i++)for(let j=i+2;j<n;j++)other.push([i,j]);
  for(let w=0;w<1<<other.length;w++){
    const rows=Array(n).fill(0);
    for(let i=0;i+1<n;i++)rows[i]|=1<<(i+1);
    other.forEach(([i,j],k)=>{if(w>>k&1)rows[i]|=1<<j;else rows[j]|=1<<i;});
    const h=table(rows),u=h.length-1,D=h[u],Q=h.reduce((a,b,i)=>a+b*h[u^i],0);
    for(const val of [A(D,Q),B(D,Q)])if(val<=1747&&prime(val)&&!primeGlued.has(val))primeEscapes.add(val);
  }
}
same([...primeEscapes],[1747],"first prime scalar escape");
