"use strict";
const assert=(ok,why)=>{if(!ok)throw Error(why)};
const eq=(a,b,why)=>assert(JSON.stringify(a)===JSON.stringify(b),why+" "+JSON.stringify(a));
function table(rows) {
  const n=rows.length,N=1<<n,d=new Float64Array(N*n);
  for(let i=0;i<n;i++){
    assert(!(rows[i]>>i&1),"loop");
    for(let j=0;j<i;j++)assert(((rows[i]>>j&1)+(rows[j]>>i&1))===1,"orientation");
    d[(1<<i)*n+i]=1;
  }
  for(let s=1;s<N;s++)for(let i=0;i<n;i++)if(d[s*n+i]){
    for(let z=rows[i]&~s;z;z&=z-1){
      const b=z&-z,j=31-Math.clz32(b);
      d[(s|b)*n+j]+=d[s*n+i];
    }
  }
  const h=new Float64Array(N); h[0]=1;
  for(let s=1;s<N;s++)for(let i=0;i<n;i++)h[s]+=d[s*n+i];
  return h;
}
function strong(r){
  let z=1,old=0;
  while(old!==z){old=z;for(let i=0;i<r.length;i++)if(z>>i&1)z|=r[i];}
  if(z!==(1<<r.length)-1)return false;
  z=1;old=0;
  while(old!==z){old=z;for(let i=0;i<r.length;i++)if(r[i]&z)z|=1<<i;}
  return z===(1<<r.length)-1;
}
function toggle(r,e,c){
  const s=r.slice();
  e.forEach(([i,j],k)=>{if(c>>k&1){s[i]^=1<<j;s[j]^=1<<i;}});
  return s;
}
function fixture(r,e,want,allStrong=false){
  for(let c=0;c<4;c++){
    const f=toggle(r,e,c);
    eq(table(f).at(-1),want[c],"marked corner");
    if(allStrong)assert(strong(f),"need strong corner");
  }
}
fixture([202,76,9,64,239,207,132,14],[[3,4],[0,6]],[31,209,37,217]);
fixture([254,508,504,240,224,448,128,256,89],[[1,7],[3,6]],[103,411,221,721],true);
fixture([482,68,505,499,483,450,128,258,66],[[3,6],[4,7]],[149,577,343,1043],true);
function usingBoth(rows,e){
  const n=rows.length;let count=0;
  function go(s,tail,used){
    if(s===(1<<n)-1){if(used===3)count++;return;}
    for(let z=rows[tail]&~s;z;z&=z-1){
      const b=z&-z,j=31-Math.clz32(b);
      let u=used;
      e.forEach(([a,v],k)=>{if((tail===a&&j===v)||(tail===v&&j===a))u|=1<<k;});
      go(s|b,j,u);
    }
  }
  for(let i=0;i<n;i++)go(1<<i,i,0);
  return count;
}
const U=[202,76,9,64,239,207,132,14], ue=[[3,4],[0,6]];
eq([0,1,2,3].map(c=>usingBoth(toggle(U,ue,c),ue)),[0,42,0,44],"unit mixed sign");

const F={
 A:{tags:[2,0,3,0,3,3,3],block:1,marks:[[0,5],[0,6]],
    ab:[[3,1],[9,3],[15,5],[21,7]],ba:[[14,5],[48,17],[82,29],[98,35]]},
 B:{tags:[1,0,1,3,2,3,3],block:2,marks:[[0,3],[1,6]],
    ab:[[13,6],[18,9],[73,24],[98,39]],ba:[[18,9],[21,12],[100,51],[119,66]]}
};
function skeleton(tags,b,order){
  const fixed=[0,1,2,3,4,5,6].filter(i=>i!==b),out=[0,0,0];
  function has(i,j){
    if(i<7&&j<7)return i<j;
    if(i>=7&&j>=7)return i===8&&j===7;
    return i>=7 ? Boolean(tags[j]>>>(i-7)&1) : !(tags[i]>>>(j-7)&1);
  }
  for(let k=0;k<729;k++){
    let p=k;
    const base=[[],[],[]];
    for(const v of fixed){base[p%3].push(v);p=Math.floor(p/3);}
    for(let m=1;m<8;m++){
      const runs=base.map((a,i)=>m>>i&1?[...a,b].sort((x,y)=>x-y):a);
      const q=[...runs[0],order[0],...runs[1],order[1],...runs[2]];
      if(q.slice(1).every((v,i)=>has(q[i],v)))out[(m&1)+(m>>1&1)+(m>>2&1)-1]++;
    }
  }
  return out;
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
const TT=n=>Array.from({length:n},(_,i)=>(1<<n)-(1<<(i+1)));
const modules=[TT(1),TT(2),TT(4),[2,4,1],[30,4,8,18,6],[18,4,9,19,6]];
for(const [name,o] of Object.entries(F))for(let c=0;c<4;c++){
  const tags=o.tags.slice();
  o.marks.forEach(([a,j],i)=>{if(c>>i&1)tags[j]^=1<<a;});
  const ab=skeleton(tags,o.block,[7,8]), ba=skeleton(tags,o.block,[8,7]);
  eq(ab,[...o.ab[c],0],"skeleton ab");
  eq(ba,[...o.ba[c],0],"skeleton ba");
  for(const M of modules){
    const h=table(M),D=h.at(-1),Q=h.reduce((z,a,i)=>z+a*h[h.length-1-i],0);
    const R=literal(tags,o.block,M), H=table(R).at(-1);
    eq(H,(ab[0]+ba[0])*D+(ab[1]+ba[1])*(Q-2*D),"module substitution");
    if(c===3)assert(strong(R),"strong upper");
    if(name==="A"&&M===modules[4])eq(H,[409,1365,2321,2863][c],"p409");
    if(name==="B"&&M===modules[5])eq(H,[1511,2067,7753,10577][c],"p1511");
  }
}
const prime=n=>{for(let i=2;i*i<=n;i++)if(n%i===0)return false;return n>=2};
for(let p=5;p<200;p+=2)if(prime(p))for(let r=1;r<p;r++){
  let ok=false;
  for(let D=1;D<p&&!ok;D++)for(let Q=1;Q<p&&!ok;Q++)
    if((6*Q+5*D-r)%p===0&&(Q+D)%p&&(20*Q+17*D)%p&&(34*Q+29*D)%p)ok=true;
  assert(ok,"phase wheel at "+p);
}
for(const r of [1,2,4,5,7,8]){
  let count=0;
  for(let Q=0;Q<9;Q++)for(let D=0;D<9;D++)
    if((Q-2*r)%3===0&&(6*Q+5*D-r)%9===0&&(34*Q+29*D)%9!==0) {
      assert((20*Q+17*D)%3!==0&&(34*Q+29*D)%3===0,"3-adic route");
      count++;
    }
  eq(count,2,"3-adic lifts");
}
