// Exact small controls for the stochastic Wilson repair. Run with Node.
"use strict";
function assert(ok, msg) { if (!ok) throw Error(msg); }
const zero=n=>Array(n*n).fill(0);
const idx=(a,b,d)=>a*d+b;
const dot=(p,q)=>p.reduce((z,x,i)=>z+x*q[i],0);
const trans=(p,d)=>p.map((_,i)=>p[(i%d)*d+Math.floor(i/d)]);
function product(p,q,d) {
  const r=zero(d);
  for(let i=0;i<d;i++) for(let m=0;m<d;m++) if(p[idx(i,m,d)])
    for(let j=0;j<d;j++) r[idx(i,j,d)]+=p[idx(i,m,d)]*q[idx(m,j,d)];
  return r;
}
function pairs(k) {
  const ps=[];
  for(let i=0;i<k;i++) for(let j=i+1;j<k;j++) ps.push([i,j]);
  return ps;
}
function stats(chart,k,b) {
  // Inputs represent probabilities with common denominator b.
  // Return E, W, and D all multiplied by b^4.
  const d=k-1,h=b*b, ps=pairs(k), link=Array.from({length:k},()=>[]);
  ps.forEach(([i,j],l)=>{link[i][j]=chart[l];link[j][i]=trans(chart[l],d);});
  let E=0,W=0,V=0,D=0;
  for(const [i,j] of ps) {
    D+=b*b*(d*b*b-dot(link[i][j],link[i][j]));
    const rr=[link[i][j].map(p=>b*p)];
    for(let l=0;l<k;l++) if(l!==i&&l!==j) rr.push(product(link[i][l],link[l][j],d));
    const sum=zero(d);
    for(let l=0;l<rr.length;l++) {
      const p=rr[l];
      for(let x=0;x<p.length;x++) sum[x]+=p[x];
      V+=d*h*h-dot(p,p);
      for(let m=0;m<l;m++) W+=2*dot(p,rr[m]);
    }
    E+=sum.reduce((z,x)=>z+(x-h)**2,0);
  }
  assert(W===E+V,"variance identity");
  assert((k-1)*D<=V&&V<=(2*k-3)*D,"variance bracket");
  return {E,W,D};
}
function rooted(k,adj) {
  const heads=adj[0], leaves=heads.map(i=>adj[i].filter(x=>x!==0)), d=k-1;
  assert(heads.length===k&&new Set(leaves.flat()).size===k*d,"root chart");
  return pairs(k).map(([i,j])=>leaves[i].flatMap(x=>leaves[j].map(y=>+adj[x].includes(y))));
}
function petersen() {
  const X=pairs(5);
  return rooted(3,X.map(([x,y])=>X.flatMap(([a,b],i)=>
    [a,b].includes(x)||[a,b].includes(y)?[]:[i])));
}
function hoffmanSingleton() {
  const adj=Array.from({length:50},()=>[]);
  const add=(x,y)=>{adj[x].push(y);adj[y].push(x);};
  for(let i=0;i<5;i++) for(let x=0;x<5;x++) {
    add(5*i+x,5*i+(x+1)%5);
    add(25+5*i+x,25+5*i+(x+2)%5);
    for(let j=0;j<5;j++) add(5*i+x,25+5*j+(x+i*j)%5);
  }
  return rooted(7,adj);
}
function swapRows(P,d,a,b) {
  const vals=P.slice();
  for(let j=0;j<d;j++) {
    vals[idx(a,j,d)]=P[idx(b,j,d)];vals[idx(b,j,d)]=P[idx(a,j,d)];
  }
  return vals;
}
function check(k,C) {
  const d=k-1;
  assert(stats(C,k,1).W===0,"initial Moore witness");
  const unif=stats(C.map(()=>Array(d*d).fill(1)),k,d);
  assert(unif.E===0&&unif.W===C.length*d*(d-1)*d**4,"uniform false zero");
  let count=0;
  for(let l=0;l<C.length;l++) for(const [a,b] of pairs(d)) {
    const S=swapRows(C[l],d,a,b), b0=19;
    const cp=C.map((p,i)=>i===l?S:p);
    const cost=stats(cp,k,1).W;
    assert(cost===4*(2*k-3),"transposition cost");
    const mix=C.map((p,i)=>p.map((v,j)=>i===l?(b0-1)*v+S[j]:b0*v));
    assert(stats(mix,k,b0).W===cost*b0**3,"one-link affine slope");
    count++;
  }
  const n=Math.min(5,C.length);
  const switched=C.map((p,l)=>swapRows(p,d,l%d,(l+1)%d));
  const mix=C.map((p,l)=>p.map((x,j)=>l<n?x+switched[l][j]:2*x));
  let ans=0;
  for(let m=0;m<(1<<n);m++)
    ans+=stats(C.map((p,l)=>m&(1<<l)&&l<n?switched[l]:p),k,1).W;
  assert(ans*16===stats(mix,k,2).W*(1<<n),"independent rounding expectation");
  const N=17;
  const test=C.map((p,l)=>p.map((v,i)=>l===0?(N-1)*v+switched[0][i]:N*v));
  assert(stats(test,k,N).W<4*N**4,"strict rounding certificate");
  const best=[C[0],switched[0]].map(p=>[p,stats([p,...C.slice(1)],k,1).W])
    .sort((a,b)=>a[1]-b[1])[0];
  assert(best[1]===0&&best[0].every((v,i)=>!v||test[0][i]>0),"support rounding");
  console.log(JSON.stringify({valency:k,switches:count,blendVertices:1<<n,
    uniformAction:unif.W/d**4,oneLinkCertificate:(4*(2*k-3))+"/"+N,roundedAction:best[1]}));
}
check(3,petersen());
check(7,hoffmanSingleton());
