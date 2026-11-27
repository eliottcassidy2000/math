// Exact census. Run: node replay.js
const assert=require("node:assert/strict");
const perms=A=>A.length ? A.flatMap((x,i)=>perms(A.filter((_,j)=>i!==j))
    .map(p=>[x,...p])) : [[]];
const pop=x=>x===0?0:1+pop(x&(x-1));
function Aof(n,k){
  let i=0;const A=Array.from({length:n},()=>Array(n).fill(0));
  for(let u=0;u<n;u++)for(let v=u+1;v<n;v++){
    A[u][v]=(k>>i++&1)?1:-1;A[v][u]=-A[u][v];}
  return A;
}
const arc=(A,u,v,m)=>A[u][v]*(((m>>u^m>>v)&1)?-1:1);
const allowed=(m,i)=>((m>>i)&1)!==0;
const swmask=[0,3,5,6];
const norm=(v)=>v>>3?v^15:v;
const ranks={1:[1,1,5,5],3:[3,3,3,3]};
const P4=perms([0,1,2,3]), seen=new Set(), counts=[0,0];
function h(A,m,P,profile){
  const slots=[];
  let k=0;
  for(const q of profile){
    for(let j=k+1;j<k+q;j++)slots.push(j); k+=q;
  }
  return P.filter(p=>slots.every(i=>arc(A,p[i-1],p[i],m)>0)).length;
}
function code(A,m) {
  let r=0,i=0;
  for(let u=0;u<4;u++)for(let v=u+1;v<4;v++)
    r|=(arc(A,u,v,m)>0?1:0)<<i++;
  return r;
}
const pft=(A,m)=>[arc(A,0,1,m)*arc(A,2,3,m),-arc(A,0,2,m)*arc(A,1,3,m),arc(A,0,3,m)*arc(A,1,2,m)];
const eq=(x,y)=>assert.deepEqual(x,y);
for(let k=0;k<64;k++){
  const A=Aof(4,k),seed=pft(A,0);
  const pf=Math.abs(seed.reduce((a,b)=>a+b,0)), all=swmask.map(m=>h(A,m,P4,[4]));
  eq([...all].sort(),ranks[pf]);
  eq(swmask.map(m=>pft(A,m)),Array(4).fill(seed));
  if(pf===1){
    const nt=swmask.filter(t=>t!==0&&swmask.every((m,i)=>all[i]===all[swmask.indexOf(norm(m^t))]));
    eq(nt.length,1);
  }
  if(!seen.has(k)){counts[pf===1?0:1]++;for(const m of swmask)seen.add(code(A,m));}
}
eq(counts,[12,4]);
console.log("balanced K4 orbits: 12 with (1,1,5,5); 4 with (3,3,3,3)");

const cases=[[3,[3]],[4,[4]],[4,[2,2]],[5,[5]],[5,[2,3]],[6,[6]],[6,[3,3]],[6,[2,2,2]],[7,[7]],[7,[1,3,3]]];
let checked=0;
for(const [n,profile] of cases){
  const V=Array.from({length:n},(_,i)=>i), P=perms(V);
  const degree=2*profile.reduce((s,c)=>s+Math.floor(c/3),0);
  const size=1<<(n-1);
  for(const k of [0, 713, 15491]){
    const A=Aof(n,k);
    const hs=Array.from({length:size},(_,m)=>h(A,m,P,profile));
    for(let S=0;S<(1<<n);S++){
      const w=pop(S);
      if((w&1)|| (w<=degree&&S!==0))continue;
      const s=hs.reduce((a,b,m)=>a+((pop(S&m)&1)?-b:b),0);
      if(S)eq(s,0);else eq(s*(1<<(n-profile.length)),P.length*size);
    }
    const cy=Array.from({length:size},(_,m)=>P.filter(p=>p.every((u,i)=>arc(A,u,p[(i+1)%n],m)>0)).length/n);
    const v=P.reduce((s,p)=>s+p.reduce((w,u,i)=>w*A[u][p[(i+1)%n]],1),0);
    eq(n*(1<<n)*cy.reduce((a,b)=>a+b,0),size*(P.length+v));
    if(n%2===0)eq(cy.reduce((s,c,m)=>s+((pop(m)&1)?-c:c),0),0);
    checked++;
  }
}
console.log("all Fourier support and cycle-mean checks passed on "+checked+" cases");
