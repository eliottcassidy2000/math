// node replay_h5_shifts.js
"use strict";
const rows=[8,1,3,6], h=5, f=32;
function chk(x,msg){if(!x)throw Error(msg);}
function profiles(){
 const n=4,full=15;
 function ends(R){
  const E=Array.from({length:16},()=>Array(4).fill(0));
  for(let i=0;i<4;i++)E[1<<i][i]=1;
  for(let s=1;s<16;s++)for(let j=0;j<4;j++)if(s>>j&1)
   for(let i=0;i<4;i++)if(R[i]>>j&1)E[s][j]+=E[s^(1<<j)][i];
  return E;
 }
 const E=ends(rows),A=ends(rows.map((r,i)=>15^r^(1<<i)));
 const H=E.map((r,i)=>i?r.reduce((p,x)=>p+x,0):1);
 chk(H[15]===h&&H.reduce((p,x,i)=>p+x*H[15^i],0)===f,"core");
 return H.map((_,U)=>{
  let p=0,q=0,S=0,T=0,L=0;
  for(let i=0;i<n;i++)if(U>>i&1){S+=E[15][i];
   for(let a=1;a<16;a++){
    p+=A[a][i]*H[15^a];
    for(let j=0;j<n;j++)if(U>>j&1)L+=E[15^a][j]*A[a][i];
   }
  }else{T+=A[15][i];
   for(let a=1;a<16;a++)q+=E[a][i]*H[15^a];
  }
  return {p,q,S,T,K:h+p-S-L};
 });
}
function graph(){
 const G=new Map(),B=new Map();
 for(let c=-h;c<=f;c++)for(let e=-h;e<=f;e++){
  const s=c+","+e, adj=[];
  for(let x=0;x<=1;x++)for(let z=0;z<=1;z++){
   const a=c+f*x-h*z,b=2*e-f*z+h*x;
   if(a%2||b< -h||b>f)continue;
   const j=a/2+","+b;adj.push([j,x,z]);
   B.set(j,[...(B.get(j)||[]),s]);
  }
  G.set(s,adj);
 }
 return {G,B};
}
function canReach(B,end){
 const Q=[end],set=new Set(Q);
 for(let j=0;j<Q.length;j++)for(const u of B.get(Q[j])||[]){
  if(!set.has(u)){Q.push(u);set.add(u);}
 }
 return set;
}
const abs=x=>x<0n?-x:x;
function gcd(x,y){x=abs(x);y=abs(y);while(y){[x,y]=[y,x%y];}return x;}
function independent(basis,v){
 let a=v.slice();
 for(const b of basis){
  const j=b.findIndex(x=>x!==0n);
  if(a[j]===0n)continue;
  const z=a[j];
  a=a.map((x,i)=>b[j]*x-z*b[i]);
  const d=a.reduce(gcd,0n);if(!d)return false;
  a=a.map(x=>x/d);
 }
 const d=a.reduce(gcd,0n);
 if(!d)return false;
 // Echelon order is intentional; every older pivot is already killed.
 basis.push(a.map(x=>x/d));
 return true;
}
function read(V,x,z,px,pz){
 const [o,L,X,Z,R,Y,H,J]=V;
 x=BigInt(x);z=BigInt(z);
 return [o,2n*L,X+x*L,Z+z*L,2n*R+x*o,2n*Y+z*o,
  H+x*(L-R-BigInt(px)*o),J+z*(L-Y-BigInt(pz)*o)];
}
const {G,B}=graph(), P=profiles();
chk(new Set(P.map(a=>a.p+','+a.q)).size===16,'unique intrinsic masks');
let pairs=0, maxGap=0n, states=0,witness=null;
for(let u=0;u<16;u++)for(let v=0;v<16;v++){
 const a=P[u], b=P[v], start=b.p+","+(f-h-b.q), end=a.p+","+(f-h-a.q);
 const live=canReach(B,end);if(!live.has(start))continue;
 const bases=new Map(),q=[],values=new Map();
 function add(s,x0,z0,x1,z1,V){
  const key=[s,x0,z0,x1,z1].join(","), basis=bases.get(key)||[];
  if(independent(basis,V)){
   bases.set(key,basis);q.push([key,s,x0,z0,x1,z1,V]);
  }
 }
 add(start,-1,-1,1,1,[1n,1n,0n,0n,0n,0n,1n,1n]);
 for(let i=0;i<q.length;i++){
  const [s,key,x0,z0,x1,z1,V]=q[i];
  if(key===end&&x0!==-1){
   const [o,L,X,Z,R,Y,H,J]=V;
   const d=BigInt(b.K-a.K+x0*b.T-a.S*(1-z1))*o
        +BigInt(h)*(H-J)+BigInt(b.q)*X-BigInt(a.p)*(L-o-Y);
   let dSet=[];
   for(let e=0;e<2;e++)for(let g=0;g<2;g++){
    const t=d+o*BigInt(e*(h*(1-x1)-a.T)+g*(b.S-h*z0));
    chk(abs(t)<=2n,"gap outside ±2");
    // A fixed terminal state must have a single exact affine readout.
    const k=s+","+e+","+g;
    if(values.has(k))chk(values.get(k)===t,"nonconstant terminal readout");
    else values.set(k,t);
    dSet.push(t);
    if(abs(t)>maxGap){maxGap=abs(t);witness=[u,v,String(t)];}
   }
  }
  for(const [j,x,z] of G.get(key))if(live.has(j)){
   add(j,x0<0?x:x0,z0<0?z:z0,x,z,read(V,x,z,x1,z1));
  }
 }
 states+=bases.size;
 if(values.size)pairs++;
}
chk(pairs===8&&maxGap===2n,"terminal counts");
console.log("PASS:",pairs,"positive-width mask pairs;",states,"productive boundary states; max gap",String(maxGap),witness);
