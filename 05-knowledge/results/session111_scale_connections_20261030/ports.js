// Exact finite controls for socle-shift joint-port information.
const need=(p,m)=>{if(!p)throw Error(m);};
const mod=(x,q)=>(x%q+q)%q;
const q=25,p=5,N=q*q;
const add=(a,b)=>[mod(a[0]+b[0],q),mod(a[1]+b[1],q)];
// Coordinates x+y*mu, mu^2=mu+1.
const mul=(a,b)=>[mod(a[0]*b[0]+a[1]*b[1],q),mod(a[0]*b[1]+a[1]*b[0]+a[1]*b[1],q)];
const pow=(a,n)=>{let b=[1,0];for(;n;n=Math.floor(n/2),a=mul(a,a))if(n%2)b=mul(b,a);return b;};
const g=pow([0,1],4), gp=pow(g,5), eta=[mod(gp[0]-1,q),gp[1]];
need(eta[0]!==0||eta[1]!==0,"nonzero eta");
need(mul(eta,eta).every(x=>x===0),"square-zero eta");
const ix=a=>a[0]*q+a[1],xy=i=>[Math.floor(i/q),i%q],rho=a=>mod(a[0]+3*a[1],p);
const orbit=Array(N).fill(-1);let next=0;
for(let i=0;i<N;i++)if(orbit[i]<0){let j=i;do{orbit[j]=next;j=ix(mul(g,xy(j)));}while(j!==i);next++;}
const perms=Array.from({length:p},(_,c)=>Array.from({length:N},(_,i)=>ix(add(xy(i),rho(xy(i))===c?eta:[0,0]))));
function cycles(a){let seen=new Set(),lens=[];for(let i=0;i<a.length;i++)if(!seen.has(i)){let j=i,l=0;do{seen.add(j);l++;j=a[j];}while(j!==i);lens.push(l);}return lens;}
let arcChecks=0;
for(const a of perms){
 const lens=cycles(a);need(lens.filter(x=>x===1).length===500&&lens.filter(x=>x===5).length===25,"cycle inventory");
 for(let i=0;i<N;i++)for(let j=i+1;j<N;j++){
  const u=xy(i),v=xy(j),u1=xy(a[i]),v1=xy(a[j]);
  const d=[mod(v[0]-u[0],q),mod(v[1]-u[1],q)],d1=[mod(v1[0]-u1[0],q),mod(v1[1]-u1[1],q)];
  need(orbit[ix(d)]===orbit[ix(d1)],"literal difference orbit");arcChecks++;
 }
}
function rank(rows){const basis=new Map();for(let v of rows){while(v){let k=v.toString(2).length-1;if(basis.has(k))v^=basis.get(k);else{basis.set(k,v);break;}}}return basis.size;}
const a=perms[0];let ranks=[];
for(let r=1;r<=3;r++){let rows=[];for(let j=0;j<r;j++)for(let i=0;i<N;i++)if(a[i]!==i)rows.push((1n<<BigInt(j*N+i))^(1n<<BigInt(j*N+a[i])));const got=rank(rows);need(got===100*r,"tensor rank");ranks.push(got);}
const hist=(support,coords)=>{const h=new Map();for(const w of support){let z=0;coords.forEach((c,i)=>z|=((w>>c)&1)<<i);h.set(z,(h.get(z)||0)+1);}return [...h].sort((a,b)=>a[0]-b[0]);};
const full=[0,1,2,3,4,5,6,7],parity=[0,3,5,6];
for(const pair of [[0,1],[0,2],[1,2]]){const x=hist(full,pair),y=hist(parity,pair);need(JSON.stringify(x.map(([a,b])=>[a,b/8]))===JSON.stringify(y.map(([a,b])=>[a,b/4])),"pairwise marginal");}
console.log(JSON.stringify({status:"PASS",ring:"(Z/25)[mu]/(mu^2-mu-1)",g,eta,vertices:N,clockOrbitsIncludingZero:next,generatorCount:5,literalPairOrbitChecks:arcChecks,oneFibre:{moved:125,fiveCycles:25,index:100},binaryPortConstraintRanks:ranks,independentTwoPorts:{inputBits:250,freeBits:50,constraintBits:200},correlatedTwoPorts:{effectiveInputBits:125,freeBits:25,constraintBits:100},pairwiseIdenticalTripleControls:{independentExponent:300,evenParityExponent:200},scope:"Verifies universal difference-orbit preservation and constraint ranks; does not enumerate full graph automorphisms."},null,2));
