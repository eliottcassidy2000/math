"use strict";
const eq=(a,b)=>{if(JSON.stringify(a)!==JSON.stringify(b))throw Error(JSON.stringify([a,b]));};
const adj=Array.from({length:50},()=>new Set());
const add=(a,b)=>{adj[a].add(b);adj[b].add(a);};
for(let i=0;i<5;i++)for(let x=0;x<5;x++){
  add(5*i+x,5*i+(x+1)%5);
  add(25+5*i+x,25+5*i+(x+2)%5);
  for(let j=0;j<5;j++)add(5*i+x,25+5*j+(x+i*j)%5);
}
const C=[0,1,2,3,28,8,7,6,44,15,16,17,30];
const D=[4,29,26,11,42,40,23,22,48,46,13,37,39];
const E=[5,9,36,38,35,19,32,10,41,24,45,47,49];
const V=adj.map((_,i)=>i), addBit=i=>1n<<BigInt(i);
const mask=xs=>xs.reduce((m,i)=>m|addBit(i),0n);
const edges=(a,b)=>a.reduce((n,i)=>n+b.filter(j=>adj[i].has(j)).length,0);
for(const X of [C,D,E])for(let i=0;i<13;i++)
  eq(X.filter(j=>adj[X[i]].has(j)).sort((a,b)=>a-b),[X[(i+12)%13],X[(i+1)%13]].sort((a,b)=>a-b));
const paths=[];
const inC=new Set(C);
for(const root of V)if(!inC.has(root)){
  const visit=(X,used)=>{
    const last=X.at(-1);
    for(const s of adj[last]){
      if(s<=root || used.has(s) || inC.has(s))continue;
      if(X.slice(1,-1).some(i=>adj[i].has(s)))continue;
      const closing=X.length>=2&&adj[root].has(s);
      if(X.length===12){if(closing&&X[1]<s)paths.push([...X,s]);}
      else if(!closing){used.add(s);visit([...X,s],used);used.delete(s);}
    }
  };
  visit([root],new Set([root]));
}
const masks=paths.map(mask), pairs=[];
eq(paths.length,968);
const degrees=paths.map(()=>0);
for(let i=0;i<masks.length;i++)for(let j=0;j<i;j++)if(!(masks[i]&masks[j])){
  pairs.push([j,i]);degrees[i]++;degrees[j]++;
}
eq(pairs.length,207);
const hist=X=>Object.entries(X.reduce((out,x)=>(out[x]=(out[x]||0)+1,out),{}));
eq(hist(degrees),[["0",712],["1",176],["2",42],["3",18],["4",12],["5",4],["8",4]]);
function refine(X){
  let c=V.map(()=>0);X.forEach((x,i)=>c[x]=i+1);
  while(new Set(c).size<50){
    const s=V.map(i=>[c[i],...Array.from(adj[i],j=>c[j]).sort((a,b)=>a-b)].join(","));
    const labels=new Map(Array.from(new Set(s)).sort().map((x,i)=>[x,i]));
    const d=s.map(x=>labels.get(x));
    if(labels.size===new Set(c).size)throw Error("nondiscrete refinement");
    c=d;
  }
  return c;
}
const base=refine(C), inv=Array(50);V.forEach(i=>{inv[base[i]]=i;});
function mapsToC(X){
  const maps=[];
  for(const sign of [1,-1])for(let s=0;s<13;s++){
    const src=X.map((_,i)=>X[(s+sign*i+13)%13]);
    const p=refine(src).map(c=>inv[c]);
    if(V.every(i=>Array.from(adj[i]).every(j=>adj[p[i]].has(p[j]))))maps.push(p);
  }
  eq(maps.length,2);
  return maps;
}
const maps=[C,...paths].map(mapsToC);
const cycles=[C,...paths];
function key(pack){
  return pack.flatMap(i=>maps[i].map(p=>pack.filter(j=>j!==i)
    .map(j=>mask(cycles[j].map(v=>p[v])).toString(16)).sort().join(":"))).sort()[0];
}
const pairKeys=new Map(), tripleKeys=new Map();
function take(M,k,X){const v=M.get(k);if(v)v.count++;else M.set(k,{count:1,profile:X});}
for(let i=1;i<cycles.length;i++)take(pairKeys,key([0,i]),[edges(C,cycles[i])]);
for(const [i,j] of pairs){const X=cycles[i+1],Y=cycles[j+1];take(tripleKeys,key([0,i+1,j+1]),
 [edges(C,X),edges(C,Y),edges(X,Y)].sort((a,b)=>a-b));}
eq(hist(Array.from(pairKeys.values(),x=>x.count)),[["2",6],["4",239]]);
eq(hist(Array.from(tripleKeys.values(),x=>x.count)),[["3",1],["6",34]]);
const profiles=Array.from(tripleKeys.values(),x=>x.profile);
eq(new Set(profiles.map(x=>x.join("/"))).size,29);
console.log("Complement edges by triple class",hist(profiles.map(x=>x.reduce((a,b)=>a+b,0)-59)));
const refl=maps[0].find(p=>p.some((x,i)=>x!==i));
eq(V.filter(i=>refl[i]===i),[1,14,18,31,33,34]);
eq(mask(D.map(i=>refl[i])).toString(),mask(E).toString());
console.log("HS pairs",pairKeys.size,"triples",tripleKeys.size,"triple cross-edge profiles",hist(Array.from(tripleKeys.values(),x=>x.profile.join("/"))));
console.log("reflection profile",edges(C,D),edges(C,E),edges(D,E));

const p=61n, n=Number(2n*p), target=(1n<<(p-1n))*((1n<<p)-1n);
const want=[target-2n,target];
const F=[0n,2n];for(let i=2;i<=n;i++)F[i]=1n<<BigInt(i-1);
const low=[];for(let i=2;i<n;i++) low[i]=F[i]-F[n-i]; // 2^(i-1) - f(n-i)
let tested=0, comparisons=0;
const A=[];
function test(h){tested++;if(want.includes(h))throw Error("endpoint: "+A);}
function scan(first,h){
  test(h);
  if(A.length===4){
    for(const t of want){
      let lo=first,hi=n-1;
      while(lo<=hi){
        const a=(lo+hi)>>1;
        const val=h+low[a]-A.reduce((acc,b)=>acc+F[a-b],0n);
        comparisons++;
        if(val===t)throw Error("endpoint: "+A+","+a);
        if(val<t)lo=a+1;else hi=a-1;
      }
    }
    return;
  }
  for(let a=first;a<n;a++){
    const next=h+low[a]-A.reduce((acc,b)=>acc+F[a-b],0n);
    A.push(a);scan(a+1,next);A.pop();
  }
}
scan(2,F[n]); // n is the mandatory last loss
eq(tested,8502671);
console.log("Sparse test",tested,"short candidates;",comparisons,"exact bisection comparisons");
const alphabet=[];
for(let a=1n;(1n<<a)+1n<=target+1n;a+=2n){
  if((target+1n)%((1n<<a)+1n)===0n)throw Error("upper odd block");
  if((target-1n)%((1n<<a)+1n)===0n)alphabet.push(Number(a+2n));
}
eq(alphabet,[3,5,7,17,63]);
console.log("Exact odd blocks",alphabet);

const inds=[0,1,2], mod=x=>(x+26)%13, off=(i,j)=>3*i+j;
function check(steps,S,T,U){
 const a=Array(9);
 for (let i of inds)a[off(i,i)]=[steps[i],13-steps[i]];
 for(const [i,j,D] of [[0,1,S],[0,2,T],[1,2,U]]){
  a[off(i,j)]=D;
  a[off(j,i)]=D.map(d=>mod(-d));
 }
 for(const i of inds)for(let j=i;j<3;j++){
  let seen=0;
  for(const d of a[off(i,j)])seen|=1<<d;
  for(const k of inds)for(const u of a[off(i,k)])for(const v of a[off(k,j)]){
      const d=mod(u+v);
      if(i===j && d===0)continue;
      if(seen>>d&1)return false;
      seen|=1<<d;
  }
 }
 return true;
}
const rows=[];
for(let a=1;a<=6;a++)for(let b=1;b<=6;b++)for(let c=b+1;c<=6;c++){
 if(a===b||a===c)continue;
 let count0=0,count1=0;
 for(let x=1;x<13;x++)for(let y=1;y<13;y++){
  const s=[a,b,c],S=[0,x],T=[0,y];
  if(!check(s,S,T,[]))continue;
  count0++;
  for(let u=0;u<13;u++)if(check(s,S,T,[u]))count1++;
 }
 rows.push([a,b,c,count0,count1]);
}
const hot=rows.filter(row=>row[4]).map(r=>r.slice(0,3));
eq(hot,[[1,4,6],[2,1,5],[3,1,5],[4,2,3],[5,4,6],[6,2,3]]);
const m=rows.find(r=>r[0]===1&&r[1]===4&&r[2]===6);
eq(m,[1,4,6,8,4]);
const seed=[];
for(let x=1;x<13;x++)for(let y=1;y<13;y++)for(let z=0;z<13;z++)
 if(check([1,4,6],[0,x],[0,y],[z]))seed.push([x,y,z]);
console.log("2-2-1 step classes",hot,"normalized seed shifts",seed);
