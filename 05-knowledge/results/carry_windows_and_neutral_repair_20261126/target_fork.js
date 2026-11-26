
const {cores,F}=require("../overnight_chanpak_moore_20261126/commuting_unit_forks");
const {roots,inside,inverseDouble,inversePrefix}=require("../overnight_chanpak_moore_20261126/nested_fork_inverse");
const all=[...cores,
 [[2,72,3,21,39,15,61],5,[4,6],59,87,11],
 [[2,68,17,7,43,15,61],6,[3,4],71,72,14],
 [[34,20,9,67,13,30,55],3,[4,5],103,129,11],
 [[36,65,42,19,7,90,29],6,[0,4],143,173,23]
];
const targets=[13,27,289,889,1489,9011,47184047,53084231,34379738367];
const data=Object.fromEntries(targets.map(x=>[x,{trueCenter:0,corner:0,examples:[]}]));
const G=(i,c,R,A,D)=>R*c+A*F(c-1n)-D;
const seen=new Set();
for(const [n,T] of roots(Number(process.argv[2]||10))){
 const [h,e,r,V,Q]=inside(T);
 for(let j=0;j<all.length;j++){
  const [,,,f,g,k]=all[j],H=BigInt(f),R=BigInt(g),E=BigInt(k);
  const A=H*h,B=H*r+R*Q,D=H*e+R*V+E,D0=R*V+E;
  const key=[j,A,B,D,D0].join();if(seen.has(key))continue;seen.add(key);
  for(const N of targets)for(const c0 of [N,N-2,N+2])for(const s of [1n,-1n]){
   if(c0<1)continue;
   const want=BigInt(c0),shift=s===1n?0n:2n;
   const z0=inverseDouble(want+D+shift+s,A,B), z1=inversePrefix(want+D0+shift+s,A,B);
   for(const row of [...z0.out,...z1.out]){
    if(BigInt(row[2])!==s)continue;
    const [cRaw,i,_,suffix,k0]=row, c=BigInt(cRaw), L=1n<<BigInt(i), B1=B<<BigInt(k0||0), D1=(k0?D0:D)+shift;
    if(G(i,c,B1,A,D1)!==want||G(i,c+L,B1,A,D1)!==want+2n*s)throw Error("fork control");
    const y=data[N];if(c0===N)y.trueCenter++;else y.corner++;
    if(y.examples.length<8)y.examples.push({c0,n,T,j,c:cRaw,k:k0||0,i,s:String(s),A:String(A),B:String(B1),D:String(D1),order:all[j][0].length+n+(k0||0)+1+(c-1n).toString(2).length});
   }
  }
 }
}
console.log(JSON.stringify({distinct:seen.size,data},null,2));
