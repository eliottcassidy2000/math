const F=x=>(x*x*x-x)/6;
function Q(v){let q=0,p=0; for(const bit of v.toString(2)){if(bit==="1")q+=Math.ceil(p/2);p=2*p+Number(bit);}return q;}
function H(v){return 1662*v-58*Q(v)+427;}
const high=17023575;
let B=1;while (2*F(B+1)+1<=high)B++;
const p3=3*(high-1)+4*B;
let cMax=0;while(16n*BigInt(cMax+1)**3n<=27n*BigInt(p3)**2n)cMax++;
let triples=new Map(), shapes=[];
for(let a=1;3*a*a<=cMax;a+=2)for(let b=a;a*a+2*b*b<=cMax;b+=2)for(let c=b;a*a+b*b+c*c<=cMax;c+=2){
 const n=a*a+b*b+c*c,p=a*b*c, key=1e8*n+p;
 for(const e of [-2,2]){const ts=triples.get(key+e)||[];for(const tr of ts)shapes.push([n,...(e<0?[tr,[a,b,c]]:[[a,b,c],tr])]);}
 // An equal-product triple need not be unique.
 if(!triples.has(key))triples.set(key,[]);triples.get(key).push([a,b,c]);
}
shapes.sort((a,b)=>a[0]-b[0]);
function walsh(s,[u,v,w]){return [[u,v,w],[u,-v,-w],[-u,v,-w],[-u,-v,w]].map(p=>(s+p[0]+p[1]+p[2])/2);}
const targets=new Map();
for(let v=1;1604*v+427<=high;v++)if(v.toString(2).startsWith("1010")&&H(v)+2<=high){
 let h=H(v);targets.set(h,[...(targets.get(h)||[]),v]);
}
const out=[];
for(let [norm,a,b] of shapes)for(let sign of [1,-1]){
 let [a1,b1]=sign===1?[a,b]:[b.map(y=>-y),a.map(y=>-y)];
 for(let s=1;s<=2*B;s+=2){
  let A=walsh(s,a1), B1=walsh(s,b1);
  if(Math.min(...A,...B1)<2)continue;
  let n=A.reduce((x,y)=>x+F(y),0);
  if(2*(n+1)+3>high)break;
  for(let z=2;2*(n+F(z))+3<=high;z++){
   const h=2*(n+F(z))+1;
   if(targets.has(h)){if(B1.reduce((x,y)=>x+F(y),0)!==n+1)throw Error("gap");out.push([h,norm,targets.get(h),[...A,z],[...B1,z]]);}
  }
 }
}
out.sort((a,b)=>a[0]-b[0]);
if(cMax!==163884||shapes.length!==2||out.map(t=>t[0]).join(",")!=="8541511,16818421,17021085,17023573")throw Error("unexpected scan");
console.log({high,B,cMax,distinctChartHeights:targets.size,shapeCount:shapes.length,smallShapes:shapes.filter(t=>t[0]<=45219),hits:out});
