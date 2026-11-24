function* part(n,b=1) {
  if (!n) {yield []; return;}
  for (let a=b; a<=n; a++) for (const t of part(n-a,a)) yield [a,...t];
}
const ff=(k,x)=>(1n<<BigInt(k*x-2))+1n;
const verify=x=>{if(!x)throw Error("verification failed");};
for(const [p,n,f,end,k,want] of [[11,5,0,false,5,36],[5,11,0,false,5,685],[5,10,5,false,5,444],[7,8,0,true,3,74]]) {
  const S=new Map();
  for(const a of part(n)) for(let t=0;t<2**a.length;t++) {
    const z=new Array(f).fill(1);
    a.forEach((x,i)=>z.push(...((t>>i&1)?[p*x]:new Array(p).fill(x))));
    if(end&&z.includes(1))continue;
    z.sort((a,b)=>a-b);
    S.set(z.join(","),z.reduce((h,x)=>h*ff(k,x),1n));
  }
  const H=new Set(S.values());
  verify(S.size===want&&H.size===want);
  const sorted=[...H].sort((a,b)=>a<b?-1:(a>b?1:0));
  if(p===11)verify(sorted.slice(0,3).join()===[9n**55n,9n**44n*ff(5,11),9n**33n*257n**11n].join());
  console.log(p,n,f,S.size,H.size);
}
const mod=(x,p)=>(x%p+p)%p;
function rank(A,p) {
  A=A.map(a=>a.map(x=>mod(x,p)));
  let z=0;
  for(let j=0;j<A.length;j++) {
    let t=z;
    while(t<A.length&&!A[t][j])t++;
    if(t===A.length)continue;
    [A[z],A[t]]=[A[t],A[z]];
    let v=1;
    while(mod(v*A[z][j],p)!==1)v++;
    for(let i=z+1;i<A.length;i++) {
      const c=A[i][j]*v;
      for(let k=j;k<A.length;k++)A[i][k]=mod(A[i][k]-c*A[z][k],p);
    }
    z++;
  }
  return z;
}
const arr=(f,edge,d,e,j)=>Array.from({length:f},(_,u)=>Array.from({length:f},(_,v)=>d*(u===v)+e*edge(u,v)+j));
const st=(i,j)=>i!==j&&(i===0||j===0), pent=(i,j)=>Math.abs(i-j)===1||Math.abs(i-j)===4;
verify(rank(arr(56,st,455,-65,1),2)===54);
verify((7*56+50-112)/15===22&&(7*50+50-700)/15===-20);
console.log("p=2 star −8 rank 54, ordinary trace 22; HS ordinary trace −20");
for(const [p,f,edge,r] of [...Array.from({length:9},(_,i)=>[7,2+7*i,st,7*i]),[11,5,pent,2],[13,0,pent,0],[19,1,()=>false,0]]) {
  verify(rank(arr(f,edge,400,50,-1),p)===r);
  const a=15*r-8*f+65,m=(p-1)*a/(2*p);
  verify(Number.isInteger(m));
  console.log(p,f,r,a,m);
}
