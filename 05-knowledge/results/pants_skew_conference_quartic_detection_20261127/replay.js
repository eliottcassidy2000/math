// Exact small controls; run with Node.js.
const assert=require("node:assert/strict");
const pop=x=>x===0?0:1+pop(x&(x-1));
const fac=n=>n<2?1:n*fac(n-1);
const sum=a=>a.reduce((x,y)=>x+y,0);
function paley(q) {
  const squares=new Set(Array.from({length:q-1},(_,i)=>(i+1)**2%q));
  return Array.from({length:q+1},(_,i)=>Array.from({length:q+1},(_,j)=>
    i===j?0:i===q?1:j===q?-1:squares.has((j-i+q)%q)?1:-1));
}
const vertices=(n,s)=>Array.from({length:n},(_,i)=>i).filter(i=>s>>i&1);
function ham(A,V,flip) {
  const l=V.length,N=1<<l,M=new Uint16Array(l*N);
  for(let i=0;i<l;i++) M[(1<<i)*l+i]=1;
  for(let s=1;s<N;s++)for(let v=0;v<l;v++){
    const c=M[s*l+v];if(!c)continue;
    for(let w=0;w<l;w++)if(!(s>>w&1)&&
      A[V[v]][V[w]]===((flip>>v^flip>>w)&1?-1:1))
      M[(s|(1<<w))*l+w]+=c;
  }
  return sum([...M.subarray((N-1)*l)]);
}
function walsh(array) {
  for(let w=1;w<array.length;w*=2)for(let a=0;a<array.length;a+=2*w)
    for(let j=a;j<a+w;j++){
      const x=array[j],y=array[j+w];array[j]=x+y;array[j+w]=x-y;
    }
  return array;
}
function check(n) {
  const A=paley(n-1), N=1<<n;
  for(let i=0;i<n;i++)for(let j=0;j<n;j++)
    assert.equal(sum(A[i].map((x,k)=>x*A[j][k])),i===j?n-1:0);
  const quartic=new Map();let mass=0,pc=0;
  for(let s=0;s<N;s++)if(pop(s)===4){
    const V=vertices(n,s),[a,b,c,d]=V;
    const p=Math.abs(A[a][b]*A[c][d]-A[a][c]*A[b][d]+A[a][d]*A[b][c]);
    const q=sum(A.map((_,x)=>V.reduce((w,i)=>w*A[i][x],1)));
    pc+=+(p===3);mass+=q*q;quartic.set(s,q);
    const r=n-4+(p===1?4:0);
    assert.equal((q-r)%8||0,0);
    assert(Math.abs(q)<=n-(p===1?8:4));
    assert.equal(Math.abs(q),((n===8&&p===3)||(n===12&&p===1))?4:0);
  }
  assert.equal(pc,n*n*(n-1)*(n-2)/96);
  assert.equal(mass,n*n*(n-1)*(n-2)*(n-4)/24);
  for(let l=3;l<=6;l++){
    const h=1<<(l-1), F=new Array(N).fill(0);
    for(let s=0;s<N;s++)if(pop(s)===l){
      const V=vertices(n,s);
      const W=walsh(Array.from({length:h},(_,flip)=>ham(A,V,flip)));
      for(let i=0;i<h;i++){
        let key=0;
        for(let j=0;j<l-1;j++)if(i>>j&1)key|=1<<V[j];
        if(pop(i)&1)key|=1<<V[l-1];
        F[key]+=fac(n-l)*W[i];
      }
    }
    assert.equal(F[0],fac(n));
    for(let s=1;s<N;s++)assert.equal(4*F[s],
      l===6&&quartic.get(s)?-3*fac(n-6)*h*quartic.get(s):0);
  }
  console.log("order",n,"Pfaffian-three four-sets",pc,
              "quartic squared mass",mass,"all run spectra checked");
}
check(8); check(12);
