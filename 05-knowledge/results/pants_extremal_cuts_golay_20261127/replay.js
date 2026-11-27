// Exact replay; run with Node.js.
const assert=require("node:assert/strict");
const range=n=>Array.from({length:n},(_,i)=>i);
const pop=x=>x===0?0:1+pop(x&(x-1));
const sum=a=>a.reduce((x,y)=>x+y,0);
const fact=n=>n<2?1:n*fact(n-1);
function data(p) {
  const sq=new Set(range(p-1).map(i=>(i+1)**2%p));
  const mod=a=>((a%p)+p)%p;
  const chi=a=>sq.has(mod(a))?1:-1;
  const A=range(p+1).map(i=>range(p+1).map(j=>i===j?0:i===p?1:j===p?-1:chi(j-i)));
  return {A,chi,mod};
}
const mr=(A,i)=>A[i].reduce((v,a,j)=>a<0?v|1<<j:v,0);
function pathHist(A) {
  const n=A.length,N=1<<(n-1),tot=new Int32Array(N);
  for(let s=0;s<1<<n;s++)if(pop(s)===6){
    const V=range(n).filter(i=>s>>i&1);
    const count=[];
    for(let x=0;x<32;x++){
      const D=new Uint16Array(64*6);
      for(let i=0;i<6;i++)D[(1<<i)*6+i]=1;
      for(let S=1;S<64;S++)for(let v=0;v<6;v++){
        const a=D[S*6+v];if(!a)continue;
        for(let w=0;w<6;w++)if(!(S>>w&1)&&
          A[V[v]][V[w]]===(((x>>v^x>>w)&1)?-1:1))
          D[(S|1<<w)*6+w]+=a;
      }
      count.push(sum([...D.subarray(378)]));
    }
    for(let x=0;x<N;x++){
      let w=0;for(let i=0;i<5;i++)w|=((x>>V[i]^x>>V[5])&1)<<i;
      tot[x]+=count[w];
    }
  }
  return tot;
}
for(const n of [8,12]){
  const {A}=data(n-1),N=1<<(n-1), all=(1<<n)-1;
  const rows=range(n).map(i=>mr(A,i)),actual=pathHist(A),hist={};
  const occ=new Array(3).fill(0),L=n*(n-2);
  for(let x=0;x<N;x++){
    const y=rows.map((r,i)=>(n-1-2*pop((x^r)&(all^1<<i)))**2-1);
    assert.equal(sum(y),L);
    const s=sum(y.map(v=>v*v)),T=s+n+2*L;
    const d=(n-1-Math.sqrt(1+Math.max(...y)))/2;
    const F=(fact(n)/fact(n-6)+n*(3*(n-1)**2-2*(n-1))-T)/32;
    assert.equal(F,actual[x]);
    if(d<=2){
      occ[d]++;
      assert.equal(s,[L*L,((n-2)*(n-4))**2+32*(n-2),
        ((n-4)*(n-6))**2+144*(n-4)+16*n][d]);
    }
    const val=String(F*fact(n-6));hist[val]=(hist[val]??0)+1;
  }
  console.log("order",n,"shell sizes",occ.join(","),"full histogram",JSON.stringify(hist));
}
const {A,chi,mod}=data(23), full=(1<<24)-1;
const rows=range(24).map(i=>mr(A,i));
assert.equal(rows[0],0xfacca0);
assert(rows.every(r=>pop(r)%2===0));
const M=new Set([...rows,...rows.map(r=>r^full)]);
assert.equal(M.size,48);
const basis=[];
for(let x of [full,...rows]){
  while(x){ const i=31-Math.clz32(x);
    if(basis[i])x^=basis[i];else {basis[i]=x;break;}}
}
let C=[0];
for(const x of basis.filter(x=>x))C.push(...C.map(y=>x^y));
assert.equal(C.length,4096);
const set=new Set(C), W={};
for(const x of C)W[pop(x)]=(W[pop(x)]??0)+1;
assert.deepEqual(W,{0:1,8:759,12:2576,16:759,24:1});
for(let i=0;i<24;i++)for(let j=i+1;j<24;j++)for(let k=j+1;k<24;k++)
  for(let l=k+1;l<24;l++){
    const S=[i,j,k,l];
    const Q=sum(A.map((_,x)=>S.reduce((p,y)=>p*A[y][x],1)));
    const pp=Math.abs(A[i][j]*A[k][l]-A[i][k]*A[j][l]+A[i][l]*A[j][k]);
    // Golay has no weight four, so the intrinsic weight-four test cannot trigger.
    const h=A[i][j]*A[i][k]*A[i][l];
    assert(!(pp===3&&Q===20*h));
  }
for(const name of ["translation","inversion"]){
  const inv=x=>range(23).find(y=>mod(x*y)===1);
  const p=range(24).map(x=>name==="translation"?x===23?23:mod(x+1):
    x===0?23:x===23?0:mod(-inv(x)));
  const P=x=>range(24).reduce((m,i)=>x>>i&1?m|1<<p[i]:m,0);
  const t=range(24).reduce((m,i)=>name==="inversion"&&i!==23&&chi(i===0?-1:-i)<0?
    m|1<<p[i]:m,0);
  assert(name==="translation"?t===0:t===0xfacca0);
  assert(set.has(t));
  assert([...M].every(x=>M.has(P(x)^t)));
  assert(C.every(x=>set.has(P(x))));
  console.log(name,"Golay action, affine shift",t.toString(16));
}
console.log("Golay weight distribution",JSON.stringify(W));
