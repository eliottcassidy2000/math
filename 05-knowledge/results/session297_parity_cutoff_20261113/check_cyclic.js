// Direct, exhaustive chord-mask action. Run: node check_cyclic.js
// Uses only integer operations and does not use the fixed-set formula.
const expected = [[1,1],[1,1],[8,8],[88,16],[2344,256],[131072,1024]];
for (let n=3;n<=8;n++) {
  const E=[];
  for (let i=0;i<n;i++) for (let j=i+1;j<n;j++)
    if (j-i!==1 && j-i!==n-1) E.push([i,j]);
  const D=E.length;
  const P=new Map(E.map(([i,j],k)=>[i*n+j,k]));
  function action(s,reflection) {
    return E.map(([i,j])=>{
      const a=(s+(reflection?-i:i)+n)%n, b=(s+(reflection?-j:j)+n)%n;
      return [P.get(Math.min(a,b)*n+Math.max(a,b)), +(a>b)^+reflection];
    });
  }
  const rotations=Array.from({length:n},(_,i)=>action(i,false));
  const mirror=action(n-1,true);
  const fixed=(x,t)=>t.every(([j,flip],i)=>(((x>>i)&1)^((x>>j)&1))===flip);
  const U={},B={};
  let total=0,self=0;
  const b=n%2?1:2;
  for (let x=0;x<1<<D;x++) {
    const r=rotations.filter(t=>fixed(x,t)).length;
    U[r]=(U[r]||0)+1;
    if (fixed(x,mirror)) B[r]=(B[r]||0)+1;
  }
  for (const r in U) {
    if (U[r]*r%n || (B[r]||0)%b) throw Error('nonintegral orbit count');
    total+=U[r]*r/n;
    self+=(B[r]||0)/b;
  }
  if (total!==expected[n-3][0] || self!==expected[n-3][1]) throw Error('count mismatch');
  console.log(JSON.stringify({n,U,B,total,self}));
}
