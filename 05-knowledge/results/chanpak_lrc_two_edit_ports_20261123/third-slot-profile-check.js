function paths(R){const n=R.length,N=1<<n,E=new Uint32Array(N*n),H=new Uint32Array(N);H[0]=1;for(let i=0;i<n;i++)E[(1<<i)*n+i]=1;
 for(let m=1;m<N;m++)for(let i=0;i<n;i++){const q=E[m*n+i];if(q){H[m]+=q;for(let a=R[i]&~m;a;a&=a-1){const b=a&-a;E[(m|b)*n+31-Math.clz32(b)]+=q;}}}return H;
}
function p(R){const H=paths(R),F=H.length-1;const Q=new Uint32Array(F+1);let r=0;for(let S=0;S<=F;S++){for(let A=S;;A=(A-1)&S){Q[S]+=H[A]*H[S^A];if(A===0)break;}r+=Q[S]*H[F^S];}return [H[F],Q[F],r];}
function host(tags,m,R){const f=R.length,k=tags.length-1,n=k+f,a=n,b=n+1, A=Array(n+2).fill(0);
 for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)if(!(m<=i&&j<m+f))A[i]|=1<<j;
 for(let i=0;i<f;i++)A[m+i]|=R[i]<<m;
 for(let i=0;i<n;i++){const t=tags[i<m?i:i<m+f?m:i-f+1];for(let j=0;j<2;j++)if(t>>j&1)A[n+j]|=1<<i;else A[i]|=1<<(n+j);}
 A[b]=A[b]|1<<a;return A;
}
const entries=[['000210230',6,[1,5],[[155,102,20],[935,606,116],[335,222,44],[1085,714,140]]],['0003311301',7,[0,2],[[675,524,122],[3491,2672,610],[2083,1598,366],[4725,3668,854]]]];
function eq(a,b){if(JSON.stringify(a)!==JSON.stringify(b))throw Error(JSON.stringify([a,b]));}
function audit(log){
 for(const [tags,m,marks,expected] of entries){const W=[...tags].map(Number),out=[];
  for(let c=0;c<4;c++){let B=W.slice();for(let i=0;i<2;i++)if(c>>i&1)B[marks[i]>>1]^=1<<(marks[i]&1);
   const P=[1,2,3].map(r=>paths(host(B,m,Array(r).fill(0))).at(-1)/[1,1,2,6][r]);out.push(P);}
  eq(out,expected);log('independent skeleton slot patterns',out);
 }
 const F=[[[86,284,392,401,484,79,398,291,33],[2393,38272,317403],[7981807,46862827,17465251,55872649]],[[6,244,344,51,417,69,153,301,107],[2393,38272,317643],[7986607,46890667,17475811,55906249]]];
 for(const [R,profile,want] of F){eq(p(R),profile);const [tags,m,marks]=entries[0],W=[...tags].map(Number),A=host(W,m,R), out=[];
  for(let c=0;c<4;c++){const B=A.slice();for(let i=0;i<2;i++)if(c>>i&1){const u=marks[i]>>1,v=B.length-1;B[u]^=1<<v;B[v]^=1<<u;}out.push(paths(B).at(-1));}
  for(let j=2;j*j<=out[0];j++)if(out[0]%j===0)throw Error('composite');
  eq(out,want);log('physical profile',profile,out,(out[3]-out[1]-out[2]+out[0])/2);
 }
 eq(p([94,60,184,112,32,193,150,27]),[345,5640,45225]);eq(p([30,92,216,48,32,199,153,27]),[345,5640,44505]);
}

function canonical(R){
 const n=R.length;
 function vClass(v,C){return C.map(G=>G.reduce((x,j)=>x+(R[v]>>j&1),0)).join(',');}
 let C=[R.map((_,i)=>i)];
 for(;;){const Q=[];for(const G of C){const M=new Map();for(const v of G){const x=vClass(v,C);if(!M.has(x))M.set(x,[]);M.get(x).push(v);}for(const k of [...M.keys()].sort())Q.push(M.get(k));}
  if(Q.length===C.length)break;C=Q;
 }
 let best=2**(n*(n-1)/2),W=[];
 function group(i){
  if(i===C.length){let k=0;for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)k=2*k+(R[W[i]]>>W[j]&1);best=Math.min(best,k);return;}
  function perm(U){if(!U.length){group(i+1);return;}for(let j=0;j<U.length;j++){W.push(U[j]);perm(U.slice(0,j).concat(U.slice(j+1)));W.pop();}}
  perm(C[i]);
 }
 group(0);return best;
}
function profiles(log){
 let A=[[0]],set=new Set(),triples=new Map(), collisions=[];
 for(let n=1;n<=8;n++){
  const tally=new Set();
  for(const v of A){const [d,q,r]=p(v),key=[d,q].join();tally.add(key);if(!triples.has(key))triples.set(key,new Set());triples.get(key).add(r);}
  eq(A.length,[0,1,1,2,4,12,56,456,6880][n]);if(n>=7)eq(tally.size,[156,1727][n-7]);log(n,A.length,tally.size);
  if(n===8)break;
  const M=new Map();
  for(const v of A)for(let m=0;m<1<<n;m++){const R=v.map((a,i)=>a|((!(m>>i&1))<<n)).concat(m),key=canonical(R);if(!M.has(key))M.set(key,R);}
  A=[...M.values()];
 }
 for(const [key,ss] of triples)if(ss.size>1)collisions.push([key,[...ss].sort((a,b)=>a-b)]);
 eq(collisions,[['345,5640',[44505,45225]]]);log('all cross-order collisions',collisions);
}

audit(console.log);profiles(console.log);
