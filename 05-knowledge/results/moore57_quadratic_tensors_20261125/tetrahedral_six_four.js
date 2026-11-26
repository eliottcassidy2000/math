// Tetrahedral Wilson action, small tournaments and prime-cut replay.

const I=[0,1,2],Ps=[[0,1,2],[0,2,1],[1,0,2],[2,1,0],[1,2,0],[2,0,1]];
function mul(f,g){return f.map(x=>g[x]);}
function inv(f){return I.map(x=>f.indexOf(x));}
function eq(f,g){return f.filter((x,i)=>x===g[i]).length;}
function charts(){
 const count={},geom={};let best=Infinity,wit=[];
 for(const a of Ps)for(const b of Ps)for(const c of Ps){
   const P=Array.from({length:4},()=>Array(4));
   for(let i=1;i<4;i++)P[0][i]=P[i][0]=I;
   for(const [i,j,v] of [[1,2,a],[1,3,b],[2,3,c]]){P[i][j]=v;P[j][i]=inv(v);}
   let t=0,s=0;
   for(let i=0;i<4;i++)for(let j=i+1;j<4;j++){
    const [k,l]=[0,1,2,3].filter(x=>x!==i&&x!==j),
          x=mul(P[i][k],P[k][j]),y=mul(P[i][l],P[l][j]);
    t+=eq(P[i][j],x)+eq(P[i][j],y);s+=eq(x,y);
   }
   if(t%3||s%2)throw Error("cycle");
   t/=3;s/=2;
   const E=6*t+4*s;
   count[E]=(count[E]||0)+1;
   const ts=[t,s].join(",");
   geom[ts]=(geom[ts]||0)+1;
   if(E<best){best=E;wit=[];}if(E===best)wit.push([a,b,c,t,s,[eq(a,I),eq(b,I),eq(c,I),eq(mul(mul(a,c),inv(b)),I)].findIndex(z=>z===3)]);
 }
 return {count,geom,best,wit};
}
function primeControl(){

function perms(a){if(!a.length)return [[]];return a.flatMap((x,i)=>perms(a.filter((_,j)=>i!==j)).map(r=>[x,...r]));}
const pairs=[];
for(let i=0;i<5;i++)for(let j=i+1;j<5;j++)pairs.push([i,j]);
const adj=(a,b)=>!a.some(v=>b.includes(v));
const E=pairs.filter(x=>!x.includes(4)),F=pairs.filter(x=>x.includes(4));
const p=[3,13,23,2,7,29,5,11,17,19], S=new Set([4,9,16,25,36,49]);
const Q=[];for(let i=0;i<10;i++)for(let j=i+1;j<10;j++)if(S.has(p[i]+p[j]))Q.push([i,j]);
let best=-1,N=0;
for(const e of perms(E))for(const f of perms(F)){
 const V=[...e,...f];
 const c=Q.filter(([i,j])=>adj(V[i],V[j])).length;
 if(c>best){best=c;N=0;}if(c===best)N++;
}
let example=null;
function embed(V, unused){
 if(!unused.length)return V;
 const i=V.length;
 for(const pair of unused){
  if(!Q.every(([a,b])=>b!==i||adj(V[a],pair)))continue;
  const w=embed([...V,pair],unused.filter(e=>e!==pair));
  if(w)return w;
 }
 return null;
}
example=embed([],pairs);
const T={};
for(let m=0;m<64;m++){
 const e=[]; for(let i=0;i<4;i++)for(let j=i+1;j<4;j++)e.push([i,j]);
 const A=[[],[],[],[]];
 for(let i=0;i<6;i++){let [a,b]=e[i];if(m>>i&1)[a,b]=[b,a];A[a].push(b);}
 let c=0,H=0;
 for(let a=0;a<4;a++)for(const b of A[a])for(const q of A[b])if(A[q].includes(a))c++;
 c/=3;
 for(const w of perms([0,1,2,3]))if(w.slice(1).every((x,i)=>A[w[i]].includes(x)))H++;
 if(H!==1+2*c)throw Error("H");
 const t=[A.map(r=>r.length).sort().join(""),c,H].join(",");
 T[t]=(T[t]||0)+1;
}
return {best,N,Q,example,T};

}
const W=charts();
const S=primeControl();

const expect={"18":8,"20":36,"24":12,"32":36,"36":18,"38":72,"48":12,"60":3,"68":18,"108":1};
if(JSON.stringify(W.count)!==JSON.stringify(expect))throw Error('energy census');
if(W.best!==18||W.wit.length!==8)throw Error('minimizers');
const V=[0,0,0,0];
for(const x of W.wit)V[x.at(-1)]++;
if(V.join()!=="2,2,2,2")throw Error('face support');
if(S.best!==3||S.N!==1152||!S.example||S.Q.length!==7)throw Error('prime split');
if(Object.values(S.T).sort((a,b)=>a-b).join()!=="8,8,24,24")throw Error('tournaments');
console.log("K4 rooted chart energy counts:",W.count,"; minimum bad-face multiplicities:",V);
console.log("4-vertex tournaments keyed by score / directed triangles / Hamilton paths:",S.T);
console.log("Prime split-preserving maximum:",S.best,"; maximizers:",S.N,
            "; unrestricted embedding:",S.example);
