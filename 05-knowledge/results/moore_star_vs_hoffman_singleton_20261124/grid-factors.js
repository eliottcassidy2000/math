// node grid-factors.js
const m=(x,n=5)=>(x%n+n)%n;
const pairs=(n)=>Array.from({length:n},(_,i)=>Array.from({length:n-i-1},(_,k)=>[i,i+k+1])).flat();
const hs=Array.from({length:50},()=>new Set());
function add(A,i,j){A[i].add(j);A[j].add(i)}
for(let a=0;a<5;a++)for(let x=0;x<5;x++){
 add(hs,5*a+x,5*a+m(x+1));add(hs,25+5*a+x,25+5*a+m(x+2));
 for(let b=0;b<5;b++)add(hs,5*a+x,25+5*b+m(x+a*b));
}
function offsets(i,j){
 const d=Math.min(Math.abs(i-j),10-Math.abs(i-j));
 if(d===1){const s=1+(j===(i+1)%10?i:j)%2;return [s,m(-s)]}
 return d===5||d===(i%2===0?4:2)?[0]:[];
}
const right=Array.from({length:50},()=>new Set());
for(let i=0;i<10;i++)for(let j=0;j<10;j++)for(const a of offsets(i,j))for(let k=0;k<5;k++)add(right,5*i+k,5*j+m(k+a));
function moore(A){
 if(A.some(s=>s.size!==7))throw Error('not seven-regular');
 for(const [i,j] of pairs(A.length)){
  const c=[...A[i]].filter(x=>A[j].has(x)).length;
  if(c!==1-Number(A[i].has(j)))throw Error('not Moore');
 }
}
moore(hs);moore(right);
const E=pairs(50).filter(([i,j])=>!hs[i].has(j)).map(([i,j])=>[i,j,0]);
const F=pairs(10).flatMap(([i,j])=>[0,1,2,3,4].filter(a=>!offsets(i,j).includes(a)).map(a=>[i,j,a]));
for(let i=0;i<10;i++)for(let a=1;a<=2;a++)F.push([i,i,a]);
function factor(N,ed){
 const inc=Array.from({length:N},()=>[]), out=inc.map(()=>[]),seen=new Set();
 ed.forEach(([i,j],k)=>{inc[i].push(k);inc[j].push(k)});
 for(let i=0;i<N;i++){
  let u=i;
  while(true){
   while(inc[u].length&&seen.has(inc[u].at(-1)))inc[u].pop();
   if(!inc[u].length)break;
   const k=inc[u].pop();seen.add(k);
   const [a,b,r]=ed[k];
   if(a===u){out[u].push([b,r]);u=b}else{out[u].push([a,m(-r)]);u=a}
  }
  if(u!==i)throw Error('Euler trail unclosed');
 }
 if(seen.size!==ed.length||out.some(a=>a.length!==21))throw Error('orientation');
 const fs=[];
 for(let k=0;k<21;k++){
  const owner=Array(N).fill(-1), sel=Array(N).fill(-1);
  function match(i,seen){
   for(let j=0;j<out[i].length;j++){
    const v=out[i][j][0];if(seen[v])continue;seen[v]=true;
    if(owner[v]===-1||match(owner[v],seen)){owner[v]=i;sel[i]=j;return true}
   }
   return false;
  }
  for(let i=0;i<N;i++)if(!match(i,Array(N).fill(false)))throw Error('Hall');
  fs.push(out.map((v,i)=>v.splice(sel[i],1)[0]));
 }
 return fs;
}
const a=factor(50,E), b=factor(10,F);
const C=Array.from({length:2500},()=>new Set());
for(let i=0;i<21;i++)for(let f=0;f<50;f++)for(let h=0;h<50;h++){
 const [j,v]=b[i][Math.floor(h/5)];
 add(C,50*f+h,50*a[i][f][0]+5*j+m(h+v));
}
const t=x=>x-x%5+m(x+1);
for(let x=0;x<2500;x++){
 const f=Math.floor(x/50),h=x%50,L=new Set(),R=new Set();
 for(const y of C[x]){
  L.add(Math.floor(y/50));R.add(y%50);
  if(!C[t(x)].has(t(y)))throw Error('equivariance');
 }
 if(C[x].size!==42||L.size!==42||R.size!==42||L.has(f)||R.has(h)||[...L].some(y=>hs[f].has(y))||[...R].some(y=>right[h].has(y)))throw Error('projection');
}
let tri=0, square=0;
for(let x=0;x<2500;x++){
 const c=new Map();
 for(const y of C[x])for(const z of C[y])if(z>x)c.set(z,1+(c.get(z)||0));
 for(const [y,n] of c){if(C[x].has(y))tri+=n;square+=n*(n-1)/2}
}
let good=0;
for(let r=0;r<10;r+=2)for(let a=0;a<5;a++){
 if(right.every((s,h)=>{const k=5*m(Math.floor(h/5)+r,10)+m(h+a);return k!==h&&!s.has(k)}))good++;
}
if(good!==20)throw Error('automorphism filter');
console.log('2500 rows, 52500 edges, both coordinate bijections and C5 symmetry verified');
console.log('triangles',tri/3,'squares',square/2,'admissible automorphisms',good);

// The smaller model: two of the five anticomplete pentagons are the grid
// coordinates and the other three will be reconstructed from their blocks.
const Dg=hs.slice(0,25).map(s=>new Set([...s].filter(x=>x<25)));
const blocks=Array.from({length:15},(_,i)=>new Set([...hs[35+i]].filter(x=>x<25)));
const has=(i,j)=>blocks[i].has(j);
const intersection=(a,b)=>[...a].filter(x=>b.has(x)).length;
const ds=Array.from({length:15},(_,i)=>new Set(Array.from({length:15},(_,j)=>j)
 .filter(j=>j!==i&&!intersection(blocks[i],blocks[j]))));
for(const [x,y] of pairs(25)){
 const xy=(x%5)!==(y%5)&&m(x%5+Math.floor(x/5))!==m(y%5+Math.floor(y/5));
 const q=Number(xy&&!Dg[x].has(y)&&!intersection(Dg[x],Dg[y]));
 if(blocks.filter(b=>b.has(x)&&b.has(y)).length!==q)throw Error('missing pairs');
}
function solve(A){
 const p=101,n=A.length;
 for(let i=0;i<n;i++){
  const j=A.findIndex((r,j)=>j>=i&&m(r[i],p));if(j<0)throw Error('singular');
  [A[i],A[j]]=[A[j],A[i]];
  const c=Array.from({length:100},(_,j)=>j+1).find(x=>m(x*A[i][i],p)===1);
  A[i]=A[i].map(x=>m(x*c,p));
  for(let k=0;k<n;k++)if(k!==i){const s=A[k][i];A[k]=A[k].map((x,j)=>m(x-s*A[i][j],p))}
 }
 return A.map(r=>r[n]);
}
for(let j=0;j<15;j++){
 const comp=[...ds[j]]; // each disjointness component is K5 here
 const Y=Array.from({length:25},(_,x)=>1-Number(has(j,x))-[...Dg[x]].filter(y=>has(j,y)).length);
 const sol=solve(comp.map(i=>[...comp.map(k=>i===k?4:-Number(ds[i].has(k))),
  [...blocks[i]].reduce((a,x)=>a+Y[x],0)-2]));
 for(let i=0;i<15;i++)if((sol[comp.indexOf(i)]||0)!==Number(hs[35+i].has(35+j)))throw Error('root');
}
console.log('all 15 pinned Laplacian solves recover the three leftover Hoffman–Singleton pentagons');
