// node minimum_trades.js ../moore57_quadratic_tensors_20261125/period_fourteen_site.txt
const need=(c,m)=>{if(!c)throw Error(m)};
const N=56, mod=x=>((x%N)+N)%N, ident=Array.from({length:N},(_,x)=>x);
const slots=ident.filter(h=>h!==0&&h!==28), centers=[56,0,28];
const triples=[[9,10,19],[8,12,20],[7,11,18],[6,15,21],[5,22,27],[4,13,17],[3,23,26],[2,14,16],[1,24,25]];
const inverse=x=>{const a=[];x.forEach((y,i)=>a[y]=i);return a;};
const t=x=>ident.map(y=>mod(y+x));
function chart(trades) {
 const s=[], p=[],q=[], C=Array.from({length:57},()=>[]);
 for(let [a,b,c] of triples){const d=[a,b,-c];for(let i=0;i<3;i++)for(let m of [-1,1])s[mod(m*d[i])]=mod(-m*d[(i+1)%3]);}
 for(let h of slots){p[h]=t(-h);q[h]=t(-s[h]);}
 for(let [h,k,x] of trades){
  let I=inverse(p[k]),J=inverse(q[k]),seen=new Set([x]),stack=[x];
  while(stack.length){let z=stack.pop();for(let y of [I[p[h][z]],J[q[h][z]]])if(!seen.has(y)){seen.add(y);stack.push(y);}}
  need(seen.size<N,"whole trade orbit");
  for(let j of seen){[p[h][j],p[k][j]]=[p[k][j],p[h][j]];[q[h][j],q[k][j]]=[q[k][j],q[h][j]];}
 }
 for(let h of ident)C[56][h]=C[h][56]=ident;
 C[0][28]=C[28][0]=t(28);
 for(let h of slots){C[h][0]=p[h];C[0][h]=inverse(p[h]);C[h][28]=q[h];C[28][h]=inverse(q[h]);}
 return C;
}
function boards(file){
 const a=file.trim().split(/\r?\n/), nt=Number(a[0].split(/\s+/)[1]);
 const rows=l=>l.trim().split(/\s+/).map(Number);
 const trades=a.slice(1,nt+1).map(rows);
 const seed=[a[nt+1].trim().split(/\s+/).slice(1).map(Number),...a.slice(nt+2,nt+55).map(rows)];
 const before=chart(trades),after=chart([...trades,[10,45,2],[27,49,0]]);
 const pos=(C,h,x)=>centers.map(c=>C[h][c][x]);
 const key=x=>x.join(","), lookup=new Map();
 for(let h of slots)for(let x of ident)lookup.set(key(pos(after,h,x)),[h,x]);
 const placed=[], read=seed.map(([h,x])=>pos(before,h,x));
 for(let d=0;d<N;d+=4){
  const image=read.map(z=>lookup.get(key(z.map(y=>mod(y+d)))));
  need(image.every(Boolean)&&new Set(image.map(([h])=>h)).size===54,"planted arrangement");
  placed.push(image);
 }
 need(placed.length===14,"fixture size");
 const full=[],marked=[], marker=pos(after,2,7),x=1,b=1;
 const reserved=centers.map((c,j)=>new Set([after[b][c][x],...centers.filter(k=>k!==c).map(k=>after[k][c][after[b][k][x]]),marker[j]]));
 const used=Array.from({length:56},()=>new Set());
 for(let star of placed)for(let [h,y] of star.slice(1))used[h].add(y);
 for(let h of slots)for(let y of ident) {
  let z=pos(after,h,y); full.push([h,...z]);
  if(h!==1&&h!==2&&!used[h].has(y)&&z.every((v,i)=>!reserved[i].has(v)))marked.push([h,...z]);
 }
 need(full.length===3024&&marked.length===1741,"board sizes");
 return {full,marked};
}

function catalogue(rows){
 const n=rows.length, stride=56*56;
 const ids=new Int32Array(16*stride).fill(-1);
 const idx=(i,j,x,y)=>(4*i+j)*stride+56*x+y;
 const get=(i,j,x,y)=>ids[idx(i,j,x,y)];
 for(let a=0;a<n;a++)for(let i=0;i<4;i++)for(let j=0;j<4;j++)if(i!==j){
  const p=idx(i,j,rows[a][i],rows[a][j]);need(ids[p]<0,"not linear");ids[p]=a;
 }
 const ds=[];
 for(let i=1;i<4;i++)for(let j=0;j<4;j++)if(j!==1&&j!==i)
 for(let k=0;k<4;k++)if(k!==2&&k!==j&&k!==i){
  let l=6-i-j-k;if(l!==3)ds.push([i,j,k,l]);
 }
 need(ds.length===9,"derangements");
 const M=new Int32Array(16), L=new Int32Array(4),R=new Int32Array(4);
 const qs=new Map(), coverage=new Set(), examples=[];
 let candidates=0,total=0;
 for(let a=0;a<n;a++)for(let b=a+1;b<n;b++) {
  const X=rows[a],Y=rows[b];
  if(X[0]===Y[0]||X[1]===Y[1]||X[2]===Y[2]||X[3]===Y[3])continue;
  for(let i=0;i<4;i++)for(let j=0;j<4;j++)if(i!==j)M[4*i+j]=get(i,j,X[i],Y[j]);
  for(let d of ds){
   const w=M[d[0]],x=M[4+d[1]],y=M[8+d[2]],z=M[12+d[3]];
   if(w<=a||x<=a||y<=a||z<=a)continue;
   candidates++;
   const A=rows[w],B=rows[x],C=rows[y],D=rows[z];
   let ok=true;
   for(let i=0;i<4;i++){
    let u=A[i], v=B[i], k=C[i], l=D[i];
    if(u===v||u===k||u===l||v===k||v===l||k===l){ok=false;break;}
    let hit=0;
    for(let s of [u,v,k,l])if(s!==X[i]&&s!==Y[i]){
     if(hit++===0)L[i]=s;else R[i]=s;
    }
    if(hit!==2){ok=false;break;}
   }
   if(!ok)continue;
   for(let mask=0;mask<8;mask++){
    const p0=L[0],p1=mask&1?R[1]:L[1],p2=mask&2?R[2]:L[2],p3=mask&4?R[3]:L[3];
    const q0=R[0],q1=mask&1?L[1]:R[1],q2=mask&2?L[2]:R[2],q3=mask&4?L[3]:R[3];
    const s=get(0,1,p0,p1),t=get(0,1,q0,q1);
    if(s<=b||t<=b||rows[s][2]!==p2||rows[s][3]!==p3||rows[t][2]!==q2||rows[t][3]!==q3)continue;
    const old=[a,b,s,t], nw=[w,x,y,z], perm=[];
    for(let i of old){
     const A=rows[i],j=nw.find(x=>rows[x][0]===A[0]),B=rows[j];
     need([1,2,3].every(k=>mod(A[k]+28)===B[k]),"not half turn");
    }
    // Relative to the branch correspondence the coordinate maps are double transpositions.
    for(let k=1;k<4;k++){
      const R0=old.map(i=>old.findIndex(j=>rows[j][0]===rows[nw.find(t=>rows[t][k]===rows[i][k])][0]));
      need(R0.every((v,i)=>v!==i&&R0[v]===i),"not double transposition");
      perm.push(R0);
    }
    total++;for(let q of [...old,...nw])coverage.add(q);
    const branches=old.map(i=>rows[i][0]).sort((a,b)=>a-b);
    const k=branches.join(",");qs.set(k,(qs.get(k)||0)+1);
    examples.push(branches);
   }
  }
 }
 const max=Math.max(...qs.values()),extremal=[...qs].filter(([k,v])=>v===max).map(([k])=>k).sort();
 const five=new Map();
 for(const t of qs.keys())for(let h of slots) {
  const q=t.split(",").map(Number);if(q.includes(h))continue;
  const k=[...q,h].sort((a,b)=>a-b).join(",");five.set(k,(five.get(k)||0)+1);
 }
 return {n,candidates,total,coverage:coverage.size,quartets:qs.size,max,extremal,fiveMax:Math.max(...five.values()),table:qs};
}

const seed=require("node:fs").readFileSync(process.argv[2],"utf8");
const B=boards(seed);
for(const [key,want] of [["full",[8916508,1057,2478,152,14,2]],["marked",[515656,83,484,33,5,2]]]){
 const o=catalogue(B[key]), z=[o.candidates,o.total,o.coverage,o.quartets,o.max,o.fiveMax];
 need(JSON.stringify(z)===JSON.stringify(want),"unexpected catalogue");
 if(key==="full")need(o.extremal.join(";")==="1,2,36,37;1,31,44,54;2,23,38,54","wrong extremal quartets");
 console.log(key,JSON.stringify({...o,table:undefined}));
}
const lanes=[];
for(let h=0;h<=22;h++)for(let s=6;s<=8;s++)if(97-4*h<=2*s){
 const pp=[];
 for(let p=2;p<=28;p+=2)if((54-h-s)*s*p<=(8-s)*(s+7)*(56-p))pp.push(p);
 if(pp.length)lanes.push([h,s,pp]);
}
need(JSON.stringify(lanes)==="[[21,7,[2,4]],[22,6,[2,4,6,8]],[22,7,[2,4]]]","unexpected lanes");
console.log(JSON.stringify(lanes));
