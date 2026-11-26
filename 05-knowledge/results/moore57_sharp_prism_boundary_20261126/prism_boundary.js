// Node.js: node prism_boundary.js
// Labels 0..5 are the six hexagon parents. Ordinary sheets are 0..53.
const n=54, L=n;
const id=n=>Array.from({length:n},(_,i)=>i);
const comp=(f,g)=>g.map(i=>f[i]);
function inv(f){const r=[]; f.forEach((v,i)=>r[v]=i);return r;}
function assert(p,msg="certificate failure"){if(!p)throw Error(msg);}
function cycles(f){
  const used=new Set(), r=[];
  for(const i of f)if(!used.has(i)){
    let k=i,m=0; do{used.add(k);m++;k=f[k];}while(k!==i);
    r.push(m);
  }
  return r.sort((a,b)=>a-b);
}
function splice(equal,m=n){
  let F=[[1,2,3,0,5,4],[1,0,3,2,5,4],[1,0,3,2,5,4],
         equal?[1,3,4,5,2,0]:[2,4,0,5,1,3]];
  for(let k=6;k<m;k+=4){
    F=F.map((p,j)=>p.concat((j?[3,2,1,0]:[2,0,3,1]).map(x=>k+x)));
    const t=id(k+4);t[0]=k;t[k]=0;
    const j=equal?3:1;F[j]=comp(t,comp(F[j],t));
  }
  return F;
}
const eq=(a,b)=>JSON.stringify(a)===JSON.stringify(b);
function prism(equal){
  const [T,...R]=splice(equal);
  const U=comp(T,comp(R[2],comp(R[1],R[0])));
  const typ=[T,U,...R].map(cycles);
  const rep=(k,x)=>Array(k).fill(x);
  const hyp=equal
    ?[[2,...rep(13,4)],[2,...rep(13,4)],rep(27,2),rep(27,2),[...rep(25,2),4]]
    :[[2,...rep(13,4)],[...rep(12,4),6],rep(27,2),rep(27,2),rep(27,2)];
  assert(eq(typ,hyp),JSON.stringify(typ));
  const top=[0,2,4],bottom=[3,5,1];
  const g=Array.from({length:6*n},()=>new Set());
  const edge=(a,b,f)=>f.forEach((j,i)=>{let x=n*a+i,y=n*b+j; g[x].add(y);g[y].add(x);});
  const I=id(n);
  for(let j=0;j<3;j++){
    const p=j===2?T:I;
    edge(top[j],top[(j+1)%3],p);
    edge(bottom[j],bottom[(j+1)%3],comp(p,R[j]));
    edge(top[j],bottom[j],I);
  }
  const reached=new Set([0]), pending=[0];
  for(let i=0;i<pending.length;i++)for(const w of g[pending[i]])
    if(!reached.has(w)){reached.add(w);pending.push(w);}
  assert(reached.size===6*n);
  const short=g.map((N,v)=>{
    assert(N.size===3);
    const W=new Set([...N,v]);
    for(const u of N)for(const w of g[u])if(w!==v){
      assert(!W.has(w));W.add(w);
    }
    return W;
  });
  const code=(i,x,j,y)=>((6*i+j)*(n+1)+x)*(n+1)+y;
  function allowed(i,x,j,y){
    return x===L||y===L||!short[n*i+x].has(n*j+y);
  }
  const profiles={};
  for(let i=0;i<6;i++)for(let j=i+1;j<6;j++){
    const d=Math.min(j-i,6-j+i);
    const f=Array.from({length:2*n},()=>[]);
    for(let x=0;x<n;x++)for(let y=0;y<n;y++)if(!allowed(i,x,j,y)){
      f[x].push(n+y);f[n+y].push(x);
    }
    assert(f.every(N=>N.length===(d===3?1:2)));
    const seen=new Set(),typ=[];
    for(let k=0;k<2*n;k++)if(!seen.has(k)){
      const Q=[k];seen.add(k);
      for(let t=0;t<Q.length;t++)for(const u of f[Q[t]])if(!seen.has(u)){seen.add(u);Q.push(u);}
      typ.push(Q.length);
    }
    const c={}; for(const x of typ)c[x]=(c[x]||0)+1;
    const k=d+":"+Object.entries(c).map(([l,q])=>l+"^"+q).join(" ");
    profiles[k]=(profiles[k]||0)+1;
  }
  return {allowed,code,profiles};
}
function rng(seed){
  return m=>{seed^=seed<<13;seed^=seed>>>17;seed^=seed<<5;return (seed>>>0)%m;};
}
function shuffle(x,rand){x=x.slice();for(let i=x.length-1;i>0;i--){const j=rand(i+1);[x[i],x[j]]=[x[j],x[i]];}return x;}
// Rows with a distinguished block. Cases can be solved independently.
function boundary({allowed,code},concurrent,seed){
  const rand=rng(seed);
  const special=i=>[i,i+3];
  const all=id(6);
  const root=()=>Array(6).fill(L);
  const prefix=[];
  if(concurrent)prefix.push(root());
  else for(let k=0;k<3;k++){
    const row=root();
    // Ordinary labels on the two unconstrained columns of each pair-only row.
    row[k]=0;row[k+3]=id(n).find(y=>allowed(k,0,k+3,y));
    assert(row[k+3]!==undefined);
    prefix.push(row);
  }
  const seen=new Set();
  function pairs(row,uses){
    for(let i=0;i<6;i++)if(row[i]!==L)for(let j=i+1;j<6;j++)if(row[j]!==L)
      uses(i,row[i],j,row[j]);
  }
  function add(row){pairs(row,(i,x,j,y)=>{const q=code(i,x,j,y);assert(!seen.has(q));seen.add(q);});}
  for(const r of prefix)add(r);
  const output=prefix.slice(),info=[];
  for(let k=0;k<3;k++){
    const cs=all.filter(i=>!special(k).includes(i)),len=n-(concurrent?0:1);
    const pools=cs.map(i=>id(n).filter(x=>concurrent||x!==prefix[i%3][i]));
    assert(pools.every(p=>p.length===len));
    const ok=(i,x,j,y)=>allowed(i,x,j,y)&&!seen.has(i<j?code(i,x,j,y):code(j,y,i,x));
    function score(a){
      let bad=0;
      for(let i=0;i<4;i++)for(let j=i+1;j<4;j++)if(!ok(cs[i],a[i],cs[j],a[j]))bad++;
      return bad;
    }
    let done=false,flips=0;
    for(let attempt=0;attempt<80&&!done;attempt++){
      const cols=pools.map(x=>shuffle(x,rand));
      const rows=id(len).map(i=>cols.map(p=>p[i]));
      let scores=rows.map(score);
      let total=scores.reduce((a,b)=>a+b,0);
      for(let t=0;t<40000&&total;t++){
        const bad=id(len).filter(i=>scores[i]);
        const i=bad[rand(bad.length)],j=rand(len),k=1+rand(3);
        if(i===j)continue;
        const a=rows[i],b=rows[j];
        [a[k],b[k]]=[b[k],a[k]];
        const p=score(a),q=score(b),change=p+q-scores[i]-scores[j];
        if(change<=0||rand(50)===0){total+=change;scores[i]=p;scores[j]=q;flips++;}
        else [a[k],b[k]]=[b[k],a[k]];
      }
      if(total===0){
        for(const a of rows){
          const row=root();cs.forEach((j,i)=>row[j]=a[i]);add(row);output.push(row);
        }
        done=true;
      }
    }
    assert(done,"boundary search did not find a witness");
    info.push(flips);
  }
  const counts=new Map();
  for(const row of output){
    assert(row.length===6&&row.every(a=>Number.isInteger(a)&&a>=0&&a<=L));
    assert(id(3).every(i=>(row[i]===L)===(row[i+3]===L)));
    for(let i=0;i<6;i++)for(let j=i+1;j<6;j++){
      if(i+3===j&&row[i]===L)continue; // It is the same large block.
      const key=code(i,row[i],j,row[j]);
      counts.set(key,1+(counts.get(key)||0));
      assert(counts.get(key)===1&&allowed(i,row[i],j,row[j]),"invalid joint intersection");
    }
  }
  for(let i=0;i<6;i++)for(let j=i+1;j<6;j++)if(j!==i+3)
    for(let x=0;x<=L;x++)for(let y=0;y<=L;y++)if(x===L||y===L)
      assert(counts.get(code(i,x,j,y))===1,"lost large-block pair");
  assert(output.length===(concurrent?163:162));
  assert(id(6).every(i=>output.filter(r=>r[i]===L).length===55));
  const m=concurrent?2754:2755;
  const countsExpected=(i,x)=>!concurrent&&x===prefix[i%3][i]?52:51;
  for(let i=0;i<6;i++)for(let j=i+1;j<6;j++){
    const a=id(n).map(()=>0),b=a.slice();
    for(let x=0;x<n;x++)for(let y=0;y<n;y++)
      if(allowed(i,x,j,y)&&!counts.has(code(i,x,j,y))){a[x]++;b[y]++;}
    assert(a.every((v,x)=>v===countsExpected(i,x)));
    assert(b.every((v,x)=>v===countsExpected(j,x)));
    assert(a.reduce((u,v)=>u+v,0)===m);
  }
  return {output,info};
}
function show(m){return Object.entries(m).map(([k,v])=>k+" x "+v).join("; ");}
for(const e of [true,false]){
  const c=prism(e);
  console.log(e?"Equal triangles":"All squares swap", show(c.profiles));
  for(const t of [true,false]){
    const {output,info}=boundary(c,t,e?(t?401:411):(t?601:611));
    console.log(t?"  Concurrent":"  Separated",output.length,"validated U-tuples",info.join(","));
    if(typeof process!=="undefined"&&process.argv.includes("--tuples"))
      console.log(JSON.stringify(output));
  }
}

function primeForests(){
  const P=id(102).filter(p=>p>=2&&id(Math.floor(Math.sqrt(p))+1).slice(2).every(q=>p%q));
  for(const V of [P.slice(1,-1),P.slice(0,-1),P.slice(1),P]){
    const E=[];
    const square=(i,j)=>Number.isInteger(Math.sqrt(V[i]+V[j]));
    for(let i=0;i<V.length;i++)for(let j=i+1;j<V.length;j++)if(square(i,j))E.push([i,j]);
    const d=id(V.length).map(()=>0), parent=id(V.length), chosen=[];
    const root=i=>{while(parent[i]!==i)i=parent[i];return i;};
    let best=0, sols=[];
    function visit(k){
      if(chosen.length+E.length-k<best)return;
      if(k===E.length){
        if(chosen.length>best){best=chosen.length;sols=[];}
        sols.push(chosen.slice());return;
      }
      const [i,j]=E[k],a=root(i),b=root(j);
      if(d[i]<2&&d[j]<2&&a!==b){
        d[i]++;d[j]++;parent[a]=b;chosen.push(k);
        visit(k+1);
        chosen.pop();parent[a]=a;d[j]--;d[i]--;
      }
      visit(k+1);
    }
    visit(0);
    const names=ks=>ks.map(k=>E[k].map(i=>V[i]).join("-"));
    const fixed=sols[0].filter(k=>sols.every(f=>f.includes(k)));
    const select=sols.map(s=>names(s.filter(k=>!fixed.includes(k))).join(",")).sort();
    assert(eq(select,["3-13,3-61","3-13,3-97","3-61,3-97"]));
    let oriented=0;
    for(const s of sols){
      const labels=new Set(s.flatMap(k=>E[k]));
      const k=V.length-s.length,r=labels.size-s.length;
      oriented+=Array.from({length:k},(_,i)=>1+i).reduce((a,b)=>a*b,1)*2**r;
    }
    const f=names(fixed);
    assert(f.includes("2-7")===V.includes(2)&&f.includes("2-79")===V.includes(2)
           &&f.includes("43-101")===V.includes(101));
    console.log("Prime window",V[0]+".."+V.at(-1),"max square edges",best,
      "forests",sols.length,"orders up to reversal",oriented/2,"fixed",f.join(","));
  }
}
primeForests();
