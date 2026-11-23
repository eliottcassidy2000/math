function H(r){let n=r.length,N=1<<n,d=new Uint32Array(n*N);for(let i=0;i<n;i++)d[(1<<i)*n+i]=1;for(let s=1;s<N;s++)for(let i=0;i<n;i++){let h=d[s*n+i];if(h)for(let w=r[i]&~s;w;w&=w-1){let b=w&-w;d[(s|b)*n+31-Math.clz32(b)]+=h;}}let h=0;for(let i=0;i<n;i++)h+=d[(N-1)*n+i];return h;}
function anchor(n,z){let r=Array.from({length:n},(_,i)=>i<n-1?1<<(i+1):0);for(let i=0,k=0;i<n;i++)for(let j=i+2;j<n;j++,k++)if(z>>k&1)r[i]|=1<<j;else r[j]|=1<<i;return r;}
function strong(r){let S=1;for(let j=0;j<r.length;j++)for(let i=0;i<r.length;i++)if(r[i]&S)S|=1<<i;return S===(1<<r.length)-1;}
function iso(r,s){const n=r.length;let pop=(x)=>x===0?0:1+pop(x&(x-1));const R=r.map(pop),S=s.map(pop);const p=[];function rec(i,vis){if(i===n)return true;for(let j=0;j<n;j++)if(!(vis>>j&1)&&R[i]===S[j]&&p.every((q,k)=>((r[i]>>k&1)===(s[j]>>q&1)))){p.push(j);if(rec(i+1,vis|1<<j))return true;p.pop();}return false;}return rec(0,0);}
function join(B){let R=[];for(let r of B){let n=R.length;R=R.map(x=>x|(((1<<r.length)-1)<<n)).concat(r.map(x=>x<<n));}return R;}
function reps(n,h){let V=[];for(let z=0;z<(1<<((n-1)*(n-2)/2));z++){let R=anchor(n,z);if(strong(R)&&H(R)===h&&!V.some(S=>iso(R,S)))V.push(R);}return V;}
function shapes(n,parts){let V=[];for(let B of parts){let spare=n-B.reduce((x,y)=>x+y.length,0);for(let a=0;a<=spare;a++)if(B.length===1)V.push(join([...Array(a).fill([0]),B[0],...Array(spare-a).fill([0])]));else for(let b=0;b<=spare-a;b++)V.push(join([...Array(a).fill([0]),B[0],...Array(b).fill([0]),B[1],...Array(spare-a-b).fill([0])]));}return V;}
function flip(r,i,j){r[i]^=1<<j;r[j]^=1<<i;}
function tests(V,t,two){let count=0,hits=0;for(let r of V){let E=[];for(let i=0;i<r.length;i++)for(let j=i+1;j<r.length;j++)E.push([i,j]);for(let i=0;i<E.length;i++){flip(r,...E[i]);if(two){for(let j=i+1;j<E.length;j++){flip(r,...E[j]);count++;hits+=H(r)===t;flip(r,...E[j]);}}else {count++;hits+=H(r)===t;}flip(r,...E[i]);}}return [V.length,count,hits];}

function Hs(r){let n=r.length,R=(1<<n)-1,S=r.map((x,i)=>x|(1<<i));for(let k=0;k<n;k++)for(let i=0;i<n;i++)if(S[i]>>k&1)S[i]|=S[k];const T=[];while(R){const b=R&-R,i=31-Math.clz32(b),C=[];for(let j=0;j<n;j++)if((R&S[i])>>j&1&&S[j]>>i&1){C.push(j);R^=1<<j;}T.push(C);}let h=1;for(const C of T)if(C.length>1)h*=H(C.map(i=>C.reduce((z,j,k)=>z|((r[i]>>j&1)<<k),0)));return h;}
function replay(log){
 const A=[1022,1020,1016,1008,992,448,640,768,576,32], B=A.slice(),ed=[[1,3],[0,4],[0,1],[1,3],[0,1],[0,4]], hs=[H(B)],core=[H(B.slice(0,5).map(t=>t&31))];
 for(let e of ed){flip(B,...e);hs.push(H(B));core.push(H(B.slice(0,5).map(t=>t&31)));}
 if(String(A)!==String(B)||String(hs)!=='15,45,225,195,75,135,15')throw Error('literal loop');
 log('literal',hs,core);
 const C=[14,28,24,16,1],nb=[];for(let i=0;i<5;i++)for(let j=i+1;j<5;j++){flip(C,i,j);nb.push(H(C));flip(C,i,j);}log('near-transitive H9 neighbors',nb);
 const c15=[...reps(5,15),...reps(6,15)],c5=reps(4,5),tri=[2,4,1];const bank=[...c15.map(x=>[x]),...c5.flatMap(x=>[[x,tri],[tri,x]])];
 const S=shapes(11,bank);log('15 to 55',tests(S,55,false));
 const E=[];for(let i=0;i<11;i++)for(let j=i+1;j<11;j++)E.push([i,j]);
 let stats=[0,0,0,0],off=[0,0,0,0],closures=0;
 const targets=[195,135,15];
 for(const r of S){
  const orig=r.slice();
  function walk(k,prev,valid){
   if(k===targets.length){closures+=iso(orig,r);return;}
   for(let a=0;a<E.length;a++){flip(r,...E[a]);if(Hs(r)===targets[k]){stats[k+1]++; if(valid&&a!==prev)off[k+1]++;walk(k+1,a,valid&&a!==prev);}flip(r,...E[a]);}
  }
  for(let a=0;a<E.length;a++){flip(r,...E[a]);for(let b=a+1;b<E.length;b++){
   flip(r,...E[b]);if(Hs(r)===225){stats[0]++;off[0]++;walk(0,b,true);}flip(r,...E[b]);
  }flip(r,...E[a]);}
 }
 log('unrestricted',stats,'no immediate undos',off,'closures',closures);
}

function min9(log){
 const n=9,c15=[...reps(5,15),...reps(6,15)],c5=reps(4,5),tri=[2,4,1],bank=[...c15.map(x=>[x]),...c5.flatMap(x=>[[x,tri],[tri,x]])],V=shapes(n,bank);
 const E=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)E.push([i,j]);const L=E.length;
 let found=0,stat=[0,0];
 for(const orig of V){
 const cache=new Map(),r=orig.slice();
 function HH(cur){const key=cur.toString();if(cache.has(key))return cache.get(key);for(const a of cur)flip(r,...E[a]);const h=Hs(r);for(const a of cur)flip(r,...E[a]);cache.set(key,h);return h;}
 function tog(A,a){return A.includes(a)?A.filter(x=>x!==a):[...A,a].sort((x,y)=>x-y);}
 function suffix(cur,support,seq,nt){
  if(seq.length===6){if(cur.length===0&&nt===4){found++;if(found===1)log('found',orig,seq.map(x=>E[x]));}return;}
  const left=6-seq.length;
  if(cur.length>left||4-nt>left)return;
  const t=[225,195,135,15];
  for(let a=0;a<L;a++){
   if(!support.includes(a)&&support.length===3)continue;
   const next=tog(cur,a);
   if(next.length>=left)continue;
   const h=HH(next),nv=nt+(h===t[nt]);
   if(4-nv>left-1)continue;
   suffix(next,support.includes(a)?support:[...support,a],[...seq,a],nv);
  }
 }
 for(let a=0;a<L;a++)for(let b=a+1;b<L;b++)if(HH([a,b])===225){stat[0]++;suffix([a,b],[a,b],[a,b],1);}
 const last=[];for(let a=0;a<L;a++)if(HH([a])===135)last.push(a);
 for(const a of last)for(let b=0;b<L;b++)if(b!==a)for(let c=b+1;c<L;c++)if(c!==a){
   const tri=[a,b,c].sort((x,y)=>x-y);
   if(HH(tri)!==225)continue;stat[1]++;
   for(const pair of [[a,b],[a,c]])if(HH(pair.sort((x,y)=>x-y))===195){found++;if(found===1)log('found delayed',orig,tri.map(x=>E[x]),pair.map(x=>E[x]),E[a]);}
 }
 }
 log('labelled at order',n,'sources',V.length,'cases',stat,'loops',found);
}

function paths(r){const n=r.length,N=1<<n,d=Array.from({length:N},()=>Array(n).fill(0));for(let i=0;i<n;i++)d[1<<i][i]=1;for(let S=1;S<N;S++)for(let i=0;i<n;i++)if(d[S][i])for(let m=r[i]&~S;m;m&=m-1){let b=m&-m;d[S|b][31-Math.clz32(b)]+=d[S][i];}return d;}
function shell(log){const A=[82,52,121,115,96,65,2],I=66,F=127,D=paths(A),R=paths(A.map((r,i)=>r^F^(1<<i)));let G=0,L=0;for(let s=1;s<=F;s++){const t=F^s,h=D[s].reduce((u,x,i)=>u+(I>>i&1)*x,0);G+=h*(t?R[t].reduce((a,b)=>a+b,0):1);L+=h*R[t].reduce((u,x,i)=>u+(I>>i&1)*x,0);}if(String([H(A),G,L])!=='49,272,34')throw Error('shell');log('signed shell',H(A),G,L);
 for(let k=5;k<=8;k++){const s=k-4,f=((1<<s)-1)<<8,r=A.map((x,i)=>x|((I>>i&1)<<7));r.push((F^I)|f);for(let i=8;i<8+s;i++)r.push(F|(f&~((1<<(i+1))-1)));if(H(r)!==17*2**k+15)throw Error('shell literal');const d=2**k+1,B=17*d-1,C=d*(B-1);log('lift',k,B,[1,3,15,13,5,9,1].map(a=>a*C));}}

function checked(...args) {
  if(args[0]==="unrestricted" && (String(args[1])!=="348,324,954,0" || args[5]!==0))
    throw Error("five-edit census mismatch");
  if(args[0]==="15 to 55" && String(args[1])!=="56,3080,0")
    throw Error("normalized first jump census mismatch");
  if(args[0]==="labelled at order" &&
    (args[1]!==9 || args[3]!==30 || String(args[5])!=="44,10" || args[7]!==0))
    throw Error("labelled six-edit census mismatch");
  console.log(...args);
}
replay(checked); min9(checked); shell(checked);

function shellStats(r,tag) {
  const n=r.length, F=(1<<n)-1;
  function dp(rows) {
    const d=Array.from({length:F+1},()=>Array(n).fill(0));
    for(let i=0;i<n;i++)d[1<<i][i]=1;
    for(let s=1;s<=F;s++) for(let i=0;i<n;i++) if(d[s][i])
      for(let j=0;j<n;j++)if((rows[i]>>j&1)&&!(s>>j&1))d[s|(1<<j)][j]+=d[s][i];
    return d;
  }
  const end=dp(r), start=dp(r.map((x,i)=>F^(1<<i)^x));
  let D=end[F].reduce((a,b)=>a+b,0), G=0,L=0;
  for(let s=1;s<=F;s++) {
    const a=end[s].reduce((z,v,i)=>z+((tag>>i&1)?v:0),0);
    const b=F^s,c=b?start[b].reduce((z,v)=>z+v,0):1;
    const e=start[b].reduce((z,v,i)=>z+((tag>>i&1)?v:0),0);
    G+=a*c; L+=a*e;
  }
  return {D,G,L};
}

const core4lift=[14,28,24,16,1], tag4lift=24;
const out4lift=shellStats(core4lift,tag4lift);
if([out4lift.D,out4lift.G,out4lift.L].join()!=="9,40,6")
  throw Error("base-4 lifted shell coefficients");
for(let w=1;w<=4;w++){
  const tr=Array.from({length:w},(_,i)=>((1<<w)-1)^((1<<(i+1))-1));
  const r=join([[0],tr,core4lift]);
  for(let i=0;i<core4lift.length;i++)if(tag4lift>>i&1)flip(r,0,1+w+i);
  if(H(r)!==40*2**w+3)throw Error("base-4 lifted shell");
}
console.log("base-4 lifted shell verified");
