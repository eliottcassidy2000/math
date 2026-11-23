const text = x => console.log(JSON.stringify(x));

function H(r){let n=r.length,full=(1<<n)-1,f=new Uint32Array((full+1)*n);for(let i=0;i<n;i++)f[(1<<i)*n+i]=1;for(let s=1;s<=full;s++)for(let b=s;b;b&=b-1){let i=31-Math.clz32(b&-b),v=f[s*n+i];if(v)for(let a=r[i]&~s;a;a&=a-1){let j=31-Math.clz32(a&-a);f[(s|(a&-a))*n+j]+=v;}}let h=0;for(let i=0;i<n;i++)h+=f[full*n+i];return h;}
function pairs(n){let e=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)e.push([i,j]);return e;}
function flip(r,[i,j]){r[i]^=1<<j;r[j]^=1<<i;}
function join(C){let n=C.reduce((t,c)=>t+c.length,0),r=[],s=0;for(let c of C){for(let x of c)r.push((x<<s)|((1<<n)-(1<<(s+c.length))));s+=c.length;}return r;}
function canonical(r){let n=r.length,m=Infinity;function f(p,s){if(p.length===n){let mask=0,k=0;for(let i=0;i<n;i++)for(let j=i+1;j<n;j++,k++)if(r[p[i]]>>p[j]&1)mask|=1<<k;m=Math.min(m,mask);return;}for(let z=s;z;z&=z-1){p.push(31-Math.clz32(z&-z));f(p,s^(z&-z));p.pop();}}f([],(1<<n)-1);return m;}
let coreBank=new Map(),free=pairs(5).filter(([i,j])=>j>i+1);for(let mask=0;mask<1<<free.length;mask++){let r=Array(5).fill(0);for(let i=0;i<4;i++)r[i]|=1<<(i+1);for(let t=0;t<free.length;t++){let [i,j]=free[t];if(mask>>t&1)r[i]|=1<<j;else r[j]|=1<<i;}if(H(r)===11)coreBank.set(canonical(r),r);}
if(coreBank.size!==1)throw Error("wrong core bank");
let core=[...coreBank.values()][0];text({core});
let schedules=new Map();for(let b of [[2,1,1,1,2,2],[2,1,1,2,1,2]])for(let i=0;i<b.length;i++){let s=[...b];s[i]++;schedules.set(s.join(""),s);}
let checkpoints=[121,99,55,33,77,11];
function test(n,exact){
 let starts=[];for(let i=0;i<=n-5;i++)starts.push(join([...Array(i).fill([0]),core,...Array(n-5-i).fill([0])]));
 let es=pairs(n),receipts={};
 for(let [key,s] of schedules){
   let stats=Array(6).fill(0);receipts[key]=stats;
   for(let start of starts){
    let r=[...start],odd=new Uint8Array(es.length);
    function phase(p,remaining,parity){
     if(p===6)return;
     function select(k,lo,left,sz){
      if(exact&&sz>left)return;
      if(!k){
       if(H(r)!==checkpoints[p])return;
       stats[p]++;phase(p+1,left,sz);return;
      }
      for(let i=lo;i<=es.length-k;i++){
       flip(r,es[i]);odd[i]^=1;
       select(k-1,i+1,left-1,sz+(odd[i]?1:-1));
       odd[i]^=1;flip(r,es[i]);
      }
     }
     select(s[p],0,remaining,parity);
    }
    phase(0,10,0);
   }
 }
 return {n,starts:starts.length,exact,receipts};
}
function witness(start,edits,map){let r=[...start],path=[H(r)];for(let e of edits){flip(r,e);path.push(H(r));}for(let [i,j] of pairs(r.length))if((start[i]>>j&1)!==(r[map[i]]>>map[j]&1))throw Error("invalid endpoint map");return {r,path,identical:r.every((x,i)=>x===start[i]),map};}
text(witness([1022,972,984,976,994,974,896,768,512,0],[[3,7],[3,9],[3,6],[1,3],[3,9],[1,3],[0,3],[3,7],[0,3],[3,6]],[0,1,2,3,4,5,6,7,8,9]));
text(witness([486,492,488,497,487,448,384,256,0],[[2,6],[2,8],[2,5],[0,2],[2,8],[0,2],[1,4],[2,7],[1,5],[2,3]],[4,0,5,1,3,6,7,2,8]));

const small=test(8,false),fixed=test(9,true);
if(small.starts!==4||Object.keys(small.receipts).length!==11||
   Object.values(small.receipts).some(a=>a[2]!==0))
  throw Error("order-eight numerical obstruction");
if(fixed.starts!==5||Object.keys(fixed.receipts).length!==11||
   Object.values(fixed.receipts).some(a=>a[5]!==0)||
   String(fixed.receipts["212122"].slice(2))!=="6,4,0,0"||
   String(fixed.receipts["212212"].slice(2))!=="6,16,0,0")
  throw Error("order-nine labelled obstruction");
text(small);text(fixed);


function shellCore(A,tag,w){
 const trans=Array.from({length:w},(_,i)=>((1<<w)-1)^((1<<(i+1))-1));
 const r=join([[0],trans,A]);
 for(let i=0;i<A.length;i++)if(tag>>i&1)flip(r,[0,w+1+i]);
 return r;
}
for(let k=5;k<=7;k++){
 const x=2**k-2,d=11*x+1,e=13*x+1,B=143*x+12;
 if(H(shellCore([14,28,24,16,33,15],28,k-4))!==d
   ||H(shellCore([6,4,8,19,39,15],37,k-3))!==e
   ||d*(B-1)!==11*d*e)throw Error("base-12 shell lifts");
 text({k,B,d,e});
}