
const {cores,dp,F}=require("./commuting_unit_forks");
const {inside,inverseDouble,inversePrefix}=require("./nested_fork_inverse");
const assert=require("node:assert/strict");
const MAX=22, target=(1n<<60n)*((1n<<61n)-1n);
const log2=x=>x<4294967296?31-Math.clz32(x):63-Math.clz32(Math.floor(x/4294967296));
function gcd(x,y){while(y){const z=x%y;x=y;y=z;}return x;}
function read(n,v,h,r,Q){
 if(n>8)return;
 const A=new Array(n+1).fill(0);
 for(let j=0;j<n;j++){
  A[j]|=(1<<j)-1;
  if(v>>j&1)A[n]|=1<<j;else A[j]|=1<<n;
 }
 const [a,b,c,d,e]=inside(A);
 assert.deepEqual([a,b,c,d,e],[h,1,r,v+1,Q].map(BigInt));
 assert.equal(F(BigInt(v)),BigInt(h));
}
function scan(max, screen=true){
 const totals={v:0,fullControls:0,passes:0,checks:0,lanesP:0,splitsP:0,readsP:0,reads:0,hits:[]};
 const p2=Array.from({length:max+8},(_,i)=>2**i);
 function visit(j,v,h,r,Q,u1,u3,w,prev){
  if(j===max)return;
  const t=p2[j];
  for(const x of [0,1]){
   const vv=v+x*t, hh=h+x*(w+(prev===0?1:0)),
    rr=2*r+(1-x), QQ=2*Q+x*2*(u1+u3);
   if(x){
    totals.v++;
    if(j<8){read(j+1,vv,hh,rr+1,QQ);totals.fullControls++;}
    if(screen)for(let g=0;g<cores.length;g++){
     const [, , ,H,R,E]=cores[g],a=H*hh,b=H*(rr+1)+R*QQ;
     assert(Number.isSafeInteger(a)&&Number.isSafeInteger(b));
     if(gcd(a,b)!==1)continue;
     const k=log2(a),m=log2(b);
     const jLo=Math.max(1,Math.floor((k+119)/2-3-m)+1), jHi=Math.ceil((k+121)/2+1-m)-1;
     const A=BigInt(a),B=BigInt(b);
     let test=false;
     for(let J=jLo;J<=jHi;J++){
      const u=B<<BigInt(J),q=u%A;if(q!==2n&&q!==A-2n)continue;
      totals.passes++;test=true;
     }
     if(!test)continue;
     totals.checks++;
     for(const shift of [0,2]){
      const D=BigInt(H+R*(vv+1)+E+shift),D0=BigInt(R*(vv+1)+E+shift);
      const i=inverseDouble(target+D,A,B),p=inversePrefix(target+D0,A,B);
      totals.reads+=i.reads;totals.lanesP+=p.lanes;totals.splitsP+=p.splits;totals.readsP+=p.reads;
      for(const y of [...i.out,...p.out])totals.hits.push([vv,g,shift,y]);
     }
    }
   }
   visit(j+1,vv,hh,rr,QQ,u1+(1-x)*t,3*u3+(1-x)*t,2*w+1-x,x);
  }
 }
 read(0,0,1,1,2);
 // The low bit is zero. Initialise j=1 after appending that bit.
 visit(1,0,1,1,4,1,1,1,0);
 return totals;
}
const r=scan(Number(process.argv[2]||MAX));
assert.equal(r.fullControls,127);
if(!process.argv[2]){
 assert.equal(r.v,2097151);assert(!r.hits.length);assert.equal(r.reads+r.readsP,0);
 assert.deepEqual([r.passes,r.lanesP,r.splitsP],[15,4,0]);
}
console.log(JSON.stringify(r,null,2));
