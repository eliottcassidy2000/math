"use strict";
const {cores,flip,dp,F,reverse,seam,expect}=require("./commuting_unit_forks");

function* roots(maxN){
 const strong=[null,[[0]]];yield [1,[0]];
 for(let n=2;n<=maxN;n++){
  const own=[], root=n-1, rt=1<<root;
  for(let c=1;c<n;c++)for(const core of strong[c]){
   const q=n-1-c, C=(rt-1)^((1<<q)-1);
   const pref=Array.from({length:q},(_,i)=>(rt-1)^((1<<(i+1))-1));
   const y=core.map((a,i)=>(a<<q)|(i===c-1?rt:0));
   for(let m=0;m<(1<<q);m++){
    const a=pref.map((v,i)=>v|(m>>i&1?0:rt)).concat(y);
    a.push((C^(1<<(q+c-1)))|m);
    if(q?(m&1):c>1)own.push(a);
    yield [n,a];
   }
  }
  strong[n]=own;
 }
}
function inside(A){
 const n=A.length, r=n-1, [P,D]=dp(A), [A0]=dp(A,r), M=P.length-1;
 let V=0,Q=0, R=0;
 for(let s=0;s<=M;s++){
  Q+=P[s]*P[M^s];
  if(s>>r&1){V+=A0[s]*P[M^s];R+=D[s*n+r]*P[M^s];}
 }
 expect(D[M*n+r]===1);
 return [P[M],A0[M],R,V,Q].map(BigInt);
}
function sqrt(x){
 if(x<2n)return x;let t=1n<<BigInt(Math.ceil(x.toString(2).length/2));
 for(;;){const y=(t+x/t)>>1n;if(y>=t)return t;t=y;}
}
function inverseDouble(E,H,R){
 const b=sqrt(E*H/2n), lo=b>12n*H?b-12n*H:1n, hi=b+12n*H+1n;
 let ell=1n,i=0, k=hi/R, reads=0;while(2n*ell<=k){ell*=2n;i++;}
 const out=[];
 for(;i>=1;i--,ell>>=1n){
  const u=R*ell;if(u+H*(ell-1n)<lo)break;
  const z0=u>=lo?0n:(lo-u+H-1n)/H, zhi=(hi-u)/H, z1=zhi<ell-1n?zhi:ell-1n;
  for(const s of [1n,-1n]){
   if((u-2n*s)%H)continue;
   for(let z=z0;z<=z1;z++){
    const r=(ell-1n)^reverse(z,i);if(r&1n)continue;
    const b=(u-2n*s)/H+z+(z&1n);expect(!(b&1n));
    for(const q of [b,b-1n]){
     if(q<=0n)continue;
     reads++;const c=2n*ell*q+r+1n;
     if(R*c+H*F(c-1n)+s===E){
      expect(R*(c+ell)+H*F(c+ell-1n)===E+s);
      out.push([String(c),i,String(s),String(r)]);
     }
    }
   }
  }
 }
 return {out,reads};
}
function nested(maxN,exponents){
 const ex=exponents.map(p=>[p,(1n<<BigInt(p-1))*((1n<<BigInt(p))-1n)]);
 const totals=Object.fromEntries(ex.map(([p])=>[p,{reads:0,lanesP:0,splitsP:0,readsP:0,witnesses:[]}]));
 const classCounts={}, distinct=new Set(), distinctPrefix=new Set();let candidates=0;
 for(const [n,A] of roots(maxN)){
  classCounts[n]=(classCounts[n]||0)+1;
  const [h,e,R,V,Q]=inside(A);
  for(let j=0;j<cores.length;j++){
   const [, , ,H0,R0,El]=cores[j],H=BigInt(H0),L=BigInt(R0);
   const a=H*h,b=H*R+L*Q,d=H*e+L*V+BigInt(El);
   if(!distinct.has([a,b,d].join())){
   distinct.add([a,b,d].join());candidates++;
   expect(a<=b);
   for(const [p,T] of ex)for(const k of [0n,2n]){
    const r=inverseDouble(T+d+k,a,b);
    totals[p].reads+=r.reads;
    for(const w of r.out)totals[p].witnesses.push({n,A,j,profile:[a,b,d+k].map(String),w});
   }
   }
   const dp=L*V+BigInt(El), kP=[a,b,dp].join();
   if(!distinctPrefix.has(kP)){
    distinctPrefix.add(kP);
    for(const [p,T] of ex)for(const k of [0n,2n]){
     const r=inversePrefix(T+dp+k,a,b);
     totals[p].readsP+=r.reads;totals[p].lanesP+=r.lanes;totals[p].splitsP+=r.splits;
     for(const w of r.out)totals[p].witnesses.push({n,A,j,profile:[a,b,dp+k].map(String),w});
    }
   }
  }
 }
 return {classCounts,candidates,prefixCandidates:distinctPrefix.size,totals};
}

function join(A,u,B,v){
 const n=A.length,r=n+v,T=A.concat(B.map(x=>x<<n));
 for(let y=n;y<T.length;y++)for(let x=0;x<n;x++){
  if(x===u&&y===r)T[x]|=1<<y;else T[y]|=1<<x;
 }return T;
}
function nestedControls(){
 const small=[];
 for(const [n,C] of roots(6)){
  const [h,e,R,V,Q]=inside(C), [A,u,[x,y],h0,r0,e0]=cores[0], H=BigInt(h0), L=BigInt(r0);
  const a=H*h, b=H*R+L*Q,d=H*e+L*V+BigInt(e0);
  for(let i=1;i<=40;i++)for(const s of [1n,-1n]){
   const ell=1n<<BigInt(i);
   if((b*ell-2n*s)%a)continue;
   const t=(b*ell-2n*s)/a+(ell-1n)+1n;
   const c=2n*ell*t+1n;
   const M=b*c+a*F(c-1n)+s;
   const rr=inverseDouble(M,a,b);
   expect(rr.out.some(w=>w[0]===String(c)&&w[1]===i&&w[2]===String(s)));
   if(c<32n){
    const W=seam(C,n-1,Number(c)), Z=join(A,u,W,W.length-1);
    const [P]=dp(Z);expect(BigInt(P[P.length-1])===M-s-d);
    const ZZ=flip(Z,A.length+n+i,Z.length-1);
    const [PP]=dp(ZZ);expect(BigInt(PP[PP.length-1])===M+s-d);
   }
   if(small.length<16)small.push([n,i,String(s),String(c),String(M-d)]);
  }
 }
 return small;
}

function inversePrefix(E,H,B){
 const b=sqrt(E*H/2n),lo=b>12n*H?b-12n*H:1n,hi=b+12n*H+1n;
 let u=B, t=0, reads=0, splits=0, lanes=0;
 while(2n*u<=hi){u*=2n;t++;}
 const out=[];
 for(;t>=2&&3n*u>=2n*lo;t--,u>>=1n){
  for(const s of [1n,-1n]){
   if((u-2n*s)%H)continue;lanes++;
   for(let i=1,ell=2n;i<t;i++,ell<<=1n){
    if(u+H*(ell-1n)<lo)continue;
    splits++;
    const R=u/ell,z0=u>=lo?0n:(lo-u+H-1n)/H,zhi=(hi-u)/H,z1=zhi<ell-1n?zhi:ell-1n;
    for(let z=z0;z<=z1;z++){
     const r=(ell-1n)^reverse(z,i);if(r&1n)continue;
     const beta=(u-2n*s)/H+z+(z&1n);expect(!(beta&1n));
     for(const q of [beta,beta-1n]){
      if(q<=0n)continue;reads++;
      const c=2n*ell*q+r+1n;
      if(R*c+H*F(c-1n)+s===E){
       expect(R*(c+ell)+H*F(c+ell-1n)===E+s);
       out.push([String(c),i,String(s),String(r),t-i]);
      }
     }
    }
   }
  }
 }
 return {out,reads,splits,lanes};
}
function literalNested(){
 let tested=0;
 const [A,u,,h0,r0,e0]=cores[0], H=BigInt(h0),R=BigInt(r0);
 for(const [n,C] of roots(5))for(const c of [1,2,3,5,9]){
  const [h,e,x,V,Q]=inside(C), a=H*h,b=H*x+R*Q,d=H*e*BigInt(c%2)+R*V+BigInt(e0);
  const B=seam(C,n-1,c),Z=join(A,u,B,B.length-1), [P]=dp(Z);
  expect(BigInt(P[P.length-1])===a*F(BigInt(c-1))+b*BigInt(c)-d);
  tested++;
 }
 for(const k of [1,2,3,5,11]){
  const a=23n,B=71n;
  for(let i=1;i<=24;i++){
   const ell=1n<<BigInt(i),R=B<<BigInt(k),u=R*ell;
   if((u-2n)%a)continue;
   const q=(u-2n)/a+ell;
   const c=2n*ell*q+1n;
   const rr=inversePrefix(R*c+a*F(c-1n)+1n,a,B).out;
   expect(rr.some(w=>w[0]===String(c)&&w[1]===i&&w[2]==="1"&&w[4]===k));
  }
 }
 return tested;
}

function eligiblePerfect(maxP){
 const out=[], tests={eligible:0,reads:0};
 for(let p=5;p<=maxP;p+=22){
  const T=(1n<<BigInt(p-1))*((1n<<BigInt(p))-1n);
  tests.eligible++;
  for(const e of [8n,10n]){
   const rr=inverseDouble(T+e,23n,24n);
   tests.reads+=rr.reads;
   if(rr.out.length)out.push([p,String(e),rr.out]);
  }
 }
 return {maxP,tests,out};
}

function explicit(){
 const C=[126,60,56,48,96,64,14], [h,e,R,V,Q]=inside(C);
 const a=25n*h,b=25n*R+28n*Q,d=25n*e+28n*V+8n;
 const c=253n,N=a*F(c-1n)+b*c-d,W=a*F(c+1n)+b*(c+2n)-d;
 const counts=[N,N-2n,W,W-2n].map(String);
 expect(counts.join()==="4618983,4618981,4618985,4618983");
 return {core:C,profile:[h,e,R,V,Q].map(String),counts};
}
if(require.main!==module){module.exports={inside,inverseDouble,inversePrefix};}else{
const n=Number(process.argv[2]||10);
expect(n>=1&&n<=13);
const out={literal:literalNested(),example:explicit(),inverseControls:nestedControls(),search:nested(n,[19,31,61,89,107,127])};
expect(out.literal===110);
if(n===12){
 expect(out.search.candidates===53103);expect(out.search.prefixCandidates===53097);
 expect(Object.values(out.search.totals).every(x=>!x.witnesses.length));
}
if(process.argv.includes("--perfect")){
 out.perfect=eligiblePerfect(5000);
 expect(out.perfect.tests.eligible===228&&out.perfect.tests.reads===21860&&!out.perfect.out.length);
 for(const p of [1,3,5,7,9,11]){
  const T=(1n<<BigInt(p-1))*((1n<<BigInt(p))-1n);
  for(const row of cores)for(const d of [0,2])
   expect(!inverseDouble(T+BigInt(row[5]+d),BigInt(row[3]),BigInt(row[4])).out.length);
 }
}
console.log(JSON.stringify(out,null,2));

}
