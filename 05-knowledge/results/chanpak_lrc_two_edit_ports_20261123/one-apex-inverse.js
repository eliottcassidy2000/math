"use strict";
// O is bit 1, contact index 0 is the least significant bit. The apex
// is not included in the word length. No floating-point admission tests.
const pow = n => 1n << BigInt(n);
const min = (a,b) => a < b ? a : b;
const max = (a,b) => a > b ? a : b;
const bit = (w,j) => (w >> BigInt(j)) & 1n;
function height(w,n) {
  let H=1n, sum=0n, prev=0n;
  for(let j=0;j<n;j++){
    const b=bit(w,j);
    if(!b) H+=sum+prev;
    sum=2n*sum+b;
    prev=b;
  }
  return H;
}
function floor(a,b){const q=a/b; return a%b<0n?q-1n:q;}
const ceil = (a,b) => -floor(-a,b);
function sqrt(a) {
  if(a<0n) throw Error("negative square root");
  if(a<2n) return a;
  let t=pow(Math.ceil(a.toString(2).length/2));
  for(;;){const u=(t+a/t)/2n;if(u>=t)return t;t=u;}
}
function band(a,b,c,e) {
  const discr=b*b-4n*a*c, hi=discr+4n*a*e, lo=discr-4n*a*e;
  if(hi<0n)return [];
  const u=sqrt(hi);
  let intervals=[[-u,u]];
  if(lo>0n) {
    const v0=sqrt(lo), v=v0*v0<lo?v0+1n:v0;
    intervals=v>u ? [] : [[-u,-v],[v,u]];
  }
  return intervals.map(([x,y])=>[ceil(x-b,2n*a),floor(y-b,2n*a)]);
}
function reverse(w,n) {
  let v=0n;
  for(let i=0;i<n;i++)v=2n*v+bit(w,i);
  return v;
}
const key = z => [z.n,z.w,z.i,z.j].join(":");
function inverse(p) {
  p=BigInt(p);
  if(p<5n||!(p&1n))throw Error("use odd starting heights >= 5");
  const ans=new Map();
  const add = (w,n,i,j) => {
    if(height(w,n)!==p || height(w^pow(i)^pow(j),n)!==7n*p)return;
    const a={w,n,i,j};
    const b={w:reverse(w,n)^(pow(n)-1n),n,i:n-1-j,j:n-1-i};
    ans.set(key(a),a);ans.set(key(b),b);
  };
  const len=p.toString(2).length;
  for(let m=len;m<=len+1;m++) {
    const T=pow(m), E=BigInt(m+1), bL=pow(Math.ceil(m/2)), bR=pow(Math.floor(m/2));
    const eta = s => s===1 ? 0n : pow(s-1);
    function test(a,b,c,e,emit,deriv=null) {
      const seen=new Set();
      for(let [l,u] of band(a,b,c,e)) {
        l=max(l,T/2n);u=min(u,T-1n);
        if(deriv!==null)l=max(l,ceil(deriv-b,2n*a));
        if(l>u)continue;
        for(let lb=l/bL;lb<=u/bL;lb++){
          const lmin=max(l,lb*bL), lmax=min(u,(lb+1n)*bL-1n);
          // Avoid negative lower numerator before using monotonicity.
          const rmin=max(T/2n,ceil(max(0n,p-E)*T,lmax));
          const rmax=min(T-1n,floor((p+E)*T,lmin));
          if(rmin>rmax)continue;
          for(let rb=rmin/bR;rb<=rmax/bR;rb++){
            const word=((T-1n)^(lb*bL)) & (T-bL) | reverse(rb,Math.ceil(m/2));
            if((word&1n)===0n||bit(word,m-1)||seen.has(word))continue;
            seen.add(word);emit(word);
          }
        }
      }
    }
    // Two same-side extensions: exact binary inverse.
    if((6n*p-2n)%5n===0n) {
      const L=(6n*p-2n)/5n;
      if(L>=T/2n&&L<T&&!(L&1n))add(((T-1n)^L)<<3n,m+3,0,2);
    }
    // Opposite-side exterior marks.
    for(let s=1;s<=3;s++)for(let t=1;t<=4-s;t++){
      const ks=pow(s-1),kt=pow(t-1),C=pow(s+t-2)*T+eta(s)+eta(t);
      test(ks,C-6n*p,kt*p*T,kt*E*T,word=>{
        const full=(word<<BigInt(s))|((pow(t)-1n)<<BigInt(m+s));
        add(full,m+s+t,0,m+s+t-1);
      });
    }
    // Mixed: exterior on the left, then include each converse.
    for(let s=1;s<=6;s++)for(let j=0;j<m;j++)for(let sign of [1n,-1n]){
      const l=sign*pow(j), r=-sign*pow(m-1-j), A=pow(s-1)*T+r;
      const q=s===1 ? (j===0&&sign===1n?1n:0n) : eta(s);
      test(A,A*l+(q-6n*p)*T,l*p*T,2n*E*T*T,(word)=>{
        if(bit(word,j)!==(sign===1n?1n:0n))return;
        add(word<<BigInt(s),m+s,0,s+j);
      },m>=9 ? T*T/16n : null);
    }
  }
  return [...ans.values()];
}
function replay() {
  let first={}, tests=0;
  const hits=new Map();
  const cap=4095n;
  // Independent enumeration over active words and reduced pad placements.
  // Both marks are contacts; a reduced host has no unused exterior pad.
  function scan(w,n,i,j,p,kind) {
    if(height(w^pow(i)^pow(j),n)!==7n*p)return;
    const a={w,n,i,j}, b={w:reverse(w,n)^(pow(n)-1n),n,i:n-1-j,j:n-1-i};
    if(!hits.has(p))hits.set(p,new Set());
    hits.get(p).add(key(a));hits.get(p).add(key(b));
    let prime=p!==7n;
    for(let d=3n;d*d<=p;d+=2n)if(p%d===0n){prime=false;break;}
    if(prime&&(!first[kind]||n<first[kind].n))first[kind]={w,n,i,j,p};
  }
  for(let m=2;m<=cap.toString(2).length+1;m++)for(let w=1n;w<pow(m-1);w+=2n){
    const p=height(w,m);if(p<5n||p>cap)continue;
    for(let s=2;s<=4;s++)for(let t=1;t<s;t++)
      scan(w<<BigInt(s),m+s,0,s-t,p,"same");
    for(let s=1;s<=3;s++)for(let t=1;t<=4-s;t++)
      scan((w<<BigInt(s))|((pow(t)-1n)<<BigInt(m+s)),m+s+t,0,m+s+t-1,p,"opposite");
    for(let s=1;s<=6;s++)for(let j=0;j<m;j++)
      scan(w<<BigInt(s),m+s,0,s+j,p,"mixed");
  }
  let count=0, byKind=[0,0,0];
  for(let p=5n;p<=cap;p+=2n) {
    const got=inverse(p), want=hits.get(p)??new Set();
    if(got.length!==want.size||got.some(x=>!want.has(key(x))))throw Error("p="+p);
    tests++;count+=got.length;
    for(const {w,n,i,j} of got) {
      let lo=0,hi=n-1;
      while(!bit(w,lo))lo++;
      while(bit(w,hi))hi--;
      const outs=[i,j].filter(k=>k<lo||k>hi);
      byKind[outs.length===1 ? 2 : (outs[0]<lo)!==(outs[1]<lo) ? 1 : 0]++;
    }
  }
  for(const [k,p,n] of [["same",37n,9],["opposite",97n,10],["mixed",103n,11]]) {
    if(first[k].p!==p||first[k].n!==n)throw Error(k);
  }
  const packet = (digits,i,j) => {
    const w=[...digits].reduce((z,b,k)=>z+(b==="1"?pow(k):0n),0n), n=digits.length;
    return [w,w^pow(i),w^pow(j),w^pow(i)^pow(j)].map(v=>height(v,n)).join(",");
  };
  if(packet("0110100011",0,9)!=="97,213,307,679")throw Error("97 square");
  if(packet("00010100110",0,5)!=="103,723,85,721")throw Error("103 square");
  for(const p of [103n,227n,1399n]) {
    const a=inverse(p);
    if(a.length!==2||a.some(({w,n,i,j})=>[i,j].some(k=>{
      const h=height(w^pow(k),n);return p<=h&&h<=7n*p;
    }))) throw Error("unique nonmonotone packet " + p);
  }
  if(byKind.join(",")!=="6,2,32")throw Error("packet partition");
  return {tests,count,byKind,first};
}
if(typeof module!=="undefined"){
  module.exports={inverse,height,replay};
  if(require.main===module) console.log(replay());
}
