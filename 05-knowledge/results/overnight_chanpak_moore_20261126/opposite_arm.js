// Usage: node opposite_arm.js [address cap, default 33554432].
function run(N) {
if(!Number.isInteger(N)||N<2||N>2**26)throw Error("scope");
const f=new Uint32Array(N+1);f[0]=1;
for(let v=1;v<=N;v++) {let q=v>>1; f[v]=f[q]+((v&1)?0:q+(q&1));}
const F=v=>v<=N?f[v]:F(v>>1)+((v&1)?0:(v>>1)+((v>>1)&1));
const H=(b,c)=>b*c+F(b-1)*F(c-1)-(b&1)*(c&1);
const hit=new Map();let readings=0;
function test(b,c,ell,s) {
  if(b<c||b>N)return;
  ++readings;
  const h=H(b,c);
  if(H(1+((b-1)^ell),c)-h!==s)throw Error("bad inverse");
  for(let m=1;m<N;m*=2) if(H(b,1+((c-1)^m))-h===-s) {
    const row=[b,c,ell,m,s,h];
    if(b===c&&m<ell) continue;
    hit.set(String(row),row);
  }
}
for(let c=1;c<=N;c++) {
  const g=f[c-1];
  for(let ell=1;ell*c<=g*(Math.floor(N/(2*ell))+1)+2;ell*=2) {
    const gamma=ell===1?(c&1):0;
    for(let s=-2;s<=2;s+=4) {
      const A=c*ell+gamma-s;
      if(A%g)continue;
      const u=A/g;
      // z loop updates the bit reversal by an MSB-first increment.
      for(let z=0,rev=0;z<ell;z++) {
        const beta=u+z+(ell>1?(z&1):0);
        if(beta>=0&&!(beta&1)){
          const r=(ell-1)^rev;
          for(let q=Math.max(0,beta-1);q<=beta;q++){
            if(q+(q&1)!==beta)continue;
            const b=2*ell*q+r+1;
            test(b,c,ell,s);test(b+ell,c,ell,-s);
          }
        }
        let i=ell>>1;
        while(i&&(rev&i)){rev^=i;i>>=1;}
        rev^=i;
      }
    }
  }
}
const rows=[...hit.values()].sort((x,y)=>x[0]-y[0]||x[1]-y[1]);
return {N,readings,rows};
}

// Materialize both arms and the actual reversed contacts; count independently.
function literal(rows) {
  const rb=(w,k)=>w+1>0&&k>=0 ? Math.max(1,Math.ceil(Math.log2(w+1)),k+1) : 1;
  function graph(b,c,l,r,flip) {
    const i=Math.log2(l),j=Math.log2(r);
    const x=rb(b-1,i), y=x+1+rb(c-1,j), n=y+1;
    const a=new Uint32Array(n);
    const put=(u,v)=>{a[u]|=1<<v;};
    const u=(b-1)^(flip===1?l:0),v=(c-1)^(flip===2?r:0);
    for(let p=0;p<x;p++){
      for(let q=p+1;q<x;q++)put(p,q);
      if((u>>p)&1)put(p,x);else put(x,p);
    }
    for(let p=x+1;p<y;p++){
      for(let q=p+1;q<y;q++)put(q,p);
      if((v>>(p-x-1))&1)put(y,p);else put(p,y);
    }
    for(let p=0;p<=x;p++)for(let q=x+1;q<=y;q++) {
      if(p===x&&q===y)put(p,q);else put(q,p);
    }
    return a;
  }
  function dp(a){
    const n=a.length, full=(1<<n)-1;
    const D=new Float64Array(n*(full+1));
    for(let i=0;i<n;i++)D[(1<<i)*n+i]=1;
    for(let s=1;s<full;s++)for(let i=0;i<n;i++){
      const k=D[s*n+i];if(!k)continue;
      for(let m=a[i]&~s;m;m&=m-1){
        const z=m&-m,j=31-Math.clz32(z);
        D[(s|z)*n+j]+=k;
      }
    }
    let sum=0;for(let i=0;i<n;i++)sum+=D[full*n+i];
    return sum;
  }
  return rows.map(([b,c,l,r,sign,h])=>{
    const a=graph(b,c,l,r,0),z=graph(b,c,l,r,1),w=graph(b,c,l,r,2);
    for(const g of [a,z,w])
      for(let i=0;i<g.length;i++)for(let j=i+1;j<g.length;j++)
       if(((g[i]>>j)&1)+((g[j]>>i)&1)!==1)throw Error("not tournament");
    const Y=[a,z,w].map(dp);if(Y[0]!==h||Y[1]!==h+sign||Y[2]!==h-sign)throw Error("wrong literal counts");
    for(const g of [z,w])if(g.reduce((s,t,i)=>s+(t!==a[i]),0)!==2)throw Error("not one flip");
    return [a.length,...Y];
  });
}

const r=run(process.argv[2]?Number(process.argv[2]):2**25);
console.log(JSON.stringify(r));
console.log(JSON.stringify(literal(r.rows)));
