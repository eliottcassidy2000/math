"use strict";
const assert=(ok,msg)=>{if(!ok)throw Error(msg||"failed check")};
const bits=x=>x.toString(2).length;
function choose(n,k){
  if(k<0||k>n)return 0n;
  k=Math.min(k,n-k);
  let b=1n;
  for(let j=1;j<=k;j++)b=b*BigInt(n-k+j)/BigInt(j);
  return b;
}
function fact(n){let a=1n;for(let k=2;k<=n;k++)a*=BigInt(k);return a;}
function bounds(){
  const r=1520,s=1730,i0=17024,j=57*r-2*i0;
  const u=r*(r-1)/2,v=r*s;
  assert(j===52592 && u===1154440 && v===2629600);
  assert(32n*BigInt(u-i0)*BigInt(j)*BigInt(j-1)<BigInt(i0+1)*BigInt(v+1-j)*BigInt(v+2-j));
  const first=choose(u,i0)*choose(v,j);
  assert(bits(first)===499859 && 32n*first<31n*(1n<<499859n));
  const labels=bits(choose(r+s,r)-1n);
  assert(labels===3235);
  const lambda=36n;
  let z=0n;
  for(let d=0;d<=57;d++)z+=choose(r-1,d)*lambda**BigInt(d)*choose(s,57-d)**2n;
  const denom=lambda**BigInt(i0);
  const fin=z**760n;
  assert(fin<denom*(1n<<493538n) && fin>denom*(1n<<493537n));
  const cap=(n,shim)=>n*(n-1)/2-(shim?n-1:0);
  for(const [p,n,spine] of [[499859,1001,false],[499859,1002,true],[493538,995,false],[493538,996,true],
         [499859+labels,1004,false],[499859+labels,1005,true],[493538+labels,998,false],[493538+labels,999,true]]){
    assert(cap(n,spine)>=p && cap(n-1,spine)<p);
    if(spine && p===499859)assert(cap(n,true)-p===641);
  }
  for(const [n,p,e] of [[1002,499859,7563],[996,493538,7509]]){
    const m=n/2,a=fact(m),b=fact(m+1),c=fact(n);
    assert(cap(n,true)-p>=Math.floor(n/3)+1);
    assert(b**2n>(1n<<BigInt(15*(m+1))));
    assert(c**5n<(1n<<BigInt(43*n)));
    assert((a*b)**10n*2048n<(1n<<BigInt(10*e)));
  }
  console.log("Exact bits 499859/493538; spine-height upper exponents 7563/7509.");
}
const bit=i=>1n<<BigInt(i);
function control(d){
  const n=d*d+1;
  const adj=Array(n).fill(0n);
  const add=(i,j)=>{adj[i]|=bit(j);adj[j]|=bit(i)};
  if(d===3){
    const pairs=[];
    for(let i=0;i<5;i++)for(let j=i+1;j<5;j++)pairs.push([i,j]);
    pairs.forEach((q,i)=>pairs.forEach((r,j)=>{if(q.every(x=>!r.includes(x)))add(i,j)}));
  }else for(let i=0;i<5;i++)for(let x=0;x<5;x++){
    add(5*i+x,5*i+(x+1)%5); add(25+5*i+x,25+5*i+(x+2)%5);
    for(let j=0;j<5;j++) add(5*i+x,25+5*j+(x+i*j)%5);
  }
  const all=(1n<<BigInt(n))-1n;
  // P at degree 3; I+P at degree 7.
  const mat=adj.map((x,i)=>x^all^(d===3?bit(i):0n));
  const basis=[],rows=new Map();
  for(let i=0;i<n;i++){
    let x=mat[i];
    while(x){
      const j=bits(x)-1;
      if(rows.has(j))x^=rows.get(j);
      else {rows.set(j,x);basis.push(i);break}
    }
  }
  const r=basis.length;
  assert(r===(d===3?4:22));
  const restrict=x=>basis.reduce((y,p,j)=>x&bit(p)?y|bit(j):y,0n);
  const inv=basis.map((p,j)=>restrict(mat[p])|(bit(j)<<BigInt(r)));
  for(let j=0;j<r;j++){
    const i=inv.findIndex((x,i)=>i>=j&&Boolean(x&bit(j)));
    assert(i>=0);
    [inv[i],inv[j]]=[inv[j],inv[i]];
    for(let k=0;k<r;k++)if(k!==j&&inv[k]&bit(j))inv[k]^=inv[j];
  }
  for(let i=0;i<n;i++){
    const y=restrict(mat[i]);let c=0n,recovered=0n;
    for(let j=0;j<r;j++)if(y&bit(j))c^=inv[j]>>BigInt(r);
    for(let j=0;j<r;j++)if(c&bit(j))recovered^=mat[basis[j]];
    assert(recovered===mat[i]);
  }
  console.log("Degree",d,"projector basis rank",r);
}
bounds();control(3);control(7);
