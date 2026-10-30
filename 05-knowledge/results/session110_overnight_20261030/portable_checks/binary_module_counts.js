function verifyBinaryModuleCounts() {
  const requireTrue = (ok, msg) => { if (!ok) throw new Error(msg); };
  const factorial = [1n];
  for (let n = 1; n <= 12; n++) factorial[n] = factorial[n-1] * BigInt(n);
  const gcd = (a,b) => b ? gcd(b,a%b) : a;
  function choose(n,k) { let x=1n; for(let i=1;i<=k;i++) x=x*BigInt(n-k+i)/BigInt(i); return x; }
  function gaussian(n,k,q) {
    if(k<0 || k>n) return 0n;
    let x=1n; const Q=BigInt(q);
    for(let i=1;i<=k;i++) {
      const num=x*(Q**BigInt(n-k+i)-1n), den=Q**BigInt(i)-1n;
      requireTrue(num%den===0n,"Gaussian integral"); x=num/den;
    }
    return x;
  }
  function partitions(n,lo=1,prefix=[],out=[]) {
    if(n===0) { out.push(prefix); return out; }
    for(let l=lo;l<=n;l+=2) partitions(n-l,l,[...prefix,l],out);
    return out;
  }
  function groups(parts) {
    const out=[];
    for(const l of parts) {
      if(out.length && out[out.length-1][0]===l) out[out.length-1][1]++;
      else out.push([l,1]);
    }
    return out;
  }
  function edgeOrbits(parts) {
    let e=parts.reduce((s,l)=>s+(l-1)/2,0);
    for(let i=0;i<parts.length;i++) for(let j=0;j<i;j++) e+=gcd(parts[i],parts[j]);
    return e;
  }
  const fixedT = parts => 1n << BigInt(edgeOrbits(parts));
  const strongMemo = new Map([["",0n]]);
  function strong(parts) {
    const key=parts.join(",");
    if(strongMemo.has(key)) return strongMemo.get(key);
    let out=fixedT(parts), gs=groups(parts);
    function visit(i,left,right,mult) {
      if(i===gs.length) {
        if(left.length && right.length) out-=mult*strong(left)*fixedT(right);
        return;
      }
      const [l,m]=gs[i];
      for(let k=0;k<=m;k++) visit(i+1,[...left,...Array(k).fill(l)],[...right,...Array(m-k).fill(l)],mult*choose(m,k));
    }
    visit(0,[],[],1n); requireTrue(out>=0n,"Strong fixed count");
    strongMemo.set(key,out); return out;
  }
  function weight(parts,r,k=0,affine=false) {
    if(!parts.length) return 1n;
    const c1=parts.length-(affine?1:0), c3=parts.filter(l=>l%3===0).length, c7=parts.filter(l=>l%7===0).length;
    const other=Array(r+1).fill(0n); other[0]=1n;
    if(r>=2) other[2]=gaussian(c3,1,4);
    if(r>=3) other[3]=2n*gaussian(c7,1,8);
    let out=0n;
    for(let j=0;j<=Math.min(r,c1);j++) for(let v=0;j+v<=r;v++) {
      if(!other[v]) continue;
      let wordCount=0n;
      for(let t=0;t<=Math.min(j,k);t++) if(j+v-t<=r-k) {
        let fullRank=1n;
        for(let i=0;i<t;i++) fullRank*= (1n<<BigInt(k))-(1n<<BigInt(i));
        wordCount+=gaussian(j,t,2)*fullRank;
      }
      out+=gaussian(c1,j,2)*other[v]*wordCount;
    }
    return out;
  }
  const matrices=[];
  for(let a=0;a<8;a++) for(let b=0;b<8;b++) for(let c=0;c<8;c++) {
    const map=Array.from({length:8},(_,x)=>((x&1)?a:0)^((x&2)?b:0)^((x&4)?c:0));
    if(new Set(map).size===8) matrices.push({cols:[a,b,c],map});
  }
  requireTrue(matrices.length===168,"GL3 order");
  const paletteGroups={
    GL3:matrices.map(x=>x.map),
    AGL3:matrices.flatMap(x=>Array.from({length:8},(_,b)=>x.map.map(v=>v^b))),
    PORT1:matrices.filter(x=>(x.cols[0]&1)===1 && !(x.cols[1]&1) && !(x.cols[2]&1)).map(x=>x.map),
    NAMED3:[Array.from({length:8},(_,i)=>i)]
  };
  function cycleLengths(perm) {
    const used=new Set(), out=[];
    for(let x=0;x<perm.length;x++) if(!used.has(x)) {
      let y=x,l=0; do { used.add(y); y=perm[y]; l++; } while(y!==x); out.push(l);
    }
    return out;
  }
  const cycles=Object.fromEntries(Object.entries(paletteGroups).map(([k,v])=>[k,v.map(cycleLengths)]));
  function explicitPalette(parts,name) {
    let out=0n;
    for(const cs of cycles[name]) {
      let n=1n;
      for(const l of parts) n*=BigInt(cs.filter(c=>l%c===0).reduce((a,b)=>a+b,0));
      out+=n;
    }
    const d=BigInt(cycles[name].length);
    requireTrue(out%d===0n,"Twisted color average integral");
    return out/d;
  }
  let comparisons=0;
  const all=[], strongCounts=[];
  const gauges=[["UNCOLORED",0,0,false],["GL3",3,0,false],["AGL3",3,0,true],["PORT1",3,1,false],["NAMED3",3,3,false]];
  for(let n=0;n<=12;n++) {
    const row={n}, sr={n};
    for(const [name,r,k,aff] of gauges) {
      let total=0n, st=0n;
      for(const parts of partitions(n)) {
        const w=weight(parts,r,k,aff);
        if(name!=="UNCOLORED") { requireTrue(w===explicitPalette(parts,name),"Module versus explicit palette group"); comparisons++; }
        let z=1n; for(const [l,m] of groups(parts)) z*=BigInt(l)**BigInt(m)*factorial[m];
        requireTrue(factorial[n]%z===0n,"Cycle class integral");
        const classes=factorial[n]/z;
        total+=classes*fixedT(parts)*w; st+=classes*strong(parts)*w;
      }
      requireTrue(total%factorial[n]===0n && st%factorial[n]===0n,"Burnside integral");
      row[name]=(total/factorial[n]).toString(); sr[name]=(st/factorial[n]).toString();
    }
    all.push(row); strongCounts.push(sr);
  }
  requireTrue(all[4].GL3==="188" && all[4].AGL3==="48" && all[4].PORT1==="756","Inherited controls");
  requireTrue(JSON.stringify(strongCounts.slice(1,6).map(x=>x.GL3))===JSON.stringify(["2","0","8","66","1568"]),"Strong controls");
  requireTrue(all[12].UNCOLORED==="154108311168","A000568 order12 control");
  return {status:"PASS",field:2,rankLimit:3,method:"Semisimple invariant subspaces, independently compared with complete palette groups",
    explicitPaletteOrders:Object.fromEntries(Object.entries(paletteGroups).map(([k,v])=>[k,v.length])),
    termwiseComparisons:comparisons,total:all,strong:strongCounts,
    scope:"Labelled palette action quotients combined with tournament isomorphism; no timing claim"};
}

console.log(JSON.stringify(verifyBinaryModuleCounts(), null, 2));
