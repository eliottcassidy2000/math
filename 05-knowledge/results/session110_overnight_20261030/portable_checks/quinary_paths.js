function countTower(p) {
  const need=(ok,msg)=>{if(!ok)throw new Error(msg);};
  const maskLimit=(1<<p)-1,half=(p-1)/2,vertices=p*p;
  const rotate=(m,v)=>((m>>>v)|(m<<(p-v)))&maskLimit;
  const canonical=Array.from({length:1<<p},(_,m)=>Math.min(...Array.from({length:p},(_,v)=>rotate(m,v))));
  const reps=[...new Set(canonical)].sort((a,b)=>a-b),base=reps.length;
  const position=new Map(reps.map((x,i)=>[x,i]));
  const encode=a=>a.slice(1).reduce((k,m)=>k*base+position.get(m),a[0]>>>1);
  function decode(k) {const a=Array(p);for(let j=p-1;j>0;j--){a[j]=reps[k%base];k=Math.floor(k/base);}a[0]=2*k+1;return a;}
  const states=2**(p-1)*base**(p-1),levels=Array.from({length:vertices+1},()=>[]),transitions=Array(states),mass=Array(states).fill(0n);
  const pop=m=>m.toString(2).split("1").length-1;
  for(let k=0;k<states;k++) {
    const a=decode(k);levels[a.reduce((s,m)=>s+pop(m),0)].push(k);
    const out=new Map(),add=b=>{const j=encode(b);out.set(j,(out.get(j)||0)+1);};
    for(let v=1;v<=half;v++) if(!(a[0]>>v&1)) {const b=a.slice();b[0]=rotate(a[0]|1<<v,v);add(b);}
    for(let j=1;j<=half;j++) for(let v=0;v<p;v++) if(!(a[j]>>v&1))
      add(Array.from({length:p},(_,i)=>i===0?rotate(a[j]|1<<v,v):canonical[a[(i+j)%p]]));
    transitions[k]=[...out];
  }
  const start=encode([1,...Array(p-1).fill(0)]);
  mass[start]=BigInt(vertices);
  for(let level=1;level<vertices;level++)for(const k of levels[level])for(const [j,w]of transitions[k])mass[j]+=mass[k]*BigInt(w);
  const H=levels[vertices].reduce((s,k)=>s+mass[k],0n);
  const colors=Array(states),intern=new Map(),completion=Array(states).fill(0n),supports=[];
  for(let level=vertices;level>=1;level--) {
    let support=0;
    for(const k of levels[level]) {
      const out=new Map();let value=level===vertices?1n:0n;
      for(const [j,w]of transitions[k]) {out.set(colors[j],(out.get(colors[j])||0)+w);value+=BigInt(w)*completion[j];}
      completion[k]=value;if(value>0n)support++;
      const signature=JSON.stringify([level,[...out].sort((a,b)=>a[0]-b[0])]);
      if(!intern.has(signature))intern.set(signature,intern.size);
      colors[k]=intern.get(signature);
    }
    need(support>0,"Each level has a nonzero completion observation");supports.push(support);
  }
  need(BigInt(vertices)*completion[start]===H,"Forward/backward count");
  function rawNormalize(mask,last) {
    const child=Math.floor(last/p),v=last%p;
    const a=Array.from({length:p},(_,i)=>{let m=0;for(let j=0;j<p;j++)if(mask&(1<<(((child+i)%p)*p+j)))m|=1<<j;return i===0?rotate(m,v):canonical[m];});
    return encode(a);
  }
  function arc(x,y) {
    const a=Math.floor(x/p),b=Math.floor(y/p),d=a===b?(y%p-x%p+p)%p:(b-a+p)%p;
    return d>0&&d<=half;
  }
  function checkRaw(mask,last) {
    const expected=new Map();
    for(let v=0;v<vertices;v++)if(!(mask&(1<<v))&&arc(last,v)){
      const key=rawNormalize(mask|1<<v,v);expected.set(key,(expected.get(key)||0)+1);
    }
    const sort=x=>JSON.stringify([...x].sort((a,b)=>a[0]-b[0]));
    need(sort(expected)===sort(transitions[rawNormalize(mask,last)]),"Raw transition versus orbit quotient");
  }
  let rawChecks=0;
  if(p===3) {
    for(let mask=1;mask<(1<<vertices);mask++)for(let last=0;last<vertices;last++)if(mask&(1<<last)){checkRaw(mask,last);rawChecks++;}
  } else {
    let seed=110;const random=()=>{seed^=seed<<13;seed^=seed>>>17;seed^=seed<<5;return seed>>>0;};
    for(let i=0;i<1000;i++){const last=random()%vertices,mask=(random()&((1<<vertices)-1))|(1<<last);checkRaw(mask,last);rawChecks++;}
  }
  need(rawChecks===(p===3?2304:1000),"Exact raw-validation coverage");
  let rawH=null;
  if(p===3) {
    const dp=Array.from({length:1<<vertices},()=>Array(vertices).fill(0n));
    for(let v=0;v<vertices;v++)dp[1<<v][v]=1n;
    for(let mask=1;mask<(1<<vertices);mask++)for(let v=0;v<vertices;v++)if(dp[mask][v])
      for(let w=0;w<vertices;w++)if(!(mask&(1<<w))&&arc(v,w))dp[mask|1<<w][w]+=dp[mask][v];
    rawH=dp[(1<<vertices)-1].reduce((a,b)=>a+b,0n);
    need(rawH===H&&H===3159n,"Independent raw9 count");
  }
  const reachable=mass.filter(x=>x>0n).length,minimalReachable=new Set(colors.filter((c,k)=>mass[k]>0n)).size;
  if(p===5) {
    need(H===2398512090439546875n&&H%15625n===0n,"Inherited25 count and free automorphism action");
    need(states===65536&&reachable===63488&&intern.size===25437&&minimalReachable===24297,"Quinary carrier census");
  } else need(states===64&&reachable===48&&intern.size===57&&minimalReachable===42,"Ternary carrier census");
  return {p,vertices,H:H.toString(),pointedOrbitStates:states,reachableStates:reachable,
    coarsestFullWeightedLevelPartition:intern.size,coarsestReachableWeightedLevelPartition:minimalReachable,
    rawTransitionControls:rawChecks,independentRawHamiltonCount:rawH===null?null:rawH.toString(),
    acceptanceOnlyLinearRank:vertices,completionSupportSizesFromTop:supports};
}
function verifyQuinaryPaths() {
  return {status:"PASS",cases:[countTower(3),countTower(5)],
    attribution:"The25-vertex Hamilton count is inherited from lattice_bridges_20261021/circ25_H.txt; the orbit quotient and partition counts are the current contribution.",
    boundary:"These64 ternary states are not the earlier tau-pair folded64 carrier. Acceptance-only linear rank does not retain explicit levels on dead states."};
}

console.log(JSON.stringify(verifyQuinaryPaths(), null, 2));
