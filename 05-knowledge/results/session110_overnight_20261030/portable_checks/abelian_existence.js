function verifyAbelianExistenceControls() {
 const need=(ok,msg)=>{if(!ok)throw new Error(msg);},results=[];
 const specs=[[9],[15],[21],[25],[5,5],[27],[9,3],[3,3,3],[33],[35],[3,3]];
 const expectedBounds=[4,16,32,16,96,64,256,832,128,64,16];
 for(let caseIndex=0;caseIndex<specs.length;caseIndex++){
  const mod=specs[caseIndex],N=mod.reduce((a,b)=>a*b,1);
  const coords=Array.from({length:N},(_,n)=>mod.map(m=>{const r=n%m;n=Math.floor(n/m);return r;}));
  const index=c=>{let n=0,b=1;c.forEach((v,i)=>{n+=v*b;b*=mod[i];});return n;};
  const add=(a,b)=>index(coords[a].map((v,i)=>(v+coords[b][i])%mod[i]));
  const subgroups=new Map([["0",[0]]]),queue=[[0]];
  for(let at=0;at<queue.length;at++){
   const H=queue[at];
   for(let g=1;g<N;g++){
    if(H.includes(g))continue;
    const S=new Set(H);let x=g;
    while(!H.includes(x)){for(const h of H)S.add(add(h,x));x=add(x,g);}
    const L=[...S].sort((a,b)=>a-b),key=L.join(",");
    if(!subgroups.has(key)){subgroups.set(key,L);queue.push(L);}
   }
  }
  const counts={};let bound=0n;
  for(const H of subgroups.values())if(H.length>1&&H.length<N){
   counts[H.length]=(counts[H.length]||0)+1;
   bound+=1n<<BigInt((H.length+N/H.length-2)/2);
  }
  const total=1n<<BigInt((N-1)/2);
  need(bound===BigInt(expectedBounds[caseIndex]),"Subgroup union-bound control");
  let exceptionMasks=null;
  if(mod.length===2&&mod[0]===3&&mod[1]===3){
   const lines=[...subgroups.values()].filter(H=>H.length===3),pairs=[],seen=new Set([0]);
   for(let x=1;x<N;x++)if(!seen.has(x)){const y=index(coords[x].map(v=>(3-v)%3));pairs.push([x,y]);seen.add(x);seen.add(y);}
   let checked=0;
   for(let mask=0;mask<16;mask++){
    const S=new Set(pairs.map((pair,i)=>pair[(mask>>i)&1]));
    const cuts=lines.filter(H=>coords.every((_,x)=>H.includes(x)||H.every(h=>S.has(x)===S.has(add(x,h)))));
    need(cuts.length===1,"Exactly one line module in every exceptional state");checked++;
   }
   exceptionMasks=checked;need(bound===total,"Exceptional coverage");
  }else need(bound<total,"Small nonexceptional existence");
  results.push({cyclicFactors:mod,subgroupCounts:counts,compositeUnionBound:bound.toString(),allStates:total.toString(),exceptionMasks});
 }
 need(1369n**3n<2n**34n,"Exact N37 boundary");
 return {status:"PASS",cases:results,boundaryCertificate:{left:(1369n**3n).toString(),right:(2n**34n).toString()},
  scope:"Exact subgroup-count union bounds, not exact prime counts except the separately checked exceptional coverage"};
}

console.log(JSON.stringify(verifyAbelianExistenceControls(), null, 2));
