function verifyRamifiedPeriods() {
  const need=(ok,msg)=>{if(!ok)throw new Error(msg);};
  const reports=[];
  for(const [p,a,k] of [[5,1,1],[5,1,2],[5,1,3],[13,3,1],[13,3,2],[17,8,1],[17,8,2]]) {
    const q=p**k,N=q*q,D=a*a+4,mod=x=>(x%q+q)%q;
    need(D%p===0 && D%(p*p)!==0,"Simple ramification");
    const multiply=(x,y)=>[mod(x[0]*y[0]+D*x[1]*y[1]),mod(x[0]*y[1]+x[1]*y[0])];
    const mu=[mod(a*((q+1)/2)),(q+1)/2];
    let g=multiply(mu,mu);g=multiply(g,g);
    need(mod(g[0]*g[0]-D*g[1]*g[1])===1,"Norm one");
    const index=(x,y)=>x*q+y,point=x=>[Math.floor(x/q),x%q];
    const step=x=>{const z=multiply(g,point(x));return index(z[0],z[1]);};
    const valuation=x=>{if(x===0)return k;let v=0;while(x%p===0){x/=p;v++;}return v;};
    const seen=new Uint8Array(N);let orbitCount=0,checks=0;
    const maxima=[0,0];
    for(let start=1;start<N;start++) {
      if(seen[start])continue;
      const orbit=[];let x=start;
      do {need(!seen[x],"Orbit partition");seen[x]=1;orbit.push(x);x=step(x);} while(x!==start);
      const z=point(start),j=Math.min(2*valuation(z[0]),2*valuation(z[1])+1);
      need(j<2*k && orbit.length===p**Math.floor((2*k-j)/2),"Orbit valuation length");
      orbitCount++;
      for(let e=0;e<=1;e++) {
        const frequencies=new Map();
        for(const y of orbit){const c=point(y)[e===0?1:0];frequencies.set(c,(frequencies.get(c)||0)+1);}
        const correlation=new Int32Array(q);
        for(const [u,cu] of frequencies)for(const [v,cv]of frequencies)correlation[mod(u-v)]+=cu*cv;
        let squared;
        if(e+j>=2*k) squared=orbit.length*orbit.length;
        else if((e+j)%2===0) squared=0;
        else squared=p**(k+1+Math.floor((e+j)/2)-2*Math.ceil(j/2));
        need(Number.isInteger(squared) && squared>=0 && squared<=p**(k+1),"Period bound");
        correlation[0]-=squared;
        const width=q/p;
        for(let r=0;r<width;r++)for(let h=1;h<p;h++)
          need(correlation[r]===correlation[r+h*width],"Exact cyclotomic norm identity");
        maxima[e]=Math.max(maxima[e],squared);checks++;
      }
    }
    need(orbitCount===2*(p**k-1),"Orbit census");
    reports.push({p,a,k,vertices:N,nonzeroOrbits:orbitCount,exactCyclotomicIdentities:checks,
      maximumPeriodNormSquared:maxima,upperBound:p**(k+1)});
  }
  return {status:"PASS",method:"Integer cyclic autocorrelation modulo the prime-power cyclotomic polynomial; no floating tolerance",
    characters:"Coefficient-of-delta character and its delta multiple, representing the two new conductor types",
    scope:"Finite controls for the analytic all-depth proof",reports};
}

console.log(JSON.stringify(verifyRamifiedPeriods(), null, 2));
