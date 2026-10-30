function verifyArithmeticProbes() {
  const need=(ok,msg)=>{if(!ok)throw new Error(msg);};
  const safe=n=>{need(Number.isSafeInteger(n),"Exact integer range");return n;};
  function factors(n) {const out=[];for(let d=2;d*d<=n;d++)if(n%d===0){let e=0;while(n%d===0){n/=d;e++;}out.push([d,e]);}if(n>1)out.push([n,1]);return out;}
  function sigma(n) {let out=1;for(const [p,e]of factors(n)){let r=1,t=1;for(let j=1;j<=e;j++){t=safe(t*p);r=safe(r+t);}out=safe(out*r);}return out;}
  function divisorSigma(n) {let out=0;for(let d=1;d*d<=n;d++)if(n%d===0)out+=d*d===n?d:d+n/d;return out;}
  const gcd=(a,b)=>b?gcd(b,a%b):a;
  const R=(p,e)=>{let r=1,t=1;for(let j=1;j<=e;j++){t=safe(t*p);r=safe(r+t);}return r;};
  function primePowerInverse(p,t,a,b) {
    need(b>a,"Positive difference");let e=0,d=b-a;
    while(d%p===0){d/=p;e++;}
    const denominator=R(p,t-1),numerator=b-p**t*a;
    need(numerator%denominator===0,"Integral core sigma");
    return safe(R(p,e)*(numerator/denominator)-a);
  }
  let reconstructionChecks=0;
  for(let n=1;n<=1000;n++)for(const p of [2,3,5,7])for(const t of [1,2,3]){
    const a=sigma(n)-n,b=sigma(n*p**t)-n*p**t;
    need(primePowerInverse(p,t,a,b)===n,"Universal prime-power inverse");reconstructionChecks++;
  }
  function inverse(m,a,b) {
    const fs=factors(m),D=b-m*a;
    need(fs.length>=2,"Composite-port inverse domain");
    if(D<=0)return {sources:[],signatures:0};
    const choices=fs.map(([p,t])=>{const es=[];for(let value=p**t,e=0;value<=D;value*=p,e++)es.push(e);return es;});
    let signatures=0;const found=new Set();
    function visit(i,h) {
      if(i<fs.length){for(const e of choices[i])visit(i+1,safe(h*fs[i][0]**e));return;}
      signatures++;const K=sigma(m*h)-m*sigma(h);
      need(K>0,"Positive overlap coefficient");
      if(D%K)return;
      const c=D/K,N=safe(sigma(h)*c-a);
      if(N<=0 || N%h)return;
      const u=N/h;
      if(gcd(u,m)!==1 || sigma(u)!==c)return;
      need(sigma(N)-N===a && sigma(m*N)-m*N===b,"Inverse sufficiency");
      found.add(N);
    }
    visit(0,1);return {sources:[...found].sort((x,y)=>x-y),signatures};
  }
  const witnesses=[];
  for(const [m,n1,n2]of [[10,387152,437510],[15,563472,759888]]) {
    const a=sigma(n1)-n1,b=sigma(m*n1)-m*n1;
    need(sigma(n2)-n2===a && sigma(m*n2)-m*n2===b,"Composite collision");
    for(const n of [n1,n2,m*n1,m*n2])need(sigma(n)===divisorSigma(n),"Independent divisor-pair sum");
    const fibre=inverse(m,a,b);
    need(JSON.stringify(fibre.sources)===JSON.stringify([n1,n2]),"Complete two-point fibre");
    witnesses.push({m,a,b,...fibre});
  }
  let positiveControls=0;
  for(const m of [6,10,12,15])for(let n=1;n<=60;n++){
    need(inverse(m,sigma(n)-n,sigma(m*n)-m*n).sources.includes(n),"Inverse retained source");positiveControls++;
  }
  function isqrt(n) {if(n<2n)return n;let x=n,y=(x+1n)/2n;while(y<x){x=y;y=(x+n/x)/2n;}return x;}
  function codeBank(primes) {
    const support=[...new Set(primes.flatMap(p=>factors(p+1).filter(x=>x[1]%2).map(x=>x[0])))].sort((a,b)=>a-b);
    const columns=primes.map(p=>factors(p+1).reduce((v,[q,e])=>e%2?v^(1<<support.indexOf(q)):v,0));
    const observed=columns.map(c=>[2,3,5].reduce((v,q,i)=>support.includes(q)&&((c>>support.indexOf(q))&1)?v^(1<<i):v,0));
    const rank=cs=>{const pivots=[];let r=0;for(let x of cs){for(let b=31;b>=0;b--)if(x>>>b&1){if(pivots[b])x^=pivots[b];else{pivots[b]=x;r++;break;}}}return r;};
    let projected=0,actual=0;
    for(let mask=0;mask<128;mask++){
      let full=0,shadow=0,product=1n;
      for(let j=0;j<7;j++)if(mask>>j&1){full^=columns[j];shadow^=observed[j];product*=BigInt(primes[j]+1);}
      if(shadow===0)projected++;
      const square=isqrt(product)**2n===product;
      need(square===(full===0),"Square test versus full syndrome");
      if(square)actual++;
    }
    return {primes,support,fullRank:rank(columns),observedRank:rank(observed),projectedSquareSelections:projected,actualSquareSelections:actual};
  }
  const falseFano=codeBank([7,2,41,19,109,1019,389]),genuineFano=codeBank([7,2,19,5,89,59,29]);
  need(falseFano.fullRank===7 && falseFano.observedRank===3 && falseFano.projectedSquareSelections===16 && falseFano.actualSquareSelections===1,"False Fano firewall");
  need(genuineFano.fullRank===3 && genuineFano.actualSquareSelections===16,"Inherited genuine Fano control");
  return {status:"PASS",reconstructionChecks,compositePositiveControls:positiveControls,witnesses,falseFano,genuineFano,
    scope:"Exact finite controls for the proved universal inverse and projection theorems; no minimal-collision claim"};
}

console.log(JSON.stringify(verifyArithmeticProbes(), null, 2));
