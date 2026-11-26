// Exact lower bounds, all integer arithmetic. D is n times the exterior cut surplus.
function separation(iterations=8) {
  const n=3024n, ceil=(a,b)=>(a+b-1n)/b;
  const baseSupport=a=>Math.max(a*a+1,2*a+1,2*a*(57-a)-(a-1)-a*(a-1)/2-(a*a-4*a),a===4?421:0);
  const results=[];
  for(let ai=4;ai<=53;ai++){
    const a=BigInt(ai);
    let min=Infinity,arg=null;
    for(let ui=baseSupport(ai);ui<=Number(n);ui++){
      const u=BigInt(ui);
      const base=46n*u*(n-u), KN=(50n-a)*n-46n*u;
      const marked=KN>0n?2n*a:0n;
      const K=KN>0n?KN:0n;
      const g=(L)=>L+((u-L)&1n);
      let cut=g(ceil(base,n));
      const cutMax=53n*(u<n-u?u:n-u);
      if(u<n){
        const l=u-marked, r=n-u;
        if(l<=0n)throw Error("missing residual interior");
        // min of squared sum over a common lattice coset, with a lower
        // bound on the nonnegative mean; interpolate consecutive squares.
        const sqBound=(D,k,res)=>{
          if(2n*res<=n && D<k*res) D=k*res;
          const t=D-k*res, d=k*n;
          const p=res+n*(t>=0n?t/d:-ceil(-t,d));
          return k*p*p+(D-k*p)*(2n*p+n);
        };
        const res=46n*u%n, inv=(n-res)%n;
        for(let i=0;i<iterations && cut<=cutMax;i++){
          const D=cut*n-base;
          const S=marked*K*K+sqBound(D+marked*K,l,res)+sqBound(D,r,inv);
          const next=g(ceil(15n*n*base+S,15n*n*n));
          if(next<=cut)break;
          cut=next;
        }
      }
      if(cut>cutMax)continue;
      const M=Number(ceil(a*(53n*u+cut),2n*(49n+a)));
      if(M<min){min=M;arg=[ui,Number(cut)]}
    }
    if(!arg)throw Error("no numerical support candidate");
    results.push({a:ai,M:min,u:arg[0],cut:arg[1]});
  }
  return results;
}

const r=separation(12);
if(Math.min(...r.map(x=>x.M))!==1523)throw Error("unexpected minimum");
console.log(r.slice(0,8));
