const assert=require("node:assert/strict");
for(let m=1;m<=4;m++){
 const k=2**m,p=2**k+1; let M=1;for(let j=1;j<m;j++)M*=2**(2**j)-1;
 const sq=new Set(Array.from({length:(p-1)/2},(_,i)=>(i+1)**2%p)), best=[0,0],a0=[0,0],groups=new Map;
 for(let a=1;a<p;a++){
  let x=a,d=p;for(let j=0;j<k;j++){d=Math.min(d,x,p-x);x=2*x%p;}
  const c=sq.has(a)?0:1;if(d>best[c])best[c]=d,a0[c]=a;
  if(!groups.has(d))groups.set(d,[0,0]);groups.get(d)[c]++;
 }
 assert.equal(Math.max(...best),M);assert.equal(groups.get(M).reduce((a,b)=>a+b),2*k);
 if(m>=2)assert.equal(groups.get(M)[0],0);
 console.log(p,M,JSON.stringify(best),JSON.stringify(a0),JSON.stringify([...groups].sort((a,b)=>b[0]-a[0]).slice(0,4)));
}
for(let m=1;m<=10;m++){
 let k=2**m,p=(1n<<BigInt(k))+1n,M=1n;
 for(let j=1;j<m;j++)M*=((1n<<BigInt(2**j))-1n);
 let z=M,obs=p;const orb=new Set;
 for(let j=0;j<2*k;j++){orb.add(z);const y=z<p-z?z:p-z;if(y<obs)obs=y;z=2n*z%p;}
 assert.equal(z,M);assert.equal(orb.size,2*k);assert.equal(obs,M);
 if(m<3){if(M*BigInt(k+1)>=p)throw Error("grid unexpectedly enough");}
 if(m>=3){if(M*BigInt(k+1)<=p)throw Error("grid unexpectedly inadequate");}
}
console.log("constructive orbit checks through m=10, including composite moduli");
