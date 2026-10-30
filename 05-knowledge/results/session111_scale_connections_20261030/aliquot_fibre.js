const need=(x,m)=>{if(!x)throw Error(m);};
function sigma(n){let s=0;for(let d=1;d*d<=n;d++)if(n%d===0)s+=d+(d*d===n?0:n/d);return s;}
function factors(n){const out=[];for(let p=2;p*p<=n;p++){let e=0;while(n%p===0){n/=p;e++;}if(e)out.push([p,e]);}if(n>1)out.push([n,1]);return out;}
const target=343,bound=(target-1)**2,s=new Uint32Array(bound+1);
for(let d=1;d<=bound/2;d++)for(let n=2*d;n<=bound;n+=d)s[n]+=d;
const fibre=[];for(let n=1;n<=bound;n++)if(s[n]===target)fibre.push(n);
const direct=[];for(let n=1;n<=bound;n++)if(sigma(n)-n===target)direct.push(n);
need(JSON.stringify(fibre)===JSON.stringify(direct),"independent divisor-pair enumeration");
need(fibre.length===19,"fibre count");need(fibre.filter(n=>n%2===0).join()==="578","unique even source");
const semiprimes=fibre.filter(n=>{let f=factors(n);return f.length===2&&f.every(x=>x[1]===1);});
need(semiprimes.length===17,"semiprime classification");for(const n of semiprimes){let f=factors(n);need(f[0][0]+f[1][0]===342,"Goldbach sum");}
need(fibre.filter(n=>!semiprimes.includes(n)).join()==="578,1001","exceptional types");
need(sigma(578)-578===343&&sigma(1001)-1001===343,"merge");
need(sigma(343)-343===57&&343+57===20**2,"square edge");
need(17**2+3*17+3===7**3,"cyclotomic identity");
console.log(JSON.stringify({status:"PASS",target,completeBound:bound,fibre,classified:fibre.map(n=>({n,factors:factors(n)})),uniqueEvenSource:578,semiprimeCount:17,orbitMerge:[{source:578,path:[578,343,57]},{source:1001,path:[1001,343,57]}],scope:"Complete finite fibre from proved elementary bound; global perfect-power uniqueness is a cited theorem, not this computation."},null,2));