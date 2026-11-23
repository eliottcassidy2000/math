function Hdp(r){const n=r.length,F=(1<<n)-1,D=new Int32Array((F+1)*n);
for(let i=0;i<n;i++)D[(1<<i)*n+i]=1;
for(let s=1;s<=F;s++)for(let i=0;i<n;i++)if(D[s*n+i])for(let b=r[i]&(F^s);b;b&=b-1){const t=b&-b;D[(s|t)*n+31-Math.clz32(t)]+=D[s*n+i];}
return D;}
function H(r){const n=r.length,D=Hdp(r);return D.slice(((1<<n)-1)*n).reduce((x,y)=>x+y,0);}
function strong(r){let s=1,z=0;while(z!==s){z=s;for(let i=0;i<r.length;i++)if(r[i]&s)s|=1<<i;}return s===(1<<r.length)-1;}
function rows(n,t){const R=Array(n).fill(0);let k=0;for(let j=1;j<n;j++)for(let i=0;i<j;i++){if(j===i+1||(t>>k++)&1)R[i]|=1<<j;else R[j]|=1<<i;}return R;}
function rooted(r){const n=r.length,F=(1<<n)-1,D=Hdp(r),S=Hdp(r.map((a,i)=>a^F^(1<<i))),hs=new Int32Array(F+1);hs[0]=1;const u=Array(n).fill(0);for(let m=1;m<=F;m++)for(let i=0;i<n;i++)hs[m]+=D[m*n+i];for(let m=1;m<=F;m++)for(let i=0;i<n;i++)u[i]+=hs[F^m]*S[m*n+i];return [hs[F],u,D.slice(F*n)];}
function prime(x){if(x<2)return false;for(let i=2;i*i<=x;i++)if(x%i===0)return false;return true;}


const expected=[5,11,13,17,19,23,29,31,37,41,43];
const stats=new Map(expected.map(p=>[p,[0,0,0,0]])),max={};
let primeFrames=0, rootedHits=0;
for(let n=3;n<=7;n++) for(let t=0;t<1<<((n-1)*(n-2)/2);t++){
 const r=rows(n,t); if(!strong(r))continue;
 const h=H(r);max[n]=Math.max(max[n]??0,h);
 if(h<5||h===7||!prime(h))continue;
 const [p,u,x]=rooted(r), a=stats.get(p);
 primeFrames++;
 if(a)a[0]++;
 for(let i=0;i<n;i++){
  if(u[i]-x[i]===6*p){rootedHits++;if(a)a[1]++;}
  if(u[i]===3*p){rootedHits++;if(a)a[2]++;}
  if(a)for(let j=i+1;j<n;j++){
   r[i]^=1<<j;r[j]^=1<<i;
   if(H(r)===7*p)a[3]++;
   r[i]^=1<<j;r[j]^=1<<i;
  }
 }
}
if(rootedHits||[...stats.values()].some(a=>!a[0]||a.slice(1).some(Boolean)))throw Error("low-order census");
console.log({stats:[...stats],max,primeFrames,rootedHits});
for(const [p,C,i,k,expectedU,expectedE] of [
 [163,[118,60,88,209,32,12,50,119],7,1,983,5],
 [73,[252,213,248,242,32,130,48,80],1,2,219,null]
]){
 const [h,U,E]=rooted(C),n=C.length+k,F=(1<<n)-1;
 if(h!==p||U[i]!==expectedU||(expectedE!==null&&E[i]!==expectedE))throw Error("profile");
 const r=C.map(a=>a|(F^((1<<C.length)-1)));
 for(let j=C.length;j<n;j++)r.push(F^((1<<(j+1))-1));
 if(H(r)!==p)throw Error("old height");
 r[i]^=1<<(n-1);r[n-1]^=1<<i;
 if(H(r)!==7*p)throw Error("literal factor seven");
 // Forward and backward reachability give strength.
 let t=1,q=0;
 while(t!==q){q=t;for(let j=0;j<n;j++)if(t>>j&1)t|=r[j];}
 let z=1;q=0;
 while(z!==q){q=z;for(let j=0;j<n;j++)if(r[j]&z)z|=1<<j;}
 if(t!==F||z!==F)throw Error("terminal strength");
 console.log({p,q:7*p,n,U:U[i],ending:E[i],r});
}
