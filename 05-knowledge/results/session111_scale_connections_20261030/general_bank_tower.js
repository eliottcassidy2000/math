const need=(x,m)=>{if(!x)throw Error(m);};
function check(p,metal,k){
 const q=p**k,N=q*q,mod=x=>(x%q+q)%q,ix=([x,y])=>x*q+y,xy=i=>[Math.floor(i/q),i%q];
 const mul=(a,b)=>[mod(a[0]*b[0]+a[1]*b[1]),mod(a[0]*b[1]+a[1]*b[0]+metal*a[1]*b[1])],mu2=mul([0,1],[0,1]),g=mul(mu2,mu2),neg=a=>a.map(x=>mod(-x));
 need((metal*metal+4)%p===0&&(metal*metal+4)%(p*p)!==0,"simple ramification");
 const orb=Array(N).fill(-1),os=[];orb[0]=-2;
 for(let i=1;i<N;i++)if(orb[i]===-1){let o=[],j=i;do{o.push(j);orb[j]=os.length;j=ix(mul(g,xy(j)));need(o.length<=q,"clock order");}while(j!==i);os.push(o);}
 const bit=Array(os.length).fill(-1),sign=Array(os.length).fill(0),reps=[];
 for(let i=0;i<os.length;i++)if(bit[i]===-1){let j=orb[ix(neg(xy(os[i][0])))];need(j!==i,"antipodes");bit[i]=bit[j]=reps.length;sign[j]=1;reps.push(os[i][0]);}
 need(reps.length===q-1,"bits");
 const n=p**(k-1),b0=1+n,inv2=(q+1)/2,h=(p-1)/2;
 function transform(v,a,b){let scalar=1;for(let i=0;i<a;i++)scalar=scalar*(1+p)%q;let bb=1;for(let i=0;i<b;i++)bb=bb*b0%q;return[mod(scalar*(v[0]+metal*v[1]*(1-bb)*inv2)),mod(scalar*v[1]*bb)];}
 const gcd=(a,b)=>b?gcd(b,a%b):a,hist={};let masks=0n;
 for(let a=0;a<n;a++)for(let b=0;b<p;b++){
  const m=reps.map(i=>{const o=orb[ix(transform(xy(i),a,b))];return[bit[o],sign[o]];});
  const used=new Set();let cycles=0;for(let i=0;i<reps.length;i++)if(!used.has(i)){let j=i,s=0;do{used.add(j);s^=m[j][1];j=m[j][0];}while(j!==i);need(!s,"odd group cycle sign");cycles++;}
  let predicted=h*(gcd(a,n)+gcd((a+(n/p)*b)%n,n));for(let j=0;j<=k-2;j++)predicted+=(p-1)*gcd(a,p**j);
  need(cycles===predicted,"actual joint cycles");hist[cycles]=(hist[cycles]||0)+1;masks+=2n**BigInt(cycles);
 }
 let expected={};const put=(e,c)=>expected[e]=(expected[e]||0)+c;
 put(p*n-1,1);put(n-1+h*(n+n/p),2*(p-1));put(n-1+2*h*n/p,(p-1)**2);
 for(let v=0;v<=k-3;v++)put(((p-1)*(k-v)+1)*p**v-1,(p-1)*p**(k-1-v));
 need(JSON.stringify(Object.entries(hist).sort())===JSON.stringify(Object.entries(expected).sort()),"cycle sum");
 need(masks%BigInt(p**k)===0n,"Burnside divisibility");
 return{p,metal,k,vertices:N,bankBits:q-1,oddEffectiveGroup:p**k,cycleHistogram:hist,binaryBankOrbits:(masks/BigInt(p**k)).toString()};
}
console.log(JSON.stringify({status:"PASS",cases:[[5,1,2],[5,1,3],[13,3,2],[17,8,2]].map(x=>check(...x)),recurrence:"B_(k+1)(Q)=(Q^(p-1) B_k(Q^p)+(p-1)Q^((p-1)(k+1)))/p",scope:"Exhaustive actual signed action in four specified simple-ramification cases; general result is proved separately."},null,2));