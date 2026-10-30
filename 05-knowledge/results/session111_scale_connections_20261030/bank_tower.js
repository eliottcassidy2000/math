const need=(x,m)=>{if(!x)throw Error(m);};
function check(k){
 const q=5**k,N=q*q,mod=x=>(x%q+q)%q,ix=([x,y])=>x*q+y,xy=i=>[Math.floor(i/q),i%q];
 const mul=(a,b)=>[mod(a[0]*b[0]+a[1]*b[1]),mod(a[0]*b[1]+a[1]*b[0]+a[1]*b[1])],g=[2,3],neg=a=>a.map(x=>mod(-x));
 const orb=Array(N).fill(-1),os=[];orb[0]=-2;
 for(let i=1;i<N;i++)if(orb[i]===-1){let o=[],j=i;do{o.push(j);orb[j]=os.length;j=ix(mul(g,xy(j)));}while(j!==i);os.push(o);}
 const bit=Array(os.length).fill(-1),sign=Array(os.length).fill(0),reps=[];
 for(let i=0;i<os.length;i++)if(bit[i]===-1){let j=orb[ix(neg(xy(os[i][0])))];bit[i]=bit[j]=reps.length;sign[j]=1;reps.push(os[i][0]);}
 need(reps.length===q-1,"bits");
 const n=5**(k-1),b0=1+n,inv2=(q+1)/2;
 function transform(v,a,b){let scalar=1;for(let i=0;i<a;i++)scalar=scalar*6%q;let bb=1;for(let i=0;i<b;i++)bb=bb*b0%q;return[mod(scalar*(v[0]+v[1]*(1-bb)*inv2)),mod(scalar*v[1]*bb)];}
 const gcd=(a,b)=>b?gcd(b,a%b):a,hist={};let masks=0n;
 for(let a=0;a<n;a++)for(let b=0;b<5;b++){
  const m=reps.map(i=>{const o=orb[ix(transform(xy(i),a,b))];return[bit[o],sign[o]];});
  const used=new Set();let cycles=0;for(let i=0;i<reps.length;i++)if(!used.has(i)){let j=i,s=0;do{used.add(j);s^=m[j][1];j=m[j][0];}while(j!==i);need(!s,"odd group cycle sign");cycles++;}
  let predicted=2*gcd(a,n)+2*gcd((a+(n/5)*b)%n,n);for(let j=0;j<=k-2;j++)predicted+=4*gcd(a,5**j);
  need(cycles===predicted,"actual versus joint quotient cycles");hist[cycles]=(hist[cycles]||0)+1;masks+=2n**BigInt(cycles);
 }
 let expected={};const put=(e,c)=>expected[e]=(expected[e]||0)+c;
 put(5*n-1,1);put(3*n+2*n/5-1,8);put(n+4*n/5-1,16);
 for(let v=0;v<=k-3;v++)put((4*k+1-4*v)*5**v-1,20*5**(k-2-v));
 need(JSON.stringify(Object.entries(hist).sort())===JSON.stringify(Object.entries(expected).sort()),"closed cycle sum");
 return{k,vertices:N,bankBits:q-1,oddEffectiveGroup:5**k,cycleHistogram:hist,binaryBankOrbits:(masks/BigInt(5**k)).toString()};
}
const cases=[check(2),check(3)];
function cyclepoly(k){if(k===2)return{24:1n,16:8n,8:16n};const out={};for(const[e,c]of Object.entries(cyclepoly(k-1)))out[5*Number(e)+4]=c;const e=4*k;out[e]=(out[e]||0n)+20n*5n**BigInt(k-2);return out;}
for(let k=2;k<=8;k++){const h=cyclepoly(k);need(Object.values(h).reduce((a,b)=>a+b,0n)===5n**BigInt(k),"mass");for(const Q of [2n,3n])need(Object.entries(h).reduce((s,[e,c])=>s+c*Q**BigInt(e),0n)%5n**BigInt(k)===0n,"integral orbit count");}
console.log(JSON.stringify({status:"PASS",cases,recurrence:"B_(k+1)(Q)=(Q^4 B_k(Q^5)+4Q^(4k+4))/5 for k>=2",base:"B_2(Q)=(Q^24+8Q^16+16Q^8)/25",scope:"Actual signed actions exhausted at k2,k3. Formal recurrence mass and integrality checked k2..8; all-k result requires proof."},null,2));