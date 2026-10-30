const need=(x,m)=>{if(!x)throw Error(m);},q=25,N=625,mod=x=>(x%q+q)%q,ix=([x,y])=>x*q+y,xy=i=>[Math.floor(i/q),i%q];
const mul=(a,b)=>[mod(a[0]*b[0]+a[1]*b[1]),mod(a[0]*b[1]+a[1]*b[0]+a[1]*b[1])],g=[2,3],neg=a=>a.map(x=>mod(-x));
const orb=Array(N).fill(-1),os=[];orb[0]=-2;
for(let i=1;i<N;i++)if(orb[i]===-1){let o=[],j=i;do{o.push(j);orb[j]=os.length;j=ix(mul(g,xy(j)));}while(j!==i);os.push(o);}
const bit=Array(os.length).fill(-1),sign=Array(os.length).fill(0),reps=[];
for(let i=0;i<os.length;i++)if(bit[i]===-1){let j=orb[ix(neg(xy(os[i][0])))];bit[i]=bit[j]=reps.length;sign[j]=1;reps.push(os[i][0]);}
need(reps.length===24,"bank bits");
const scalar=a=>a.map(x=>mod(6*x)),alpha=([x,y])=>[mod(x+10*y),mod(6*y)];
function power(f,a,e){while(e--)a=f(a);return a;}

const units=Array.from({length:25},(_,i)=>i).filter(x=>x%5!==0);
const autos=units.filter(b=>5*(b*b-1)%25===0),signed=new Set(),unsigned=new Set(),kernel=[];
for(const u of units)for(const b of autos){
 const m=reps.map(i=>{const [x,y]=xy(i),v=[mod(u*(x+y*(1-b)*13)),mod(u*y*b)],o=orb[ix(v)];return[bit[o],sign[o]];});
 signed.add(JSON.stringify(m));unsigned.add(JSON.stringify(m.map(x=>x[0])));
 if(m.every((v,i)=>v[0]===i))kernel.push({scalar:u,deltaMultiplier:b,signPattern:m.map(x=>x[1])});
}
need(signed.size===200&&unsigned.size===50,"signed/unsigned image sizes");
need(kernel.length===4,"four-element effective kernel");
need(kernel.map(x=>x.scalar+","+x.deltaMultiplier).sort().join("|")==="1,1|1,24|24,1|24,24","kernel generators");
const rank=(rows)=>{const b={};for(let v of rows){while(v){let j=v.toString(2).length-1;if(b[j])v^=b[j];else{b[j]=v;break;}}}return Object.keys(b).length;};
need(rank(kernel.map(x=>x.signPattern.reduce((v,b,i)=>v|(BigInt(b)<<BigInt(i)),0n)))===2,"two independent lost signs");
console.log(JSON.stringify({status:"PASS",normalizerOrder:5000,clockKernel:25,signedImage:signed.size,unsignedKernel:100,unsignedImage:unsigned.size,effectiveForgetfulKernel:kernel,forgottenSignRank:2,scope:"Exact finite action on all24 coordinates; unsigned means forgetting flips, not identifying tournament graph classes."},null,2));