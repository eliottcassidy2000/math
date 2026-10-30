const need=(x,m)=>{if(!x)throw Error(m);},q=25,N=625,mod=x=>(x%q+q)%q,ix=([x,y])=>x*q+y,xy=i=>[Math.floor(i/q),i%q];
const mul=(a,b)=>[mod(a[0]*b[0]+a[1]*b[1]),mod(a[0]*b[1]+a[1]*b[0]+a[1]*b[1])],g=[2,3],neg=a=>a.map(x=>mod(-x));
const orb=Array(N).fill(-1),os=[];orb[0]=-2;
for(let i=1;i<N;i++)if(orb[i]===-1){let o=[],j=i;do{o.push(j);orb[j]=os.length;j=ix(mul(g,xy(j)));}while(j!==i);os.push(o);}
const bit=Array(os.length).fill(-1),sign=Array(os.length).fill(0),reps=[];
for(let i=0;i<os.length;i++)if(bit[i]===-1){let j=orb[ix(neg(xy(os[i][0])))];bit[i]=bit[j]=reps.length;sign[j]=1;reps.push(os[i][0]);}
need(reps.length===24,"bank bits");
const scalar=a=>a.map(x=>mod(6*x)),alpha=([x,y])=>[mod(x+10*y),mod(6*y)];
function power(f,a,e){while(e--)a=f(a);return a;}
const maps=[];
for(let a=0;a<5;a++)for(let b=0;b<5;b++)maps.push(reps.map(i=>{let v=power(alpha,power(scalar,xy(i),a),b),o=orb[ix(v)];return[bit[o],sign[o]];}));
const joint=[],seen=new Set();
for(let i=0;i<24;i++)if(!seen.has(i)){const o=[...new Set(maps.map(m=>m[i][0]))].sort((a,b)=>a-b);o.forEach(x=>seen.add(x));const kernel=[];maps.forEach((m,j)=>{if(o.every(x=>m[x][0]===x))kernel.push([Math.floor(j/5),j%5]);});joint.push({size:o.length,kernel});}
const fixedHist={};let total=0n;
for(const m of maps){let used=new Set(),cycles=0,valid=true;for(let i=0;i<24;i++)if(!used.has(i)){let j=i,s=0;do{used.add(j);s^=m[j][1];j=m[j][0];}while(j!==i);if(s)valid=false;cycles++;}
 const f=valid?2n**BigInt(cycles):0n;fixedHist[f.toString()]=(fixedHist[f.toString()]||0)+1;total+=f;}
need(fixedHist["16777216"]===1&&fixedHist["65536"]===8&&fixedHist["256"]===16,"fixed histogram");
need(total/25n===692224n,"Burnside");
need(joint.map(x=>x.size).sort((a,b)=>a-b).join()==="1,1,1,1,5,5,5,5","joint orbit structure");
const nontrivial=joint.filter(x=>x.size===5),kernels={};for(const x of nontrivial){let key=JSON.stringify(x.kernel);kernels[key]=(kernels[key]||0)+1;}need(Object.values(kernels).sort().join()==="2,2","two repeated kernels");
// A joint fixed signed assignment exists by orbit averaging over an odd group.
let gauge=Array(24).fill(0);for(const m of maps)for(let i=0;i<24;i++)gauge[m[i][0]]^=m[i][1];
for(const m of maps)for(let i=0;i<24;i++)need((gauge[i]^m[i][1])===gauge[m[i][0]],"simultaneous sign gauge");
function binrow(n){const a=[1n];for(let k=1;k<=n;k++)a.push(a[k-1]*BigInt(n-k+1)/BigInt(k));return a;}
function polymul(a,b){const c=Array(a.length+b.length-1).fill(0n);for(let i=0;i<a.length;i++)for(let j=0;j<b.length;j++)c[i+j]+=a[i]*b[j];return c;}
const B=binrow(10).map((x,k)=>(x+(k%5===0?4n*binrow(2)[k/5]:0n))/5n),F=polymul(binrow(4),polymul(B,B));
let Q=[F.at(-1)];for(let k=F.length-2;k>=0;k--){const R=Array(Q.length+1).fill(0n);Q.forEach((x,j)=>{R[j]-=x;R[j+1]+=x;});R[0]+=F[k];Q=R;}
const depths=[Q.findIndex(x=>x!==0n),Q.findIndex(x=>x%5n!==0n),Q.findIndex(x=>x%2n!==0n)];
need(depths.join()==="8,16,24","weighted depths");need(F.reduce((a,b)=>a+b,0n)===692224n,"weighted sum");
console.log(JSON.stringify({weightedPolynomial:F.map(String),depths_Q_5_2:depths,status:"PASS",group:"C5 x C5",bits:24,jointOrbits:joint,fixedHist,bankOrbits:(total/25n).toString(),factorization:"16 * ((2^10 + 4*2^2)/5)^2",signGauge:gauge,scope:"Actual signed normalizer subgroup action; graph classes are not counted."},null,2));