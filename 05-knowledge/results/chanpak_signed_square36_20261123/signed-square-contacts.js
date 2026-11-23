const bit=(v,w)=>{let a=Math.min(v,w),b=Math.max(v,w);return a===1&&b===3?1:a===7&&b===9?2:0};
const sq=x=>Number.isInteger(Math.sqrt(x));
function contact(N){let out=[];for(let x=1;x<=N;x+=2)for(let s of [-1,1]){
  let y=3*x+s,a=0;while(!(y%2)){y/=2;a++}
  if(y!==x&&y<=N&&sq(x+y))out.push([x,y,s,a]);
}return out;}
function census(N){
 const A=Array.from({length:N+1},(_,x)=>Array.from({length:N},(_,i)=>i+1).filter(y=>x!==y&&sq(x+y)));
 const c=[0,0,0,0], witness={};
 if(A.slice(1).some(a=>a.length<2))return [c,witness];
 const r=A.findIndex((a,i)=>i>0&&a.length===2)>0?A.findIndex((a,i)=>i>0&&a.length===2):1,V=new Uint8Array(N+1),P=[r];V[r]=1;
 function dfs(v,b){
   if(P.length===N){if(A[v].includes(r)&&P[1]<v){let z=b|bit(v,r);c[z]++;witness[z]??=P.slice()}return;}
   for(let i=1;i<=N;i++)if(!V[i]&&A[i].filter(w=>!V[w]||w===v||w===r).length<2)return;
   for(let w of A[v])if(!V[w]){V[w]=1;P.push(w);dfs(w,b|bit(v,w));P.pop();V[w]=0;}
 }
 dfs(r,0);return [c,witness];
}

const eq=(a,b)=>JSON.stringify(a)===JSON.stringify(b);
const assert=(x,m)=>{if(!x)throw Error(m)};
assert(eq(contact(36),[[3,1,-1,3],[9,7,1,2]]),"Q36 contacts");
const counts=[[1,0,0,0],[1,0,0,0],[1,10,0,0],[29,10,15,3],[15,6,9,1]];
for(let n=3;n<=36;n++){
 const [c,w]=census(n);assert(eq(c,n<32?[0,0,0,0]:counts[n-32]),"cycle census "+n);
 if(n>=31)console.log(n,c,w[3]||null);
}
const sigma=x=>(x&2)?-1:1, rho=x=>(3*x+sigma(x))/4;
const odd=Array.from({length:18},(_,i)=>2*i+1);
assert(eq(odd.map(rho),Array.from({length:27},(_,i)=>i+1).filter(i=>i%3)),"odd chart");
let edge=0,equal=0;
function guard(x){let v=3*x+sigma(x),a=0;while(!(v&1)){v/=2;a++;}return a;}
for(let a=1;a<=36;a++)for(let b=a+1;b<=36;b++){
 if((a&1)&&(b&1)){
  assert(sq(a+b)===((rho(a)+rho(b))%3===0&&sq((rho(a)+rho(b))/3)),"chart graph");
  if(sq(a+b)){edge++;equal+=+(guard(a)===guard(b));}
 }
 if(sq(a+b))assert((sigma(a)!==sigma(b))===Boolean((a&1)&&(b&1)),"sign edge");
}
assert(edge===16&&equal===5,"equal depth count");
const fixedMinus=[[1],[5,7],[17,25,37,55,41,61,91]];
for(const orb of fixedMinus)for(let i=0;i<orb.length;i++){
 let v=3*orb[i]-1;while(v%2===0)v/=2;
 assert(v===orb[(i+1)%orb.length],"known minus cycles");
}

const H=(r,e=[])=>{r=r.slice();for(let [i,j] of e){r[i]^=1<<j;r[j]^=1<<i}let n=r.length,N=1<<n, d=new Float64Array(N*n);for(let i=0;i<n;i++)d[(1<<i)*n+i]=1; for(let s=1;s<N;s++)for(let z=s;z;z&=z-1){let i=31-Math.clz32(z&-z),v=d[s*n+i];if(v)for(let t=r[i]&~s;t;t&=t-1){let j=31-Math.clz32(t&-t);d[(s|1<<j)*n+j]+=v}}return d.subarray((N-1)*n).reduce((a,b)=>a+b,0)};
const D=[[57,43,[32,57,43,113,101,64,7],[3,4]],[139,185,[8,141,9,160,111,71,143,53],[2,7]],[329,247,[218,136,139,32,110,199,14,88],[0,5]],[157,419,[58,80,107,2,172,74,25,111],[3,7]],[439,585,[164,137,74,17,7,94,27,124],[3,5]]];
for(const [x,y,A,e] of D){for(let i=0;i<A.length;i++)for(let j=0;j<A.length;j++)if(((A[i]>>j&1)+(A[j]>>i&1))!==+(i!==j))throw Error('orientation');
 if(H(A)!==x||H(A,[e])!==y||!sq(x+y))throw Error('seam');}
if(JSON.stringify(contact(56))!==JSON.stringify([[3,1,-1,3],[9,7,1,2]]))throw Error('first contact');
let mx=0;
for(let z=0;z<1024;z++){let A=Array(6).fill(0),k=0;for(let i=0;i<6;i++)for(let j=i+1;j<6;j++){
if(j===i+1||(z>>k++)&1)A[i]|=1<<j;else A[j]|=1<<i;}mx=Math.max(mx,H(A));}
if(mx!==45)throw Error('six-vertex max');

const strong=r=>{
 let n=r.length,s=1,t=0;
 while(t!==s){t=s;for(let i=0;i<n;i++)if((s>>i)&1)s|=r[i];}
 let p=1; t=0;
 while(t!==p){t=p;for(let i=0;i<n;i++)if(r[i]&p)p|=1<<i;}
 return s===(1<<n)-1&&p===s;
};
const T=D[0][2].slice(),[a,b]=D[0][3];
assert(strong(T),"57 strong");T[a]^=1<<b;T[b]^=1<<a;assert(strong(T),"43 strong");
assert(eq(contact(57),[[3,1,-1,3],[9,7,1,2],[57,43,1,2]]),"tail threshold");
