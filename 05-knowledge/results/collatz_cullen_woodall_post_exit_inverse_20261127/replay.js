const A=x=>{if(!x)throw Error("replay failed")};
function mod(x,m){return (x%m+m)%m;}
function pow(a,n,m){let r=1n;for(;n;n>>=1n,a=a*a%m)if(n&1n)r=r*a%m;return r;}
function T(x,sign=1n){return x%2n===0n?x/2n:(3n*x+sign)/2n;}
function word(x,k,sign=1n){let w="";for(let i=0;i<k;i++){w+=String(x%2n);x=T(x,sign)}return w;}
function inverseWord(w,sign=1n){
 let x=0n,M=1n;
 for(let j=0;j<w.length;j++,M*=2n)if(word(x,j+1,sign).at(-1)!==w[j])x+=M;
 A(word(x,w.length,sign)===w);return x;
}
function inverseCenter(y,k){
 let n=1n,M=2n;
 for(let r=1;r<k;r++,M*=2n)if(mod(n*pow(3n,n,2n*M)-y,2n*M)!==0n)n+=M;
 A(mod(n*pow(3n,n,2n**BigInt(k))-y,2n**BigInt(k))===0n);return n;
}
function indexClass(w,s,eps){
 const r=inverseWord(w,-eps),k=w.length+s;
 return [inverseCenter(2n**BigInt(s)*r-eps,k),2n**BigInt(k)];
}
function post(n,s,eps){
 const y=n*3n**n+eps,h=2n**BigInt(s);
 A(y%(2n*h)===h);return y/h;
}
let checked=0;
for(const eps of [-1n,1n])for(let s=1;s<=4;s++)for(let k=1;k<=10;k++)for(let q=0;q<1<<(k-1);q++){
 const w="1"+q.toString(2).padStart(k-1,"0").slice(0,k-1);
 const [n,M]=indexClass(w,s,eps),N=n+M;
 const mm=2n**BigInt(s+k),X=mod(N*pow(3n,N,mm)+eps,mm)/2n**BigInt(s);
 A(word(X,k,-eps)===w);checked++;
}
console.log("Post-exit index inverse checks:",checked);
const w="111010",s=1,eps=-1n,[n,M]=indexClass(w,s,eps);
let p=n;while(mod(p*pow(2n,p,5n)+1n,5n)!==0n||mod(p*pow(2n,p,7n)-1n,7n)!==0n)p+=M;
A(word(post(p,s,eps),w.length,1n)===w);
console.log("Composite pair with prescribed Woodall exit:",[w,String(n),String(M),String(p)]);
function decomp(n){let k=n,e=n;while(k%2n===0n){e++;k/=2n}return [k,e]}
function cullen(n,a){const N=(n<<n)+1n,[k,e]=decomp(n);let R=pow(a,k,N);for(let j=1n;j<e;j++)R=R*R%N;return R===N-1n}
function woodall(n,x){
 const N=(n<<n)-1n,[k,e]=decomp(n),D=x*x-1n;
 const mul=([a,b],[c,d])=>[mod(a*c+D*b*d,N),mod(a*d+b*c,N)];
 let r=[1n,0n],z=[x,1n];
 for(let j=k;j;j>>=1n,z=mul(z,z))if(j&1n)r=mul(r,z);
 let t=mod(2n*r[0],N);
 for(let j=1n;j<e;j++)t=mod(t*t-2n,N);
 return t===N-2n;
}
for(const n of [141n,4713n])A(cullen(n,5n));
for(const [n,x] of [[3n,4n],[75n,6n],[123n,4n]])A(woodall(n,x));
console.log("Half-exponent countdown machines replayed on five prime terms.");
