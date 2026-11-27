const A=x=>{if(!x)throw Error("failed")};
function mod(x,n){return (x%n+n)%n}
function gcd(a,b){while(b)[a,b]=[b,a%b];return a<0n?-a:a}
function pow(a,n,m){let r=1n;for(;n;n>>=1n,a=a*a%m)if(n&1n)r=r*a%m;return r;}
function inv(a,m){
 let [b,u,v]=[m,1n,0n];a=mod(a,m);
 while(b){const q=a/b;[a,b]=[b,a%b];[u,v]=[v,u-q*v]}
 if(a!==1n)throw Error("no inverse");return mod(u,m);
}
function value(n,sign){return (n<<n)+BigInt(sign);}
function plusWitness(n,a){const N=value(n,1);return pow(a,(N-1n)/2n,N)===N-1n;}
function ringPow(z,e,D,N){
 const mul=([a,b],[c,d])=>[mod(a*c+D*b*d,N),mod(a*d+b*c,N)];
 let r=[1n,0n];for(;e;e>>=1n,z=mul(z,z))if(e&1n)r=mul(r,z);return r;
}
function minusWitness(n,D,x,y){
 const N=value(n,-1);
 if(n<2n||mod(x*x-D*y*y,N)!==1n)return false;
 const [a,b]=ringPow([x,y],(N+1n)/2n,D,N);return a===N-1n&&b===0n;
}
function traceWitness(n,x){
 if(n<2n)return false;
 const N=value(n,-1),D=mod(x*x-1n,N);
 const [a]=ringPow([x,1n],(N+1n)/2n,D,N);
 return a===N-1n;
}
function findMinus(n){
 const N=value(n,-1);
 for(let x=0n;x<100n&&x<N;x++)if(traceWitness(n,x))return x;
 return null;
}
for(const n of [1n,141n,4713n]){
 let a=2n;while(a<100n&&!plusWitness(n,a))a++;
 A(a<100n);console.log("Cullen prime certificate: n and base",String(n),String(a))
}
for(const n of [2n,3n,6n,30n,75n,81n,115n,123n]){
 const x=findMinus(n);A(x!==null); console.log("Woodall prime certificate: n and trace seed",String(n),String(x));
}
function certificate(n,sign){
 const N=value(n,sign);
 if(N<2n)return {prime:false};
 for(let x=2n;x<N;x++){
   if((sign===1?plusWitness(n,x):traceWitness(n,x)))return {prime:true,seed:x};
   if(x*x<=N&&N%x===0n)return {prime:false,factor:x};
 }
 throw Error("candidate dichotomy failed");
}
for(const sg of [1,-1])for(let n=1n;n<=24n;n++){
 const result=certificate(n,sg);A(result.prime===isprime(value(n,sg)));
}
console.log("Complete witness-or-factor machine verified at 48 small terms.");
function isprime(n){if(n<2n)return false;for(let a=2n;a*a<=n;a++)if(n%a===0n)return false;return true;}
function lcm(a,b){return a/gcd(a,b)*b;}
function plant(a,D,H,lower=0n){
 let ps=[],B=D>20n?D:20n;
 // Any fixed finite list can be found by a prime-in-(B,2B) interval scan
 // and widening B until this dyadic interval holds 2H primes.
 while(ps.length<2*H){
   B*=2n;ps=[];
   for(let p=B+1n;p<2n*B&&ps.length<2*H;p+=2n)if(isprime(p))ps.push(p);
 }
 let L=D;for(const p of ps)L=lcm(L,p-1n);
 let x=mod(a,L), M=L, out=[];
 for(let i=0;i<H;i++)for(let s of [1,-1]){
   const p=ps[out.length],off=BigInt(i)*D,res=mod(-off-BigInt(s)*inv(pow(2n,a+off,p),p),p);
   x+=M*mod((res-x)*inv(M,p),p);M*=p;out.push([off,s,p]);
 }
 const min=ps.at(-1)>lower?ps.at(-1):lower;
 if(x<=min)x+=M*((min-x)/M+1n);
 A(x%D===a%D);
 for(const [off,s,p] of out)A(mod((x+off)*pow(2n,x+off,p)+BigInt(s),p)===0n&&x+off>p);
 return [x,M,out];
}
for(const [a,D,H] of [[0n,1n,5],[0n,30n,3],[11n,30n,2]]){
 const [x,M,out]=plant(a,D,H,1000n);
 console.log("Plant: a,D,H,start,period,primes",...[a,D,H,x,M].map(String),out.map(v=>String(v[2])));
}

