function assert(x){if(!x)throw Error("replay failed")}
function mod(x,m){return (x%m+m)%m}
function pow(a,n,m){let r=1n;for(;n;n>>=1n,a=a*a%m)if(n&1n)r=r*a%m;return r}
function value(m,power){return mod(m*pow(2n,m,power),power)}
function lift(p,L,j,r,target){
 assert(pow(2n,L,p)===1n&&j>=0n&&j<L);
 let v=0, P=1n, d=L;
 while(d%p===0n){v++;d/=p;P*=p}
 if(value(j,p**BigInt(Math.min(r,v)))!==mod(target,p**BigInt(Math.min(r,v))))return null;
 if(r<=v)return [0n,1n];
 const inv=pow(d*pow(2n,j,p),p-2n,p);
 let t=0n,Q=1n;
 for(let level=v;level<r;level++){
   const val=value(j+L*t,P*p);
   const delta=mod(target-val,P*p);
   assert(delta%P===0n);
   t+=Q*mod(delta/P*inv,p);
   P*=p;Q*=p;
 }
 assert(value(j+L*t,p**BigInt(r))===mod(target,p**BigInt(r)));
 return [t,Q];
}
let passes=0;
for(let p of [3n,5n,7n,11n,13n])for(let m=1;m<=3;m++){
 const L=(p-1n)*3n;
 for(let j=0n;j<L;j++)for(let eps of [-1n,1n]){
   const ans=lift(p,L,j,m,eps);
   let choices=[];for(let t=0n;t<p**BigInt(m);t++)if(value(j+L*t,p**BigInt(m))===mod(eps,p**BigInt(m)))choices.push(t);
   const sample=ans?choices.filter(x=>x%ans[1]===ans[0]):[];
   assert(sample.length===choices.length&&(!ans||choices.length===Number(p**BigInt(m)/ans[1])));
   passes++;
 }
}
console.log("Clock/phase/depth/sign combinations:",passes);
const a=lift(3n,6n,2n,3,-1n), b=lift(7n,6n,2n,2,1n);
assert(lift(3n,6n,0n,1,-1n)===null && lift(3n,6n,2n,1,1n)===null);
let t=0n;while(t%a[1]!==a[0]||t%b[1]!==b[0])t++;
const m=2n+6n*t;
console.log("C divisible by 27 and W divisible by 49 in phase 2 mod 6:",String(m),String(6n*a[1]*b[1]));
assert(value(m,27n)===26n&&value(m,49n)===1n);
function pair(n){
 let y=n<<n,x=y-1n,z=y+1n,rows=[[x,z]],steps=0;
 while((x&1n)&&(z&1n)){
   x=(3n*x+1n)/2n;
   z=(3n*z-1n)/2n;
   assert(z-x===2n);
   steps++;rows.push([x,z]);
 }
 let v=0,m=n;while(m%2n===0n){v++;m/=2n}
 assert(steps===Number(n)+v&&x%2n===0n&&z%2n===0n&&(x%4n===0n)!==(z%4n===0n));
 return rows.map(r=>r.map(String));
}
for(let n=1n;n<=100n;n++)pair(n);
console.log("Signed Collatz for n=3:",pair(3n));
