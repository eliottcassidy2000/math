const A=x=>{if(!x)throw Error("failed replay")};
function gcd(a,b){while(b)[a,b]=[b,a%b];return a}
function mod(a,b){return (a%b+b)%b}
function pow(x,n,m){let z=1;for(;n;n>>=1,x=x*x%m)if(n&1)z=z*x%m;return z;}
function beta(n,j){const d=gcd(n,j);return n/2%d?2n**BigInt(d/2):0n}
function sum(n,f){let s=0n;for(let j=0;j<n;j++)s+=f(j);A(s%BigInt(n)===0n);return s/BigInt(n)}
function sqCount(p,h){
 const n=(p-1)/h,c=sum(n,j=>beta(n,j)),b=sum(n,j=>beta(n,j)**2n);
 const a=sum(n*p,j=>beta(n*p,j)*beta(n,j));return [a,b,c*c,a-b+c*c]
}
const C=sqCount(23,11);A(C[3]===364723n);console.log("529 linked, saturated-linked, saturated-independent, actual",C.map(String));
function paley(p){const Q=new Set();for(let i=1;i<p;i++)Q.add(i*i%p);return Q;}
function calc(m,D){
 const H={};
 for(let t=1;t<m;t++){let r=D.has(t)?1:0;for(const x of D)if(D.has(mod(t-x,m)))r++;H[r]=(H[r]||0)+1}
 return H;
}
const Q=paley(23), D=new Set();
for(let x=1;x<529;x++)if(Q.has(x%23||x/23))D.add(x);
const B=new Set(D),g=399, O=[];
for(let z=1;;){O.push(z);B.delete(z);B.add(529-z);z=z*g%529;if(z===1)break;}
for(let x=1;x<529;x++)A(B.has(x)===B.has(x*g%529)&&B.has(x)!==B.has(529-x));
console.log("Paley product, one-orbit flip",calc(529,D),calc(529,B),O);
function C5paths(){
 function adj(a,b){return mod(b-a,5)<=2&&a!==b}
 const w=[0,0,0,0,0,0];
 function perm(xs){
  if(xs.length===5){
    for(let cuts=0;cuts<16;cuts++){
      let j=1,ok=true;
      for(let i=0;i<4;i++)if(cuts>>i&1)j++;else if(!adj(xs[i],xs[i+1]))ok=false;
      if(ok)w[j]++;
    }return;
  }
  for(let i=0;i<5;i++)if(!xs.includes(i))perm(xs.concat(i));
 }
 perm([]);
 const F=Array.from({length:6**5},()=>Array(5).fill(0n)), step=[1,6,36,216,1296];for(let i=0;i<5;i++)F[step[i]][i]=1n;
 let H=0n;
 for(let M=1;M<F.length;M++){
  const count=step.map(d=>Math.floor(M/d)%6), weight=count.reduce((p,k)=>p*BigInt(w[k]),1n);
  for(let i=0;i<5;i++)if(F[M][i]){
    H+=F[M][i]*weight;
    for(let j=0;j<5;j++)if(adj(i,j)&&count[j]<5)F[M+step[j]][j]+=F[M][i];
  }
 }
 A(H===2398512090439546875n);console.log("C5[C5] per-part cover counts and exact HP",w.slice(1),String(H));
}
C5paths();
for(const [p,h] of [[5,1],[7,3],[11,5],[19,3],[23,11],[31,15]])console.log("p, order, residue at p^2",p,h,mod(pow(9,h,p*p)-pow(4,h,p*p),p*p));
