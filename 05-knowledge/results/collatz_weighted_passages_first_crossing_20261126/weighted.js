// Node.js: node weighted.js. All mathematical comparisons use BigInt.
const S=10n**51n, N=1n<<71n, K=137500000000n;
function assert(v,msg='certificate failed'){if(!v)throw Error(msg);}
const ceil=(a,b)=> (a+b-1n)/b;
function logBounds(n,d){
  const a=n-d,b=n+d;
  assert(0n<=a&&a<b);
  let lo=0n,hi=0n,x=a,y=b;
  for(let j=0n;j<100n;j++){
    lo+=(2n*S*x)/(y*(2n*j+1n));
    hi+=ceil(2n*S*x,y*(2n*j+1n));
    x*=a*a;y*=b*b;
  }
  // Upper bound for the omitted tail via a geometric series.
  return [lo,hi+ceil(2n*S*x*b*b,y*201n*(b*b-a*a))];
}
const [L2,U2]=logBounds(2n,1n), [L3,U3]=logBounds(3n,1n);
const al=L3*S/U2, ah=ceil(U3*S,L2);
const Gnum=158496251n,Gden=100000000n;
assert(ah*Gden<Gnum*S&&133n*46057n*S<9n*1000000n*L2);
const bLog=K.toString(2).length-1, [,extraLog]=logBounds(K,1n<<BigInt(bLog));
const logK=BigInt(bLog)*S+ceil(extraLog*S,L2);
const Gcap=ceil(Gnum*S,Gden);
const powerOfTwo=x=>x>1n?BigInt((x-1n).toString(2).length):0n;
function library(num,den){
  let C=1n,L=0n,cut=0n,mp=0n,mq=1n;
  for(let b=1n;b<=32n;b++)for(let a=1n;a<=b;a++){
    const p=3n**a,q=2n**b,d=p-q;
    if(p<=q||a*den>num*b)continue;
    const B=3n**(a-1n)+2n**(b-a)*(p-2n**a-3n**(a-1n));
    if(ceil(B,d)>C)C=ceil(B,d);
    if(powerOfTwo(d)>L)L=powerOfTwo(d);
    if(a>=2n&&powerOfTwo(ceil(2n**(b-a+1n),3n))>cut)cut=powerOfTwo(ceil(2n**(b-a+1n),3n));
    if(a*p*mq>mp*d){mp=a*p;mq=d;}
  }
  return {C,L,cut,mp,mq};
}
const general=library(1n,1n);
assert((4n*(N+general.C))**200n>2n**83n*(3n*N+1n+4n*general.C)**200n);
const caps=[[1n,1n],[2n,3n],[7n,11n]];
function periodUpper(r,over,cap,lib){
  let P=S, add=0n;
  for(let i=0;i<r;i++){
    const [a,b]=i<over?caps[0]:cap;
    add+=ceil(a*P,b);
    P=ceil(P*a*Gcap,b*S);
  }
  add+=ceil((P-S)*200n,83n);
  const {cut,L,mp,mq}=lib, logPrice=powerOfTwo(ceil(BigInt(r)*(mp+5n*mq),mq));
  const e=ceil(143n*logK,10n), f=logPrice*S+ceil(133n*logK,10n);
  // n_min < 2^8 times the smaller of the two envelopes, source cut adds 'cut'.
  const H=(cut+L+9n)*S+(e<f?e:f);
  return ceil(H*add,S*S);
}
for(const cap of caps){
  const lib=library(...cap);
  let r=0;while(periodUpper(r+1,0,cap,lib)<K)r++;
  console.log('uniform',cap.join('/'),lib.C.toString(),lib.L.toString(),lib.cut.toString(),r,periodUpper(r,0,cap,lib).toString(),periodUpper(r+1,0,cap,lib).toString());
}
for(const q of [1,2,5,10]){
 let r=q;while(periodUpper(r+1,q,caps[2],general)<K)r++;
 console.log('unrestricted',q,'excluded total',r,periodUpper(r,q,caps[2],general).toString(),periodUpper(r+1,q,caps[2],general).toString());
}

// First-crossing Farey barriers. Differences and all thresholds are directed.
function farey(l,h){
 let a=0n,b=1n,c=1n,d=0n;
 for(let steps=0;steps<100000;steps++){
  assert(c*b-a*d===1n);
  const p=a+c,q=b+d;
  if(p*S<l*q){a=p;b=q;}
  else if(p*S>h*q){c=p;d=q;}
  else return [p,q,a,b,c,d];
 }
 throw Error('unresolved Farey interval');
}
const baseErr=ceil(S*S,3n*N*L2);
const early=farey(al,ah+baseErr);
const q=190537n,p=301994n;
assert(q*q*al>q*p*S-S&&q*q*ah<q*p*S+S);
const mean=ceil(S*S,2n*L2)+ceil(3n*S,2n*q)+ceil(q*S,early[1]);
assert(mean*10000n<7214n*S);
console.log('old Farey',early.map(String).join(' '));
console.log('mean scale upper',mean.toString(),'/',S.toString());
for(const [m,v] of [[1n,1n],[133163n,100000n],[4n,3n],[2n,1n],[4n,1n],[8n,1n],[16n,1n]]){
 const err=ceil(7214n*S*S*v,10000n*3n*N*L2*m);
 const f=farey(al,ah+err);
 assert((f[0]-1n)*S < f[1]*al);
 console.log('source/N',String(m)+'/'+String(v),'first numerator/denominator',f.slice(0,2).map(String).join('/'),'safe through',String(f[0]-1n));
}
const target=farey(al,ah+ceil(7214n*S*S*250n,10000n*3n*N*L2*333n));
assert(target[0]===217976794617n && (target[0]-1n)*S < target[1]*al);
const critical=early[0]*S-early[1]*ah;
assert(critical>0n);
// If n/N exceeds this, the first previously admissible denominator is excluded.
const needed=ceil(7214n*S*S*S*early[1],10000n*3n*N*L2*critical);
console.log('sufficient multiplier for excluding the first mediant by itself',needed.toString(),'/',S.toString());

// Exhaustive small first-crossing control against the exact positional extremizer.
for(let j=1;j<=10;j++){
 const z=BigInt(j);
 const pmax=Number(ceil(al*z,S));
 assert(BigInt(pmax)===ceil(ah*z,S));
 const iopt=Array.from({length:j},(_,k)=>{
  const a=al*BigInt(k)/S,b=ah*BigInt(k)/S;assert(a===b);return Number(a);
 });
 const best=iopt.reduce((sum,i,k)=>sum+2n**BigInt(i)*3n**BigInt(j-k-1),0n);
 let seen=0,max=0n,winners=0;
 function go(i,m,A){
  if(i===pmax){
    if(m===j&&3n**BigInt(j)<2n**BigInt(i)){
      seen++;
      if(A>max){max=A;winners=1;}else if(A===max)winners++;
    }
    return;
  }
  if(3n**BigInt(m)<2n**BigInt(i))return;
  go(i+1,m,A);
  if(m<j)go(i+1,m+1,3n*A+2n**BigInt(i));
 }
 go(0,0,0n);
 assert(max===best&&winners===1);
}
console.log('position-extremizer exhaustive controls: j=1 through 10');
