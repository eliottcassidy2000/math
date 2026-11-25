const A=require('node:assert/strict');
const p=61n,E=(1n<<(p-1n))*((1n<<p)-1n);
for(const [h,a,k,expected] of [[29n,640n,111n,[[37,29],[37,30],[37,29]]],[651n,17616n,2313n,[[36,27],[36,27],[36,27]]]]){
 let pairs=[];
 for(const d of [-3n,-1n,1n]){
  let first,balanced;
  for(let r=1;r<=59;r++){
   const T=1n<<BigInt(r);
   let inv=1n;
   for(let i=0;(1n<<BigInt(i))<T;i++)inv=inv*(2n-h*inv)%T;
   const Q=((inv*((k-d)/2n))%T+T)%T,H=E+d;
   A.equal(((h*Q-(k-d)/2n)%T+T)%T,0n);
   if(first===undefined&&2n*a*T*T*(2n*Q-1n)>H)first=r;
   if(balanced===undefined&&16n*a*T*T*Q*Q>9n*H)balanced=r;
  }
  pairs.push([first,balanced]);
 }
 A.deepEqual(pairs,expected);console.log(String(h),pairs);
}
const profiles=new Set();
for(let i=0;i<=6;i++)for(let j=0;j<=6-i;j++){
 const k=6-i-j;
 if((i-j)%3===0&&(i-k)%3===0)profiles.add([i,j,k].sort((a,b)=>a-b).join(','));
}
A.deepEqual([...profiles].sort(),['0,0,6','0,3,3','1,1,4','2,2,2']);
console.log([...profiles].sort());