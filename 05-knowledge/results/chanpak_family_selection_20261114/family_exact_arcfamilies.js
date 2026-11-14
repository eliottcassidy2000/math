function prof(rows){
  const n=rows.length,T=1<<n,F=T-1;
  const dp=new Float64Array(n*T),H=new Float64Array(T); H[0]=1;
  for (let i=0;i<n;i++)dp[n*(1<<i)+i]=1;
  for(let s=1;s<T;s++){
   let z=s;
   while(z){let b=z&-z,i=31-Math.clz32(b); z^=b; let cur=dp[n*s+i];let t=rows[i]&~s;
     while(t){let b=t&-t,j=31-Math.clz32(b);t^=b;dp[n*(s|b)+j]+=cur;}
     H[s]+=cur;
   }
  }
  const F2=new Float64Array(T);
  for(let s=0;s<T;s++)for(let u=s;;u=(u-1)&s){F2[s]+=H[u]*H[s^u];if(!u)break}
  let f3=0;for(let u=0;u<T;u++)f3+=H[u]*F2[F^u];
  return [H[F],F2[F],f3]
}
function seed(a,b,c){
 let R=Array(a+b+c).fill(0),g=Array.from({length:a+b+c},(_,i)=>i<a?0:i<a+b?1:2);
 for(let i=0;i<R.length;i++)for(let j=0;j<R.length;j++)if((i<j&&g[i]===g[j])||(g[j]-g[i]+3)%3===1)R[i]|=1<<j;
 return R;
}
for(let n=1;n<=11;n++)for(let j=1;j<=n;j++)for(let z=0;z<2;z++){
 const R=seed(1,2,n),a=j-1,b=n-j;
 R[1+z]^=1<<(2+j);R[2+j]^=1<<(1+z);
 const d=z ? 2**b*(2*3**a-3) - Number(b===0) : 2**(b+1)*(3**a-1)+2**a+Number(a===0)-2*Number(b===0);
 const N=3*2**n-1+d;
 if(prof(R)[0]!==N)throw Error([n,j,z,N,prof(R)[0]]);
}
console.log('PASS finite exact arc representatives through n=11')
