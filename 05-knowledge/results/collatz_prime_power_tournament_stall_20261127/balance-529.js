function balance529(){
 const m=529,g=399,bs=23,N=1<<bs, coeff=new Float64Array(N), used=new Int16Array(m), reps=[];
 const Q=new Set(Array.from({length:22},(_,i)=>(i+1)**2%23));
 for(let i=1;i<23;i++)used[i*23]=Q.has(i)?24:-24;
 for(let a=1;a<m;a++)if(!used[a]){
   reps.push(a);const idx=reps.length;
   for(let x=a;;){used[x]=idx;used[m-x]=-idx;x=x*g%m;if(x===a)break;}
 }
 if(reps.length!==bs)throw Error("orbit count");
 reps.push(23);
 for(const t of reps){
   const poly=new Map([[0,-1]]);
   for(let x=1;x<m;x++){
     const y=(t-x+m)%m;if(y===0)continue;
     const a=used[x],b=used[y], mask=((1<<(Math.abs(a)-1))^(1<<(Math.abs(b)-1)))&(N-1);
     poly.set(mask,(poly.get(mask)||0)+Math.sign(a)*Math.sign(b));
   }
   for(const [a,v] of poly)if(v)for(const [b,w] of poly)if(w)coeff[a^b]+=22*v*w;
 }
 for(let bit=1;bit<N;bit*=2)for(let i=0;i<N;i+=2*bit)for(let j=i;j<i+bit;j++){const a=coeff[j],b=coeff[j+bit];coeff[j]=a+b;coeff[j+bit]=a-b;}
 let best=Infinity,count=0,mask=0,next=Infinity;
 for(let i=0;i<N;i++){
   const v=coeff[i];
   if(v<best){next=best;best=v;count=0;mask=i}
   if(v===best)count++;
   else if(v<next)next=v;
 }
 
 function chosen(mask,i){const v=used[i];return Math.sign(v)===(mask>>(Math.abs(v)-1)&1?-1:1);}
 function canonical(mask){
   let min=mask;
   for(let k=1;k<23;k++){
     let a=0;for(let j=0;j<bs;j++)if(!chosen(mask,reps[j]*(1+23*k)%m))a|=1<<j;
     min=Math.min(min,a);
   }
   return min;
 }
 const canon=new Set();
 for(let i=0;i<N;i++)if(coeff[i]===best)canon.add(canonical(i));
 if(best!==89408||count!==46||canon.size!==2)throw Error("minimizer census");
 const D=new Set();for(let i=1;i<m;i++)if(chosen(mask,i))D.add(i);
 const H={};
 for(let t=1;t<m;t++){let r=D.has(t)?1:0;for(const x of D)if(D.has((t-x+m)%m))r++;H[r]=(H[r]||0)+1}
 if(Object.entries(H).reduce((s,[v,n])=>s+(4*(+v)-530)**2*n,0)!==best)throw Error("direct balance replay");
 const square=new Set(Array.from({length:22},(_,i)=>(i+1)**2%23));
 let prod=0,flip=0;
 for(let j=0;j<bs;j++)if(!square.has(reps[j]%23))prod|=1<<j;
 flip=prod^1;
 console.log({best:best/16, count, mask, next:next/16, hist:H, products:coeff[prod]/16, flipped:coeff[flip]/16, canon:[...canon], reps});
}
balance529();
