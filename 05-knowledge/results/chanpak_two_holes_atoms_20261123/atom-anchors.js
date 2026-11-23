
function H(r){const n=r.length,F=(1<<n)-1,d=new Uint32Array((F+1)*n);for(let i=0;i<n;i++)d[(1<<i)*n+i]=1;for(let s=1;s<=F;s++)for(let a=s;a;a&=a-1){let i=31-Math.clz32(a&-a),q=d[s*n+i];if(q)for(let b=r[i]&(F^s);b;b&=b-1){const t=b&-b,j=31-Math.clz32(t);d[(s|t)*n+j]+=q;}}let h=0;for(let i=0;i<n;i++)h+=d[F*n+i];return h;}

for(const [h,r] of [
 [49,[114,120,3,117,100,68,4]],
 [63,[238,236,248,240,163,192,16,64]],
 [343,[178,168,227,165,46,192,155,16]]
]){
 if(H(r)!==h) throw Error("height "+h);
 const n=r.length,F=(1<<n)-1;
 for(let i=0;i<n;i++)for(let j=0;j<n;j++)
  if(((r[i]>>j&1)+(r[j]>>i&1))!==Number(i!==j))throw Error("orientation");
 for(let i=0;i<n;i++){
  let reach=1<<i;
  for(let j=0;j<n;j++)for(let k=0;k<n;k++)if(reach>>k&1)reach|=r[k];
  if(reach!==F)throw Error("not strong");
 }
 console.log(h,n);
}
