// Run with Node.js. All comparisons are exact integers.
const assert = require("node:assert/strict");
function C(n,k) {
  if (k<0 || k>n) return 0n;
  let z=1n;
  for (let i=1;i<=k;i++) z=z*BigInt(n-i+1)/BigInt(i);
  return z;
}
function search(r) {
  const R=BigInt(r);
  let tested=0, maxN=0,maxH=0,min=null,zeros=[], ct=[];
  for(let s=r-1;s<5*r-11;s++){
    for(let h=1n+C(s,r-2);h<4n*C(s,r-3);h++){
      for(let n=s+1;C(n,r-1)<=h*C(s,r-1);n++){
        const m=C(n,r-1)-C(s,r-1),H=h-1n;
        if(m<H)continue;
        let l=r-1;
        while(l<s&&H*C(l+1,r-1)<=m)l++;
        const den=C(l,r-2);
        const bNum=den*(C(s,r-3)+H*C(l,r-3)-BigInt(r-3)*C(n,r-2))
          + C(l,r-4)*(m-H*C(l,r-1));
        const intPart=h*h+(2n*R-8n)*h+R*(6n-R);
        const num=den*intPart-4n*bNum;
        assert(num>=0n);
        if(n>s+1 && num<=(r===4?28n:60n)*den) {
          // In rank four on six points, 4-element hyperplanes
          // have disjoint complementary pairs, hence at most three.
          const impossible = r===4&&s===4&&n===6&&
            ((20n-h)%3n!==0n||(20n-h)/3n>3n);
          if(!impossible)ct.push([s,Number(h),n,String(num/den)]);
        }
        if(num===0n)zeros.push([s,Number(h),n]);
        else if(min===null || num*min[1]<den*min[0])min=[num,den];
        tested++; maxN=Math.max(maxN,n); maxH=Math.max(maxH,Number(h));
      }
    }
  }
  assert.deepEqual(zeros,[[r-1,r,r]]);
  assert.equal(tested,r===4?214:6738);
  assert.deepEqual(ct,r===4?[[3,10,5,"28"]]:[[4,15,6,"60"]]);
  console.log(JSON.stringify({rank:r,tested,maxN,maxH,zeros, connectedSharp:ct,
    ceilMinPositive:String((min[0]+min[1]-1n)/min[1])}));
}
search(4);
search(5);
