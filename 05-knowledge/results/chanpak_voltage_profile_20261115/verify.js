// Run with Node. Check digit-split drift and constant phase separation.
const need=(ok,z)=>{if(!ok)throw Error(JSON.stringify(z));};
const bits=(x,r)=>x.toString(2).padStart(r,"0");
const rev=(x,r)=>parseInt([...bits(x,r)].reverse().join(""),2);
function p(x){let y=0;while(x){let u=Math.floor(x/2);if(x%2===0)y+=u+u%2;x=u;}return y;}
function gcd(a,b){while(b){[a,b]=[b,a%b];}return Math.abs(a);}
let orders=0, drifts=0;
for(let m=1;m<=256;m++){
 let a=0;while(m%2**(a+1)===0)a++;
 const r=2+Math.floor(Math.log2(2*m)), tests=[];
 for(const h of [r,r+1]){
  const b=2**(h-1)-1;
  const d=c=>rev(c+2*m,h)-rev(c,h);
  need(d(b+1)-d(b)===3*(2**(h-a-2)-1),[m,h]);
  for(const c of [b,b+1])tests.push([h,c,d(c)]);
 }
 need(gcd(tests[1][2]-tests[0][2],tests[3][2]-tests[2][2])===3,m);
 for(let q=2;q<=m;q++)if(m%q===0 && q!==3){
  const z=tests.find(([, ,d])=>d%q!==0);
  need(z,[m,q]);
  const [h,c,d]=z;
  let min=1, sumCos=0,sumSin=0;
  for(let v=0;v<m;v++){
   const i=2**h*2*v+c;
   const diff=(p(i+2*m)-p(i))/2;
   need(diff===(p(c+2*m)-p(c))/2-v*d,[m,q,v,diff]);
   const phase=2*Math.PI*((diff%q)/q);
   min=Math.min(min,Math.cos(phase));
   sumCos+=Math.cos(phase);sumSin+=Math.sin(phase);
   drifts++;
  }
  need(min<1e-10 && Math.abs(sumCos)+Math.abs(sumSin)<1e-9,[m,q,min]);
  orders++;
 }
}
console.log({orders,drifts});
