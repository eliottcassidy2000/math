const assert=require('node:assert/strict');
function sq(n,sum){
 if(n===0)return sum===0?0:Infinity;
 if(sum<0)return Infinity;
 const q=Math.floor(sum/n),r=sum-q*n;
 return n*q*q+r*(2*q+1);
}
function cases(f,g,enhanced){
 const M=Math.max(f,g),u=2+f+g,t=M+2,R=f*g,A=f*(56-g),B=(56-f)*g,C=(56-f)*(56-g);
 let count=0,minX=null,check=0, witness=null;
 const ms=[];
 if(enhanced)for(let m=0;m<=Math.min((Math.min(f,g)-1)*R,Math.floor((50*R*R+7*3250*R)/3250));m+=2){
  if(m*3250<65*R*R-8*3250*R)continue;
  let lhs=sq(R,m)+sq(A,(f-1)*R-m)+sq(B,(g-1)*R-m);
  const sm=(57-f-g)*R+m,budget=(56-f-g)*R+R*R-m-lhs;
  if(budget>=sq(C,sm))ms.push([m,sm,budget]);
 }
 for(let x=t*t+1;x<=C;x++){
  const l=Math.max(f+g+M-53,u-Math.floor((x-1)/t));
  const z=6500*l-115*R;
  if(z>0&&z*z*x>225*R*(3250-x)*(3250-R))continue;
  check++;
  if(enhanced){
   let found=null;
   for(let [m,sm,budget] of ms){
    if(sm<x*l)continue;
    const lower=sm>=C*l?sq(C,sm):x*l*l+sq(C-x,sm-x*l);
    if(lower<=budget){found=m;break}
   }
   if(found===null)continue;
   if(witness===null)witness=[x,found,l];
  }
  count++;if(minX===null)minX=x;
 }
 return {f,g,count,minX,check,witness};
}
const frontier=[[22,19,490],[23,16,545],[24,13,596],[25,9,674],[26,5,578]];
for(const [f,g,k] of frontier){
 const r=cases(f,g,true);
 assert.equal(r.count,0);assert.equal(r.check,k);
 assert.equal(cases(g,f,true).count,0);
 assert(cases(f,g-1,true).count>0);
}
console.log('rectangle frontier:',frontier);
function q(n){
 let r=0n;while(n>0n){const a=n>>1n;if(n&1n)r+=(a+1n)>>1n;n=a}
 return r;
}
function H(v){return 1662n*v-58n*q(v)+427n}
function reverse(u,k){let v=0n;for(let i=0;i<k;i++){v=(v<<1n)|(u&1n);u>>=1n}return v;}
function inverse(h, prefix='1010') {
 const low=h>=427n?(h-427n+1661n)/1662n:0n,high=(h-427n)/1604n;
 if(high<low||high<0n)return [];
 const out=[];
 for(let L=Math.max(1,low.toString(2).length);L<=high.toString(2).length;L++){
  if(L<prefix.length)continue;
  const t=Math.max(0,Math.floor((L-7)/2)), pow=1n<<BigInt(t);
  const first=1n<<BigInt(L-1), last=(first<<1n)-1n;
  const lo=low>first?low:first,hi=high<last?high:last;
  if(lo>hi)continue;
  for(let u=lo/pow;u<=hi/pow;u++){
   const base=u*pow;
   if(prefix&&!base.toString(2).startsWith(prefix))continue;
   let a=0n,b=pow;
   while(a<b){const m=(a+b)>>1n;if(H(base+reverse(m,t))>h)a=m+1n;else b=m}
   if(a<pow){const v=base+reverse(a,t);if(v>=lo&&v<=hi&&H(v)===h)out.push(v)}
  }
 }
 return out;
}
let x=inverse(546598033n).map(v=>2n*v);
assert(x.includes(659056n)&&x.includes(666074n));
for(let v=160n;v<2560n;v++)if(v.toString(2).startsWith('1010'))assert(inverse(H(v)).includes(v));
console.log('explicit collision preimages:', x.map(String));
function step(state,b,M){return [(2*state[0]+b)%(2*M),(state[1]+b*Math.ceil(state[0]/2))%M]}
function run(w,M,state=[0,0]){for(let b of w)state=step(state,Number(b),M);return state}
for(let M=1;M<=128;M++){
 const stack=[[0,0]],seen=new Set(['0,0']),ind=new Int32Array(2*M*M);
 for(let i=0;i<stack.length;i++)
  for(let b=0;b<=1;b++){
   const r=step(stack[i],b,M),s=r.join(',');
   ind[r[0]*M+r[1]]++;
   if(!seen.has(s)){seen.add(s);stack.push(r)}
  }
 assert.equal(seen.size,2*M*M);assert(ind.every(a=>a===2))
}
for(let a=1;a<9;a++){
 const M=1<<a;
 for(let r=0;r<M;r++){
  let w='';
  for(let i=0;i<a;i++)if(r>>i&1)w+='1'+'0'.repeat(i+1)+'1'+'0'.repeat(a+1);
  assert.deepEqual(run(w,M),[0,r])
 }
}
console.log('mixed residue state checks through M=128; dyadic routing through 256');