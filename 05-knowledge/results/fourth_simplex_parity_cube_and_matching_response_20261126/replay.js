"use strict";
const eq=(x,y)=>{if(JSON.stringify(x)!==JSON.stringify(y)) throw Error(JSON.stringify([x,y]));};
const pc=x=>{let n=0;for(;x;x&=x-1)n++;return n;};
function bitrows(n,edge){return Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i!==j&&edge(i,j)?1<<j:0).reduce((a,b)=>a|b,0));}
function hs(out){
 const n=out.length,K=1<<n, H=Array(K).fill(0), dp=Array.from({length:K},()=>Array(n).fill(0));
 for(let i=0;i<n;i++)dp[1<<i][i]=1;
 for(let m=1;m<K;m++)for(let i=0;i<n;i++)if(dp[m][i]){
   H[m]+=dp[m][i];
   for(let rem=out[i]&~m;rem;rem&=rem-1){const b=rem&-rem,j=31-Math.clz32(b);dp[m|b][j]+=dp[m][i];}
 }
 return H;
}
function cover(out){
 const H=hs(out),K=H.length,n=out.length, dp=Array.from({length:K},()=>Array(n+1).fill(0));
 dp[0][0]=1;
 for(let m=1;m<K;m++)for(let u=m;u;u=(u-1)&m)for(let k=1;k<=pc(m);k++)dp[m][k]+=H[u]*dp[m^u][k-1];
 return dp[K-1];
}
function strong(o){
 for(let s=0;s<o.length;s++){let m=1<<s,old;do{old=m;for(let i=0;i<o.length;i++)if(m>>i&1)m|=o[i];}while(old!==m);if(m!==(1<<o.length)-1)return false;}
 return true;
}
function sample(mask,k,cyclic=false){
 return bitrows(8+k,(i,j)=>{
  if(i>=8&&j>=8)return cyclic?(j-i+3)%3===1:i<j;
  if(i>=8)return j<4;
  if(j>=8)return i>=4;
  if((i<4)!==(j<4))return (i<4)===((i^j)%2===0);
  if((j-i+4)%4===1)return true;
  if((i-j+4)%4===1)return false;
  return (i%4<2)!==Boolean(mask>>(2*(i>=4)+(i%2))&1);
 });
}
function runs(out){
 const c=Array(6).fill(0),s=8;
 const arc=(i,j)=>Boolean(out[i]>>j&1);
 function rec(p,used){
   if(p.length<8){for(let i=0;i<8;i++)if(!(used>>i&1))rec([...p,i],used|(1<<i));return;}
   let r=0,f=0;
   for(let i=0;i<=8;i++){
     const ok=(!i||arc(p[i-1],s))&&(i===8||arc(s,p[i]));
     const need=i>0&&i<8&&!arc(p[i-1],p[i]);
     if(need&&!ok)return;
     if(need)r++;else if(ok)f++;
   }
   let coeff=1;
   for(let j=0;j<=f;j++){c[r+j]+=coeff;coeff=coeff*(f-j)/(j+1);}
 }
 rec([],0);return c.slice(1);
}
eq(runs(sample(0,1)),[2917,4638,3274,984,96]);
eq(runs(sample(1,1)),[2919,4639,3274,984,96]);
const common=[38622,232422,821040,1857000,2751840,2610720,1451520,362880];
for(let m=0;m<16;m++){
 const a=sample(m,1);
 eq(strong(a),true);eq(cover(a),[0,2917+2*(pc(m)%2),...common]);
}
function simplex(x){return x*(x-1)*(x-2)*(x-3)/24;}
for(const [k,cyclic,t,d,common,L,g] of [
 [1,false,7,[1,1,1,1],[4,4,9,10],2917,2],
 [3,true,2,[1,1,1,6],[4,30],56223,12],
 [4,false,7,[1,2,2,2],[24,40],209329,16]
]){
 const sides=[[...common],[...common]];
 for(let m=0;m<16;m++)sides[pc(m)%2].push(d.reduce((a,x,i)=>a+((m>>i&1)?x:0),t));
 const sum=(y,p)=>y.reduce((a,b)=>a+b**p,0);
 for(let p=0;p<4;p++)eq(sum(sides[0],p),sum(sides[1],p));
 const f=sides.map(y=>y.reduce((a,b)=>a+simplex(b),0)).sort((a,b)=>a-b);
 eq(f,[(L-1)/2,(L+g-1)/2]);
 for(let m=0;m<16;m++){const o=sample(m,k,cyclic);eq(strong(o),true);eq(hs(o).at(-1),L+g*(pc(m)%2));}
 console.log("literal 16-state cube:",L,L+g,"simplex:",f.join(","));
}
function mutate(o,e,mask){
 const a=o.slice();
 e.forEach(([i,j],k)=>{if(mask>>k&1){a[i]^=1<<j;a[j]^=1<<i;}});return a;
}
function contrast(o,e,stat){
 let val=0;for(let m=0;m<(1<<e.length);m++)val+=(-1)**pc(m)*stat(mutate(o,e,m));return val;
}
function principal(o,e,S,u){
 const vs=S.flatMap(i=>e[i]).concat(u);
 const ee=S.map((_,i)=>[2*i,2*i+1]);
 return contrast(bitrows(vs.length,(i,j)=>o[vs[i]]>>vs[j]&1),ee,b=>hs(b).at(-1));
}
function partitions(n){
 const all=[];
 function rec(i,p){if(i===n){all.push(p);return;}
  rec(i+1,[...p,[i]]);
  for(let j=0;j<p.length;j++)rec(i+1,p.map((a,k)=>k===j?[...a,i]:a));
 }
 rec(0,[]);return all;
}
const fac=j=>j<=1?1:j*fac(j-1);
let state=81411;
const random=()=>{state=(Math.imul(state,1664525)+1013904223)>>>0;return state>>>31;};
for(const [n,r] of [[4,2],[5,2],[6,2],[7,2],[8,2],[7,3],[8,3],[9,3],[10,3]]){
 let nonzero=0;
 for(let it=0;it<12;it++){
   const edges=Array.from({length:r},(_,i)=>[2*i,2*i+1]);
   const out=Array(n).fill(0);
   for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)if(random())out[i]|=1<<j;else out[j]|=1<<i;
   const d=n-2*r;
   const deltas=Array(n+1).fill(0);
   for(let m=0;m<1<<r;m++)cover(mutate(out,edges,m)).forEach((a,i)=>deltas[i]+=(-1)**pc(m)*a);
   eq(deltas.slice(d+1),Array(n-d).fill(0));
   let val=0;
   for(const p of partitions(r)){
    function rec(j,used,v){if(j===p.length){val+=v;return;}
      for(let u=2*r;u<n;u++)if(!(used>>u&1))rec(j+1,used|(1<<u),v*principal(out,edges,p[j],u));
    }
    rec(0,0,1);
   }
   eq(deltas[d],val*fac(d));
   if(val) nonzero++;
 }
 console.log("top-layer link formula",n,r,"nonzero controls",nonzero,"of 12");
}
