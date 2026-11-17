function need(x,msg){if(!x)throw Error(msg);}
function pathTable(rows,root=-1){
  const n=rows.length,F=(1<<n)-1;
  for(let i=0;i<n;i++){need(!(rows[i]>>i&1),"loop");for(let j=0;j<i;j++)need((rows[i]>>j&1)+(rows[j]>>i&1)===1,"not tournament");}
  const d=new Float64Array((F+1)*n), h=new Float64Array(F+1);h[0]=1;
  for(let i=0;i<n;i++)if(root<0||root===i)d[(1<<i)*n+i]=1;
  for(let s=1;s<=F;s++)for(let i=0;i<n;i++){
    const v=d[s*n+i];if(!v)continue;h[s]+=v;
    for(let u=rows[i]&~s;u;u&=u-1){const b=u&-u,j=31-Math.clz32(b);d[(s|b)*n+j]+=v;}
  }
  return {d,h,F,n};
}
function inflate(D,q,B){
  const n=D.length,k=B.length;
  return [...D.map((a,i)=>a|(i===q?0:((1<<k)-1)<<n)),...B.map(a=>(a<<n)|(1<<q))];
}
function trans(n){return Array.from({length:n},(_,j)=>((1<<n)-1)^((1<<(j+1))-1));}
function ordinal(X,Y){return [...X.map(a=>a|(((1<<Y.length)-1)<<X.length)),...Y.map(b=>b<<X.length)];}
function full(rows){const H=pathTable(rows);return H.h[H.F];}
function ports(D,q){
  const s=pathTable(D),t=pathTable(D,q);
  const G=s.h[s.F],J=s.h[s.F^(1<<q)],e=s.d[s.F*s.n+q],L=G-e;
  let W=-J;
  for(let m=1;m<=s.F;m++)W+=t.h[m]*s.h[s.F^m];
  return [G,J,e,L,W];
}



function converse(a){return a.map((_,i)=>a.reduce((m,x,j)=>m|(((x>>i)&1)<<j),0));}
function led(R){
 const n=R.length, left=pathTable(R), right=pathTable(converse(R)),F=left.F;
 const w=Array(n).fill(0), C=Array.from({length:n},()=>Array(n).fill(0));
 for(let S=1;S<=F;S++)for(let y=0;y<n;y++){
   const z=right.d[S*n+y];if(!z)continue;w[y]+=z*left.h[F^S];
   for(let x=0;x<n;x++)C[x][y]+=z*left.d[(F^S)*n+x];
 }
 return [left.h[F],w,C];
}
function attach_q(R,Q){const n=R.length;return [...R.map((x,i)=>x|(Q>>i&1?0:1<<n)),Q];}
function strong(a){
 const F=(1<<a.length)-1;
 for(let i=0;i<a.length;i++){let r=1<<i;for(let j=0;j<a.length;j++)for(let k=0;k<a.length;k++)if(r>>k&1)r|=a[k];if(r!==F)return false;}
 return true;
}
function groups(A,limit){
 if(!A.length)return [[]];
 const ans=[];
 function run(prefix,vs){if(!vs.length){ans.push(prefix);return;}
 for(let t=0;t<=limit;t++)if(prefix[1]+t<=limit)run([ordinal(ordinal(prefix[0],trans(t)),vs[0]),prefix[1]+t],vs.slice(1));}
 run([[],0],A);
 let result=[];for(const [R,m] of ans)for(let k=0;m+k<=limit;k++)result.push(ordinal(R,trans(k)));
 return result;
}
function scan(R,g,K,out){
 const [h,w,C]=led(R);need(h===g,"bank");
 let F=(1<<R.length)-1;
 for(let Q=1;Q<=F;Q++){
  let W=0,U=0;
  for(let y=0;y<R.length;y++)if(Q>>y&1){W+=w[y];for(let x=0;x<R.length;x++)if(Q>>x&1)U+=C[x][y];}
  if(U!==g-K||W%7)continue;
  while(!(W%2)&&W)W>>=1;
  if(!out[W])out[W]={R,Q};
 }
}
const triangle=[2,4,1];

function all(n,fun){for(let code=0;code<1<<(n*(n-1)/2);code++){
 let a=Array(n).fill(0),k=0;for(let i=0;i<n;i++)for(let j=i+1;j<n;j++,k++)if(code>>k&1)a[i]|=1<<j;else a[j]|=1<<i;fun(a);
}}
function key(a){let best=null;const p=a.map((_,i)=>i);
 function go(i){if(i===p.length){let x=0,k=0;for(let i=0;i<p.length;i++)for(let j=i+1;j<p.length;j++,k++)if(a[p[i]]>>p[j]&1)x|=1<<k;if(best===null||x<best)best=x;return;}
 for(let j=i;j<p.length;j++){[p[i],p[j]]=[p[j],p[i]];go(i+1);[p[i],p[j]]=[p[j],p[i]];}}
 go(0);return a.length+":"+best;
}

const cores={3:{},5:{},9:{},11:{},13:{},15:{}};
for(let n=3;n<=6;n++)all(n,a=>{const h=full(a);if(cores[h]&&strong(a))cores[h][key(a)]=a;});
need(JSON.stringify(Object.values(cores).map(x=>Object.keys(x).length))==="[1,1,2,1,1,4]","core catalogue");
const out={};
for(const g of [1,3,5,9,11,13,15]){
 const seqs=g===1?[[]]:Object.values(cores[g]).map(a=>[a]);
 if(g===9)seqs.push([triangle,triangle]);
 if(g===15){const five=Object.values(cores[5])[0];seqs.push([triangle,five],[five,triangle]);}
 out[g]={};
 for(const K of [-7,7]){
  if(g<K)continue;
  const allowed={},limit=2+Math.floor(Math.log2(g-K));
  for(const seq of seqs){
   const bank=seq.length?groups(seq,limit):Array.from({length:limit},(_,i)=>trans(i+1));
   for(const r of bank)scan(r,g,K,allowed);
  }
  out[g][K]=Object.keys(allowed).map(Number).sort((a,b)=>a-b);
 }
}
const expected={"1":{"-7":[]},"3":{"-7":[]},"5":{"-7":[]},"9":{"7":[],"-7":[]},"11":{"7":[],"-7":[]},"13":{"7":[],"-7":[49,119]},"15":{"7":[21],"-7":[]}};
need(JSON.stringify(out)===JSON.stringify(expected),"bounded deletion passport");
for(const [R,Q,g,W,K,N] of [[[44,53,56,50,33,0],36,13,49,-7,91],
                            [[62,24,42,48,36,2],3,13,119,-7,231],
                            [[6,12,24,17,3],5,15,42,7,91]]){
 const D=attach_q(R,Q),[G,J,e,L,T]=ports(D,R.length);
 need(J===g&&T===W&&L+J-T===K,"literal profile");
 const child=inflate(D,R.length,trans(1));
 need(full(child)===N&&strong(child),"small sharpness case");
}
// Independent comparison between the contact ledger and the rooted profiles.
for(let n=1;n<=4;n++)all(n,R=>{
 const [g,w,C]=led(R);
 for(let Q=0;Q<1<<n;Q++){
  let W=0,U=0;
  for(let y=0;y<n;y++)if(Q>>y&1){W+=w[y];for(let x=0;x<n;x++)if(Q>>x&1)U+=C[x][y];}
  const [G,J,e,L,W0]=ports(attach_q(R,Q),n);
  need(J===g&&W===W0&&U===W-L,"ledger identity");
 }
});
globalThis.console?.log(JSON.stringify(out));
