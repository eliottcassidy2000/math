
function Q(m,a){const j=Math.floor(m/a);return a*j*j+(m%a)*(2*j+1);}
function gap(p,s,caps){
 const m=p*s;
 return 56*caps.reduce((sum,[n,a])=>sum+n*Q(m,a),0)
  -caps.reduce((n,x)=>n+x[0],0)*m*m-(8-s)*(s+7)*m*(56-p);
}
function lanes(t){
 const out=[];
 for(let h=0;h<=2*t;h++)for(let s=Math.max(1,54-h-Math.floor(t/2));s<=7;s++){
   const L=54-h-s, B=Math.max(0,L-Math.floor((15+2*t-h)/2));
   const caps=[[B,14],[L-B,28],[57-s-L,56]];
   const ps=[];
   for(let p=2;p<=(B?14:28);p+=2)if(p>=s-1&&gap(p,s,caps)<=0)ps.push(p);
   if(ps.length)out.push({h,s,B,p:ps});
 }
 return out;
}

function spineCheck(){
 const F=(u)=>{let v=1;while(u){const q=u>>1;if(!(u&1))v+=q+(q&1);u=q;}return v;};
 const H=(b,c)=>b*c+F(b-1)*F(c-1)-(b&1)*(c&1);
 function graph(b,c,start,t,flip){
  const l=Math.max(start+t,32-Math.clz32(b-1)),r=Math.max(1,32-Math.clz32(c-1));
  const x=l,y=l+r+1, n=y+1, A=new Uint32Array(n);
  const add=(i,j)=>A[i]|=1<<j;
  for(let i=0;i<l;i++){
   for(let j=i+1;j<l;j++){
    if(flip&&i===start&&j===start+t-1)add(j,i);else add(i,j);
   }
   if(((b-1)>>i)&1)add(i,x);else add(x,i);
  }
  for(let i=x+1;i<y;i++){
   for(let j=i+1;j<y;j++)add(j,i);
   if(((c-1)>>(i-x-1))&1)add(y,i);else add(i,y);
  }
  for(let i=0;i<=x;i++)for(let j=x+1;j<=y;j++)if(i===x&&j===y)add(i,j);else add(j,i);
  return A;
 }
 function dp(A){
  const n=A.length,M=(1<<n)-1,D=new Float64Array((M+1)*n);
  for(let i=0;i<n;i++)D[(1<<i)*n+i]=1;
  for(let s=1;s<M;s++)for(let i=0;i<n;i++){
   const z=D[s*n+i];if(!z)continue;
   for(let m=A[i]&~s;m;m&=m-1){const j=m&-m;D[(s|j)*n+31-Math.clz32(j)]+=z;}
  }
  let k=0;for(let i=0;i<n;i++)k+=D[M*n+i];return k;
 }
 let tested=0, positive=0, unit=0,twoRuns=0;
 for(let t=3;t<=6;t++)for(let start=0;start<=2;start++)for(let above=0;above<4;above++)for(let type=0;type<2;type++)for(let c=1;c<=7;c+=2){
  const below=(1<<start)-1;
  const b=(above<<(start+t))+(type*((1<<t)-1)<<start)+below+1;
  const a=k=>(above<<(start+k))+(type*((1<<k)-1)<<start)+below+1;
  if(above===0&&type===0&&c>3)continue; // discard redundant high-zero padding controls
  const H1=H(a(1),c),H2=H(a(2),c);
  const d=2**(t-2)*H1+(3**(t-2)-1-2**(t-2))*(H2-H1);
  if(H2>H1)twoRuns++;
  const old=dp(graph(b,c,start,t,false)),now=dp(graph(b,c,start,t,true));
  if(old!==H(b,c)||now!==old+d)throw Error("module identity");
  tested++;positive+=d>0;if(d===2){unit++;if(old!==1)throw Error("unexpected unit");}
 }
 return {tested,twoRuns,positive,unit};
}

const r={deltas:Array.from({length:9},(_,i)=>gap(i+6,7,[[2,14],[7,28],[38,56]])),at127:lanes(19),at128:lanes(20)};
if(r.deltas.some(a=>a<=0)||r.at127.length!==0||r.at128.length!==10)throw Error("closure check");
const sp=spineCheck();
if(sp.tested!==360||sp.twoRuns!==336||sp.unit!==1)throw Error("spine check");
console.log(JSON.stringify(r));
console.log(JSON.stringify(sp));
