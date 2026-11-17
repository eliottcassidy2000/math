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

const D=[270,140,264,272,167,463,31,333,82], q=8;
const ordinary=[270,268,264,256,239,463,271,335,16];
const order=[8,4,5,7,6,0,1,2,3];
const flips=[[1,7],[1,8],[3,4],[6,8],[4,6]];
const coefficients=[[128,1],[256,5],[276,-9],[475,109],[823,13],[896,7]];
let core=[...ordinary];
for(let j=0;j<coefficients.length;j++){
  need(order.slice(1).every((v,i)=>core[order[i]]>>v&1),"lost shared frame");
  const [G,J,e,L,W]=ports(core,q);
  need(W===coefficients[j][0]&&L+J-W===coefficients[j][1],"incorrect profile");
  if(j===flips.length)break;
  const [u,v]=flips[j];
  core[u]^=1<<v;core[v]^=1<<u;
}
need(core.toString()===D.toString()&&ports(D,q).toString()==="877,127,101,776,896","terminal profile");

function replace(rows,slot,B){
  const parent=rows.flatMap((_,i)=>i===slot?B.map((_,j)=>[slot,j]):[[i,0]]);
  return parent.map(([i,u])=>parent.reduce((a,[j,v],l)=>a|
      (((i===slot&&j===slot?B[u]>>v:rows[i]>>j)&1)<<l),0));
}
function apex(w){
  const n=w.length, A=trans(n);
  let q=0;
  for(let i=0;i<n;i++)if(w[i]==="1")q|=1<<i;else A[i]|=1<<n;
  return [...A,q];
}
const cycle=[2,4,1];
const small=[
 [35,[110,108,120,48,99,64,8]],
 [63,apex("1001100")],
 [119,replace([250,228,249,242,194,16,32,96],5,trans(2))],
 [231,ordinal([74,32,3,6,15,29,62],cycle)],
 [455,ordinal([74,32,3,22,7,29,62],[8,1,3,6])],
 [903,ordinal([52,5,136,227,110,6,39,115],cycle)],
];
for(const [want,T] of small)need(full(T)===want,"small plus case "+want);
for(let m=1;m<=6;m++)
  need(full(inflate(D,q,trans(m)))===7*(2**(m+7)+1),"plus rail "+m);

function negative(B){
  const left=replace([4,5,72,19,39,79,155,63],0,ordinal(B,[0]));
  const n=left.length, v=B.length+5, right=trans(3).map((r,i)=>(r<<n)|((1<<n)-1));
  left[v]|=1<<n;right[0]^=1<<v;
  return [...left,...right];
}
for(const [B,h,P] of [[trans(1),1,2],[trans(3),1,8],[cycle,3,12],[trans(4),1,16]]){
  const plus=full(inflate(D,q,B)), minus=full(negative(B));
  need(plus===7*(128*P+h)&&minus===7*(128*P-h),"signed window");
  need((plus/7)%3===0||(minus/7)%3===0,"mod3 guard");
}

// Independent port profile check for c=10 P in the old negative scaffold.
for(let m=1;m<=2;m++){
  const rows=replace([4,5,72,19,39,79,155,63],0,trans(m));
  const root=m+4;
  const s=pathTable(rows),t=pathTable(rows,root);
  let term=0;
  for(let M=1;M<=s.F;M++)term+=s.d[M*s.n+root]*s.h[s.F^M];
  need(s.h[s.F]===64*2**m+1&&term===96*2**m-2&&t.h[t.F]===10*2**m,"negative triple");
}
globalThis.console?.log("Six small values, six positive heights, five monotone edits, and signed profiles: PASS");
