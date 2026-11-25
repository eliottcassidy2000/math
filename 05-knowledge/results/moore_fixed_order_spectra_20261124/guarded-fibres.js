function census(n){
const N=1<<n, k=N-1; const tiles=[];
for(let i=0;i<n;i++)for(let j=i+2;j<n;j++)tiles.push([i,j]);
const work=new Int32Array(N*n*2), p=new Int32Array(N), end=new Int32Array(N), start=new Int32Array(N), outs=new Int32Array(n), inns=new Int32Array(n);
const data=new Int32Array((1<<tiles.length)*N*5);
for(let m=0;m<(1<<tiles.length);m++){
 outs.fill(0);inns.fill(0);work.fill(0);p.fill(0);p[0]=1;
 for(let i=0;i<n-1;i++)outs[i]|=1<<(i+1);
 tiles.forEach(([i,j],b)=>{if((m>>b)&1)outs[j]|=1<<i;else outs[i]|=1<<j});
 for(let i=0;i<n;i++){for(let j=0;j<n;j++)if((outs[j]>>i)&1)inns[i]|=1<<j;work[(1<<i)*2*n+i]=work[(1<<i)*2*n+n+i]=1;}
 for(let a=1;a<N;a++)for(let e=a;e;e&=e-1) {
  let i=31-Math.clz32(e&-e), j0=a^(1<<i), E=0,S=0;
  if(!j0){p[a]=1;break;}
  for(let b=inns[i]&j0;b;b&=b-1)E+=work[j0*2*n+31-Math.clz32(b&-b)];
  for(let b=outs[i]&j0;b;b&=b-1)S+=work[j0*2*n+n+31-Math.clz32(b&-b)];
  work[a*2*n+i]=E;work[a*2*n+n+i]=S;p[a]+=E;
 }
 let Q=0;for(let a=0;a<N;a++)Q+=p[a]*p[k^a];
 for(let o=0;o<N;o++){
  end[0]=start[0]=1;
  for(let a=1;a<N;a++){
   let E=0,S=0;
   for(let b=a&~o;b;b&=b-1)E+=work[a*2*n+31-Math.clz32(b&-b)];
   for(let b=a&o;b;b&=b-1)S+=work[a*2*n+n+31-Math.clz32(b&-b)];
   end[a]=E;start[a]=S;
  }
  let H=0,R=0,J=0;
  for(let a=0;a<N;a++){let b=k^a;H+=end[a]*start[b];R+=p[a]*start[b];if(a)J+=end[a]*p[b];}
  let z=(m*N+o)*5;data[z]=p[k];data[z+1]=Q;data[z+2]=R;data[z+3]=J;data[z+4]=H;
 }
}
let edges=0;const best=[];
for(let o=0;o<N;o++){
 const up=new Map(), length=new Map(),prev=new Map();
 for(let m=0;m<(1<<tiles.length);m++){
  let z=(m*N+o)*5;
  let H=data[z+4];
  for(let b=0;b<tiles.length;b++){
   let t=m^(1<<b), w=(t*N+o)*5;
   if(data[w+4]!=H+2)continue;
   if([0,1,2,3].some(x=>data[z+x]!=data[w+x]))continue;
   if(!up.has(m))up.set(m,[]);
   up.get(m).push(t);edges++;
  }
 }
 const entries=[...up.keys()].sort((a,b)=>data[(a*N+o)*5+4]-data[(b*N+o)*5+4]);
 for(let m of entries){
  let x=length.get(m)??1;
  for(let t of up.get(m)){
   if((length.get(t)??1)<x+1){length.set(t,x+1);prev.set(t,m)}
   if(x+1>=(best[0]?.len??2)){
    if(x+1>(best[0]?.len??2))best.length=0;
    const path=[t];let w=t;while(prev.has(w)){w=prev.get(w);path.push(w);}path.reverse();
    const z=(t*N+o)*5;
    best.push({len:x+1,o,path,vals:path.map(q=>data[(q*N+o)*5+4]),signature:[...data.slice(z,z+4)]})
   }
  }
 }
}
return {n,edges,longest:best[0]?.len,count:best.length,guardHist:n<=6?best.reduce((a,b)=>{const k=b.signature.join(',')+';'+b.vals.join(',');a[k]=(a[k]??0)+1;return a},{}):undefined}
};

const assert=require('node:assert/strict');
const upto=Number(process.argv[2]??6);
for(let n=2;n<=upto;n++){
 const r=census(n);
 assert.equal(r.edges,n<6?0:n==6?48:96);
 if(n>=6){assert.equal(r.longest,2);assert.equal(r.count,r.edges)}
 if(n==6){
   assert.deepEqual(r.guardHist,{'29,320,158,162;111,113':24,'29,320,191,129;111,113':24})
 }
 console.log(r);
}
function h(v){if(!v)return 1;const q=v>>1;return h(q)+(v&1?0:q+(q&1))}
function H(v){return 802*v+29*(h(v)-1)+427}
if(H(659056)!==546598033||H(666074)!==546598033)throw Error('collision');
const seen=new Set();
for(let B=10*(1<<19);B<11*(1<<19);B+=2)seen.add(H(B));
if(seen.size!==262012)throw Error('width-23 support');
console.log('width 23 distinct:',seen.size);