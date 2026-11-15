
function* words(k,a,b,p=[]){if(!a&&!b){yield p;return;}if(a)yield* words(k+1,a-1,b,[...p,k]);if(b)yield* words(k+1,a,b-1,p);}
function floor(n){return n%3===0?3*5**((n-3)/3):n%3===1?5**((n-1)/3):9*5**((n-5)/3);}
function* shapes(n, pref=[], c=1){if(!n){yield [pref,c];return;}for(const m of [3,4,5])if(n>=m-1)yield* shapes(n-m+1,[...pref,m],c*(2**(m-2)+1));}
function chain(s){let n=1+s.reduce((x,y)=>x+y-1,0),bits=Array.from({length:n},(_,i)=>((1<<n)-1)^((1<<(i+1))-1)),end=[0],inner=[];let i=0;for(const m of s){const j=i+m-1;bits[i]^=1<<j;bits[j]^=1<<i;end.push(j);for(let k=i+1;k<j;k++)inner.push(k);i=j;}return {n,bits,cyc:[...end.reverse(),...inner]};}
function order(c,v){const i=c.indexOf(v);return [...c.slice(i+1),...c.slice(0,i)];}
function* perm(xs){if(!xs.length){yield [];return;}for(let i=0;i<xs.length;i++)for(const y of perm(xs.filter((_,j)=>i!==j)))yield [xs[i],...y];}
function h(e){const n=e.length,limit=1<<n,d=new Float64Array(n*limit);for(let i=0;i<n;i++)d[(1<<i)*n+i]=1;for(let s=1;s<limit;s++){for(let z=s;z;z&=z-1){let k=31-Math.clz32(z&-z),x=d[s*n+k];if(!x)continue;for(let b=e[k]&~s;b;b&=b-1){let j=31-Math.clz32(b&-b);d[(s|(1<<j))*n+j]+=x;}}}return d.subarray((limit-1)*n).reduce((a,b)=>a+b,0);}
function between(n,m){let target=floor(n)+4,L=n-m+1,checked=0,badcases=0,win=[];
for(let [sa,ha] of shapes(m-1))for(let [sb,hb] of shapes(L-1)){if(ha*hb!==floor(n)||ha!==floor(m)||hb!==floor(L))continue;let A=chain(sa),B=chain(sb);
for(let v=0;v<m;v++)for(let w=0;w<L;w++){let ao=order(A.cyc,v),bo=order(B.cyc,w),p=ao.length,q=bo.length;let bmap=Array(L).fill(v);for(let j=0;j<q;j++)bmap[bo[j]]=m+j;
let base=[...A.bits,...Array(q).fill(0)];for(let i=0;i<L;i++)for(let j=0;j<L;j++)if(B.bits[i]>>j&1)base[bmap[i]]|=1<<bmap[j];
const opts=m===4?[0,1,2]:[0];
const patterns=[];for(let arr of perm(ao)){let cyc=[v,...arr],bad=[];for(let i=0;i<m;i++)if(!(A.bits[cyc[i]]>>cyc[(i+1)%m]&1))bad.push(i);if(bad.length!==1)continue;const g=bad[0];for(const opt of opts){let gs=opt===0?[-1]:Array.from({length:m},(_,i)=>i).filter(j=>opt===1?j<g:j>g);for(const j of gs){const P=[];for(let k=0;k<m;k++){P.push(cyc[k]);if(k===g)for(let i=(opt===1?1:0);i<q-(opt===2?1:0);i++)P.push(m+i);if(k===j)P.push(m+(opt===1?0:q-1));}patterns.push(P);}}}
if(!patterns.length)continue;
for(const apos of words(0,p,q)){checked++;const e=[...base],chosen=new Set(apos),bpos=[];for(let j=0;j<p+q;j++)if(!chosen.has(j))bpos.push(j);for(let i=0;i<p;i++)for(let j=0;j<q;j++){if(apos[i]<bpos[j])e[ao[i]]|=1<<(m+j);else e[m+j]|=1<<ao[i];}
if(!patterns.some(P=>P.every((x,j)=>e[x]>>P[(j+1)%n]&1)))continue;badcases++;
if(h(e)===target)win.push({sa,sb,v,w,apos,e});}
}}
return {n,m,checked,badcases,solutions:win.length,ex:win[0]};}

let rooted=0,candidateMatrices=0;
for(let n=9;n<=14;n++)for(let m=4;m<=5;m++)if(n-m+1>=6){
  const result=between(n,m);
  console.log(JSON.stringify(result));
  if(result.solutions)throw Error("Unexpected second-gap witness");
  rooted+=result.checked;candidateMatrices+=result.badcases;
}
console.log({rooted,candidateMatrices});
if(rooted!==250316 || candidateMatrices!==54476)throw Error('Unreconciled enumeration totals');
if(h(chain([4,2,4]).bits)!==49)throw Error('Order-eight boundary failed');
