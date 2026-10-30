const need=(x,m)=>{if(!x)throw Error(m);};
function binrow(n){let a=[1n];for(let k=1;k<=n;k++)a.push(a[k-1]*BigInt(n-k+1)/BigInt(k));return a;}
function mul(a,b){let c=Array(a.length+b.length-1).fill(0n);for(let i=0;i<a.length;i++)for(let j=0;j<b.length;j++)c[i+j]+=a[i]*b[j];return c;}
function inspect(p,a,r){
 const M=p**a,n=M*r,B=binrow(n),D=binrow(n/p),mm=BigInt(M);
 const numer=B.map((x,k)=>x+(k%p===0?BigInt(M-1)*D[k/p]:0n));
 need(numer.every(x=>x%mm===0n),"Burnside integrality");
 const P=numer.map(x=>x/mm);
 let Q=[P[n]];
 for(let k=n-1;k>=0;k--){let next=Array(Q.length+1).fill(0n);for(let j=0;j<Q.length;j++){next[j]-=Q[j];next[j+1]+=Q[j];}next[0]+=P[k];Q=next;}
 const bp=binrow(p),base=bp.map((x,j)=>j===0?0n:x*BigInt((p-j)%2?-1:1));
 let direct=[1n];for(let k=0;k<n/p;k++)direct=mul(direct,base);
 need(Q.every((x,k)=>mm*x===BigInt(M-1)*direct[k]+(k===n?1n:0n)),"independent expansion");
 const dep=q=>Q.findIndex(x=>q?x%BigInt(q)!==0n:x!==0n);
 need(dep(0)===n/p,"rational depth");need(dep(p)===(r%p?n-p+1:n),"modular depth");
 const other=[];for(let l=2;l<=M-1;l++){if((M-1)%l)continue;let prime=true;for(let d=2;d*d<=l;d++)if(l%d===0)prime=false;if(prime){need(dep(l)===n,"other characteristic");other.push(l);}}
 return {p,a,r,gaps:n,rationalDepth:dep(0),characteristicPDepth:dep(p),fullDepthCharacteristics:other};
}
const cases=[];for(const [p,a] of [[3,1],[3,2],[3,3],[5,1],[5,2],[7,1],[7,2]])for(let r=1;r<=p+1;r++)cases.push(inspect(p,a,r));
for(const r of [7,15,31,63])cases.push(inspect(5,2,r));
const doubling=Array.from({length:20},(_,i)=>{const j=i+1,r=2**j-1;return{j,vertices:25*2**j,depth:25*r-(j%4?4:0)};});
for(let j=0;j+5<doubling.length;j++)need(doubling[j+5].depth===2*doubling[j+4].depth+doubling[j+1].depth-2*doubling[j].depth,"order-five recurrence");
const S=new Set(['0,1','0,2','1,0','1,1','1,2','1,3','1,4','2,0','2,1','2,2','2,3','3,1']),n=25;
const A=Array.from({length:n},(_,a)=>Array.from({length:n},(_,b)=>S.has(`${(Math.floor(b/5)-Math.floor(a/5)+5)%5},${(b%5-a%5+5)%5}`)));
let homogeneousPairs=0;
for(let a=0;a<n;a++)for(let b=a+1;b<n;b++){need(A[a][b]!==A[b][a],"tournament");if(!Array.from({length:n},(_,x)=>x).some(x=>x!==a&&x!==b&&A[a][x]!==A[b][x]))homogeneousPairs++;}
need(homogeneousPairs===0,"seed homogeneous pairs");
let m=Array(n).fill(-1),used=Array(n).fill(false),nodes=0,count=0;m[0]=0;used[0]=true;
function dfs(depth){nodes++;if(depth===n){count++;return;}let best=-1,opts=null;for(let a=0;a<n;a++)if(m[a]<0){let c=[];for(let b=0;b<n;b++)if(!used[b]){let ok=true;for(let x=0;x<n;x++)if(m[x]>=0&&A[a][x]!==A[b][m[x]]){ok=false;break;}if(ok)c.push(b);}if(c.length===0)return;if(opts===null||c.length<opts.length){best=a;opts=c;}}for(const b of opts){m[best]=b;used[b]=true;dfs(depth+1);m[best]=-1;used[b]=false;}}
dfs(1);need(count===1,"seed full group");
console.log(JSON.stringify({status:"PASS",cases,doublingFirstEight:doubling.slice(0,8),seed:{homogeneousPairs,fixedZeroAutomorphisms:count,dfsNodes:nodes,fullAut:25*count},scope:"Exact coefficient expansions; first20 doubling recurrence terms from formula. Seed full group exhaustively verified; all-r tournament construction uses proof."},null,2));