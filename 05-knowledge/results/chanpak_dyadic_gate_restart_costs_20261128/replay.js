
function masks(n,m){let A=new Uint32Array(n);for(let i=0,k=0;i<n;i++)for(let j=i+1;j<n;j++,k++)if((m>>>k)&1)A[i]|=1<<j;else A[j]|=1<<i;return A}
function count(A){let n=A.length,N=1<<n,I=A.map((_,i)=>A.reduce((x,r,j)=>x|((r>>i&1)<<j),0));let H=new Float64Array(N),E=new Float64Array(n*N),S=new Float64Array(n*N);H[0]=1;for(let T=1;T<N;T++)for(let z=T;z;z&=z-1){let j=31-Math.clz32(z&-z),R=T^(1<<j),x=+!R,y=x;for(let q=R&I[j];q;q&=q-1)x+=E[R*n+31-Math.clz32(q&-q)];for(let q=R&A[j];q;q&=q-1)y+=S[R*n+31-Math.clz32(q&-q)];E[T*n+j]=x;S[T*n+j]=y;H[T]+=x;}return {H,E,S,b:H[N-1]}}
function data(A){let n=A.length,F=(1<<n)-1,{H,E,S,b}=count(A);let w=Array.from({length:n},(_,i)=>E[F*n+i]);for(let t=1;t<=F;t++)for(let i=0;i<n;i++)w[i]+=S[t*n+i]*H[F^t];return {w,b,Q:w.reduce((a,b)=>a+b,0)}}
function pointed(A,c){let n=A.length;let B=Array.from(A,(x,i)=>x+((c>>i&1)?0:1<<n));B.push(c);let rev=B.map((_,i)=>B.reduce((v,d,j)=>v|((d>>i&1)<<j),0));let {H,E,S}=count(B),D=count(rev),F=(1<<(n+1))-1;let u=0,v=0;for(let t=1;t<=F;t++){u+=H[F^t]*S[t*(n+1)+n];v+=H[F^t]*D.S[t*(n+1)+n];}return [H[F],E[F*(n+1)+n],D.E[F*(n+1)+n],u,v]}
function find(n,A,U){let d=data(A);for(let p=1;p<=U*d.Q/d.b;p*=2){let t=p*d.b-(U-1)*d.Q;if(t<0||t>d.Q)continue;
for(let c=0;c<1<<n;c++){let z=0;for(let i=0;i<n;i++)if(c>>i&1)z+=d.w[i];if(z===t)return [p,c,d,pointed(A,c)];}}return null}

function canon(A){const n=A.length,L=[];for(let k=0;k<n;k++)L.push([]);for(let i=0;i<n;i++)L[A[i].toString(2).split("1").length-1].push(i);let cells=L.filter(x=>x.length),q=[],min=Infinity;
 function cell(t){if(t==cells.length){let z=0;for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)z=2*z+(A[q[i]]>>q[j]&1);min=Math.min(z,min);return;}
 const C=cells[t];function fill(p){if(p===C.length){cell(t+1);return;}for(let j=p;j<C.length;j++){[C[j],C[p]]=[C[p],C[j]];q.push(C[p]);fill(p+1);q.pop();[C[j],C[p]]=[C[p],C[j]];}}fill(0);}
 cell(0);return min;}
function strong(A){for(const D of [A,A.map((_,i)=>A.reduce((x,r,j)=>x|((r>>i&1)<<j),0))]){let s=1,v=0;while(s!==v){v=s;for(let i=0;i<D.length;i++)if(v>>i&1)s|=D[i];}if(s!==(1<<D.length)-1)return false;}return true;}
function catalog(maxn,cap){let U=[[0]], sizes=[1], rows=new Map([[1,new Set([1])]]);
 for(let n=2;n<=maxn;n++){let X=new Map();for(const A of U)for(let c=0;c<1<<(n-1);c++){let D=A.map((r,i)=>r+((c>>i&1)?1<<(n-1):0));D.push(((1<<(n-1))-1)^c);let z=canon(D);if(!X.has(z))X.set(z,D);}
 U=[...X.values()];sizes.push(U.length);
 if(n>=3)for(const A of U)if(strong(A)){let d=data(A);if(d.b>cap)continue;if(!rows.has(d.b))rows.set(d.b,new Set());rows.get(d.b).add(d.Q);}
 }
 let all=new Map();for(let b=1;b<=cap;b+=2){let S=new Set(rows.get(b)??[]);for(let a=3;a<b;a+=2)if(b%a===0)for(const x of (all.get(a)??[]))for(const y of (all.get(b/a)??[]))S.add(x*y);all.set(b,S);}
 let candidates=[];for(const [b,S] of all)for(const Q of S)for(let p=1;p*b<=63*Q;p*=2)if(62*Q<=p*b&&p*b<=63*Q)candidates.push([b,Q,p]);
 return {sizes,strong:[...rows].sort((a,b)=>a[0]-b[0]).map(([b,s])=>[b,[...s].sort((x,y)=>x-y)]),candidates};
}

let R=catalog(7,43);
function check(v){if(!v)throw Error("verification failed");}
check(R.sizes.join()=="1,1,2,4,12,56,456"&&R.candidates.length===0);
console.log("isomorphism class catalogue",JSON.stringify(R.sizes));
console.log("strong (b,Q) profiles through 43",JSON.stringify(R.strong));
console.log("upper branch: possible powers in all ordinal joins through b=43",R.candidates.length);
for(const [n,m,c,U,target,b0] of [[5,10,29,55,512,9],[8,238981866,127,63,2048,45]]){
 const A=masks(n,m),{w,b,Q}=data(A),[j,e,s,u,v]=pointed(A,c);
 const sum=w.reduce((y,w,i)=>y+((c>>i&1)?w:0),0);
 const h=(303-9)*b+9*s+(9*U-9)*(v-b)+9*(j-s),z=(9*U-9)*Q+9*u;
 check(b===b0&&sum===u-e&&z-9*e===9*b*target);
 console.log("suffix",JSON.stringify({n,m,c,U,b,Q,w,sum,C:[j,e,s,u,v],state:[h,9*e,z,9*b],guard:target}));
}
