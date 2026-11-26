function pairs(n){let r=[];for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)r.push([i,j]);return r}
function cyc(n,m){let a=Array(n).fill(0),k=0;for(const [i,j] of pairs(n)){let b=i===0&&j===n-1;if(!b&&j!==i+1)b=!!(m>>k++&1);if(b)a[j]|=1<<i;else a[i]|=1<<j}return a}
function paths(a){let n=a.length,u=(1<<n)-1,h=new Int32Array(u+1),d=new Int32Array((u+1)*n);h[0]=1;for(let s=1;s<=u;s++)for(let b=s;b;b&=b-1){let v=31-Math.clz32(b&-b),t=s^(b&-b),x=t?0:1;for(let c=t;c;c&=c-1){let j=31-Math.clz32(c&-c);if(a[j]>>v&1)x+=d[t*n+j]}h[s]+=d[s*n+v]=x}return h}
function split(a,v){let n=a.length,b=[...a,0];for(let i=0;i<n;i++)if(i===v)b[v]|=1<<n;else if(a[v]>>i&1)b[n]|=1<<i;else b[i]|=1<<n;return b}
function fill(A,B){const R=[10,4,9,18,39,15],M=[[0],A,[0],[0],[0],B],I=M.flatMap((m,i)=>m.map((_,j)=>[i,j]));return I.map(([i,j])=>I.reduce((h,[k,l],v)=>h|(((i===k?M[i][j]>>l:R[i]>>k)&1)<<v),0))}
function chain(n){return Array.from({length:n},(_,i)=>(1<<n)-(1<<(i+1)))}
function DQ(a){const h=paths(a),u=h.length-1;return [h[u],h.reduce((s,x,i)=>s+x*h[u^i],0)]}
const bank={};
for(const [n,m] of [[5,0],[5,1],[6,0],[6,11],[5,4],[6,1],[6,4],[6,9]]){
 let a=cyc(n,m);if(n===5&&m===4)a=[62,...a.map(t=>t<<1)];bank[n+","+m]=a;
 const [D,Q]=DQ(a);let u=(1<<a.length)-1;
 const cl=a.map((_,v)=>paths(split(a,v)).at(-1));
 for(let w=1;w<=4;w++)if(paths(fill(a,chain(w))).at(-1)!==(4*Q+D)*2**w+3*Q-D)throw Error("fill mismatch");
 console.log(JSON.stringify({n,m,a,D,Q,c:4*Q+D,d:3*Q-D,del:a.map((_,v)=>paths(a)[u^(1<<v)]),clone:cl}));}

function root(n){let x=1n<<BigInt(Math.ceil(n.toString(2).length/2));for(;;){let y=(x+n/x)/2n;if(y>=x)return x;x=y;}}
const hits=[];
for(let u=1;u<=256;u++)for(let v=1;v<=256;v++){
 const a=1n<<BigInt(u),b=1n<<BigInt(v),z=4n*a*b+3n*a+b+1n,r=root(z);
 if(r*r===z&&u!==v)hits.push([u,v,String(r)]);
}
if(JSON.stringify(hits)!=='[[7,3,"67"]]')throw Error(JSON.stringify(hits));
console.log("off-diagonal through 256: "+JSON.stringify(hits));
function oddminus(x){let y=3*x-1;while(y%2===0)y/=2;return y;}
for(const s of [7,9,15,17,21,23]){
 let p=[s];while(!p.slice(0,-1).includes(p.at(-1)))p.push(oddminus(p.at(-1)));console.log('minus: '+p.join(' '));
}

for(const [u,v] of [[0,0],[0,1],[1,0],[0,2],[2,0],[0,3],[3,0],[0,4],[4,0]]){
 const a=chain(u),b=chain(v);console.log('empty '+JSON.stringify({u,v,actual:paths(fill(a,b)).at(-1),formal:4*2**(u+v)+3*2**u+2**v-1}));
}

const R=[10,4,9,18,39,15];
function skeleton(i,j){const L=[1,i,1,1,1,j];let answer=0;
 function f(last,d){if(!d){answer++;return;}for(let k=0;k<6;k++)if(L[k]&&(last===-1||(R[last]>>k&1))){L[k]--;f(k,d-1);L[k]++;}}
 f(-1,4+i+j);return answer;
}
const st=[0,1,2].map(i=>[0,1,2].map(j=>skeleton(i,j)));
if(JSON.stringify(st)!=='[[5,9,4],[13,23,9],[6,11,4]]')throw Error(JSON.stringify(st));
console.log("skeleton table "+JSON.stringify(st));
for(const a of [[],[0],chain(2),bank["5,1"],bank["6,11"]])for(const b of [[],[0],chain(2),bank["5,0"]]){
 function weights(a){if(!a.length)return [1,0,0];const [d,q]=DQ(a);return [0,d,q-2*d];}
 const A=weights(a),B=weights(b);const h=A.reduce((v,x,i)=>v+B.reduce((w,y,j)=>w+x*y*st[i][j],0),0);
 if(h!==paths(fill(a,b)).at(-1))throw Error("occupancy mismatch");
}
console.log("20 occupancy controls passed");

for(let s=1;s<63;s++)if(s&(s-1)){
 if(R.every((r,i)=>(s>>i&1)||(r&s)===0||(r&s)===s))throw Error("proper quotient module");
}
for(const [n,m] of [[5,0],[5,1],[6,0],[6,11],[5,4],[6,1]]){
 const a=bank[n+","+m],b=bank[n+","+(n===5?(m^1):11)];
 if(n===5&&m===0&&a.reduce((x,r,i)=>x+(r^b[i]).toString(2).split('1').length-1,0)!==2)throw Error('not one flip');
}
console.log("seed has no nontrivial proper module; height-nine pair differs by one arc");
