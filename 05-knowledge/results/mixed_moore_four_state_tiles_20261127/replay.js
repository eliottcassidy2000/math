function verify(adj,t){
 const n=adj.length, out=Array.from({length:n},()=>0), seen=new Uint32Array(n);
 for(let u=0;u<n;u++){
   const tag=u+1;
   for(const v of adj[u]){if(v===u||seen[v]===tag)return false;seen[v]=tag;out[v]++}
   let reach=adj[u].length,t0=0;
   for(const v of adj[u])for(const w of adj[v]){
     if(w===u){t0++;continue}
     if(seen[w]===tag)return false;
     seen[w]=tag;reach++;
   }
   if(reach!==n-1||(t!==undefined&&t!==t0))return false;
 }
 return true;
}
function kautz(q){
 const V=[];for(let a=0;a<q;a++)for(let b=0;b<q;b++)if(a!==b)V.push([a,b]);
 const id=(a,b)=>a*(q-1)+b-(b>a?1:0);
 return V.map(([a,b])=>Array.from({length:q},(_,c)=>c).filter(c=>c!==b).map(c=>id(b,c)));
}
function induced(A,S){const id=new Map(S.map((x,i)=>[x,i]));return S.map(x=>A[x].filter(y=>id.has(y)).map(y=>id.get(y)))}
function assert(x){if(!x)throw Error("replay failed")}
const census={};
for(let n=1;n<=4;n++){
 const m=n*(n-1),counts={exact:0,relaxed:0,stars:0};
 for(let z=0;z<(1<<m);z++){
   const A=Array.from({length:n},()=>[]);let k=0;
   for(let a=0;a<n;a++)for(let b=0;b<n;b++)if(a!==b){if(z>>k&1)A[a].push(b);k++}
   if(verify(A,0))counts.exact++;
   if(verify(A)){counts.relaxed++;if(new Set(A.map(x=>x.length)).size>1)counts.stars++}
 }
 census[n]=counts;
}
assert(JSON.stringify(Object.values(census))===JSON.stringify([{exact:1,relaxed:1,stars:0},{exact:0,relaxed:1,stars:0},{exact:2,relaxed:5,stars:3},{exact:0,relaxed:4,stars:4}]));
console.log(census);
const K=kautz(5), counts={};
for(let mask=1;mask<1<<K.length;mask++){
 const S=[];for(let i=0;i<K.length;i++)if(mask>>i&1)S.push(i);
 const A=induced(K,S);
 if(verify(A))counts[S.length]=(counts[S.length]||0)+1;
}
assert(JSON.stringify(counts)===JSON.stringify({1:20,2:10,3:20,6:10,12:5,20:1}));
console.log("Induced perfect subsets of the five-letter Kautz graph:",counts);
const big=kautz(58);assert(big.length===3306&&verify(big,1));
console.log("Verified all ordered-pair entries for the 3306-vertex Kautz graph.");
