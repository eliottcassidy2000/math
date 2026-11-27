const A=x=>{if(!x)throw Error("replay failed")};
function next(x,c,b){return [2*x+b,c+(1-b)*Math.ceil(x/2)];}
function suffix(x,c,w){for(const b of w)[x,c]=next(x,c,Number(b));return c;}
for(let m=2;m<=128;m++){
 const n=2*m*m, seen=new Uint8Array(n),Q=[1];
 seen[1]=1;const keys=new Set();
 for(let i=0;i<Q.length;i++){
   const x=Q[i]%(2*m),c=Math.floor(Q[i]/(2*m));
   const sig=["0","10","110"].map(w=>suffix(x,c,w)%m).join();
   A(!keys.has(sig));keys.add(sig);
   for(const b of [0,1]){
     const [y,d]=next(x,c,b),key=y%(2*m)+2*m*(d%m);
     if(!seen[key]){Q.push(key);seen[key]=1;}
   }
 }
 A(Q.length===n);
}
console.log("Reachability and 3-probe distinguishability verified moduli 2..128.");
function pct(b){
 const C=3n**BigInt(1-b);
 return [[11n*5n**4n*C,7n*4n],[7n,11n],[13n*25n*C,14n],[13n,7n],
         [17n*5n**BigInt(b),13n],[38n,85n],[17n,19n],[23n,17n]];
}
function fappend(x,c,b){
 const P=pct(b);let z=7n*2n**BigInt(x)*3n**BigInt(c),i=0;
 while(true){const f=P.find(([a,d])=>z*a%d===0n);if(!f)break;z=z*f[0]/f[1];A(++i<1000);}
 const [y,h]=next(x,c,b);A(z===23n*2n**BigInt(y)*3n**BigInt(h));return [y,h];
}
for(let x=0;x<=128;x++)for(let c=0;c<5;c++)for(let b=0;b<2;b++)fappend(x,c,b);
console.log("FRACTRAN microprogram: 1290 port inputs.");
function dp(word){
 const n=word.length+1, out=Array.from({length:n},()=>[]);
 for(let a=0;a<n-1;a++){for(let b=a+1;b<n-1;b++)out[a].push(b);if(word[a]==="1")out[n-1].push(a);else out[a].push(n-1);}
 const D=Array.from({length:1<<n},()=>Array(n).fill(0));for(let i=0;i<n;i++)D[1<<i][i]=1;
 for(let S=1;S<1<<n;S++)for(let j=0;j<n;j++)if(D[S][j])for(let k of out[j])if(!(S>>k&1))D[S|1<<k][k]+=D[S][j];
 return D.at(-1).reduce((x,y)=>x+y,0);
}
const rows=["1000","1110","100010","111010"].map(w=>[w,dp(w)]);A(rows.map(x=>x[1]).join()==="9,9,27,39");
console.log("Literal tournament collision and split:",rows);
