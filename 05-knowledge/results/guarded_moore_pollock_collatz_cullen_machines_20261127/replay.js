const A=(x)=>{if(!x)throw Error("replay failed")};
function F(x){return x*(x*x-1)/6}
function pollock(N){
 const M=(1n<<BigInt(N+1))-1n, b=[1n,0n,0n,0n,0n,0n];
 for(let f=1,g=3,h=3; f<=N;[f,g,h]=[f+g,g+h,h+1]){
   for(let j=1;j<=5;j++)b[j]|=(b[j-1]<<BigInt(f))&M;
 }
 return b;
}
function bit(b,n){return Boolean(b>>BigInt(n)&1n)}
const X=1000000,P=pollock(X);
A(P.reduce((a,b)=>a|b,0n)===(1n<<BigInt(X+1))-1n);
const fewer=P.slice(0,5).reduce((a,b)=>a|b,0n);
A(!bit(fewer,343867)&&bit(P[5],343867));
let defect=((1n<<BigInt(X+1))-1n)^fewer,cnt=0;
for(const nibble of defect.toString(16))cnt+=parseInt(nibble,16).toString(2).replace(/0/g,"").length;
A(cnt===241);
console.log("Pollock: all 0..1000000; four-term defects",cnt);
function* permutations(n){
 const V=Array.from({length:n},(_,i)=>i);
 function* rec(i){if(i===n){yield V.slice();return}for(let j=i;j<n;j++){[V[i],V[j]]=[V[j],V[i]];yield* rec(i+1);[V[i],V[j]]=[V[j],V[i]];}}
 yield* rec(0);
}
function moore(d,budget=Infinity){
 const b=d-1, I=Array.from({length:b},(_,i)=>i), P=I.map(()=>Array(d).fill(null));
 P.push(Array(d).fill(null));
 for(let j=1;j<d;j++)P[0][j]=P[j][0]=I;
 function check(){
   for(let i=0;i<d;i++)for(let j=i+1;j<d;j++)for(let x=0;x<b;x++){
     const used=Array(b).fill(false);
     if(P[i][j])used[P[i][j][x]]=true;
     for(let h=0;h<d;h++)if(h!==i&&h!==j&&P[i][h]&&P[h][j]){
       const y=P[h][j][P[i][h][x]];if(used[y])return false;used[y]=true;
     }
   }
   return true;
 }
 const jobs=[];for(let i=1;i<d;i++)for(let j=i+1;j<d;j++)jobs.push([i,j]);
 let tested=0;
 function visit(pos){
   if(pos===jobs.length)return P.map(L=>L.slice());
   const [i,j]=jobs[pos];
   for(const p of permutations(b)){
     if(tested++>=budget)return "budget";
     const q=[];p.forEach((y,x)=>q[y]=x);
     P[i][j]=p;P[j][i]=q;
     if(check()){const ans=visit(pos+1);if(ans)return ans;}
   }
   P[i][j]=P[j][i]=null;
   return null;
 }
 return {bank:visit(0),tested};
}
function expand(P){
 const d=P.length,b=d-1, V=1+d+d*b, adj=Array.from({length:V},()=>[]);
 function add(i,j){adj[i].push(j);adj[j].push(i)}
 const leaf=(i,x)=>1+d+i*b+x;
 for(let i=0;i<d;i++){add(0,i+1);for(let x=0;x<b;x++)add(i+1,leaf(i,x));}
 for(let i=0;i<d;i++)for(let j=i+1;j<d;j++)for(let x=0;x<b;x++)add(leaf(i,x),leaf(j,P[i][j][x]));
 return adj;
}
for(const d of [2,3,7]){
 const {bank,tested}=moore(d,d===7?2000:Infinity);
 if(d===7){A(bank==="budget");console.log("degree 7 bounded smoke test, no claim",tested-1);continue}
 const adj=expand(bank);A(adj.length===d*d+1);
 for(let i=0;i<adj.length;i++){let all=[...adj[i],...adj[i].flatMap(j=>adj[j].filter(k=>k!==i))];A(all.length===adj.length-1&&new Set(all).size===all.length)}
 console.log("Moore verified",d,adj.length,tested);
}
function pow(a,m,n){let r=1n;for(;m;m>>=1n,a=a*a%n)if(m&1n)r=r*a%n;return r}
function C(m,a,M){return (m*pow(a,m,M)+1n)%M}
function ord(a,p){let k=1n,z=a%p;while(z!==1n){z=z*a%p;k++}return k;}
function lift(p,a,j,R,target=0n){
 const d=ord(a,p), slope=d*pow(a,j,p)%p, inv=pow(slope,p-2n,p);
 let m=j,P=1n,out=[];
 for(let r=0;r<R;r++){
   const residue=C(m,a,P*p), rep=target%(P*p), E=(rep-residue)/P;
   A((rep-residue)%P===0n);
   m+=d*P*((E*inv%p+p)%p);P*=p;A(C(m,a,P)===target%P);out.push(m.toString());
 }
 return out;
}
for(const v of [3n,5n,7n,13n,29n,563n]){
 const a=2n,d=ord(a,v), branch=(v-1n)%d;
 console.log("Cullen base 2, phase zero, p",String(v),lift(v,a,branch,5));
 let fact=1n,p2=v*v,u=0n;
 for(let a=1n;a<v;a++){fact=fact*a%p2;u+=C(v-1n,a,p2)/v;}
 const W=(fact+1n)/v%v;
 A((u+1n+W)%v===0n);
}
for(let n=1n;n<1000n;n++){
 let x=n;
 for(let r=1n;r<=4n;r++){
   if(x%4n!==3n){A((n+1n)%4n**r!==0n);break}
   x=(9n*x+5n)/4n;A((n+1n)%4n**r===0n&&4n**r*(x+1n)===9n**r*(n+1n));
 }
}
console.log("Two-odd Collatz guard replay: starts 1..999");
function hoffman(){
 const adj=Array.from({length:50},()=>[]);
 const L=(a,b)=>5*a+b, R=(a,b)=>25+5*a+b;
 const add=(u,v)=>{adj[u].push(v);adj[v].push(u)};
 for(let i=0;i<5;i++)for(let x=0;x<5;x++){
   add(L(i,x),L(i,(x+1)%5));add(R(i,x),R(i,(x+2)%5));
   for(let j=0;j<5;j++)add(L(i,x),R(j,(x+i*j)%5));
 }
 return adj;
}
function extract(adj,r){
 const parents=adj[r], d=parents.length, b=d-1;
 const leaves=parents.map(p=>adj[p].filter(v=>v!==r));
 for(let j=1;j<d;j++)leaves[j]=leaves[0].map(v=>leaves[j].find(y=>adj[v].includes(y)));
 const P=leaves.map((L,i)=>leaves.map((M,j)=>i===j?null:L.map(x=>M.findIndex(y=>adj[x].includes(y)))));
 for(let i=0;i<d;i++)for(let j=i+1;j<d;j++)for(let x=0;x<b;x++){
   const all=[P[i][j][x]];
   for(let h=0;h<d;h++)if(h!==i&&h!==j)all.push(P[h][j][P[i][h][x]]);
   A(new Set(all).size===b&&all.every(y=>y>=0&&y<b));
 }
 const I=Array.from({length:b},(_,i)=>i), G=new Set([I.join()]), Q=[I], gens=P.flat().filter(Boolean);
 for(let i=0;i<Q.length;i++)for(const p of gens){
   const q=Q[i].map(x=>p[x]), key=q.join();if(!G.has(key)){G.add(key);Q.push(q)}
 }
 return G.size;
}
const hg=extract(hoffman(),0);A(hg===720);console.log("Hoffman-Singleton: valid bank, transport group size",hg);
function fermat(a,p){return (pow(a,p-1n,p*p)-1n)/p;}
const plateau=[];
for(let p=5n;p<10000n;p+=2n){
 let is=true;for(let d=3n;d*d<=p;d+=2n)if(p%d===0n){is=false;break}
 if(!is)continue;
 const g=9n*pow(4n,p*(p-1n)-1n,p*p)%(p*p);
 const D2=C(p-1n,2n,p*p)/p, D3=C(p-1n,3n,p*p)/p;
 A((fermat(g,p)-2n*(D2-D3))%p===0n);
 if(D2===D3)plateau.push(String(p));
}
console.log("First-order Collatz clock plateaus for 5<=p<10000:",plateau);
const orders=[];
for(const M of [23n,23n**2n,23n**3n]){
 const g=9n*pow(4n,(M/23n)*22n-1n,M)%M;
 orders.push(String(ord(g,M)));
}
A(orders.join(",")==="11,11,253"); console.log("Orders at 23,23^2,23^3:",orders);

