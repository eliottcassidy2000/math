function H(r){let n=r.length,N=1<<n,d=new Uint32Array(n*N);for(let i=0;i<n;i++)d[(1<<i)*n+i]=1;for(let s=1;s<N;s++)for(let i=0;i<n;i++){let h=d[s*n+i];if(h)for(let w=r[i]&~s;w;w&=w-1){let b=w&-w;d[(s|b)*n+31-Math.clz32(b)]+=h;}}let h=0;for(let i=0;i<n;i++)h+=d[(N-1)*n+i];return h;}
function flip(r,i,j){r[i]^=1<<j;r[j]^=1<<i;}

const p=[
 [10,[9,81,63,27,45,9],[486,484,488,499,487,448,384,256,0],[[2,8],[2,4],[7,8],[6,8],[2,4],[6,8],[2,8],[7,8]]],
 [12,[11,121,99,55,33,77,11],[1022,972,984,976,994,974,896,768,512,0],[[3,7],[3,9],[3,6],[1,3],[3,9],[1,3],[0,3],[3,7],[0,3],[3,6]]],
 [14,[13,169,143,91,13],[234,244,249,242,225,192,128,0],[[3,7],[4,6],[1,3],[1,4],[1,3],[4,6],[1,4],[3,7]]],
 [14,[39,117,39],[430,508,504,496,481,448,385,256,0],[[6,8],[6,8]]],
 [14,[65,65],[18,20,9,19,36,79,31],[]],
 [16,[15,225,195,135,15],[1022,1020,1016,1008,992,448,640,768,576,32],[[1,3],[0,4],[0,1],[1,3],[0,1],[0,4]]],
 [16,[45,165,75,105,45],[250,252,249,112,160,192,144,8],[[0,3],[1,7],[1,2],[0,3],[1,2],[1,7]]],
];
const constants=[[139,[10,108,73,112,39,69,17]],[787,[478,284,296,304,36,323,158,318,80]],
 [307,[222,124,40,240,36,65,148,54]],[1069,[478,316,408,48,352,325,142,314,72]]];
function assertion(x){if(!x)throw Error('audit assertion failed');}
function gcd(a,b){return b===0n?a:gcd(b,a%b);}
function run(log){
 for(const [height,r] of constants){
  assertion(r.every((v,i)=>!(v>>i&1)&&v>>r.length===0));
  for(let i=0;i<r.length;i++)for(let j=0;j<i;j++)assertion((r[i]>>j&1)+(r[j]>>i&1)===1);
  assertion(H(r)===height);log('constant',height,r.length);
 }
 const S=[[6,12,8,1],[10,20,25,18,1],[40,53,57,50,1,16]], sh=()=>S.reduce((a,r)=>a*H(r),1),program=[sh()];
 assertion(String(S.map(H))==='5,13,17');
 for(const [i,u,v] of [[1,0,1],[0,0,1],[1,2,4],[2,4,5],[2,0,1]]){flip(S[i],u,v);program.push(sh());}
 assertion(String(program)==='1105,765,459,561,297,429');log('selector',program);
 let L=1n;
 for(const b of [10n,12n,14n,16n]){const x=b*b-1n;L=L*x/gcd(L,x);}
 assertion(L===109395n);
 for(const y of [1n,3n]){
  const M=L*y,B=M-1n,h=M-2n,product=y===1n?139n*787n:307n*1069n;
  assertion(product===h);
  const occupied=new Set();
  for(const [base,phases,A,edits] of p){
   const r=A.slice(),trace=[H(r)];for(const e of edits){flip(r,...e);trace.push(H(r));}assertion(String(A)===String(r));
   let i=0;for(const x of trace)if(x===phases[i])i++;assertion(i===phases.length||edits.length===0&&i===1);
   const b=BigInt(base),d=M/(b+1n),den=b-1n,K=h*d/den;assertion(K*den===h*d&&K===(L/((b+1n)*den))*y*product);
   const target=phases.map(x=>BigInt(x)*K),go=g=>2n*g>M?2n*g-M:M-2n*g;
   for(let j=0;j+1<target.length;j++){
    assertion(target[j]%h===0n&&go(target[j]/h)*h===target[j+1]);
    const key=String(target[j]);assertion(!occupied.has(key));occupied.add(key);
   }
   const rden=M/gcd(M,target[0]/h);
   log('base',String(B),'source',base,'K',String(K),'reduced denominator',String(rden),'orbit',target.map(String),'scaled literal trace',trace.map(x=>String(BigInt(x)*K)));
  }
 }
}

run(console.log);

const compact=[444,401,138,322,76,414,423,24,148];
const program=[[1,2],[1,4],[3,6],[0,5],[7,8]];
const counts=[H(compact)];
for(const [i,j] of program){flip(compact,i,j);counts.push(H(compact));}
assertion(String(counts)==="1105,765,955,561,405,429");
for(const [i,j] of program.slice().reverse())flip(compact,i,j);
assertion(String(compact)==="444,401,138,322,76,414,423,24,148");
console.log("compact selector",counts);
// Count every positive gap cycle, not just the seven transported ones.
function gapCycles(M) {
 const seen=new Uint8Array(M),by={};
 for(let x=1;x<M;x+=2)if(!seen[x]){
  let y=x,k=0;do{seen[y]=1;k++;y=Math.abs(2*y-M);}while(y!==x);
  by[k]=(by[k]||0)+1;
 }
 return by;
}
for(const [M,n] of [[109395,701],[328185,1142]]){
 const by=gapCycles(M);
 assertion(Object.values(by).reduce((a,b)=>a+b,0)===n);
 console.log("positive gap cycle totals",M-1,n,by);
}
