function verifyRegularAugmentation() {
 const check=b=>{if(!b)throw Error("Verification failed");};
 const rot=x=>((x<<1)&31)|(x>>4),reps=[];
 for(let x=0;x<32;x++){
  let y=x,min=x,stabilizer=0;
  for(let j=0;j<5;j++){stabilizer+=+(y===x);min=Math.min(min,y);y=rot(y);}
  if(x===min)reps.push({mask:x,e:+(stabilizer===5)});
 }
 check(JSON.stringify(reps.map(x=>x.mask))===JSON.stringify([0,1,3,5,7,11,15,31]));
 const seen=new Set(),coefficients=Array(7).fill(0);
 for(let n=0;n<8**5;n++){
  let q=n;const a=Array.from({length:5},()=>{const x=q%8;q=Math.floor(q/8);return x;});
  const codes=a.map((_,j)=>a.slice(j).concat(a.slice(0,j)).reduce((s,x)=>8*s+x,0));
  const canonical=Math.min(...codes);
  if(seen.has(canonical))continue;seen.add(canonical);
  const e=a.reduce((s,x)=>s+reps[x].e,0)+ +(a.every(x=>x===a[0]));
  coefficients[e]++;
 }
 check(JSON.stringify(coefficients)===JSON.stringify([1554,2598,1728,576,96,6,2]));
 check(seen.size===6560);
 const permutations=[];
 function permute(a,remaining){
  if(!remaining.length){permutations.push(a);return;}
  remaining.forEach((x,i)=>permute(a.concat(x),remaining.filter((_,j)=>j!==i)));
 }
 permute([],[0,1,2,3,4,5]);
 const classes=new Map();
 const encode=(A,p)=>{let s="";for(let i=0;i<6;i++)for(let j=i+1;j<6;j++)s+=A[p[i]][p[j]];return s;};
 for(let mask=0;mask<32;mask++){
  const A=Array.from({length:6},()=>Array(6).fill(0));
  for(let i=0;i<5;i++)for(let j=0;j<5;j++)if(i!==j)A[i][j]=+((j-i+5)%5<=2);
  for(let i=0;i<5;i++){A[5][i]=(mask>>i)&1;A[i][5]=1-A[5][i];}
  const base=encode(A,[0,1,2,3,4,5]);let best=null,aut=0;
  for(const p of permutations){const s=encode(A,p);aut+=+(s===base);if(best===null||s<best)best=s;}
  if(classes.has(best))check(classes.get(best)===aut);classes.set(best,aut);
 }
 const orders=[...classes.values()].sort((a,b)=>a-b);
 check(JSON.stringify(orders)===JSON.stringify([1,1,1,1,1,5,5]));
 return {status:"PASS",baseSubsetRepresentatives:reps,subsetStabilizerPolynomial:coefficients,
  marked26:seen.size,unmarked26:seen.size-1,rigid26:coefficients[0],
  unmarked6:classes.size,automorphismOrders6:orders,
  scope:"Extensions of the specified regular parent, not all26-vertex tournaments"};
}

console.log(JSON.stringify(verifyRegularAugmentation(), null, 2));
