function verifyDepth5Hostile(){
 const req=(x,m)=>{if(!x)throw Error(m)},U=125,V=25,mod=(x,n)=>(x%n+n)%n,idx=(u,v)=>mod(u,U)*V+mod(v,V),xy=x=>[Math.floor(x/V),x%V],g=x=>{let[u,v]=xy(x);return idx(66*u+70*v,14*u+16*v)},N=U*V;
 const seen=new Set(),orbits=[],lookup=new Map();
 for(let x=1;x<N;x++){if(seen.has(x))continue;let o=[],y=x;do{o.push(y);seen.add(y);y=g(y)}while(y!==x);o.forEach(z=>lookup.set(z,orbits.length));orbits.push(o)}
 const used=new Set(),vars=new Map(),pairs=[];
 for(let i=0;i<orbits.length;i++){if(used.has(i))continue;let[u,v]=xy(orbits[i][0]),j=lookup.get(idx(-u,-v)),bit=pairs.length;req(i!==j,"no selfnegative");used.add(i);used.add(j);pairs.push([i,j]);orbits[i].forEach(x=>vars.set(x,[bit,0]));orbits[j].forEach(x=>vars.set(x,[bit,1]))}
 req(pairs.length===74,"74 bits");
 req((26**5)%125===1 && 26%125!==1 && (26*26)%125===51,"extra scalar order and nontrivial norm");
 req(lookup.get(idx(26,0))!==lookup.get(idx(1,0)),"extra scalar outside H");
 const action=pairs.map(([i])=>{let[u,v]=xy(orbits[i][0]);return vars.get(idx(26*u,26*v))});
 const eqs=[null];
 for(let r=1;r<=4;r++){let a=5**Math.ceil(r/2),b=5**Math.floor(r/2),gens=[[a,0],[0,b]].filter(([u,v])=>idx(u,v)!==0),eq=new Map();
 for(let x=1;x<N;x++){let[u,v]=xy(x);if(u%a===0&&v%b===0)continue;let[i,e]=vars.get(x);for(let[c,d]of gens){let[j,f]=vars.get(idx(u+c,v+d)),row=[i,j,e^f];eq.set(row.join(","),row)}}eqs.push([...eq.values()])}
 function solve(cuts){
 let parent=Array.from({length:74},(_,i)=>i),parity=Array(74).fill(0),components=74;
 function find(x){if(parent[x]!==x){let p=parent[x];parent[x]=find(p);parity[x]^=parity[p]}return parent[x]}
 function join(i,j,e){let a=find(i),b=find(j),p=parity[i]^parity[j]^e;if(a===b)return p===0;parent[a]=b;parity[a]=p;components--;return true}
 for(let i=0;i<74;i++)if(!join(i,action[i][0],action[i][1]))return {count:0n};
 for(let r=1;r<=4;r++)if(cuts>>(r-1)&1)for(let[i,j,e]of eqs[r])if(!join(i,j,e))return {count:0n};
 let mask=0n;for(let i=0;i<74;i++){find(i);if(parity[i])mask|=1n<<BigInt(i)}
 return {count:1n<<BigInt(components),mask,components};
 }
 let primeFixed=0n,fixedCounts=[];
 for(let c=0;c<16;c++){let f=solve(c).count,par=0;for(let z=c;z;z>>=1)par^=z&1;primeFixed+=(par?-1n:1n)*f;fixedCounts.push(f.toString())}
 req(primeFixed===17045913600n,"extra scalar fixes prime states");
 let solution=solve(0),mask=solution.mask;
 const cuts=[];for(let r=1;r<=4;r++)if(eqs[r].every(([i,j,e])=>Number(((mask>>BigInt(i))^(mask>>BigInt(j)))&1n)===e))cuts.push(r);
 req(cuts.length===0,"explicit fixed solution is prime by exact ideal criterion");
 req(solution.count===(1n<<34n),"fixed dimension34");
 return {status:"PASS",uniformizer_depth:5,vertices:3125,bank_bits:74,extra_scalar:26,extra_scalar_order:5,all_scalar_fixed_states:solution.count.toString(),prime_scalar_fixed_states:primeFixed.toString(),fixed_cut_intersections:fixedCounts,prime_witness_mask_hex:mask.toString(16),prime_witness_cuts:cuts,meaning:"Extra odd normalizer elements can fix prime banks beyond depth4; do not extrapolate stabilizer H."};
}

console.log(JSON.stringify(verifyDepth5Hostile(), null, 2));
