function verifyDualNumberFamily(){
 const req=(x,m)=>{if(!x)throw Error(m)},out=[];
 for(const p of [3,5,7]){
 const q=p*p,N=q*q,mod=x=>(x%q+q)%q,idx=(u,v)=>mod(u)*q+mod(v),xy=x=>[Math.floor(x/q),x%q],add=(x,y)=>{let[u,v]=xy(x),[a,b]=xy(y);return idx(u+a,v+b)},g=x=>{let[u,v]=xy(x);return idx(u,v+u)},delta=x=>idx(0,Math.floor(x/q)),key=H=>[...H].sort((a,b)=>a-b).join(","),groups=new Map();
 function put(H){let S=new Set(H);groups.set(key(S),S)}
 put([0]);put(Array.from({length:N},(_,i)=>i));
 for(let x=1;x<N;x++){let[u,v]=xy(x),order=(u%p||v%p)?q:p;put(Array.from({length:order},(_,j)=>idx(j*u,j*v)))}
 put(Array.from({length:q},(_,i)=>idx(p*Math.floor(i/p),p*(i%p))));
 for(let a=0;a<p;a++)for(let b=0;b<p;b++){if(a+b===0)continue;put(Array.from({length:N},(_,i)=>i).filter(x=>{let[u,v]=xy(x);return (a*u+b*v)%p===0}))}
 req(groups.size===p*p+3*p+5,"complete additive subgroup classification");
 const ideals=[...groups.values()].filter(H=>H.size>1&&H.size<N&&[...H].every(x=>H.has(delta(x))));
 const expected=new Set();
 expected.add(key(new Set(Array.from({length:p},(_,j)=>idx(0,p*j)))));
 expected.add(key(new Set(Array.from({length:q},(_,j)=>idx(p*Math.floor(j/p),p*(j%p))))));
 expected.add(key(new Set(Array.from({length:p*q},(_,j)=>idx(p*Math.floor(j/q),j%q)))));
 for(let t=0;t<p;t++)expected.add(key(new Set(Array.from({length:q},(_,j)=>idx(p*t*j,j)))));
 req(ideals.length===p+3&&expected.size===p+3&&ideals.every(H=>expected.has(key(H))),"explicit proper ideal list");
 const seen=new Set(),O=[],lookup=new Map();
 for(let x=1;x<N;x++){if(seen.has(x))continue;let y=x,o=[];do{o.push(y);seen.add(y);y=g(y)}while(y!==x);o.forEach(z=>lookup.set(z,O.length));O.push(o)}
 let bits=0,vars=new Map(),used=new Set();
 for(let i=0;i<O.length;i++){if(used.has(i))continue;let[u,v]=xy(O[i][0]),j=lookup.get(idx(-u,-v)),bit=bits++;req(i!==j,"skew orbit pairs");used.add(i);used.add(j);O[i].forEach(x=>vars.set(x,[bit,0]));O[j].forEach(x=>vars.set(x,[bit,1]))}
 req(bits===(3*p*p-2*p-1)/2,"orbit exponent formula");
 const eqs=ideals.map(H=>{let E=new Map();for(let x=1;x<N;x++){if(H.has(x))continue;let[i,a]=vars.get(x);for(let h of H){let[j,b]=vars.get(add(x,h)),row=[i,j,a^b];E.set(row.join(","),row)}}return [...E.values()]});
 function count(mask){let parent=Array.from({length:bits},(_,i)=>i),par=Array(bits).fill(0),components=bits;
 function find(x){if(parent[x]!==x){let y=parent[x];parent[x]=find(y);par[x]^=par[y]}return parent[x]}
 for(let r=0;r<eqs.length;r++)if(mask>>r&1)for(let[i,j,e]of eqs[r]){let a=find(i),b=find(j),v=par[i]^par[j]^e;if(a===b){if(v)return 0n}else{parent[a]=b;par[a]=v;components--}}
 return 1n<<BigInt(components)}
 let prime=0n;for(let m=0;m<2**ideals.length;m++){let f=count(m),par=0;for(let z=m;z;z>>=1)par^=z&1;prime+=(par?-1n:1n)*f}
 const t=1n<<BigInt((p-1)/2),B=1n<<BigInt((p*p-1)/2),C=1n<<BigInt((2*p*p-p-1)/2),D=1n<<BigInt(bits),formula=D-2n*t*C-B*B+3n*t*t*B-t**4n;
 req(prime===formula,"typed general formula versus complete ideal intersection census");
 out.push({p,vertices:N,additive_subgroups:groups.size,proper_ideals:ideals.length,orbit_pair_bits:bits,all_states:D.toString(),prime_states:prime.toString(),composite_states:(D-prime).toString(),ideal_subsets_checked:2**ideals.length});
 }
 return {status:"PASS",ring:"(Z/p^2)[epsilon]/epsilon^2",clock:"1+epsilon",formula:"P=D-2tC-B^2+3t^2B-t^4; t=2^((p-1)/2), B=2^((p^2-1)/2), C=2^((2p^2-p-1)/2), D=2^((3p^2-2p-1)/2)",cases:out};
}

console.log(JSON.stringify(verifyDualNumberFamily(), null, 2));
