function verify625(){
 const req=(x,m)=>{if(!x)throw Error(m)};
 const mod=x=>(x%25+25)%25, idx=(u,v)=>mod(u)*25+mod(v), xy=x=>[Math.floor(x/25),x%25];
 const mul=(x,y)=>{let[a,b]=xy(x),[c,d]=xy(y);return idx(a*c+5*b*d,a*d+b*c)};
 const g=idx(16,14), gp=[idx(1,0)]; for(let j=1;j<=25;j++)gp.push(mul(gp[j-1],g));req(gp[25]===25&&new Set(gp.slice(0,25)).size===25,"g order25");
 const seen=new Set(),orbits=[],lookup=new Map();
 for(let x=1;x<625;x++){if(seen.has(x))continue;let o=[],y=x;do{o.push(y);seen.add(y);y=mul(g,y)}while(y!==x);o.forEach(z=>lookup.set(z,orbits.length));orbits.push(o)}
 const vars=new Map(),pairs=[],used=new Set();
 for(let i=0;i<orbits.length;i++){if(used.has(i))continue;let[u,v]=xy(orbits[i][0]),j=lookup.get(idx(-u,-v));req(j!==i,"no selfnegative");let bit=pairs.length;pairs.push([i,j]);used.add(i);used.add(j);orbits[i].forEach(x=>vars.set(x,[bit,0]));orbits[j].forEach(x=>vars.set(x,[bit,1]))}
 req(pairs.length===24,"24 bits");
 const action=(unit,b)=>pairs.map(([i])=>{let[r,s]=xy(orbits[i][0]);let[j,e]=vars.get(mul(unit,idx(r,b*s)));return j+32*e});
 const acts=new Map(),matrices=new Set(),exponents={};
 for(let u=0;u<25;u++)for(let v=0;v<25;v++){if(u%5===0)continue;let unit=idx(u,v);for(let b=0;b<25;b++){if(b%5!==1&&b%5!==4)continue;
 let gg=idx(16,14*b),j=gp.indexOf(gg);req(j>0&&j<25&&j%5!==0,"clock normalized");exponents[j]=(exponents[j]||0)+1;
 matrices.add([u,mod(5*v*b),v,mod(u*b)].join(","));let a=action(unit,b);acts.set(a.join(","),a);
 }}
 req(matrices.size===5000&&acts.size===200,"normalizer and effective size");
 const eqs=[null];
 for(let r=1;r<=3;r++){let equations=new Map(),U=5**Math.ceil(r/2),V=5**Math.floor(r/2),gens=r===1?[[0,1],[5,0]]:r===2?[[5,0],[0,5]]:[[0,5]];
 for(let x=1;x<625;x++){let[u,v]=xy(x);if(u%U===0&&v%V===0)continue;let[i,e]=vars.get(x);for(let[a,b]of gens){let[j,f]=vars.get(idx(u+a,v+b));let q=[i,j,e^f];equations.set(q.join(","),q)}}eqs.push([...equations.values()])}
 function count(a,cuts){
 let parent=Array.from({length:24},(_,i)=>i),parity=Array(24).fill(0),components=24;
 function find(x){if(parent[x]!==x){let p=parent[x];parent[x]=find(p);parity[x]^=parity[p]}return parent[x]}
 function join(i,j,e){let ri=find(i),rj=find(j),v=parity[i]^parity[j]^e;if(ri===rj)return v===0;parent[ri]=rj;parity[ri]=v;components--;return true}
 for(let i=0;i<24;i++)if(!join(i,a[i]%32,a[i]>>5))return 0;
 for(let r=1;r<=3;r++)if(cuts>>(r-1)&1)for(let[i,j,e]of eqs[r])if(!join(i,j,e))return 0;
 return 2**components;
 }
 let hist={},primeHist={},sum=0,sumPrime=0,rows=[];
 for(let a of acts.values()){let fixed=[],pf=0;for(let c=0;c<8;c++){let f=count(a,c);fixed.push(f);let parity=0;for(let z=c;z;z>>=1)parity^=z&1;pf+=(parity?-1:1)*f}
 hist[fixed[0]]=(hist[fixed[0]]||0)+1;primeHist[pf]=(primeHist[pf]||0)+1;sum+=fixed[0];sumPrime+=pf;rows.push({action:a,fixed,prime_fixed:pf})}
 req(sum%200===0&&sumPrime%200===0,"Burnside integral");
 const tau=action(idx(1,0),6),unitOnly=action(idx(6,0),21),scalar=action(idx(6,0),1);
 req(count(tau,0)===65536&&count(tau,4)===65536,"tau fixed iff cut3");
 req(count(unitOnly,0)===65536&&count(unitOnly,1)===65536,"unit-only fixed iff cut1");
 req(count(scalar,0)===256&&count(scalar,7)===256,"scalar fixed coherent");
 req(primeHist[16646400]===1&&primeHist[0]===199,"all prime effective stabilizers trivial");
 req(sum/200===86528&&sumPrime/200===83232,"corrected orbit counts");
 return {status:"PASS",normalizer_order:5000,effective_order:200,normalizer_exponents:exponents,fixed_histogram:hist,prime_fixed_histogram:primeHist,normalizer_all_orbits:sum/200,prime_full_affine_classes:sumPrime/200,prime_affine_aut_order:15625,composite_normalizer_orbits:(sum-sumPrime)/200,action_rows:rows};
}
function verify125Composite(){
 const req=(v,s)=>{if(!v)throw Error(s)},mod=(x,n)=>(x%n+n)%n,idx=(u,v)=>mod(u,25)*5+mod(v,5),xy=x=>[Math.floor(x/5),x%5];
 function sign(x,bits){let e=x>2?1:0,i=Math.min(x,5-x)-1;return ((bits>>i)&1)^e}
 const states=[],lookup=new Map();
 for(let mask=0;mask<64;mask++){let S=[];for(let x=1;x<125;x++){let[u,v]=xy(x),a=u%5;let yes=a?sign(a,mask&3):v?sign(v,(mask>>2)&3):sign(u/5,(mask>>4)&3);if(yes)S.push(x)}req(S.length===62,"skew state size");let key=S.reduce((n,x)=>n|(1n<<BigInt(x)),0n).toString();lookup.set(key,mask);states.push(S)}
 const act=(M,x)=>{let[a,b,c,d]=M,[u,v]=xy(x);return idx(a*u+5*b*v,c*u+d*v)},g=x=>{let[u,v]=xy(x);return idx(16*u+20*v,4*u+v)};
 const autos=[],normal=[];
 for(let a=0;a<25;a++)for(let b=0;b<5;b++)for(let c=0;c<5;c++)for(let d=1;d<5;d++){if(a%5===0)continue;let M=[a,b,c,d];autos.push(M);let found=false;for(let j=1;j<5;j++){let ok=true;for(let x of [idx(1,0),idx(0,1)]){let y=act(M,x);for(let k=0;k<j;k++)y=g(y);if(act(M,g(x))!==y)ok=false}if(ok)found=true}if(found)normal.push(M)}
 req(autos.length===2000&&normal.length===200,"groups");
 function orbits(group){let unseen=new Set(Array.from({length:64},(_,i)=>i)),out=[];while(unseen.size){let m=Math.min(...unseen),o=new Set();for(let M of group){let key=states[m].reduce((n,x)=>n|(1n<<BigInt(act(M,x))),0n).toString();req(lookup.has(key),"all additive images stay in coherent bank");o.add(lookup.get(key))}for(let x of o){req(unseen.has(x),"disjoint orbit");unseen.delete(x)}out.push({representative:m,size:o.size,masks:[...o].sort((a,b)=>a-b)})}return out}
 let full=orbits(autos),norm=orbits(normal);req(full.length===4&&full.every(o=>o.size===16),"four full affine classes");req(norm.length===8&&norm.every(o=>o.size===8),"eight normalizer classes");
 return {status:"PASS",universe:"64 orientation triples on characteristic flag pR subset R[5] subset R",full_additive_group_order:2000,normalizer_order:200,composite_full_affine_classes:4,composite_normalizer_classes:8,composite_additive_stabilizer:125,composite_affine_aut_order:15625,full_affine_orbits:full,normalizer_orbits:norm};
}

console.log(JSON.stringify({normalizer625:verify625(), composite125:verify125Composite()}, null, 2));
