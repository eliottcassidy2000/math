function verify125Graphs(){
 const req=(v,s)=>{if(!v)throw Error(s)},mod=(x,n)=>(x%n+n)%n,idx=(u,v)=>mod(u,25)*5+mod(v,5),xy=x=>[Math.floor(x/5),x%5],add=(x,y)=>{let[a,b]=xy(x),[c,d]=xy(y);return idx(a+c,b+d)},g=x=>{let[u,v]=xy(x);return idx(16*u+20*v,4*u+v)};
 const seen=new Set(),orbits=[],lookup=new Map();
 for(let x=1;x<125;x++){if(seen.has(x))continue;let o=[],y=x;do{o.push(y);seen.add(y);y=g(y)}while(y!==x);o.forEach(z=>lookup.set(z,orbits.length));orbits.push(o)}
 const vars=new Map(),pairs=[],used=new Set();
 for(let i=0;i<orbits.length;i++){if(used.has(i))continue;let[u,v]=xy(orbits[i][0]),j=lookup.get(idx(-u,-v)),bit=pairs.length;req(i!==j,"no selfnegative");pairs.push([i,j]);used.add(i);used.add(j);orbits[i].forEach(x=>vars.set(x,[bit,0]));orbits[j].forEach(x=>vars.set(x,[bit,1]))}
 req(pairs.length===14,"14 bits");
 const act=(M,x)=>{let[a,b,c,d]=M,[u,v]=xy(x);return idx(a*u+5*b*v,c*u+d*v)},normal=[];
 for(let a=0;a<25;a++)for(let b=0;b<5;b++)for(let c=0;c<5;c++)for(let d=1;d<5;d++){if(a%5===0)continue;let M=[a,b,c,d];for(let j=1;j<5;j++){let ok=true;for(let x of [idx(1,0),idx(0,1)]){let y=act(M,x);for(let k=0;k<j;k++)y=g(y);if(act(M,g(x))!==y)ok=false}if(ok)normal.push(M)}}
 req(normal.length===200,"normalizer");
 const actions=normal.map(M=>pairs.map(([i])=>vars.get(act(M,orbits[i][0]))));
 const move=(m,A)=>A.reduce((out,[j,e],i)=>out|((((m>>i)&1)^e)<<j),0);
 const eqs=[];
 for(let r=1;r<=2;r++){let U=5**Math.ceil(r/2),V=5**Math.floor(r/2),eq=new Map(),gens=[[U,0],[0,V]].filter(([u,v])=>idx(u,v)!==0);
 for(let x=1;x<125;x++){let[u,v]=xy(x);if(u%U===0&&v%V===0)continue;let[i,e]=vars.get(x);for(let[a,b]of gens){let[j,f]=vars.get(idx(u+a,v+b)),q=[i,j,e^f];eq.set(q.join(","),q)}}eqs.push([...eq.values()])}
 const hasCut=m=>eqs.some(E=>E.every(([i,j,e])=>(((m>>i)^(m>>j))&1)===e));
 const visited=new Uint8Array(16384),prime=[],composite=[];
 for(let m=0;m<16384;m++){if(visited[m])continue;let O=new Set(actions.map(A=>move(m,A)));for(let x of O){req(!visited[x],"disjoint normalizer orbits");visited[x]=1}let out=hasCut(m)?composite:prime;out.push({mask:m,size:O.size})}
 req(prime.length===408&&prime.every(e=>e.size===40)&&composite.length===8&&composite.every(e=>e.size===8),"complete prime/composite orbit census");
 const pop=x=>{x-=(x>>>1)&0x55555555;x=(x&0x33333333)+((x>>>2)&0x33333333);return (((x+(x>>>4))&0x0f0f0f0f)*0x01010101)>>>24};
 const bit=(row,j)=>(row[j>>>5]>>>(j&31))&1;
 function rowsFor(m){let S=[...vars].filter(([x,[i,e]])=>((m>>i)&1)^e).map(([x])=>x),rows=Array.from({length:125},()=>new Uint32Array(4));req(S.length===62,"connection set");for(let x=0;x<125;x++)for(let s of S){let y=add(x,s);rows[x][y>>>5]|=1<<(y&31)}return rows}
 function refine(rows,parts){
 while(true){let masks=parts.map(C=>{let M=new Uint32Array(4);for(let v of C)M[v>>>5]|=1<<(v&31);return M}),next=[];
 for(let C of parts){let groups=new Map();for(let v of C){let signature=masks.map(M=>pop(rows[v][0]&M[0])+pop(rows[v][1]&M[1])+pop(rows[v][2]&M[2])+pop(rows[v][3]&M[3])),key=signature.join(",");if(!groups.has(key))groups.set(key,{signature,vertices:[]});groups.get(key).vertices.push(v)}
 let ordered=[...groups.values()].sort((a,b)=>{for(let j=0;j<a.signature.length;j++)if(a.signature[j]!==b.signature[j])return a.signature[j]-b.signature[j];return 0});next.push(...ordered.map(x=>x.vertices))}
 if(next.length===parts.length)return next;parts=next;
 }}
 function certificate(rows,h){
 let parts=refine(rows,[[0],Array.from({length:124},(_,i)=>i+1)]),i=-1;
 for(let j=0;j<parts.length;j++)if(parts[j].length>1&&(i<0||parts[j].length<parts[i].length))i=j;
 req(i>=0,"nontrivial root orbit");let C=parts[i],v=Math.min(...C),orbit=new Set([v]),w=h[v];while(w!==v){orbit.add(w);w=h[w]}
 req(C.length===5&&orbit.size===5&&C.every(x=>orbit.has(x)),"chosen cell exactly H orbit");
 let finer=refine(rows,[...parts.slice(0,i),[v],C.filter(x=>x!==v),...parts.slice(i+1)]);req(finer.length===125,"discrete after one individualization");let order=finer.map(C=>C[0]),bytes=[];
 for(let x of order)for(let q=0;q<125;q+=8){let b=0;for(let j=0;j<8&&q+j<125;j++)b|=bit(rows[x],order[q+j])<<j;bytes.push(b)}
 return {key:String.fromCharCode(...bytes),chosen:v,order,initial_cell_sizes:parts.map(C=>C.length).sort((a,b)=>a-b)};
 }
 const h=Array.from({length:125},(_,i)=>g(i)),keys=new Set(),certificates=[];
 let rng=110;const random=()=>{rng^=rng<<13;rng^=rng>>>17;rng^=rng<<5;return rng>>>0};
 for(let pos=0;pos<prime.length;pos++){let m=prime[pos].mask,rows=rowsFor(m),cert=certificate(rows,h);req(!keys.has(cert.key),"full graph keys distinct");keys.add(cert.key);req(cert.initial_cell_sizes.filter(x=>x===1).length===5&&cert.initial_cell_sizes.filter(x=>x===5).length===24,"root partition profile");
 if(pos<3){let perm=Array.from({length:125},(_,i)=>i);for(let i=124;i>1;i--){let j=1+random()%i;[perm[i],perm[j]]=[perm[j],perm[i]]}let rr=Array.from({length:125},()=>new Uint32Array(4)),hh=Array(125);
 for(let x=0;x<125;x++){hh[perm[x]]=perm[h[x]];for(let y=0;y<125;y++)if(bit(rows[x],y))rr[perm[x]][perm[y]>>>5]|=1<<(perm[y]&31)}
 req(certificate(rr,hh).key===cert.key,"nonaffine rooted relabelling control");
 }
 certificates.push({mask:m,chosen_vertex:cert.chosen,discrete_order:cert.order});
 }
 let rejected=false;
 try{certificate(rowsFor(composite[0].mask),h)}catch(error){req(["chosen cell exactly H orbit","discrete after one individualization"].includes(error.message),"unexpected negative-control failure");rejected=true}
 req(rejected,"composite graph must not receive the prime rigidity certificate");
 req(keys.size===408,"complete full graph census");
 return {status:"PASS",composite_certificate_rejected:true,bank_states:16384,prime_states:16320,composite_states:64,prime_normalizer_classes:408,composite_normalizer_classes:8,prime_full_graph_classes:408,prime_root_automorphism_order:5,prime_full_automorphism_order:625,composite_graph_classes_by_product_theorem:1,total_full_graph_classes:409,normalizer_representatives:prime,certificates};
}

console.log(JSON.stringify(verify125Graphs(), null, 2));
