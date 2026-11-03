import sys, itertools, collections, json
sys.path.insert(0, '/tmp/s168_cd')
import check as c

def require(x, msg):
    if not x: raise RuntimeError(msg)

results={}
# Build permutations from dense algebra, never calling sign or full_permutation_group.
for n in range(1,5):
    N=1<<n
    units=[tuple(int(i==j) for i in range(N)) for j in range(N)]
    gens=[]
    for a in range(1,N):
        base=[]
        for x in range(N):
            p=c.multiply(units[a],units[x])
            inds=[j for j,z in enumerate(p) if z]
            require(len(inds)==1,'basis closure')
            j=inds[0]; require(abs(p[j])==1,'unit coefficient')
            base.append(j+N*(p[j]<0))
        gens.append(tuple(base+[q^N for q in base]))
    found={tuple(range(2*N))}; todo=list(found)
    while todo:
        p=todo.pop()
        for g in gens:
            q=tuple(p[g[i]] for i in range(2*N))
            if q not in found: found.add(q); todo.append(q)
    if n>=3:
        actual=set()
        for p in found:
            if p[0]%N: continue
            require(all(p[x]%N==x for x in range(N)), 'kernel fixes every index')
            actual.add(tuple(p[x]//N for x in range(N)))
        half=N//2
        expected=set()
        for bits in itertools.product(range(2),repeat=half):
            if sum(bits)%2: continue
            for top in range(2): expected.add(bits+tuple(b^top for b in bits))
        require(actual==expected,'full actual kernel equals claimed functions')
    results[f'dense_group_n{n}']=len(found)
# Independent coefficient formula: ANF coefficient at M is XOR over subsets of M.
# Physically run all nonzero generators, two signs, and every index through n=7.
gates=0; macros=0
for n in range(2,8):
    N=1<<n
    arc_positions=[]
    order=6*(n+1)
    for block in range(n+1):
        u,v=6*block+1,6*block+2
        pos=sum(order-i-1 for i in range(u))+(v-u-1)
        arc_positions.append(pos)
    zero=c.tournament_bits([0]*(n+1))[1]
    def encode(s,x):
        value=zero^(s<<arc_positions[0])
        for j in range(n): value^=((x>>j)&1)<<arc_positions[j+1]
        return value
    for a in range(1,N):
        row=[c.sign(n,a,x) for x in range(N)]
        coeff=[]
        for mask in range(N):
            val=0; subset=mask
            while True:
                val^=row[subset]
                if subset==0:break
                subset=(subset-1)&mask
            if val:coeff.append(mask)
        require(coeff==c.anf(row),'independent subset coefficient')
        for s in range(2):
            for x in range(N):
                physical=encode(s,x)
                for m in coeff:
                    if x&m==m:
                        before=physical; physical^=1<<arc_positions[0]
                        require((before^physical).bit_count()==1,'one arc sign gate');gates+=1
                for j in range(n):
                    if (a>>j)&1:
                        before=physical;physical^=1<<arc_positions[j+1]
                        require((before^physical).bit_count()==1,'one arc index gate');gates+=1
                require(physical==encode(s^row[x],x^a),'full normal-form physical endpoint')
                macros+=1
results['full_anf_physical_macros']=macros
results['full_anf_physical_arc_flips']=gates
# Fully classify the two-block physical single-flip neighbors by SCCs and seed orbits.
orbits={}
for bit,seed in enumerate([153,185]):
    for perm in itertools.permutations(range(6)):orbits[c.relabel(seed,perm)]=bit

def ordered_sccs(out):
    N=len(out)
    reach=[]
    for v in range(N):
        seen={v}; todo=[v]
        while todo:
            u=todo.pop()
            for w in range(N):
                if (out[u]>>w)&1 and w not in seen:seen.add(w);todo.append(w)
        reach.append(seen)
    remain=set(range(N)); comps=[]
    while remain:
        v=min(remain); block={w for w in remain if w in reach[v] and v in reach[w]}
        comps.append(block);remain-=block
    comps.sort(key=lambda block:-sum(1 for other in comps if other is not block and (out[min(block)]>>min(other))&1))
    return comps

def classify(out):
    blocks=ordered_sccs(out)
    if len(blocks)!=2 or any(len(b)!=6 for b in blocks):return None
    answer=[]
    for block in blocks:
        verts=sorted(block)
        mask=sum(((out[verts[i]]>>verts[j])&1)<<k for k,(i,j) in enumerate(itertools.combinations(range(6),2)))
        if mask not in orbits:return None
        answer.append(orbits[mask])
    return tuple(answer)
flips=0
for bits in itertools.product(range(2),repeat=2):
    out,_=c.tournament_bits(bits)
    destinations=set()
    for u,v in itertools.combinations(range(12),2):
        changed=list(out);changed[u]^=1<<v;changed[v]^=1<<u
        word=classify(changed);flips+=1
        if word is not None:destinations.add(word)
    require(destinations-{bits}=={(bits[0]^1,bits[1]),(bits[0],bits[1]^1)},'induced two-block class cube')
results['physical_two_block_single_flips_classified']=flips
# Direct polynomial evaluation of F from p_j for a range of integer slot counts.
p=[29,131,162,77,15,1]
for r in range(16):
    falling=1;total=0
    for j,pj in enumerate(p,1):falling*=r-j+1;total+=pj*falling
    require(total==r**6+12*r**4+16*r*r,'F polynomial')
results['status']='all independent audit controls passed'
print(json.dumps(results,sort_keys=True,indent=2))
