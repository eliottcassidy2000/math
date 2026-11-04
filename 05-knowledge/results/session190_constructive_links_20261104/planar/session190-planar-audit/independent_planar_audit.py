"""Independent finite checks of planar interfaces, probes, and graph counts.
No imports from the producer's code. Python standard library only.
"""
import itertools as it
from collections import Counter
from functools import lru_cache
from pathlib import Path
import json, math, hashlib

HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'session190-planar'
report={}

def nc(xs):
    if not xs:
        yield (); return
    a=xs[0]
    # The two intervals between paired endpoints must each have even size.
    for j in range(1,len(xs),2):
        for left in nc(xs[1:j]):
            for right in nc(xs[j+1:]):
                yield tuple(sorted(((a,xs[j]),)+left+right))

def states(k):
    out=[]
    for m in range(1,k//2+1):
        for D in it.combinations(range(k),2*m):
            rest=tuple(i for i in range(k) if i not in D)
            for values in it.product((0,2),repeat=len(rest)):
                d=[1]*k
                for i,x in zip(rest,values): d[i]=x
                out.extend((tuple(d),M) for M in nc(D))
    return sorted(out)

def compatible(s,t):
    if tuple(2-x for x in s[0])!=t[0]: return 0
    D=tuple(i for i,x in enumerate(s[0]) if x==1)
    parents={i:i for i in D}
    def root(i):
        while parents[i]!=i:i=parents[i]
        return i
    for a,b in s[1]+t[1]:parents[root(a)]=root(b)
    return int(len({root(i) for i in D})==1)

def gf2rank(matrix):
    basis={}
    for row in matrix:
        bits=sum(x<<j for j,x in enumerate(row))
        while bits:
            p=bits.bit_length()-1
            if p not in basis: basis[p]=bits;break
            bits^=basis[p]
    return len(basis)

def determinant_small(A):
    n=len(A);total=0
    for p in it.permutations(range(n)):
        inversions=sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        total+=(-1 if inversions%2 else 1)*math.prod(A[i][p[i]] for i in range(n))
    return total

def build_probe(s):
    d,M=s;k=len(d);D=[i for i in range(k) if d[i]==1]
    branches={}
    for a,b in zip(D,D[1:]+D[:1]):
        branch=[a];i=(a+1)%k
        while i!=b:
            if d[i]==2:branch.append(i)
            i=(i+1)%k
        branches[a]=branch
    E=[];v=k
    for a,b in M:
        path=branches[a]+branches[b][::-1]
        for x,y in zip(path,path[1:]):
            E.extend(((x,v),(y,v)));v+=1
    return v,E

def enumerate_profile(n,E,k):
    out=Counter()
    for mask in range(1<<len(E)):
        deg=[0]*n;parent=list(range(n));cycles=False
        def root(i):
            while parent[i]!=i:i=parent[i]
            return i
        for j,(a,b) in enumerate(E):
            if not mask>>j&1:continue
            deg[a]+=1;deg[b]+=1
            u,v=root(a),root(b)
            if u==v:cycles=True;break
            parent[u]=v
        if cycles or any(x>2 for x in deg[:k]) or any(x!=2 for x in deg[k:]):continue
        endpoint_groups={}
        for i in range(k):
            if deg[i]==1:endpoint_groups.setdefault(root(i),[]).append(i)
        M=tuple(sorted(tuple(v) for v in endpoint_groups.values()))
        if not M or any(len(p)!=2 for p in M):continue
        out[(tuple(deg[:k]),M)]+=1
    return out

rows=[];probe_count=0
for k in range(2,7):
    S=states(k);K=[[compatible(s,t) for t in S] for s in S]
    distribution=Counter(map(sum,K))
    expected=[1,6,26,100,365][k-2]
    if len(S)!=expected:raise RuntimeError('state count')
    if k<=5 and set(distribution)!={1}:raise RuntimeError('not permutation')
    if gf2rank(K)!=(expected if k<=5 else expected-1):raise RuntimeError('wrong rank')
    for s in S:
        n,E=build_probe(s)
        if n>2*k-1 or n<=k or any(a<k and b<k for a,b in E):raise RuntimeError('probe bound')
        if enumerate_profile(n,E,k)!=Counter({s:1}):raise RuntimeError('probe profile')
        probe_count+=1
    rows.append({'k':k,'states':len(S),'row_sum_distribution':dict(distribution),'gf2_rank':gf2rank(K)})
M=[((0,1),(2,3),(4,5)),((0,1),(2,5),(3,4)),((0,3),(1,2),(4,5)),((0,5),(1,2),(3,4)),((0,5),(1,4),(2,3))]
K6=[[compatible(((1,)*6,a),((1,)*6,b)) for b in M] for a in M]
if determinant_small(K6)!=-2:raise RuntimeError('determinant')
if [sum(K6[i][j] for j in (1,2,4))%2 for i in range(5)]!=[0]*5:raise RuntimeError('kernel')
six_states=states(6)
full_six=[[compatible(s,t) for t in six_states] for s in six_states]
block_ids={i for i,s in enumerate(six_states) if s[0]==(1,)*6}
outside=[i for i in range(len(six_states)) if i not in block_ids]
permutation={i:next(j for j,x in enumerate(full_six[i]) if x) for i in outside}
if set(permutation.values())!=set(outside) or any(sum(full_six[i])!=1 for i in outside):raise RuntimeError('permutation complement')
unseen=set(outside);cycles_count=0
while unseen:
    seed=next(iter(unseen));i=seed;cycles_count+=1
    while i in unseen:
        unseen.remove(i);i=permutation[i]
    if i!=seed:raise RuntimeError('permutation cycle')
complement_det=(-1)**(len(outside)-cycles_count)
full_det=complement_det*determinant_small(K6)
report['interfaces']={'rows':rows,'independent_probes':probe_count,'six_block':K6,'six_block_det':determinant_small(K6),'full_six_det':full_det,'permutation_complement_size':len(outside),'permutation_complement_cycles':cycles_count,'full_six_rational_rank':len(outside)+len(K6),'full_six_gf2_rank':gf2rank(full_six)}

# Held-Karp recurrence on visited subsets, independent of vertex-order backtracking.
def cycle_count(vertices,edges):
    V=tuple(sorted(vertices));n=len(V)
    if n<3:return 0
    ix={v:i for i,v in enumerate(V)};adj=[0]*n
    for a,b in edges:
        i,j=ix[a],ix[b];adj[i]|=1<<j;adj[j]|=1<<i
    full=(1<<n)-1
    @lru_cache(None)
    def paths(mask,last):
        if mask==full:return int(bool(adj[last]&1))
        ways=0;choices=adj[last]&~mask
        while choices:
            b=choices&-choices;choices-=b
            ways+=paths(mask|b,b.bit_length()-1)
        return ways
    result=paths(1,0)
    if result%2:raise RuntimeError('cycle orientation parity')
    return result//2

data=json.loads((SOURCE/'RESULTS.json').read_text())
for case in data['physical_gluing_cases']:
    V=set(case['A']['vertices'])|set(case['B']['vertices'])
    E={tuple(e) for e in case['A']['edges']+case['B']['edges']}
    if cycle_count(V,E)!=case['cycles']:raise RuntimeError('physical cycle mismatch')
report['independent_physical_cycle_counts']=len(data['physical_gluing_cases'])

# All-start subset DP records the set of endpoints for every reachable visited set.
def traceable(V,E):
    V=tuple(sorted(V));n=len(V);ix={v:i for i,v in enumerate(V)};adj=[0]*n
    for a,b in E:
        i,j=ix[a],ix[b];adj[i]|=1<<j;adj[j]|=1<<i
    current={1<<i:1<<i for i in range(n)}
    for length in range(1,n):
        following={}
        for mask,endpoints in current.items():
            while endpoints:
                bit=endpoints&-endpoints;endpoints-=bit
                choices=adj[bit.bit_length()-1]&~mask
                while choices:
                    nxt=choices&-choices;choices-=nxt
                    nm=mask|nxt;following[nm]=following.get(nm,0)|nxt
        if not following:return False
        current=following
    return bool(current)

square=[]
for row in data['square_deletion_profiles']:
    n=row['n'];V=set(range(1,n+1));E={(i,j) for i in V for j in V if i<j and math.isqrt(i+j)**2==i+j}
    if traceable(V,E):raise RuntimeError('square base')
    restored=[]
    for v in V:
        if traceable(V-{v},{e for e in E if v not in e}):restored.append(v)
    if set(restored)!=set(map(int,row['restoring_deletions'])):raise RuntimeError('square cards')
    square.append({'n':n,'restoring_deletions':restored})
report['independent_square_subset_dp']={'base_graphs':6,'deletion_cards':124,'results':square}

st=json.loads((SOURCE/'STELLATION.json').read_text())
V=set(st['vertices']);E={tuple(e) for e in st['edges']}
actual={str(v):cycle_count(V-{v},{e for e in E if v not in e}) for v in V}
if cycle_count(V,E)!=0 or actual!=st['deletion_cycle_counts']:raise RuntimeError('stellation counts')
report['stellation']={'cycles':0,'deletion_cycles':actual,'traceable':traceable(V,E)}
report['source_hashes']={p:hashlib.sha256((SOURCE/p).read_bytes()).hexdigest() for p in ('PROOF.md','verify.py','independent_kernel.py','stellation.py','RESULTS.json','KERNEL_REPLAY.json','STELLATION.json')}
report['status']='PASS'
(HERE/'INDEPENDENT_RESULTS.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
