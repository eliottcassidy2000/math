"""Independent bounded order-seven audit. No imports from producer code.
Uses a Gray sweep of all 2^21 orientations and path-edge predicate counting.
"""
from itertools import combinations, permutations
from collections import Counter, deque
from pathlib import Path
from fractions import Fraction
import hashlib, json, subprocess, sys

ROOT=Path(__file__).resolve().parent.parent
SOURCE=ROOT/'session23b-fano'
def check(value, message):
    if not value: raise RuntimeError(message)

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
manifest=json.loads((SOURCE/'manifest.json').read_text())
check(all(sha(SOURCE/p)==h for p,h in manifest['files'].items()), 'frozen manifest hashes')

V=tuple(range(7)); PAIRS=list(combinations(V,2)); EDGE={p:i for i,p in enumerate(PAIRS)}
ALL=(1<<21)-1

def beats(m,u,v):
    if u<v:return bool(m & (1<<EDGE[u,v]))
    return not bool(m & (1<<EDGE[v,u]))

def encode_arcs(arcs):
    arcs=set(arcs)
    check(len(arcs)==21,'21 arcs')
    check(all(((u,v) in arcs)^((v,u) in arcs) for u,v in PAIRS),'tournament arcs')
    return sum(1<<i for i,(u,v) in enumerate(PAIRS) if (u,v) in arcs)

# A different exhaustive generator: traverse every possible orientation once.
scores=list(range(7)); bad=sum(d!=3 for d in scores); regular=[]; old=0
for k in range(1<<21):
    m=k^(k>>1)
    if k:
        changed=m^old; i=changed.bit_length()-1; u,v=PAIRS[i]
        bad-=(scores[u]!=3)+(scores[v]!=3)
        d=1 if m&changed else -1
        scores[u]+=d; scores[v]-=d
        bad+=(scores[u]!=3)+(scores[v]!=3)
    if bad==0: regular.append(m)
    old=m
check(len(regular)==2640 and len(set(regular))==2640,'Gray regular census')
regular=set(regular)

# Hamilton counts by 5040 literal path predicates, not dynamic programming.
path_conditions=[]
for p in permutations(V):
    required=0; selected=0
    for u,v in zip(p,p[1:]):
        bit=1<<EDGE[tuple(sorted((u,v)))]; selected|=bit
        if u<v:required|=bit
    path_conditions.append((selected,required))
h={m:sum((m&s)==r for s,r in path_conditions) for m in regular}
check(Counter(h.values())=={189:240,171:1680,175:720},'Hamilton classes')

def relabel(m,p):
    return encode_arcs((p[u],p[v]) if beats(m,u,v) else (p[v],p[u]) for u,v in PAIRS)
reps={H:min(m for m in regular if h[m]==H) for H in (189,171,175)}
orbits={H:{relabel(m,p) for p in permutations(V)} for H,m in reps.items()}
check(all(O=={m for m in regular if h[m]==H} for H,O in orbits.items()),'three full labelled isomorphism classes')

# Directly construct the full simple triangle-reversal graph.
TRIS=list(combinations(V,3)); tri_masks={t:sum(1<<EDGE[p] for p in combinations(t,2)) for t in TRIS}
graph={}
for m in regular:
    cyclic=[t for t in TRIS if beats(m,t[0],t[1])==beats(m,t[1],t[2])==beats(m,t[2],t[0])]
    graph[m]=tuple(m^tri_masks[t] for t in cyclic)
    check(len(cyclic)==14 and set(graph[m])<=regular,'14 regular triangle neighbors')
order=(189,171,175)
rows={H:set(tuple(Counter(h[n] for n in graph[m])[K] for K in order) for m in regular if h[m]==H) for H in order}
Q=((0,14,0),(2,9,3),(0,7,7))
check(all(rows[H]=={Q[i]} for i,H in enumerate(order)), 'equitable full quotient')
check(all(m in graph[n] for m in regular for n in graph[m]),'undirected full graph')

# Construct Fano by XOR-zero triples, then associate dual normals.
X=tuple(range(1,8))
fano={frozenset(t) for t in combinations(X,3) if t[0]^t[1]^t[2]==0}
lines={a:tuple(sorted(next(L for L in fano if all((a&x).bit_count()%2==0 for x in L)))) for a in X}
check(len(fano)==7 and Counter(p for L in fano for p in combinations(sorted(L),2))==Counter(combinations(X,2)),'Fano exact STS')

def cube_tournament(word):
    arcs=[]
    for a,L in lines.items():
        x,y,z=(v-1 for v in L)
        cycle=[(x,y),(y,z),(z,x)]
        arcs.extend(cycle if word&(1<<(a-1)) else [(v,u) for u,v in cycle])
    return encode_arcs(arcs)
cube={w:cube_tournament(w) for w in range(128)}; reverse_cube={m:w for w,m in cube.items()}
check(len(reverse_cube)==128 and set(cube.values())<=regular,'128 embedded presentations')

def syn(word):
    z=0
    for a in X:
        if (word^96)&(1<<(a-1)):z^=a
    return z
check(all(h[m]==(189 if syn(w)==0 else 171) for w,m in cube.items()),'affine Hamming Paley locus')
check(Counter(syn(w) for w in cube)=={s:16 for s in range(8)},'16 fibers')
check(Counter((syn(w),w.bit_count()%2) for w in cube)=={(s,p):8 for s in range(8) for p in (0,1)},'8 parity fibers')
for w,m in cube.items():
    internal=set(graph[m]) & set(reverse_cube)
    check(internal=={cube[w^(1<<i)] for i in range(7)},'induced cube and 7+7 neighbors')
    check({(syn(w^(1<<i)),(w^(1<<i)).bit_count()%2) for i in range(7)}=={(s,1-w.bit_count()%2) for s in range(8) if s!=syn(w)},'crown graph quotient')
    dist={m:0}; queue=deque([m])
    while queue:
        u=queue.popleft()
        for v in graph[u]:
            if v not in dist:dist[v]=dist[u]+1;queue.append(v)
    check(len(dist)==2640,'connected full graph')
    check(all(dist[n]==(w^v).bit_count() for v,n in cube.items()),'all 16384 cube-pair distances')
volts=set()
for a,b in combinations(X,2):
    c=a^b; d=(1<<(a-1))|(1<<(b-1))|(1<<(c-1))
    check(d.bit_count()==3 and syn(96^d)==0,'nonzero triangle holonomy')
    volts.add(d)
check(len(volts)==7,'seven holonomy supports')

# Tetrahedral edge-face incidence and explicit oriented boundary.
E=list(combinations(range(4),2))
evec={e:sum(1<<v for v in e if v!=3) for e in E}
faces=[(1,2,3),(0,3,2),(0,1,3),(0,2,1)]
facial=[]; fixed={}
for face in faces:
    cyc=tuple(evec[tuple(sorted((u,v)))] for u,v in zip(face,face[1:]+face[:1]))
    facial.append(frozenset(cyc)); a=next(a for a,L in lines.items() if set(L)==set(cyc))
    x,y,z=lines[a];fixed[a]=int(cyc in ((x,y,z),(y,z,x),(z,x,y)))
opp=[(evec[e],evec[f]) for e,f in combinations(E,2) if set(e).isdisjoint(f)]
check(set(facial)|{frozenset((7,*p)) for p in opp}==fano,'tetrahedral Fano map')
check(fixed=={1:0,2:1,4:0,7:1},'coherent orientation bits')
coherent=[w for w in cube if all((w>>(a-1))&1==b for a,b in fixed.items())]
check(coherent==[66,70,82,86,98,102,114,118],'eight coherent masks')
check(all(h[cube[w]]==171 for w in coherent),'coherent H171')
for w in coherent:
    m=cube[w]
    for x,y in opp:
        others=set(range(1,7))-{x,y}
        check(all(beats(m,x-1,z-1)==beats(m,y-1,z-1) for z in others),'all-completion opposite twins')
        check(sum(beats(m,x-1,z) and beats(m,y-1,z) for z in V if z not in (x-1,y-1))==2,'pair has two common outneighbors')
repairs={a:sorted(w^(1<<(a-1)) for w in coherent if h[cube[w^(1<<(a-1))]]==189) for a in fixed}
check(repairs=={1:[83,103],2:[84,96],4:[74,126],7:[6,50]},'one-face seam repairs')

# Group and sheet audit, using an independently constructed arithmetic Paley7.
paley=encode_arcs((u,v) if (v-u)%7 in (1,2,4) else (v,u) for u,v in PAIRS)
aut=[p for p in permutations(V) if relabel(paley,p)==paley]
anti=[p for p in permutations(V) if relabel(paley,p)==(paley^ALL)]
check(len(aut)==len(anti)==21,'21 automorphisms and 21 anti-automorphisms')

def sign(p):return (-1)**sum(p[i]>p[j] for i,j in PAIRS)
fano0={frozenset(x-1 for x in L) for L in fano}
fano_aut=[p for p in permutations(V) if {frozenset(p[x] for x in L) for L in fano0}==fano0]
check(len(fano_aut)==168 and all(sign(p)==1 for p in fano_aut),'Fano 168 even group')
check(all(sign(p)==-1 for p in anti),'all Paley anti odd')
check(sum(p[6]==6 for p in fano_aut)==24,'tetrahedral point stabilizer24')
cyclic=[t for t in TRIS if beats(paley,t[0],t[1])==beats(paley,t[1],t[2])==beats(paley,t[2],t[0])]
sheets=[]
for blocks in combinations(cyclic,7):
    if Counter(p for B in blocks for p in combinations(B,2))==Counter(PAIRS):sheets.append(frozenset(frozenset(B) for B in blocks))
check(len(sheets)==2 and sheets[0].isdisjoint(sheets[1]),'two disjoint Paley Fano sheets')
check(all(frozenset(frozenset(p[x] for x in B) for B in sheets[i])==sheets[1-i] for p in anti for i in (0,1)),'every anti exchanges sheets')

subsets={}
for name,G in [('aut',aut),('extended',aut+anti),('fano',fano_aut)]:
    profile={}
    for k in range(8):
        profile[k]=dict(Counter(sum({g[x] for x in S}==set(S) for g in G) for S in combinations(V,k)))
    subsets[name]=profile
check([subsets['aut'][k].get(1,0) for k in range(8)]==[0,0,21,21,21,21,0,0],'84 ordinary regular subsets')
check(all(1 not in subsets[n][k] for n in ('extended','fano') for k in range(8)),'no extended or Fano regular subsets')
paths=[]
for k in range(1,8):
    P=[p for p in permutations(V,k) if all(beats(paley,u,v) for u,v in zip(p,p[1:]))]
    canonical=lambda p:min(tuple(g[x] for x in p) for g in aut)
    O={canonical(p) for p in P}
    inv=anti[0]
    fixed_orbits=sum(canonical(tuple(inv[x] for x in p[::-1]))==p for p in O)
    paths.append([k,len(P),len(O),fixed_orbits,(len(O)+fixed_orbits)//2])
check(paths==[[1,7,1,1,1],[2,21,1,1,1],[3,63,3,1,2],[4,147,7,3,5],[5,273,13,3,8],[6,315,15,3,9],[7,189,9,3,6]],'all marked-path counts')

# Actual graph walks, not merely quotient-matrix powers.
full_weights={paley:1}; weighted=[]
for t in range(13):
    weighted.append(sum(h[m]*c for m,c in full_weights.items()))
    nxt=Counter()
    for m,c in full_weights.items():
        for n in graph[m]:nxt[n]+=c
    full_weights=nxt
num=[189,-630,-1470]; den=[1,-16,14,196]
check(all(sum(den[j]*weighted[t-j] for j in range(min(3,t)+1))==(num[t] if t<3 else 0) for t in range(13)),'actual graph weighted rational series through12')
weights={96:1}; cube_returns=[]
for t in range(13):
    locus=sum(c for w,c in weights.items() if syn(w)==0)
    check(8*locus==7**t+7*(-1)**t,'actual cube locus return formula')
    cube_returns.append([t,locus,weights.get(96,0)])
    nxt=Counter()
    for w,c in weights.items():
        for i in range(7):nxt[w^(1<<i)]+=c
    weights=nxt
check(cube_returns[3]==[3,42,0],'42 locus returns but zero exact returns')
# Exact characteristic coefficients from trace, trace-square, determinant.
trace=sum(Q[i][i] for i in range(3))
trace2=sum(Q[i][j]*Q[j][i] for i in range(3) for j in range(3))
det=sum((-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))*Q[0][p[0]]*Q[1][p[1]]*Q[2][p[2]] for p in permutations(range(3)))
check([1,-trace,(trace*trace-trace2)//2,-det]==[1,-16,14,196],'characteristic polynomial independent coefficients')
pi=[Fraction(len(orbits[H]),2640) for H in order]
check(all(sum(pi[i]*Q[i][j] for i in range(3))==14*pi[j] for j in range(3)),'class stationary law')
check(sum(pi[i]*order[i] for i in range(3))==Fraction(1911,11),'stationary mean')
# Solve 5a-3b=14 and -7a+7b=14 by exact elimination.
hitA=Fraction(14+3*2,2);hitB=hitA+2
check((hitA,hitB)==(10,12),'first hit means')

# Producer replay is supplemental and done only after independent checks above.
replays={}
for script,receipt in [('check.py','receipt.json'),('global_triangle_check.py','global-receipt.json')]:
    expected=(SOURCE/receipt).read_bytes()
    for optimize in (False,True):
        command=[sys.executable]+(['-O'] if optimize else [])+[str(SOURCE/script)]
        proc=subprocess.run(command,check=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        check(proc.stdout==expected, f'replay bytes {script} optimized={optimize}')
        replays[f'{script}:optimized={optimize}']=hashlib.sha256(proc.stdout).hexdigest()
check(all(sha(SOURCE/p)==h for p,h in manifest['files'].items()), 'frozen hashes unchanged after audit')
result={
 'status':'PASS', 'scope':'Exact labelled n=7 only; all 2^21 orientations, all 2640 regular tournaments; no larger-order inference.',
 'independence':'Gray-code global generation; direct 5040 Hamilton-path predicates per regular tournament; all full-graph BFS distances for 128 cube sources; actual-state walk distributions.',
 'source_manifest_sha256':sha(SOURCE/'manifest.json'), 'verified_source_hashes':manifest['files'],
 'all_orientation_count':1<<21,'regular_count':len(regular),'H_class_sizes':dict(Counter(h.values())),
 'class_order':order,'neighbor_matrix':Q,'cube_distance_pairs_checked':128*128,'full_graph_directed_edges':sum(map(len,graph.values())),
 'tetrahedron_boundary_masks':coherent,'single_face_repairs':repairs,'triangle_holonomy_supports':sorted(volts),
 'fano_sheet_count_in_Paley':len(sheets),'paley_aut_order':len(aut),'paley_anti_count':len(anti),'fano_aut_order':len(fano_aut),
 'path_profile_columns':['k','labelled_paths','Aut_orbits','converse_fixed','converse_merged'],'path_profile':paths,
 'weighted_actual_graph_walks':weighted,'cube_returns_columns':['t','to_Paley_locus','to_same_presentation'],'cube_returns':cube_returns,
 'characteristic_polynomial_high_first':[1,-16,14,196],'stationary_probabilities':list(map(str,pi)),'stationary_mean_H':'1911/11','first_hit_means':[10,12],
 'producer_replays':replays, 'findings':'No mathematical error found. Frozen artifacts reproduced; inheritance and operation-relative limitations retained.'}
print(json.dumps(result,sort_keys=True,indent=2))
