from itertools import combinations, permutations
from collections import Counter, defaultdict
import json
from fractions import Fraction


def require(q, why):
    if not q: raise AssertionError(why)

V=list(range(1,8))
lines={a:tuple(x for x in V if (a&x).bit_count()%2==0) for a in V}
pairs=list(combinations(V,2))
require(Counter(p for L in lines.values() for p in combinations(L,2))==Counter(pairs),'Fano STS')

def tournament(mask):
    A={u:set() for u in V}
    for a,(x,y,z) in lines.items():
        cyc=(x,y,z) if mask>>(a-1)&1 else (z,y,x)
        for i in range(3): A[cyc[i]].add(cyc[(i+1)%3])
    return A

def hp(A):
    d={(1<<(v-1),v):1 for v in V}
    for s in range(1,128):
        for v in V:
            c=d.get((s,v),0)
            if c:
                for w in A[v]:
                    if not s>>(w-1)&1:
                        d[(s|1<<(w-1),w)]=d.get((s|1<<(w-1),w),0)+c
    return sum(d.get((127,v),0) for v in V)

def syndrome(mask):
    s=0
    for a in V:
        if (mask^96)>>(a-1)&1:s^=a
    return s

rows=[]
for mask in range(128):
    A=tournament(mask);H=hp(A); s=syndrome(mask)
    require(all(len(A[v])==3 for v in V),'regularity')
    require(H==(189 if s==0 else 171),'inherited H census')
    brute=sum(all(p[i+1] in A[p[i]] for i in range(6)) for p in permutations(V))
    require(H==brute,'independent HP')
    neighbor_s=[syndrome(mask^(1<<a)) for a in range(7)]
    require(set(neighbor_s)==set(range(8))-{s},'syndrome K8 covering')
    b=Counter(int(x!=0) for x in neighbor_s)
    require(b==({1:7} if s==0 else {0:1,1:6}),'equitable quotient')
    charsum=sum((-1)**((u&s).bit_count()%2) for u in V)
    require(4*H==693+9*charsum,'Fourier identity')
    rows.append(dict(mask=mask,syndrome=s,H=H,neighbors_by_type=dict(b)))

# Full vertex permutation automorphism/anti groups for Paley, not assumed affine.
Q={1,2,4};P=list(range(7));A={u:{v for v in P if (v-u)%7 in Q} for u in P}
G=[];anti=[]
for p in permutations(P):
    if all(((v in A[u])==(p[v] in A[p[u]])) for u,v in combinations(P,2)):G.append(p)
    if all(((v in A[u])!=(p[v] in A[p[u]])) for u,v in combinations(P,2)):anti.append(p)
require(len(G)==len(anti)==21,'Paley automorphism group')
require(set(G)=={tuple((a*x+b)%7 for x in P) for a in Q for b in P},'affine Aut')
require(set(anti)=={tuple((a*x+b)%7 for x in P) for a in {3,5,6} for b in P},'affine anti')

def subset_profile(group,vertices):
    result={}
    for k in range(8):
        hist=Counter()
        for S in combinations(vertices,k):
            ss=frozenset(S)
            stab=sum(frozenset(p[x] for x in S)==ss for p in group)
            hist[stab]+=1
        result[k]=dict(sorted(hist.items()))
    return result
paley_subset=subset_profile(G,P);extended_subset=subset_profile(G+anti,P)
require([paley_subset[k].get(1,0) for k in range(8)]==[0,0,21,21,21,21,0,0],'84 rigid subsets')
require(all(h.get(1,0)==0 for h in extended_subset.values()),'no converse-rigid subset')
F=[]
line_set=set(map(frozenset,lines.values()))
for seq in permutations(V):
    p=dict(zip(V,seq))
    if {frozenset(p[x] for x in L) for L in lines.values()}==line_set:F.append(p)
require(len(F)==168,'Fano full automorphism group')
fano_subset=subset_profile(F,V)
require(all(h.get(1,0)==0 for h in fano_subset.values()),'no Fano rigid subset')

# Marked directed paths in this one fixed Paley tournament.
path_profile=[]
for k in range(1,8):
    paths=[p for p in permutations(P,k) if all(p[i+1] in A[p[i]] for i in range(k-1))]
    unseen=set(paths);orbits=[]
    while unseen:
        p=min(unseen);O={tuple(g[x] for x in p) for g in G};unseen-=O;orbits.append(O)
    inv=anti[0]
    B=sum({tuple(inv[x] for x in p[::-1]) for p in O}==O for O in orbits)
    fixes=[sum(all(g[p[i]]==p[k-1-i] for i in range(k)) for p in paths) for g in anti]
    require(sum(fixes)==21*B,'twisted Burnside')
    path_profile.append(dict(k=k,labelled_paths=len(paths),ordinary_orbits=len(orbits),converse_fixed=B,merged=(len(orbits)+B)//2,anti_fix_hist=dict(Counter(fixes))))

# Marked tetrahedron edges plus infinity => Fano.
E=list(combinations(range(4),2));edge_vector={e:sum(1<<x for x in e if x<3) for e in E}
require(set(edge_vector.values())==set(range(1,7)),'six K4 column labels')
facial=[frozenset(edge_vector[e] for e in combinations(face,2)) for face in combinations(range(4),3)]
opposition=[frozenset({7,edge_vector[e],edge_vector[f]}) for e,f in combinations(E,2) if set(e).isdisjoint(f)]
require(set(facial+opposition)==line_set,'tetrahedral Fano incidence')
require(sum(p[7]==7 for p in F)==24,'S4 marked-point stabilizer')

# Coherently boundary-oriented tetrahedron faces: induce a cycle on their three edge vertices.
# Oriented simplex boundary uses (-1)^i times increasing face excluding i.
boundarybits={}
for i in range(4):
    face=[v for v in range(4) if v!=i]
    if i%2:face[0],face[1]=face[1],face[0]
    edgecycle=[edge_vector[tuple(sorted((face[j],face[(j+1)%3])))] for j in range(3)]
    normal=next(a for a,L in lines.items() if set(L)==set(edgecycle))
    x,y,z=lines[normal]
    bit=int(any(edgecycle==list(t) for t in [(x,y,z),(y,z,x),(z,x,y)]))
    boundarybits[normal]=bit
consistent=[]
for mask in range(128):
    if all((mask>>(a-1)&1)==b for a,b in boundarybits.items()):
        consistent.append(dict(mask=mask,syndrome=syndrome(mask),H=hp(tournament(mask))))
require(len(consistent)==8,'boundary completion count')

# Direct seam mechanism: opposite edge pairs are twins outside their pair and the apex.
T=tournament(consistent[0]['mask'])
opp_pairs=[(edge_vector[e],edge_vector[f]) for e,f in combinations(E,2) if set(e).isdisjoint(f)]
for x,y in opp_pairs:
    rest=set(range(1,7))-{x,y}
    require(T[x]&rest==T[y]&rest,'opposite-edge twins')
    require(len(T[x]&T[y])==2,'two common out-neighbors obstruct Paley')
seam_repairs={}
for a in boundarybits:
    bits=boundarybits.copy();bits[a]^=1
    candidates=[m for m in range(128) if all((m>>(b-1)&1)==c for b,c in bits.items())]
    repaired=[m for m in candidates if syndrome(m)==0]
    require(len(candidates)==8 and len(repaired)==2,'one-face repair gives two Paley completions')
    seam_repairs[a]=repaired

# Holonomy and count-return certificates. No claim of full-state mixing.
weights={0:1};return_counts=[]
for t in range(13):
    to_locus=sum(c for m,c in weights.items() if syndrome(m^96)==0)
    to_self=weights.get(0,0)
    require(8*to_locus==7**t+7*(-1)**t,'walk return to affine code')
    if t%2: require(to_self==0,'hypercube has no odd closed walk')
    return_counts.append(dict(t=t,return_to_Paley_locus=to_locus,return_to_exact_mask=to_self))
    nxt=Counter()
    for m,c in weights.items():
        for a in range(7):nxt[m^(1<<a)]+=c
    weights=nxt
triangle_voltages=[]
for a,b in combinations(V,2):
    c=a^b
    voltage=(1<<(a-1))^(1<<(b-1))^(1<<(c-1))
    require(voltage and voltage.bit_count()==3,'triangle voltage is nonzero weight three')
    require(syndrome(voltage^96)==0,'triangle voltage belongs to linear Hamming code')
    triangle_voltages.append(voltage)
require(len(set(triangle_voltages))==7,'seven Fano triangle voltages')

def permutation_sign(p,domain):
    vals=[p[x] for x in domain]
    return (-1)**sum(vals[i]>vals[j] for i in range(len(vals)) for j in range(i+1,len(vals)))
require(all(permutation_sign(p,V)==1 for p in F),'Fano permutations even')
require(all(permutation_sign(p,P)==-1 for p in anti),'Paley anti permutations odd')

receipt=dict(status='FINITE-EXACT',universe='All 128 orientations of fixed labelled Fano lines; all 7! vertex permutations for groups and independent HP checks; all subsets and paths of fixed Paley7.',lines=lines,orientation_H_counts=dict(Counter(r['H'] for r in rows)),orientation_rows=rows,paley_aut_size=len(G),paley_anti_size=len(anti),paley_subset_stabilizers=paley_subset,extended_subset_stabilizers=extended_subset,fano_subset_stabilizers=fano_subset,paley_marked_paths=path_profile,tetrahedron_edge_vectors={str(k):v for k,v in edge_vector.items()},tetrahedron_boundary_bits=boundarybits,tetrahedron_boundary_completions=consistent,tetrahedron_single_face_repairs=seam_repairs,line_flip_return_counts=return_counts,distinct_triangle_voltages=sorted(set(triangle_voltages)),all_Fano_automorphisms_even=True,all_Paley_anti_automorphisms_odd=True)
print(json.dumps(receipt,indent=2,sort_keys=True))
