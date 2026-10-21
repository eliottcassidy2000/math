#!/usr/bin/env python3
"""Exact, optimization-safe controls for Fano supports and flow decorations.
No repository writes, network use, or third-party packages.
"""
from collections import Counter, defaultdict
from itertools import combinations, permutations, product
import json


def require(ok, label):
    if not ok:
        raise RuntimeError(label)


def edges(n):
    return tuple(combinations(range(n), 2))


def orientation(n, mask):
    return tuple((j,i) if mask>>k&1 else (i,j) for k,(i,j) in enumerate(edges(n)))


def directed_cycles(n, arcs):
    es=set(arcs); index={frozenset(e):i for i,e in enumerate(edges(n))}
    out=[]
    for k in range(3,n+1):
        for seq in permutations(range(n),k):
            if seq[0]!=min(seq):
                continue
            cyc=tuple((seq[i],seq[(i+1)%k]) for i in range(k))
            if all(e in es for e in cyc):
                out.append((sum(1<<index[frozenset(e)] for e in cyc),cyc))
    return out


def transport_mask(mask,p,n=4):
    es=edges(n); idx={e:i for i,e in enumerate(es)}
    return sum(1<<idx[tuple(sorted((p[a],p[b])))] for j,(a,b) in enumerate(es) if mask>>j&1)


def permutations_of(arcs,n,anti=False):
    es=set(arcs)
    return [p for p in permutations(range(n)) if {(p[b],p[a]) if anti else (p[a],p[b]) for a,b in arcs}==es]


def strong(n,arcs):
    reach=[{i} for i in range(n)]
    for a,b in arcs: reach[a].add(b)
    for k in range(n):
        for i in range(n):
            if k in reach[i]: reach[i] |= reach[k]
    return all(len(s)==n for s in reach)


def all_vectors(m,bound):
    """All nonnegative vectors with total at most bound, independently of cycles."""
    if m==0:
        yield ()
        return
    for x in range(bound+1):
        for tail in all_vectors(m-1,bound-x):
            yield (x,)+tail


def balanced(arcs,v,n):
    div=[0]*n
    for (a,b),x in zip(arcs,v): div[a]+=x;div[b]-=x
    return not any(div)


def transform_flow(arcs,v,p,anti=False):
    idx={e:i for i,e in enumerate(arcs)}
    out=[0]*len(arcs)
    for (a,b),x in zip(arcs,v):
        e=(p[b],p[a]) if anti else (p[a],p[b])
        out[idx[e]]=x
    return tuple(out)


# All seven binary-cycle-space supports and their Fano incidence.
es4=edges(4)
points=[]
for mask in range(1,64):
    deg=[0]*4
    for k,(a,b) in enumerate(es4):
        if mask>>k&1: deg[a]+=1;deg[b]+=1
    if all(d%2==0 for d in deg): points.append(mask)
lines={tuple(sorted((a,b,a^b))) for a,b in combinations(points,2)}
tri={s for s in points if s.bit_count()==3}; quad=set(points)-tri
require((len(points),len(lines),len(tri),len(quad))==(7,7,4,3),'Fano counts')
require(all(sum(p in l for l in lines)==3 for p in points),'Fano degree')
base=next((a,b,c) for a,b,c in combinations(points,3) if a^b!=c)
coordinates={0:0}
for bits in range(1,8):
    coordinates[bits]=0
    for i in range(3):
        if bits>>i&1: coordinates[bits]^=base[i]
gl=[]
for a,b,c in permutations(points,3):
    if a^b==c: continue
    mp={}
    for bits,src in coordinates.items():
        val=0
        for i,t in enumerate((a,b,c)):
            if bits>>i&1: val^=t
        mp[src]=val
    gl.append(mp)
require(len(gl)==168,'full Fano group')
require(sum({m[p] for p in tri}==tri for m in gl)==24,'tetrahedral subgroup')

census=Counter();patterns=Counter(); strong_records=[]
for mask in range(64):
    arcs=orientation(4,mask); cycs=directed_cycles(4,arcs)
    supp=tuple(sorted(s for s,c in cycs))
    if not supp:
        census['transitive']+=1
        require(not strong(4,arcs),'transitive strongly connected')
    elif len(supp)==1:
        require(supp[0] in tri,'one directed cycle must be triangle')
        out=[sum(a==i for a,b in arcs) for i in range(4)]
        kind='source_C3' if 3 in out else 'sink_C3'
        require((3 in out) != (0 in out),'one-cycle source/sink distinction')
        census[kind]+=1
    else:
        require(strong(4,arcs),'three-cycle type is strong')
        require(len(supp)==3 and sum(s in tri for s in supp)==2,'strong support profile')
        a,b=[s for s in supp if s in tri];c=next(s for s in supp if s in quad)
        require(a^b!=c,'strong cycles must be Fano noncollinear')
        au=permutations_of(arcs,4); an=permutations_of(arcs,4,True)
        stabilizer=[p for p in permutations(range(4)) if tuple(sorted(transport_mask(s,p) for s in supp))==supp]
        require(len(au)==len(an)==1 and len(stabilizer)==2,'orientation/support stabilizers')
        census['strong']+=1;patterns[supp]+=1
        strong_records.append((mask,supp))
require(census==Counter(transitive=24,source_C3=8,sink_C3=8,strong=24),'64-tournament census')
require(len(patterns)==12 and set(patterns.values())=={2},'converse support ambiguity')

arcs=((0,1),(1,2),(2,0),(1,3),(3,0),(2,3))
t4_anti=(1,0,3,2)
require(permutations_of(arcs,4,True)==[t4_anti],'T4 anti map')
raw=Counter(); fixed=Counter(); tested_vectors=0
for v in all_vectors(6,18):
    tested_vectors+=1
    if not balanced(arcs,v,4): continue
    a,b,c=v[2],v[3],v[5]
    require(v==(a+b+c,a+c,a,b,b+c,c),'integer cone coordinates')
    require(transform_flow(arcs,v,t4_anti,True)==(a+b+c,b+c,b,a,a+c,c),'anti swaps a and b')
    raw[sum(v)]+=1
    if transform_flow(arcs,v,t4_anti,True)==v: fixed[sum(v)]+=1
pred=Counter(); pred_fixed=Counter()
for a,b,c in product(range(7),repeat=3):
    w=3*a+3*b+4*c
    if w<=18: pred[w]+=1
    if a==b and w<=18:pred_fixed[w]+=1
require(raw==pred and fixed==pred_fixed,'T4 ordinary/self-converse rational series')
paths=[p for p in permutations(range(4)) if all((p[i],p[i+1])in set(arcs) for i in range(3))]
require(len(paths)==5 and sum(tuple(t4_anti[v] for v in p[::-1])==p for p in paths)==3,'T4 underlying h,b')

# Independently enumerate all decorated SCC words with the small strong alphabet.
N,Q=10,12
# Letter (type,parameter1,parameter2,parameter3), types K1,C3,T4 are 0,1,2.
letters=[((0,0,0,0),1,0)]
letters += [((1,a,0,0),3,3*a) for a in range(Q//3+1)]
letters += [((2,a,b,c),4,3*a+3*b+4*c) for a,b,c in product(range(Q//3+1),repeat=3) if 3*a+3*b+4*c<=Q]

def opposite(letter):
    t,a,b,c=letter
    return (t,b,a,c) if t==2 else letter

def reverse_word(word):
    return tuple(opposite(s) for s in word[::-1])

def underlying_degree(word):
    ts=tuple(t[0] for t in word)
    r=ts.count(2);h=5**r
    if ts!=ts[::-1]: return h
    half=ts[:len(ts)//2];x=5**half.count(2)
    center=ts[len(ts)//2] if len(ts)%2 else -1
    a,b=(5,3) if center==2 else (1,1)
    return x*(x*a-b)//2

words=defaultdict(list)
def visit(word,n,m):
    words[(n,m)].append(word)
    for let,d,w in letters:
        if n+d<=N and m+w<=Q:visit(word+(let,),n+d,m+w)
visit((),0,0)
ordinary=Counter({nm:len(ws) for nm,ws in words.items()})
fixed_words=Counter();under_sc=Counter();merged=Counter();bydegree=Counter()
for nm,ws in words.items():
    reps=set()
    for w in ws:
        op=reverse_word(w)
        if w==op:fixed_words[nm]+=1
        ts=tuple(s[0] for s in w)
        if ts==ts[::-1]:under_sc[nm]+=1
        reps.add(min(w,op))
    merged[nm]=len(reps)
    for w in reps: bydegree[(underlying_degree(w),)+nm]+=1

# Sparse bivariate series operations, with finite truncation only in tests.
def add(*series):
    o=Counter()
    for d in series:o.update(d)
    return +o

def mul(a,b):
    o=Counter()
    for (n,m),x in a.items():
        for (r,s),y in b.items():
            if n+r<=N and m+s<=Q:o[n+r,m+s]+=x*y
    return +o

def renewal(a):
    out=Counter({(0,0):1});term=out
    for i in range(N):
        term=mul(term,a);out=add(out,term)
    return out

atom=[Counter() for i in range(3)];fixatom=[Counter() for i in range(3)]
for let,n,m in letters:
    atom[let[0]][n,m]+=1
    if opposite(let)==let:fixatom[let[0]][n,m]+=1
A=add(*atom)
W=renewal(A)
W2=Counter({(2*n,2*m):v for (n,m),v in W.items() if 2*n<=N and 2*m<=Q})
J=mul(W2,add(Counter({(0,0):1}),*fixatom))
require(W==ordinary,'word rational denominator')
require(J==fixed_words,'self-converse word factorization')
require(all(2*merged[nm]==ordinary[nm]+fixed_words[nm] for nm in ordinary),'converse quotient')
paired_atoms=[]
for p in atom:
    pair=Counter()
    for (n,m),v in p.items():
        for (r,s),w in p.items():
            require(n==r,'fixed atom vertex order')
            if 2*n<=N and m+s<=Q:pair[2*n,m+s]+=v*w
    paired_atoms.append(pair)
U=mul(renewal(add(*paired_atoms)),add(Counter({(0,0):1}),*atom))
require(U==under_sc,'independent paired decoration denominator')
require(under_sc[6,3]==4 and fixed_words[6,3]==0,'underlying vs decorated self-converse hostile')
neutral_hostile=((1,0,0,0),(1,1,0,0))
require(neutral_hostile!=reverse_word(neutral_hostile) and underlying_degree(neutral_hostile)==0,'explicit C3 pair hostile')

# General fixed-underlying-degree formula, in the restricted test alphabet.
# Half-word h=5^r, so weight index r suffices; track two independently constructed carriers.
ordinary_by_r=defaultdict(Counter);paired_by_r=defaultdict(Counter)
ordinary_by_r[0][0,0]=1;paired_by_r[0][0,0]=1
for n in range(N+1):
  for m in range(Q+1):
    for r in range(N//4+1):
      for t,ap in enumerate(atom):
        v=ordinary_by_r[r][n,m]
        if v:
          for (dn,dm),w in ap.items():
            if n+dn<=N and m+dm<=Q:ordinary_by_r[r+(t==2)][n+dn,m+dm]+=v*w
      for t,ap in enumerate(paired_atoms):
        v=paired_by_r[r][n,m]
        if v:
          for (dn,dm),w in ap.items():
            if n+dn<=N and m+dm<=Q:paired_by_r[r+(t==2)][n+dn,m+dm]+=v*w
ordinary_by_r={r:+v for r,v in ordinary_by_r.items()}
paired_by_r={r:+v for r,v in paired_by_r.items()}
U_by_r=defaultdict(Counter);L=defaultdict(Counter);M=defaultdict(Counter)
centers=[(-1,Counter({(0,0):1}),Counter({(0,0):1}))]+[(t,atom[t],fixatom[t]) for t in range(3)]
for r,G in paired_by_r.items():
    x=5**r
    for t,ca,cf in centers:
        h,b=(5,3)if t==2 else(1,1);d=x*(x*h-b)
        prod=mul(G,ca); U_by_r[2*r+(t==2)].update(prod);L[d].update(prod)
for r,F in ordinary_by_r.items():
    x=5**r
    mirrored=Counter({(2*n,2*m):v for(n,m),v in F.items() if 2*n<=N and 2*m<=Q})
    for t,ca,cf in centers:
        h,b=(5,3)if t==2 else(1,1);d=x*(x*h-b)
        M[d].update(mul(mirrored,cf))
pred_degree=Counter()
possible=set(k for k,n,m in bydegree)
for k in possible:
    for n in range(N+1):
      for m in range(Q+1):
        numerator=L[2*k][n,m]+M[2*k][n,m]
        if k%2:
          for r,F in ordinary_by_r.items():
            if 5**r==k:numerator+=F[n,m]-U_by_r[r][n,m]
        require(numerator%2==0,'fixed-degree integral quotient')
        if numerator:pred_degree[k,n,m]=numerator//2
require(pred_degree==bydegree,'corrected fixed-underlying-degree rational construction')

# Generic nontrivial Aut and anti-coset Burnside control on cyclic R5.
r5=tuple((i,(i+d)%5) for i in range(5) for d in (1,2))
au5=permutations_of(r5,5);an5=permutations_of(r5,5,True)
require(len(au5)==len(an5)==5,'R5 symmetry groups')
require(sum(all((p[i],p[i+1]) in set(r5) for i in range(4)) for p in permutations(range(5)))==15,'R5 Hamilton path count')
cyc5=directed_cycles(5,r5)
require(Counter(len(c) for _,c in cyc5)==Counter({3:5,4:5,5:2}),'R5 12 primitive cycle rays')
idx={e:i for i,e in enumerate(r5)}
rays=[tuple(int(e in c) for e in r5) for _,c in cyc5]
bound=14
seen={tuple(0 for e in r5)}; todo=list(seen)
for v in todo:
    for c in rays:
        w=tuple(a+b for a,b in zip(v,c))
        if sum(w)<=bound and w not in seen:seen.add(w);todo.append(w)
brute={v for v in all_vectors(10,8) if balanced(r5,v,5)}
require({v for v in seen if sum(v)<=8}==brute,'R5 cycle decomposition exhaustive control')
for v in seen:
    x=v[0::2];y=v[1::2];ts={x[i]+y[i]+y[(i-1)%5] for i in range(5)}
    require(len(ts)==1,'R5 fractional stable-set cone parametrization')
    t=next(iter(ts));require(sum(v)==5*t-sum(y),'R5 edge-use weight')
orbits={min(transform_flow(r5,v,g) for g in au5) for v in seen}
raw5=Counter(sum(v) for v in seen);orbit5=Counter(sum(v) for v in orbits)
fixed5=Counter()
for v in orbits:
    op=transform_flow(r5,v,an5[0],True)
    if min(transform_flow(r5,op,g)for g in au5)==v:fixed5[sum(v)]+=1
b_orb=Counter();b_fix=Counter()
for v in seen:
    for g in au5:
        if transform_flow(r5,v,g)==v:b_orb[sum(v)]+=1
    for g in an5:
        if transform_flow(r5,v,g,True)==v:b_fix[sum(v)]+=1
require(all(b_orb[i]==5*orbit5[i]and b_fix[i]==5*fixed5[i]for i in range(bound+1)),'ordinary and anti-coset Burnside')
# One directed triangle breaks the order-five automorphism group.
one_triangle=next(c for c in rays if sum(c)==3)
require(sum(transform_flow(r5,one_triangle,g)==one_triangle for g in au5)==1,'decoration changes Aut')

out={
 'status':'PASS',
 'fano':{'point_count':7,'line_count':7,'full_group_order':168,'tetrahedral_subgroup_order':24,'tournament_census':dict(census),'strong_support_patterns':len(patterns),'strong_support_fiber_size':2},
 'T4':{'arc_order':arcs,'anti':t4_anti,'all_edge_vectors_tested':tested_vectors,'flow_weight_bound':18,'ordinary_coefficients':[raw[i]for i in range(19)],'self_converse_coefficients':[fixed[i]for i in range(19)],'underlying_h':5,'underlying_b':3},
 'word_controls':{'vertex_bound':N,'edge_use_bound':Q,'decorated_words':sum(ordinary.values()),'converse_orbits':sum(merged.values()),'actual_self_converse':sum(fixed_words.values()),'all_over_underlying_self_converse':sum(under_sc.values()),'underlying_degree_slices':sorted(possible),'underlying_vs_decorated_hostile':{'vertices':6,'edge_use':3,'ordinary_over_underlying_self_converse':under_sc[6,3],'actual_self_converse':fixed_words[6,3]}},
 'R5':{'flow_weight_bound':bound,'raw_flows':len(seen),'flow_orbits':len(orbits),'simple_cycle_lengths':dict(sorted(Counter(len(c)for _,c in cyc5).items())),'raw_coefficients':[raw5[i]for i in range(bound+1)],'orbit_coefficients':[orbit5[i]for i in range(bound+1)],'self_converse_coefficients':[fixed5[i]for i in range(bound+1)],'one_triangle_stabilizer_order':1,'underlying_Aut_order':5,'underlying_h':3,'same_H_over_decorated_Aut':15},
 'scope':'Finite controls only. All-order claims are proved separately. No Euler tours or new metagraph adjacency are enumerated.'
}
print(json.dumps(out,sort_keys=True,indent=2))
