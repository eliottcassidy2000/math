#!/usr/bin/env python3
"""External audit: vector quaternion arithmetic, signed matrices, graph fibers."""
import importlib.util, itertools, json
spec=importlib.util.spec_from_file_location('target', '/tmp/s167_knots/check.py')
m=importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
def need(ok, why):
    if not ok: raise RuntimeError(why)
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def neg(a): return tuple(-x for x in a)
def conj(a): return (a[0],)+tuple(-x for x in a[1:])
def qmul(a,b):
    w,x,y,z=a; t,u,v,s=b
    return (w*t-x*u-y*v-z*s, w*u+x*t+y*s-z*v,
            w*v-x*s+y*t+z*u, w*s+x*v-y*u+z*t)
def om(a,b):
    p,q=a[:4],a[4:]; r,s=b[:4],b[4:]
    return add(qmul(p,r),neg(qmul(conj(s),q)))+add(qmul(s,p),qmul(q,conj(r)))
E=[tuple(int(i==j) for i in range(8)) for j in range(8)]
units=E+[neg(e) for e in E]; lookup={e:i for i,e in enumerate(units)}
Ls=[tuple(lookup[om(e,u)] for u in units) for e in E]
need(Ls==m.L, 'all left maps match independent vector arithmetic')
for a,b in itertools.product(range(8),repeat=2):
    s,k=m.basis_mul(a,b)
    need(om(E[a],E[b])==tuple(s*x for x in E[k]), 'all 64 vector basis products')
def compose(a,b): return tuple(a[i] for i in b)
I=tuple(range(16)); Z=tuple((i+8)%16 for i in range(16))
# Clifford normal forms from subsets of the first six operators.
forms={}
for mask in range(64):
    p=I
    for bit in range(6):
        if mask>>bit&1: p=compose(p,Ls[bit+1])
    forms[mask,0]=p; forms[mask,1]=compose(Z,p)
G=set(forms.values()); need(len(G)==128,'128 signed Clifford subset forms')
need(all(compose(g,l) in G for g in G for l in Ls),'closed under all 7 maps')
need(len({g[0] for g in G})==16,'transitivity')
# Independent centralizer argument using G itself, not generator propagation.
centralizers=[]
for y in range(16):
    pairs={(g[0],g[y]) for g in G}
    if len(pairs)!=16: continue
    p=tuple(dict(pairs)[x] for x in range(16))
    need(len(set(p))==16,'equivariant self-map bijective')
    need(all(compose(p,g)==compose(g,p) for g in G),'full-group commute')
    centralizers.append(p)
need(set(centralizers)=={I,Z},'centralizer exactly +/-identity')
# Verify tower deck generators by explicit vertex maps and edge sets.
base=list(itertools.combinations(range(5),2)); chords=[e for e in base if e[0]!=0]
v16={e:I for e in base}; v8={e:tuple(range(8)) for e in base}
for a,e in enumerate(chords,1):
    v16[e]=Ls[a]; v8[e]=tuple(x^a for x in range(8))
def edge_set(volts,d):
    return {tuple(sorted((v*d+x,w*d+volts[v,w][x]))) for v,w in base for x in range(d)}
A=edge_set(v16,16); B=edge_set(v8,8)
need(len(A)==160 and len(B)==80,'edge cardinality')
need({tuple(sorted((16*(a//16)+Z[a%16],16*(b//16)+Z[b%16]))) for a,b in A}==A,'sign deck on E')
for shift in range(8):
    need({tuple(sorted((8*(a//8)+((a%8)^shift),8*(b//8)+((b%8)^shift)))) for a,b in B}==B,'translation deck on C')
project=lambda x:8*(x//16)+(x%8)
need({tuple(sorted((project(a),project(b)))) for a,b in A}==B,'global edge projection')
need(all(sum(tuple(sorted((project(a),project(b))))==e for a,b in A)==2 for e in B),'edge double fibers')
print(json.dumps({'status':'PASS','independent_vector_basis_products':64,
 'independent_left_maps':8,'signed_Clifford_subset_forms':128,
 'full_group_centralizer':2,'explicit_sign_deck_verified':True,
 'all_eight_intermediate_decks_verified':True,'edge_double_fibers_verified':True},indent=2))
