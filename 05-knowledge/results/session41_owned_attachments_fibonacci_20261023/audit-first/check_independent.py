#!/usr/bin/env python3
"""Independent small certificates; checks remain active with python -O."""
import itertools,json
from pathlib import Path

def need(p,m):
    if not p: raise ValueError(m)

def subsets(weights):
    d={}
    for bits in itertools.product((0,1),repeat=len(weights)):
        if any(a*b for a,b in zip(bits,bits[1:])):continue
        # zero slot is kept as a forbidden-adjacency spacer, never selected
        if any(b and w==0 for b,w in zip(bits,weights)):continue
        s=tuple(i for i,b in enumerate(bits) if b)
        d.setdefault(sum(weights[i] for i in s),[]).append(s)
    return d
f=[1,2]
while f[-1]<100:f.append(sum(f[-2:]))
a=subsets([1,1]+f); b=subsets([1,0,1]+f)
need(len(a[1])==len(b[1])==3,'three units')
need(len(a[2])==2 and len(b[2])==3,'zero changes adjacency')
def r(n):
    if n<0:return 0
    for w in reversed(f[1:]):
        if n>=w:n-=w
    return n
for n in range(1,101):
    need(len(a[n])==3-r(n-1),f'compressed fibre {n}')
    expected=3 if n==1 else 4-r(n-1)-r(n-2)
    need(len(b[n])==expected,f'zero fibre {n}')
# Two closed affine walks, each admitting an owner, but no common owner.
reflections=[tuple((-x+c)%3 for x in range(3)) for c in (0,1)]
fixed=[{x for x in range(3) if p[x]==x} for p in reflections]
need(all(fixed) and not set.intersection(*fixed),'joint owner hostile')
# Every edge of a triangle admits a gauge, but its three translations cannot
# be simultaneously flattened. Composition has nonzero translation.
edge_shifts=(0,0,1)
gauges=[g for g in itertools.product(range(3),repeat=3)
        if all((g[(i+1)%3]-g[i])%3==edge_shifts[i] for i in range(3))]
need(not gauges,'joint gauge hostile')
# All 64 affine F3 maps triples: common fixed points are exact section test.
affine=[tuple((e*x+c)%3 for x in range(3)) for e in (1,2) for c in range(3)]
counts={0:0,1:0,3:0}
for p,q in itertools.product(affine,repeat=2):
    owners=[x for x in range(3) if p[x]==x and q[x]==x]
    counts[len(owners)]+=1
need(counts=={0:26,1:9,3:1},'affine common fixed census')
# Parity in the inserted sequence really is O,E,O,O,E,O,O,...;
# calling the three consecutive leading ones globally Fibonacci is false.
seq=[1,0,1,1,2,3,5,8,13,21]
need(all(seq[i+2]==seq[i+1]+seq[i] for i in range(len(seq)-2)), 'extended recurrence')
need([x%2 for x in seq]==[1,0,1,1,0,1,1,0,1,1],'parity')
need(1+1!=1,'three consecutive units break recurrence')
result={'status':'PASS','independent_subset_bound':100,'first_adjacency_difference':{'N':2,'compressed_count':len(a[2]),'zero_spacer_count':len(b[2])},'reflection_fixed_sets':[sorted(x) for x in fixed],'triangle_shift':list(edge_shifts),'common_fixed_pair_census':counts,'slot_counts_1_to_20':{'compressed':[len(a[n]) for n in range(1,21)],'zero_spacer':[len(b[n]) for n in range(1,21)]}}
# Additional independent shared-frame controls.
I=(1,0,0,1); S=(0,1,1,0)
pts=[x+y for x,y in itertools.product((I,S),repeat=2) if x[0]==y[0]]
need(len(pts)==2,'shared table exact fibre')
x=pts[0]
def mat_kernel(g):
    return (g[0]+g[1],g[2]+g[3],g[0]+g[2],g[1]+g[3],g[4]+g[5],g[6]+g[7],g[4]+g[6],g[5]+g[7],g[0]-g[4])
kernel=[g for g in itertools.product((-1,0,1),repeat=8) if any(g) and not any(mat_kernel(g))]
need(len(kernel)==2,'two coupled atoms')
need(all(sum(abs(t) for t in g)==8 for g in kernel),'coupled move length')
# Commutation within root rotations, while owner sections remain empty.
rot=[tuple((v+k)%3 for v in range(3)) for k in range(3)]
t=rot[1];s=(0,2,1)
def comp(p,q):return tuple(p[q[x]] for x in range(3))
central=[g for g in rot if comp(g,t)==comp(t,g)]
central_full=[g for g in central if comp(g,s)==comp(s,g)]
need(len(central)==3 and len(central_full)==1,'centralizer vs section')
need(not [v for v in range(3) if t[v]==v],'empty section')
result.update({'coupled_tables':'PASS','kernel_atoms':kernel,'rotational_centralizers':[len(central),len(central_full)],'both_section_counts':0})
print(json.dumps(result,indent=2))
