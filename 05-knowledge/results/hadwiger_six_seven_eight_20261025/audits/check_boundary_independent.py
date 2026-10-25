"""Independent proof controls for cube_boundary; producer code not imported."""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import hashlib
import json
import sys


def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)


points=tuple(product((0,1),repeat=3))
signs=tuple((-1)**sum(x) for x in points)
pairs=tuple(combinations(range(3),2))


def marginal(w):
    return tuple(sum(w[i] for i,x in enumerate(points) if all(x[k]==v for k,v in zip(pair,vals)))
                 for pair in pairs for vals in product((0,1),repeat=2))


def walsh(w):
    return tuple(sum(w[i]*(-1)**sum(a*b for a,b in zip(x,s)) for i,x in enumerate(points))
                 for s in points)


hist=Counter()
fibers=defaultdict(list)
full_records=[]
for w in product((0,1),repeat=8):
    sums=marginal(w)
    if all(sums):
        hist[sum(w)]+=1
        full_records.append(w)
    fibers[sums].append(w)
    missing=[i for i in range(8) if not w[i]]
    independent=all(sum(a!=b for a,b in zip(points[i],points[j]))!=1 for i,j in combinations(missing,2))
    require(all(sums)==independent,"support projection versus cube independent set")
require(hist=={4:2,5:8,6:16,7:8,8:1},"35 support relations")
ambiguous=[v for v in fibers.values() if len(v)>1]
parities={tuple(int(sum(x)%2==eta) for x in points) for eta in (0,1)}
require(len(ambiguous)==1 and set(ambiguous[0])==parities,"only Boolean counting collision")

weighted=defaultdict(list)
for w in product(range(3),repeat=8):
    weighted[marginal(w)].append(w)
weighted_hist=Counter()
for group in weighted.values():
    base=group[0]
    weighted_hist[len(group)]+=1
    for w in group:
        lam=w[0]-base[0]
        require(all(w[i]-base[i]==lam*signs[i] for i in range(8)),"weighted marginal kernel")
    allowable=[]
    for lam in range(-2,3):
        test=tuple(base[i]+lam*signs[i] for i in range(8))
        if all(0<=v<=2 for v in test):
            allowable.append(test)
    require(set(group)==set(allowable),"complete bounded integer fiber")

# Direct summation over points and coefficient convolution for all Boolean pairs.
all_tables=tuple(product((0,1),repeat=8))
all_walsh=[walsh(w) for w in all_tables]
for i,a in enumerate(all_tables):
    for j,b in enumerate(all_tables):
        w=tuple(x*y for x,y in zip(a,b))
        actual=walsh(w)
        for s in range(8):
            convolution=sum(all_walsh[i][t]*all_walsh[j][s^t] for t in range(8))
            require(convolution==8*actual[s],"complete gluing convolution")

root=Path(__file__).resolve().parents[1]
source=root/'cube_boundary'/'receipt.json'
packet=json.loads(source.read_text(encoding='utf-8'))
geometry_results={}
for name,data in packet['unit_graphs'].items():
    coords=[tuple(Fraction(v) for v in xy) for xy in data['points']]
    require(len(set(coords))==len(coords),"injective coordinates")
    def distance_squared(i,j):
        a=coords[i][0]-coords[j][0]
        b=coords[i][1]-coords[j][1]
        return a*a+a*b+b*b
    edges={p for p in combinations(range(len(coords)),2) if distance_squared(*p)==1}
    declared={tuple(e) for e in data['edges']}
    require(edges==declared,"independent all-pairs unit graph "+name)
    fixed={int(v):c for v,c in data['fixed_colors'].items()}
    ports=tuple(data['ports'])
    neighbors={i:{v if u==i else u for u,v in edges if i in (u,v)} for i in range(len(coords))}
    require(all(len(neighbors[v])==1 for v in fixed),"all prescribed added vertices are leaves")
    weights=[]
    for bits in points:
        pinned=dict(fixed)
        pinned.update(dict(zip(ports,bits)))
        remaining=tuple(i for i in range(len(coords)) if i not in pinned)
        # Brute-force all 6^4 assignments for F, rather than producer's menu recursion.
        count=0
        for values in product(range(6),repeat=len(remaining)):
            colors=dict(pinned)
            colors.update(dict(zip(remaining,values)))
            if all(colors[u]!=colors[v] for u,v in edges):
                count+=1
        weights.append(count)
    want=tuple(0 if sum(x) in (0,3) else 4 for x in points) if name=='F' else (3,)*8
    require(tuple(weights)==want,"independent six-palette extension count")
    require(len(coords)==(32 if name=='F' else 19) and len(edges)==(31 if name=='F' else 15),"geometry size")
    W=walsh(weights)
    require(sum(weights)==24 and all(W[s]==0 for s in (1,2,4,7)),"shared coarse Fourier data")
    require(all(marginal(weights)),"full existential pair shadows")
    geometry_results[name]={"vertices":len(coords),"edges":len(edges),"weights":weights,"walsh":W,
                             "prescribed_leaves":len(fixed),"coordinate_sha256":hashlib.sha256(json.dumps(data['points']).encode()).hexdigest()}
require(geometry_results['F']['walsh']!=geometry_results['G']['walsh'],"pair modes retain difference")

# Matrix is independently reconstructed from tuple Hamming distance.
M=[[-2 if i==j else 2 if sum(a!=b for a,b in zip(x,y))==1 else 0 for j,y in enumerate(points)] for i,x in enumerate(points)]
for s in points:
    char=[(-1)**sum(a*b for a,b in zip(x,s)) for x in points]
    eig=4*(1-sum(s))
    require(all(sum(M[i][j]*char[j] for j in range(8))==eig*char[i] for i in range(8)),"K3 cube eigenvalue")

data={"status":"FINITE-EXACT independent verification",
      "source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      "producer_proof_sha256":hashlib.sha256((root/'cube_boundary'/'RESULT.md').read_bytes()).hexdigest(),
      "producer_check_sha256":hashlib.sha256((root/'cube_boundary'/'check.py').read_bytes()).hexdigest(),
      "producer_receipt_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
      "support_relations":256,"full_pair_shadow_histogram":dict(sorted(hist.items())),
      "weighted_tables":6561,"weighted_fiber_histogram":dict(sorted(weighted_hist.items())),
      "gluing_pairs":65536,"unit_graphs":geometry_results}
encoded=json.dumps(data,indent=2,sort_keys=True)+'\n'
if len(sys.argv)>1:
    Path(sys.argv[1]).write_text(encoded,encoding='utf-8',newline='\n')
else:
    print(encoded,end='')
