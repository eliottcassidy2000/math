from itertools import product
from collections import Counter
from pathlib import Path
import json, hashlib
ROOT=Path(__file__).resolve().parent.parent/'session27-orthogonal'

def control(r):
    n=2*r
    zero=(0,)*n
    def Q(v): return sum(v[i]*v[i+r] for i in range(r))%3
    def bil(v,w): return sum(v[i]*w[i+r]+v[i+r]*w[i] for i in range(r))%3
    def canon(v):
        z=next((x for x in v if x),0)
        return tuple(x*z%3 for x in v) if z else zero
    pts=sorted({canon(v) for v in product(range(3),repeat=n) if any(v) and Q(v)==0})
    def extend(U,v):
        return frozenset(tuple((x+a*y)%3 for x,y in zip(u,v)) for u in U for a in range(3))
    levels=[{frozenset([zero])}]
    # Enumerate spaces by actual vector-set closure, not canonical matrices.
    for d in range(r):
        nxt=set()
        for U in levels[-1]:
            for v in pts:
                if v not in U and all(bil(v,u)==0 for u in U): nxt.add(extend(U,v))
        levels.append(nxt)
    maxs=levels[-1]
    assert len(maxs)==2*__import__('math').prod(3**i+1 for i in range(1,r))
    for d, spaces in enumerate(levels):
        den=__import__('math').prod(3**j+1 for j in range(r-d,r))
        for U in spaces: assert sum(U<=L for L in maxs)*den==len(maxs)
    result={'r':r,'isotropic_points':len(pts),'isotropic_spaces_by_dimension':[len(s) for s in levels]}
    if r==2:
        masks=[sum(1<<i for i,v in enumerate(pts) if v in U) for U in maxs]
        # Span closure ranks and pairwise bilinear tests supply the IE weights.
        spans=[frozenset([zero])]*(1<<len(pts)); weights=[0]*len(spans)
        for s in range(len(spans)):
            if s:
                bit=s&-s; i=bit.bit_length()-1; t=s^bit
                if spans[t] is None or any(bil(pts[i],u) for u in spans[t]): spans[s]=None
                else: spans[s]=extend(spans[t],pts[i])
            U=spans[s]
            if U is not None:
                rank={1:0,3:1,9:2}[len(U)]
                weights[s]=(-1)**s.bit_count()*{0:8,1:2,2:1}[rank]
        for i in range(len(pts)):
            for s in range(len(weights)):
                if s&(1<<i): weights[s]+=weights[s^(1<<i)]
        direct=[sum(not (s&m) for m in masks) for s in range(len(weights))]
        assert weights==direct
        result.update(subsets_checked=len(weights),histogram=dict(sorted(Counter(direct).items())))
        expected=json.loads((ROOT/'results.json').read_text())['avoidance_count_histogram']
        assert {str(k):v for k,v in result['histogram'].items()}==expected
    if r==3:
        e=[tuple(int(j==i) for j in range(n)) for i in range(3)]
        ep=tuple((a+b)%3 for a,b in zip(e[0],e[1]))
        dep=e[:2]+[ep]; ind=e
        assert all(bil(v,w)==0 for B in [dep,ind] for v in B for w in B)
        a=sum(not any(v in L for v in dep) for L in maxs)
        b=sum(not any(v in L for v in ind) for L in maxs)
        assert (a,b)==(60,61)
        result.update(dependent_K3_avoiders=a,independent_K3_avoiders=b)
    return result
out={'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.iterdir() if p.is_file()},'controls':[control(r) for r in [1,2,3]],'threshold':7**23+1,'status':'PASS'}
assert out['threshold']==27368747340080916344
Path(__file__).with_name('independent_results.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
