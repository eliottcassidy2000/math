#!/usr/bin/env python3
"""Exact independent finite audit of the proved saturation increment lemma.
Standard library only; diagnostic checks do not replace the analytic proof.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from random import Random
import hashlib,json


def require(test,message):
    if not test:raise RuntimeError(message)


def close(pred):
    p=list(pred)
    for k in range(len(p)):
        for j in range(len(p)):
            if p[j]>>k&1:p[j]|=p[k]
    require(all(not(p[j]>>j&1) for j in range(len(p))),'Cycle')
    return tuple(p)


def degrees(p):
    return tuple(len(p)-1-p[v].bit_count()-sum(bool(p[u]>>v&1) for u in range(len(p))) for v in range(len(p)))


def count(p):
    full=(1<<len(p))-1
    @lru_cache(None)
    def f(mask):
        if mask==full:return 1
        return sum(f(mask|1<<v) for v in range(len(p)) if not(mask>>v&1) and not(p[v]&~mask))
    return f(0)


def words(p):
    def rec(mask,w):
        if len(w)==len(p):yield w;return
        for v in range(len(p)):
            if not(mask>>v&1) and not(p[v]&~mask):yield from rec(mask|1<<v,w+(v,))
    return list(rec(0,()))


def induced(p,vs,ordered=()):
    at={v:i for i,v in enumerate(vs)}
    q=[sum(1<<at[u] for u in vs if p[v]>>u&1) for v in vs]
    for a,b in zip(ordered,ordered[1:]):q[at[b]]|=1<<at[a]
    return tuple(q)


def main():
    rng=Random(202610190746)
    pool=[(0,4,0,0,8,0,32,1,1),(0,4,0,0,8,24,56,1,1)]
    edges=[(2,z) for z in range(3,7)]+[(2,a) for a in (7,8)]+[(z,a) for z in range(3,7) for a in (7,8)]+[(7,8)]+list(combinations(range(3,7),2))
    while len(pool)<4:
        p=[0,4,0,0,0,0,0,1,1]
        for a,b in edges:
            if rng.random()<.3:p[b]|=1<<a
        p=close(p);d=degrees(p)
        if max(d)<=7 and any(d[v]==7 for v in (2,3,4,5,6,7,8)) and p not in pool:pool.append(p)
    data={}
    for p in pool:
        ww=words(p);gens=Counter();pre={()};suf={()}
        for o in ww:
            h=min(o.index(0),o.index(1));k=max(o.index(0),o.index(1))+1
            for i in range(h+1):
                pre.add(o[:i])
                for j in range(k,10):gens[o[:i],o[j:]]+=1
            for j in range(k,10):suf.add(o[j:])
        data[p]=(ww,gens,pre,suf)
    checked=attempts=zero_u=zero_v=maximum_n=0;digest=hashlib.sha256()
    while checked<80:
        attempts+=1;require(attempts<=10000,'Sampling bound exhausted')
        core=pool[checked%len(pool)];saturated={v for v,d in enumerate(degrees(core)) if d==7}
        na=rng.randrange(4);nb=rng.randrange(4);A=list(range(9,9+na));B=list(range(9+na,9+na+nb));p=list(core)+[0]*(na+nb)
        for side in (A,B):
            for a,b in combinations(side,2):
                if rng.random()<.4:p[b]|=1<<a
        for a in A:
            for v in range(9):
                if v in (0,1) or v in saturated or rng.random()<.8:p[v]|=1<<a
            for b in B:p[b]|=1<<a
        for b in B:
            for v in range(9):
                if v in (0,1) or v in saturated or rng.random()<.8:p[b]|=1<<v
        p=close(p)
        require(tuple(v&511 for v in p[:9])==core,'Core changed')
        if max(degrees(p))>7:continue
        for s in saturated:
            require(all(p[s]>>a&1 for a in A) and all(p[b]>>s&1 for b in B),'Forced direction missing')
        ww,gens,pre,suf=data[core];eA=count(induced(p,A));eB=count(induced(p,B));lc={():Fraction(1)};rc={():Fraction(1)}
        def L(u):
            if u not in lc:lc[u]=Fraction(count(induced(p,A+list(u),u)),eA)
            return lc[u]
        def R(v):
            if v not in rc:rc[v]=Fraction(count(induced(p,list(v)+B,v)),eB)
            return rc[v]
        U={u:L(u)-L(u[:-1]) if u else Fraction(1) for u in pre}
        V={v:R(v)-R(v[1:]) if v else Fraction(1) for v in suf}
        require(min(U.values())>=0 and min(V.values())>=0,'Negative tree increment')
        for u in pre:
            if saturated.intersection(u):require(U[u]==0,'Nonzero saturated prefix increment');zero_u+=1
        for v in suf:
            if saturated.intersection(v):require(V[v]==0,'Nonzero saturated suffix increment');zero_v+=1
        fibre=sum(L(o[:min(o.index(0),o.index(1))])*R(o[max(o.index(0),o.index(1))+1:]) for o in ww)
        full=sum(E*U[u]*V[v] for (u,v),E in gens.items())
        pruned=sum(E*U[u]*V[v] for (u,v),E in gens.items() if not(saturated.intersection(u) or saturated.intersection(v)))
        E=count(p)
        require(fibre==full==pruned==Fraction(E,eA*eB),'Fibre or saturation-pruned expansion mismatch')
        checked+=1;maximum_n=max(maximum_n,len(p));digest.update((json.dumps([list(p),E,str(pruned)],separators=(',',':'))+'\n').encode())
    result={'status':'PASS','scope':'Exact finite audit of a separately proved zero-increment lemma; no c=4 majority exclusion theorem.','actual_seven_thin_posets_checked':checked,'maximum_order':maximum_n,'distinct_core_shapes':len(pool),'saturated_prefix_increments_checked':zero_u,'saturated_suffix_increments_checked':zero_v,'full_and_pruned_fibre_formulas_match_whole_poset_DP':True,'audit_stream_sha256':digest.hexdigest(),'core_shapes':[{'predecessor_masks':list(p),'incomparability_degrees':list(degrees(p))} for p in pool],'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    out=Path(__file__).with_name('saturation_increment_lemma_controls_20261019.json');out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
