#!/usr/bin/env python3
"""Independent audit of mixed (one lower, two upper) core completeness
and its exceptional weighted comparison certificate. Standard library only.
"""
from itertools import combinations, permutations
from functools import lru_cache
from fractions import Fraction as F
from random import Random
from pathlib import Path
import json

PRED=(0,4,0,4,13,5)
PAIRS=((0,3,3),(3,1,4),(1,5,1),(5,4,1))

def close(p):
    p=list(p)
    for k in range(len(p)):
        for j in range(len(p)):
            if p[j]>>k&1:p[j]|=p[k]
    assert all(not(p[j]>>j&1) for j in range(len(p)))
    return tuple(p)

def words(p):
    n=len(p)
    def rec(mask,o):
        if len(o)==n:yield o;return
        for i in range(n):
            if not(mask>>i&1) and not(p[i]&~mask):yield from rec(mask|1<<i,o+(i,))
    return list(rec(0,()))

def count(p):
    n=len(p)
    @lru_cache(None)
    def f(mask):
        if mask==(1<<n)-1:return 1
        return sum(f(mask|1<<i) for i in range(n) if not(mask>>i&1) and not(p[i]&~mask))
    return f(0)

def induced(p,vs,ordered=()):
    at={v:i for i,v in enumerate(vs)}
    q=[sum(1<<at[u] for u in vs if p[v]>>u&1) for v in vs]
    for a,b in zip(ordered,ordered[1:]):q[at[b]]|=1<<at[a]
    return q

def completeness_counts():
    results=[]
    for c in range(3):
        n=c+5;Z=list(range(3,3+c));U=[c+3,c+4]
        edges=[(2,z) for z in Z]+[(2,a) for a in U]+[(z,a) for z in Z for a in U]+[tuple(U)]+list(combinations(Z,2))
        seen=set()
        for mask in range(1<<len(edges)):
            p=[0,4,0]+[0]*c+[1,1]
            for j,(a,b) in enumerate(edges):
                if mask>>j&1:p[b]|=1<<a
            seen.add(close(p))
        word_count=admissible=0
        for p in seen:
            inc=lambda t:{u for u in range(n) if u!=t and not(p[t]>>u&1 or p[u]>>t&1)}
            assert inc(0)=={1,2,*Z};assert inc(1)=={0,*Z,*U}
            pairs=[(a,b) for a,b in combinations(range(n),2) if b in inc(a)]
            ww=words(p);word_count+=len(ww)
            for M in ww:
                at={v:i for i,v in enumerate(M)};span=at[1]-at[0]
                if span>0 and all(abs(at[a]-at[b])<=span for a,b in pairs):admissible+=1
        results.append({'common_count':c,'distinct_cores_with_natural_common_and_upper_labels':len(seen),'legal_core_words':word_count,'admissible_maxspan_orders':admissible})
    assert [(r['distinct_cores_with_natural_common_and_upper_labels'],r['legal_core_words'],r['admissible_maxspan_orders']) for r in results]==[(7,92,60),(40,2123,722),(357,79681,15899)]
    return results

def exception_certificate():
    ww=words(PRED)
    # Independent direct permutation enumeration agrees with recursive generation.
    independent=[]
    for o in permutations(range(6)):
        at={v:i for i,v in enumerate(o)}
        if all(not(PRED[v]>>u&1) or at[u]<at[v] for u in range(6) for v in range(6)):independent.append(o)
    assert set(ww)==set(independent)
    prefixes={()};suffixes={()}
    for o in ww:
        left=min(o.index(0),o.index(1));right=max(o.index(0),o.index(1))+1
        prefixes.update(o[:j] for j in range(1,left+1))
        suffixes.update(o[j:] for j in range(right,6))
    prefixes=sorted(prefixes,key=lambda u:(len(u),u));suffixes=sorted(suffixes,key=lambda v:(len(v),v))
    table=[]
    for u in prefixes:
        row=[]
        for v in suffixes:
            chosen=[o for o in ww if (not u or o[:len(u)]==u) and(not v or o[-len(v):]==v)]
            coefficient=sum(6-sum(w*(o.index(a)<o.index(b)) for a,b,w in PAIRS) for o in chosen)
            assert coefficient>=0;row.append(coefficient)
        table.append(row)
    expected=[[0,1,6,7,9,1,6,3,3],[10,4,8,4,9,3,4,2,2],[7,1,3,0,2,0,0,0,0]]
    assert table==expected
    counts=[sum(o.index(a)<o.index(b) for o in ww) for a,b,_ in PAIRS]
    assert len(ww)==37 and counts==[27,24,22,23]
    M=(2,0,3,1,5,4);at={v:i for i,v in enumerate(M)}
    for a,b,_ in PAIRS:
        assert not(PRED[a]>>b&1 or PRED[b]>>a&1) and at[a]<at[b]
    return prefixes,suffixes,table,{'extensions':37,'comparison_counts':counts,'weighted_numerator':222,'prefixes':[list(u) for u in prefixes],'suffixes':[list(v) for v in suffixes],'coefficient_table_for_6E_minus_weighted_N':table}

def exterior_audit(prefixes,suffixes,table):
    rng=Random(202610190703);maximum_order=0;mingap=None
    for trial in range(200):
        na=rng.randrange(4);nb=rng.randrange(4);A=list(range(6,6+na));B=list(range(6+na,6+na+nb));n=6+na+nb
        p=list(PRED)+[0]*(na+nb);maximum_order=max(maximum_order,n)
        def edge(a,b):p[b]|=1<<a
        for side in (A,B):
            for a,b in combinations(side,2):
                if rng.random()<.4:edge(a,b)
        for a in A:
            for v in (0,1,4,5):edge(a,v)
            for v in (2,3):
                if rng.random()<.5:edge(a,v)
            for b in B:edge(a,b)
        for b in B:
            for v in (0,1,2):edge(v,b)
            for v in (3,4,5):
                if rng.random()<.5:edge(v,b)
        p=close(p);E=count(p);weighted=0
        for a,b,w in PAIRS:
            q=list(p);q[b]|=1<<a;weighted+=w*count(q)
        delta=6*E-weighted;assert delta>=0
        eA=count(induced(p,A));eB=count(induced(p,B))
        lc={():F(1)};rc={():F(1)}
        def L(u):
            if u not in lc:lc[u]=F(count(induced(p,A+list(u),u)),eA)
            return lc[u]
        def R(v):
            if v not in rc:rc[v]=F(count(induced(p,list(v)+B,v)),eB)
            return rc[v]
        evaluated=0
        for i,u in enumerate(prefixes):
            U=L(u)-L(u[:-1]) if u else 1;assert U>=0
            for j,v in enumerate(suffixes):
                V=R(v)-R(v[1:]) if v else 1;assert V>=0
                evaluated+=table[i][j]*U*V
        assert evaluated==F(delta,eA*eB)
        gap=F(delta,E);mingap=gap if mingap is None else min(mingap,gap)
    return {'actual_posets_checked':200,'maximum_order':maximum_order,'exact_tree_evaluation_matches':True,'minimum_normalized_deficit':str(mingap)}

def main():
    pre,suf,tab,certificate=exception_certificate()
    result={'status':'PASS','core_completeness_counts':completeness_counts(),'exceptional_weighted_certificate':certificate,'actual_exterior_audit':exterior_audit(pre,suf,tab)}
    out=Path(__file__).with_name('mixed12_exception_independent_audit_20261019.json');out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
