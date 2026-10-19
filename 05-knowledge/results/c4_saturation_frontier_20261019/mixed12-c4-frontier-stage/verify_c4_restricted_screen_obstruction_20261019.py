#!/usr/bin/env python3
"""Exact self-contained check of a limited-screen obstruction and its repair.
The rational mixture is a generator-cone point, not an asserted actual exterior.
"""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations
from pathlib import Path
import hashlib,json

P=(0,4,0,0,4,0,60,29,45)
M=(2,0,3,4,5,6,1,8,7)
PAIRS=((0,4),(4,1),(1,8),(8,7))
SUPPORT=(((),(8,7),374647),((3,),(8,),21897),((3,5),(),631448),((3,5),(6,),1153771))
DEN=2181763


def require(test,message):
    if not test:raise RuntimeError(message)


def words():
    def rec(mask,o):
        if len(o)==9:yield o;return
        for v in range(9):
            if not(mask>>v&1) and not(P[v]&~mask):yield from rec(mask|1<<v,o+(v,))
    return list(rec(0,()))


@lru_cache(None)
def extension_count(mask=0):
    if mask==511:return 1
    return sum(extension_count(mask|1<<v) for v in range(9) if not(mask>>v&1) and not(P[v]&~mask))


def main():
    ww=words();require(len(ww)==len(set(ww))==extension_count()==3184,'Core enumeration mismatch')
    inc=[(a,b) for a,b in combinations(range(9),2) if not(P[b]>>a&1 or P[a]>>b&1)]
    degree=[sum(v in pair for pair in inc) for v in range(9)]
    require(degree==[6,7,3,5,5,6,4,4,4],'Degree mismatch')
    require(M in set(ww),'Candidate order not a core extension')
    r={v:i for i,v in enumerate(M)};require(all(abs(r[a]-r[b])<=r[1]-r[0] for a,b in inc),'Root not maximum span')
    require(all(tuple(sorted(pair)) in inc and r[pair[0]]<r[pair[1]] for pair in PAIRS),'Comparison not incomparable/forward')
    gens=Counter()
    for o in ww:
        h=min(o.index(0),o.index(1));k=max(o.index(0),o.index(1))+1
        for i in range(h+1):
            for j in range(k,10):gens[o[:i],o[j:]]+=1
    require(len(gens)==1204,'Generator count mismatch')
    mixture=[];probs=[Q(0)]*4
    for u,v,numerator in SUPPORT:
        require((u,v) in gens and 1 not in u+v,'Invalid or saturation-deleted generator')
        chosen=[o for o in ww if (not u or o[:len(u)]==u) and(not v or o[-len(v):]==v)]
        E=len(chosen);N=[sum(o.index(a)<o.index(b) for o in chosen) for a,b in PAIRS]
        weight=Q(numerator,DEN)
        for i in range(4):probs[i]+=weight*Q(N[i],E)
        mixture.append({'prefix':list(u),'suffix':list(v),'E':E,'N':N,'convex_weight':str(weight),'generator_weight':str(weight/E)})
    require(sum(Q(num,DEN) for _,_,num in SUPPORT)==1,'Mixture not normalized')
    require(probs==[Q(1458653,2181763)]*4 and min(probs)>Q(2,3),'Restricted-screen obstruction failed')
    # Check the repair against every generator, not merely the four mixture terms.
    bounds=[];stream=hashlib.sha256()
    for (u,v),E in sorted(gens.items()):
        chosen=[o for o in ww if (not u or o[:len(u)]==u) and(not v or o[-len(v):]==v)]
        require(len(chosen)==E,'Generator count disagreement')
        N=sum((o.index(0)<o.index(3))+(o.index(3)<o.index(1)) for o in chosen)
        require(158*N<=203*E,'Repair inequality failed')
        bounds.append(Q(N,E));stream.update((json.dumps([u,v,E,N],separators=(',',':'))+'\n').encode())
    require(max(bounds)==Q(203,158)<Q(4,3),'Repair cone bound mismatch')
    result={'status':'PASS','scope':'Only the old four-comparison screen fails in a generator-cone relaxation. The full screen is repaired by a different common neighbor. No actual exterior or conjecture counterexample is claimed.','predecessor_masks':list(P),'candidate_order':list(M),'core_degrees':degree,'core_extensions':len(ww),'generator_columns':len(gens),'four_comparisons':[list(p) for p in PAIRS],'positive_mixture':mixture,'four_probabilities':[str(p) for p in probs],'strict_majority_margin_3p_minus_2':str(3*probs[0]-2),'repair_comparisons':[[0,3],[3,1]],'universal_repair_upper_bound':'203/158','repair_column_audit_sha256':stream.hexdigest(),'verifier_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    Path(__file__).with_name('c4_restricted_screen_obstruction_controls_20261019.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
