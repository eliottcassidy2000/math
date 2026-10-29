"""Exact reflection trace of good-time components and finite residue sets."""
from fractions import Fraction as F
from itertools import combinations
from functools import reduce
from math import gcd
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
CHECKS=0
def need(b,msg):
    global CHECKS
    CHECKS+=1
    if not b:raise RuntimeError(msg)

def good_components(speeds,delta):
    intervals=[(F(0),F(1))]
    for v in speeds:
        allowed=[((F(k)+delta)/v,(F(k+1)-delta)/v) for k in range(v)]
        new=[]
        for a,b in intervals:
            for c,d in allowed:
                lo,hi=max(a,c),min(b,d)
                if lo<=hi:new.append((lo,hi))
        merged=[]
        for a,b in sorted(new):
            if merged and a<=merged[-1][1]:merged[-1]=(merged[-1][0],max(b,merged[-1][1]))
            else:merged.append((a,b))
        intervals=merged
    return intervals

def good_residues(speeds,delta,q):
    return [k for k in range(q) if all(delta.denominator*min(v*k%q,(-v*k)%q)>=delta.numerator*q for v in speeds)]

def main():
    cases=0;grids=0;largest=0
    for size in range(1,5):
        for speeds in combinations(range(1,11),size):
            if reduce(gcd,speeds)!=1:continue
            for delta in sorted({F(1,size+1),F(1,4),F(2,5),F(1,2)}):
                components=good_components(speeds,delta);fixed=sum(a+b==1 for a,b in components)
                need(set(components)=={(1-b,1-a) for a,b in components},'literal reflected components')
                need(fixed==int(all(v%2 for v in speeds)),'exact reflection trace')
                need(len(components)%2==fixed,'component parity')
                for a,b in components:
                    need(0<a<=b<1,'proper components exclude origin')
                for q in range(1,41):
                    good=good_residues(speeds,delta,q)
                    need(set(good)=={(-k)%q for k in good},'residue reflection')
                    need(len(good)%2==int(q%2==0 and all(v%2 for v in speeds)),'residue parity')
                    grids+=1
                cases+=1;largest=max(largest,len(components))
    # A composite LRC boundary and a positive-length mixed-parity witness.
    examples=[]
    for speeds,delta in [((1,2),F(1,3)),((1,2),F(1,4)),((1,2),F(2,5)),
                         ((1,3,5),F(1,4)),(tuple(range(1,14)),F(1,14))]:
        pieces=good_components(speeds,delta)
        examples.append({'speeds':speeds,'delta':str(delta),'components':[[str(a),str(b)] for a,b in pieces],
                         'component_count':len(pieces),'fixed_components':sum(a+b==1 for a,b in pieces)})
    need(len(good_components(tuple(range(1,14)),F(1,14)))==6,'AP13 exact isolated boundary')
    need(all(a==b for a,b in good_components(tuple(range(1,14)),F(1,14))),'AP13 no interval interior')
    # The unique reflection preserving every individual distance function.
    reflection_cases=0
    for speeds in ((1,2),(2,3),(3,5),(2,5,7)):
        for den in range(1,21):
            for num in range(den):
                phase=F(num,den)
                necessary=all((v*phase).denominator==1 for v in speeds)
                need(necessary==(phase==0),'primitive phase-preserving reflection centre')
                reflection_cases+=1
    report={'status':'FINITE-EXACT','primitive_speed_max':10,'speed_cardinalities':[1,2,3,4],
            'thresholds':'unique of 1/(k+1), 1/4, 2/5, 1/2','component_cases':cases,
            'residue_moduli':[1,40],'residue_cases':grids,'max_components_seen':largest,
            'reflection_phase_cases':reflection_cases,'examples':examples,'checks':CHECKS}
    (HERE/'lrc_reflection.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report))
if __name__=='__main__':main()
