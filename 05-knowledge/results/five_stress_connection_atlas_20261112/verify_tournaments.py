"""Complete finite round-family and stress-sampling censuses, with controls."""
from pathlib import Path
from itertools import combinations,product
from collections import Counter
from math import comb
import json,sys
from fractions import Fraction
from tournament_tools import *

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'modular_centers_collatz_pell_20261111'))
import tournament_checks as tc

def need(b,msg):
    if not b:raise RuntimeError(msg)

def shift(w):return w[1:]+(1-w[0],)
def orbit(w):
    out=[];v=w
    while v not in out:out.append(v);v=shift(v)
    return out
def representatives(n):
    seen=set();out=[]
    for w in product((0,1),repeat=n):
        if w in seen:continue
        oo=orbit(w);seen.update(oo);out.append(min(oo))
    return out
def matrix(w):return switched_tt(len(w),sum(v<<i for i,v in enumerate(w)))
def changes(w):return sum(a!=b for a,b in zip(w,w[1:]))
def aut_shifts(w):
    v=w;out=[]
    for k in range(2*len(w)):
        if v==w:out.append(k%len(w))
        v=shift(v)
    return out
def move(mask,k,n):return sum(1<<((i-k)%n) for i in range(n) if mask>>i&1)
def composition_orbits(n,k):
    # k here is an odd prime, 3 or 5.
    if n<k:return 0
    return (comb(n-1,k-1)+(k-1)*(n%k==0))//k

def curve_points(w):
    # Every <=3-change word on t=1..n can be produced by one circle
    # polynomial A(t^4+t^2)+b*t+d*t^2+C. Use roots between changes,
    # and a negative fourth root to cancel the cubic coefficient.
    roots=[Fraction(2*i+3,2) for i in range(len(w)-1) if w[i]!=w[i+1]]
    if len(roots)==3:roots.append(-sum(roots))
    # Degrees 0,1,2 already lie in the circle space.
    coeff=[Fraction(1)]
    for r in roots:
        new=[Fraction(0)]*(len(coeff)+1)
        for j,a in enumerate(coeff):new[j]-=r*a;new[j+1]+=a
        coeff=new
    def ev(t):return sum(a*t**j for j,a in enumerate(coeff))
    factor=(1 if w[0]==0 else -1)*(1 if ev(1)>0 else -1)
    coeff=[factor*a for a in coeff]
    need(all((ev(i+1)<0)==bool(v) for i,v in enumerate(w)),'parabola sign witness')
    need(len(coeff)<4 or coeff[3]==0,'no cubic coefficient')
    return coeff

def run(max_n=10):
    levels=[];records=[];literal=0;all_sign_controls=0
    for n in range(1,max_n+1):
        words=representatives(n)
        # Four-dimensional central arrangement, quotiented by global sign.
        parwords=[w for w in product((0,1),repeat=n) if w[0]==0 and changes(w)<=3]
        for w in parwords:curve_points(w);all_sign_controls+=1
        expected=sum(comb(n-1,j) for j in range(min(3,n-1)+1))
        need(len(parwords)==expected,'parabola chamber count')
        parorbits={min(orbit(w)) for w in parwords}
        cubicorbits={min(orbit(w)) for w in product((0,1),repeat=n) if changes(w)<=4}
        need(len(parorbits)==1+composition_orbits(n,3),'parabola type formula')
        need(len(cubicorbits)==1+composition_orbits(n,3)+composition_orbits(n,5),'cubic type formula')
        hist=Counter();shist=Counter();uh=Counter();unit_orbits=0;unit_presentations=0
        if n<=6:
            cs=[canonical(matrix(w)) for w in words]
            need(len(set(cs))==len(words),'word quotient vs graph quotient')
        for w in words:
            rows=matrix(w);need(locally_transitive(rows),'round local transitivity')
            H,U,W,_,_,_=tc.path_tables(rows)
            need(hp(rows)==H,'forward/reverse DP count')
            if n<=7:
                hh,uu,ww=tc.literal_cover_counts(rows)
                need((H,U,W)==(hh,uu,ww),'independent permutation cover count');literal+=1
            hist[H]+=1
            if strong(rows):shist[H]+=1
            auto=aut_shifts(w)
            if n<=7:need(len(auto)==tc.automorphism_count(rows),'shift stabilizer vs automorphisms')
            selected=set()
            for mask in range(1,1<<n):
                verts=[v for v in range(n) if mask>>v&1]
                defect=H-sum(W[a][b] for a,b in combinations(verts,2))
                if defect!=1:continue
                unit_presentations+=1
                canonical_mask=min(move(mask,k,n) for k in auto)
                if canonical_mask in selected:continue
                selected.add(canonical_mask);unit_orbits+=1
                value=sum(U[v] for v in verts);odd=value
                valuation=0
                while odd%2==0:odd//=2;valuation+=1
                uh[odd]+=1
                records.append({'n':n,'word':''.join(map(str,w)),'rows':rows,'H':H,'selected':canonical_mask,'U':value,'odd_part':odd,'v2_U':valuation,'parabola_core':w in parorbits,'cubic_core':w in cubicorbits,'aut_order':len(auto)})
        levels.append({'n':n,'round_isomorphism_types':len(words),'parabola_types':len(parorbits),'cubic_types':len(cubicorbits),'parabola_sign_pairs':len(parwords),'all_switch_sign_pairs':2**(n-1),'H_histogram':dict(sorted(hist.items())),'strong_H_min':min(shist) if shist else None,'unit_mark_presentations_on_unmarked_reps':unit_presentations,'unit_marked_isomorphism_types':unit_orbits,'odd_seed_parts':sorted(uh),'odd_seed_part_count':len(uh)})
        print(json.dumps({k:v for k,v in levels[-1].items() if k not in ('H_histogram','odd_seed_parts')}),flush=True)
    # Full labelled universe through order5 checks that the intrinsic property
    # and sign-word construction agree; no theorem inferred from a sample.
    fullcontrols=[]
    for n in range(1,6):
        pairs=list(combinations(range(n),2));local=set();total=0
        for mask in range(1<<len(pairs)):
            rows=[0]*n
            for bit,(a,b) in enumerate(pairs):
                u,v=(a,b) if mask>>bit&1 else (b,a);rows[u]|=1<<v
            if locally_transitive(rows):local.add(canonical(rows));total+=1
        need(local=={canonical(matrix(w)) for w in representatives(n)},'complete labelled control')
        fullcontrols.append({'n':n,'labelled_total':1<<len(pairs),'labelled_local':total,'local_classes':len(local)})
    targets=[1,3,5,7,9,11,13,19,27,33,1039]
    witnesses={str(c):next((r for r in records if r['odd_part']==c),None) for c in targets}
    allparts={r['odd_part'] for r in records}
    out={'scope':f'Every locally transitive isomorphism class and every nonempty marked set through order {max_n}; marked sets quotiented by actual automorphisms. No arbitrary tournaments beyond controls.','levels':levels,'full_labelled_controls':fullcontrols,'permutation_cover_controls':literal,'exact_parabola_sign_witnesses':all_sign_controls,'seed_records':records,'requested_seed_part_witnesses':witnesses,'first_absent_odd_parts_le_201':[c for c in range(1,202,2) if c not in allparts]}
    dest=HERE/('tournament_receipt.optimized.json' if not __debug__ else 'tournament_receipt.json')
    dest.write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print(json.dumps({'output':dest.name,'marked_orbits':len(records),'requested_covered':[c for c in targets if str(c) in witnesses and witnesses[str(c)] is not None],'missing_le201':out['first_absent_odd_parts_le_201']}),flush=True)

if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 10)
