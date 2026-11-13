"""Search actual odd-slope seeds from the whole70-family diagonal tensors."""
from context import *
from joint_nine import expo
from collections import defaultdict,Counter
from fractions import Fraction

def run():
    query=defaultdict(list)
    for r in load('even_valuations.json.gz')['records']:
        for d in range(1,13):
            for shift in range(13):
                p=[Fraction(c*(i+1)**shift,2**d) for i,c in enumerate(r['exponential_coefficients'])]
                if any(x.denominator!=1 for x in p):continue
                p=tuple(map(int,p))
                # Sufficient and here exact on the positive clocks: parity of
                # sum of odd-base coefficients is constant.
                if sum(p[::2])%2==0:continue
                query[trim(list(p))].append({'birth_id':r['birth_id'],'port':r['port'],'division':d,'old_shift':shift})
    bank=load('synchronized_descent.json.gz',HERE);hits=[];stats=Counter()
    for r in bank['records']+bank['other_profiles']:
        U={tuple(x[:3]):x[3] for x in r['slope_runs']};E=expo(U)
        for m in range(1,21):
            e=defaultdict(int)
            for (b,c,z),v in E.items():e[z]+=v*(b*c)**m
            for shift in range(13):
                p=trim([e[b]*b**shift for b in range(1,max(e)+1)])
                stats['profiles']+=1
                if p in query:
                    hits.extend({**q,'family':r['family'],'clock':r['clock'],'coupled_length':m,'new_shift':shift,'profile':p} for q in query[p])
    out={'scope':'All78 exact synchronized seed tensors; equal coupled lengths1..20, free-clock shifts0..12; old478 axes with division exponents1..12 and shifts0..12. A match must be a whole odd U profile, so it would pay every exponent for that shifted coefficient family; finite omitted prefixes would still require witnesses.','queries':sum(map(len,query.values())),'distinct_queries':len(query),'statistics':dict(stats),'hits':hits,'matched_axes':sorted({(h['birth_id'],h['port']) for h in hits})}
    save('synchronized_repairs.json.gz',out);print(out if len(str(out))<2500 else {k:v for k,v in out.items() if k!='hits'})
if __name__=='__main__':run()
