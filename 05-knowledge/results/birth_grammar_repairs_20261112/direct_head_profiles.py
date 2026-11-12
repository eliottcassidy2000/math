"""Reproducible full-profile search of the complete pointed bank through9."""
from common import *
from hamilton_clocks import trim
from fractions import Fraction
from math import gcd

def run():
    bank=read('hamilton_clocks_9.json.gz')['records'];lookup={tuple(r['exponential_coefficients']):r for r in bank};hits=[]
    for r in read('even_valuations.json.gz')['records']:
        for d in range(9):
            for shift in range(7):
                target=[Fraction(c*(i+1)**shift,2**d) for i,c in enumerate(r['exponential_coefficients'])];target[0]+=1
                if any(c.denominator!=1 for c in target):continue
                target=[int(c) for c in target];G=abs(gcd(*target))
                for factor in sorted({1,G}):
                    p=lookup.get(trim([c//factor for c in target]))
                    if p is not None:hits.append({'birth_id':r['birth_id'],'port':r['port'],'division_exponent':d,'shift':shift,'factor':factor,'profile':p})
    need(len({(h['birth_id'],h['port']) for h in hits})==3 and all(h['factor']==1 for h in hits),'declared complete search controls')
    save('direct_head_profiles.json.gz',{'scope':'All478 even axes; all pointed tournaments through9; parameter shifts0..6; divisions2^d,d0..8; factor1 or full profile content. Whole exponential equality, no finite interpolation. Positive factors all1.','hits':hits});print(json.dumps(hits),flush=True)
if __name__=='__main__':run()
