#!/usr/bin/env python3
"""Independent quotient-residue checks, including common exponent extensions."""
from pathlib import Path
import runpy,math
m=runpy.run_path(str(Path(__file__).with_name('check.py')))
def require(condition):
    if not condition: raise AssertionError('Independent audit failed')
def rem(p,a,modulus):
    numerator=pow(p,a+1,modulus*(p-1))-1
    require(numerator%(p-1)==0)
    return (numerator//(p-1))%modulus
count=common_count=valuation_count=0
for n in range(1,6):
    for mask in range(1<<(n*(n-1)//2)):
        ps,exponents=m['labels'](n,mask)
        base=math.lcm(*(a+1 for a in exponents))
        # This global multiplier safely imposes H=2 for every target,
        # also when p == 1 (mod r). Extra valuations are permitted.
        high=base*math.prod(r*r for r in ps)
        for i,p in enumerate(ps):
            require(exponents[i]>0 and exponents[i]%2==0)
            for j,r in enumerate(ps):
                if i==j: continue
                outgoing=pow(p,(r-1)//2,r)==1
                require((rem(p,exponents[i],r)==0)==outgoing)
                count+=1
                for desired in (0,2):
                    length=base if (base-1)%4==desired else 3*base
                    require((length-1)%4==desired)
                    require((rem(p,length-1,r)==0)==outgoing)
                    common_count+=1
                    length=high if (high-1)%4==desired else 3*high
                    require((length-1)%4==desired)
                    require((rem(p,length-1,r*r)==0)==outgoing)
                    valuation_count+=1
require((count,common_count,valuation_count)==(21300,42600,42600))
print('PASS: per-source',count,'common/mod4',common_count,'H=2/mod4',valuation_count)
