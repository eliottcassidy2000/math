from math import isqrt
from pathlib import Path
root=Path(__file__).resolve().parent.parent
def prime(x): return x>=2 and all(x%d for d in range(2,isqrt(x)+1))
for label,lo,hi,base in [('frontier',84261,90001,84259),('100001',90003,100001,90001)]:
    rows=[tuple(map(int,l.split())) for l in (root/f'certificates_{label}'/'witnesses.txt').read_text().splitlines()]
    hs=[h for h,_,_ in rows]; s=set(hs)
    if len(hs)!=len(s) or any(not (lo<=h<=hi and h%2) for h in hs): raise ValueError('bad witness targets')
    splits={}
    for n in range(lo,hi+1,2):
        if n in s: continue
        pair=next(((d,n//d) for d in range(3,isqrt(n)+1,2) if n%d==0 and d not in (7,21) and n//d not in (7,21) and n//d<=base),None)
        if pair is None: raise ValueError(('uncovered',n))
        splits[n]=pair
    required={n for n in range(lo,hi+1,2) if prime(n) or n%7==0 and prime(n//7)}
    if s!=required: raise ValueError(('lane mismatch',sorted(s-required),sorted(required-s)))
    print(f'PASS {label}: witnesses={len(s)}, products={len(splits)}, endpoint={hi}, endpoint_product={splits.get(hi)}')
