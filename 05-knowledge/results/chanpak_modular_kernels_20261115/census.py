#!/usr/bin/env python3
from collections import defaultdict
from pathlib import Path

def table():
    return [tuple(map(int,(n,h,p,c)))+(None if sc=='-' else int(sc),)
        for n,h,p,c,sc in map(str.split,Path(__file__).with_name('profiles.tsv').read_text().splitlines())]

def inspect(rows):
    for n in range(1,10):
        a=[(h,p,c,sc) for i,h,p,c,sc in rows if i==n]
        assert all(p>=(n-1)*h and p+2*h>=2**n and (p+2*h-(-1)**n*h)%3==0 for h,p,c,sc in a)
        s=[(h,(p+2*h).bit_length()-1,sc) for h,p,c,sc in a
            if sc is not None and not ((p+2*h)&(p+2*h-1))]
        print(f'order={n} profiles={len(a)} strong_dyadic={s}')
    profiles={(h,p+2*h) for _,h,p,_,_ in rows}
    for h0,q0 in ((1,2),(3,12),(5,32),(19,256),(217,4096),(223,4096),(599,8192)):
        image={3*q*q0-h*h0 for h,q in profiles}
        windows=defaultdict(set)
        for x in image:
            k=x.bit_length()-1
            if 3*2**k < 2*x:
                windows[k].add(x)
        print(f'fixed h={h0} Q={q0}: distinct={len(image)} leading_11={sum(map(len,windows.values()))}')
        print(f'  band histogram (k,#): {sorted((k,len(v))for k,v in windows.items())}')
        print(f'  first leading_11: {sorted(set.union(*windows.values()))[:16]}')
    print(f'distinct across orders={len(profiles)}')

if __name__=='__main__':inspect(table())
