#!/usr/bin/env python3
"""Check a finite extension, conditional only on the stated inherited range."""
from verify import verify
w={}
for line in open('frontier_witnesses.txt'):
 h,c,cut=map(int,line.split())
 if h in w: raise ValueError('duplicate witness')
 verify(h,c,cut)
 w[h]=(c,cut)
exceptions={80423:11489,80437:11491,80479:11497}
for n,p in exceptions.items():
 if 7*p!=n or p>11519 or p<2 or any(p%d==0 for d in range(2,1+__import__('math').isqrt(p))):raise ValueError('bad inherited lane')
for n in range(80407,80489,2):
 if n in w or n in exceptions:continue
 split=next(((a,n//a) for a in range(3,int(n**0.5)+1,2) if n%a==0 and a not in (7,21) and n//a not in (7,21) and n//a<=80405),None)
 if not split:raise ValueError(('no certified split',n))
 print(n,*split)
if set(w)!={80407,80429,80447,80449,80471,80473}:raise ValueError('unexpected witnesses')
print('PASS: all odd n in [80407,80487] attained; six strong witnesses; three inherited seven-lane values')
