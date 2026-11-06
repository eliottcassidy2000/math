"""Finite regression checks only; all-order arguments are in ROOT_OBSTRUCTION.md."""
import json
from functools import lru_cache
from math import log
@lru_cache(None)
def Q(m):
    return 0 if not m else Q(m//2)+(m//2+(m//2)%2 if m%2 else 0)
def T(m):
    n=m.bit_length();a=m-(1<<(n-1));r=0
    for _ in range(n-1):r=2*r+(a&1);a>>=1
    return (1<<n)-1-r
def require(condition, message):
    if not condition:raise RuntimeError(message)
def J(s):
    h,R,e,d=s;return (h-e+R,h+R,h-e,h)
def O(s):
    h,R,e,d=s;return (h,2*R,0,d)
for m in range(1,1<<17):
    t=T(m)
    require(T(t)==m, f'involution {m}')
    require(t.bit_length()==m.bit_length(), f'dyadic block {m}')
    require(Q(t)-Q(m)==2*(t-m), f'scalar identity {m}')
for n in range(1,18):
    expected=0 if (n-1)%2 else 1<<((n-1)//2)
    require(sum(T(m)==m for m in range(1<<(n-1),1<<n))==expected, f'fixed points {n}')
s=O(J(J((1,2,0,1))))
require(s==(5,12,0,3),'physical seed arithmetic')
require(s[0]+2*s[1]-s[3]*Q(2)==29,'new scalar 29')
require(all(1+2*m-Q(m)!=29 for m in range(28)),'29 is not in root image')
coefficient=log(2)/2+log(4/3)/3
require(coefficient<.5,'mixed-bank strict inequality')
result={'status':'PASS','checks':'Finite regression only; proof is all-order in companion Markdown','involution_max_m':(1<<17)-1,'fixed_point_max_bit_length':17,'root_image_upper_coefficient':log(2)/2,'mixed_bank_seed':s,'mixed_bank_upper_coefficient':coefficient,'added_scalar':29}
print(json.dumps(result,indent=2))
