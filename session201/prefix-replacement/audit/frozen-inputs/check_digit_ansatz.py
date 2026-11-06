"""Exact controls for the separate-index formal doubling obstruction."""
from itertools import product
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import hashlib,json
from verify import step,need
@lru_cache(None)
def g(u):
    if u==1:return 1
    if u%2==0:return g(u//2)
    h=(u+1)//2;return g(h)+h-h%2

def index_op(u,c):
    if c=='O':return 2*u
    if c=='Z':return 2*u-1
    if c=='S':return u-1 if u%2==0 else u+1
    raise ValueError(c)
def affine(word,r):
    # New index=a*u+b, new g=A+c*u+d, fixed original u mod4.
    a,b,c,d=1,0,F(0),F(0)
    for letter in word:
        m=(a*r+b)%4
        if letter=='O':a,b=2*a,2*b
        elif letter=='Z':c+=a;d+=b-m%2;a,b=2*a,2*b-1
        elif m%2==0:
            c+=F(a,2);d+=F(b,2)-m//2;b-=1
        else:
            c-=F(a,2);d-=F(b+1,2);d+=(m+1)//2%2;b+=1
    return a,b,c,d
central=0
for length in range(1,9):
    words=['J']if length==1 else ['J'+''.join(x)+'J'for x in product('OZJ',repeat=length-2)]
    for word in words:
        s=(1,0,0,0)
        for letter in word:s=step(s,letter)
        need(s[0]%2==s[3]%2==1 and (s[2]-s[1])%2==0,'formal parity invariant')
        need(s[3]%2==1,'odd AB coefficient');central+=1
checks=0
for length in range(4):
    for letters in product('OZS',repeat=length):
        word=''.join(letters)
        for u in range(2,301):
            v=u;legal=True
            for c in word:
                v=index_op(v,c)
                if v<2:legal=False;break
            if not legal:continue
            a,b,c,d=affine(word,u%4)
            need(v==a*u+b and g(v)==g(u)+c*u+d,'physical index affine formula');checks+=1
# Formal empty expansion has coefficients 1 on A,B before and after (1).
# The nonempty coefficient is unchanged because A',B' have coefficients1
# and u',v' contain neither A nor B; the written proof is all-order.
print(json.dumps(dict(status='PASS',central_word_length_bound=8,central_words=central,
    affine_composition_length_bound=3,index_upper_bound=300,exact_affine_controls=checks,
    verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    scope='Coefficient-wise separate-index doubling ansatz only; no functional or occupancy exclusion.'),indent=2,sort_keys=True))
