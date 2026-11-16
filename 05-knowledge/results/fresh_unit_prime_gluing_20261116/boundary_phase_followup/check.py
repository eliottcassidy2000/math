from functools import cache
from itertools import product

seeds = [
    ([228,29,304,501,801,770,822,374,3,399],(6015,5934,6283,729),1,12217),
    ([242,156,465,821,640,406,186,8,723,231],(6017,5931,6286,729),1,12217),
    ([314,32,387,710,526,412,439,19,666,231],(6017,7883,8219,1185),1,16102),
    ([344,849,579,918,548,463,24,855,596,97],(6019,3661,4390,1185),2,16102),
]

def profile(a, r=9):
    n=len(a)
    F=(1<<n)-1
    for i in range(n):
        assert not (a[i]>>i)&1 and 0<=a[i]<=F
        for j in range(i):
            assert ((a[i]>>j)&1)+((a[j]>>i)&1)==1
    def bits(s):
        while s:
            b=s&-s
            yield b.bit_length()-1
            s^=b
    @cache
    def end(s,j):
        rem=s^(1<<j)
        return sum(end(rem,i) for i in bits(rem) if a[i]>>j&1) if rem else 1
    def h(s):
        return sum(end(s,j) for j in bits(s)) if s else 1
    E=O=0
    def rooted(s,j):
        nonlocal E,O
        w=h(F^s)
        if s.bit_count()%2:
            O+=w
        else:
            E+=w
        for i in bits(a[j]&(F^s)):
            rooted(s|(1<<i),i)
    rooted(1<<r,r)
    assert O-E==end(F,r)
    return h(F),E,O,h(F^(1<<r))

for a,p,s,Q in seeds:
    assert profile(a)==p
    for _ in range(s):
        a=[x|(1<<len(a)) for x in a]+[0]
    assert profile(a)==(p[0],Q,Q,p[3])
    print('padded',p[0],Q,Q,p[3])

def test(u,v,d):
    w=(abs(d)&-abs(d)).bit_length()-1
    cap=max(u,v+1)*(1<<w)
    hs=[(a,b) for a,b in product(range(10),repeat=2) if (u<<a)-(v<<b)==d]
    assert len(hs)<=2
    assert all((u<<a)<=cap if a==w or a==b else (u<<a)<= (v<<w)+abs(d) for a,b in hs)

for u,v in product(range(1,20,2),repeat=2):
    for d in range(-60,61):
        if d: test(u,v,d)
print('dyadic crossing controls passed')

from math import factorial
@cache
def max_product(g):
    ans=1
    for n in range(3,10):
        q=next(i for i in range(1,101) if i**3>=5**(n-1))
        if q<=g:
            ans=max(ans,factorial(n+1)*max_product(g//q))
    return ans
assert max_product(100)==factorial(10)
B=(1<<10)*max_product(100)
assert 756*2**30 <= (B+1)*231 < 756*2**31
assert 131+2*B < 756*2**31
print('quantitative cutoff controls passed')
