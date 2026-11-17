from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import comb, factorial, gcd

def oddpart(k, lo=1):
    if not k:
        yield ()
    for l in range(lo,k+1,2):
        for rest in oddpart(k-l,l):
            yield (l,)+rest

def sc(n):
    k=n//2
    ans=Fraction()
    for part in oddpart(k):
        lengths=part*2+(1,)*(n%2)
        e=sum((x-1)//2 for x in lengths)+sum(gcd(x,y) for x,y in combinations(lengths,2))
        d=(e+len(part))//2
        if 2*d!=e+len(part):raise ArithmeticError(part)
        term=Fraction(2**d,2**k)
        for l,a in Counter(part).items():
            term*=Fraction(2**(l-1),l)**a/factorial(a)
        ans+=term
    if ans.denominator!=1:raise ArithmeticError(ans)
    return ans.numerator

known=(1,1,1,2,2,8,12,88,176,2752)
if tuple(sc(n) for n in range(len(known)))!=known:
    raise ArithmeticError('small class totals')
for n in range(0,17):
    k=n//2;E=comb(n,2)
    base=Fraction(2**((E+k)//2),2**k*factorial(k))
    main=1+Fraction(8*comb(k,3),2**max(2*n-6,0))+Fraction(384*comb(k,5),2**max(4*n-20,0))+Fraction(640*comb(k,6),2**max(4*n-24,0))
    rem=Fraction(sc(n))-base*main
    if rem<0 or (n<14 and rem):raise ArithmeticError((n,rem))
    print(n,sc(n),rem)
