from collections import Counter
from fractions import Fraction
def rev(x,n):
    return int(f'{x:0{n}b}'[::-1],2)
for n in range(1,15):
    C=(1<<n)-1
    freq=Counter(Fraction(C-rev(t,n),t) for t in range(1,C))
    assert max(freq.values(),default=0)<=2**((n+1)//2)
    assert sum(k*k for k in freq.values())+2 <= (1<<n)*2**((n+1)//2)+2
assert 2**3<3**2
print('All-ratio fibers through width 14 checked against the exact count bound.')
for n in range(1,10):
    C=(1<<n)-1
    r=[C-rev(t,n) for t in range(C+1)]
    determinants=Counter(a*r[b]-b*r[a] for a in range(C+1) for b in range(C+1))
    assert max(determinants.values(),default=0)<=2*(1<<n)*2**((n+1)//2)
print('All signed determinant fibers through width 9 checked.')
