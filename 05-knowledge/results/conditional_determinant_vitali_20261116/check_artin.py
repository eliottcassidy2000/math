from math import gcd
def orbit(m):
    start=1;a=start;z=[]
    while True:
        z.append(a);a=(2*a)%m
        if a==start:return z
def union(m):
    z=orbit(m)
    a={(34*t-1)%m for t in z}
    b={(32*t+5)%m for t in z}
    return a,b
for m in range(3,258,2):
    p=next(d for d in range(2,m+1) if m%d==0)
    k=m
    while k%p==0:k//=p
    claim=k==1 and p not in (3,17) and len(orbit(m))==m-m//p
    a,b=union(m)
    assert (len(a|b)==m)==claim
assert len(orbit(67))==66 and pow(2,66,67**2)==671
assert len(orbit(69))==22 and len(set.union(*union(69)))==32
print('All odd moduli through 257 and the 67-adic lifting witness passed.')
for n in range(1,13):
    lo=32*2**n+3;hi=lo+2
    term,p=(hi,3) if n%2 else (lo,5) if n%4==0 else {
        2:(hi,7),6:(lo,7),10:(hi,13)
    }[n%12]
    assert term>p and term%p==0 and pow(2,12,p)==1
print('Period-12 twin-prime exclusion verified.')
from math import isqrt
matrix=((67,81,24),(32,86,42),(0,22,22))
def slots(n):
    return 1,2**n-2,3**n-3*2**n+3
tt=sum(a*c*b for a,line in zip(slots(8),matrix) for c,b in zip(line,slots(3)))
assert tt==1734041 and all(x%d for x in (tt,tt+2) for d in range(2,isqrt(x)+1))
print('Physical two-axis twin-prime escape verified.')
