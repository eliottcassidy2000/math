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
