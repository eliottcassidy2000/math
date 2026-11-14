from functools import cache
def check(x):
    if not x:raise AssertionError('cycle-nine enumerator')
@cache
def g(x):
    if x==0:return 1
    v,b=divmod(x,2)
    return g(v) if b else g(v)+v+v%2
for ell in range(1,16):
    M=1<<(ell-1)
    values=[]
    for x in range(M,2*M):
        p=x%2
        a=84*x*x+75*x+18*g(x)+9*p
        b=84*x*x+87*x+18*g(x)+7*p
        small=[a+11,a+13,b+11,b+13,b+15]
        check(len(set(small))==5 and min(small)>84*x*x and max(small)<84*(x+1)**2)
        values+=small
    check(len(set(values))==5*M)
    check(sum(values)==(1023 if ell==1 else 980*M**3+54*M*M-8*M))
for n in range(10,60):
    count=sum(5*2**(ell-2)*2**(n-6-ell-2) for ell in range(2,n-7))
    check(count==5*(n-9)*2**(n-10))
print('Exact diagonal class and path totals through width 15; independent tape series coefficients pass.')
from itertools import permutations
R=(10,28,25,16,1)
def group(a):
    return [p for p in permutations(range(5)) if all(
        (((R[i]>>j)&1)^((R[p[i]]>>p[j])&1))==a
        for i in range(5) for j in range(i+1,5))]
check(group(0)==[(0,1,2,3,4)] and group(1)==[(0,4,3,2,1)])
p=group(1)[0]
masks=(12,13,15,23,14)
check(tuple(sum(1<<p[i] for i in range(5) if not (u>>i)&1)
            for u in masks)==(19,18,2,4,3))
print('Unique core antiautomorphism takes every selected complemented mask outside the selected bank.')
