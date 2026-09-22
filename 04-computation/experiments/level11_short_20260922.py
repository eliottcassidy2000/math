"""Exact eta-product and dyadic recurrence controls. No removable asserts."""
import json
from pathlib import Path
from math import gcd

checks=0
def check(ok):
    global checks
    checks+=1
    if not ok: raise ArithmeticError('exact control failed')

def product(N,bound):
    a=[1]+[0]*(bound-1)
    for step in (1,N):
        for k in range(step,bound,step):
            for unused in range(2):
                for j in range(bound-1,k-1,-1): a[j]-=a[j-k]
    return [0]+a

a=product(11,512)
check(a[1:11]==[1,-2,-1,2,1,2,-2,0,-2,-2])
b=[1,-2]
for r in range(2,101): b.append(-2*b[-1]-2*b[-2])
for r in range(97): check(b[r+4]==-4*b[r])
for r in range(101): check((b[r]==0)==(r%4==3))
for r in range(10): check(a[2**r]==b[r])
for n in range(1,513):
    r=(n&-n).bit_length()-1; m=n>>r
    check(a[n]==b[r]*a[m])
for N in range(8,31): check(product(N,8)[8]==0)
bare=product(1000,12)
check(a[1:12]==bare[1:12] and a[12]==bare[12]-2)
for p in (2,3,5,7,13,17,19,23,29,31):
    count=1+sum((y*y+y-x*x*x+x*x+10*x+20)%p==0 for x in range(p) for y in range(p))
    check(a[p]==p+1-count)
def oddcore(n):
    while n%2==0: n//=2
    return n
cycles=[[1],[5,7],[17,25,37,55,41,61,91]]
for cycle in cycles:
    for i,n in enumerate(cycle):
        check(oddcore(3*n-1)==cycle[(i+1)%len(cycle)])
        check(oddcore(3*(-n)+1)==-cycle[(i+1)%len(cycle)])
for n in range(1,100,2):
    check(oddcore(-3*n-1)==-oddcore(3*n+1))
result={'status':'PASS','checks':checks,'bound':512,'coefficients_1_to_32':a[1:33],
        'dyadic_coefficients_r_0_to_12':b[:13], 'forced_zero_density':'1/15',
        'scope':'Finite checks of product, recurrence, multiplicativity and good-prime point counts; infinite claims proved in companion note.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
print(json.dumps(result))
