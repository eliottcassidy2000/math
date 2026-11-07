#!/usr/bin/env python3
"""Compare general inverse-phi recursion with the claimed supplier formula."""
import json, math
from pathlib import Path

def prime(n):
    return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))
def divisors(n):
    ds=[]
    for d in range(1,math.isqrt(n)+1):
        if n%d==0:
            ds.append(d)
            if d*d!=n: ds.append(n//d)
    return sorted(ds)
def inverse(m):
    ps=[d+1 for d in divisors(m) if d+1>2 and prime(d+1)]
    out=set()
    def rec(i,rem,n):
        if i==len(ps):
            if rem==1: out.update((n,2*n))
            elif rem&(rem-1)==0: out.add(n*(2*rem))
            return
        p=ps[i];rec(i+1,rem,n)
        f=p-1;t=p
        while rem%f==0:
            rec(i+1,rem//f,n*t)
            f*=p;t*=p
    rec(0,m,1)
    return sorted(out)
def dyadic_coeff(a,omit=None):
    if a<0:return 0
    weights=[1<<i for i in range(4) if (1<<i)<=a and (1<<(1<<i))+1!=omit]
    # Checks use a<=12; F4 weight16 is immaterial.
    sums=[0]
    for w in weights:sums += [s+w for s in sums]
    return sum(2 if s==a else 1 for s in sums if s<=a)
def formula(a,q):
    c=sum(dyadic_coeff(a-j) for j in range(1,a+1) if prime((1<<j)*q+1))
    w=q-1
    if w&(w-1)==0:
        j=w.bit_length()-1
        c+=dyadic_coeff(a-j,omit=q)
    return c
rows=[]
for q in [3,5,7,11,13,17,19,23,29,31]:
    for a in range(13):
        m=(1<<a)*q
        actual=inverse(m);expect=formula(a,q)
        assert len(actual)==expect,(m,actual,expect)
        rows.append({'q':q,'a':a,'m':m,'count':len(actual)})
receipt={'status':'PASS','universe':'q in [3,5,7,11,13,17,19,23,29,31], 0<=a<=12','targets':len(rows),'method':'general candidate-prime support recursion, compared with supplier-convolution formula','rows':rows}
Path(__file__).with_name('supplier_independent_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print('PASS:',len(rows),'complete inverse fibers match supplier convolution')
