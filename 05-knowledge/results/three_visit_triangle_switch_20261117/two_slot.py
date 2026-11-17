"""Exact exponential coefficients by fourier-like inclusion/exclusion on visits."""
from collections import Counter
from math import comb
from replay import rows, subset_count

def append(a, ports):
    n = len(a)
    return [v | (0 if i in ports else 1<<n) for i,v in enumerate(a)]+[sum(1<<i for i in ports)]

def coefficients(first, pulse):
    # spine a,M,b,c; then apices x,y; then V; optionally final N
    a = append(rows((first,1,3,0)), set() if pulse else {0,4})
    if pulse: a = append(a,{0,4})
    M,V = 1,6
    U = set(range(len(a))) - {M,V}
    def traces(prev, remain, r, s):
        if not remain and r and s:
            yield r,s
        for v in sorted(remain|{M,V}):
            if prev is None or a[prev]>>v&1:
                yield from traces(v,remain-{v},r+(v==M),s+(v==V))
    C = Counter()
    for r,s in traces(None,U,0,0):
        for i in range(1,r+1):
            for j in range(1,s+1):
                C[i,j] += (-1)**(r+s-i-j)*comb(r,i)*comb(s,j)
    return {k:v for k,v in C.items() if v}

if __name__ == '__main__':
    for pulse, expect in [
        (False,((71,-80,22),(-116,187,-76),(42,-104,66))),
        (True,((-5,-15),(3,29),(4,28))),
    ]:
        for first in (2,3):
            C = coefficients(first,pulse)
            known = {(i+1,j+1):c+2*(i==j==0 and first==2)
                     for i,row in enumerate(expect) for j,c in enumerate(row)}
            if C != known: raise AssertionError((C,known))
            for m in range(1,5):
                for k in range(1,5):
                    a = rows((first,)+(1,)*m+(3,0))
                    for _ in range(k): a = append(a,set() if pulse else {0,m+3})
                    if pulse: a = append(a,{0,m+3})
                    exp = sum(c*i**m*j**k for (i,j),c in C.items())
                    actual = subset_count(a)
                    if exp != actual: raise AssertionError((first,m,k,exp,actual))
        print('one pulse' if pulse else 'uniform return', expect)
