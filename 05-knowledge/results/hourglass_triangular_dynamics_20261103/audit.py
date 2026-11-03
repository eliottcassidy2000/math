"""Separate exact verification: does not import verify.py."""
from collections import Counter
from math import gcd, isqrt, prod
import json


def check(b, msg):
    if not b:raise RuntimeError(msg)


def triangular(n):
    return n>=0 and isqrt(8*n+1)**2==8*n+1


def expected_count(k):
    return 4*(triangular(k)+triangular(k-2))+2*(k in [1,2,5,6])+8*(k==4)


def orbit_test(k):
    # An independent deliberately loose interval. Outside it, magnitude escapes.
    B=4*abs(k)+20
    state={}
    for start in range(-B,B+1):
        path=[];at={};b=start
        while -B<=b<=B and b not in state and b not in at:
            at[b]=len(path);path.append(b);b=b*(b+1)//2-k
        good=(b in at) or state.get(b,False)
        for x in path:state[x]=good
    return sorted(b for b,good in state.items() if good)


def main():
    tested=0
    for k in range(-20,1501):
        S=orbit_test(k)
        check(len(S)==expected_count(k),('orbit count',k,S))
        tested+=1
    # At fixed denominator 16 all other odd numerators immediately lose denominator 4.
    den_checks=0
    from fractions import Fraction
    for c in range(-99,101,2):
        if c%8==5:continue
        for a in range(-41,42,2):
            check((Fraction(a,4)**2-Fraction(c,16)).denominator!=4,('denominator',c,a))
            den_checks+=1
    # All modulus-3 line words on the six-vertex layer graph, with a fixed root.
    from itertools import product
    E=[(0,1),(0,2),(1,3),(1,4),(1,5),(2,3),(2,4),(2,5)]
    grad={tuple((x[v]-x[u])%3 for u,v in E)
          for tail in product(range(3),repeat=5) for x in [(0,)+tail]}
    words=0
    for w in product(range(3),repeat=8):
        rectangles=(w[6]-w[5]-w[3]+w[2],w[7]-w[5]-w[4]+w[2])
        hourglass=w[1]+w[5]-w[0]-w[2]
        zero=all(t%3==0 for t in rectangles+(hourglass,))
        check(zero==(w in grad),('mod3 lines',w));words+=1
    # Verify all-size wheel formulas independently from maximal zero-free blocks.
    wheel_count=0
    for m in list(range(4,1002,2))+[1160,1890,2310,5510,30030]:
        flags=[gcd(2*i+1,m)==1 for i in range(m//2)]
        if all(flags):continue
        start=flags.index(False)
        blocks=[];length=0
        for t in range(1,len(flags)+1):
            if flags[(start+t)%len(flags)]:length+=1
            elif length:blocks.append(length);length=0
        check(length==0,'cyclic break')
        observed=Counter(blocks)
        P=[p for p in range(3,m+1,2) if m%p==0 and all(p%d for d in range(2,isqrt(p)+1))]
        p0=min(P);L=m//(2*prod(P))
        s=lambda r:L*prod(p-r for p in P) if r<p0 else 0
        pred={r:s(r)-2*s(r+1)+s(r+2) for r in range(1,p0)}
        check(observed==Counter({r:v for r,v in pred.items() if v}),('blocks',m))
        wheel_count+=1
    # Exact transfer matrix, independent column convention.
    import sympy as sp
    arbitrary_layers=0
    for m in range(2,5):
        for sizes in product(range(1,4),repeat=m):
            es=[((i,u),(i+1,v)) for i in range(m-1)
                for u in range(sizes[i]) for v in range(sizes[i+1])]
            tree={((i,0),(i+1,0)) for i in range(m-1)}
            tree.update(((0,u),(1,0)) for u in range(1,sizes[0]))
            tree.update(((i-1,0),(i,u)) for i in range(1,m) for u in range(1,sizes[i]))
            non=[e for e in es if e not in tree]
            rows=[]
            for i in range(1,m-1):
                for u in range(1,sizes[i]):
                    d={((i-1,0),(i,u)):1,((i,u),(i+1,0)):1,
                       ((i-1,0),(i,0)):-1,((i,0),(i+1,0)):-1}
                    rows.append([d.get(e,0) for e in non])
            for i in range(m-1):
                for u in range(1,sizes[i]):
                    for v in range(1,sizes[i+1]):
                        d={((i,u),(i+1,v)):1,((i,u),(i+1,0)):-1,
                           ((i,0),(i+1,v)):-1,((i,0),(i+1,0)):1}
                        rows.append([d.get(e,0) for e in non])
            check(len(rows)==len(non),('square minor',sizes))
            check(not rows or abs(sp.Matrix(rows).det())==1,('unimodularity',sizes))
            arbitrary_layers+=1
    F=[6,3,1,0,0,1,3,6];A=sp.zeros(8)
    for j,k in enumerate(F):A[k,j]=1
    z=sp.symbols('z')
    check(sp.expand((sp.eye(8)-z*A).det())==1-z**3,'transfer determinant')
    check([int((A**r).trace()) for r in range(1,13)]==[0,0,3]*4,'return clock')
    print(json.dumps({'independent_integer_orbit_parameters':tested,
       'parameters':[-20,1500],'denominator_tests':den_checks,
       'ternary_line_words':words,'arbitrary_layer_sizes':arbitrary_layers,'wheel_moduli':wheel_count,
       'transfer_determinant':str(1-z**3),'return_traces':[0,0,3]*4},indent=2,sort_keys=True))


if __name__=='__main__':main()
