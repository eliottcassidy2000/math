#!/usr/bin/env python3
"""Small exact controls for the algebraic-boundary note; standard library only."""
from fractions import Fraction
from itertools import permutations
import json


def mul(A, B):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]


def power(A, n):
    out=[[int(i==j) for j in range(len(A))] for i in range(len(A))]
    while n:
        if n&1: out=mul(out,A)
        A=mul(A,A);n//=2
    return out


def addpoly(a,b):
    c=[0]*max(len(a),len(b))
    for i,x in enumerate(a):c[i]+=x
    for i,x in enumerate(b):c[i]+=x
    return c


def mulpoly(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return c


def charpoly(A):
    n=len(A);out=[0]*(n+1)
    for p in permutations(range(n)):
        sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
        term=[sign]
        for i,j in enumerate(p):term=mulpoly(term,[-A[i][j],int(i==j)])
        out=addpoly(out,term)
    return out


def ideals_by_rank(p):
    rows=[{0}]
    for _ in p:
        rows.append({I|1<<x for I in rows[-1] for x in range(len(p))
                     if not I>>x&1 and p[x]&I==p[x]})
    return rows


def normalized_key(I,k,n):
    return (tuple(x-k for x in range(k) if not I>>x&1),
            tuple(x-k for x in range(k,n) if I>>x&1))


def sign(x):return (x>0)-(x<0)


def run():
    # The natural distance-two incomparability poset, at an interior rank.
    n=16;p=[(1<<max(0,x-2))-1 for x in range(n)];rows=ideals_by_rank(p);k=6
    left=sorted(rows[k],key=lambda I:normalized_key(I,k,n))
    right=sorted(rows[k+1],key=lambda I:normalized_key(I,k+1,n))
    keys=[normalized_key(I,k,n) for I in left]
    assert keys==[normalized_key(J,k+1,n) for J in right]
    T=[[int(I&J==I) for J in right] for I in left]
    assert T==[[1,0,1,1],[1,0,0,0],[1,1,0,0],[0,1,1,0]]
    cp=charpoly(T);assert cp==[-1,-3,-1,-1,1]
    T3=power(T,3);T4=power(T,4)
    assert any(x==0 for row in T3 for x in row) and all(x>0 for row in T4 for x in row)
    # Mod-2 quartic has no linear factor; the only irreducible quadratic leaves x+1.
    assert sum(cp)%2==1 and cp[0]%2==1
    rem=[x%2 for x in cp]
    for deg in range(4,1,-1):
        if rem[deg]:
            for j in range(3):rem[deg-2+j]^=1
    assert rem[:2]==[1,1]
    cr=max(Fraction(T4[i][j]*T4[h][l],T4[i][l]*T4[h][j])
           for i in range(4) for h in range(4) for j in range(4) for l in range(4))
    assert cr==2

    # Exact two-sided stationary-equality hostile; no poset event realization claimed.
    A=[[2,1],[1,1]];K=[[0,1],[-1,0]];N=[[6,3],[3,3]];E=[[2,2],[0,1]]
    assert all(0<=E[i][j]<=N[i][j] for i in range(2) for j in range(2))
    skew=[]
    for m in range(8):
        for r in range(8):
            numerator=mul(mul(power(A,m),K),power(A,r))[0][0]
            assert sign(numerator)==sign(r-m)
            if m<3 and r<3:skew.append([m,r,numerator])

    # Positive rational transfer with complex subdominant modes.
    P=[[0,0,1],[1,0,0],[0,1,0]]
    M=[[1+P[i][j] for j in range(3)] for i in range(3)]
    f=[[3,3,3]];c=[[2,1,0]];b=[[1],[0],[0]];osc=[]
    for t in range(21):
        h=mul(power(M,t),b)
        pt=Fraction(mul(c,h)[0][0],mul(f,h)[0][0])
        expected=Fraction(1,3)+Fraction((1,0,-1)[t%3],3*4**t)
        assert pt==expected
        if t<9:osc.append(str(pt))
    assert charpoly(M)==[-4, -3, -3, 1]

    # Conservative D=7, epsilon=1/100 buffered margin certificate.
    assert Fraction(17,18)**105 < Fraction(1,400)
    print(json.dumps({'status':'all exact fixed controls passed',
        'quartic':{'n':n,'rank':k,'states':keys,'T':T,'charpoly_low_to_high':cp,
                   'T4':T4,'T4_max_crossratio':str(cr)},
        'two_sided_skew_controls':skew,'complex_mode_matrix':M,
        'complex_mode_probabilities_t0_to_t8':osc,
        'degree7_buffer_blocks':106,'degree7_buffer_ranks_each_side':742},indent=2))


if __name__=='__main__':run()
