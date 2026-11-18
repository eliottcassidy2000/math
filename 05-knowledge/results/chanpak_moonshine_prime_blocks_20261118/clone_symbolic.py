"""Extract a literal smallest-total clone witness using exact rational series.

Requires SymPy. Run this file alongside clone_replay.py and its imports.
"""
from functools import lru_cache
from itertools import combinations_with_replacement
from math import prod, factorial
from random import Random
from gap_replay import rand_t
import sympy as sp
from clone_replay import inflate, read, TT
from frame_replay import cut, add, WITNESSES

def rational(A, marks, q, z):
    """Include words missing some colors; the caller imposes full support."""
    s=len(marks)
    X=sp.diag(*z)
    adjacency=sp.Matrix([[(A[i]>>j)&1 for j in marks] for i in marks])
    M=sp.eye(s)-adjacency*X
    det=sp.expand(M.det())
    trans=X*M.adjugate()
    outer=tuple(i for i in range(len(A)) if i not in marks)
    e=len(outer)

    @lru_cache(None)
    def stretch(a,b):
        l=sp.Matrix([[1 if a is None else (A[a]>>i)&1 for i in marks]])
        r=sp.Matrix([1 if b is None else (A[i]>>b)&1 for i in marks])
        empty=a is None or b is None or (A[a]>>b)&1
        return sp.expand((l*trans*r)[0]+int(empty)*det)

    @lru_cache(None)
    def go(todo, a, slots):
        # End current slot, or append a single as-yet-unused ordinary vertex.
        if not slots:
            return sp.S.One if not todo else sp.S.Zero
        result=stretch(a,None)*go(todo,None,slots-1)
        for b in todo:
            result+=stretch(a,b)*go(tuple(i for i in todo if i!=b),b,slots)
        return sp.expand(result)
    return go(outer,None,q),det,e+q

def derivative(P,D,g,z):
    """Numerator for mixed derivative, with D^(g+len(z)) denominator."""
    for h,x in enumerate(z,g):
        P=sp.expand(D*sp.diff(P,x)-h*P*sp.diff(D,x))
    return P

def smallest_witness(A,B,marks,q):
    """Return (widths, difference B minus A), or None for universal equality."""
    z=sp.symbols('t:'+str(len(marks)))
    P,D,g=rational(A,marks,q,z)
    Q,E,g2=rational(B,marks,q,z)
    assert g==g2
    h=g+len(z)
    ab,bb=sp.Poly(D,z),sp.Poly(E,z)
    common=sp.lcm(ab,bb)
    # Normalize so the common denominator has constant one.
    common=common.mul_ground(1/common.eval(dict.fromkeys(z,0)))
    a=derivative(P,D,g,z)*common.exquo(ab).as_expr()**h
    b=derivative(Q,E,g,z)*common.exquo(bb).as_expr()**h
    num=sp.Poly(b-a,z)
    if num.is_zero:
        return None
    u,c=min(num.terms(),key=lambda term:(sum(term[0]),term[0]))
    widths=tuple(i+1 for i in u)
    # A surviving raw run coefficient is c/prod(widths); literal
    # width-r TT modules fill r nonempty runs in r! ways.
    assert c%prod(widths)==0
    assert sum(widths)<=2*len(z)*(g+len(z))
    return widths, int(c)*prod(factorial(i-1) for i in widths)

def run():
    A=(10,4,9,2);B=cut(A,8);marks=(0,1,2)
    assert B==(2,12,1,5)
    R=(2,1,0,3)
    assert all(((A[j]>>i)&1)==((B[R[i]]>>R[j])&1)
               for i in range(4) for j in range(4))
    for q in (1,2,3):
        k,w=smallest_witness(A,B,marks,q)
        Q={i:TT(j) for i,j in zip(marks,k)}
        assert read(inflate(B,Q),q)-read(inflate(A,Q),q)==w
        for r in range(1,sum(k)):
            for cuts in combinations_with_replacement(range(r+1),len(k)-1):
                Ks=tuple(b-a for a,b in zip((0,)+cuts,cuts+(r,)))
                if min(Ks)<=0:continue
                T={i:TT(j) for i,j in zip(marks,Ks)}
                assert read(inflate(B,T),q)==read(inflate(A,T),q)
        print(f'q={q}: smallest total width {sum(k)}, example {k}, difference {w}')
    assert smallest_witness(A,A,marks,2) is None
    Q={0:TT(2),1:TT(1),2:TT(1)}
    assert tuple(read(inflate(T,Q),2) for T in (A,B))==(84,100)

    U,I,_,_,_=WITNESSES[0]; A=add(U,I); B=cut(A,256)
    m=(0,3,8); z=sp.symbols('x y z')
    assert read(A,2)==read(B,2)==38344
    assert read(inflate(A,{0:TT(2)}),2)==225088
    assert read(inflate(B,{0:TT(2)}),2)==222280
    for T, D in ((A,1-prod(z)),(B,1)):
        C=sp.Matrix([[(T[i]>>j)&1 for j in m] for i in m])
        assert (sp.eye(3)-C*sp.diag(*z)).det()==D
    print('equal-height elementary pair splits its marked pole and fails at width four')

    rng=Random(117122)
    for trial in range(18):
        n=rng.randrange(3,6)
        A=rand_t(n,rng); B=cut(A,rng.randrange(1<<n))
        marks=tuple(rng.sample(range(n),rng.randrange(1,n+1)))
        q=1+trial%2
        w=smallest_witness(A,B,marks,q)
        if w:
            r, delta=w
            Q={i:TT(j) for i,j in zip(marks,r)}
            assert read(inflate(B,Q),q)-read(inflate(A,Q),q)==delta
        else:
            for j in range(3):
                Q={i:rand_t(rng.randrange(1,4),rng) for i in marks}
                assert read(inflate(B,Q),q)==read(inflate(A,Q),q)
    print('18 additional symbolic versus literal substitution checks')

if __name__=='__main__':run()
