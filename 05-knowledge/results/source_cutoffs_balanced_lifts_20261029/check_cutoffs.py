"""Finite signed-residual LP duality, with rational post-verification."""
from fractions import Fraction as F
from itertools import product
import json
import numpy as np
from scipy.optimize import linprog
import sympy as sp


def require(ok,msg):
    if not ok:raise RuntimeError(msg)


def truth(T,z):
    n=len(T)
    M=sp.eye(n)
    for i,j in enumerate(T):
        if j>=0:M[i,j]-=sp.Rational(z.numerator,z.denominator)
    result=M.inv()*sp.ones(n,1)*sp.Rational((1-z).numerator,(1-z).denominator)
    return [F(v) for v in result]


def box_vertices(T,z,lo,hi,indices):
    n=len(T);r=len(indices);rows=[];rhs=[]
    for i in indices:
        row=[F(0)]*n;row[i]=1
        if T[i]>=0:row[T[i]]-=z
        rows.append(row);rhs.append(1-z)
    best=None
    # Vertices may have redundant cube bounds; trying all n-r chosen bounds suffices.
    from itertools import combinations
    for chosen in combinations(range(n),n-r):
        for bits in product((0,1),repeat=n-r):
            A=[row[:] for row in rows];b=rhs[:]
            for i,bit in zip(chosen,bits):
                row=[F(0)]*n;row[i]=1;A.append(row);b.append((lo,hi)[bit][i])
            mat=sp.Matrix(A)
            if mat.det()==0:continue
            sol=[F(v) for v in mat.inv()*sp.Matrix(b)]
            if any(sol[i]<lo[i] or sol[i]>hi[i] for i in range(n)):continue
            yield sol


def score(w,T,z,v,lo,hi):
    e=w[:]
    for i,t in enumerate(v):
        e[i]-=t
        if T[i]>=0:e[T[i]]+=z*t
    return (1-z)*sum(v)+sum(t*(hi[i] if t>=0 else lo[i]) for i,t in enumerate(e)),e


def signed_proposal(w,T,z,lo,hi,Fset):
    n=len(T);m=len(Fset)
    # e=w-(I-zQ)v = e+ - e-. Both nonnegative parts are explicit LP variables.
    A=np.zeros((n,m+2*n))
    for c,i in enumerate(Fset):
        A[i,c]=1
        if T[i]>=0:A[T[i],c]-=float(z)
    A[:,m:m+n]=np.eye(n);A[:,m+n:]=-np.eye(n)
    cost=np.r_[np.full(m,float(1-z)),list(map(float,hi)),-np.array(list(map(float,lo)))]
    result=linprog(cost,A_eq=A,b_eq=list(map(float,w)),
                   bounds=[(None,None)]*m+[(0,None)]*(2*n),method='highs')
    require(result.success,'finite signed LP proposal')
    v=[F(0)]*n
    for c,i in enumerate(Fset):v[i]=F(float(result.x[c])).limit_denominator(10**8)
    value,e=score(w,T,z,v,lo,hi)
    return value,v,e


def finite_controls():
    z=F(1,2);w=[F(1,6),F(1,12),F(1,20)];count=0;negative=0;coarse=0;nonsuper=0
    for T in product(range(-1,3),repeat=3):
        actual=truth(T,z)
        for Fset in [(0,),(0,1),(0,1,2)]:
            for mode in range(3):
                lo=[F(0)]*3 if mode==0 else [x/2 for x in actual]
                # Asymmetric authentic information tests genuine non-supersolution boxes.
                hi=[F(1)]*3 if mode<2 else [actual[0],F(1),F(1)]
                if any(hi[i] < 1-z+(z*hi[T[i]] if T[i]>=0 else 0) for i in Fset):nonsuper+=1
                optimum=max(sum(x*y for x,y in zip(w,a)) for a in box_vertices(T,z,lo,hi,Fset))
                val,v,e=signed_proposal(w,T,z,lo,hi,Fset)
                require(val==optimum,'exact rational signed primal / independent box-vertex dual')
                require(val>=sum(x*y for x,y in zip(w,actual)),'bound dominates genuine source value')
                if any(x<0 for x in v):negative+=1
                # Coarse supersolutions have an optimal nonnegative stopped occupation flow.
                if mode<2:
                    M=sp.eye(len(Fset));index={v:i for i,v in enumerate(Fset)}
                    for c,i in enumerate(Fset):
                        if T[i] in index:M[index[T[i]],c]-=sp.Rational(1,2)
                    sol=M.inv()*sp.Matrix([w[i] for i in Fset]);positive=[F(0)]*3
                    for i,x in zip(Fset,sol):positive[i]=F(x)
                    val2,e2=score(w,T,z,positive,lo,hi)
                    require(val2==optimum and all(x>=0 for x in positive+e2),'nonnegative finite-section optimum')
                    coarse+=1
                count+=1
    # Actual Collatz 3->5->1; F={3}, old upper information is stored at3.
    w=[F(1,6),F(1,12)];T=(1,-1);lo=[F(0),F(0)];hi=[F(3,4),F(1)]
    val,v,e=signed_proposal(w,T,z,lo,hi,(0,))
    # Unmodelled original active sources have exactly total mass1/4 and bound1.
    require(val+F(1,4)==F(5,12) and v==[-F(1,6),F(0)],'native reverse-credit receipt')
    return {'rational_duality_cases':count,'coarse_nonnegative_optima':coarse,
            'boxes_failing_supersolution_on_support':nonsuper,
            'proposals_using_negative_flow':negative,
            'native_3_to_5':{'z':'1/2','original_weights':['1/6','1/12'],
               'old_bound':'11/24','new_bound':'5/12','flow':list(map(str,v)),
               'residual':list(map(str,e)),'outside_original_mass':'1/4'}}


if __name__=='__main__':
    print(json.dumps(finite_controls(),indent=2,sort_keys=True))
