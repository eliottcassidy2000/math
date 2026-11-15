#!/usr/bin/env python3
"""Exact suffix-cone replay. Needs NumPy for object arrays and SymPy for one rational inverse."""
from collections import Counter
from pathlib import Path
from fractions import Fraction
import json
import sys

import numpy as np
from sympy import Matrix

NAMES = '0.2 1.2 0.123 2.13 0.1.3 0.3.1 0.13.2 0.1.123 0.3.12 0.12.3 1.2.3 1.12.3 0.1.23'.split()
PREF = (1,0,2,3,4,7,9)
SUFF = (0,1,2,3,5,6,7)


def fragments(state, a):
    # State lists ordered (head,tail) ages of the current directed paths.
    # -1 is old; -2 is the freshly added vertex before ages advance.
    new = -2
    def rev(age):
        return 0 <= age < 3 and bool(a & (1 << age))
    def finish(pieces):
        def age(i):
            return 0 if i == new else i+1 if 0 <= i < 2 else -1
        w = tuple(sorted((age(h),age(t)) for h,t in pieces))
        if sum(h == -1 for h,t in w) <= 1:
            yield w
    yield from finish(state+((new,new),))
    for i,(h,t) in enumerate(state):
        if not rev(t):
            yield from finish(state[:i]+((h,new),)+state[i+1:])
            for j,(hh,tt) in enumerate(state):
                if j != i and rev(hh):
                    yield from finish(tuple(p for k,p in enumerate(state) if k not in (i,j))+((h,tt),))
        if rev(h):
            yield from finish(state[:i]+((new,t),)+state[i+1:])


def raw_system():
    states=[()]
    index={():0}
    transitions=[{} for _ in range(8)]
    for u in states:
        for a in range(8):
            for w,c in Counter(fragments(u,a)).items():
                if w not in index:
                    index[w]=len(states)
                    states.append(w)
                transitions[a][index[w],index[u]]=c
    n=len(states)
    assert n==119
    M=[]
    for T in transitions:
        A=np.zeros((n,n),dtype=np.int64)
        for (i,j),c in T.items(): A[i,j]=c
        M.append(A)
    h=np.array([int(len(s)==1) for s in states],dtype=np.int64)
    return M,h


def multiply(ops,n):
    X=np.eye(n,dtype=np.int64)
    for A in ops: X=A@X
    return X


def rank_gte_two(A):
    ix=np.argwhere(A != 0)
    if not len(ix): return False
    i,j=ix[0]
    return bool(np.any(A[i,j]*A != A[:,j:j+1]*A[i:i+1,:]))


def transfer():
    M,h=raw_system()
    n=len(h)
    unitid=lambda w: sum(1<<(int(c)-1) for c in w if c!='0')
    U=[multiply((M[unitid(w)] for w in name.split('.')),n) for name in NAMES]
    init=M[0][:,0]
    assert sum(init)==1 and init[1]==1
    P=np.column_stack([init]+[U[i]@init for i in PREF])
    Q=np.vstack([h]+[h@U[i] for i in SUFF])
    bound=max(int(np.max(A.sum(axis=0))) for A in U)
    assert bound**3 <= np.iinfo(np.int64).max
    X=Matrix(Q@P)
    assert X.det() != 0
    Xi=X.inv()
    R=[]
    for A in U:
        twice=2*Matrix(Q@A@P)*Xi
        assert all(x.q == 1 for x in twice)
        B=np.array(twice.tolist(),dtype=object)
        assert not np.any(A[0])
        assert (B@Q[:,1:] == 2*(Q@A)[:,1:]).all()
        R.append(B)
    return Q,U,R,init,bound


def replay(path):
    Q,U,R,init,bound=transfer()
    n=Q.shape[1]
    G=json.loads(Path(path).read_text())
    p,q,den=G['n'],G['d'],G['den']
    assert p>q>0 and den>0
    I={tuple(w):np.array(vals,dtype=object) for w,vals in G['phi']}
    assert all(len(v)==8 for v in I.values())
    # In the fixed-memory certificate every tag has length three. The
    # adaptive certificate remembers one extra predecessor precisely when
    # the final triple is listed in the promotion table.
    k=3
    promote={tuple(w) for w in G.get('promote', [])}
    assert all(len(w)==k for w in promote)
    max_k=k+bool(promote)
    assert all(len(w)==(k+1 if w[-k:] in promote else k) for w in I)
    # All raw products have at most max_k+2 unit factors. Their columns
    # have l1 norm at most bound at each step. Rule out int64 overflow.
    assert bound**(max_k+2) <= np.iinfo(np.int64).max
    alpha=tuple(sorted({a for w in I for a in w}))
    assert alpha and all(0<=a<len(U) for a in alpha)
    prefixes={():np.eye(8,dtype=object)}
    for _ in range(k):
        prefixes={w+(j,):B for w,A in prefixes.items() for j in alpha if rank_gte_two(B:=R[j]@A)}
    assert promote<=set(prefixes)
    for w in promote:
        A=prefixes.pop(w)
        for j in alpha:
            B=A@R[j]
            if rank_gte_two(B):
                prefixes[(j,)+w]=B
    assert set(prefixes)==set(I)
    edges=total=0
    for w,A in prefixes.items():
        starting=multiply((U[j] for j in w),n)
        J=Q@starting
        V=I[w]@J
        assert (2000*V>=den*J[0]).all() and (V<=2*den*J[0]).all(), ('compare',w)
        for j in alpha:
            if not rank_gte_two(R[j]@A): continue
            tail=w+(j,)
            t=tail[-k-1:] if tail[-k:] in promote else tail[-k:]
            assert t in I
            Y=I[t]@(Q@U[j]@starting)
            cost=len(NAMES[j].split('.'))
            assert (q**cost*Y>=p**cost*V).all(), ('expand',w,j)
            edges+=1; total+=n
    print('exact suffix cone passed:', 'units=',tuple(NAMES[j] for j in alpha),
          'rate=',p/q, 'states=',len(I),'edges=',edges,'coordinate inequalities=',total)


def obstruction(path):
    # Exact Farkas witness: a positive sum of actual depth-three inequalities
    # cancels every coefficient of every functional, but has negative RHS.
    Z,U,R,_,bound=transfer()
    n=Z.shape[1]
    assert bound**5 <= np.iinfo(np.int64).max
    data=json.loads(Path(path).read_text())
    total={}
    rhs=Fraction(0)
    rate=Fraction(17699,10000)
    for row in data:
        w=tuple(row['w'])
        assert len(w)==3
        P=multiply((U[j] for j in w),n)
        V=Z@P
        assert rank_gte_two(multiply((R[j] for j in w),8))
        a=row['a']
        weight=Fraction(row['weight'])
        assert weight>0
        if 't' not in row:
            assert any(tuple(-v)==tuple(a) for v in V.T)
            assert a[0]<0
            rhs+=weight*Fraction(a[0],2000)
            pairs=[(w,a)]
        else:
            t=tuple(row['t']); j=t[-1]
            assert t==w[1:]+(j,)
            assert rank_gte_two(multiply((R[i] for i in w+(j,)),8))
            Y=Z@U[j]@P
            b=row['b']
            assert any(tuple(v)==tuple(a) and tuple(-y)==tuple(b) for v,y in zip(V.T,Y.T))
            f=rate**len(NAMES[j].split('.'))
            pairs=[(w,[f*x for x in a]),(t,b)]
        for tag,vals in pairs:
            for i,x in enumerate(vals):
                total[tag,i]=total.get((tag,i),Fraction(0))+weight*x
    assert rhs<0 and all(x==0 for x in total.values())
    print('exact depth-three Farkas obstruction passed:',len(data),'rows')


if __name__=='__main__':
    if len(sys.argv)==3 and sys.argv[1]=='--obstruction':
        obstruction(sys.argv[2])
    else:
        replay(sys.argv[1])
