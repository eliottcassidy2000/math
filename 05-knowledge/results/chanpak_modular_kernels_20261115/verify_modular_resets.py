#!/usr/bin/env python3
"""Standalone exact replay for characteristic polynomials, modular resets and strong return blocks."""
from collections import Counter

import numpy as np
from sympy import Matrix

NAMES = '0.2 1.2 0.123 2.13 0.1.3 0.3.1 0.13.2 0.1.123 0.3.12 0.12.3 1.2.3 1.12.3 0.1.23'.split()
PREF = (1,0,2,3,4,7,9)
SUFF = (0,1,2,3,5,6,7)


def fragments(state, a, old_heads=1):
    # State lists ordered (head,tail) ages of the current directed paths.
    # -1 is old; -2 is the freshly added vertex before ages advance.
    new = -2
    def rev(age):
        return 0 <= age < 3 and bool(a & (1 << age))
    def finish(pieces):
        def age(i):
            return 0 if i == new else i+1 if 0 <= i < 2 else -1
        w = tuple(sorted((age(h),age(t)) for h,t in pieces))
        if sum(h == -1 for h,t in w) <= old_heads:
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


def rows(word):
    out=[0]*len(word)
    for j,b in enumerate(word):
        for i in range(j):
            if 1 <= j-i <= 3 and str(j-i) in b: out[j]|=1<<i
            else: out[i]|=1<<j
    return out

def is_strong(out):
    n=len(out)
    rev=[sum(((out[j]>>i)&1)<<j for j in range(n)) for i in range(n)]
    for arcs in (out,rev):
        s=1
        while True:
            t=s
            for i,v in enumerate(arcs):
                if s>>i&1: t|=v
            if t==s: break
            s=t
        if t!=(1<<n)-1: return False
    return True

def ham(out):
    dp=[{} for _ in range(1<<len(out))]
    top=len(dp)-1
    for i in range(len(out)):dp[1<<i][i]=1
    for s,d in enumerate(dp):
        for i,c in d.items():
            b=out[i]&~s
            while b:
                j=b&-b; b^=j
                D=dp[s|j]; k=j.bit_length()-1
                D[k]=D.get(k,0)+c
    return sum(dp[top].values())

def literal(w):
    return rows(['0']+[a for j in w for a in NAMES[j].split('.')])

def check_two_ports():
    from functools import lru_cache
    from collections import defaultdict
    @lru_cache(None)
    def step(s, k):
        return tuple(Counter(fragments(s,k,2)).items())
    S=[()]; seen={()}
    for s in S:
        for k in range(8):
            for t,_ in step(s,k):
                if t not in seen: S.append(t); seen.add(t)
    n=len(S)
    assert n==206
    tests=[
      (3,(7,9),[1,0]),
      (5,(1,5),[1,3,3,4]),
      (7,(0,9,3),[1,5,2,1]),
      (11,(12,9,8),[1,1,3,8]),
      (13,(0,4,3),[1,8,7,12]),
    ]
    for p,w,coeff in tests:
        v={((0,0),):1}
        z=[]
        seq=[sum(1<<(int(c)-1) for c in a if c!='0')
             for j in w for a in NAMES[j].split('.')]
        d=len(coeff)-1
        for i in range(n+d):
            z.append(2*sum(c for s,c in v.items() if len(s) in (1,2))%p)
            for k in seq:
                out=defaultdict(int)
                for s,c in v.items():
                    for t,b in step(s,k): out[t]+=b*c
                v={s:x%p for s,x in out.items() if x%p}
        assert z[0]==2
        # The 206 residuals imply every future one by Cayley-Hamilton
        # on the 206-dimensional linear period transfer.
        assert all(sum(c*z[i-j] for j,c in enumerate(coeff))%p==0
                   for i in range(d,n+d))
        print('Pi2 mod',p,'recurrence',coeff,'begins',z[:8])


def check():
    pds=[
      (3,(7,9),39,-30,4,(4,5)),
      (5,(1,5),25,10,3,(0,)),
      (7,(0,9,3),161,42,3,(2,)),
      (11,(12,9,8),473,550,3,(4,)),
      (13,(0,4,3),117,-52,3,(2,)),
    ]
    Q,U,Rs,v,_=transfer()
    Rs=[Matrix(A.tolist())/2 for A in Rs]
    v=Matrix((Q@v).tolist())
    def evalm(w):
        A=Matrix.eye(8)
        for j in w: A=Rs[j]*A
        return A
    def pzero(A,p):
        return all(x.p%p==0 and x.q%p!=0 for x in A)
    for p,w,t,d,k,alt in pds:
        A=evalm(w)
        assert all(x.q in (1,2) for x in A)
        assert A.charpoly().all_coeffs()==[1,-t,d]+[0]*6
        assert t%p==0 and d%p==0 and d%(p*p)!=0
        Z0=w*k; Z1=alt+Z0[len(alt):]
        assert len(literal(Z0))==len(literal(Z1))
        T=A**k; B=evalm(Z1)
        E=Matrix.hstack(*T.columnspace())
        assert E.cols==2 and pzero(T,p) and not pzero(A**(k-1),p)
        assert Matrix.hstack(E,B).rank()==2 and (B*E).rank()==2
        assert is_strong(literal(w*2))
        assert is_strong(literal(Z0)) and is_strong(literal(Z1))
        assert (T*v)[0]!=(B*v)[0]
        for j in (1,2) if len(literal(w*2))<=17 else (1,):
            assert (A**j*v)[0]==ham(literal(w*j))
        print(p,NAMES[w[0]],'period',w,'reset',k,
              'block length',len(literal(Z0))-1,'always reset',pzero(B,p),
              'first counts',int((T*v)[0]),int((B*v)[0]))
    print('Five exact modular resets and literal subset-DP controls passed.')

if __name__=='__main__':
    check()
    check_two_ports()
