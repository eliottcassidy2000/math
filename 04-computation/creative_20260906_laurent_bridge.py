#!/usr/bin/env python3
"""Exact coupled-window observer gates; no repository producer imports."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement
from math import factorial

GATES = 0
TRACE = sha256()

def gate(ok, label, data=None):
    global GATES
    if not ok:
        raise RuntimeError(f"FAILED {label}: {data}")
    GATES += 1
    TRACE.update((label + ':' + repr(data) + '\n').encode())

def mul(a,b):
    out = [F(0)] * (len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): out[i+j] += x*y
    return out

def deriv(a, r):
    for _ in range(r): a = [i*a[i] for i in range(1,len(a))]
    return a

def falling(a,r):
    if a < r or a < 0: return 0
    return factorial(a)//factorial(a-r)

def window(a,r,s):
    # Fixed ambient reversal retains zero endpoint slots.
    return deriv(list(reversed(deriv(list(reversed(a)),s))),r)

def coeff(a,j): return a[j] if 0<=j<len(a) else F(0)

def quadratic(a,k,r,s):
    n=len(a)-1
    return sum((falling(j,r)*falling(2*k-j,r)*falling(n-j,s)*
                falling(n-2*k+j,s)*coeff(a,j)*coeff(a,2*k-j)
                for j in range(n+1)),F(0))

def verify(a,k,label):
    n=len(a)-1
    gate(a[0] != 0 and a[-1] != 0 and a[k] == 0,'carrier',label)
    gate(a[k-1]*a[k+1]<0,'neighbor-sign',label)
    rows=[]
    for r in range(k):
        for s in range(n-k):
            b=window(a,r,s)
            q=quadratic(a,k,r,s)
            gate(coeff(b,k-r)==0,'preserved-zero',(label,r,s))
            gate(coeff(mul(b,b),2*(k-r))==q,'independent-convolution',(label,r,s,q))
            gate(q<0,'negative-window',(label,r,s,q))
            for j in range(n+1):
                d=j-k
                w=1
                for i in range(r): w*=((k-i)**2-d*d)
                for i in range(s): w*=((n-k-i)**2-d*d)
                wf=(falling(j,r)*falling(2*k-j,r)*falling(n-j,s)*falling(n-2*k+j,s))
                # Outside the actual paired support the coefficient is zero.
                gate((w-wf)*coeff(a,j)*coeff(a,2*k-j)==0,'centered-kernel',(label,r,s,j))
            rows.append(q)
    gate(sum((i+1)*q for i,q in enumerate(rows))<0,'positive-cone',label)
    # The next derivative reaches a coefficient boundary; strictness ends.
    gate(coeff(mul(window(a,k,0),window(a,k,0)),0)==0,'left-boundary',label)
    right=window(a,0,n-k)
    gate(coeff(mul(right,right),2*k)==0,'right-boundary',label)
    return rows

cases=0
for n in range(3,9):
    for roots in combinations_with_replacement((-3,-1,1,2),n-1):
        p=[F(1)]
        for root in roots: p=mul(p,[F(1),F(root)])
        for k in range(1,n):
            if not p[k-1] or not p[k]: continue
            last=-p[k]/p[k-1]
            a=mul(p,[F(1),last])
            verify(a,k,(n,roots,k,str(last)))
            cases+=1

# A real-rooted carrier with an individually positive opposite product.
a=[F(1),F(-22,5),F(0),F(54,5),F(27,5)]
rows=verify(a,2,'mixed-sign-hostile')
gate(a[0]*a[4]>0,'positive-outer-pair',a)
print('MIXED_PAIR carrier=',','.join(map(str,a)),'window_values=',','.join(map(str,rows)))

# Full-support path control: F=(1+u)^4, G=4+u, rho=-1.
a=[F(1)]
for _ in range(4): a=mul(a,[F(1),F(1)])
a=mul(a,[F(-1),F(4)])
rows=verify(a,1,'path-common-carrier')
print('PATH_CONTROL support=(-3,1,9) rho=-1 carrier=',','.join(map(str,a)),
      'window_values=',','.join(map(str,rows)))

# Preserving the zero alone is insufficient: Euler shift destroys real roots.
a=[F(1),F(0),F(-1)]
k=1
bad=[(j-k)*c for j,c in enumerate(a)]
gate(bad[k]==0 and coeff(mul(bad,bad),2*k)==2,'zero-preserving-Euler-hostile',bad)
# Real-rootedness cannot be dropped.
bad=[F(1),F(0),F(1)]
gate(coeff(mul(bad,bad),2)==2,'non-real-core-hostile',bad)
print('HOSTILES arbitrary_zero_preserving_Euler_square=2 nonreal_core_square=2')

def basis(M,t,d):
    out=1
    for a in range(t): out*=((M-a)**2-d*d)
    return out

def cone_expansion(M,N,r,s):
    c=[0]*M
    c[r]=1
    for b in range(s):
        cc=[0]*M
        for t,v in enumerate(c):
            if not v: continue
            gate(t>=max(r,M-N+b),'cone-support-invariant',(M,N,r,s,b,t))
            multiplier=(N-b)**2-(M-t)**2
            gate(multiplier>=0,'cone-nonnegative-step',(M,N,r,s,b,t))
            cc[t]+=v*multiplier
            if t+1<M: cc[t+1]+=v
        c=cc
    return c

def back_substitute(M,values):
    c=[]
    for t in range(M):
        d=M-t
        c.append((values[d]-sum(c[i]*basis(M,i,d) for i in range(t)))/F(basis(M,t,d)))
    return c

cone_cases=0
for M in range(1,9):
    for N in range(M,21):
        for r in range(M):
            for s in range(N):
                c=cone_expansion(M,N,r,s)
                gate(all(v>=0 for v in c) and any(c),'cone-coefficients',(M,N,r,s,c))
                vals={d:basis(M,r,d)*basis(N,s,d) for d in range(1,M+1)}
                for d,value in vals.items():
                    gate(value==sum(v*basis(M,t,d) for t,v in enumerate(c)),
                         'cone-kernel-identity',(M,N,r,s,d))
                gate(back_substitute(M,vals)==c,'cone-independent-triangular-inverse',(M,N,r,s))
                cone_cases+=1
outside=back_substitute(2,{1:F(0),2:F(1)})
gate(outside==[F(1),F(-1,3)],'nonnegative-kernel-outside-cone',outside)
print('CONE universe=1<=M<=8,M<=N<=20,0<=r<M,0<=s<N cases=',cone_cases)
print('OUTSIDE_CONE W(1)=0 W(2)=1 coefficients=',','.join(map(str,outside)),
      'mixed_pair_response=54/5')

print(f'UNIVERSE degree=3..8 roots_prefix_multisets=(-3,-1,1,2) tuned_nonzero_last_root cases={cases}')
print(f'PASS explicit_gates={GATES} semantic_sha256={TRACE.hexdigest()}')
