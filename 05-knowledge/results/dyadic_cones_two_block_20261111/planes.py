"""Invariant planes P = span(p, q) of words W of the O/Z/J semigroup with W p = lam p (lam >= 2 integer),
W q = q + b p: on the affine lines beta q + Z p the word acts as alpha -> lam alpha + b beta.
Group words by the plane; report planes shared by words with different primitive roots, with the affine maps,
and the reachable states (short words) lying in the plane."""
import itertools, sys, collections
import sympy as sp
from semigroup import GEN
L=int(sys.argv[1])
x=sp.symbols('x')
def primitive(w):
    n=len(w)
    for d in range(1,n):
        if n%d==0 and w==w[:d]*(n//d): return w[:d]
    return w
def plucker(p,q):
    M=sp.Matrix([list(p),list(q)])
    co=[]
    for i in range(4):
        for j in range(i+1,4):
            co.append(M[0,i]*M[1,j]-M[0,j]*M[1,i])
    den=sp.ilcm(*[sp.Rational(c).q for c in co]); co=[int(sp.Rational(c)*den) for c in co]
    g=0
    from math import gcd
    for c in co: g=gcd(g,abs(c))
    if g==0: return None
    co=[c//g for c in co]
    for c in co:
        if c!=0:
            if c<0: co=[-t for t in co]
            break
    return tuple(co)
planes=collections.defaultdict(list)
words=[]
for n in range(1,L+1):
    for w in itertools.product('10r',repeat=n):
        W=sp.eye(4)
        for ch in w: W=GEN[ch]*W
        words.append((''.join(w),W))
for name,W in words:
    if W.rank()<2: continue
    roots=sp.roots(W.charpoly(x))
    if 1 not in roots: continue
    for lam in roots:
        if not (lam.is_integer and lam>=2): continue
        ns=(W-lam*sp.eye(4)).nullspace()
        if len(ns)!=1: continue
        p=ns[0]
        # q with (W - I) q = b p : solve for (q, b) : [W - I | -p] (q; b) = 0
        A=(W-sp.eye(4)).row_join(-p)
        for sol in A.nullspace():
            q=sol[:4,0]; b=sol[4]
            if sp.Matrix([list(p),list(q)]).rank()<2: continue
            key=plucker(p,q)
            if key: planes[key].append((name,int(lam),p,q,b))
print("words:",len(words),"; invariant planes of the required form:",len(planes))
# reachable states
states={}
v0=sp.Matrix([1,0,1,1])
for name,W in [('',sp.eye(4))]+words:
    states[name]=W*v0
shared=[(k,v) for k,v in planes.items() if len({primitive(t[0]) for t in v})>=2]
print("planes shared by words with at least two primitive roots:",len(shared))
for key,v in sorted(shared,key=lambda t:-len(t[1]))[:12]:
    p=v[0][2]; q=v[0][3]
    P=sp.Matrix.hstack(p,q)
    inplane=[(nm,list(s)) for nm,s in states.items() if sp.Matrix.hstack(P,s).rank()==2 and len(nm)<=5]
    print("plane",key,"words:",sorted({(t[0],t[1]) for t in v},key=lambda t:(len(t[0]),t[0]))[:10])
    print("    reachable states in the plane (words of length <= 5):",inplane[:6])
