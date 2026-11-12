"""Search for invariant affine lines of the O/Z/J linear semigroup: a direction p that is an eigenvector
(integer eigenvalue >= 2) of several words."""
import itertools, sys, collections
import sympy as sp
from semigroup import GEN
L=int(sys.argv[1])
x=sp.symbols('x')
byp=collections.defaultdict(list)
def norm(p):
    from math import gcd
    v=[sp.Rational(t) for t in p]
    den=sp.ilcm(*[t.q for t in v]); v=[int(t*den) for t in v]
    g=0
    for t in v: g=gcd(g,abs(t))
    if g==0: return None
    v=[t//g for t in v]
    for t in v:
        if t!=0:
            if t<0: v=[-s for s in v]
            break
    return tuple(v)
nw=0
for n in range(1,L+1):
    for w in itertools.product('10r',repeat=n):
        W=sp.eye(4)
        for ch in w: W=GEN[ch]*W
        nw+=1
        cp=W.charpoly(x).as_expr()
        for lam,mult in sp.roots(sp.Poly(cp,x)).items():
            if lam.is_integer and lam>=2:
                ns=(W-lam*sp.eye(4)).nullspace()
                if len(ns)==1:
                    p=norm(list(ns[0]))
                    if p: byp[p].append((''.join(w),int(lam)))
print("words:",nw,"; eigen-directions with integer eigenvalue >= 2:",len(byp))
def primitive(w):
    n=len(w)
    for d in range(1,n):
        if n%d==0 and w==w[:d]*(n//d): return w[:d]
    return w
shared=[]
for p,ws in byp.items():
    roots={primitive(w) for w,_ in ws}
    if len(roots)>=2: shared.append((p,ws,roots))
print("directions shared by words with at least two distinct primitive roots:",len(shared))
for p,ws,roots in sorted(shared,key=lambda t:-len(t[2]))[:25]:
    print(" p=%s  roots=%d  examples: %s" % (p,len(roots),sorted(ws,key=lambda t:(len(t[0]),t[0]))[:8]))
