"""Automorphism groups of the unit-circle Cayley graphs G_q = Cay(F_{q^2}, mu_{q+1}) via nauty (dreadnaut).
Affine reference group: F_{q^2} x| (mu_{q+1} x| Gal(F_{q^2}/F_p)), order q^2 (q+1) 2e, q = p^e.
Session mac-mini-2026-10-12-logic4.
"""
import itertools, subprocess, sys
from math import gcd
def is_prime(n): return n>1 and all(n%k for k in range(2,int(n**.5)+1))
def factor_pp(q):
    for p in range(2,q+1):
        if q%p==0:
            e=0
            while q%p==0: q//=p; e+=1
            return (p,e) if q==1 else None
def polymul(f,g,p):
    r=[0]*(len(f)+len(g)-1)
    for i,fi in enumerate(f):
        if fi:
            for j,gj in enumerate(g): r[i+j]=(r[i+j]+fi*gj)%p
    return r
def polymod(f,m,p):
    f=f[:]; d=len(m)-1; inv=pow(m[-1],-1,p)
    while len(f)-1>=d and any(f):
        if f[-1]==0: f.pop(); continue
        coef=f[-1]*inv%p; sh=len(f)-1-d
        for i in range(len(m)): f[sh+i]=(f[sh+i]-coef*m[i])%p
        f.pop()
    while len(f)<d: f.append(0)
    return f[:d]
def irreducible(m,p):
    # brute: m monic degree d over F_p irreducible iff no root-extension factor; use Rabin-free brute force: check gcd with x^{p^i}-x not needed for d<=10 small p: brute by trial division
    d=len(m)-1
    for e in range(1,d//2+1):
        for coeffs in itertools.product(range(p),repeat=e):
            g=list(coeffs)+[1]
            # division check
            f=m[:]
            # polynomial remainder of f by g
            r=f[:]
            while len(r)-1>=e:
                if r[-1]==0: r.pop(); continue
                coef=r[-1]; sh=len(r)-1-e
                for i in range(len(g)): r[sh+i]=(r[sh+i]-coef*g[i])%p
                r.pop()
            if not any(r): return False
    return True
def field(p,d):
    """F_{p^d} as tuples; returns (elements list, mul, add, pow)."""
    for coeffs in itertools.product(range(p),repeat=d):
        m=list(coeffs)+[1]
        if irreducible(m,p): break
    elems=[tuple(t) for t in itertools.product(range(p),repeat=d)]
    idx={e:i for i,e in enumerate(elems)}
    def mul(f,g): return tuple(polymod(polymul(list(f),list(g),p),m,p))
    def add(f,g): return tuple((fi+gi)%p for fi,gi in zip(f,g))
    def pw(f,n):
        r=tuple([1]+[0]*(d-1)); base=f
        while n:
            if n&1: r=mul(r,base)
            base=mul(base,base); n>>=1
        return r
    return elems,idx,mul,add,pw
def gq(q):
    p,e=factor_pp(q)
    elems,idx,mul,add,pw=field(p,2*e)
    one=tuple([1]+[0]*(2*e-1))
    mu=[t for t in elems if pw(t,q+1)==one]
    assert len(mu)==q+1
    adj={i:set() for i in range(len(elems))}
    for i,t in enumerate(elems):
        for m_ in mu:
            adj[i].add(idx[add(t,m_)])
    return adj,p,e,mu,elems,idx,add
def nauty_grpsize(adj):
    n=len(adj)
    lines=["n=%d g"%n]
    for i in range(n):
        lines.append(" ".join(str(j) for j in sorted(adj[i]) if j>i)+(";" if i<n-1 else "."))
    inp="\n".join(lines)+"\nx o\nq\n"
    out=subprocess.run(["dreadnaut"],input=inp,capture_output=True,text=True).stdout
    import re
    m=re.findall(r"grpsize=([0-9.e+]+)",out)
    orb=re.findall(r"(\d+) orbits?;",out)
    return (m[-1] if m else None), (orb[-1] if orb else None), out
def nauty_stab_orbits(adj,v0):
    """orbits of the stabiliser of vertex v0 (dreadnaut: fix v0 via partition)"""
    n=len(adj)
    lines=["n=%d g"%n]
    for i in range(n):
        lines.append(" ".join(str(j) for j in sorted(adj[i]) if j>i)+(";" if i<n-1 else "."))
    inp="\n".join(lines)+"\nf=[%d]\nx o\nq\n"%v0
    out=subprocess.run(["dreadnaut"],input=inp,capture_output=True,text=True).stdout
    import re
    m=re.findall(r"grpsize=([0-9.e+]+)",out); orb=re.findall(r"(\d+) orbits?;",out)
    return (m[-1] if m else None), (orb[-1] if orb else None)
qs=[q for q in range(2,82) if factor_pp(q)]
print("q  p^e   |V|  deg  |Aut|        |Aff|=q^2(q+1)2e  ratio   orbits(Aut)  |Stab_0|  orbits(Stab_0)")
for q in qs:
    adj,p,e,mu,elems,idx,add=gq(q)
    g,orb,out=nauty_grpsize(adj)
    aff=q*q*(q+1)*2*e
    gs,orbs=nauty_stab_orbits(adj,0)
    G=float(g.replace('e','e')) if g else float('nan')
    print("%-3d %d^%d  %4d  %3d  %-12s %-15d %-7.4g %-12s %-9s %s"%(q,p,e,q*q,q+1,g,aff,G/aff,orb,gs,orbs))
    sys.stdout.flush()
