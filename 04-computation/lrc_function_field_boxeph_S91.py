#!/usr/bin/env python3
"""
The function-field Lonely Runner: setting it up honestly  (boxeph-2026-07-18-S91)
=================================================================================
Test whether the ultrametric over F_p[t] makes the packing / INV UNCONDITIONAL.
Polynomials over F_p as coefficient tuples (low->high, no trailing zeros).
Distance ||v * a/Q|| = p^(deg((v*a) mod Q) - deg Q)  (higher deg remainder = farther).
We work in LOG units: logdist(v; a,Q) = deg((v*a)%Q) - deg Q  (in {-degQ,...,-1}, or -inf).
M(V) = max over rational times a/Q of min_v logdist  (the max-min; -1 is the farthest possible).
Goal: find near-tight families, inspect the residues at the maximizer, test 'core = subspace/AP'.
"""
from itertools import product
from functools import reduce
from math import gcd as igcd

P = 3  # base field F_P

def norm(a):  # strip trailing zeros
    a=list(a)
    while a and a[-1]%P==0: a.pop()
    return tuple(x%P for x in a)
def deg(a):
    a=norm(a); return len(a)-1  # deg 0 poly -> -1 (the zero poly)
def padd(a,b):
    n=max(len(a),len(b)); return norm(tuple((a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(n)))
def psub(a,b):
    n=max(len(a),len(b)); return norm(tuple((a[i] if i<len(a) else 0)-(b[i] if i<len(b) else 0) for i in range(n)))
def pmul(a,b):
    a=norm(a); b=norm(b)
    if not a or not b: return ()
    r=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): r[i+j]+=x*y
    return norm(r)
def _degl(a):  # degree of a coefficient list (may have trailing zeros); -1 if zero
    i=len(a)-1
    while i>=0 and a[i]%P==0: i-=1
    return i
def pdivmod(a,b):
    a=list(norm(a)); b=norm(b); assert b, "div by 0"
    lb=len(b); q=[0]*(max(0,len(a)-lb+1))
    inv=pow(b[-1],P-2,P)
    while _degl(a)>=lb-1:
        da=_degl(a); d=da-(lb-1); c=(a[da]*inv)%P; q[d]=c
        for i,y in enumerate(b): a[d+i]=(a[d+i]-c*y)%P
    return norm(q), norm(a)
def pmod(a,b): return pdivmod(a,b)[1]
def pgcd(a,b):
    a=norm(a); b=norm(b)
    while b: a,b=b,pmod(a,b)
    return a
def monic(a):
    a=norm(a)
    if not a: return a
    inv=pow(a[-1],P-2,P); return norm(tuple((x*inv)%P for x in a))

def polys_upto_deg(d):
    out=[]
    for length in range(1,d+2):
        for coeffs in product(range(P),repeat=length):
            a=norm(coeffs)
            if a and len(a)==length: out.append(a)
    # dedup
    seen=set(); res=[]
    for a in out:
        if a not in seen: seen.add(a); res.append(a)
    return res

def M_logdist(V, maxdegQ=3):
    """max over a/Q (deg a<deg Q, gcd(a,Q)=1, Q monic) of min_v logdist; returns (val, a, Q, residues)."""
    best=None; bestinfo=None
    Qs=[q for q in polys_upto_deg(maxdegQ) if deg(q)>=1 and q==monic(q)]
    for Q in Qs:
        dQ=deg(Q)
        for a in polys_upto_deg(dQ-1)+[()]:
            if a==(): continue
            if deg(a)>=dQ: continue
            if pgcd(a,Q)!=(1,): continue
            # min over runners of deg((v a) mod Q) - dQ ; if any remainder 0 -> collision (skip: min=-inf)
            worst=0; ok=True; res={}
            for v in V:
                r=pmod(pmul(v,a),Q); res[v]=r
                if not r: ok=False; break
                worst=min(worst, deg(r)-dQ)
            if not ok: continue
            if best is None or worst>best:
                best=worst; bestinfo=(a,Q,dict(res))
    return best, bestinfo

if __name__=="__main__":
    print(f"Function-field Lonely Runner over F_{P}[t].  logdist -1 = farthest (dist 1/{P}).")
    print("="*80)
    one=(1,); t=(0,1)
    # candidate families: subspaces (difference-closed = F_p-subspace) + a 'far' element
    fams = {
        "F3-constants {1,2}":            [(1,),(2,)],
        "deg<=1 nonzero (subspace-ish)": [p for p in polys_upto_deg(1)],   # 1,2,t,t+1,2t,... 8 elts
        "{1,2} + far t^2":               [(1,),(2,),(0,0,1)],
        "monomials {1,t,t^2}":           [(1,),(0,1),(0,0,1)],
        "{1, t, t+1} (deg<=1 partial)":  [(1,),(0,1),(1,1)],
    }
    for name,V in fams.items():
        val,info=M_logdist(V, maxdegQ=3)
        if info is None:
            print(f"  {name}: no lonely time found (all collide?)"); continue
        a,Q,res=info
        print(f"  {name}  (|V|={len(V)}): M(log)={val}  at a={a} Q={Q}")
        print(f"     residues (v: (v*a mod Q)): "+", ".join(f"{v}:{res[v]}" for v in V))
    print()
    print("READING: inspect whether near-tight families have their 'core' residues forced to fill")
    print("a subspace / coset (the F_p[t] analog of the AP {1..12}=F_13^*), UNCONDITIONALLY.")

# ============ THE FUNCTION-FIELD DEEP WELL and the UNCONDITIONAL packing ============
def vanishing_poly():  # t^P - t = prod_{c in F_P}(t-c), covers all points
    r=[0]*(P+1); r[P]=1; r[1]=(r[1]-1)%P; return norm(r)
def evalp(a,c):
    s=0
    for x in reversed(norm(a)): s=(s*c+x)%P
    return s
def roots(a): return [c for c in range(P) if evalp(a,c)==0]
def is_covering_ff(V):  # roots of speeds cover all of F_P
    covered=set()
    for v in V: covered|=set(roots(v))
    return covered==set(range(P))

def deep_well_ff():
    core=[(c,) for c in range(1,P)]  # F_P^* as constants = the 'AP' core (p-1 of them)
    return core+[vanishing_poly()]

if __name__=="__main__" or True:
    print("\n"+"="*80)
    print(f"THE FUNCTION-FIELD DEEP WELL over F_{P}[t]:  V = F_{P}^* (constants) + (t^{P}-t)")
    print("="*80)
    V=deep_well_ff()
    print(f"  core = F_{P}^* = {[c[0] for c in V[:-1]]} (constants);  killer = t^{P}-t = {vanishing_poly()}")
    print(f"  covering (roots cover F_{P})? {is_covering_ff(V)}")
    # no level-1 (deg Q=1) lonely time since covering => must use deg Q=2
    val1,info1=M_logdist(V, maxdegQ=1)
    print(f"  best at deg Q<=1 (level 1): {val1}  (None/-inf => no level-1 loneliness, = covering)")
    val,info=M_logdist(V, maxdegQ=2)
    if info:
        a,Q,res=info
        print(f"  M(log) at deg Q<=2: {val}  at a={a}, Q={Q} (deg {deg(Q)})")
        # top coefficient (degree deg Q - 1) of each core residue
        dQ=deg(Q)
        tops=[]
        for v in V[:-1]:  # core
            r=res[v]; top = r[dQ-1] if len(r)>=dQ else 0
            tops.append(top)
        killer_r=res[V[-1]]; killer_top=killer_r[dQ-1] if len(killer_r)>=dQ else 0
        print(f"  CORE residue top-coeffs (deg {dQ-1}): {sorted(tops)}  == F_{P}^*={list(range(1,P))}? {sorted(tops)==list(range(1,P))}")
        print(f"  killer residue={killer_r}, top-coeff={killer_top}  (0 => killer is the 'anomaly', collides in top-coeff)")
    print()
    print("UNCONDITIONAL PACKING (the point): p speeds at a level-2 lonely time have top-coeffs in")
    print(f"F_{P}^* (only {P-1} values), so the {P-1}-element difference-closed core FILLS F_{P}^* exactly.")
    print("No packing slack, no sub-case B, no q=13val+1 ambiguity -- the ultrametric makes it EXACT.")
