#!/usr/bin/env python3
r"""
lrc_chic_gw_quasiperiodic_kps_S76.py  (kind-pasteur-2026-07-07-S76)

chi_c(G_GW) in (13,14]?  Search QUASI-PERIODIC circular colorings of G(Z,GW):
    C(x) = c(x mod T) + (x div T)*Delta   (mod p),   c: {0..T-1} -> Z_p.
This is an INFINITE coloring of G(Z,GW) (quasi-periodicity => all constraints reduce to
x in {0..T-1}).  The ROTATION coloring is the T=1 case (gives exactly 1/M = 14); a T>=2
quasi-periodic coloring with p/q < 14 would prove chi_c(G_GW) < 14 -- i.e. LINEARIZATION
FAILS at the GW tight instance (the circular rung is blind, like chi_f=13).

A valid (p,q)-coloring here (min edge-gap >= q, 13 < p/q < 14) => chi_c(G_GW) <= p/q < 14.
Tiny SAT per (p,q,T,Delta): T variables in Z_p, T*13 gap constraints.  Sweep T, Delta.
"""
from math import gcd
from pysat.solvers import Cadical153 as Sat
from pysat.card import CardEnc, EncType

GW = [1,2,3,4,5,6,7,8,9,10,11,13,24]

def circ_dist(a, b, p):
    d = abs(a - b) % p
    return min(d, p - d)

def quasi_color(p, q, T, Delta):
    """SAT: exists c:{0..T-1}->Z_p with C(x)=c(x%T)+(x//T)*Delta a (p,q)-coloring of G(Z,GW)?
    Returns the color list c or None."""
    nv = T * p
    def V(r, a): return r * p + a + 1
    cls = []
    for r in range(T):
        cls.append([V(r, a) for a in range(p)])                 # >=1 color
        for a in range(p):
            for b in range(a+1, p):
                cls.append([-V(r, a), -V(r, b)])                # <=1 color
    # gap constraints: for base r in {0..T-1}, s in GW: circ_dist(C(r+s),C(r)) >= q
    for r in range(T):
        for s in GW:
            j = r + s
            carry = (j // T) * Delta                            # r//T = 0 since r<T
            jmod = j % T
            for a in range(p):          # c(r) = a  => C(r) = a
                for b in range(p):      # c(jmod) = b => C(j) = b + carry
                    if circ_dist((b + carry) % p, a, p) < q:
                        cls.append([-V(r, a), -V(jmod, b)])
    s = Sat(bootstrap_with=cls)
    ok = s.solve()
    col = None
    if ok:
        m = set(v for v in s.get_model() if v > 0)
        col = []
        for r in range(T):
            for a in range(p):
                if V(r, a) in m: col.append(a); break
    s.delete()
    return col

def verify(col, p, q, T, Delta, span=4000):
    def C(x): return (col[x % T] + (x // T) * Delta) % p
    return all(circ_dist(C(x+s), C(x), p) >= q for x in range(span) for s in GW)

# candidate p/q in (13,14), ascending p (small = fast)
CANDS = []
for q in range(2, 9):
    for p in range(13*q+1, 14*q):
        if gcd(p, q) == 1:
            CANDS.append((p, q, p/q))
CANDS.sort(key=lambda t: (t[0], t[1]))

print("="*82)
print("chi_c(G_GW): quasi-periodic (p,q)-coloring search, 13 < p/q < 14")
print("  rotation (T=1) gives 14; any T>=2 hit with p/q<14 => chi_c < 14 (linearization FAILS)")
print("="*82)
found = None
for (p, q, r) in CANDS:
    hit = None
    for T in range(2, 27):
        for Delta in range(1, p):
            if gcd(Delta, p) == 0: continue
            col = quasi_color(p, q, T, Delta)
            if col is not None:
                if verify(col, p, q, T, Delta):
                    hit = (T, Delta, col); break
        if hit: break
    if hit:
        T, Delta, col = hit
        print(f"  p/q={p}/{q}={r:.4f}: *** SAT ***  T={T} Delta={Delta}  c={col}")
        found = (p, q, r, T, Delta, col); break
    else:
        print(f"  p/q={p}/{q}={r:.4f}: no quasi-periodic coloring (T=2..26)")
print()
if found:
    p,q,r,T,Delta,col = found
    print(f"=> chi_c(G_GW) <= {p}/{q} = {r:.5f} < 14  ==>  LINEARIZATION FAILS AT GW.")
    print(f"   The circular chromatic rung is BLIND to GW's tightness (like chi_f=13);")
    print(f"   only the INTEGER rung chi=14 sees it (opus THM-652).  Non-rotation witness:")
    print(f"   C(x) = c(x mod {T}) + (x div {T})*{Delta} mod {p},  c = {col}")
else:
    print("=> NO quasi-periodic sub-14 coloring found (T<=26): strong evidence chi_c(G_GW)=14")
    print("   (rotation optimal among quasi-periodic; linearization HOLDS at GW).")
    print("   Next: general (aperiodic) SAT with large budget, or an obstruction proof.")
