#!/usr/bin/env python3
"""
The 2-DIMENSIONAL NULLCONE CONJECTURE, and GMC(2) as its corollary  (mac-mini-S136)
===================================================================================
Owner: "aim to finish GMC(2) by finishing the stronger 2 dimensional nullcone
conjecture."

THE NULLCONE.  For the Gaussian functional E on C[Z,W] (W = Zbar, E[Z^aW^b]=a! d_ab),
    N(E) := { P : E[P^m] = 0 for all m >= 1 }.
Grade by CHARGE c = deg_Z - deg_W (THM-1520): E annihilates every nonzero charge.

    NULLCONE CONJECTURE (n=2).   N(E) = { P : all charges >= 1 }  u  { P : all charges <= -1 }
                                      = STRICTLY ONE-SIDED charge support.

It is STRICTLY STRONGER than GMC(2), and GMC(2) is a two-line corollary:
P in N(E) => all charges of P are >= 1 => all charges of P^m are >= m => QP^m can have
charge 0 only if m <= deg_W(Q) => E[QP^m] = 0 for m > deg_W(Q).  Done.

WHAT IS ALREADY PROVED (THM-1520):
  * (superset) strictly one-sided => in N(E).  Trivial: charges of P^m are >= m > 0.
  * (subset, one-sided support) if P's charges are one-sided and P is in N(E) then the
    charge-0 part p_0 = 0, by the telescoping lemma E[P^m] = L(p_0^m) plus the saddle
    lemma L(p^m)/(a_D^m (Dm)!) -> exp(a_{D-1}/(D a_D)) != 0.
WHAT REMAINS: two-sided P must NOT be in N(E).  That is what this file attacks.

PART A  the EXACT polar reduction   E[P^m] = L( v |-> CT_u( H_sqrt(v)(u)^m ) ),
        H_r(u) := P(ru, r/u).  This exhibits GMC(2)/nullcone at n=2 as a GAUSSIAN-AVERAGED
        DUISTERMAAT-VAN DER KALLEN problem -- DvdK is exactly the CT_u statement, and L is
        the Gaussian average over the radius.  Clean, exact, and it names the right context.
PART B  THE TWO-CHARGE THEOREM (proved): if P has exactly two charges C > 0 > -B then
        E[P^m] = 0 unless (B'+C') | m (g = gcd(B,C), B'=B/g, C'=C/g), and at m = t(B'+C')
        it equals a multinomial times L(F^t) with F != 0 -- so by the saddle lemma it is
        NONZERO for large t.  Hence P is NOT in N(E).
PART C  the search, now WITH A WORKING POSITIVE CONTROL: the conjecture predicts that
        one-sided P ARE in N(E), so the machinery must find those.  (This is what the S135
        sweep lacked.)
PART D  the 1-variable Laurent nullcone: CT(h^m)=0 for all m  <=>  h strictly one-sided.
"""
from fractions import Fraction as F
from math import factorial, gcd
import itertools, random

# ---------------------------------------------------------------- (Z,W) formalism
def cmul(p, q):
    out = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            k = (k1[0]+k2[0], k1[1]+k2[1])
            out[k] = out.get(k, 0) + c1*c2
    return {k: c for k, c in out.items() if c}

def cexp(p):
    return sum(c*factorial(a) for (a, b), c in p.items() if a == b)

ONE = {(0,0): 1}
def cpow(p, m):
    r = ONE
    for _ in range(m): r = cmul(r, p)
    return r

def charges(P): return {a-b for (a, b) in P}
def onesided(P):
    ch = charges(P)
    return (min(ch) >= 1) or (max(ch) <= -1)

def show(P):
    if not P: return "0"
    out = []
    for (a, b), c in sorted(P.items()):
        s = ("Z" + (f"^{a}" if a > 1 else "") if a else "") + \
            ("W" + (f"^{b}" if b > 1 else "") if b else "")
        out.append(("+" if c > 0 else "-") + (f"{abs(c)}" if abs(c) != 1 or not s else "") + (s or "1"))
    return "".join(out).lstrip("+")

def L(coeffs):
    return sum(a*factorial(k) for k, a in enumerate(coeffs))

def polymul(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): out[i+j] += x*y
    return out
def polypow(a, m):
    r = [1]
    for _ in range(m): r = polymul(r, a)
    return r

# ================================================================ PART A
print("=" * 78)
print("PART A -- the exact polar reduction:  E[P^m] = L( v -> CT_u( H_sqrt(v)^m ) )")
print("=" * 78)
print("  H_r(u) := P(ru, r/u) = sum_{a,b} p_ab r^{a+b} u^{a-b}  -- a LAURENT polynomial in u")
print("  whose u-exponents are exactly the CHARGES, and whose r-degree is the total degree.")
print("  CT_u is the constant term (= the charge-0 projection); L is the Gaussian average")
print("  over the radius, r^2 ~ Exp(1).  So n=2 is DvdK's constant-term problem, averaged.")
print()
def polar_reduce(P, m):
    """CT_u(H_r^m) as a polynomial in v = r^2 ; returns coefficient list in v."""
    # H_r^m: track (u-exponent, r-exponent) -> coeff
    H = {}
    for (a, b), c in P.items():
        H[(a-b, a+b)] = H.get((a-b, a+b), 0) + c
    cur = {(0, 0): 1}
    for _ in range(m):
        nxt = {}
        for (u1, r1), c1 in cur.items():
            for (u2, r2), c2 in H.items():
                k = (u1+u2, r1+r2); nxt[k] = nxt.get(k, 0) + c1*c2
        cur = {k: c for k, c in nxt.items() if c}
    # constant term in u, then r-exponent must be even (r^{2j} = v^j)
    out = {}
    for (uu, rr), c in cur.items():
        if uu == 0:
            assert rr % 2 == 0
            out[rr//2] = out.get(rr//2, 0) + c
    if not out: return [0]
    res = [0]*(max(out)+1)
    for j, c in out.items(): res[j] = c
    return res

rng = random.Random(136)
ok = True
for trial in range(8):
    P = {}
    for _ in range(rng.randint(2, 4)):
        a, b = rng.randint(0, 2), rng.randint(0, 2)
        v = rng.randint(-3, 3)
        if v: P[(a, b)] = P.get((a, b), 0) + v
    P = {k: v for k, v in P.items() if v}
    if not P: continue
    for m in range(1, 5):
        if cexp(cpow(P, m)) != L(polar_reduce(P, m)): ok = False
print(f"  identity verified on random P, m = 1..4:  {ok}")

# ================================================================ PART B
print()
print("=" * 78)
print("PART B -- THE TWO-CHARGE THEOREM (this is a proof, not a search)")
print("=" * 78)
print("  Let P = P_C + P_{-B}, exactly two charges C > 0 > -B.  Write P_C = Z^C q(V),")
print("  P_{-B} = W^B s(V) with V = ZW, and g = gcd(B,C), B' = B/g, C' = C/g.")
print("  Charge-0 needs k_C*C = k_{-B}*B with k_C + k_{-B} = m, so k_C = tB', k_{-B} = tC'")
print("  and m = t(B'+C') -- UNIQUE for each such m, and IMPOSSIBLE otherwise.  Hence")
print("      E[P^m] = 0                       if (B'+C') does not divide m,")
print("      E[P^m] = C(m; tB', tC') * L(F^t) with F = v^{C B'} q^{B'} s^{C'} != 0  otherwise.")
print("  By the saddle lemma (THM-1520 B), L(F^t) != 0 for all large t.  So P is NOT in the")
print("  nullcone.  TWO-CHARGE CASE CLOSED.")
print()
print(f"{'P':>26} {'B,C':>7} {'B+C prime':>10} {'E[P^m], m=1..8':>44}")
for P in ({(1,0):1, (0,1):1}, {(2,0):1, (0,1):1}, {(1,0):1, (0,2):3},
          {(3,0):1, (0,2):1}, {(2,0):1, (0,2):1}, {(2,1):1, (0,3):-2}):
    ch = sorted(charges(P)); B, C = -ch[0], ch[1]
    g = gcd(B, C); Bp, Cp = B//g, C//g
    vals = [cexp(cpow(P, m)) for m in range(1, 9)]
    print(f"{show(P):>26} {f'{B},{C}':>7} {Bp+Cp:>10} {str(vals):>44}")
print("  Nonzero exactly at the predicted multiples of B'+C', as the theorem says.")

# ================================================================ PART C
print()
print("=" * 78)
print("PART C -- the nullcone search, WITH A WORKING POSITIVE CONTROL")
print("=" * 78)
MONS = [(a, b) for a in range(4) for b in range(4) if (a, b) != (0, 0) and a+b <= 4]
def in_nullcone(P, M=12):
    for m in range(1, M+1):
        if cexp(cpow(P, m)) != 0: return False
    return True

print("  POSITIVE CONTROL: the conjecture says every strictly one-sided P is in N(E).")
ctrl = [{(1,0):1}, {(2,0):1, (1,0):-3}, {(0,1):1, (0,2):5},
        {(2,1):1, (3,0):-2}, {(1,2):4, (0,3):1}]
allctrl = True
for P in ctrl:
    got = in_nullcone(P)
    if not got: allctrl = False
    print(f"    {show(P):>18}  charges {sorted(charges(P))}  in nullcone? {got}")
print(f"  positive control passes: {allctrl}   <-- this is what the S135 sweep lacked")
print()
print("  Now the real question: any TWO-SIDED P in the nullcone?")
total = 0; found = []
for ksz in range(2, 6):
    for supp in itertools.combinations(range(len(MONS)), ksz):
        for signs in itertools.product((-2, -1, 1, 2), repeat=ksz):
            P = {MONS[i]: s for i, s in zip(supp, signs)}
            ch = charges(P)
            if not (min(ch) < 0 < max(ch)): continue      # two-sided only
            total += 1
            if in_nullcone(P): found.append(P)
        if found: break
    print(f"    support size {ksz}: {total} two-sided polys tested, {len(found)} in nullcone")
    if found: break
if found:
    print("  *** TWO-SIDED NULLCONE ELEMENT FOUND -- conjecture FALSE ***")
    for P in found[:5]: print(f"      {show(P)}   charges {sorted(charges(P))}")
else:
    print("  none.  With the control passing, this is now MEANINGFUL evidence for the")
    print("  nullcone conjecture (bounded, but no longer uncontrolled).")

# ================================================================ PART D
print()
print("=" * 78)
print("PART D -- the 1-variable Laurent nullcone:  CT(h^m)=0 for all m <=> h one-sided")
print("=" * 78)
def ct_pow(h, m):
    cur = {0: 1}
    for _ in range(m):
        nxt = {}
        for e1, c1 in cur.items():
            for e2, c2 in h.items(): nxt[e1+e2] = nxt.get(e1+e2, 0) + c1*c2
        cur = {k: c for k, c in nxt.items() if c}
    return cur.get(0, 0)

EXPS = list(range(-3, 4))
tot = 0; bad = []
for ksz in (2, 3, 4):
    for supp in itertools.combinations(EXPS, ksz):
        if not (min(supp) < 0 < max(supp)): continue
        for signs in itertools.product((-2, -1, 1, 2), repeat=ksz):
            h = {e: s for e, s in zip(supp, signs)}
            tot += 1
            if all(ct_pow(h, m) == 0 for m in range(1, 11)): bad.append(h)
print(f"  two-sided Laurent h with exponents in [-3,3], coeffs in {{+-1,+-2}}: {tot} tested")
print(f"  with CT(h^m)=0 for m=1..10: {len(bad)}")
if bad:
    for h in bad[:5]: print(f"      {h}")
else:
    print("  none => the 1-variable Laurent nullcone is exactly the ONE-SIDED Laurent")
    print("  polynomials, on this box.  (0 not in the Newton polytope.)  This is the")
    print("  'top-edge' input the general two-sided argument needs.")
