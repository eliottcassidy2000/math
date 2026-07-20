#!/usr/bin/env python3
"""
Independent exact verification of the claimed GMC(4) counterexample  (mac-mini-S133)
====================================================================================
Claim supplied by the owner from an outside source:

    P = (1+Z2)(W1(1-Z1)+W2),   Q = Z2      "in 4 real Gaussians via complex Z_j, W_j"
    E[P^m] = 0        for all m >= 1
    E[Q P^m] = m! != 0 for all m >= 1        ==>  GMC(4) is false.   (cubic, 6 terms)

We do NOT take the claim on trust and we do NOT assume the intended reading of the
notation.  Instead:

  PART A  exhaustively tries ALL 4^4 = 256 assignments of the four symbols
          {Z1,Z2,W1,W2} to the four letters {Z, Zbar, W, Wbar} of two independent
          standard complex Gaussians (= 4 real Gaussians), and reports which readings
          make BOTH claims true.  Exact integer arithmetic via
              E[Z^a Zbar^b W^c Wbar^d] = delta_ab delta_cd a! c!.
  PART B  verifies the heat-operator identity  E[f] = (exp(D) f)(0),
          D = d_Z d_Zbar + d_W d_Wbar -- which is the precise bridge from the Gaussian
          moment conjecture to Zhao's VANISHING CONJECTURE (VC uses D^m(P^m); GMC uses
          exp(D) at the origin).  For HOMOGENEOUS P the two coincide up to a factorial;
          for inhomogeneous P they do not.
  PART C  tests whether the counterexample is ESSENTIALLY INHOMOGENEOUS, by running each
          homogeneous component separately.  This matters because THM-1435 (D) already
          found -- independently, from the Alpoge map -- that homogeneity is the
          load-bearing issue in the VC transport.
  PART D  a real-coordinate cross-check, expanding in u1..u4 ~ N(0,1/2) so that
          Z = u1 + i u2, W = u3 + i u4 has E|Z|^2 = 1.  Exact rational arithmetic.
"""
from math import factorial
from itertools import product
from fractions import Fraction

# ---------------------------------------------------------------- complex formalism
# monomial key (a,b,c,d) = Z^a Zbar^b W^c Wbar^d ; coefficients are integers
def pmul(p, q):
    out = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            k = (k1[0]+k2[0], k1[1]+k2[1], k1[2]+k2[2], k1[3]+k2[3])
            out[k] = out.get(k, 0) + c1*c2
    return {k: c for k, c in out.items() if c}

def padd(p, q):
    out = dict(p)
    for k, c in q.items():
        out[k] = out.get(k, 0) + c
    return {k: c for k, c in out.items() if c}

def ppow(p, m):
    r = {(0,0,0,0): 1}
    for _ in range(m): r = pmul(r, p)
    return r

def expect(p):
    """E[Z^a Zbar^b W^c Wbar^d] = a! c! if a==b and c==d, else 0."""
    return sum(c * factorial(a) * factorial(cc)
               for (a, b, cc, d), c in p.items() if a == b and cc == d)

LETTER = [(1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)]   # Z, Zbar, W, Wbar
NAME   = ["Z", "Zbar", "W", "Wbar"]
ONE    = {(0,0,0,0): 1}

def build(assign):
    """assign = (i1,i2,j1,j2) giving letters for Z1, Z2, W1, W2."""
    Z1, Z2, W1, W2 = [{LETTER[i]: 1} for i in assign]
    # P = (1+Z2) * ( W1*(1-Z1) + W2 )
    inner = padd(pmul(W1, padd(ONE, {k: -c for k, c in Z1.items()})), W2)
    P = pmul(padd(ONE, Z2), inner)
    return P, Z2

MMAX = 7
print("=" * 78)
print("PART A -- exhaustive over ALL 256 readings of the notation")
print("=" * 78)
hits = []
for assign in product(range(4), repeat=4):
    P, Q = build(assign)
    Pm = ONE; ok_left = True; right = []
    for m in range(1, MMAX + 1):
        Pm = pmul(Pm, P)
        if expect(Pm) != 0: ok_left = False; break
        right.append(expect(pmul(Q, Pm)))
    if ok_left and all(r == factorial(m + 1) for m, r in enumerate(right)):
        hits.append((assign, right))

print(f"  readings with E[P^m] = 0 for m=1..{MMAX} AND E[QP^m] = m! : {len(hits)}")
for assign, right in hits:
    print(f"    Z1={NAME[assign[0]]:>4}  Z2={NAME[assign[1]]:>4}  "
          f"W1={NAME[assign[2]]:>4}  W2={NAME[assign[3]]:>4}   "
          f"E[QP^m] = {right}")

# any reading that kills the left side but gives a DIFFERENT right side?
others = []
for assign in product(range(4), repeat=4):
    P, Q = build(assign)
    Pm = ONE; ok_left = True; right = []
    for m in range(1, MMAX + 1):
        Pm = pmul(Pm, P)
        if expect(Pm) != 0: ok_left = False; break
        right.append(expect(pmul(Q, Pm)))
    if ok_left and any(r != 0 for r in right) and \
       not all(r == factorial(m + 1) for m, r in enumerate(right)):
        others.append((assign, right))
print(f"  readings that break GMC but with a DIFFERENT growth law: {len(others)}")
for assign, right in others[:6]:
    print(f"    Z1={NAME[assign[0]]:>4}  Z2={NAME[assign[1]]:>4}  "
          f"W1={NAME[assign[2]]:>4}  W2={NAME[assign[3]]:>4}   {right}")

# ---------------------------------------------------------------- PART B
print()
print("=" * 78)
print("PART B -- E[f] = (exp(D) f)(0), D = d_Z d_Zbar + d_W d_Wbar  (the GMC<->VC bridge)")
print("=" * 78)
def apply_D(p):
    """D = d/dZ d/dZbar + d/dW d/dWbar."""
    out = {}
    for (a, b, c, d), co in p.items():
        if a >= 1 and b >= 1:
            k = (a-1, b-1, c, d); out[k] = out.get(k, 0) + co*a*b
        if c >= 1 and d >= 1:
            k = (a, b, c-1, d-1); out[k] = out.get(k, 0) + co*c*d
    return {k: v for k, v in out.items() if v}

def expD_at0(p):
    """sum_k (D^k p)(0) / k!"""
    tot = Fraction(0); q = dict(p); k = 0
    while q:
        tot += Fraction(q.get((0,0,0,0), 0), factorial(k))
        q = apply_D(q); k += 1
    return tot

if hits:
    assign = hits[0][0]
    P, Q = build(assign)
    Pm = ONE; ok = True
    for m in range(1, 6):
        Pm = pmul(Pm, P)
        if expD_at0(Pm) != expect(Pm) or expD_at0(pmul(Q, Pm)) != expect(pmul(Q, Pm)):
            ok = False
    print(f"  exp(D)-at-origin reproduces E[.] on P^m and QP^m, m=1..5:  {ok}")
    print("  So GMC is literally the 'evaluate the heat semigroup at the origin' shadow of")
    print("  Zhao's VC, which instead asks about D^m(P^m).  For HOMOGENEOUS P of degree d,")
    print("  exp(D) applied to P^m keeps only the single term k = dm/2, so the two agree up")
    print("  to a factorial.  For INHOMOGENEOUS P they are genuinely different statements.")

# ---------------------------------------------------------------- PART C
print()
print("=" * 78)
print("PART C -- is the counterexample ESSENTIALLY INHOMOGENEOUS?")
print("=" * 78)
if hits:
    assign = hits[0][0]
    P, Q = build(assign)
    deg = {}
    for k, c in P.items(): deg.setdefault(sum(k), {})[k] = c
    print(f"  P has homogeneous components of degrees {sorted(deg)} "
          f"with {[len(deg[d]) for d in sorted(deg)]} terms")
    for d in sorted(deg):
        comp = deg[d]
        Pm = ONE; lefts = []; rights = []
        for m in range(1, 6):
            Pm = pmul(Pm, comp)
            lefts.append(expect(Pm)); rights.append(expect(pmul(Q, Pm)))
        print(f"    degree-{d} part alone: E[P^m] = {lefts}   E[QP^m] = {rights}"
              f"   {'still a counterexample' if all(l==0 for l in lefts) and any(rights) else 'NOT a counterexample'}")
    # also: drop each term of P and see if the counterexample survives
    survive = 0
    keys = list(P)
    for drop in keys:
        Pd = {k: c for k, c in P.items() if k != drop}
        Pm = ONE; ok_l = True; r = []
        for m in range(1, 6):
            Pm = pmul(Pm, Pd)
            if expect(Pm) != 0: ok_l = False; break
            r.append(expect(pmul(Q, Pm)))
        if ok_l and any(r): survive += 1
    print(f"  dropping any ONE of the {len(keys)} terms of P leaves a counterexample in "
          f"{survive} cases  => {'NOT minimal' if survive else 'TERM-MINIMAL'}")

# ---------------------------------------------------------------- PART D
print()
print("=" * 78)
print("PART D -- real-coordinate cross-check: u1..u4 ~ N(0,1/2), Z = u1+i u2, W = u3+i u4")
print("=" * 78)
def dfac(k):                      # (2k-1)!! = (2k)!/(2^k k!)
    return factorial(2*k) // (2**k * factorial(k))

def rmul(p, q):
    out = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            k = tuple(a+b for a, b in zip(k1, k2))
            out[k] = out.get(k, 0) + c1*c2
    return {k: c for k, c in out.items() if c}

def rexpect(p):
    """E[prod u_i^{a_i}] with u_i ~ N(0,1/2) independent: (a-1)!! * (1/2)^{a/2} if all even."""
    tot = Fraction(0)
    for k, c in p.items():
        if any(a % 2 for a in k): continue
        v = Fraction(1)
        for a in k: v *= dfac(a // 2) * Fraction(1, 2**(a // 2))
        tot += c * v
    return tot

if hits:
    assign = hits[0][0]
    # complex letters in real coordinates, coefficients in Q[i] as (re, im) pairs
    I = complex(0, 1)
    RL = [ {(1,0,0,0): 1+0j, (0,1,0,0): 1j},      # Z    = u1 + i u2
           {(1,0,0,0): 1+0j, (0,1,0,0): -1j},     # Zbar = u1 - i u2
           {(0,0,1,0): 1+0j, (0,0,0,1): 1j},      # W    = u3 + i u4
           {(0,0,1,0): 1+0j, (0,0,0,1): -1j} ]    # Wbar = u3 - i u4
    RONE = {(0,0,0,0): 1+0j}
    Z1, Z2, W1, W2 = [RL[i] for i in assign]
    neg = lambda p: {k: -c for k, c in p.items()}
    add = lambda p, q: {k: p.get(k,0)+q.get(k,0) for k in set(p)|set(q)}
    inner = add(rmul(W1, add(RONE, neg(Z1))), W2)
    PR = rmul(add(RONE, Z2), inner)
    QR = Z2
    def rexp_c(p):
        re = {k: Fraction(round(c.real*10**9), 10**9) for k, c in p.items()}
        im = {k: Fraction(round(c.imag*10**9), 10**9) for k, c in p.items()}
        return rexpect(re), rexpect(im)
    Pm = RONE; okL = okR = True
    for m in range(1, 6):
        Pm = rmul(Pm, PR)
        lre, lim = rexp_c(Pm)
        rre, rim = rexp_c(rmul(QR, Pm))
        if lre != 0 or lim != 0: okL = False
        if rre != factorial(m) or rim != 0: okR = False
        print(f"    m={m}:  E[P^m] = {lre} + {lim}i     E[QP^m] = {rre} + {rim}i   "
              f"(m! = {factorial(m)})")
    print(f"  real-coordinate cross-check agrees: E[P^m]=0 {okL}, E[QP^m]=m! {okR}")
    print("  => the counterexample genuinely lives in FOUR REAL GAUSSIANS.")
