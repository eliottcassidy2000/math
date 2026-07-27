#!/usr/bin/env python3
# -*- coding: ascii -*-
"""
msg2153_trichotomy_hostile_probe_opus.py

Adversarial refutation probe of MSG-2153's two-colour trichotomy
(composition of THM-2323 Sections 6-7 with THM-2326, at gcd colour 91).

The trichotomy under attack (7|a branch, a=c2/g, g=gcd(c2,c3)):
  take an adjacent coherent-fibre pair A=c2*q_z0, B=c2*q_z1 from the
  THM-2323 Section 7.2 fibre (whose construction does NOT use hypothesis
  (31); (31) is used only for the final gcd(m,91)=1 of (39)), so
  B-A=m2*c3 with 7|m2 and 13 not | m2; apply THM-2326 at A to get
  C=A+m1*c3 with 7 not | m1; if 13 not | m1 the edge A-C is 91-unit,
  else the pivot edge C-B=(m1-m2)*c3 is 91-unit.

Attack modes implemented (all exact integer / cyclotomic-integer
arithmetic; floats only in one convergence sanity check):

  PART A  exhaustive residue stress of the (m1,m2) split: search for any
          pair with 7 not|m1, 7|m2, 13 not|m2 where NEITHER A-C NOR C-B
          is 91-unit, or where C=B (m1=m2), plus grade/root-character
          retention sweeps for all three vertices.
  PART B  coherent-fibre degeneracy stress: K_z distinctness mod N,
          retained-fibre sizes (K_7xK_13 vs K_6xK_13), edge existence,
          bracket units for every edge and gauge offset, over several
          (a,d') configurations including 7|a.
  PART C  the MSG-2152 sharp boundary hostile, EXTENDED adversarially to
          carry a (7|m2, 13-unit) pair atom: real step H vanishing on a
          1/7 arc, spectrum in 7Z u 13Z, Hhat(91)!=0, every incident
          multiplier m with 7 not|m divisible by 13 (blocks
          disjointness+Prony alone), PLUS Hhat(98)!=0 (pair atom, m2=7).
          Exact cyclotomic verification that the extension yields the
          pivot unit edge, i.e. CANNOT kill the composed trichotomy.
  PART D  concrete 7|a canonical-shaped instance (c2,c3)=(91,169):
          a=7, d'=13, b=1, c=2; full symbolic composition with integer
          witnesses and retention checks.

Every load-bearing check raises AssertionError on failure.
"""

import cmath
from fractions import Fraction
from functools import reduce
from math import gcd, pi, sin


def check(label, cond):
    print(("[PASS] " if cond else "[FAIL] ") + label)
    if not cond:
        raise AssertionError(label)


def nu13(n):
    n = abs(n)
    assert n != 0
    k = 0
    while n % 13 == 0:
        n //= 13
        k += 1
    return k


def lcm2(x, y):
    return x * y // gcd(x, y)


# ---------------------------------------------------------------------------
# exact cyclotomic machinery: zero-testing sums of roots of unity in Z[zeta_L]
# ---------------------------------------------------------------------------

def divisors(n):
    ds = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            ds.append(i)
            if i != n // i:
                ds.append(n // i)
        i += 1
    return sorted(ds)


def poly_divmod(a, b):
    # b monic; integer coefficients; ascending order
    a = list(a)
    db = len(b) - 1
    assert b[-1] == 1
    q = [0] * max(len(a) - db, 1)
    for i in range(len(a) - 1, db - 1, -1):
        coef = a[i]
        if coef:
            q[i - db] = coef
            for j in range(db + 1):
                a[i - db + j] -= coef * b[j]
    return q, a[:db]


_cyc = {}


def cyclotomic(n):
    if n in _cyc:
        return _cyc[n]
    p = [-1] + [0] * (n - 1) + [1]  # X^n - 1
    for d in divisors(n):
        if d < n:
            q, r = poly_divmod(p, cyclotomic(d))
            assert all(c == 0 for c in r), "cyclotomic division not exact"
            while len(q) > 1 and q[-1] == 0:
                q.pop()
            p = q
    _cyc[n] = p
    return p


def build_pow_table(L):
    """pow_table[t] = X^t mod Phi_L(X) as an integer vector of length D."""
    Phi = cyclotomic(L)
    D = len(Phi) - 1
    tbl = []
    cur = [0] * D
    cur[0] = 1
    for _ in range(L):
        tbl.append(cur)
        nf = [0] + cur
        top = nf[D]
        if top:
            nf = [nf[j] - top * Phi[j] for j in range(D + 1)]
        cur = nf[:D]
    return tbl, D


# ---------------------------------------------------------------------------
# exact rational step functions on the circle
# ---------------------------------------------------------------------------

def step_value(pieces, x):
    v = pieces[-1][1]
    for s, val in pieces:
        if s <= x:
            v = val
        else:
            break
    return v


def step_jumps(pieces):
    js = []
    n = len(pieces)
    for i in range(n):
        d = pieces[i][1] - pieces[i - 1][1]
        if d != 0:
            js.append((pieces[i][0], d))
    return js


def hat_nonzero_exact(jumps, n, L, tbl, D):
    """Exact test (1/(2*pi*i*n)) * sum Delta_j zeta^{-n t_j} != 0, n != 0."""
    assert n != 0
    accd = {}
    for pos, d in jumps:
        t = pos * L
        assert t.denominator == 1
        idx = (-n * int(t)) % L
        accd[idx] = accd.get(idx, 0) + d
    acc = [0] * D
    for t, cval in accd.items():
        if cval:
            row = tbl[t]
            for j in range(D):
                acc[j] += cval * row[j]
    return any(v != 0 for v in acc)


def hat_float(jumps, n):
    s = sum(d * cmath.exp(-2j * pi * n * float(pos)) for pos, d in jumps)
    return s / (2j * pi * n)


# ---------------------------------------------------------------------------
print("=" * 72)
print("PART A -- exhaustive residue stress of the trichotomy (m1,m2) split")
print("=" * 72)

M1MAX = 1400
M2MAX = 1400
pairs = 0
m1_eq_m2 = 0
no_unit_edge = 0
for m1 in range(1, M1MAX + 1):
    if m1 % 7 == 0:
        continue
    unit1 = (gcd(m1, 91) == 1)
    # case coverage: 7 not | m1  =>  (13 not|m1 <=> gcd(m1,91)=1)
    assert unit1 == (m1 % 13 != 0)
    for m2 in range(-M2MAX, M2MAX + 1):
        if m2 == 0 or m2 % 7 != 0 or m2 % 13 == 0:
            continue
        pairs += 1
        if m1 == m2:
            m1_eq_m2 += 1
        # C = B  <=>  m1 = m2, impossible since m1-m2 != 0 mod 7
        assert (m1 - m2) % 7 != 0
        # trichotomy: does at least one of A-C, C-B carry a 91-unit?
        pivot_unit = (gcd(m1 - m2, 91) == 1)
        if m1 % 13 == 0:
            # case 13|m1: the pivot MUST be unit
            assert pivot_unit, (m1, m2)
        if not (unit1 or pivot_unit):
            no_unit_edge += 1
check("A1: swept %d (m1,m2) pairs, 7 not|m1, 7|m2, 13 not|m2" % pairs,
      pairs > 300000)
check("A1: m1 = m2 (degenerate C=B) occurrences: %d" % m1_eq_m2,
      m1_eq_m2 == 0)
check("A1: pairs where NEITHER A-C nor C-B is a 91-unit edge: %d"
      % no_unit_edge, no_unit_edge == 0)

# A2: grade and root-character retention for A, B=A+m2*c3, C=A+m1*c3
ret_rows = 0
for b in (0, 1, 2):
    for c in (b + 1, b + 2):
        for u2 in (1, 2, 7, 14, 49):        # 13 not | u2
            for u3 in (1, 2, 7, 11):        # 13 not | u3
                c3 = 13 ** c * u3
                for q in (1, 2, 3, 5, 17, 275, 3576):
                    if q % 13 == 0:
                        continue
                    A = (13 ** b * u2) * q   # grade-b atom, char u2*q mod 13
                    charA = (A // 13 ** b) % 13
                    for m in list(range(1, 60)) + [-13, -26, 91, 1925, -1925]:
                        if m == 0:
                            continue
                        V = A + m * c3
                        assert V != 0
                        assert nu13(V) == b
                        assert (V // 13 ** b) % 13 == charA
                        ret_rows += 1
check("A2: %d grade/root retention rows: nu13 preserved (=b) and "
      "root character preserved for every c3-translate" % ret_rows,
      ret_rows > 30000)
print("A verdict: the residue split is airtight -- m1=m2 is impossible mod 7,")
print("and when 13|m1 the pivot difference m1-m2 is automatically a 91-unit;")
print("no translate can move grade or root character since nu13(c3)=c>b.")

# ---------------------------------------------------------------------------
print()
print("=" * 72)
print("PART B -- coherent-fibre degeneracy stress (THM-2323 S7.2 fibre)")
print("=" * 72)


def fibre_probe(a, dprime, kappa):
    assert gcd(a, dprime) == 1 and dprime % 13 == 0
    if a % 7 == 0:
        assert dprime % 7 != 0  # coprimality: 7|a and 7|d' impossible
    N = 91 * dprime
    K0 = None
    for k in range(1, dprime):
        if gcd(k, dprime) == 1 and k % 13 == kappa:
            K0 = k
            break
    assert K0 is not None, "CRT representative exists"
    Ks = [K0 + dprime * z for z in range(91)]
    # distinctness mod N: kills any q_z0 = q_z1 degeneracy at the source
    assert len(set(k % N for k in Ks)) == 91
    retained = [z for z in range(91) if gcd(Ks[z], N) == 1]
    expect = 91 if dprime % 7 == 0 else 78
    assert len(retained) == expect, (a, dprime, len(retained))
    # all retained K_z share the prescribed residue mod 13
    assert all(Ks[z] % 13 == kappa for z in retained)
    edges = [(z1, z2) for z1 in retained for z2 in retained
             if z1 < z2 and (z2 - z1) % 7 != 0 and (z2 - z1) % 13 != 0]
    assert len(edges) > 0, "fibre graph edgeless"
    checked = 0
    for (z1, z2) in edges:
        for dh in range(-6, 7):
            br = (z2 - z1) + 91 * dh
            assert gcd(br, 91) == 1          # bracket is a 91-unit
            m2 = a * br
            assert m2 % 13 != 0              # 13 never divides m2
            if a % 7 == 0:
                assert m2 % 7 == 0           # 7|a forces 7|m2
            checked += 1
    return len(retained), len(edges), checked


configs = [(7, 13, 1), (7, 143, 2), (14, 13, 5), (7, 2431, 4),
           (1, 91, 3), (2, 91, 1)]
for (a, dp, kp) in configs:
    nret, nedge, nchk = fibre_probe(a, dp, kp)
    graph = "K_7xK_13" if dp % 7 == 0 else "K_6xK_13"
    print("  a=%-3d d'=%-5d %s: retained=%d edges=%d bracket-checks=%d"
          % (a, dp, graph, nret, nedge, nchk))
check("B1: all fibre configs: K_z distinct mod N, retained fibre full "
      "(78 or 91), graph has edges, every edge+offset bracket is a 91-unit",
      True)
print("B verdict: the adjacent pair can never be degenerate (adjacency needs")
print("both CRT coordinates to differ, so z0!=z1, and K_z are distinct mod N,")
print("so q_z0 != q_z1 and A != B); the fibre construction never uses (31),")
print("so it survives 7|a with multiplier m2 = a*(91-unit bracket).")

# ---------------------------------------------------------------------------
print()
print("=" * 72)
print("PART C -- MSG-2152 boundary hostile, extended with a pair atom")
print("=" * 72)

# H = 1_F - per7(1_F restricted to [0,1/7)), F = union_r [r/13,(r+1/4)/13)
w = Fraction(1, 4)
piecesF = []
for r in range(13):
    piecesF.append((Fraction(r, 13), 1))
    piecesF.append((Fraction(r, 13) + w / 13, 0))

arc = Fraction(1, 7)
sub = [(s, v) for (s, v) in piecesF if s < arc]
assert sub[0][0] == 0
piecesP = []
for k in range(7):
    for (s, v) in sub:
        piecesP.append((Fraction(k, 7) + s, v))
piecesP.sort()

starts = sorted(set(s for s, _ in piecesF) | set(s for s, _ in piecesP))
piecesH = [(s, step_value(piecesF, s) - step_value(piecesP, s))
           for s in starts]

# C1: H is a real integer-valued step function vanishing on the arc [0,1/7)
check("C1: H real integer-valued step function",
      all(isinstance(v, int) for _, v in piecesH))
check("C1: H vanishes identically on the 1/7 arc [0,1/7)",
      all(v == 0 for s, v in piecesH if s < arc))

jumpsH = step_jumps(piecesH)
L = reduce(lcm2, [s.denominator for s, _ in jumpsH], 1)
check("C1: jump positions live on the exact 1/%d lattice, %d jumps"
      % (L, len(jumpsH)), L == 364)
tbl, D = build_pow_table(L)
check("C1: cyclotomic reduction ready, deg Phi_%d = %d" % (L, D), D == 144)

# C2: spectrum contained in 7Z u 13Z (exact, range check)
bad = [n for n in range(1, 501)
       if n % 7 != 0 and n % 13 != 0
       and hat_nonzero_exact(jumpsH, n, L, tbl, D)]
check("C2: Hhat(n)=0 exactly for every n in [1,500] with 7 not|n, 13 not|n "
      "(spectrum in 7Z u 13Z)", not bad)

# C3: marked atom and adversarial pair atom
check("C3: Hhat(91) != 0 (marked atom A=91, exact)",
      hat_nonzero_exact(jumpsH, 91, L, tbl, D))
check("C3: Hhat(98) != 0 (pair atom B=98=A+7, m2=7: 7|m2, 13 not|m2, exact)",
      hat_nonzero_exact(jumpsH, 98, L, tbl, D))

# C4: the hostile blocking property (this is what kills disjointness+Prony)
incident = [m for m in range(-260, 261)
            if m != 0 and m % 7 != 0
            and (91 + m) != 0
            and hat_nonzero_exact(jumpsH, 91 + m, L, tbl, D)]
check("C4: every incident multiplier m with 7 not|m in [-260,260] is "
      "divisible by 13 (found %d, e.g. %s...)"
      % (len(incident), incident[:6]),
      len(incident) > 0 and all(m % 13 == 0 for m in incident)),

# C5: THM-2326 landing exists despite the hostile (identity forces some m1)
pos_m1 = [m for m in incident if m > 0]
check("C5: a positive m1 with 7 not|m1 and Hhat(91+m1)!=0 exists "
      "(smallest: %s)" % (min(pos_m1) if pos_m1 else None),
      len(pos_m1) > 0)

# C6: the trichotomy pivot on the extended hostile: for EVERY admissible m1
#     the edge C-B is a 91-unit, so the extension cannot kill the trichotomy
m2 = 7
pivot_fail = [m1 for m1 in incident if gcd(m1 - m2, 91) != 1]
check("C6: for EVERY incident m1 (all 13-divisible), the pivot edge "
      "C-B = (m1-7)*c has gcd(m1-7,91)=1 -- failures: %d" % len(pivot_fail),
      not pivot_fail)
m1w = min(pos_m1)
Cw = 91 + m1w
check("C6: explicit unit edge inside the hostile itself: C=%d, B=98, "
      "C-B=%d, gcd(%d,91)=1, Hhat(C)!=0, Hhat(B)!=0"
      % (Cw, Cw - 98, Cw - 98),
      gcd(Cw - 98, 91) == 1
      and hat_nonzero_exact(jumpsH, Cw, L, tbl, D)
      and hat_nonzero_exact(jumpsH, 98, L, tbl, D))

# C7: float sanity of the modulated disjointness identity (THM-2326 (9)),
#     applied to H shifted by 1/14 so that it vanishes on D_1
def T(n):
    return hat_float(jumpsH, n) * cmath.exp(2j * pi * n / 14.0)

total = T(91) / 7.0
MTR = 40000
for m in range(-MTR, MTR + 1):
    if m == 0 or m % 7 == 0 or m % 13 != 0 or (91 + m) == 0:
        continue
    total += T(91 + m) * sin(pi * m / 7.0) / (pi * m)
check("C7: modulated-disjointness identity residual |sum| = %.2e < 1e-4 "
      "(truncation |m|<=%d)" % (abs(total), MTR), abs(total) < 1e-4)

print("C verdict: the sharp MSG-2152 hostile CAN be extended to carry the")
print("(7|m2, 13 not|m2) pair atom while keeping spectrum in 7Z u 13Z and")
print("blocking every 7-unit incident multiplier into 13Z -- and precisely")
print("then the trichotomy pivot fires: every admissible m1 is 13-divisible,")
print("so C-B = (m1-m2)*c is automatically a 91-unit edge between two atoms")
print("of the hostile itself. The extension CANNOT break the composition;")
print("the blocking property and the pair are arithmetically incompatible")
print("with the absence of a unit edge.")

# ---------------------------------------------------------------------------
print()
print("=" * 72)
print("PART D -- concrete 7|a canonical-shaped instance (c2,c3)=(91,169)")
print("=" * 72)

c2, c3 = 91, 169
b, c = nu13(c2), nu13(c3)
g = gcd(c2, c3)
a, dprime = c2 // g, c3 // g
check("D1: b=%d < c=%d, g=%d, a=%d (7|a), d'=%d, gcd(a,d')=1, 13|d'"
      % (b, c, g, a, dprime),
      b < c and g == 13 and a == 7 and dprime == 13
      and gcd(a, dprime) == 1 and dprime % 13 == 0)
check("D1: gcd(a,91)=%d > 1, so THM-2323 (31) FAILS and (40) blocks any "
      "unit edge inside the c2-multiple lattice" % gcd(a, 91),
      gcd(a, 91) == 7)

N = 91 * dprime
K0 = 1
z0, z1, h0, h1 = 0, 2, 0, 3
Kz0, Kz1 = K0 + dprime * z0, K0 + dprime * z1
assert gcd(Kz0, N) == 1 and gcd(Kz1, N) == 1
assert (z1 - z0) % 7 != 0 and (z1 - z0) % 13 != 0   # adjacency
q0, q1 = Kz0 + N * h0, Kz1 + N * h1
A, B = c2 * q0, c2 * q1
m2 = (B - A) // c3
check("D2: pair A=%d, B=%d, B-A=m2*c3 with m2=%d = a*bracket, "
      "bracket=%d unit mod 91, 7|m2, 13 not|m2"
      % (A, B, m2, (z1 - z0) + 91 * (h1 - h0)),
      B - A == m2 * c3 and m2 == a * ((z1 - z0) + 91 * (h1 - h0))
      and gcd((z1 - z0) + 91 * (h1 - h0), 91) == 1
      and m2 % 7 == 0 and m2 % 13 != 0)
charA = (A // 13 ** b) % 13
check("D2: A,B common-atom shape: nu13(A)=nu13(B)=%d=b, shared root "
      "character %d mod 13" % (b, charA),
      nu13(A) == b and nu13(B) == b
      and (B // 13 ** b) % 13 == charA)

# sweep every admissible THM-2326 multiplier m1 (7 not|m1, 1<=m1<=14S-1;
# take a generous range) and verify the trichotomy output on this instance
tri_ok = 0
for m1 in range(1, 400):
    if m1 % 7 == 0:
        continue
    C = A + m1 * c3
    assert nu13(C) == b and (C // 13 ** b) % 13 == charA
    if m1 % 13 != 0:
        assert gcd(m1, 91) == 1            # branch (i): A-C unit, A word-marked
    else:
        assert C != B
        assert gcd(m1 - m2, 91) == 1       # branch (ii): C-B unit, B word-marked
    tri_ok += 1
check("D3: all %d admissible m1 in [1,399], 7 not|m1: trichotomy delivers a "
      "grade/char-retentive 91-unit c3 edge with a word-marked endpoint"
      % tri_ok, tri_ok == 342)

print("D verdict: on the concrete 7|a instance the composition never fails;")
print("C = A + m1*c3 is the non-c2-multiple escape demanded by THM-2323 (40)")
print("exactly when 13|m1 forces the pivot (C is not divisible by d'=13 times")
print("the lattice step, i.e. C sits outside the c2-multiple affine coset).")

# ---------------------------------------------------------------------------
print()
print("=" * 72)
print("FINAL VERDICT")
print("=" * 72)
print("REFUTATION FAILED on all attack modes: the MSG-2153 two-colour")
print("trichotomy SURVIVES, subject to statement repairs (scope = positive-")
print("return middle-owner word strata of the 150 strict rows; the 7|a branch")
print("needs the observation that THM-2323 S7.2's fibre construction is")
print("unconditional and only (39)'s unit conclusion uses (31); the unit edge")
print("has exactly ONE word-marked endpoint in the 7|a branch, not two).")
print("ALL CHECKS PASSED")
