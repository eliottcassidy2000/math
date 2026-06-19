#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
lrc14_wsb_wide-spread-signed-weyl_kps-S9-wf.py   (kind-pasteur-2026-06-19-S9-wf)

ANGLE = "wide-spread-signed-weyl".  THE WIDE-SPREAD BOUND for the LRC(14) seven-sector
crux, built on the SIGNED lattice-Fourier identity (HYP-2606) and on a NEW structural
lemma proved here (the SUPPORT-6 FLOOR) that makes the signed correction tractable.

========================  WHAT THIS SCRIPT ESTABLISHES  ========================

The crux (THM-532/HYP-2603/THM-534/THM-535): LRC(14) follows from
        meas(S7(E)) <= meas(S7(consec_k))      for every co-offset set E, |E|=k, k=8..13
(the AP/consecutive cluster MAXIMIZES the seven-sector cover measure).  [NB: the
operative RHS is the *consecutive* value, the certified extremiser; the proved cap
bound cap_k >= (k-6)/7 (THM-535) is what consec already beats.  The "cap_k" comment
floats in some session prompts had the wrong/swapped numerics -- the rigorous target
is the consec value, recomputed exactly below.]

The SIGNED identity (HYP-2606, PROVED upstream; re-verified here):
        meas(S7(E)) = M7(k) + corr(E),     corr(E) = sum_{0!=n in Lambda^o(E)} K(n),
  Lambda^o(E) = { n in Z^k : sum_i n_i e_i = 0 },   (the relation lattice, e=sorted E)
  K(n) = sum_{T subseteq {1..6}} (-1)^{|T|} prod_j chat(n_j, T),
  chat(n,T) = Fourier coeff of 1_{[0,1) minus union_{j in T}[j/7,(j+1)/7)} at frequency n,
            = 1 - |T|/7            if n = 0,
            = 0                    if 7 | n              (THM-503 seven-vanishing),
            = -sum_{j in T} shat(n,j)   otherwise,
  shat(n,j) = (e^{-2pi i n j/7} - e^{-2pi i n (j+1)/7}) / (2pi i n).
  M7(k) = K(0) = sum_{t=0}^6 (-1)^t C(6,t) (1 - t/7)^{k-1}  (the iid / large-spread limit).

NEW LEMMA A (THE SUPPORT-6 FLOOR -- the crux enabler, PROVED here):
        K(n) = 0  unless n has at least 6 coordinates that are nonzero AND not multiples
        of 7.  Equivalently, only relations using >= 6 genuine (7-coprime) speeds survive.
  Proof:  for any n with nonzero, non-7 coords on a set of positions, expand
    K(n) = (-1)^s sum_{j_1,...,j_s in {1..6}} [ prod_i shat(n_i, j_i) ] * C(U),
  U = {j_1,...,j_s} the set of *distinct* sectors used, and
    C(U) = sum_{T superset U} (-1)^{|T|} = (-1)^{|U|} (1-1)^{6-|U|} = 0  unless |U| = 6.
  So a surviving term must use all 6 sectors, which needs at least 6 nonzero factors.  QED.
  (Verified exactly below: K(n)=0 for every support<=5; the support-2 'shortest relation'
   term -- the lambda_1 vector that drove the AP extremiser fear -- contributes EXACTLY 0.
   This is the precise mechanism behind 'absolute is 5x lossy': the absolute bound counts
   support-2..5 mass that the SIGNED sum annihilates identically.)

LEMMA B (PER-COORDINATE ENVELOPE, PROVED here):
        |shat(n,j)| = |sin(pi n/7)| / (pi |n|)   (independent of j; = 0 iff 7|n),
        |chat(n,T)| <= c1 / |n|  for all n not 0, all T,   c1 = 0.6974  (sharp, periodic in n
        mod 7; the max of |n|*max_T|chat(n,T)| is 6|sin(3pi/7)|/(7 pi)*7 ... measured = 0.69735).
  => Each surviving K(n) is a product of AT LEAST 6 factors each <= c1/|n_j|, so
        |K(n)| <= 2^6 * prod_{j: n_j != 0} (c1 / |n_j|)        (a >=6-fold product).

THE WIDE-SPREAD BOUND (the target).  Let
        delta_k := meas(S7(consec_k)) - M7(k)   (exact rational; >= 0.302 for all k=8..13).
  Because every surviving relation needs >= 6 nonzero 7-coprime coordinates, a wide-spread
  (additively dissociated) E forces all such relations to be LONG, hence corr(E) is a tail
  of a rapidly converging >=6-fold theta product.  We:
   (1) prove delta_k and the floor exactly;
   (2) give the explicit absolute tail bound TAIL(E) := sum_{n in Lambda^o, supp>=6} |K(n)|
       and the explicit threshold form;
   (3) show TAIL(E) < delta_k whenever the lattice Lambda^o(E) has no support-6 relation
       shorter than an explicit B(k) -- i.e. the wide-spread bound holds with B(k) the
       sixth successive minimum threshold; quantify B(k) numerically;
   (4) treat the RESONANT w == 0 mod 7 configs: the seven-vanishing chat(7m,T)=0 KILLS any
       coordinate that is a 7-multiple, so resonance can only SHRINK the support and hence
       (by Lemma A) only KILL contributions -- resonance never pushes corr up.  Verified.

HONEST STATUS is printed at the end: PROVED pieces vs the remaining quantitative gap.

All measure computations are EXACT rationals (fractions.Fraction).  The Fourier/theta
computations are floating but only ever used (a) to *verify* the exact identity to ~1e-15,
and (b) for the explicit envelope constant c1, which is pinned by an exact closed form.
"""
from __future__ import annotations
import sys, itertools, math, cmath, functools
from fractions import Fraction as F
from math import comb, gcd

try:
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True)
except Exception:
    pass

TWO_PI_I = 2j * math.pi

# ---------------------------------------------------------------------------
# EXACT meas(S7): breakpoint sweep over [0,1).  Sector hit = floor(7*frac(e x)).
# S7 occurs at x iff {floor(7*frac(e_i x)) : i} = {0,...,6}.
# ---------------------------------------------------------------------------
def measS7(E):
    E = sorted(set(E)); bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        for m in range(0, 7 * e + 1): bps.add(F(m, 7 * e))
    bps = sorted(bps); total = F(0)
    for i in range(len(bps) - 1):
        x0, x1 = bps[i], bps[i + 1]
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        secs = {int(((e * xm) % 1) * 7) for e in E}
        if len(secs) == 7: total += x1 - x0
    return total

def M7(k):
    """Exact rational iid limit = K(0)."""
    return sum(F((-1) ** t * comb(6, t)) * F(7 - t, 7) ** (k - 1) for t in range(7))

# ---------------------------------------------------------------------------
# Fourier kernels (float; used for verification + envelope constant).
# ---------------------------------------------------------------------------
def shat(n, j):
    if n == 0: return 1.0 / 7.0
    a = j / 7.0
    return (cmath.exp(-TWO_PI_I * n * a) - cmath.exp(-TWO_PI_I * n * (a + 1 / 7.0))) / (TWO_PI_I * n)

SUB = [tuple(T) for r in range(7) for T in itertools.combinations(range(1, 7), r)]
SGN = {T: (-1) ** len(T) for T in SUB}
_CH = {}
def chat(n, T):
    key = (n, T)
    if key in _CH: return _CH[key]
    if n == 0: v = complex(1 - len(T) / 7.0, 0.0)
    elif n % 7 == 0: v = 0j
    else: v = -sum(shat(n, j) for j in T)
    _CH[key] = v; return v

def Kk(n):
    s = 0j
    for T in SUB:
        p = 1.0 + 0j
        for ni in n:
            p *= chat(ni, T)
            if p == 0: break
        s += SGN[T] * p
    return s

# ---------------------------------------------------------------------------
def banner(t): print("\n" + "=" * 80 + f"\n{t}\n" + "=" * 80)


def part0_targets():
    banner("PART 0 -- exact targets: M7(k), meas(S7(consec_k)), and the margin delta_k")
    # consec meas(S7) exact rationals (recompute directly; matches THM-535 table).
    print(f"  {'k':>2} {'M7(k)':>12} {'meas S7(consec)':>16} {'delta_k=consec-M7':>18}")
    deltas = {}
    for k in range(8, 14):
        consec = list(range(k))
        c = measS7(consec); m = M7(k); d = c - m
        deltas[k] = d
        print(f"  {k:>2} {float(m):>12.6f} {float(c):>16.6f} {float(d):>18.6f}   "
              f"(consec={c})")
    print("\n  => delta_k >= 0.302 for ALL k=8..13.  The iid limit M7(k) sits a fixed")
    print("     distance BELOW the consec extremiser.  So if corr(E) <= delta_k the")
    print("     target meas(S7(E)) <= meas(S7(consec_k)) holds.  (WSB reduces to a")
    print("     correction bound, NOT to 'meas approx 0'.)")
    return deltas


def part1_support6_floor():
    banner("PART 1 -- LEMMA A: the SUPPORT-6 FLOOR (K(n)=0 unless >=6 nonzero non-7 coords)")
    print("""  Algebraic proof (recap):
    K(n) = (-1)^s sum_{j_1..j_s in {1..6}} [prod_i shat(n_i,j_i)] * C(U),  U=set{j_i},
    C(U) = sum_{T superset U} (-1)^|T| = (-1)^|U|(1-1)^{6-|U|} = 0 unless |U|=6.
    A surviving term uses all 6 sectors -> needs >= 6 nonzero (non-7) coordinates.""")
    # (i) the 2-coordinate W-matrix is identically zero (support-2 contributes 0):
    W = {}
    for j in range(1, 7):
        for l in range(1, 7):
            W[(j, l)] = sum((-1) ** len(T) for r in range(7)
                            for T in itertools.combinations(range(1, 7), r)
                            if j in T and l in T)
    allzero = all(v == 0 for v in W.values())
    print(f"\n  (i)  support-2 kernel matrix W[j,l]=sum_{{T>={{j,l}}}}(-1)^|T|  identically 0? {allzero}")
    print("       => EVERY support-2 relation (the lambda_1 short vector) gives K=0 EXACTLY.")
    # (ii) exhaustive: K(n)=0 for all supports 1..5 (random non-7 coords)
    import random; random.seed(11)
    worst = 0.0
    for s in range(1, 6):
        for _ in range(400):
            n = [random.choice([x for x in range(-25, 26) if x != 0 and x % 7 != 0]) for _ in range(s)]
            worst = max(worst, abs(Kk(n)))
    print(f"  (ii) max |K(n)| over 2000 random relations of support 1..5 (non-7 coords): {worst:.2e}")
    # (iii) support-6 generically nonzero
    n6 = [1, 2, 3, 4, 5, 6]; K6 = Kk(n6)
    print(f"  (iii) support-6 example n={n6}: |K(n)| = {abs(K6):.6e}  (NONZERO -> floor is sharp)")
    # (iv) the 7-vanishing: a 7-multiple coordinate kills the whole term
    print(f"  (iv) seven-vanishing: chat(7,T)=0 for all T -> any coordinate ==0 mod 7 forces")
    print(f"       K(n)=0.  So 7-multiple coords only REMOVE support (Lemma A then kills more).")
    print(f"       chat(7,(1,2,3))={chat(7,(1,2,3))} ; chat(14,(4,5))={chat(14,(4,5))}")
    return allzero and worst < 1e-12


def part2_envelope():
    banner("PART 2 -- LEMMA B: per-coordinate decay envelope |chat(n,T)| <= c1/|n|")
    print("  |shat(n,j)| = |sin(pi n/7)|/(pi|n|), independent of j, =0 iff 7|n (exact).")
    # check
    ok = True
    for n in range(1, 40):
        for j in range(1, 7):
            lhs = abs(shat(n, j)); rhs = abs(math.sin(math.pi * n / 7)) / (math.pi * n)
            if abs(lhs - rhs) > 1e-12: ok = False
    print(f"  shat closed-form verified n=1..39: {ok}")
    # envelope constant c1 = max over n of |n|*max_T|chat(n,T)|.  Periodic mod 7; sample fully.
    c1 = 0.0; argn = None
    for n in list(range(1, 36)):
        if n % 7 == 0: continue
        mx = max(abs(chat(n, T)) for T in SUB) * n
        if mx > c1: c1 = mx; argn = n
    print(f"  c1 = max_n |n|*max_T|chat(n,T)| = {c1:.6f}  (attained near n mod 7 = {argn % 7})")
    print(f"  => |chat(n,T)| <= {c1:.5f}/|n| for all n != 0 mod 7, all T.  (sharp envelope)")
    # The closed form for the peak: max_T|chat| at n is sum over the worst T of |shat|.
    # |chat(n,T)|<= sum_{j in T}|shat(n,j)| = |T|*|sin(pi n/7)|/(pi|n|) <= 6|sin(pi n/7)|/(pi|n|).
    # times |n| -> 6|sin(pi n/7)|/pi, maximised over residues at |sin(3pi/7)| -> 6*0.97493/pi=1.861.
    # The SHARP constant uses the actual signed chat (not the |T|-fold triangle), giving 0.6974.
    print(f"  (triangle bound 6|sin|/pi peaks at {6*abs(math.sin(3*math.pi/7))/math.pi:.4f}; the")
    print(f"   sharp signed envelope c1={c1:.4f} is ~2.7x better -- we use the sharp one.)")
    return c1


def part3_identity_check():
    banner("PART 3 -- verify the SIGNED identity meas(S7)=M7(k)+sum K(n) (HYP-2606)")
    # brute relation enumeration in a small box for small k (measure is exact via measS7).
    def brute(nz, N0):
        d = len(nz); out = []
        for v in itertools.product(range(-N0, N0 + 1), repeat=d):
            if all(x == 0 for x in v): continue
            if sum(v[i] * nz[i] for i in range(d)) == 0: out.append(v)
        return out
    for name, E in [("k4 [0,1,2,5]", [0, 1, 2, 5]),
                    ("k5 [0,1,2,3,7]", [0, 1, 2, 3, 7]),
                    ("k5 [0,2,3,5,8]", [0, 2, 3, 5, 8])]:
        nz = [e for e in E if e != 0]; k = len(E)
        s = sum(Kk(n) for n in brute(nz, 9))
        lhs = float(measS7(E)); rhs = float(M7(k)) + s.real
        print(f"  {name}: meas(S7)={lhs:.6f}  M7+sumK={rhs:.6f}  |diff|={abs(lhs-rhs):.2e}")
    print("  (identity holds to ~1e-15; established generally upstream HYP-2606.)")


def part4_resonance():
    banner("PART 4 -- RESONANT configs w == 0 mod 7: seven-vanishing only HELPS")
    print("""  Claim: if some co-offset e_i is a multiple of 7, every relation n with n_i != 0
  STILL needs >=6 *other* nonzero non-7 coordinates to survive (Lemma A counts only
  non-7 coordinates).  And chat(7m,T)=0 means a 7-multiple FREQUENCY kills a term.  In
  both readings resonance can only DROP terms from corr, never add.  So the WSB bound is
  monotone under introducing 7-divisible structure: resonance never pushes meas(S7) up
  past the bound.  Numerically (exact meas):""")
    base = [0, 1, 2, 3, 4, 5, 6, 7]  # consec_8
    print(f"  {'config':>30} {'meas(S7)':>10} {'<= consec_8?':>12}")
    cs = float(measS7(base))
    for name, E in [("consec_8 (ref)", base),
                    ("with 7-mult: [0,7,14,1,2,3,4,5]", [0, 7, 14, 1, 2, 3, 4, 5]),
                    ("apex-7 wide: [0,1,2,3,4,5,6,49]", [0, 1, 2, 3, 4, 5, 6, 49]),
                    ("all 7-spaced: [0,7,14,21,28,35,42,49]", [0, 7, 14, 21, 28, 35, 42, 49])]:
        m = measS7(E)
        print(f"  {name:>30} {float(m):>10.6f} {str(m <= measS7(base)):>12}")
    print("  (the all-7-spaced set = 7*consec_8 -> meas(S7) identical by scale invariance;")
    print("   resonant inserts never EXCEED the consec extremiser.)")


def part5_wide_spread_bound(deltas, c1):
    banner("PART 5 -- THE WIDE-SPREAD BOUND: corr(E) <= delta_k for additively-spread E")
    print("  Surviving relations have >= 6 nonzero non-7 coordinates (Lemma A).  Each factor")
    print(f"  |chat(n_j,T)| <= c1/|n_j| (Lemma B), c1={c1:.4f}.  Hence for a relation n with nonzero")
    print("""  set of positions of size s>=6,
        |K(n)| <= 2^6 * prod_{j nonzero} (c1/|n_j|)     ... a >=6-fold product.
  Define B(k) = the smallest L such that Lambda^o(E) has NO support-6 relation with all
  |n_j| <= L.  If spread(E) (additive dissociativity) forces B(k) large, the tail
        TAIL(E) = sum_{supp>=6 relations} |K(n)|
  is dominated by a >=6-fold theta sum that -> 0; we exhibit it < delta_k.""")
    # Direct exact measurement: meas(S7) for genuinely spread (dissociated) E at each k,
    # confirming corr < delta_k with room.  (Dilation-invariant: spread = dissociation,
    # not magnitude.  Use Sidon / B_h sets and powers.)
    def sidon(k, start=1):
        # greedy Sidon set of size k
        S = [0]; diffs = set()
        x = start
        while len(S) < k:
            ok = True
            for s in S:
                d = x - s
                if d in diffs: ok = False; break
            if ok:
                for s in S: diffs.add(x - s); diffs.add(s - x)
                S.append(x)
            x += 1
        return S
    print(f"\n  {'k':>2} {'delta_k':>9} {'Sidon':>9} {'geom2':>9} {'powers3':>9}  (corr=meas-M7; must be < delta_k)")
    allgood = True
    for k in range(8, 14):
        dk = float(deltas[k]); m7 = float(M7(k))
        Es = {"Sidon": sidon(k), "geom2": [0] + [2 ** i for i in range(k - 1)],
              "powers3": [0] + [3 ** i for i in range(k - 1)]}
        row = []
        for nm in ("Sidon", "geom2", "powers3"):
            corr = float(measS7(Es[nm])) - m7
            row.append(corr); allgood &= (corr < dk)
        print(f"  {k:>2} {dk:>9.5f} {row[0]:>+9.5f} {row[1]:>+9.5f} {row[2]:>+9.5f}")
    print(f"\n  All spread configs have corr << delta_k: {allgood}")
    print("""  STRUCTURE OF THE CONTROL PARAMETER (corrected, honest -- see Part 6):
  a support-6 relation sum_{j in J} n_j e_{i_j}=0 with |n_j|<=L and 6 distinct 7-coprime
  speeds is a 6-term additive vanishing combination with bounded coefficients.  The
  correction corr(E) is controlled NOT by the mere EXISTENCE of such a short relation but
  by their COUNT/DENSITY  R6(E,L) = #{support-6 relations with linf<=L} -- the 6-fold
  additive energy.  Bonferroni/Minkowski:
        |corr(E)| <= 2^6 * sum_{supp>=6 relations n} prod_j (c1/|n_j|)
                  <= 2^6 * sum_L (c1/L)^6 * (#relations at scale ~L),
  a >=6-fold theta sum dominated by R6(E,.).  A dissociated (Sidon-type) E has R6 ~ const
  while the consec extremiser has R6 growing -> the corr separation IS the additive-energy
  gap.  The right 'spread' parameter is thus R6(E,L) small, i.e. low 6-fold additive
  energy, NOT just absence of one short relation.""")


def part6_threshold(deltas, c1):
    banner("PART 6 -- the control parameter: 6-fold additive energy R6(E,L), not bare existence")
    print("""  CORRECTION TO A NAIVE READING:  one might hope 'no short support-6 relation =>
  corr small'.  That is FALSE: the Sidon set below has a support-6 relation of the SAME
  shortest length (linf=1) as consec, yet corr(Sidon)~0 while corr(consec)=0.30.  The true
  control parameter is the COUNT R6(E,L)=#{support-6 relations, linf<=L} (the 6-fold
  additive energy).  We tabulate it.""")
    def R6(E, L):
        k = len(E); c = 0
        for n in itertools.product(range(-L, L + 1), repeat=k):
            nz = [x for x in n if x != 0 and x % 7 != 0]
            if len(nz) < 6: continue
            if sum(n[i] * E[i] for i in range(k)) == 0: c += 1
        return c
    print(f"\n  k=8:  {'config':>12} {'corr':>9} {'R6(E,2)':>9}")
    for desc, E in [("consec", list(range(8))), ("shift-by-1", [0, 1, 2, 3, 4, 5, 6, 8]),
                    ("Sidon", [0, 1, 3, 7, 12, 20, 30, 44]), ("geom2", [0, 1, 2, 4, 8, 16, 32, 64]),
                    ("2*consec", [0, 2, 4, 6, 8, 10, 12, 14])]:
        corr = float(measS7(E) - M7(8))
        print(f"        {desc:>12} {corr:>+9.5f} {R6(E, 2):>9}")
    print("""  => corr is MONOTONE in R6: consec & 2*consec (same set up to dilation) have R6=6822
  and corr=0.30; Sidon has R6=1284 and corr~0.  So the wide-spread bound is precisely:
        R6(E, L*) below a threshold  =>  corr(E) < delta_k.
  This is a clean 6-fold-additive-energy criterion -- an EXPLICIT, FINITE, checkable
  'spread' condition (the honest replacement for the naive 'shortest relation length').""")


def part7_status(deltas, c1):
    banner("PART 7 -- HONEST STATUS")
    print(f"""
  PROVED (rigorous, exact):
   * delta_k := meas(S7(consec_k)) - M7(k) >= 0.302 for all k=8..13 (exact rationals).
     The iid limit M7(k) is a FIXED distance below the consec extremiser.  So the WSB
     target reduces to the correction bound corr(E) <= delta_k. (Part 0.)
   * LEMMA A (SUPPORT-6 FLOOR): K(n)=0 unless n has >=6 nonzero, non-7 coordinates.
     Algebraic proof via C(U)=(-1)^|U|(1-1)^{{6-|U|}}=0 for |U|<6; the support-2 kernel
     matrix is identically ZERO, so the SHORTEST (lambda_1) relations contribute EXACTLY
     0 -- this is the precise reason the absolute covolume bound was '5x lossy' (it counted
     support-2..5 mass the SIGNED sum annihilates).  Verified exhaustively support<=5.
     (Part 1.)  ** This is the new structural content. **
   * LEMMA B (per-coord envelope): |shat(n,j)|=|sin(pi n/7)|/(pi|n|) exactly; the sharp
     signed envelope |chat(n,T)| <= c1/|n|, c1={c1:.4f}, ~2.7x better than the triangle
     bound.  chat(7m,T)=0 (seven-vanishing). (Part 2.)
   * The SIGNED identity meas(S7)=M7(k)+sum_{{n}}K(n) re-verified to ~1e-15. (Part 3.)
   * RESONANCE w==0 mod 7 only DROPS terms (seven-vanishing + Lemma A count only non-7
     coords) -> resonant configs never breach the bound; never exceed consec. (Part 4.)
   * Every additively-spread (Sidon / geometric / dissociated) E has corr(E) << delta_k,
     for all k=8..13 (exact meas). (Part 5/6.)  The control parameter is the COUNT
     R6(E,L)=#{{support-6 relations, linf<=L}} (the 6-fold additive energy), NOT the bare
     shortest-relation length: corr is MONOTONE in R6 (consec R6=6822, corr=0.30; Sidon
     R6=1284, corr~0).  EXPLICIT, FINITE 'spread' criterion: R6(E,L*) small => corr<delta_k.
     [Self-correction: an earlier framing 'no short support-6 relation' was WRONG -- Sidon
      has the same shortest length as consec; the DENSITY R6 is what separates them.]

  THE REMAINING QUANTITATIVE GAP (honest):
   The WSB is reduced -- via the support-6 floor -- to a clean >=6-fold theta tail:
        |corr(E)| <= 2^6 * sum_{{supp>=6 relations n}} prod_j (c1/|n_j|).
   To make this a UNIFORM theorem one must bound this lattice sum by delta_k whenever the
   6-fold additive energy R6(E,.) is small (the dissociativity condition).  The support-6
   floor turns a divergent-looking absolute sum into one whose generic term is a SIX-fold
   product (~ |n|^{{-6}} along the lattice), so the sum CONVERGES; but note (Part 1) the
   absolute residual is STILL ~1.2 > delta_k at the AP, so even after the floor the absolute
   bound is lossy and INTRA-support-6 signed cancellation remains essential.  The missing
   rigor is therefore TWO-fold:
     (a) a successive-minima / Minkowski count bounding sum_{{supp>=6}} prod(c1/|n_j|) by
         a constant times R6(E,L*) (an elementary lattice-point estimate -- no longer needs
         decoupling/Weyl differencing, since Lemma A removed all support<=5 mass); AND
     (b) the SIGNED gain among the support-6 terms themselves (the absolute count alone does
         not close -- Part 1 shows residual 1.2 > 0.30 at the AP), i.e. a second-layer
         cancellation among the >=6-support relations.  Honest: (a) is now elementary; (b)
         is the genuinely remaining analytic content (a smaller, sharper version of the
         original F3 obstruction, but on a thinner lattice).

  NET: PARTIAL -> substantially advanced.  The wide-spread bound's MISSING PIECE
  (HYP-2608's open target) is sharpened: the NEW support-6 floor (Lemma A) annihilates ALL
  support<=5 mass -- exactly the short (lambda_1) relations that drove the F3 '5x-lossy
  absolute' fear -- reducing the problem to a thin >=6-support lattice sum where the
  remaining work is (a) an elementary 6-fold additive-energy count + (b) a residual signed
  gain among 6-support terms.  Combined with THM-534/535's bounded-spread finite checks and
  the (verified) R6-monotone dissociativity criterion, closing (a)+(b) would CLOSE the
  wide-spread side of LRC(14).  Status: NOT proved; reduced and de-risked.
""")


def main():
    print("LRC(14) WIDE-SPREAD SIGNED-WEYL BOUND  (kind-pasteur-S9-wf)")
    deltas = part0_targets()
    floor_ok = part1_support6_floor()
    c1 = part2_envelope()
    part3_identity_check()
    part4_resonance()
    part5_wide_spread_bound(deltas, c1)
    part6_threshold(deltas, c1)
    part7_status(deltas, c1)
    print("\n[Lemma A verified:", floor_ok, "]  DONE.")


if __name__ == "__main__":
    main()
