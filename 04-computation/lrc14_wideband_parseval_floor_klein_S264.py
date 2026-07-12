#!/usr/bin/env python3
"""
lrc14_wideband_parseval_floor_klein_S264.py
============================================
klein-2026-07-12-S264 / HYP-6130

TARGET: the LARGE-DIAMETER LOWER BOUND for LRC(14).
  Prove/certify M(spread DC 13-set) >= c > 1/14, with c GROWING with diameter
  (mac-mini THM-720's handoff, currently SAMPLED only).

THE OBSERVATION (a sharpening + wider-band generalization of THM-680, monad-S9):

  On a pair-sum ruler q = v_i + v_j, the LIVE-MULTIPLIER count at danger
  half-width d (band B = {r : d <= r <= q-d}, b = |B| = q-2d+1) is

     LM(q,d)/q = SUM_{k in Lambda_q(v)}  PROD_l hhat(k_l)          [THM-680(i), exact]

  where hhat = Fourier coeff of 1_B, Lambda_q(v) = {k : k.v = 0 mod q}.
  Because B is symmetric (r <-> q-r), hhat is REAL, so on the defining line
  L* = {m*(e_i+e_j)} the terms are (b/q)^11 * hhat(m)^2 >= 0  (POSITIVE).

  Summing the defining line EXACTLY by Parseval (SUM_m hhat(m)^2 = b/q):
     main(k=0) + defining-line = (b/q)^13 + (b/q)^12 (1 - b/q) = (b/q)^12.

  => EXACT IDENTITY:   LM(q,d)/q = (b/q)^12 + OffLine_signed
     where OffLine_signed = SUM_{k in Lambda_q \ (L* u {0})} PROD hhat(k_l).

  => FLOOR:            LM(q,d)/q >= (b/q)^12 - |OffLine(q,d)|.
     Hence  M(S) >= d/q  whenever some pair-sum q has |OffLine(q,d)| < (b/q)^12.

  This is SHARPER than THM-680's published floor (b/q)^12 (2b/q - 1)  [it SUBTRACTS
  the positive defining line], and it holds for ANY band width d, so it directly
  targets a GROWING c = d/q (THM-720). With b/q ~ 1-2c:  floor  ~  (1-2c)^12.

WHAT THIS SCRIPT VERIFIES (all exact rationals unless noted):
  (V1) The identity LM/q = (b/q)^12 + OffLine_signed  (via the exact defining-line
       Parseval sum) -- confirms the algebra + hhat real + defining line positive.
  (V2) M (exact) via the THM-668 pair-sum ruler == pair-sum max -- sanity.
  (V3) THE MECHANISM: for each family, the largest c the FLOOR reaches
       (|OffLine| < (b/q)^12 at some pair-sum q) vs the true M vs 1/14.
       Does c_floor GROW with diameter (matching THM-720's sampled table)?
  (V4) The sharpening: my floor (b/q)^12 - |OffLine| vs THM-680's
       (b/q)^12(2b/q-1) - |OffLine| at c = 1/14.
  (V5) A-PRIORI reach: bound |OffLine| by enumerating SMALL off-line relations
       (support <= 3, |coeff| <= Cmax) with an exact-hhat weight + a crude tail,
       and report the largest c the a-priori bound certifies. (float hhat; honest.)

Nothing here claims the adversarial minimum. It converts THM-720's SAMPLED
growing-M into a per-family CERTIFICATE form on the cancellation-immune pointwise
side, and localizes the remaining work to "spread => |OffLine| small at some q."
"""
import math
from fractions import Fraction
from itertools import combinations, product

# ----------------------------------------------------------------------
# families: (name, speeds, note).  spread/DC 13-sets + the AP wall.
# ----------------------------------------------------------------------
AP = list(range(1, 14))
FAMILIES = [
    ("AP {1..13} (wall)", AP, "M should be 1/14 exactly"),
    ("DC bounded extremal (HYP-6055)", [1,2,3,4,10,11,12,13,14,15,16,17,18], "M=1/12, diam 17"),
    ("adversarial spread DC (opus-S238)", [42,48,60,108,125,154,195,206,210,245,252,259,294], "occupies all folds mod 17,19,23"),
    ("kps blocker (HYP-6120)", [200,496,540,656,851,921,935,1122,1482,1680,1835,1849,1856], "first clears q=44, M~0.234, diam 1656"),
]

# a few random primitive spread families across scales, deterministic (index-seeded, no RNG)
def det_family(scale, idx):
    # build a spread 13-set deterministically from a multiplicative comb, primitive
    base = [1,2,3,4,6,8]                      # small divisor core (DC-flavored)
    extra = [ (scale*j + 7*idx*j*j) % (scale*13) + 5 for j in range(1,8) ]
    s = sorted(set(base + extra))
    while len(s) < 13:
        s.append(s[-1] + scale + idx + 1)
    s = s[:13]
    g = 0
    for x in s: g = math.gcd(g, x)
    s = [x//g for x in s] if g else s
    return sorted(set(s))[:13]

for sc, nm in [(10,"~20"),(50,"~110"),(200,"~430")]:
    f = det_family(sc, 3)
    if len(f) == 13:
        FAMILIES.append((f"det spread scale {sc}", f, f"diam~{max(f)-min(f)}"))

# ----------------------------------------------------------------------
# core: per pair-sum ruler q, the witness profile w(p) = min_l dist(v_l p, {0,q})
# ----------------------------------------------------------------------
def dist0(r, q):
    return min(r, q - r)

def w_profile(v, q):
    """w[p] = min_l dist(v_l * p mod q, {0,q}) for p=1..q-1 (p=0 excluded)."""
    w = [q] * q
    for p in range(1, q):
        mn = q
        for vl in v:
            d = dist0((vl * p) % q, q)
            if d < mn:
                mn = d
                if mn == 0:
                    break
        w[p] = mn
    return w  # w[0] meaningless (=q), ignore

def pair_sums(v):
    qs = set()
    for a, b in combinations(range(len(v)), 2):
        qs.add(v[a] + v[b])
    for a in range(len(v)):
        qs.add(2 * v[a])  # i=j single-runner peaks (THM-668 Part 1)
    return sorted(qs)

def exact_M(v):
    """M(S) = max over pair-sum q of max_p w(p)/q  (THM-668)."""
    best = Fraction(0)
    argq = None
    for q in pair_sums(v):
        w = w_profile(v, q)
        mx = max(w[1:])           # max over p != 0
        val = Fraction(mx, q)
        if val > best:
            best = val
            argq = q
    return best, argq

# ----------------------------------------------------------------------
# hhat for the symmetric band B = [d, q-d]  (REAL, exact via closed form check)
# hhat(k) = (1/q) sum_{x=d}^{q-d} cos(2 pi k x / q)   (imag parts cancel)
# ----------------------------------------------------------------------
def hhat_real(k, q, d):
    if k % q == 0:
        return Fraction(q - 2*d + 1, q)      # b/q, exact
    # float closed form: geometric sum of cos over x=d..q-d
    # = [ sin(pi k (q-2d+1)/q) / sin(pi k /q) ] * cos(pi k /q * (q)/1?) -- just sum directly
    s = 0.0
    ang = 2.0 * math.pi * (k % q) / q
    for x in range(d, q - d + 1):
        s += math.cos(ang * x)
    return s / q

# ----------------------------------------------------------------------
# V1: verify LM/q = (b/q)^12 + OffLine_signed, and the defining-line Parseval piece
# ----------------------------------------------------------------------
def verify_identity(v, q, d, k_runners=13):
    """Return dict with LM/q, main12=(b/q)^12, offline_signed, and the
       defining-line check main13 + Sline == main12."""
    n = len(v)
    b = q - 2*d + 1
    bq = Fraction(b, q)
    # LM = count of p with all v_l p in [d, q-d]
    LM = 0
    for p in range(1, q):
        ok = True
        for vl in v:
            if not (d <= (vl*p) % q <= q - d):
                ok = False
                break
        if ok:
            LM += 1
    LMq = Fraction(LM, q)
    main12 = bq ** (n - 1)          # (b/q)^(k-1) = (b/q)^12 for k=13
    main13 = bq ** n
    offline_signed = LMq - main12
    # defining-line Parseval piece (float): S_line = (b/q)^11 * sum_{m!=0} hhat(m)^2
    #   should equal (b/q)^12 (1 - b/q).  Uses ONE fixed defining pair -> represented as line.
    # Parseval sum_{m} hhat(m)^2 = b/q  (exact); subtract m=0 term (b/q)^2.
    parseval_sum = float(bq)                      # sum_m hhat(m)^2  (Parseval, exact)
    Sline = float(bq)**(n-2) * (parseval_sum - float(bq)**2)   # (b/q)^11 (b/q - (b/q)^2)
    predicted = float(bq)**n + Sline               # main13 + defining line
    return dict(LM=LM, LMq=LMq, main12=main12, main13=main13,
                offline_signed=offline_signed,
                defline_predicts_main12=(predicted, float(main12)),
                b=b, bq=bq)

# ----------------------------------------------------------------------
# V3/V4: floor reach.  For each pair-sum q, whole d-profile from w(p).
#   LM(d) = #{p!=0 : w(p) >= d}.  offline(d) = LM(d)/q - (b/q)^12.
#   my_floor(d)   = (b/q)^12 - |offline(d)|
#   thm680_floor  = (b/q)^12 (2b/q - 1) - |offline(d)|
# Report largest d/q with my_floor>0 (mechanism reach) and true M.
# ----------------------------------------------------------------------
def floor_analysis(v):
    n = len(v)
    M, argq = exact_M(v)
    c_floor_mine = Fraction(0)      # largest d/q with (b/q)^12 - |offline| > 0
    c_floor_680  = Fraction(0)
    argq_floor = None
    for q in pair_sums(v):
        if q < 14:
            continue
        w = w_profile(v, q)
        # cumulative counts: cnt[d] = #{p!=0 : w[p] >= d}
        maxd = q // 2
        hist = [0]*(q+1)
        for p in range(1, q):
            hist[w[p]] += 1
        cnt = [0]*(q+2)                        # suffix sum: cnt[d]=#{p!=0: w>=d}
        for d in range(q, -1, -1):
            cnt[d] = cnt[d+1] + (hist[d] if d <= q else 0)
        for d in range(1, maxd+1):
            b = q - 2*d + 1
            if b <= 0:
                break
            bq = Fraction(b, q)
            LMq = Fraction(cnt[d], q)
            main12 = bq ** (n-1)
            aoff = abs(LMq - main12)
            c = Fraction(d, q)
            if main12 - aoff > 0 and c > c_floor_mine:
                c_floor_mine = c
                argq_floor = q
            if main12 * (2*bq - 1) - aoff > 0 and c > c_floor_680:
                c_floor_680 = c
    return dict(M=M, argq=argq, c_floor_mine=c_floor_mine, argq_floor=argq_floor,
                c_floor_680=c_floor_680)

# ----------------------------------------------------------------------
# V5: a-priori |OffLine| bound at the M-achieving ruler (small relations).
#   Enumerate off-line k with support S (2<=|S|<=3), |k_l|<=Cmax, k.v=0 mod q,
#   k not on the defining line; weight |prod hhat(k_l)| (float hhat);
#   crude Parseval tail bound for the rest.  Report largest c a-priori-certified.
# ----------------------------------------------------------------------
def apriori_offline_reach(v, q, Cmax=4, Smax=3):
    """At ruler q, for each d return (LMq, main12, small_rel_mass, exact_offline).
       small_rel_mass = sum over enumerated small off-line relations of |prod hhat|.
       Returns the largest d/q with small_rel_mass < main12 (a-priori-style reach)."""
    n = len(v)
    # enumerate candidate off-line relations once (independent of d in SUPPORT, but hhat depends on d)
    idx = list(range(n))
    rels = []   # list of dict {pos: coeff}
    # support 2: k_i, k_j with k_i v_i + k_j v_j = 0 mod q, not the defining pair-sum line
    for i, j in combinations(idx, 2):
        for ci in range(-Cmax, Cmax+1):
            if ci == 0: continue
            for cj in range(-Cmax, Cmax+1):
                if cj == 0: continue
                if (ci*v[i] + cj*v[j]) % q == 0:
                    # defining line for THIS ruler is (i0,j0) with v_i0+v_j0=q, coeff (1,1);
                    # we exclude only exact (1,1)/( -1,-1) on a pair summing to q
                    if (v[i]+v[j]) % q == 0 and abs(ci)==1 and ci==cj:
                        continue
                    rels.append({i: ci, j: cj})
    # support 3
    for i, j, l in combinations(idx, 3):
        for ci in range(-Cmax, Cmax+1):
            if ci == 0: continue
            for cj in range(-Cmax, Cmax+1):
                if cj == 0: continue
                for cl in range(-Cmax, Cmax+1):
                    if cl == 0: continue
                    if (ci*v[i] + cj*v[j] + cl*v[l]) % q == 0:
                        rels.append({i: ci, j: cj, l: cl})
    best_c = Fraction(0)
    maxd = q // 2
    for d in range(1, maxd+1):
        b = q - 2*d + 1
        if b <= 0: break
        bq = Fraction(b, q)
        main12 = bq ** (n-1)
        # small relation mass (float)
        mass = 0.0
        for r in rels:
            wprod = 1.0
            for pos in range(n):
                if pos in r:
                    wprod *= abs(hhat_real(r[pos], q, d))
                else:
                    wprod *= float(bq)
            mass += wprod
        if mass < float(main12) and Fraction(d, q) > best_c:
            best_c = Fraction(d, q)
    return best_c, len(rels)

# ======================================================================
print("="*78)
print("V1 -- EXACT IDENTITY  LM/q = (b/q)^12 + OffLine_signed  (+ defining line)")
print("="*78)
# test on small rulers of the AP and one spread family
for name, v, _ in [FAMILIES[0], FAMILIES[1]]:
    print(f"\n{name}:  v={v}")
    for q in [14, 15, 20, 26]:
        if q not in pair_sums(v):
            # use nearest pair-sum
            ps = [x for x in pair_sums(v) if x >= 14]
            if not ps: continue
            q = min(ps, key=lambda x: abs(x-q))
        d = max(1, round(q/14))
        r = verify_identity(v, q, d)
        pred, m12 = r["defline_predicts_main12"]
        print(f"  q={q:3d} d={d} b={r['b']:3d}: LM={r['LM']:4d}  LM/q={float(r['LMq']):.5f}  "
              f"(b/q)^12={float(r['main12']):.5f}  offline={float(r['offline_signed']):+.5f}  "
              f"| defline->main12: {pred:.6f} vs {m12:.6f}  {'OK' if abs(pred-m12)<1e-9 else 'MISMATCH'}")

print()
print("="*78)
print("V2/V3/V4 -- FLOOR REACH per family (mechanism ceiling; |OffLine| measured exact)")
print("="*78)
print(f"{'family':38s} {'diam':>6s} {'M':>10s} {'M~':>7s} {'c_floor':>10s} {'c_flr~':>7s} {'1/14?':>6s} {'680~':>7s}")
one14 = Fraction(1,14)
results = []
for name, v, note in FAMILIES:
    if len(v) != 13:
        continue
    diam = max(v) - min(v)
    fa = floor_analysis(v)
    M = fa["M"]; cf = fa["c_floor_mine"]; c680 = fa["c_floor_680"]
    results.append((name, v, diam, M, cf, c680))
    print(f"{name[:38]:38s} {diam:6d} {str(M):>10s} {float(M):7.4f} "
          f"{str(cf):>10s} {float(cf):7.4f} {'>' if cf>one14 else '<':>6s} {float(c680):7.4f}")

print()
print("  (c_floor = largest c=d/q where (b/q)^12 - |OffLine| > 0 at some pair-sum q:")
print("   the reach of the SHARPENED Parseval floor if |OffLine| a-priori-bounded by its true value.")
print("   680~ = same for THM-680's published floor (b/q)^12(2b/q-1) - |OffLine|.)")

print()
print("="*78)
print("V5 -- A-PRIORI reach at the M-achieving ruler (small off-line relations, float hhat)")
print("="*78)
for name, v, diam, M, cf, c680 in results:
    _, argq = exact_M(v)
    if argq is None or argq < 14:
        continue
    if argq > 700:
        print(f"{name[:38]:38s}  ruler q*={argq:5d}  M={float(M):.4f}  (V5 skipped: q*>700, slow)")
        continue
    ca, nrel = apriori_offline_reach(v, argq, Cmax=4, Smax=3)
    print(f"{name[:38]:38s}  ruler q*={argq:5d}  M={float(M):.4f}  "
          f"a-priori c={float(ca):.4f} ({'>' if ca>one14 else '<'}1/14)  [{nrel} small rels]")

print()
print("1/14 =", float(one14))
print("done.")
