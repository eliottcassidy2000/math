#!/usr/bin/env python3
"""
coverage_spectrum_bonferroni_corona_klein_S313.py — klein-2026-07-15-S313 (cont.4)

THE COVERAGE SPECTRUM FRAME for the LRC(14) covering case (LEM-013 candidate).

For a cluster E = {e_1..e_j} at slow-variable x, the danger arcs A_i = [frac(e_i x),
frac(e_i x) + 1/7) have coverage function C(t) = #{i : t in A_i} and SPECTRUM
mu_k(x) = meas{t : C(t) = k}.  Exact rational sweep (x rational => endpoints rational).

Sum rules and identities (all verified exactly here):
  (R0) sum_k mu_k = 1
  (R1) sum_k k mu_k = j/7                          (first moment = total arc mass)
  (R2) mu_0 = (7-j)/7 + sum_{k>=2} (k-1) mu_k      (uncovered = overlap excess - (j-7)/7)
       => at covering x (mu_0 = 0, j = 13): overlap excess EXACTLY 6/7 — the adversary's budget
  (S)  S_d = sum_k C(k,d) mu_k                     (subset-sum = spectrum moment; d<=3 cross-checked
                                                    against literal subset enumeration)
  (T)  sum_{d<=D} (-1)^d S_d = mu_0 + (-1)^D sum_{k>=1} C(k-1,D) mu_k
       => for odd D: BONFERRONI B_D = mu_0 - deficit, deficit = sum_{k>=D+1-ish} C(k-1,D) mu_k >= 0
       => **B_D is EXACT iff max multiplicity <= D** (the Kakeya-tightness criterion transferred:
          K(A5) = 15 was the multiplicity-<=2 model case where Bonferroni-2 is tight)

Roots-of-unity prediction (the solvable=tight inversion): deep multiplicity (champion overlap
points) concentrates at LOW-ORDER x = a/q; generic (high-order) x has near-binomial spectrum.
"""
from fractions import Fraction as Fr
from math import comb
from itertools import combinations
import random

W = Fr(1, 7)

def spectrum(E, x):
    """exact mu_k for arcs [frac(e x), frac(e x)+1/7) on R/Z; also max multiplicity."""
    events = []
    for e in E:
        a = (e * x) % 1
        b = a + W
        if b <= 1:
            events.append((a, 1)); events.append((b, -1))
        else:
            events.append((a, 1)); events.append((Fr(1), -1))
            events.append((Fr(0), 1)); events.append((b - 1, -1))
    pts = sorted(set([Fr(0), Fr(1)] + [p for p, _ in events]))
    mu = {}
    for i in range(len(pts) - 1):
        t0, t1 = pts[i], pts[i + 1]
        c = sum(s for p, s in events if p <= t0) if False else None
        cnt = 0
        for e in E:
            a = (e * x) % 1
            if a + W <= 1:
                if a <= t0 < a + W: cnt += 1
            else:
                if t0 >= a or t0 < a + W - 1: cnt += 1
        mu[cnt] = mu.get(cnt, 0) + (t1 - t0)
    return mu

def S_d_moment(mu, d): return sum(comb(k, d) * m for k, m in mu.items())

def S_d_literal(E, x, d):
    """intersection measure summed over d-subsets (cross-check, d small)."""
    def arc(e):
        a = (e * x) % 1
        return [(a, min(a + W, Fr(1)))] + ([(Fr(0), a + W - 1)] if a + W > 1 else [])
    tot = Fr(0)
    for T in combinations(E, d):
        segs = arc(T[0])
        for e in T[1:]:
            newsegs = []
            for s0, s1 in segs:
                for u0, u1 in arc(e):
                    lo, hi = max(s0, u0), min(s1, u1)
                    if lo < hi: newsegs.append((lo, hi))
            segs = newsegs
            if not segs: break
        tot += sum(s1 - s0 for s0, s1 in segs)
    return tot

OK = []
def check(name, cond):
    OK.append((name, bool(cond)))
    print(("PASS" if cond else "FAIL"), name)

AP13 = list(range(1, 14))
FIB13 = [1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377]
rng = random.Random(14)
RAND13 = sorted(rng.sample(range(1, 200), 13))
INSTANCES = {"tightAP{1..13}": AP13, "AP12+84": list(range(1, 13)) + [84],
             "Fibonacci13": FIB13, f"random13": RAND13}

# ---- identity suite at assorted x ----
ok_R = ok_S = ok_T = True
for name, E in INSTANCES.items():
    j = len(E)
    for x in [Fr(1, 7), Fr(1, 13), Fr(1, 14), Fr(3, 91), Fr(5, 101), Fr(377, 1009), Fr(233, 610)]:
        mu = spectrum(E, x)
        if sum(mu.values()) != 1: ok_R = False
        if sum(k * m for k, m in mu.items()) != Fr(j, 7): ok_R = False
        lhs = mu.get(0, Fr(0))
        rhs = Fr(7 - j, 7) + sum((k - 1) * m for k, m in mu.items() if k >= 2)
        if lhs != rhs: ok_R = False
        for d in (2, 3):
            if S_d_moment(mu, d) != S_d_literal(E, x, d): ok_S = False
        for D in (3, 4, 5):
            tr = sum((-1) ** d * S_d_moment(mu, d) for d in range(D + 1))
            tgt = mu.get(0, Fr(0)) + (-1) ** D * sum(comb(k - 1, D) * m for k, m in mu.items() if k >= 1)
            if tr != tgt: ok_T = False
check("(R0)(R1)(R2) sum rules exact on 4 instances x 7 rational x", ok_R)
check("(S) S_d = sum C(k,d) mu_k == literal d-subset intersection sums (d=2,3)", ok_S)
check("(T) truncation identity exact (D=3,4,5); odd-D deficit = C(k-1,D)-moment >= 0", ok_T)

# ---- covering-slack: at covering x the overlap excess is exactly 6/7 ----
cov_found = 0; slack_ok = True
for name, E in INSTANCES.items():
    for q in range(7, 120):
        x = Fr(1, q)
        mu = spectrum(E, x)
        if mu.get(0, Fr(0)) == 0:
            cov_found += 1
            if sum((k - 1) * m for k, m in mu.items() if k >= 2) != Fr(6, 7): slack_ok = False
check(f"covering-slack identity: every covering x found ({cov_found} across instances/q<120) has "
      "overlap excess EXACTLY 6/7", cov_found > 0 and slack_ok)

# ---- B_D exactness criterion + where the deficit lives ----
print()
print("max multiplicity + B5 deficit by x-order (tight AP {1..13}):")
print("  q | x=1/q | maxC | mu_0 | B5 deficit = sum C(k-1,5) mu_k")
for q in [2, 3, 5, 7, 13, 14, 28, 49, 91, 101, 187, 401, 1009]:
    x = Fr(1, q)
    mu = spectrum(AP13, x)
    maxC = max(mu)
    deficit = sum(comb(k - 1, 5) * m for k, m in mu.items() if k >= 6)
    print(f"  {q:4d} |  1/{q:<4d} | {maxC:2d} | {str(mu.get(0, Fr(0))):>8s} | {str(deficit):>10s}"
          f"{'   <- B5 EXACT' if deficit == 0 else ''}")
# THE MID-BAND LAW (prediction was backwards — this is the discovery): for the AP cluster at
# x = 1/q the max multiplicity is the window count min(13, floor(q/7) (+1 if 7∤q)); the lonely
# regime is LARGE q (single cluster, arcs stack 13-deep, exactly the far-tail/Lemma-A band);
# perfect minimal covers live in the MID-BAND q = 7..14 (cf. THM-667 mid-band realization).
ok_mid = True
for q in range(14, 201):
    mu = spectrum(AP13, Fr(1, q))
    pred = min(13, q // 7 + (1 if q % 7 else 0))
    if max(mu) != pred: ok_mid = False
check("MID-BAND LAW: maxC(AP, x=1/q) = min(13, floor(q/7)+[7∤q]) for q = 14..200 "
      "(=> clustered/lonely for q > 91; B5-exact zone = maxC <= 5 <=> q <= 35 + covering band)",
      ok_mid)

# UNIQUENESS OF THE MINIMAL COVERING SPECTRUM: covering + maxC = 2 forces mu = {1:1/7, 2:6/7}
min_cov, ok_uni = 0, True
for name, E in INSTANCES.items():
    for q in range(7, 120):
        mu = spectrum(E, Fr(1, q))
        if mu.get(0, Fr(0)) == 0 and max(mu) == 2:
            min_cov += 1
            if mu != {1: Fr(1, 7), 2: Fr(6, 7)}: ok_uni = False
check(f"RIGIDITY: every minimal-multiplicity covering x ({min_cov} found) has the UNIQUE spectrum "
      "(mu_1, mu_2) = (1/7, 6/7) — the LRC twin of the K(A5) multiplicity-<=2 Kakeya witness",
      min_cov > 0 and ok_uni)

# high-multiplicity time location: at x = 1/q the champion time carries the cluster structure
mu7 = spectrum(AP13, Fr(1, 7))
print()
print("x = 1/7 spectrum (tight AP):", {k: str(v) for k, v in sorted(mu7.items())})
mu1009 = spectrum(AP13, Fr(377, 1009))
print("x = 377/1009 spectrum      :", {k: str(v) for k, v in sorted(mu1009.items())})
bino = {k: comb(13, k) * Fr(1, 7) ** k * Fr(6, 7) ** (13 - k) for k in range(14)}
tv_low = sum(abs(mu7.get(k, Fr(0)) - bino[k]) for k in range(14)) / 2
tv_high = sum(abs(mu1009.get(k, Fr(0)) - bino[k]) for k in range(14)) / 2
check(f"generic x is near-binomial, resonant x is far (TV to Binom(13,1/7): {float(tv_high):.3f} "
      f"vs {float(tv_low):.3f})", tv_high < tv_low)

# ---- the pair moment = second-moment vein ----
# S_2 = sum over pairs of overlap. For an AP at generic x the points are EQUALLY SPACED
# (maximally correlated) so S_2 sits far from the independent value C(13,2)/49; a random
# (unstructured) cluster at generic x approaches it. Structure = correlation = the vein.
S2_ap = S_d_moment(mu1009, 2)
mu_rand = spectrum(RAND13, Fr(377, 1009))
S2_rand = S_d_moment(mu_rand, 2)
indep = Fr(comb(13, 2), 49)
check(f"second-moment vein: random cluster S_2 = {float(S2_rand):.3f} ~ independent "
      f"{float(indep):.3f}; the AP's S_2 = {float(S2_ap):.3f} deviates (equal spacing = maximal "
      "pair correlation — THM-863's rho-floor is the pointwise floor on these summands)",
      abs(S2_rand - indep) < Fr(1, 4) and abs(S2_ap - indep) > abs(S2_rand - indep))

print()
fails = [n for n, c in OK if not c]
print(f"=== {len(OK)} checks, {len(OK) - len(fails)} passed, {len(fails)} failed ===")
for f in fails: print("FAILED:", f)

# ================= (W) THE PAIR-SPREAD WITNESS CRITERION =================
# Covering x forces S_2 >= 6/7 (convexity: minimize sum C(k,2) mu_k under R0/R1, mu_0=0:
# support {1,2} => 6/7).  Contrapositive: S_2(x) < 6/7 ==> mu_0(x) > 0 ==> maxgap > 1/7:
# a WITNESS from pair data alone: S_2(x) = sum_{i<j} max(0, 1/7 - ||(e_i-e_j) x||).
def S2_pairs(E, x):
    tot = Fr(0)
    for i in range(len(E)):
        for j in range(i + 1, len(E)):
            d = ((E[j] - E[i]) * x) % 1
            dist = min(d, 1 - d)
            if dist < W: tot += W - dist
    return tot
ok_w, fires = True, {}
for name, E in INSTANCES.items():
    fires[name] = []
    best = None
    for q in range(7, 400):
        for a in range(1, min(q, 8)):
            x = Fr(a, q)
            s2 = S2_pairs(E, x)
            if best is None or s2 < best[0]: best = (s2, x)
            if s2 < Fr(6, 7):
                fires[name].append(x)
                if spectrum(E, x).get(0, Fr(0)) == 0: ok_w = False   # must be a witness
    fires[name] = (len(fires[name]), best)
check("(W) pair-spread witness criterion: S_2(x) < 6/7 ==> uncovered (verified everywhere it "
      "fires); S_2 computed from the DIFFERENCE SET alone", ok_w)
for name, (cnt, best) in fires.items():
    print(f"   {name:>16s}: fires at {cnt:4d} sampled x; min S_2 = {float(best[0]):.4f} at x = {best[1]}"
          + ("   <- NEVER fires (saturates the criterion)" if cnt == 0 else ""))
print()
fails = [n for n, c in OK if not c]
print(f"=== FINAL {len(OK)} checks, {len(OK) - len(fails)} passed ===")

# ================= (FT) THE FEJES TOTH FLOOR: why (W) is vacuous =================
# S_2 as pure configuration energy: 13 arbitrary points, kernel tent = max(0, 1/7 - dist).
# Fejes Toth (convex decreasing kernel => regular minimizes): min = 13*(1/7 - 1/13) = 6/7.
# => S_2 >= 6/7 ALWAYS: the pair criterion can never fire; pair data cannot certify
# loneliness (the precise second-moment wall). Minimal certifying moment order = 3.
rngc = random.Random(7)
ok_ft = True
regular = [Fr(i, 13) for i in range(13)]
def S2_config(pos):
    tot = Fr(0)
    for i in range(len(pos)):
        for j in range(i + 1, len(pos)):
            d = (pos[j] - pos[i]) % 1
            dist = min(d, 1 - d)
            if dist < W: tot += W - dist
    return tot
assert S2_config(regular) == Fr(6, 7)
for _ in range(4000):
    pos = sorted(Fr(rngc.randrange(1, 10007), 10007) for _ in range(13))
    if S2_config(pos) < Fr(6, 7): ok_ft = False
check("(FT) Fejes-Toth floor: S_2 >= 6/7 for ALL 13-point configurations (4000 random exact "
      "trials; regular 13-gon attains 6/7 exactly) => the pair criterion is VACUOUS: the "
      "second-moment wall, made pointwise and exact", ok_ft)
# covering refinement: S_2 - excess = sum C(k-1,2) mu_k  (the >=3-fold mass meter)
ok_meter = True
for name, E in INSTANCES.items():
    for q in (7, 13, 14, 20, 50):
        mu = spectrum(E, Fr(1, q))
        s2 = S_d_moment(mu, 2)
        excess = sum((k - 1) * m for k, m in mu.items() if k >= 2)
        if s2 - excess != sum(comb(k - 1, 2) * m for k, m in mu.items() if k >= 3): ok_meter = False
check("(M3) S_2 - overlap excess = sum C(k-1,2) mu_k: for covering x, S_2 - 6/7 measures the "
      ">=3-fold mass exactly (distance from rigidity, configurationally AND spectrally)", ok_meter)
print()
print(f"=== GRAND TOTAL {len(OK)} checks, {sum(1 for _,c in OK if c)} passed ===")
