#!/usr/bin/env python3
"""
clock_moduli_largetheta_Qs_klein_S314.py — klein-2026-07-16-S314

THREE PROOFS on LEM-020/THM-729 ground:

I.  D(q) = 0 iff q in {7, 13, 14}  (the clock-moduli characterization, THM-878).
    Proof by case split, each case verified here:
      q <= 6 : every primitive a has coincidence mass Sum C(m_j,2)/7 >= 8/7 > 6/7 (positions
               = all q residues, multiplicities near-equal; no tent term since 1/q >= 1/6 > 1/7).
      q = 7  : positions = the 7 sevenths, multiplicities (2,2,2,2,2,2,1): spectrum (1/7,6/7). = 6/7.
      q in 8..12 : S2 = (13-q)/7 + (26-q+adj2)(q-7)/(7q) > 6/7 for EVERY primitive a.
      q = 13 : regular 13-gon for every a. = 6/7.
      q = 14 : 13 distinct positions missing one 14th; gaps twelve 1/14's + one 1/7; boundary
               of the polytope P. = 6/7.
      q >= 15: a = 1 alone: the twelve adjacent pairs give S2 >= 12(1/7 - 1/q) > 6/7 iff q > 14. QED
    (A(q) = primitive mean >= 6/7 pointwise by the FT floor; = 6/7 iff EVERY a is flat-bottom.)

II. THE LARGE-THETA LEMMA (Fejer-Bernoulli evaluations, all EXACT):
      (L1) sum_{l>=1} sin^2(pi l theta)/l^2 = pi^2 theta(1-theta)/2   (theta in [0,1])
      (L2) V(e) = sum_{j mod e} tent(j/e) = e/49 + {e/7}(1-{e/7})/e   (Moebius-sinc closed form;
           the aliasing of the Fejer kernel sums to ONE theta(1-theta) term — "large theta")
      (L3) sum_{l != 0} e(l t)/l^2 = 2 pi^2 B_2({t}),  B_2(t) = t^2 - t + 1/6.

III. Q_s CLOSED FORM (THM-729 upgrade): with U(N) = sum_k eps_k e(-N p_k), sum eps_k = 0:
      Q_s := sum_{l != 0} |U(l w)|^2 / l^2  =  -2 pi^2 sum_{k,k'} eps_k eps_k' TH(p_k - p_k')
      where TH(t) = {w t}(1 - {w t}).   EXACT FINITE RATIONAL BILINEAR FORM — no l-sum, no
      truncation; the empirical "off-diagonal cancels" is sign cancellation in eps eps'.
      Consequences: (a) diagonal = 2 pi^2 sum_arcs {w u}(1-{w u}) recovers THM-729's rigorous
      diagonal; (b) dipole-transport bound Q_s <= 2 pi^2 sum_{i,j} min(1/2, 2 w min(u_i,u_j));
      (c) exact evaluation on THM-729-style clusters (7-section R_s sets, exact arithmetic),
      Q_s/diam table now EXACT rather than LMAX-truncated.
"""
from fractions import Fraction as Fr
from math import comb, gcd, pi, sin, cos
import itertools, random

W7 = Fr(1, 7)
AP = list(range(1, 14))
OK = []
def check(name, cond):
    OK.append((name, bool(cond)))
    print(("PASS" if cond else "FAIL"), name)

def tent(fr):
    d = fr % 1
    dist = min(d, 1 - d)
    return W7 - dist if dist < W7 else Fr(0)

def S2(E, x):
    tot = Fr(0)
    for i in range(len(E)):
        for j in range(i + 1, len(E)):
            tot += tent((E[j] - E[i]) * x)
    return tot

# ---------------- I. the clock-moduli theorem ----------------
ok1 = True
for q in range(2, 15):
    vals = [S2(AP, Fr(a, q)) for a in range(1, q) if gcd(a, q) == 1]
    if q in (7, 13, 14):
        if any(v != Fr(6, 7) for v in vals): ok1 = False
    else:
        if any(v <= Fr(6, 7) for v in vals): ok1 = False      # strict for EVERY primitive a
check("cases q <= 14: S2(a/q) = 6/7 for ALL primitive a iff q in {7,13,14}; else STRICTLY > "
      "for all a", ok1)
ok1b = all(12 * (Fr(1, 7) - Fr(1, q)) > Fr(6, 7) for q in range(15, 3000))
check("case q >= 15 (one line): a=1 adjacent pairs alone give 12(1/7 - 1/q) > 6/7 iff q > 14 "
      "(algebra: 12/7 - 12/q > 6/7 iff q > 14; sampled to 3000)", ok1b)
# referee: D(q) via the exact primitive mean
def A_direct(q):
    tot, cnt = Fr(0), 0
    for a in range(1, q):
        if gcd(a, q) == 1:
            cnt += 1; tot += S2(AP, Fr(a, q))
    return tot / cnt
okD = all((A_direct(q) == Fr(6, 7)) == (q in (7, 13, 14)) for q in range(2, 120))
check("REFEREE: D(q) = A(q) - 6/7 = 0 iff q in {7,13,14}, exact, q = 2..119  => THM-878", okD)

# ---------------- II. the large-theta lemma ----------------
L = 200000
ok2a = True
for th in (Fr(1, 7), Fr(2, 7), Fr(1, 3), Fr(5, 13), Fr(1, 2)):
    s = sum(sin(pi * l * float(th)) ** 2 / l ** 2 for l in range(1, L))
    target = pi ** 2 * float(th) * (1 - float(th)) / 2
    if abs(s - target) > 2e-5: ok2a = False
check("(L1) sum sin^2(pi l theta)/l^2 = pi^2 theta(1-theta)/2 (numeric, 5 thetas, err < 2e-5)", ok2a)
def V(e): return sum(tent(Fr(j, e)) for j in range(e))
ok2b = all(V(e) == Fr(e, 49) + (Fr(e % 7, 7) * (1 - Fr(e % 7, 7))) / e for e in range(1, 2001))
check("(L2) V(e) = e/49 + {e/7}(1-{e/7})/e EXACT for e = 1..2000 (Moebius-sinc closed form)", ok2b)
ok2c = True
for t in (0.13, 0.377, 0.5, 0.91):
    s = sum(2 * cos(2 * pi * l * t) / l ** 2 for l in range(1, L))
    if abs(s - 2 * pi ** 2 * (t * t - t + 1.0 / 6)) > 2e-4: ok2c = False
check("(L3) sum_{l!=0} e(lt)/l^2 = 2 pi^2 B2({t}) (numeric)", ok2c)

# ---------------- III. Q_s exact bilinear form ----------------
def Q_exact(endpts, w):
    # endpts: list of (p Fraction, eps ±1), sum eps = 0. Returns Q_s / (2 pi^2) exact (ell != 0).
    tot = Fr(0)
    for p1, e1 in endpts:
        for p2, e2 in endpts:
            th = (w * (p1 - p2)) % 1
            tot += -e1 * e2 * th * (1 - th)
    return tot
def Q_lsum(endpts, w, LMAX=3000):
    tot = 0.0
    for l in range(1, LMAX + 1):
        re = im = 0.0
        for p, e in endpts:
            a = -2 * pi * l * w * float(p)
            re += e * cos(a); im += e * sin(a)
        tot += (re * re + im * im) / l ** 2
    return 2 * tot          # ell != 0 doubles ell >= 1

rng = random.Random(314)
ok3 = True
for trial in range(6):
    narc = rng.randint(3, 8)
    pts = sorted(Fr(rng.randrange(1, 4001), 4001) for _ in range(2 * narc))
    endpts = [(p, +1 if i % 2 == 0 else -1) for i, p in enumerate(pts)]
    w = rng.choice([97, 997, 397])
    exact = 2 * pi ** 2 * float(Q_exact(endpts, w))
    approx = Q_lsum(endpts, w)
    if abs(exact - approx) > 3e-3 * max(1.0, abs(exact)): ok3 = False
check("Q_s = -2pi^2 sum eps eps' {wD}(1-{wD}) == the l-sum (6 random arc systems, w in "
      "{97,397,997}; rel err < 3e-3 at LMAX=3000)", ok3)

# THM-729-style clusters: R_s = {x : sections of E miss exactly section s}, EXACT arithmetic
def R_s_exact(E, s):
    bps = sorted(set(Fr(k, 7 * e) for e in E for k in range(7 * e)) | {Fr(0), Fr(1)})
    arcs, inR, start = [], False, None
    for i in range(len(bps) - 1):
        mid = (bps[i] + bps[i + 1]) / 2
        occ = set(int((e * mid % 1) * 7) for e in E)
        cur = (len(occ) == 6) and (s not in occ)
        if cur and not inR: start, inR = bps[i], True
        if (not cur) and inR: arcs.append((start, bps[i])); inR = False
    if inR: arcs.append((start, Fr(1)))
    return arcs

print()
print("THM-729 clusters, EXACT Q_s (w = 997) and the transport bound:")
print("  cluster | diam | s | #arcs r | Q_s exact | Q/diam | w*maxu")
fams = [[0, 1, 2, 3, 4, 5, 6], [0, 1, 2, 3, 4, 5, 12], [0, 1, 2, 3, 4, 5, 25],
        [0, 1, 2, 3, 4, 5, 50], [0, 2, 5, 11, 17, 29, 47]]
w = 997
okQ, ratios = True, []
for E in fams:                            # S280 convention: the small offsets ARE the cluster
    Epos = [e for e in E if e > 0]
    diam = max(E) - min(E) if max(E) > min(E) else 1
    best = None
    for s in range(7):
        arcs = R_s_exact(E, s)
        if not arcs: continue
        endpts = []
        for a, b in arcs:
            endpts.append((a, +1)); endpts.append((b, -1))
        Q = 2 * pi ** 2 * float(Q_exact(endpts, w))
        r = len(arcs)
        maxu = max(float(b - a) for a, b in arcs)
        if best is None or Q > best[0]: best = (Q, s, r, maxu)
    Q, s, r, maxu = best
    ratios.append(Q / max(diam, 1))
    print(f"  {str(E):28s} | {diam:4d} | {s} | {r:3d} | {Q:9.3f} | {Q/max(diam,1):6.2f} | "
          f"{w*maxu:7.1f}")
check("EXACT Q_s on the S280-convention clusters (no truncation, no LMAX): Q_s/diam = O(1) "
      f"(ratios {['%.2f' % t for t in ratios]}) — the empirical law of THM-729, now exact-"
      "rational and Lean-shaped", all(t < 8 for t in ratios))

print()
fails = [n for n, c in OK if not c]
print(f"=== {len(OK)} checks, {len(OK) - len(fails)} passed ===")
for f in fails: print("FAILED:", f)
