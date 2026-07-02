#!/usr/bin/env python3
"""
klein-2026-07-01-S88 -- HYP-3841 (overtaking-defect sieve + tangent-ladder certificate)
                        + HYP-3843 (final-window identity, (c,E) at n=20, deep-witness n=9,10)

A. THE SIEVE (HYP-3841). Lambda_S(r) is piecewise linear; kinks: gap-death (convex,
   slope jumps UP) at r = d/(v+w); overtaking (concave, slope can jump DOWN) at
   r = d/(w-v). Exposed concave defect K(S; [r0,r1]) = total downward slope jump.
   TANGENT-LADDER CERTIFICATE (THM-592 ladder, made per-set):
       Lambda(r) >= Lambda(r0) + (slope(r0+) - K) (r - r0)   for r in [r0, r1].
   If cert(1/14) >= 0 the certificate proves M(S) >= 1/14 from the single anchor r0.
   Question: does the certificate succeed on COVERING sets (which always contain a
   multiple of 14, so all shallow witnesses are dead and anchoring must be sub-critical)?
   And is K structurally small on [1/16, 1/14] (overtaking there needs w - v in
   [14d, 16d] -- a narrow difference band)?

B. FINAL-WINDOW IDENTITY (HYP-3843): Lambda_AP == Lambda_GW EXACTLY on [r*, 1/14];
   compute r* exactly (last radius where the profiles differ).

C. (c, E) at n=20: does E separate {AP_19, GW_19} the way it does at n=14?

D. Deep-witness emptiness at n=9, 10: exhaustive tight census; any tight set with a
   multiple of n? (extends the n<=8 emptiness observation)
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, random

# ---------------- exact machinery ----------------
def Lambda(S, r):
    ivs = []
    for v in S:
        rv = r / v
        for a in range(v + 1):
            c = F(a, v)
            lo, hi = c - rv, c + rv
            if hi <= 0 or lo >= 1:
                continue
            ivs.append((max(lo, F(0)), min(hi, F(1))))
    ivs.sort()
    tot = F(0); cl = ch = None
    for lo, hi in ivs:
        if ch is None: cl, ch = lo, hi
        elif lo <= ch: ch = max(ch, hi)
        else: tot += ch - cl; cl, ch = lo, hi
    if ch is not None: tot += ch - cl
    return 1 - tot

def m_at(S, t):
    best = None
    for v in S:
        x = (v * t) % 1
        d = min(x, 1 - x)
        if best is None or d < best: best = d
    return best

def M_exact(S):
    qs = set()
    for u, v in itertools.combinations(S, 2):
        qs.add(u + v); qs.add(abs(u - v))
    for v in S: qs.add(2 * v)
    qs.discard(0)
    best, wit = F(0), []
    for q in qs:
        for c in range(q + 1):
            t = F(c, q); m = m_at(S, t)
            if m > best: best, wit = m, [t]
            elif m == best: wit.append(t)
    return best, sorted(set(wit))

def breakpoints_in(S, rlo, rhi):
    """Candidate kink radii of Lambda_S in [rlo, rhi]: d/(v+w) and d/(w-v)."""
    bps = set()
    for u, v in itertools.combinations(sorted(S), 2):
        for q in (u + v, v - u):
            if q <= 0: continue
            d0 = int(rlo * q); d1 = int(rhi * q) + 1
            for d in range(max(d0, 1), d1 + 1):
                r = F(d, q)
                if rlo <= r <= rhi: bps.add(r)
    for v in S:  # same-speed tent peaks: gaps between adjacent same-v fractions die at r=1/2v...
        q = 2 * v
        d0 = int(rlo * q); d1 = int(rhi * q) + 1
        for d in range(max(d0, 1), d1 + 1):
            r = F(d, q)
            if rlo <= r <= rhi: bps.add(r)
    return sorted(bps)

def piecewise_profile(S, rlo, rhi):
    """Exact piecewise-linear profile of Lambda_S on [rlo, rhi]:
    returns list of (r_i, Lambda(r_i)) at all kinks + endpoints, with linearity verified."""
    bps = [rlo] + [b for b in breakpoints_in(S, rlo, rhi) if rlo < b < rhi] + [rhi]
    pts = [(r, Lambda(S, r)) for r in bps]
    # verify linearity between consecutive breakpoints (midpoint test)
    segs = []
    for (r0, L0), (r1, L1) in zip(pts, pts[1:]):
        rm = (r0 + r1) / 2
        Lm = Lambda(S, rm)
        assert 2 * Lm == L0 + L1, f"kink inside ({r0},{r1}) of {S}"
        slope = (L1 - L0) / (r1 - r0) if r1 != r0 else None
        segs.append((r0, r1, slope))
    return pts, segs

def defect_and_cert(S, r0, r1):
    """Exposed concave defect K on [r0, r1] and the tangent-ladder certificate at r1."""
    pts, segs = piecewise_profile(S, r0, r1)
    slopes = [s for (_, _, s) in segs if s is not None]
    K = F(0); down_kinks = []
    for s_prev, s_next, seg in zip(slopes, slopes[1:], segs[1:]):
        if s_next < s_prev:
            K += (s_prev - s_next)
            down_kinks.append((seg[0], s_prev - s_next))
    L0 = pts[0][1]
    s0 = slopes[0]
    cert = L0 + (s0 - K) * (r1 - r0)
    actual = pts[-1][1]
    return dict(L0=L0, s0=s0, K=K, cert=cert, actual=actual,
                nkinks=len(slopes) - 1, ndown=len(down_kinks), down=down_kinks[:4])

# ---------------- A. the sieve on covering sets ----------------
print("=" * 92)
print("A. OVERTAKING-DEFECT SIEVE (HYP-3841): tangent-ladder certificates, anchor r0=1/16 -> 1/14")
print("=" * 92)

def is_covering(S, n=14):
    return all(any(v % q == 0 for v in S) for q in range(2, n + 1))

random.seed(88)
def random_covering():
    """13 speeds containing a multiple of every q in 2..14."""
    while True:
        # cover 8..14 with dedicated multiples, small q usually free
        S = set()
        for q in [14, 13, 11, 9, 8]:
            S.add(q * random.randint(1, 6))
        while len(S) < 13:
            S.add(random.randint(1, 120))
        S = sorted(S)
        if len(S) == 13 and is_covering(S) and reduce(gcd, S) == 1:
            return S

CONSTR = list(range(1, 13)) + [182]
tests = [("CONSTR {1..12,182}", CONSTR)]
for i in range(8):
    tests.append((f"RANDCOV{i}", random_covering()))
# adversarial: covering with clustered large elements
tests.append(("CLUSTER", [2, 3, 4, 5, 6, 7, 9, 11, 13, 14, 100, 101, 102]))
tests.append(("TWO-SCALE", [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 28, 182]))

r0, r1 = F(1, 16), F(1, 14)
ok = 0
for name, S in tests:
    cov = is_covering(S)
    d = defect_and_cert(S, r0, r1)
    Mv, _ = M_exact(S) if max(S) <= 200 else (None, None)
    succ = d["cert"] >= 0
    ok += succ
    print(f"  {name:22s} cov={cov}  L(1/16)={float(d['L0']):.4f}  s0={float(d['s0']):+.3f}  "
          f"K={float(d['K']):.3f} (#down={d['ndown']}/{d['nkinks']})  cert(1/14)={float(d['cert']):+.5f}  "
          f"actual={float(d['actual']):.5f}  CERT-OK={succ}" + (f"  M={Mv}" if Mv else ""))
print(f"\n  certificate success: {ok}/{len(tests)}")
print("  overtaking band check: kinks at d/(w-v) in [1/16,1/14] need (w-v)/d in [14,16] -- differences")
print("  in a narrow band; K small iff few pairs have w-v in {14d..16d} with aligned fractions.")

# ---------------- B. final-window identity ----------------
print("\n" + "=" * 92)
print("B. FINAL-WINDOW IDENTITY (HYP-3843): Lambda_AP == Lambda_GW exactly on [r*, 1/14]?")
print("=" * 92)
AP = list(range(1, 14)); GW = list(range(1, 12)) + [13, 24]
# candidate distinguishing radii: union of both breakpoint sets on [1/20, 1/14]
cands = sorted(set(breakpoints_in(AP, F(1, 20), F(1, 14)) + breakpoints_in(GW, F(1, 20), F(1, 14))))
diffs = []
for r in cands + [F(1, 20), F(1, 14)]:
    dAP, dGW = Lambda(AP, r), Lambda(GW, r)
    if dAP != dGW: diffs.append((r, dGW - dAP))
diffs.sort()
if diffs:
    rlast = diffs[-1][0]
    print(f"  profiles DIFFER at {len(diffs)} candidate radii; LAST difference at r = {rlast} = {float(rlast):.6f}")
    # r* = first radius after rlast where they agree onward; verify identity on (rlast, 1/14]
    after = [r for r in cands if r > rlast]
    probe = sorted(set([ (rlast + F(1,14))/2 ] + after))
    all_eq = all(Lambda(AP, r) == Lambda(GW, r) for r in probe)
    print(f"  identity on ({rlast}, 1/14]: {all_eq}  => r* = {rlast} (exclusive)")
    print(f"  window width 1/14 - r* = {float(F(1,14) - rlast):.6f}; midpoints check exact-equal: "
          f"{Lambda(AP, (rlast + F(1,14))/2) == Lambda(GW, (rlast + F(1,14))/2)}")
else:
    print("  profiles identical on the whole [1/20, 1/14] ?!")

# ---------------- C. (c, E) at n=20 ----------------
print("\n" + "=" * 92)
print("C. (c, E) AT n=20 (HYP-3843): does mean loneliness separate {AP_19, GW_19}?")
print("=" * 92)
def candidate_times(S):
    qs = set()
    for u, v in itertools.combinations(S, 2):
        qs.add(u + v); qs.add(abs(u - v))
    for v in S: qs.add(2 * v)
    qs.discard(0)
    c = set()
    for q in qs:
        for a in range(q + 1): c.add(F(a, q))
    return c

def mean_loneliness(S):
    bps = sorted(candidate_times(S) | {F(0), F(1)})
    tot = F(0); prev = bps[0]; vprev = m_at(S, prev)
    for t in bps[1:]:
        vt = m_at(S, t)
        mid = (prev + t) / 2; vm = m_at(S, mid)
        assert vm * 2 == vprev + vt, "kink inside piece"
        tot += (t - prev) * (vprev + vt) / 2
        prev, vprev = t, vt
    return tot

AP20 = list(range(1, 20)); GW20 = list(range(1, 18)) + [19, 36]
eA, eG = mean_loneliness(AP20), mean_loneliness(GW20)
print(f"  E(AP_19) = {float(eA):.8f}  E(GW_19) = {float(eG):.8f}  GW<AP: {eG < eA}  (exact diff {float(eG-eA):+.3e})")
print(f"  (collapse rates equal by HYP-3835 family universality; E separates iff GW<AP strictly)")

# ---------------- D. deep-witness emptiness at n=9,10 ----------------
print("\n" + "=" * 92)
print("D. DEEP-WITNESS EMPTINESS (HYP-3843): tight census n=9 (v<=20), n=10 (v<=18); multiples of n?")
print("=" * 92)
def lam_float(S, r):
    ivs = []
    for v in S:
        rv = r / v; step = 1.0 / v
        for a in range(v + 1):
            c = a * step
            lo, hi = c - rv, c + rv
            if hi <= 0 or lo >= 1: continue
            ivs.append((max(lo, 0.0), min(hi, 1.0)))
    ivs.sort()
    tot, cl, ch = 0.0, None, None
    for lo, hi in ivs:
        if ch is None: cl, ch = lo, hi
        elif lo <= ch: ch = max(ch, hi)
        else: tot += ch - cl; cl, ch = lo, hi
    if ch is not None: tot += ch - cl
    return 1.0 - tot

for n, vmax in [(9, 20), (10, 18)]:
    k = n - 1; thr = 1.0 / n; thrF = F(1, n); found = []
    for c in itertools.combinations(range(1, vmax + 1), k):
        if reduce(gcd, c) != 1: continue
        if lam_float(c, thr * (1 + 1e-9)) > 1e-12: continue
        Mv, wits = M_exact(list(c))
        if Mv == thrF: found.append((c, wits))
    print(f"  n={n}: tight sets (v<={vmax}): {len(found)}")
    for S, wits in found:
        multn = [v for v in S if v % n == 0]
        print(f"    {S}  #wit={len(wits)}  multiples-of-{n}: {multn if multn else 'NONE'}")

print("\nDONE.")

# ---------------- A'. chained ladder rescue for the two failures ----------------
print("\n" + "=" * 92)
print("A'. CHAINED TANGENT-LADDER: anchors 1/16 -> 1/15 -> 1/14 for the two single-anchor failures")
print("=" * 92)
for name, S in [("CONSTR {1..12,182}", CONSTR), ("TWO-SCALE", [1,2,3,4,5,6,7,8,9,11,13,28,182])]:
    d1 = defect_and_cert(S, F(1,16), F(1,15))
    d2 = defect_and_cert(S, F(1,15), F(1,14))
    # pure certificate: step-2 starting value = step-1 certified bound (not the actual)
    cert_chain = d1["cert"] + (d2["s0"] - d2["K"]) * (F(1,14) - F(1,15))
    print(f"  {name:22s} step1 cert(1/15)={float(d1['cert']):+.5f} (actual {float(d1['actual']):.5f}, K1={float(d1['K']):.3f})")
    print(f"  {'':22s} step2 from cert1: cert(1/14)={float(cert_chain):+.5f} (K2={float(d2['K']):.3f}, s(1/15+)={float(d2['s0']):+.3f})  CHAIN-OK={cert_chain >= 0}")
