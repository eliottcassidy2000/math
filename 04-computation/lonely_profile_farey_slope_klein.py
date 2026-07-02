#!/usr/bin/env python3
"""
klein-2026-07-01-S87 -- HYP-3834 / HYP-3835 / HYP-3836

THE LONELY DISTRIBUTION FUNCTION Lambda_S(r) = meas{t : m_S(t) >= r},
where m_S(t) = min_{v in S} ||vt|| is the (codex-S179) lonely profile.

Parts:
  A. exact machinery (Fractions throughout): Lambda_S(r), M(S), witnesses, collapse rate
  B. n=14: AP {1..13} = Farey gap-sum, collapse constant 1666/6435 EXACT;
     GW {1..11,13,24} universality test; dilates; construction {1..12,182};
     the MISTAKE-075/THM-523 perturbation {1..11,13,36} -> Lambda(1/14) =? 1/1260
  C. general n: c(AP_n) = (2/n) H^x(n) for n=4..16; GW family universality at n=8, 20
  D. profile crossing (spread vs tight) + small-n exhaustive lonely envelope
  E. mean loneliness E(S) = int m_S = int Lambda (layer-cake); Farey closed form; AP-min test

All measures/values exact rationals unless noted.
"""
from fractions import Fraction as F
from math import gcd
import itertools, random

# ---------------------------------------------------------------- A. machinery

def danger_intervals(S, r):
    """All danger intervals (a/v - r/v, a/v + r/v) clipped to [0,1], as (lo,hi) Fractions."""
    ivs = []
    for v in S:
        rv = r / v
        for a in range(v + 1):
            c = F(a, v)
            lo, hi = c - rv, c + rv
            if hi <= 0 or lo >= 1:
                continue
            ivs.append((max(lo, F(0)), min(hi, F(1))))
    return ivs

def union_measure(ivs):
    ivs = sorted(ivs)
    tot = F(0)
    cur_lo, cur_hi = None, None
    for lo, hi in ivs:
        if cur_hi is None:
            cur_lo, cur_hi = lo, hi
        elif lo <= cur_hi:
            cur_hi = max(cur_hi, hi)
        else:
            tot += cur_hi - cur_lo
            cur_lo, cur_hi = lo, hi
    if cur_hi is not None:
        tot += cur_hi - cur_lo
    return tot

def Lambda(S, r):
    """Exact lonely distribution function: meas{t: min_v ||vt|| >= r}."""
    return 1 - union_measure(danger_intervals(S, r))

def m_at(S, t):
    """Exact m_S(t) = min_v ||vt||."""
    best = None
    for v in S:
        x = (v * t) % 1
        d = min(x, 1 - x)
        if best is None or d < best:
            best = d
    return best

def candidate_times(S):
    """Breakpoint candidates for max/kinks of m_S: c/q, q in {u+v, |u-v|, 2v}."""
    qs = set()
    for u, v in itertools.combinations(S, 2):
        qs.add(u + v)
        if u != v:
            qs.add(abs(u - v))
    for v in S:
        qs.add(2 * v)
    qs.discard(0)
    cands = set()
    for q in qs:
        for c in range(q + 1):
            cands.add(F(c, q))
    return cands

def M_exact(S):
    """Exact M(S) = max_t min_v ||vt|| and the exact set of maximizers."""
    best, wit = F(0), []
    for t in candidate_times(S):
        m = m_at(S, t)
        if m > best:
            best, wit = m, [t]
        elif m == best:
            wit.append(t)
    return best, sorted(set(wit))

def binding_speeds(S, t, val):
    return [v for v in S if m_at([v], t) == val]

def collapse_rate(S, thresh):
    """Exact one-sided slope of Lambda at r->thresh^- (assumes M(S)=thresh, linear regime).
    Cross-checked at two epsilons; returns (slope Fraction, consistent bool)."""
    out = []
    for k in (10**4, 2 * 10**4):
        r = thresh * (1 - F(1, k))
        L = Lambda(S, r)
        out.append(L / (thresh * F(1, k) * (1 / thresh)))  # L / (1 - r/thresh) scaled below
    # careful: 1 - r/thresh = 1/k; but the natural variable is (1 - n r) with n=1/thresh
    n = 1 / thresh
    s1 = Lambda(S, thresh * (1 - F(1, 10**4))) / (F(1, 10**4))
    s2 = Lambda(S, thresh * (1 - F(1, 2 * 10**4))) / (F(1, 2 * 10**4))
    return s1, s1 == s2

def slope_from_witnesses(S, thresh):
    """Predicted slope: sum over maximizers of (1/n)*(1/v+ + 1/v-), n = 1/thresh.
    v+/v- = binding speeds with derivative +/- at the witness (for generic double binding)."""
    Mv, wits = M_exact(S)
    assert Mv == thresh, (Mv, thresh)
    tot = F(0)
    for t in wits:
        bs = binding_speeds(S, t, thresh)
        # each binding speed contributes 1/(n*v) on one side (up or down as t moves)
        for v in bs:
            tot += thresh * F(1, v)  # thresh = 1/n
        if len(bs) == 1:
            tot += thresh * F(1, bs[0])  # single tent peak: both sides same speed
    return tot, wits

# ------------------------------------------------------------ B. n=14 verification

def farey_pairs(N):
    """Adjacent pairs in Farey_N via Stern-Brocot walk."""
    a, b, c, d = 0, 1, 1, N
    pairs = []
    while (a, b) != (1, 1):
        pairs.append(((a, b), (c, d)))
        k = (N + b) // d
        a, b, c, d = c, d, k * c - a, k * d - b
    return pairs

def farey_gap_sum(N, r):
    tot = F(0)
    for (a, q), (a2, q2) in farey_pairs(N):
        g = 1 - r * (q + q2)
        if g > 0:
            tot += g / (q * q2)
    return tot

def Hx(n):
    return sum(F(1, u) for u in range(1, n) if gcd(u, n) == 1)

print("=" * 78)
print("PART B: n=14 -- AP, GW, dilates, construction, perturbation")
print("=" * 78)

AP = list(range(1, 14))
GW = list(range(1, 12)) + [13, 24]
CONSTR = list(range(1, 13)) + [182]
PERT36 = list(range(1, 12)) + [13, 36]

n = 14
thr = F(1, n)

# B1: AP profile = Farey gap-sum (exact, several radii)
print("\nB1. AP {1..13}: Lambda vs Farey gap-sum (exact equality test)")
for r in [F(1, 100), F(1, 30), F(1, 20), F(1, 15), F(9999, 14 * 10**4)]:
    L, Fg = Lambda(AP, r), farey_gap_sum(13, r)
    print(f"  r={str(r):>14}  Lambda={str(L):>22}  Farey={str(Fg):>22}  equal={L == Fg}")

# B2: M and witnesses
for name, S in [("AP", AP), ("GW", GW), ("CONSTR", CONSTR), ("PERT36", PERT36)]:
    Mv, wits = M_exact(S)
    print(f"\nB2. {name:8s} M(S) = {Mv} = {float(Mv):.6f}   #maximizers={len(wits)}")
    if len(wits) <= 8:
        for t in wits:
            bs = binding_speeds(S, t, Mv)
            print(f"      witness t={t}   binding={bs}")

# B3: collapse rates (exact)
print("\nB3. collapse rates c(S) = lim Lambda/(1-14r)")
target = 2 * (F(1, 13) + F(1, 33) + F(1, 45))
print(f"  hand-derived target c(AP) = 2(1/13+1/33+1/45) = {target} = {float(target):.6f}")
print(f"  formula (2/n)*Hx(n)       = {F(2,14)*Hx(14)} = {float(F(2,14)*Hx(14)):.6f}")
for name, S in [("AP", AP), ("GW", GW), ("2*AP (dilate)", [2 * v for v in AP])]:
    s, ok = collapse_rate(S, thr)
    pred, wits = slope_from_witnesses(S, thr)
    print(f"  {name:14s} c = {s} = {float(s):.6f}   linear-regime={ok}   witness-formula={pred}  match={s == pred}")

# B4: construction + perturbation measures at the critical radius
print("\nB4. Lambda at r=1/14 (exact)")
for name, S in [("CONSTR {1..12,182}", CONSTR), ("PERT36 {1..11,13,36}", PERT36),
                ("AP", AP), ("GW", GW)]:
    L = Lambda(S, thr)
    print(f"  {name:22s} Lambda(1/14) = {L} = {float(L):.6f}")
print("  (THM-523 single-perturbation infimum is 1/1260 = {:.6f} -- PERT36 should match)".format(1/1260))

# ------------------------------------------------------------ C. general n
print("\n" + "=" * 78)
print("PART C: general n -- c(AP_n) = (2/n) Hx(n); GW-family universality")
print("=" * 78)
print("\nC1. AP_n = {1..n-1}")
for nn in range(4, 17):
    S = list(range(1, nn))
    Mv, wits = M_exact(S)
    s, ok = collapse_rate(S, F(1, nn))
    pred = F(2, nn) * Hx(nn)
    print(f"  n={nn:2d}  M={str(Mv):>6}  #wit={len(wits):2d} (phi={sum(1 for k in range(1,nn) if gcd(k,nn)==1)})"
          f"  c={str(s):>16} = {float(s):.6f}  (2/n)Hx={str(pred):>16}  match={s == pred}")

print("\nC2. GW family {1..m-2, m, 2(m-1)}, m=n-1=1 mod 6 -> tight at 1/n; universality c(GW_n)=c(AP_n)?")
for m in (7, 13, 19):
    nn = m + 1
    S = list(range(1, m - 1)) + [m, 2 * (m - 1)]
    Mv, wits = M_exact(S)
    tight = (Mv == F(1, nn))
    if tight:
        s, ok = collapse_rate(S, F(1, nn))
        pred = F(2, nn) * Hx(nn)
        print(f"  m={m:2d} (n={nn:2d})  S={S}  M={Mv}  tight={tight}  #wit={len(wits)}  "
              f"c={s}  =c(AP)? {s == pred}")
    else:
        print(f"  m={m:2d} (n={nn:2d})  S={S}  M={Mv}  tight={tight}  (NOT tight -- family condition fails?)")

# ------------------------------------------------------------ D. crossing + envelope
print("\n" + "=" * 78)
print("PART D: profile crossing + small-n lonely envelope")
print("=" * 78)

PRIMES = [17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67]
BLOCK100 = list(range(100, 113))
random.seed(87)
RANDOMS = []
while len(RANDOMS) < 3:
    S = sorted(random.sample(range(1, 400), 13))
    if gcd(*S) if len(S) == 1 else True:
        from functools import reduce
        if reduce(gcd, S) == 1:
            RANDOMS.append(S)

fams = [("AP", AP), ("GW", GW), ("PRIMES>13", PRIMES), ("BLOCK 100..112", BLOCK100),
        ("CONSTR", CONSTR)] + [(f"RAND{i}", S) for i, S in enumerate(RANDOMS)]
rgrid = [F(1, 1000), F(1, 100), F(1, 50), F(1, 28), F(1, 20), F(1, 16), F(1, 15),
         F(70, 1000), F(1, 14)]
print("\nD1. Lambda_S(r) table (floats for readability; exact inside)")
hdr = "r".rjust(10) + "".join(name.rjust(16) for name, _ in fams)
print(hdr)
for r in rgrid:
    row = f"{float(r):10.5f}"
    for name, S in fams:
        row += f"{float(Lambda(S, r)):16.6f}"
    print(row)
print("  union-bound floor 1-26r at r=1/1000:", float(1 - 26 * F(1, 1000)))

print("\nD2. small-n exhaustive lonely envelope: n=4 (3 speeds <= 15, primitive)")
from functools import reduce
sets4 = [list(c) for c in itertools.combinations(range(1, 16), 3) if reduce(gcd, c) == 1]
for r in [F(1, 100), F(1, 30), F(1, 12), F(1, 8), F(1, 5), F(9, 40), F(1, 4)]:
    best = min(sets4, key=lambda S: Lambda(S, r))
    bl = Lambda(best, r)
    ap = Lambda([1, 2, 3], r)
    print(f"  r={str(r):>7}  min Lambda={float(bl):.6f} at {best}   Lambda(AP)={float(ap):.6f}"
          f"   AP-is-min={bl == ap}")

print("\nD3. small-n exhaustive: n=5 (4 speeds <= 12, primitive)")
sets5 = [list(c) for c in itertools.combinations(range(1, 13), 4) if reduce(gcd, c) == 1]
for r in [F(1, 100), F(1, 30), F(1, 12), F(1, 7), F(1, 6), F(19, 100), F(1, 5)]:
    best = min(sets5, key=lambda S: Lambda(S, r))
    bl = Lambda(best, r)
    ap = Lambda([1, 2, 3, 4], r)
    print(f"  r={str(r):>7}  min Lambda={float(bl):.6f} at {best}   Lambda(AP)={float(ap):.6f}"
          f"   AP-is-min={bl == ap}")

# ------------------------------------------------------------ E. mean loneliness
print("\n" + "=" * 78)
print("PART E: mean loneliness E(S) = int_0^1 m_S dt = int_0^inf Lambda dr")
print("=" * 78)

def mean_loneliness(S):
    """Exact integral of the piecewise-linear m_S over [0,1]."""
    bps = sorted(candidate_times(S) | {F(0), F(1)})
    tot = F(0)
    prev = bps[0]
    vprev = m_at(S, prev)
    for t in bps[1:]:
        vt = m_at(S, t)
        # verify linearity on the piece via midpoint
        mid = (prev + t) / 2
        vm = m_at(S, mid)
        assert vm * 2 == vprev + vt, (prev, t, "kink inside piece!")
        tot += (t - prev) * (vprev + vt) / 2
        prev, vprev = t, vt
    return tot

def farey_E(N):
    # int_0^{1/(q+q')} (1 - r(q+q'))/(qq') dr = 1/(2 q q' (q+q'))
    return sum(F(1, 2 * q * q2 * (q + q2)) for (a, q), (a2, q2) in farey_pairs(N))

print("\nE1. E(AP_13) vs Farey closed form (1/2) sum 1/(q q' (q+q'))")
eap = mean_loneliness(AP)
efa = farey_E(13)
print(f"  E(AP) = {eap} = {float(eap):.8f}   Farey = {efa}   equal={eap == efa}")
egw = mean_loneliness(GW)
print(f"  exact comparison: E(GW) < E(AP)? {egw < eap}   (E(GW)-E(AP) = {float(egw - eap):.3e})")

print("\nE2. E over families (does AP minimize?  scale-matching predicts NO)")
for name, S in [("AP", AP), ("GW", GW), ("PRIMES>13", PRIMES), ("BLOCK 100..112", BLOCK100),
                ("CONSTR", CONSTR), ("RAND0", RANDOMS[0])]:
    e = mean_loneliness(S)
    print(f"  {name:15s} E = {float(e):.8f}   (exact {e})")

print("\nE3. exhaustive n=4 (3 speeds <= 15): min E and where AP ranks")
es = sorted((mean_loneliness(S), tuple(S)) for S in sets4)
eap4 = mean_loneliness([1, 2, 3])
rank = 1 + sum(1 for e, s in es if e < eap4)
print(f"  min E = {float(es[0][0]):.8f} at {es[0][1]};  E(AP)={float(eap4):.8f} rank {rank}/{len(es)}")
print(f"  bottom five: {[(float(e), s) for e, s in es[:5]]}")

print("\nDONE.")
