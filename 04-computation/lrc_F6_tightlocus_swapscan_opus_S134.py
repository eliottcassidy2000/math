"""
lrc_F6_tightlocus_swapscan_opus_S134.py

Stress-test of (A''-F6) on the TIGHT LOCUS (near-AP swap families), where the
anchored-tail bound has the least room (t_F6(AP) = bar exactly).

Scan (normalized): all 1-swap families {1..k}\{j} u {m}, j=1..k, m=k+1..k+27,
plus all 2-swap {1..k}\{j1,j2} u {m1,m2} with m in k+1..k+14 (sampled cap), at
k=13 and k=12. Report the minimum t_F6 margin over the scan and any family with
t_F6 < bar (also track mu for (A') itself).

Also: the F7 EXACT IDENTITY check -- max spacing of F7 is 1/7, so
{max_{a in F7} gap∋a > 1/7} == {maxgap > 1/7} pointwise. Verify numerically.
"""
import numpy as np
from fractions import Fraction
from math import gcd
import itertools, random

GOLD = 0.6180339887498949
THETA = 1.0/7.0

def farey_points(Q):
    pts = set()
    for q in range(1, Q+1):
        for p in range(q):
            pts.add(Fraction(p, q))
    return sorted(float(f) for f in pts)

F6 = farey_points(6); F7 = farey_points(7)

def farey(k):
    fr = set()
    for q in range(1, k+1):
        for p in range(0, q+1): fr.add(Fraction(p, q))
    return sorted(fr)

def mu_from_roof(k, theta=Fraction(1,7)):
    F = farey(k); tot = Fraction(0)
    for a, b in zip(F[:-1], F[1:]):
        q, qp = a.denominator, b.denominator
        L = b - a; vl, vr = Fraction(1,q), Fraction(1,qp)
        if vl <= theta and vr <= theta: continue
        if vl > theta and vr > theta: tot += L; continue
        ts = (theta - vl) * L / (vr - vl)
        tot += ts if vl > theta else L - ts
    return float(tot)

def normalize(E):
    E = sorted(set(E)); d = [e - E[0] for e in E[1:]]
    g = 0
    for x in d: g = gcd(g, x)
    if g == 0: g = 1
    return [1] + [1 + x//g for x in d]

def stats(E, xs, anchors):
    E = np.asarray(sorted(E), dtype=np.float64)
    ph = np.mod(np.outer(xs, E), 1.0); ph.sort(axis=1)
    N, k = ph.shape
    gaps = np.empty_like(ph)
    gaps[:, :-1] = ph[:, 1:] - ph[:, :-1]; gaps[:, -1] = ph[:, 0] + 1.0 - ph[:, -1]
    mg = gaps.max(axis=1)
    am = np.zeros(N); rows = np.arange(N)
    for a in anchors:
        idx = np.sum(ph < a, axis=1)
        up = np.where(idx < k, ph[rows, np.minimum(idx, k-1)], ph[:, 0] + 1.0)
        lo = np.where(idx > 0, ph[rows, np.maximum(idx-1, 0)], ph[:, -1] - 1.0)
        am = np.maximum(am, up - lo)
    return mg, am

def main():
    N = 160001
    xs = (np.arange(N)+GOLD)/N
    rng = random.Random(99)

    # --- F7 exact-identity check on a few families ---
    print("(0) F7 pointwise identity {max_a gap_at_a>1/7} == {maxgap>1/7}:")
    for E in (list(range(1,14)), [2,4,6,8,10,12,14,16,18,20,22,11,13],
              [6,11,14,16,20,26,31,34,37,38,46,47,58]):
        mg, am7 = stats(E, xs, F7)
        dis = float(((mg > THETA) != (am7 > THETA)).mean())
        print(f"    {str(E)[:44]:<46} disagreement measure = {dis:.2e}")

    for k in (13, 12):
        bar = mu_from_roof(k)
        AP = list(range(1, k+1))
        print(f"\n(1) k={k} 1-swap scan (normalized), bar={bar:.6f}")
        worst = (9.9, None); nbelow_t = 0; nbelow_mu = 0; count = 0
        for j in range(1, k+1):
            for m in range(k+1, k+28):
                E = normalize([i for i in range(1, k+1) if i != j] + [m])
                if len(E) != k: continue
                mg, am = stats(E, xs, F6)
                t = float((am > THETA).mean()); mu = float((mg > THETA).mean())
                count += 1
                if t < worst[0]: worst = (t, (j, m, tuple(E)), mu)
                if t < bar - 1.5e-3: nbelow_t += 1
                if mu < bar - 1.5e-3: nbelow_mu += 1
        print(f"    scanned {count};  t_F6 below bar: {nbelow_t};  mu below bar (A'-viol): {nbelow_mu}")
        print(f"    worst t_F6 = {worst[0]:.5f} (mu={worst[2]:.5f}) at swap j={worst[1][0]}->m={worst[1][1]}")

        print(f"(2) k={k} 2-swap sample (normalized)")
        worst2 = (9.9, None); nbelow2 = 0; c2 = 0
        pairs = list(itertools.combinations(range(1, k+1), 2))
        for (j1, j2) in pairs:
            for _ in range(6):
                m1, m2 = rng.sample(range(k+1, k+15), 2)
                E = normalize([i for i in range(1, k+1) if i not in (j1, j2)] + [m1, m2])
                if len(E) != k: continue
                mg, am = stats(E, xs, F6)
                t = float((am > THETA).mean()); mu = float((mg > THETA).mean())
                c2 += 1
                if t < worst2[0]: worst2 = (t, (j1, j2, m1, m2), mu)
                if t < bar - 1.5e-3: nbelow2 += 1
        print(f"    scanned {c2};  t_F6 below bar: {nbelow2}")
        print(f"    worst t_F6 = {worst2[0]:.5f} (mu={worst2[2]:.5f}) at {worst2[1]}")

if __name__ == "__main__":
    main()
