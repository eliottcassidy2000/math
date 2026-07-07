"""
lrc_per_q_audit_opus_S144.py   (opus-2026-07-07-S144, HYP-5137)

VALIDITY AUDIT of the per-resonance-q window reduction (kps-S72, HYP-5117), plus the
EXACT AP-side closed form it was missing.

FINDINGS TESTED HERE:
 (A) THE q=1 HOLE.  The S72 attribution (nearest p/q Voronoi, q <= 6) is a genuine
     partition, so mu = SUM_{q=1}^{6} W_q -- but the reflection states mu = sum_{q=2..6}
     and never reports W_1.  W_1 is the ORIGIN-resonance window: for the AP it is
     12/(7(k-1)) (exactly 1/7 at k=13 = 32% of mu(AP_13)), and the AP (minimal diameter)
     is the W_1 MAXIMIZER, not minimizer.  So "AP minimizes each W_q" is FALSE at q=1
     and the per-q program needs an origin-budget clause.
 (B) EXACT AP WINDOWS FROM THE ROOF (THM-637).  On the Farey-k grid the roof is affine
     per cell with node values 1/q; the 1/7-superlevel width on the cell flanking node
     p/q (flank denominator qF) is (7-q)/(7 q (qF - q)) for k >= 12 (no merging).  So
       W_q(AP_k) = SUM_{p: gcd(p,q)=1} (7-q)/(7q) * [ 1/(qL-q) + 1/(qR-q) ],
     qL,qR = Farey-k neighbor denominators of p/q.  k=13 exact values:
       W_2 = 5/77, W_3 = 3/35, W_4 = 8/147, W_5 = 23/294, W_6 = 4/245,
       and mu(AP_13) = 1/7 + (5/77+3/35+8/147+23/294+4/245) = 477/1078  (roof check).
 (C) THE LIFT-FAMILY ADVERSARY (the S65/MISTAKE-116 escape shape, absent from the S72
     test set).  E_m = {1..12, 13+60m} shares ALL residues mod 2..6 with AP_13 but the
     lifted element makes one cluster edge move at huge rate ~60m near every q<=6 node,
     collapsing ONE side of each window pair => prediction: per-q windows LOSE mass vs
     AP while total mu stays >= mu(AP) via redistributed far-scale mass.  If confirmed:
     W_q(E_m) < W_q(AP) for some q -- the per-q minimality claim as stated is REFUTED
     by the repo's own canonical adversary family (weak-adversary trap, taxonomy #1),
     and the program must be restated with rate/diameter corrections.
 (D) CLUSTER-RATE LINEARIZATION = THE ROOF, generalized: near x = p/q + delta the
     inter-cluster gap_j(delta) = G_j + delta*(min_{C_{j+1}} e - max_{C_j} e); at the AP
     the closing rate at every q<=6 node equals qF - q (Farey-flank), reproducing (B).
     Verified numerically at node 1/2.
"""
from fractions import Fraction as F
from math import gcd
import numpy as np
import time

# ---------------------------------------------------------------- exact AP-side (B)
def farey_neighbors(p, q, k):
    """Neighbors of p/q in F_k: (pL,qL), (pR,qR) with |p qX - q pX| = 1."""
    # right neighbor: smallest a/b > p/q with pb - qa = -1  => a = (p*b+1)/q
    best = None
    for b in range(1, k + 1):
        if (p * b + 1) % q == 0:
            a = (p * b + 1) // q
            f = F(a, b)
            if best is None or f < best[0]:
                best = (f, a, b)
    right = best
    best = None
    for b in range(1, k + 1):
        if (p * b - 1) % q == 0:
            a = (p * b - 1) // q
            f = F(a, b)
            if best is None or f > best[0]:
                best = (f, a, b)
    left = best
    return (left[1], left[2]), (right[1], right[2])

def W_AP_exact(k, theta=F(1, 7)):
    """Exact per-q window masses of AP_k from the roof, q = 1..6 (k >= 12: no merging).
       Returns dict q -> Fraction, plus mu = sum (checked against the known constant)."""
    Q = int(1 / theta) - 1  # 6
    W = {}
    for q in range(1, Q + 1):
        tot = F(0)
        ps = [p for p in range(0 if q == 1 else 1, q + (1 if q == 1 else 0))
              if gcd(p, q) == 1] if q > 1 else [0, 1]
        for p in ps:
            if q == 1:
                # node 0 (right side only) / node 1 (left side only)
                if p == 0:
                    (_, _), (pR, qR) = farey_neighbors(0, 1, k)
                    tot += (1 - theta) / (qR - 1) * 1  # width (1/1-1/7)/(qR-1)... see below
                else:
                    (pL, qL), (_, _) = farey_neighbors(1, 1, k)
                    tot += (1 - theta) / (qL - 1)
            else:
                (pL, qL), (pR, qR) = farey_neighbors(p, q, k)
                for qF in (qL, qR):
                    tot += (F(1, q) - theta) / (qF - q)
        W[q] = tot
    return W

# ---------------------------------------------------------------- numeric engine
def mu_and_windows_np(E, res=2_000_000, theta=1 / 7, Q=6, chunk=200_000):
    """Vectorized mu + nearest-rational-Voronoi per-q attribution (q = 1..Q)."""
    E = np.asarray(sorted(E), dtype=np.float64)
    # Voronoi: precompute all p/q nodes with denominator label (lowest terms, q<=Q)
    nodes = []
    for q in range(1, Q + 1):
        for p in range(0, q + 1):
            if gcd(p, q) == 1:
                nodes.append((p / q, q))
    nv = np.array([n[0] for n in nodes])
    nq = np.array([n[1] for n in nodes])
    Wc = np.zeros(Q + 1, dtype=np.int64)
    good_total = 0
    for lo in range(0, res, chunk):
        x = (np.arange(lo, min(lo + chunk, res)) + 0.5) / res
        ph = np.sort((x[:, None] * E[None, :]) % 1.0, axis=1)
        gaps = np.diff(ph, axis=1)
        wrap = ph[:, 0] + 1 - ph[:, -1]
        mg = np.maximum(gaps.max(axis=1), wrap)
        good = mg > theta
        good_total += int(good.sum())
        if good.any():
            xg = x[good]
            d = np.abs(xg[:, None] - nv[None, :])
            d = np.minimum(d, 1 - d)
            qq = nq[np.argmin(d, axis=1)]
            for q in range(1, Q + 1):
                Wc[q] += int((qq == q).sum())
    mu = good_total / res
    return mu, {q: Wc[q] / res for q in range(1, Q + 1)}

def report(name, E, WAPn, muAPn):
    mu, W = mu_and_windows_np(E)
    flags = []
    for q in range(1, 7):
        if W[q] < WAPn[q] - 0.002:
            flags.append(f"q{q}: {W[q]:.4f} < AP {WAPn[q]:.4f}  *** BELOW AP ***")
    wstr = " ".join(f"{W[q]:.4f}" for q in range(1, 7))
    print(f"    {name:32s} mu={mu:.4f} ({mu-muAPn:+.4f})  W_1..6: {wstr}"
          f"{'   ' + '; '.join(flags) if flags else '   per-q >= AP holds (incl q=1)'}")
    return mu, W

def main():
    t0 = time.time()
    k = 13
    print("=" * 100)
    print("(B) EXACT AP-side windows from the roof (THM-637), k=13, theta=1/7")
    print("=" * 100)
    Wx = W_AP_exact(k)
    tot = sum(Wx.values())
    print("    q :  exact W_q(AP_13)        float")
    for q in range(1, 7):
        print(f"    {q} :  {str(Wx[q]):>12}          {float(Wx[q]):.6f}")
    print(f"    sum = {tot} = {float(tot):.6f}   vs canon mu(AP_13) = 477/1078 = {477/1078:.6f}"
          f"   {'MATCH' if tot == F(477,1078) else '*** MISMATCH ***'}")
    print(f"    W_1(AP_13) = {Wx[1]} = {float(Wx[1]):.4f} = {float(Wx[1]/tot)*100:.1f}% of mu -- the unreported origin windows")

    print()
    print("=" * 100)
    print("(A)+(C) numeric audit: W_q INCLUDING q=1; the lift-family adversary vs kps-S72 families")
    print("=" * 100)
    AP = list(range(1, k + 1))
    muAP, WAP = mu_and_windows_np(AP)
    wstr = " ".join(f"{WAP[q]:.4f}" for q in range(1, 7))
    print(f"    {'AP {1..13}':32s} mu={muAP:.4f} (+0.0000)  W_1..6: {wstr}")
    print(f"    [exact check: numeric W_q vs roof closed form: "
          + " ".join(f"d{q}={abs(WAP[q]-float(Wx[q])):.4f}" for q in range(1, 7)) + "]")
    print()
    fams = {
        "lift m=1 {1..12, 73}": list(range(1, 13)) + [73],
        "lift m=5 {1..12, 313}": list(range(1, 13)) + [313],
        "lift m=20 {1..12, 1213}": list(range(1, 13)) + [1213],
        "lift420 {1..12, 433}": list(range(1, 13)) + [433],
        "double-lift {1..11,12+60,13+120}": list(range(1, 12)) + [72, 133],
        "geometric": [2 ** j for j in range(k)],
        "primes": [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41][:k],
        "all-odd (2AP-1)": [2 * j - 1 for j in range(1, k + 1)],
    }
    results = {}
    for nm, E in fams.items():
        results[nm] = report(nm, E, WAP, muAP)

    print()
    print("=" * 100)
    print("(D) cluster-rate linearization at node 1/2 (q=2): window half-widths, one per side")
    print("=" * 100)
    for nm, E in [("AP_13", AP), ("lift m=5", list(range(1, 13)) + [313])]:
        E = sorted(E)
        evens = [e for e in E if e % 2 == 0]; odds = [e for e in E if e % 2 == 1]
        # gap A: top of 0-cluster (max even) to bottom of 1/2-cluster (min odd): closes at rate max(evens)-min(odds)... signed for delta>0
        rateA = max(evens) - min(odds)      # gap = 1/2 - delta*rateA
        rateB = max(odds) - min(evens)      # wrap gap = 1/2 - delta*rateB
        wA = (1 / 2 - 1 / 7) / rateA if rateA > 0 else float("inf")
        wB = (1 / 2 - 1 / 7) / rateB if rateB > 0 else float("inf")
        # numeric: scan delta > 0 for maxgap > 1/7 window end
        deltas = np.linspace(1e-6, 0.05, 200_000)
        x = 0.5 + deltas
        Ea = np.array(E, dtype=float)
        ph = np.sort((x[:, None] * Ea[None, :]) % 1.0, axis=1)
        gaps = np.diff(ph, axis=1)
        wrap = ph[:, 0] + 1 - ph[:, -1]
        mg = np.maximum(gaps.max(axis=1), wrap)
        idx = np.nonzero(mg <= 1 / 7)[0]
        first_exit = deltas[idx[0]] if len(idx) else None
        print(f"    {nm:12s}: predicted right-window = min gap-survival = "
              f"min((1/2-1/7)/{rateA}, (1/2-1/7)/{rateB}) = {min(wA, wB):.6f}; "
              f"numeric first exit at delta = {first_exit:.6f}"
              f"   {'MATCH' if first_exit and abs(first_exit - min(wA, wB)) < 1e-3 else '(see discussion)'}")

    print(f"\n[{time.time()-t0:.0f}s]")

if __name__ == "__main__":
    main()
