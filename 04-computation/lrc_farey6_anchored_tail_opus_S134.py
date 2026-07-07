"""
lrc_farey6_anchored_tail_opus_S134.py

Refined (A'') candidate after the first experiment (lrc_anchored_tail_A2):
  * ALL observed anchored-tail failures were AFFINE IMAGES of the AP; mu is
    affine-invariant but anchored tails are not => state the lemma over
    AFFINE-NORMALIZED families N(E) = (E - min E)/gcd(diffs) + 1.
  * The principled anchor set is F6 = {p/q in [0,1): q <= 6} (12 anchors):
    WINDOW-EXACTNESS mechanism: near x = a/q + delta (q<=6, delta>0 small) the
    config clusters at positions {r/q} each spreading UPWARD by e*delta, so every
    inter-cluster gap contains the Farey point at its right cluster base (and for
    delta<0, its left base). Hence inside the q<=6 resonance windows
        max_{a in F6} gap∋a == maxgap   EXACTLY.
    Since mu_1/7(AP) is PURE window measure (Farey roof), t_F6(AP) = bar exactly.

CANDIDATE (A''-F6): for every affine-normalized k-set E (k=8..13):
    t_F6(E) := P_x( max_{a in F6} gap∋a > 1/7 )  >=  mu_1/7(AP_k)   [= (A') restricted
    to normalized E, which by affine invariance of mu implies full (A')] -- IF t_F6
    survives; otherwise measure the deficit and identify the bulk/window split.

Diagnostics per family: mu, t_F6, window-mass vs bulk-mass of {maxgap>1/7},
window-exactness violation rate, anchored bulk fraction.
"""
import numpy as np
from fractions import Fraction
from math import gcd
import random

GOLD = 0.6180339887498949
THETA = 1.0/7.0

def farey_points(Q):
    pts = set()
    for q in range(1, Q+1):
        for p in range(q):
            pts.add(Fraction(p, q))
    return sorted(float(f) for f in pts)

F6 = farey_points(6)   # 12 anchors incl 0

def farey(k):
    fr = set()
    for q in range(1, k+1):
        for p in range(0, q+1):
            fr.add(Fraction(p, q))
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
    return tot

BAR = {k: float(mu_from_roof(k)) for k in range(8,14)}

def normalize(E):
    E = sorted(set(E))
    d = [e - E[0] for e in E[1:]]
    g = 0
    for x in d: g = gcd(g, x)
    if g == 0: g = 1
    return [1] + [1 + x//g for x in d]

def gap_data(E, xs):
    E = np.asarray(sorted(E), dtype=np.float64)
    ph = np.mod(np.outer(xs, E), 1.0)
    ph.sort(axis=1)
    N, k = ph.shape
    gaps = np.empty_like(ph)
    gaps[:, :-1] = ph[:, 1:] - ph[:, :-1]
    gaps[:, -1] = ph[:, 0] + 1.0 - ph[:, -1]
    mg = gaps.max(axis=1)
    # anchored max over F6
    am = np.zeros(N)
    idxcache = np.arange(N)
    for a in F6:
        idx = np.sum(ph < a, axis=1)
        upper = np.where(idx < k, ph[idxcache, np.minimum(idx, k-1)], ph[:, 0] + 1.0)
        lower = np.where(idx > 0, ph[idxcache, np.maximum(idx-1, 0)], ph[:, -1] - 1.0)
        am = np.maximum(am, upper - lower)
    return mg, am

# window classifier: x within 2/(q*D) of some p/q (q<=6), D = diam(E)
def window_mask(E, xs):
    D = max(E) - min(E)
    m = np.zeros(len(xs), dtype=bool)
    for q in range(1, 7):
        w = 2.0/(q*D)
        # distance of x*q to nearest integer
        dq = np.abs(xs*q - np.round(xs*q))/q
        m |= (dq <= w)
    return m

def full_row(E, xs, name, bar):
    mg, am = gap_data(E, xs)
    mu = float((mg > THETA).mean())
    tF6 = float((am > THETA).mean())
    wm = window_mask(E, xs)
    win_mu = float(((mg > THETA) & wm).mean())
    blk_mu = mu - win_mu
    # window-exactness: on windows where maxgap>theta, does anchored max == maxgap?
    sel = (mg > THETA) & wm
    viol = float((np.abs(am[sel] - mg[sel]) > 1e-9).mean()) if sel.sum() else 0.0
    # anchored bulk capture: among bulk good x, fraction anchored max also > theta
    selb = (mg > THETA) & (~wm)
    cap = float((am[selb] > THETA).mean()) if selb.sum() else 1.0
    flag = ""
    if mu < bar - 1e-3: flag += " MU<BAR(A'-viol?)"
    if tF6 < bar - 1e-3: flag += " tF6<bar"
    print(f"{name:<30} mu={mu:.4f} tF6={tF6:.4f} | win={win_mu:.4f} bulk={blk_mu:.4f} "
          f"| wexact-viol={viol:.3f} bulkcap={cap:.3f}{flag}")
    return mu, tF6

def descent_tF6(seed, xs_fast, rng, iters=300, k=13):
    def ev(E):
        if len(set(E)) != k: return None
        En = normalize(E)
        mg, am = gap_data(En, xs_fast)
        return float((am > THETA).mean()), En
    cur = normalize(seed); val, _ = ev(cur)
    for _ in range(iters):
        cand = list(cur)
        i = rng.randrange(k)
        r = rng.random()
        if r < 0.5: cand[i] = max(1, cand[i] + rng.choice([-3,-2,-1,1,2,3]))
        elif r < 0.85: cand[i] = rng.randrange(1, 75)
        else:
            j = rng.randrange(k); cand[i] = max(1, cand[i]+rng.choice([-1,1])); cand[j] = max(1, cand[j]+rng.choice([-1,1]))
        cand = sorted(set(cand))
        while len(cand) < k: cand.append(rng.randrange(1,75)); cand = sorted(set(cand))
        out = ev(cand)
        if out is None: continue
        v, En = out
        if v < val: cur, val = En, v
    return cur, val

def main():
    N_fine, N_fast = 400001, 50001
    xs = (np.arange(N_fine)+GOLD)/N_fine
    xs_fast = (np.arange(N_fast)+GOLD)/N_fast
    rng = random.Random(613461346)

    print("="*112)
    print(f"(1) NORMALIZED k=13 corpus: t_F6 vs bar {BAR[13]:.6f}; window/bulk split of mu; window-exactness check")
    print("="*112)
    corpus = {
        "AP": list(range(1,14)),
        "GW": list(range(1,12))+[13,24],
        "2AP+13": [2,4,6,8,10,12,14,16,18,20,22,24,13],
        "monad rec": [2,4,6,8,10,12,14,16,18,20,22,11,13],
        "t2-breaker {1,3..13,15}": [1,3,4,5,6,7,8,9,10,11,12,13,15],
        "AP del2 ins15": [1,3,4,5,6,7,8,9,10,11,12,13,15],
        "AP del7 ins16": [1,2,3,4,5,6,8,9,10,11,12,13,16],
        "AP del{2,3} ins{15,17}": [1,4,5,6,7,8,9,10,11,12,13,15,17],
        "kps adv": [2,6,8,9,10,12,14,15,16,18,22,36,64],
        "opus spread": [6,11,14,16,20,26,31,34,37,38,46,47,58],
        "stretch": [1,3,4,5,6,7,8,9,10,11,13,18,29],
        "primes": [2,3,5,7,11,13,17,19,23,29,31,37,41],
        "deep well": list(range(1,13))+[182],
        "missadj mod3,5": None,  # fill below
    }
    allowed = [n for n in range(2,80) if n%5!=1 and n%3!=1]
    corpus["missadj mod3,5"] = allowed[:13]
    worst = (9.9, None)
    for name, E in corpus.items():
        En = normalize(E)
        mu, tF6 = full_row(En, xs, name, BAR[13])
        if tF6 < worst[0]: worst = (tF6, name)
    print(f"  worst corpus t_F6: {worst[0]:.4f} at {worst[1]}   (bar {BAR[13]:.4f})")

    print()
    print("="*112)
    print("(2) dedicated NORMALIZED descent minimizing t_F6 (k=13, then k=12, k=10)")
    print("="*112)
    for k, bar in ((13, BAR[13]), (12, BAR[12]), (10, BAR[10])):
        seeds = [list(range(1,k+1))]
        seeds.append([2*i for i in range(1,k)]+[2*k-3])
        seeds.append([1,3]+[i for i in range(4,k+2)])
        seeds.append(sorted(rng.sample(range(1,70), k)))
        seeds.append(sorted(rng.sample(range(1,40), k)))
        best, bestv = None, 9.9
        for sd in seeds:
            E, v = descent_tF6(sd, xs_fast, rng, iters=280, k=k)
            if v < bestv: best, bestv = E, v
        mg, am = gap_data(best, xs)
        tfine = float((am > THETA).mean()); mufine = float((mg > THETA).mean())
        print(f"  k={k}: min t_F6 -> {best}")
        print(f"        fine grid: t_F6={tfine:.5f}  mu={mufine:.5f}  bar={bar:.5f}  "
              f"{'*** t_F6 BELOW bar ***' if tfine < bar - 5e-4 else 'holds'}"
              f"{'  [mu below bar too!! (A -viol?)]' if mufine < bar - 5e-4 else ''}")

    print()
    print("="*112)
    print("(3) sanity: t_F6(AP_k) == bar exactly (window-only structure), k=8..13")
    print("="*112)
    for k in range(8,14):
        mg, am = gap_data(list(range(1,k+1)), xs)
        print(f"  k={k}: t_F6(AP)={float((am>THETA).mean()):.6f}  bar={BAR[k]:.6f}  mu={float((mg>THETA).mean()):.6f}")

if __name__ == "__main__":
    main()
