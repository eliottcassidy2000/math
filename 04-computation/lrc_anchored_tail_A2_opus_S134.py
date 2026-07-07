"""
lrc_anchored_tail_A2_opus_S134.py

TARGET: monad-S1's audit says the load-bearing open lemma is (A') -- per-k tail
minimality of mu_{1/7} at the AP -- NOT the mean E[maxgap] (a k=13-only sidecar
with razor margin vs T*~0.19127). This script tests a NEW candidate route to (A'):

(A'') ANCHORED TAIL DOMINATION: for a finite anchor set A (dilation acts on x, so
  any A gives a dilation-invariant statistic),
      t_A(E) := P_x( max_{a in A} gap∋a( {frac(e x): e in E} ) > 1/7 )  >=  mu_{1/7}(AP_k)
  for every k-element integer set E (min >= 1, distinct), k = 8..13.
  Since max_a gap∋a <= maxgap pointwise, (A'') ==> (A').
  At the AP itself gap∋0 = maxgap ALWAYS (Farey roof / origin-gap saturation,
  opus-S134), so with 0 in A the bound is TIGHT at the AP: t_A(AP) = mu_{1/7}(AP).

WHY anchors might suffice (mechanism): mu_{1/7}(E) = resonance windows (x near p/q,
q<=6, where the config clusters into q groups leaving gaps ~1/q > 1/7) + decorrelated
bulk. In a q-window the inter-cluster gaps sit at the q-grid positions j/q +- spread,
so anchors on a fine enough grid catch them. The adversary's counterplay: make the
big window-gap sit AWAY from all anchors (miss a residue ADJACENT to the anchor).
This script measures exactly how much that counterplay can extract.

Also reported: iid reference values, E[maxgap] (for the fleet's descent race vs
T* = 56291/294294 ~ 0.191275), and a t_A-minimizing adversarial descent.
"""
import numpy as np
from fractions import Fraction
import random, itertools, sys

GOLD = 0.6180339887498949
THETA = 1.0/7.0
TSTAR = 56291.0/294294.0

# ---------- exact AP bars from the Farey roof (opus-S134) ----------
def farey(k):
    fr = set()
    for q in range(1, k+1):
        for p in range(0, q+1):
            fr.add(Fraction(p, q))
    return sorted(fr)

def mu_from_roof(k, theta):
    F = farey(k); tot = Fraction(0)
    for a, b in zip(F[:-1], F[1:]):
        q, qp = a.denominator, b.denominator
        L = b - a; vl, vr = Fraction(1,q), Fraction(1,qp)
        if vl <= theta and vr <= theta: continue
        if vl > theta and vr > theta: tot += L; continue
        t_star = (theta - vl) * L / (vr - vl)
        tot += t_star if vl > theta else L - t_star
    return tot

BAR = {k: mu_from_roof(k, Fraction(1,7)) for k in range(8,14)}

# ---------- gap machinery (vectorized over an x-grid) ----------
def make_grid(N):
    return (np.arange(N) + GOLD) / N

def gap_stats(E, xs, anchors=(0.0, 0.5, 0.25, 0.75, 1.0/3, 2.0/3, 1.0/6, 5.0/6)):
    """Return dict with mu, E[maxgap], and anchored-gap matrix columns.
       gap∋a computed for each anchor a: the circular gap of the sorted config
       containing point a."""
    E = np.asarray(sorted(E), dtype=np.float64)
    ph = np.mod(np.outer(xs, E), 1.0)          # (N, k)
    ph.sort(axis=1)
    N, k = ph.shape
    gaps = np.empty_like(ph)
    gaps[:, :-1] = ph[:, 1:] - ph[:, :-1]
    gaps[:, -1] = ph[:, 0] + 1.0 - ph[:, -1]   # wrap gap (upper endpt = ph[:,0]+1)
    mg = gaps.max(axis=1)
    out = {"mu": float((mg > THETA).mean()), "Emg": float(mg.mean())}
    # anchored gaps: for anchor a, find gap containing a.
    anch = {}
    for a in anchors:
        # index of first phase > a  (in sorted row). gap containing a is between
        # phase[idx-1] and phase[idx] (with wrap conventions).
        idx = np.sum(ph < a, axis=1)           # 0..k
        upper = np.where(idx < k, ph[np.arange(N), np.minimum(idx, k-1)], ph[:, 0] + 1.0)
        lower = np.where(idx > 0, ph[np.arange(N), np.maximum(idx-1, 0)], ph[:, -1] - 1.0)
        anch[a] = upper - lower
    out["anch"] = anch
    return out

def tails(out, sets):
    """P(max over anchor subset > theta) for each named anchor subset."""
    res = {}
    for name, A in sets.items():
        m = None
        for a in A:
            g = out["anch"][a]
            m = g if m is None else np.maximum(m, g)
        res[name] = float((m > THETA).mean())
    return res

ANCHSETS = {
    "t1{0}":        (0.0,),
    "t2{0,1/2}":    (0.0, 0.5),
    "t3{0,1/3,2/3}":(0.0, 1.0/3, 2.0/3),
    "t4{j/4}":      (0.0, 0.25, 0.5, 0.75),
    "t6{j/6}":      (0.0, 1.0/6, 1.0/3, 0.5, 2.0/3, 5.0/6),
}

# ---------- corpus ----------
def corpus_k13():
    fams = {
        "AP {1..13}":            list(range(1,14)),
        "GW {1..11,13,24}":      list(range(1,12))+[13,24],
        "2AP+13 {2..24e,13}":    [2,4,6,8,10,12,14,16,18,20,22,24,13],
        "monad rec {2..22e,11,13}": [2,4,6,8,10,12,14,16,18,20,22,11,13],
        "3AP+2 {3..33t,13,17}":  [3,6,9,12,15,18,21,24,27,30,33,13,17],
        "3AP+2 {3..33t,11,22?dup->23}": [3,6,9,12,15,18,21,24,27,30,33,11,23],
        "kps adv {2,6,..,64}":   [2,6,8,9,10,12,14,15,16,18,22,36,64],
        "opus spread {6..58}":   [6,11,14,16,20,26,31,34,37,38,46,47,58],
        "stretch+1 {1,3..13,18,29}": [1,3,4,5,6,7,8,9,10,11,13,18,29],
        "primes<=41":            [2,3,5,7,11,13,17,19,23,29,31,37,41],
        "deep well {1..12,182}": list(range(1,13))+[182],
        "{1..12,14}":            list(range(1,13))+[14],
        "{2..14}":               list(range(2,15)),
        "odds {1,3..25}":        list(range(1,26,2)),
    }
    # evens + single odd sweep
    for m in (3,7,13,17,25):
        fams[f"2AP+{m}"] = [2*i for i in range(1,13)]+[m]
    # 11 evens + two odds sweep
    for (a,b) in ((1,13),(11,13),(7,17),(5,19),(11,15),(9,13)):
        fams[f"11e+{{{a},{b}}}"] = [2*i for i in range(1,12)]+[a,b]
    # missing-residue-adjacent adversary: no element ==1 mod 5, none ==1 mod 3
    # build from allowed residues
    allowed = [n for n in range(2,80) if n%5!=1 and n%3!=1]
    fams["missadj mod3,5"] = allowed[:13]
    return fams

def corpus_k(k):
    fams = {
        f"AP {{1..{k}}}": list(range(1,k+1)),
        f"2AP+odd": [2*i for i in range(1,k)]+[2*k-3 if (2*k-3)%2==1 else 2*k-1],
        f"spread": [3+7*i for i in range(k)],
        f"primes": [2,3,5,7,11,13,17,19,23,29,31,37,41][:k],
    }
    return fams

# ---------- descent (adversarial minimization of a chosen statistic) ----------
def descent(stat_name, seedE, xs_fast, rng, iters=250, pool=None):
    """Local search minimizing tails[stat_name] (or mu/Emg)."""
    def evalE(E):
        if len(set(E)) != len(E) or min(E) < 1: return None
        o = gap_stats(E, xs_fast)
        t = tails(o, ANCHSETS)
        t["mu"] = o["mu"]; t["Emg"] = o["Emg"]
        return t[stat_name], t
    cur = sorted(seedE); val, _ = evalE(cur)
    for it in range(iters):
        cand = list(cur)
        move = rng.random()
        i = rng.randrange(len(cand))
        if move < 0.45:
            cand[i] = max(1, cand[i] + rng.choice([-2,-1,1,2]))
        elif move < 0.8:
            cand[i] = rng.randrange(1, 70)
        else:
            j = rng.randrange(len(cand))
            cand[i] = max(1, cand[i] + rng.choice([-1,1])); cand[j] = max(1, cand[j] + rng.choice([-1,1]))
        cand = sorted(set(cand))
        while len(cand) < len(cur):
            cand.append(rng.randrange(1,70)); cand = sorted(set(cand))
        r = evalE(cand)
        if r is None: continue
        v, _ = r
        if v < val:
            cur, val = cand, v
    return cur, val

# ================= MAIN =================
def main():
    N_fast, N_fine = 40001, 400001
    xs_fine = make_grid(N_fine)
    xs_fast = make_grid(N_fast)
    rng = random.Random(134134)

    print("="*110)
    print("(1) k=13 corpus: anchored tails vs the (A') bar mu_1/7(AP_13) = 477/1078 = %.6f" % float(BAR[13]))
    print("    [also E[maxgap] vs T* = %.6f for the mean-sidecar race]" % TSTAR)
    print("="*110)
    bar13 = float(BAR[13])
    hdr = f"{'family':<28} {'mu':>7} {'t1':>7} {'t2':>7} {'t3':>7} {'t4':>7} {'t6':>7} {'Emg':>8}  verdicts"
    print(hdr); print("-"*110)
    viol = []
    for name, E in corpus_k13().items():
        if len(set(E)) != 13:
            print(f"{name:<28} SKIP (needs 13 distinct)"); continue
        o = gap_stats(E, xs_fine); t = tails(o, ANCHSETS)
        v = []
        if o["mu"] < bar13 - 1e-3: v.append("mu<BAR!! (A' viol?)")
        if t["t2{0,1/2}"] < bar13 - 1e-3: v.append("t2<bar")
        if t["t4{j/4}"] < bar13 - 1e-3: v.append("t4<bar")
        if t["t6{j/6}"] < bar13 - 1e-3: v.append("t6<bar")
        if o["Emg"] < TSTAR: v.append("Emg<T*!!")
        if v and ("mu<BAR" in v[0] or "Emg<T*!!" in v): viol.append((name, E, o, t))
        print(f"{name:<28} {o['mu']:7.4f} {t['t1{0}']:7.4f} {t['t2{0,1/2}']:7.4f} {t['t3{0,1/3,2/3}']:7.4f} "
              f"{t['t4{j/4}']:7.4f} {t['t6{j/6}']:7.4f} {o['Emg']:8.5f}  {'; '.join(v) if v else 'ok'}")

    # iid reference (simulation)
    print("-"*110)
    sim = np.random.RandomState(7).rand(200000, 13)
    sim.sort(axis=1)
    g = np.empty_like(sim); g[:, :-1] = sim[:,1:]-sim[:,:-1]; g[:,-1] = sim[:,0]+1-sim[:,-1]
    mg = g.max(axis=1)
    # anchored: gap containing 0 == wrap gap; gap containing 1/2:
    idx = (sim < 0.5).sum(axis=1); n = sim.shape[0]
    up = np.where(idx < 13, sim[np.arange(n), np.minimum(idx,12)], sim[:,0]+1)
    lo = np.where(idx > 0, sim[np.arange(n), np.maximum(idx-1,0)], sim[:,-1]-1)
    g0 = sim[:,0] + 1 - sim[:,-1]; gh = up - lo
    print(f"{'iid k=13 reference':<28} {float((mg>THETA).mean()):7.4f} {float((g0>THETA).mean()):7.4f} "
          f"{float((np.maximum(g0,gh)>THETA).mean()):7.4f} {'':>7} {'':>7} {'':>7} {float(mg.mean()):8.5f}")

    print()
    print("="*110)
    print("(2) per-k check (k=8..12): t2/t4 anchored tails vs the per-k bar mu_1/7(AP_k)")
    print("="*110)
    for k in range(8, 13):
        bark = float(BAR[k])
        print(f"  -- k={k}: bar = {BAR[k]} = {bark:.6f}")
        for name, E in corpus_k(k).items():
            if len(set(E)) != k: continue
            o = gap_stats(E, xs_fine); t = tails(o, ANCHSETS)
            flag = []
            if o["mu"] < bark - 1e-3: flag.append("mu<bar (A' VIOLATION?)")
            if t["t2{0,1/2}"] < bark - 1e-3: flag.append("t2<bar")
            if t["t4{j/4}"] < bark - 1e-3: flag.append("t4<bar")
            print(f"     {name:<22} mu={o['mu']:.4f} t1={t['t1{0}']:.4f} t2={t['t2{0,1/2}']:.4f} "
                  f"t4={t['t4{j/4}']:.4f} t6={t['t6{j/6}']:.4f}  {'; '.join(flag) if flag else 'ok'}")

    print()
    print("="*110)
    print("(3) adversarial descent minimizing each statistic (k=13, fast grid, refined on fine grid)")
    print("="*110)
    seeds = [list(range(1,14)), [2,4,6,8,10,12,14,16,18,20,22,24,13],
             [2,4,6,8,10,12,14,16,18,20,22,11,13],
             [6,11,14,16,20,26,31,34,37,38,46,47,58],
             sorted(rng.sample(range(1,60), 13))]
    for stat in ("mu", "t2{0,1/2}", "t4{j/4}", "Emg"):
        best, bestv = None, 9.9
        for sd in seeds:
            E, v = descent(stat, sd, xs_fast, rng, iters=220)
            if v < bestv: best, bestv = E, v
        o = gap_stats(best, xs_fine); t = tails(o, ANCHSETS)
        t["mu"] = o["mu"]; t["Emg"] = o["Emg"]
        ref = bar13 if stat != "Emg" else TSTAR
        print(f"  min {stat:<10} -> {best}  fine-grid value={t[stat]:.5f} (bar/T*={ref:.5f})"
              f" {'*** BELOW ***' if t[stat] < ref else '(above)'}")
        print(f"      full row: mu={o['mu']:.4f} t2={t['t2{0,1/2}']:.4f} t4={t['t4{j/4}']:.4f} Emg={o['Emg']:.5f}")

    print()
    print("READOUT: (A'') candidate = smallest anchor set whose tail column never dips below")
    print("the per-k bar. t1 alone is expected to FAIL (iid 0.427 < 0.4425). The question is")
    print("whether t2 (the observer + its antipode) or t4/t6 survive the descent.")

if __name__ == "__main__":
    main()
