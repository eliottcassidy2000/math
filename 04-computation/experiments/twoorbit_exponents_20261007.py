#!/usr/bin/env python3
"""How the two orbits' exponents correlate over long gaps (mac-mini-2026-10-07-oaimath3).

Setting (S19's lag-1 Mersenne debt pairs, HYP-9217): y odd (Haar proxy: a random odd M-bit integer), v = 2 + Geom(1/2),
x0 = 3 2^v y + 1, and the y-orbit is advanced one step, so x_t is compared with y_t := U^(t+1)(y).  Exponents
a_t = v_2(3 y_t + 1), b_t = v_2(3 x_t + 1); A_t, B_t their partial sums; L_t = L_0 + A_t - B_t (x_t = 2^(L_t) y_t + Delta_t).

The tape picture.  Both orbits read the same 2-adic digits of y; the x-orbit's head sits L_t digits behind the y-orbit's.
EXACT ALIGNMENT: B_t - B_s = L_s (the x-orbit, from time s to t, consumed exactly the L_s digits it was behind), i.e. at
time t the x-orbit starts reading exactly where the y-orbit started reading at time s.  Then
      x_t = 3^k y_s + kappa,   k = t - s,   kappa = (3^k Delta_s + c)/2^(L_s) an integer fixed by the past,
and with E := 3 kappa + 1 - 3^k, delta := v_2(E):
      b_t = a_s if a_s < delta;  b_t = delta if a_s > delta;  b_t >= delta if a_s = delta        (ultrametric coupling).
Conditionally on the past, y_s is Haar, so Cov(a_s, b_t | delta) = 2 - 3 * 2^(1 - delta), which averages to 0 iff
P(delta >= d) = 2^(1-d) (the Haar law of an even 2-adic number's valuation).

The archimedean side.  rho_t = D_t / 2^(L_t), D_t = 3 Delta_t + 1 - 2^(L_t), satisfies exactly
      rho_(t+1) = (3 / 2^(a_t)) rho_t - 1 + 2^(-L_(t+1)),
a random affine recursion driven by the y-exponents only; E[3/2^a] = 1 for a ~ Geom(1/2), so (Kesten-Goldie) its
stationary law has tail index exactly 1; E[(3/2^a)^theta] = 3^theta/(2^(theta+1) - 1) = THM-4554's Moran function rho(theta+1).

Checks:
  1. EXACT (every alignment of every sample): the identity x_t = 3^k y_s + kappa with kappa an integer, the ultrametric law,
     and past-measurability of kappa (re-run with the tape above the alignment position replaced: kappa unchanged).
  2. EXACT (every step): the rho recursion.
  3. NUMERICAL: the law of delta at alignments (overall, by k mod 4, by the misalignment class); the cross-covariance
     function Cov(a_s, b_(s+k)) by L_s-bin and lag (looking for the surfacing lag k ~ L_s/2); the aligned covariance;
     block variances of L; the tail of |rho|.
Usage: python3 twoorbit_exponents_20261007.py [NPAIRS BITS SEED NPROC]   (defaults 400 16000 2027 8)
"""
import sys, random, math, time
from fractions import Fraction
from multiprocessing import Pool
import numpy as np

FAIL = []


def check(cond, msg):
    print(('  ok   ' if cond else '  FAIL ') + msg, flush=True)
    if not cond:
        FAIL.append(msg)


def v2(n):
    return (n & -n).bit_length() - 1


args = [int(a) for a in sys.argv[1:5]]
NPAIRS, BITS, SEED, NPROC = (args + [400, 16000, 2027, 8][len(args):])[:4]
MARGIN = 96
KMAX = 160          # lags for the cross-covariance function
LBINS = [(8, 16), (16, 32), (32, 64), (64, 128)]


def run_pair(idx, perturb_at=None):
    """Exact integer orbits of the pair; returns per-step arrays and per-alignment records.
    perturb_at: (s_align, newbits) -> replace the bits of y above the tape position of y-step s_align (past-measurability test)."""
    rng = random.Random(SEED * 1000003 + idx)
    y0 = rng.getrandbits(BITS) | 1
    v = 2
    while rng.random() < 0.5:
        v += 1
    if perturb_at is not None:
        pos, newhigh = perturb_at
        y0 = (y0 & ((1 << pos) - 1)) | (newhigh << pos)
    x = 3 * (1 << v) * y0 + 1
    t3 = 3 * y0 + 1
    a0 = v2(t3)
    y = t3 >> a0
    L0 = v + a0
    L = L0
    ys, xs, a_l, b_l, L_l, rho_ok = [], [], [], [], [], True
    A = B = 0
    Apos = {}                      # A_s -> s  (y-head position at step s)
    consumed_y = a0                # bits of y0 consumed so far (for the precision budget)
    maxsteps = int((BITS - MARGIN) / 2.3)
    Delta = Fraction(x) - Fraction(2) ** L * y
    rho = (3 * Delta + 1 - Fraction(2) ** L) / Fraction(2) ** L
    aligns = []
    prev_align = None
    for t in range(maxsteps):
        if x == y:
            break
        Apos[A] = t
        ys.append(y); xs.append(x)
        a = v2(3 * y + 1)
        b = v2(3 * x + 1)
        consumed_y += a
        if consumed_y > BITS - MARGIN:
            ys.pop(); xs.pop()
            break
        # exact alignment: the x-head now (tape coordinate B - L0) equals the y-head at some s <= t
        s = Apos.get(B - L0)
        if s is not None and s < t:
            k = t - s
            kappa = x - 3 ** k * ys[s]
            E = 3 * kappa + 1 - 3 ** k
            delta = v2(E) if E != 0 else 10 ** 6
            a_s = a_l[s]
            if a_s < delta:
                law = (b == a_s)
            elif a_s > delta:
                law = (b == delta)
            else:
                law = (b >= delta)
            # continuation: (s-1, t-1) aligned and b_(t-1) = a_(s-1) (lockstep); then delta = min(prev_delta - a, nu_k) unless tie
            cont = 0
            lock_ok = True
            if prev_align is not None and prev_align[0] == s - 1 and prev_align[1] == t - 1 and b_l[t - 1] == a_l[s - 1]:
                cont = 1
                pd, pa = prev_align[2], a_l[s - 1]
                nu = v2(3 ** k - 1)
                if pd - pa != nu:
                    lock_ok = (delta == min(pd - pa, nu))
                else:
                    lock_ok = (delta > nu)
            prev_align = (s, t, delta)
            aligns.append((s, t, k, delta, a_s, b, L_l[s], kappa & ((1 << 40) - 1), law, cont, lock_ok))
        ny = (3 * y + 1) >> a
        nx = (3 * x + 1) >> b
        Lp = L + a - b
        # exact rho recursion check on the first 200 steps
        if t < 200:
            Dn = Fraction(nx) - Fraction(2) ** Lp * ny
            rho_n = (3 * Dn + 1 - Fraction(2) ** Lp) / Fraction(2) ** Lp
            if rho_n != Fraction(3, 2 ** a) * rho - 1 + Fraction(1, 1) / Fraction(2) ** Lp:
                rho_ok = False
            rho = rho_n
        a_l.append(a); b_l.append(b); L_l.append(L)
        A += a; B += b
        x, y, L = nx, ny, Lp
    return dict(a=np.array(a_l, dtype=np.int64), b=np.array(b_l, dtype=np.int64), L=np.array(L_l, dtype=np.int64),
                aligns=aligns, rho_ok=rho_ok, L0=L0, v=v, ys=None)


def past_measurability_test(idx, nmax=30):
    """For alignments in pair idx: replace the bits of y0 strictly above the alignment's tape position and check that kappa
    and delta are unchanged (kappa is a function of the past only)."""
    base = run_pair(idx)
    rng = random.Random(77 + idx)
    ok = True
    tested = 0
    # tape position of y-step s in y0-bit coordinates: bits of y0 consumed before y_s is formed = a0 + A_s (+1 parity bit)
    rng0 = random.Random(SEED * 1000003 + idx)
    y0 = rng0.getrandbits(BITS) | 1
    t3 = 3 * y0 + 1
    a0 = v2(t3)
    Acum = np.concatenate([[0], np.cumsum(base['a'])])
    for (s, t, k, delta, a_s, b, Ls, kap, law, cont, lk) in base["aligns"][:nmax]:
        if t >= len(base['a']) - 5 or delta > 60:
            continue
        pos = a0 + int(Acum[s]) + 1            # y_s depends on y0 bits >= pos only through the fresh tape
        newhigh = rng.getrandbits(BITS - pos) | 1 if BITS - pos > 0 else 0
        pert = run_pair(idx, perturb_at=(pos, newhigh))
        match = [r for r in pert['aligns'] if r[0] == s and r[1] == t]
        if not match:
            # the perturbed run may diverge only after time t; alignment (s, t) is determined by the past, so it must exist
            ok = False
            continue
        tested += 1
        if match[0][7] != kap or match[0][3] != delta:
            ok = False
    return ok, tested


def analyze(res):
    a, b, L = res['a'], res['b'], res['L']
    T = len(a)
    out = {}
    # cross-covariance sums by L-bin and lag: accumulate sum a_s b_(s+k), sum a_s, sum b_(s+k), count
    for (lo, hi) in LBINS + [(-10 ** 9, 10 ** 9)]:
        mask = (L >= lo) & (L < hi)
        S = np.zeros((KMAX + 1, 4))
        for k in range(KMAX + 1):
            if T - k <= 0:
                break
            m = mask[:T - k]
            aa = a[:T - k][m]
            bb = b[k:][m]
            S[k] = (np.dot(aa, bb), aa.sum(), bb.sum(), m.sum())
        out[(lo, hi)] = S
    # block variances of L over blocks starting at L >= 40
    bv = {}
    for K in (1, 4, 16, 64, 256):
        vals = []
        t = 0
        while t + K < T:
            if L[t] >= 40:
                vals.append(L[t + K] - L[t])
            t += K
        bv[K] = (len(vals), float(np.sum(vals)), float(np.sum(np.square(vals))))
    out['bv'] = bv
    out['T'] = T
    return out


def worker(idx):
    res = run_pair(idx)
    an = analyze(res)
    al = res['aligns']
    return idx, res['rho_ok'], an, al


if __name__ == '__main__':
    t0 = time.time()
    print(f"two-orbit exponent coupling: NPAIRS={NPAIRS} BITS={BITS} SEED={SEED} NPROC={NPROC}", flush=True)
    print("1. exact checks", flush=True)
    with Pool(NPROC) as pool:
        results = pool.map(worker, range(NPAIRS))
    rho_ok = all(r[1] for r in results)
    aligns = [x for r in results for x in r[3]]
    law_ok = all(x[8] for x in aligns)
    check(law_ok, f"ultrametric coupling law b_t = f(a_s, delta) holds at all {len(aligns)} exact alignments")
    conts = [x for x in aligns if x[9] == 1]
    check(all(x[10] for x in conts), f"lockstep recursion delta' = min(delta - a, v_2(3^k - 1)) (tie: delta' > v_2(3^k - 1)) "
                                    f"holds at all {len(conts)} continuation alignments")
    check(rho_ok, "rho_(t+1) = (3/2^a_t) rho_t - 1 + 2^(-L_(t+1)) exactly, first 200 steps of every pair")
    pm = [past_measurability_test(i, 25) for i in range(8)]
    check(all(p[0] for p in pm), f"kappa and delta are unchanged when the tape above the alignment position is replaced "
                                 f"({sum(p[1] for p in pm)} alignments in 8 pairs)")
    steps = sum(r[2]['T'] for r in results)
    print(f"   steps {steps}, alignments {len(aligns)} ({len(aligns)/steps:.3f} per step)  [{time.time()-t0:.0f}s]", flush=True)

    print("2. the law of delta at exact alignments (Haar: P(delta >= d) = 2^(1-d))", flush=True)
    al = np.array([(x[2], min(x[3], 60), x[4], x[5], x[6], x[9]) for x in aligns], dtype=np.int64)
    k_, d_, as_, bt_, Ls_, cont_ = al.T
    for name, m in [("all", np.ones(len(al), bool)), ("k odd", k_ % 2 == 1), ("k = 2 mod 4", k_ % 4 == 2),
                    ("k = 0 mod 4", k_ % 4 == 0), ("L_s in [8,32)", (Ls_ >= 8) & (Ls_ < 32)), ("L_s >= 32", Ls_ >= 32),
                    ("L_s < 8", Ls_ < 8), ("fresh, L>=8", (cont_ == 0) & (Ls_ >= 8)), ("cont., L>=8", (cont_ == 1) & (Ls_ >= 8)),
                    ("fresh k odd L>=8", (cont_ == 0) & (Ls_ >= 8) & (k_ % 2 == 1)), ("fresh k even L>=8", (cont_ == 0) & (Ls_ >= 8) & (k_ % 2 == 0))]:
        n = m.sum()
        if n == 0:
            continue
        P = [float((d_[m] >= d).mean()) for d in range(1, 9)]
        ratio = [P[i] / 2 ** (-i) for i in range(8)]
        cov = float(np.mean((as_[m] - 2) * (bt_[m] - 2)))
        pred = float(np.mean(2 - 3 * 2.0 ** (1 - np.minimum(d_[m], 60))))
        print(f"   {name:14s} n={n:8d}  P(delta>=d)/2^(1-d), d=1..8: " + " ".join(f"{r:.3f}" for r in ratio)
              + f"   Cov(a_s,b_t) = {cov:+.4f} (pred from delta {pred:+.4f})", flush=True)
    # joint law of (a_s, b_t) at alignments with L_s >= 8 against independent Geom(1/2) x Geom(1/2): chi-square on {1..6}^2 cells
    m = Ls_ >= 8
    n = int(m.sum())
    obs = np.zeros((7, 7))
    for i, j in zip(np.minimum(as_[m], 7), np.minimum(bt_[m], 7)):
        obs[i - 1, j - 1] += 1
    pg = np.array([2.0 ** -i for i in range(1, 7)] + [2.0 ** -6])
    exp_ = n * np.outer(pg, pg)
    chi2 = float(np.sum((obs - exp_) ** 2 / exp_))
    print(f"   joint law of (a_s, b_t) at alignments with L_s >= 8 vs independent Geom x Geom: chi^2 = {chi2:.1f} on 48 dof "
          f"(n = {n}); P(b_t = a_s) = {float((as_[m] == bt_[m]).mean()):.5f} (independent: 1/3)", flush=True)
    # E[2^(1-delta)] vs 2/3
    e2 = float(np.mean(2.0 ** (1 - np.minimum(d_, 60))))
    se = float(np.std(2.0 ** (1 - np.minimum(d_, 60))) / math.sqrt(len(d_)))
    print(f"   E[2^(1-delta)] = {e2:.5f} +- {se:.5f}  (Haar value 2/3; aligned covariance = 2 - 3 E = {2-3*e2:+.5f})", flush=True)
    print("   lag distribution at alignment: k / L_s mean = %.3f, median k for L_s in [30,34]: %s" % (
        float(np.mean(k_[Ls_ > 4] / Ls_[Ls_ > 4])), float(np.median(k_[(Ls_ >= 30) & (Ls_ < 34)])) if np.any((Ls_ >= 30) & (Ls_ < 34)) else None),
        flush=True)

    print("3. cross-covariance Cov(a_s, b_(s+k)) by L_s-bin (s.e. ~ 2/sqrt(n))", flush=True)
    for key in LBINS + [(-10 ** 9, 10 ** 9)]:
        S = sum(r[2][key] for r in results)
        cov = np.zeros(KMAX + 1)
        n = S[:, 3]
        with np.errstate(invalid='ignore', divide='ignore'):
            cov = S[:, 0] / n - (S[:, 1] / n) * (S[:, 2] / n)
        lo, hi = key
        name = "all L" if lo < -1000 else f"L in [{lo},{hi})"
        ks = [0, 1, 2, 4, 8, 12, 16, 24, 32, 48, 64, 96, 128, 160]
        print(f"   {name:14s} n0={int(n[0]):8d}  " + " ".join(f"k={k}:{cov[k]:+.4f}" for k in ks if k <= KMAX), flush=True)
        # window around the surfacing lag k ~ L/2
        if lo > 0:
            mid = (lo + hi) // 4
            w = range(max(1, mid - 8), min(KMAX, mid + 9))
            print(f"      window k in [{w.start},{w.stop-1}] around L/2: max |Cov| = {max(abs(cov[k]) for k in w):.4f}, "
                  f"sum = {sum(cov[k] for k in w):+.4f}", flush=True)
    print("4. block variances Var(L_(t+K) - L_t)/K, blocks starting at L >= 40", flush=True)
    for K in (1, 4, 16, 64, 256):
        n = sum(r[2]['bv'][K][0] for r in results)
        s1 = sum(r[2]['bv'][K][1] for r in results)
        s2 = sum(r[2]['bv'][K][2] for r in results)
        if n > 1:
            var = (s2 - s1 * s1 / n) / (n - 1)
            print(f"   K={K:4d}: n={n:7d}  Var/K = {var / K:.4f}", flush=True)
    print(f"[{time.time()-t0:.0f}s]  " + ("ALL CHECKS PASSED" if not FAIL else "SOME CHECK FAILED"), flush=True)
