# Independent audit of THM-4887, part 2: the strip lemma, the small-denominator count, the decomposition of
# E_n(Delta), all ordered pairs h != h', worst-case shift D.  Written from scratch (no session imports).
#
# For (n, m): C = 2^m, P = 2^(n-m-1), B = P/C (must be >= 8), h, h' in [P, 2P).
#  Z*(h,h';W) = max over integer D of #{(c,c') in [0,C)^2 : |h c - h' c' - D| <= W}, for W = C (Delta = 0) and
#  W = 2C (Delta = C).
#  beta(h,h'): the proof allows ANY Dirichlet pair (a, b), 1 <= b <= B, |b h/h' - a| < 1/B; its reduced denominator.
#  We test the lemma with the LARGEST admissible reduced denominator (strongest form, smallest bound) and also
#  record the min-b choice.  Lemma: Z*(W=2C) <= 4C^2/P + 3C/beta + 7, Z*(W=C) <= 2C^2/P + 3C/beta + 5.
#  Step 4: N_beta := #{pairs with beta(h,h') = beta} <= P (beta+1)(4C/beta+1) (checked for both choice rules),
#          sum 1/beta <= P (4C ln B + 10.58 C + B + 1 + ln B).
#  Decomposition: E_n(Delta) = P C + sum_{h != h'} actual(h,h') with actual <= Z*(W = C + Delta); same-h pairs only
#  diagonal (min same-h gap > C).
import sys, math, time, random
import numpy as np

def g_level_arrays(nmax):
    """g on X_j = [2^(j-1), 2^j) for j = 1..nmax, built level by level from g(2v+b)."""
    levels = {1: np.array([1], dtype=np.int64)}   # g(1) = 1
    for j in range(1, nmax):
        prev = levels[j]
        u = np.arange(1 << (j - 1), 1 << j, dtype=np.int64)
        nxt = np.empty(2 * len(prev), dtype=np.int64)
        nxt[0::2] = prev + u + (u & 1)    # v = 2u   : g(2u) = g(u) + u + (u mod 2)
        nxt[1::2] = prev                  # v = 2u+1 : g(2u+1) = g(u)
        levels[j + 1] = nxt
    return levels

def g_small(v):
    G, A = 1, 0
    for ch in bin(v)[2:]:
        if ch == "1": A = 2 * A + 1
        else: G += A + (A & 1); A *= 2
    return G

def beta_choices(h, hps, B):
    """For fixed h and vector hps, return (beta for the minimal admissible b, max admissible reduced denominator)."""
    bmin = np.zeros(len(hps), dtype=np.int64)
    bmax = np.zeros(len(hps), dtype=np.int64)
    for b in range(1, B + 1):
        a = (2 * b * h + hps) // (2 * hps)           # nearest integer to b h / h'
        diff = np.abs(b * h - a * hps)
        ok = B * diff < hps                          # |b h/h' - a| < 1/B
        red = b // np.gcd(a, b)
        newmin = ok & (bmin == 0)
        bmin[newmin] = red[newmin]
        bmax = np.where(ok, np.maximum(bmax, red), bmax)
    assert (bmin > 0).all() and (bmax > 0).all()
    return bmin, bmax

def run(n, m, levels, sample_pairs=None, seed=1):
    C = 1 << m; k = n - m; P = 1 << (k - 1); B = P // C
    assert P >= 8 * C, (n, m)
    Gn = levels[n]                                   # g on X_n, index v - 2^(n-1)
    Gb = Gn.reshape(P, C)                            # row h - P, column l
    cvec = np.array([sum((1 - ((l >> i) & 1)) << (m - 1 - i) for i in range(m)) for l in range(C)], dtype=np.int64)
    assert sorted(cvec.tolist()) == list(range(C))
    # identity spot check on this block: g(h 2^m + l) = g(h) + h c + gamma
    gh = levels[k]
    hs = np.arange(P, 2 * P, dtype=np.int64)
    gl = np.array([g_small(l) for l in range(C)], dtype=np.int64)
    gam = gl[None, :] - 1 + (hs[:, None] & 1) * (np.arange(C)[None, :] < C // 2)
    assert (Gb == gh[:, None] + hs[:, None] * cvec[None, :] + gam).all()
    assert gam.min() >= 0 and gam.max() <= C - 1
    # same-h separation
    srt = np.sort(Gb, axis=1)
    same_gap = int(np.diff(srt, axis=1).min())
    cc = np.arange(C, dtype=np.int64)
    pairs_h = range(P) if sample_pairs is None else None
    stats = dict(viol1=0, viol2=0, r1=0.0, r2=0.0, viol1_minb=0, viol2_minb=0, S1=0, S2=0, Sact0=0, SactC=0,
                 act_gt_Z0=0, act_gt_ZC=0, invb_min=0.0, invb_max=0.0, npairs=0)
    Nb_min = np.zeros(B + 1, dtype=np.int64); Nb_max = np.zeros(B + 1, dtype=np.int64)
    W1, W2 = C, 2 * C
    BIG = 4 * P * C + 8 * C + 10
    def process(i, js):
        h = P + i
        hps = P + js
        # worst-case strip counts
        vals = (h * cc[:, None] - hps[:, None, None] * cc[None, None, :] * 0 - 0)  # placeholder (overwritten)
        vals = (h * cc)[None, :, None] - (hps[:, None, None] * cc[None, None, :])  # shape (len(js), C, C)
        vals = np.sort(vals.reshape(len(js), C * C), axis=1)
        flat = (vals + (np.arange(len(js), dtype=np.int64) * BIG)[:, None]).ravel()
        base = np.arange(flat.size, dtype=np.int64)
        z1 = (np.searchsorted(flat, flat + 2 * W1, side="right") - base).reshape(len(js), C * C).max(axis=1)
        z2 = (np.searchsorted(flat, flat + 2 * W2, side="right") - base).reshape(len(js), C * C).max(axis=1)
        bmin, bmax = beta_choices(h, hps, B)
        b1 = 2 * C * C / P + 3 * C / bmax + 5; b2 = 4 * C * C / P + 3 * C / bmax + 7
        b1m = 2 * C * C / P + 3 * C / bmin + 5; b2m = 4 * C * C / P + 3 * C / bmin + 7
        stats["viol1"] += int((z1 > b1).sum()); stats["viol2"] += int((z2 > b2).sum())
        stats["viol1_minb"] += int((z1 > b1m).sum()); stats["viol2_minb"] += int((z2 > b2m).sum())
        stats["r1"] = max(stats["r1"], float((z1 / b1).max())); stats["r2"] = max(stats["r2"], float((z2 / b2).max()))
        stats["S1"] += int(z1.sum()); stats["S2"] += int(z2.sum())
        np.add.at(Nb_min, bmin, 1); np.add.at(Nb_max, bmax, 1)
        stats["invb_min"] += float((1.0 / bmin).sum()); stats["invb_max"] += float((1.0 / bmax).sum())
        stats["npairs"] += len(js)
        # actual pair counts for Delta = 0 and Delta = C
        d = np.abs(Gb[i][None, :, None] - Gb[js][:, None, :])     # (len(js), C, C)
        a0 = (d == 0).reshape(len(js), -1).sum(axis=1); aC = (d <= C).reshape(len(js), -1).sum(axis=1)
        stats["Sact0"] += int(a0.sum()); stats["SactC"] += int(aC.sum())
        stats["act_gt_Z0"] += int((a0 > z1).sum()); stats["act_gt_ZC"] += int((aC > z2).sum())
    t0 = time.time()
    if sample_pairs is None:
        allj = np.arange(P, dtype=np.int64)
        for i in range(P):
            process(i, allj[allj != i])
    else:
        rng = np.random.default_rng(seed)
        per = max(1, sample_pairs // P)
        for i in range(P):
            js = rng.choice(P - 1, size=min(per, P - 1), replace=False)
            js = js + (js >= i)
            process(i, js.astype(np.int64))
    # step-4 counts
    nb_viol_min = sum(1 for b in range(1, B + 1) if Nb_min[b] > P * (b + 1) * (4 * C / b + 1))
    nb_viol_max = sum(1 for b in range(1, B + 1) if Nb_max[b] > P * (b + 1) * (4 * C / b + 1))
    nb_ratio = max(Nb_max[b] / (P * (b + 1) * (4 * C / b + 1)) for b in range(1, B + 1))
    invb_bound = P * (4 * C * math.log(B) + 10.58 * C + B + 1 + math.log(B))
    # direct energies
    x = np.sort(Gn)
    _, cnt = np.unique(x, return_counts=True)
    E0 = int((cnt.astype(np.int64) ** 2).sum())
    EC = int((np.searchsorted(x, x + C, side="right") - np.searchsorted(x, x - C, side="left")).sum())
    lB = math.log(B)
    D1_0 = P * C + 8 * P * P + P * C * C * (33.74 + 12 * lB) + 3 * P * C * (1 + lB)
    D1_C = P * C + 10 * P * P + P * C * C * (35.74 + 12 * lB) + 3 * P * C * (1 + lB)
    tag = "ALL pairs" if sample_pairs is None else f"{stats['npairs']} random pairs"
    print(f"(n,m)=({n},{m}) P={P} C={C} B={B} [{tag}, {time.time()-t0:.0f}s]")
    print(f"   strip lemma, beta = largest admissible reduced denominator: W=C viol {stats['viol1']} "
          f"(max Z*/bound {stats['r1']:.3f}); W=2C viol {stats['viol2']} (max Z*/bound {stats['r2']:.3f}); "
          f"with min-b beta: viol {stats['viol1_minb']}, {stats['viol2_minb']}")
    if sample_pairs is None:
        print(f"   step 4: N_beta <= P(beta+1)(4C/beta+1) violations: min-b rule {nb_viol_min}, max-beta rule "
              f"{nb_viol_max} (max ratio {nb_ratio:.3f}); sum 1/beta = {stats['invb_min']:.1f} (min-b), "
              f"{stats['invb_max']:.1f} (max-beta) <= {invb_bound:.1f}")
        print(f"   same-h min gap {same_gap} (> Delta_max = C = {C}; P - C = {P - C})")
        print(f"   E_n(0) = {E0} = PC + sum actual = {P*C + stats['Sact0']} ; <= PC + sum Z*(W=C) = {P*C + stats['S1']}"
              f" <= (D1, Delta=0) {D1_0:.0f};  actual > Z* in {stats['act_gt_Z0']} pairs")
        print(f"   E_n(C) = {EC} = PC + sum actual = {P*C + stats['SactC']} ; <= PC + sum Z*(W=2C) = {P*C + stats['S2']}"
              f" <= (D1, Delta=C) {D1_C:.0f};  actual > Z* in {stats['act_gt_ZC']} pairs")
        assert E0 == P * C + stats["Sact0"] and EC == P * C + stats["SactC"]
    sys.stdout.flush()

if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    levels = g_level_arrays(22)
    if which in ("all", "small"):
        for (n, m) in [(10, 2), (12, 4), (13, 4), (14, 4), (13, 3), (15, 5)]:
            run(n, m, levels)
    if which in ("all", "big"):
        run(16, 5, levels)
        run(18, 6, levels, sample_pairs=300000)
        run(21, 7, levels, sample_pairs=100000)
