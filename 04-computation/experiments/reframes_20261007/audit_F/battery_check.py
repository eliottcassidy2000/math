#!/usr/bin/env python3
"""audit_F (E): independent re-verification of section 1 of coalescence_phase_diagram_20261008.md.

  mersenne : recompute (odd count o, Terras length sigma_T) of M_K = 2^K - 1 for all K <= KMAX with a 2^16-ary
             accelerated Terras map (own code) and compare with runcompress_20261007/mersenne_sigma_12800.txt;
             then (1.1) key (o, sigma_T - K) a function of o for K >= 5; identity sigma_T = o log2 3 + log2 n + eps with
             eps = sum log2(1 + 1/(3x)) at K = 101, 1001, 5001; range of eps(M_K); river offset {-o log2 3};
             jump tuning; (1.2) span table; (1.3) jump rates, river births, orphan rate by v2(K-1), braid widths.
  epsgen   : eps(n) for all n <= NMAX (does the stated range (0, 0.33) hold for general n?)
  landing  : (1.4) landing integers N_D of y* = 2/3^D - 1, 3 <= D <= DMAX, computed from the 2-adic truncation of y*
             (parity vector of an integer representative mod 2^W) and the affine parity-vector formula, not the
             author's a/3^m recursion; statistics, the bound N_D < (3/2)^D, value multiplicities.
Usage: python3 battery_check.py mersenne|epsgen|landing [args]
"""
import sys, math
from collections import defaultdict, Counter

DATA = '../../runcompress_20261007/mersenne_sigma_12800.txt'
L3 = math.log2(3)


# ---------------------------------------------------------------- Terras orbit of M_K, accelerated
KB = 16
MASK = (1 << KB) - 1
TAB_C = [0] * (1 << KB)
TAB_D = [0] * (1 << KB)
for b in range(1 << KB):
    x, c = b, 0
    for _ in range(KB):
        if x & 1:
            x = (3 * x + 1) >> 1; c += 1
        else:
            x >>= 1
    TAB_C[b] = c; TAB_D[b] = x
POW3 = [3 ** c for c in range(KB + 1)]


def orbit_counts(n):
    """(odd steps, Terras steps) from n to 1"""
    o = s = 0
    x = n
    lim = 1 << (KB + 1)
    while x >= lim:
        b = x & MASK
        c = TAB_C[b]
        x = POW3[c] * (x >> KB) + TAB_D[b]
        o += c; s += KB
    while x != 1:
        if x & 1:
            x = (3 * x + 1) >> 1; o += 1
        else:
            x >>= 1
        s += 1
    return o, s


def eps_direct(n):
    """sum over odd orbit values x > 1 (before 1) of log2(1 + 1/(3x)), plain iteration"""
    s = 0.0; x = n; o = t = 0
    while x != 1:
        if x & 1:
            s += math.log1p(1 / (3 * x)) / math.log(2); x = (3 * x + 1) >> 1; o += 1
        else:
            x >>= 1
        t += 1
    return s, o, t


def mersenne(kmax):
    import time
    t0 = time.time()
    mine = {}
    for K in range(2, kmax + 1):
        mine[K] = orbit_counts((1 << K) - 1)
    print(f"[mersenne] recomputed (o, sigma_T) for 2 <= K <= {kmax} in {time.time()-t0:.1f}s")
    data = {}
    for line in open(DATA):
        K, o, t = map(int, line.split()); data[K] = (o, t)
    diff = [K for K in range(2, kmax + 1) if data.get(K) != mine[K]]
    print(f"   data file has {len(data)} lines (K = {min(data)}..{max(data)}); mismatches with my recomputation: {len(diff)} {diff[:10]}")
    D = mine
    Ks = sorted(D)
    # (1.1) key (o, sigma_T - K) is a function of o, for K >= 5
    byo = defaultdict(set)
    for K in Ks:
        if K >= 5:
            byo[D[K][0]].add(D[K][1] - K)
    bad = {o: v for o, v in byo.items() if len(v) > 1}
    print(f"(1.1) o-levels (K >= 5): {len(byo)}; levels where sigma_T - K is not constant: {len(bad)}")
    byo_all = defaultdict(set)
    for K in Ks:
        byo_all[D[K][0]].add(D[K][1] - K)
    print(f"      including K = 2,3,4: levels with two keys: { {o: sorted(v) for o, v in byo_all.items() if len(v) > 1} }")
    for K in [k for k in (101, 1001, 5001) if k <= kmax]:
        e, o, t = eps_direct((1 << K) - 1)
        assert (o, t) == D[K]
        lhs = t - o * L3 - (K + math.log2(1 - 2.0 ** (-K)))
        print(f"      K={K}: sigma_T - o log2 3 - log2 M_K = {lhs:.7f}; sum log2(1+1/(3x)) = {e:.7f}; diff {lhs-e:.1e}")
    eps = {K: D[K][1] - D[K][0] * L3 - (K + math.log2(1 - 2.0 ** (-K))) for K in Ks}
    ev = sorted(eps.values())
    print(f"      eps(M_K), 2 <= K <= {kmax}: range [{ev[0]:.4f}, {ev[-1]:.4f}], median {ev[len(ev)//2]:.4f}")
    off = max(abs(((eps[K] + D[K][0] * L3) % 1.0 + 0.5) % 1.0 - 0.5) for K in Ks if K >= 5)
    off2 = max(abs(((eps[K] + D[K][0] * L3 + math.log2(1 - 2.0 ** (-K))) % 1.0 + 0.5) % 1.0 - 0.5) for K in Ks if K >= 5)
    off3 = max(abs(eps[K] - ((-D[K][0] * L3) % 1.0)) for K in Ks if K >= 20)
    print(f"      river offset: max_K>=5 dist(eps + o log2 3, Z) = {off:.2e} (= |log2(1-2^-5)|); after removing log2(1-2^-K): {off2:.1e}; "
          f"max_K>=20 |eps - {{-o log2 3}}| = {off3:.1e}")
    # jumps and tuning
    jumps = [(K, D[K][0] - D[K-1][0]) for K in Ks if K - 1 in D and D[K][0] != D[K-1][0]]
    def dist(x): return abs(x - round(x))
    dj = sorted(dist(abs(j) * L3) for K, j in jumps)
    deps = sorted(abs(eps[K] - eps[K-1]) for K, j in jumps)
    print(f"      jumps: {len(jumps)}; median ||dO log2 3|| = {dj[len(dj)//2]:.4f}; median |d eps| over jumps = {deps[len(deps)//2]:.4f}; "
          f"max ||dO log2 3|| = {dj[-1]:.4f} vs max eps - min eps = {ev[-1]-ev[0]:.4f}")
    big = [K for K in Ks if K >= min(1000, kmax // 2)]
    eb = sorted(eps[K] for K in big)
    print(f"      eps(M_K) for K >= {big[0]}: quartiles {eb[len(eb)//4]:.4f}, {eb[len(eb)//2]:.4f}, {eb[3*len(eb)//4]:.4f}")
    # (1.2) rivers = level sets of o; spans
    riv = defaultdict(list)
    for K in Ks:
        riv[D[K][0]].append(K)
    rows = sorted((max(v), min(v), len(v)) for v in riv.values())
    print(f"(1.2) rivers (o-level sets, all K >= 2): {len(riv)}")
    for lo, hi in ((100, 400), (400, 1600), (1600, 6400), (6400, kmax + 1)):
        sp = [(a - b) / math.sqrt(a) for a, b, n in rows if lo <= a < hi and n > 1]
        if sp:
            m = sum(sp) / len(sp)
            se = math.sqrt(sum((x - m) ** 2 for x in sp) / (len(sp) - 1) / len(sp)) if len(sp) > 1 else float('nan')
            print(f"      maxK in [{lo},{hi}): {len(sp):3d} rivers, span/sqrt(maxK) mean {m:.2f} +- {se:.2f}, max {max(sp):.2f}")
    by_span = sorted(rows, key=lambda r: -(r[0] - r[1]))[:3]
    by_norm = sorted(rows, key=lambda r: -(r[0] - r[1]) / math.sqrt(r[0]))[:3]
    by_size = sorted(rows, key=lambda r: -r[2])[:3]
    print(f"      widest by span: {[(b, a, n, a-b, round((a-b)/math.sqrt(a), 2)) for a, b, n in by_span]}")
    print(f"      widest by span/sqrt(maxK): {[(b, a, n, round((a-b)/math.sqrt(a), 2)) for a, b, n in by_norm]}")
    print(f"      largest by members: {[(b, a, n) for a, b, n in by_size]}")
    alive_end = [(a, b, n) for a, b, n in rows if a > kmax - 200]
    print(f"      rivers with a member in (KMAX-200, KMAX] (censored spans): {len(alive_end)}")
    # (1.3) jump rates, births, orphans, braid widths
    jK = [K for K, j in jumps]
    for lo, hi in ((100, 400), (400, 1600), (1600, 6400), (6400, kmax + 1)):
        n = sum(1 for K in jK if lo <= K < hi); r = n / (hi - lo)
        print(f"(1.3) jumps in [{lo},{hi}): {n:4d}, rate {r:.4f}, rate*sqrt(Kmid) {r*math.sqrt((lo+hi)/2):.2f}")
    first = {}
    for K in Ks:
        first.setdefault(D[K][0], K)
    starts = sorted(first.values())
    prev = None
    for X in (400, 1600, 6400, 12800):
        R = sum(1 for s in starts if s <= X)
        loc = '' if prev is None else f", local exponent {math.log(R/prev[1])/math.log(X/prev[0]):.3f}"
        print(f"      births R({X}) = {R}; R/X^0.41 = {R/X**0.41:.2f}{loc}")
        prev = (X, R)
    def v2(x): return (x & -x).bit_length() - 1
    seen = set(); rows2 = []
    for K in Ks:
        o = D[K][0]
        orphan = o not in seen
        if K % 2 == 1 and K >= 200:
            rows2.append((min(v2(K - 1), 5), orphan))
        seen.add(o)
    c = Counter(m for m, _ in rows2); co = Counter(m for m, x in rows2 if x)
    tot_o = sum(co.values()); tot = sum(c.values()); p0 = tot_o / tot
    chi2 = sum((co[m] - c[m] * p0) ** 2 / (c[m] * p0) for m in c)
    print(f"      orphan rate by m = v2(K-1) (odd K >= 200): " + ", ".join(f"m={m}: {co[m]}/{c[m]} = {co[m]/c[m]:.4f}" for m in sorted(c))
          + f"; chi2 (4 dof) = {chi2:.2f}")
    births_even = sum(1 for s in starts if s % 2 == 0 and s >= 200)
    births_odd = sum(1 for s in starts if s % 2 == 1 and s >= 200)
    print(f"      births at K >= 200: odd K {births_odd}, even K {births_even}")
    for W in (50, 200, 1000):
        widths = [len({D[K][0] for K in range(a, a + W)}) for a in range(1000, kmax - W, W)]
        if widths:
            print(f"      braid width W={W}: mean {sum(widths)/len(widths):.2f}, max {max(widths)}")


def epsgen(nmax):
    """eps(n), o(n), sigma(n) for all n <= nmax by memoised descent (values of T may exceed nmax: iterate until below n)"""
    eps = [0.0] * (nmax + 1)
    best = (0.0, 1)
    LN2 = math.log(2)
    for n in range(2, nmax + 1):
        x, s = n, 0.0
        while x >= n:
            if x & 1:
                s += math.log1p(1 / (3 * x)) / LN2; x = (3 * x + 1) >> 1
            else:
                x >>= 1
        eps[n] = s + eps[x]
        if eps[n] > best[0]:
            best = (eps[n], n)
    over = [n for n in range(2, nmax + 1) if eps[n] >= 0.33]
    print(f"[epsgen] eps(n) for 2 <= n <= {nmax}: max {best[0]:.4f} at n = {best[1]}; #n with eps >= 0.33: {len(over)}; "
          f"first few {over[:12]}")
    for n in (3, 7, 9, 27, best[1]):
        e, o, t = eps_direct(n)
        print(f"      n={n}: eps = {e:.4f} (o={o}, sigma_T={t}); check sigma_T - o log2 3 - log2 n = {t - o*L3 - math.log2(n):.4f}")


def landing(dmax):
    import time
    t0 = time.time()
    res = {}
    for D in range(1, dmax + 1):
        P = 3 ** D
        W = 4 * D + 80
        mod = 1 << W
        ystar = ((2 - P) * pow(P, -1, mod)) % mod      # integer representative of 2/3^D - 1 modulo 2^W
        x = ystar; rho = 0; o = 0; s = 0
        while o < D:
            if x & 1:
                rho = 3 * rho + (1 << s); x = (3 * x + 1) >> 1; o += 1
            else:
                x >>= 1
            s += 1
            assert s < W - 2, "window too small"
        num = (2 - P) + rho           # 3^D y* + rho
        assert num % (1 << s) == 0
        N = num >> s
        res[D] = (N, s)
    print(f"[landing] N_D for 1 <= D <= {dmax} via 2-adic truncation + parity-vector formula: {time.time()-t0:.1f}s")
    rng = [D for D in range(3, dmax + 1)]
    Ns = [res[D][0] for D in rng]
    S0 = [res[D][1] / D for D in rng]
    print(f"   3 <= D <= {dmax}: nonpositive {sum(1 for N in Ns if N <= 0)}; max {max(Ns)} (at D = {[D for D in rng if res[D][0] == max(Ns)][:5]}); "
          f"median {sorted(Ns)[len(Ns)//2]}; s_0/D in [{min(S0):.3f}, {max(S0):.3f}]")
    small = [(D, res[D][0]) for D in range(1, dmax + 1) if not (2 ** D) * res[D][0] < 3 ** D]   # N_D < (3/2)^D, exact
    print(f"   bound N_D < (3/2)^D violated for D in {small[:10]} (all D in 1..{dmax}); N_D for D=1..12: {[res[D][0] for D in range(1, 13)]}")
    for n in (2, 4, 8, 16, 32, 64, 128, 256, 512, 1024):
        p = sum(1 for N in Ns if N > n) / len(Ns)
        print(f"      P(N_D > {n:5d}) = {p:.4f}   n*P = {n*p:.3f}   count {sum(1 for N in Ns if N > n)}")
    cnt = Counter(Ns)
    print(f"   distinct values {len(cnt)} among {len(Ns)}; most common {cnt.most_common(12)}")
    top = sorted(cnt.items(), key=lambda kv: -kv[0])[:12]
    print(f"   largest values with multiplicity: {top}")
    same_next = sum(1 for D in range(3, dmax) if res[D][0] == res[D + 1][0])
    print(f"   N_D = N_(D+1) for {same_next} of {dmax - 3} consecutive pairs")
    big = [D for D in rng if res[D][0] > 512]
    print(f"   D with N_D > 512: {len(big)}; first 30: {big[:30]}")
    # rivers: classes of equal landing key (N_D, s_0); consecutive sharing
    for lo, hi in ((3, 4001), (3, 4500)):
        if hi - 1 > dmax:
            continue
        keys = {}
        for D in range(lo, hi):
            keys.setdefault(res[D], []).append(D)
        share = sum(1 for D in range(lo, hi - 1) if res[D] == res[D + 1])
        sizes = sorted((len(v) for v in keys.values()), reverse=True)
        Nr = [k[0] for k in keys]
        print(f"   D in [{lo},{hi}): {len(keys)} distinct landings (N_D, s_0); consecutive pairs sharing: {share}/{hi - 1 - lo}; "
              f"largest classes {sizes[:8]}")
        recs = sorted(((k[0], min(v), max(v), len(v)) for k, v in keys.items()), reverse=True)[:6]
        print(f"      records (N, minD, maxD, size): {recs}")
        for n in (4, 16, 64, 256):
            print(f"      per-class P(N > {n:3d}) = {sum(1 for x in Nr if x > n)/len(Nr):.4f}   (2.87/n = {2.87/n:.4f})")
    # odd-step excess of the orbits of small N (HYP-9240 / note: 'excess at most 9 for N <= 880')
    best = (-1e9, None)
    for N in range(1, 881):
        x, o, t = N, 0, 0
        while x != 1:
            if x & 1:
                x = (3 * x + 1) >> 1; o += 1
            else:
                x >>= 1
            t += 1
        if o - t / 2 > best[0]:
            best = (o - t / 2, N)
    print(f"   max over 1 <= N <= 880 of odd(N) - sigma_T(N)/2 = {best[0]} at N = {best[1]}")


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'mersenne':
        mersenne(int(sys.argv[2]) if len(sys.argv) > 2 else 12800)
    elif cmd == 'epsgen':
        epsgen(int(sys.argv[2]) if len(sys.argv) > 2 else 10 ** 6)
    elif cmd == 'landing':
        landing(int(sys.argv[2]) if len(sys.argv) > 2 else 4000)
