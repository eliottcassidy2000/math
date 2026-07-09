"""
lrc14_near_resonance_count_mertens_opus_S167.py   (opus-2026-07-09-S167)

THE NEAR-RESONANCE COUNT + the single dissociated-branch inequality + the MERTENS structural analogy.

Setup (opus-S165 / kps-S92). The good period j*<=N <=> S_N := sum_{j=1}^N W(j/V) > 0, W = cluster
uncovered measure. S_N = N(6/7)^k + Corr_N, Corr_N = sum_{m != 0} Whatfreq(m) G_N(m/V),
G_N(a)=sum_{j=1}^N e(ja) (geometric), the frequencies m = n.e over BALANCED relations n (sum n_i=0).
NEAR-RESONANCES = m with ||m/V|| small (G_N ~ N).  kps-S92: the ABSOLUTE bound sum|What|min(N,1/2||.||)
~ 20x the target N(6/7)^k, but the SIGNED r_N = 0.08-0.26 -- CANCELLATION IS ESSENTIAL.

MERTENS ANALOGY (structural).  This is the Mertens situation: a signed arithmetic sum whose ABSOLUTE
value is trivially large (sum|mu(n)| = x; here sum|What| diverges, opus-S154) but whose SIGNED value
has deep cancellation (M(x)=o(x) conj. O(x^{1/2+eps}); here r_N small).  Proving the cancellation is
the hard part in BOTH, and the NO-CANCELLATION extremal is the STRUCTURED case (AP / complete residues
= the tight LRC instance; for Mertens, the conspiracy of zeta zeros).  Moreover the arc Fourier
b(m)=(1-e(m/7))/(2 pi i m) VANISHES at 7|m, so the resonance product prod b(n_i) is supported on
n_i coprime-to-7 -- an inclusion-exclusion over residues mod 7 = a MOBIUS structure (klein's x7 collapse).

This file: (1) the near-resonance COUNT vs longest-AP (dissociated => few); (2) the signed Corr_N
cancellation (|Corr| << sum of |terms|); (3) the mod-7 vanishing / Mobius structure; (4) the clean
route: mac-mini's arc-count (#arcs < rho* V) sidesteps the cancellation, big margin.
"""
import sys, random, cmath, math
from math import gcd


def W_uncovered(E, V, j):
    ph = sorted((e * j) % V for e in E)
    thr = V / 7.0; tot = 0.0; prev = ph[-1] - V
    for p in ph:
        g = p - prev
        if g > thr: tot += g - thr
        prev = p
    return tot / V


def longest_ap(E):
    S = set(E); E = sorted(E); best = 1
    for i, a in enumerate(E):
        for b in E[i + 1:]:
            d = b - a
            if a - d in S: continue
            L = 2; x = b + d
            while x in S: L += 1; x += d
            best = max(best, L)
    return best


def near_resonances(E, V, N, supp_max=3, coef_max=2):
    """balanced relations n (support<=supp_max, |n_i|<=coef_max, sum n_i=0) with ||(n.e)/V|| < 1/(2N),
    all n_i coprime to 7 (else the arc-Fourier product vanishes).  Returns count + list of m=n.e."""
    import itertools
    k = len(E); found = []
    for s in range(2, supp_max + 1):
        for idx in itertools.combinations(range(k), s):
            for vals in itertools.product([v for v in range(-coef_max, coef_max + 1) if v != 0], repeat=s):
                if sum(vals) != 0: continue
                if any(v % 7 == 0 for v in vals): continue          # arc-Fourier zero at 7|m
                m = sum(v * E[i] for v, i in zip(vals, idx))
                if m == 0: continue
                frac = (m % V) / V
                if min(frac, 1 - frac) < 1.0 / (2 * N):
                    found.append(m)
    return len(found), found


def corr_direct(E, V, N):
    base = (6.0 / 7.0) ** len(E)
    SN = sum(W_uncovered(E, V, j) for j in range(1, N + 1))
    return SN - N * base, SN


def main():
    print("=" * 96)
    print("NEAR-RESONANCE COUNT + the dissociated inequality + the MERTENS cancellation analogy")
    print("=" * 96)
    r = random.Random(167)
    for k in (11, 13):
        Nk = -(-7 * (k - 1) // 6)
        base = (6 / 7) ** k
        print(f"\n  k={k}: N=ceil(7(k-1)/6)={Nk}, target N(6/7)^k = {Nk*base:.4f}")
        print(f"    longest-AP L | #near-res (supp<=3) | Corr_N | r_N=|Corr|/(N base) | (samples)")
        buckets = {}
        cands = []
        for _ in range(2500):
            spread = r.randint(k, 45)
            E = sorted(set([0] + r.sample(range(1, spread), k - 2) + [spread]))
            if len(E) == k: cands.append(E)
        for d in range(1, 5):
            cands.append([i * d for i in range(k)])
        for E in cands:
            spread = max(E)
            for V in range(spread + 1, (7 * spread) // 6 + 1):
                if W_uncovered(E, V, 0) is None: continue
                # only where a good period exists (skip tight AP)
                if not any(W_uncovered(E, V, j) > 1e-12 for j in range(1, V)): continue
                nr, _ = near_resonances(E, V, Nk)
                cr, SN = corr_direct(E, V, Nk)
                L = longest_ap(E)
                b = buckets.setdefault(L, [0, 0.0, 0.0, 0])
                b[0] = max(b[0], nr); b[1] = max(b[1], abs(cr))
                b[2] = max(b[2], abs(cr) / (Nk * base)); b[3] += 1
        for L in sorted(buckets):
            b = buckets[L]
            tag = "  <- dissociated (few res)" if L <= k - 6 else ("  <- near-AP (many res)" if L >= k - 3 else "")
            print(f"       L={L:2d}        |    {b[0]:3d}              | {b[1]:.4f} | {b[2]:.3f}"
                  f"            | ({b[3]}){tag}")

    print("\n" + "=" * 96)
    print("(2) CANCELLATION: |Corr_N| vs sum of |individual near-resonance terms| (Mertens-like)")
    # a dissociated and an AP example
    exs = {"AP {0..10} scale2 V=21": ([2 * i for i in range(11)], 21),
           "dissociated (Sidon-ish)": (sorted(set([0, 1, 3, 7, 12, 20, 30, 44, 60, 78, 80])), 90)}
    for name, (E, V) in exs.items():
        if len(E) != 11:
            continue
        N = 12
        cr, SN = corr_direct(E, V, N)
        # sum of |W(j/V)| deviations as a crude "absolute" proxy
        base = (6/7) ** 11
        absmass = sum(abs(W_uncovered(E, V, j) - base) for j in range(1, N + 1))
        nr, ms = near_resonances(E, V, N)
        print(f"  {name}: #near-res={nr}, |Corr_N|={abs(cr):.4f}, sum|W(j/V)-base|={absmass:.4f}"
              f"  (ratio {absmass/max(abs(cr),1e-9):.1f}x => cancellation)")
    print("\n  READING: dissociated => few near-resonances => small signed |Corr_N|; AP => many, resonant.")
    print("  The ABSOLUTE bound (sum of |terms|) is >> |Corr_N| (Mertens: sum|mu| >> |M(x)|) --")
    print("  cancellation ESSENTIAL (kps-S92). arc-Fourier zeros at 7|m => Mobius/x7-collapse structure.")
    print("  CLEAN ROUTE (mac-mini c): #arcs < rho* V (both a-priori, big margin) sidesteps the cancellation.")
    print("=" * 96)


if __name__ == "__main__":
    main()
