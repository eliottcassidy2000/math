#!/usr/bin/env python3
"""
lrc14_moment_layer_forensics_klein_S213.py

HYP-5770: which support layer carries the deviation from iid?

KEY IDENTITY (one histogram pass gives every moment layer EXACTLY):
    L_d(q) := Σ_{|T|=d} N_T(q) = Σ_{p≠0} C(C(p), d)
where C(p) = # killing classes at p. The B_D bounds are alternating sums of
the L_d. The iid reference per layer is C(ncl, d)·κ_q^d·(q−1), κ_q = (2c−1)/q.

MEASURE, per instance (generic corpus + @91 + the V=260 near-dilation):
 [1] per-layer relative deviation dev_d(q) = L_d(q)/(q−1) − C(ncl,d)κ^d,
     its distribution over q (mean, |mean|, max) — WHERE does the deviation live?
 [2] the SIGNED aggregated deviation per layer (mean over q) vs the ABSOLUTE
     (mean |dev|): the cancellation ratio per layer — is d=5 self-cancelling?
 [3] what B5 needs: B5 = (q−1)[Σ_d (−1)^d C κ^d] + Σ_d (−1)^d dev_d·(q−1);
     per-layer contribution to the B5 deficit; which layers must be controlled
     TWO-SIDED (d=2,4: lower; d=1,3,5: upper) and how big their actual
     deviations are vs the 0.1221 slack.
 [4] low-height relation counts per support size (h ≤ 3): the per-layer
     resonance census — correlate with per-layer |dev|.
"""

from math import gcd, comb
import random

random.seed(20260713)
QS = list(range(2, 15))


def is_covering(S):
    return all(any(s % q == 0 for s in S) for q in QS)


def coverage_hist(S, q):
    def r_safe(r):
        return 14 * r >= q and 14 * (q - r) >= q
    cov = bytearray(q)
    seen = set()
    ncl = 0
    for v in S:
        r = v % q
        key = min(r, (q - r) % q)
        if key in seen:
            continue
        seen.add(key)
        ncl += 1
        if r == 0:
            for p in range(q):
                cov[p] += 1
            continue
        g = gcd(r, q)
        rr, qq = r // g, q // g
        inv = pow(rr, -1, qq)
        for m in range(qq):
            s = m * g
            if not r_safe(s):
                p0 = (m * inv) % qq
                for t in range(g):
                    cov[p0 + t * qq] += 1
    hist = [0] * (ncl + 1)
    for p in range(1, q):
        hist[cov[p]] += 1
    return hist, ncl


def layer_devs(S):
    """per-layer signed devs dev_d(q) for d=1..5, arrays over q in (V,2V]."""
    V = max(S)
    out = {d: [] for d in range(1, 6)}
    b5s = []
    for q in range(V + 1, 2 * V + 1):
        hist, ncl = coverage_hist(S, q)
        c = -(-q // 14)
        kq = (2 * c - 1) / q
        for d in range(1, 6):
            Ld = sum(n * comb(cc, d) for cc, n in enumerate(hist))
            iid = comb(ncl, d) * kq ** d * (q - 1)
            out[d].append((Ld - iid) / (q - 1))
        b5 = sum(n * sum((-1) ** dd * comb(cc, dd) for dd in range(6))
                 for cc, n in enumerate(hist))
        b5s.append(b5 / (q - 1))
    return out, b5s


def relation_census(S, H=3):
    """# exact relations by support size (both signs), height<=H."""
    from itertools import combinations, product
    cnt = {s: 0 for s in range(2, 6)}
    coeffs = [c for c in range(-H, H + 1) if c != 0]
    for s in range(2, 6):
        for U in combinations(range(13), s):
            vs = [S[i] for i in U]
            for j in product(coeffs, repeat=s):
                if sum(a * b for a, b in zip(j, vs)) == 0:
                    cnt[s] += 1
    return cnt


def main():
    print("=" * 78)
    print("HYP-5770: moment-layer forensics (klein-S213)")
    print("=" * 78)

    corpus = [
        ([12, 33, 46, 47, 68, 73, 79, 81, 85, 87, 91, 112, 120], "adv-worst-120"),
        ([62, 66, 69, 102, 109, 118, 120, 126, 130, 136, 159, 185, 200], "adv-200"),
        ([10, 46, 48, 59, 66, 148, 177, 181, 208, 213, 236, 261, 280], "adv-280"),
        ([9, 16, 24, 33, 40, 47, 54, 62, 65, 70, 77, 84, 91], "@91-comb"),
        ([20, 41, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260], "near-dil-260"),
    ]
    for S, name in corpus:
        devs, b5s = layer_devs(S)
        rc = relation_census(S)
        V = max(S)
        print(f"\n{name} (V={V}): avg B5/q = {sum(b5s)/len(b5s):+.4f}; "
              f"#relations(H<=3) by support: {rc}")
        print(f"  {'d':>2} {'iid C(13,d)/7^d':>15} {'mean dev':>10} {'mean |dev|':>10} "
              f"{'max |dev|':>10} {'cancel ratio':>12}")
        for d in range(1, 6):
            xs = devs[d]
            iid = comb(13, d) / 7 ** d
            m = sum(xs) / len(xs)
            ma = sum(abs(x) for x in xs) / len(xs)
            mx = max(abs(x) for x in xs)
            print(f"  {d:>2} {iid:>15.4f} {m:>+10.4f} {ma:>10.4f} {mx:>10.4f} "
                  f"{(abs(m)/ma if ma else 0):>12.3f}")
        # per-layer contribution to the B5 deficit (signed, aggregated)
        contrib = [(-1) ** d * (sum(devs[d]) / len(devs[d])) for d in range(1, 6)]
        print(f"  aggregated B5-deficit decomposition by layer (d=1..5): "
              + ", ".join(f"{-c:+.4f}" for c in contrib)
              + f"  | total {-sum(contrib):+.4f}")

    print("""
READING GUIDE: 'cancel ratio' = |mean signed| / mean |abs| per layer — near 0
means the layer's deviations SELF-CANCEL across moduli (aggregation wins);
near 1 means coherent bias (must be boxed exactly). The B5-deficit
decomposition shows which layers a proof must control and in which direction
(d odd: upper bounds suffice; d even: lower bounds — the hard direction).""")
    print("DONE.")


if __name__ == '__main__':
    main()
