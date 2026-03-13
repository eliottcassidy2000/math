#!/usr/bin/env python3
"""
walsh_paley_alignment.py -- The Paley Perfect Alignment Theorem

DISCOVERY: At p=11 (p=3 mod 4), the Paley orientation sigma_P has the property
that ALL Walsh contributions at EVERY even degree are POSITIVE.

sigma_P = (1, -1, 1, 1, 1)  [+1 for QR gaps, -1 for NQR gaps]

At the Paley point:
  deg-0: +93101.25  (always positive)
  deg-2: +1402.50   (ALL 10 terms positive after sign flip)
  deg-4: +591.25    (ALL 5 terms positive after sign flip)
  Total: 95095 = H(Paley) CHECK

At the Interval point (all +1):
  deg-0: +93101.25
  deg-2: +280.50    (6 positive, 4 negative, net slightly positive)
  deg-4: -354.75    (1 positive, 4 negative, net negative)
  Total: 93027

The Paley orientation is the UNIQUE point on the hypercube where ALL Walsh
terms align. This happens because:
1. There is exactly 1 NQR among the first m gaps {1,...,m}
   (at p=11: gap 2 is NQR, all others QR)
2. Flipping this single chord negates all Walsh coefficients with odd
   multiplicity of that chord, which EXACTLY converts all negatives to positives

This script:
1. Verifies the perfect alignment at p=7, 11
2. Computes alignment quality for ALL orientations
3. Shows why no alignment exists at p=13 (p=1 mod 4)
4. Connects to the Legendre symbol structure

Author: kind-pasteur-2026-03-12-S59c
"""

from math import comb
from collections import defaultdict


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i + s) % p] = 1
    return A


def count_ham_paths(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def all_circulant_H(p):
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]
    results = {}
    for bits in range(1 << m):
        sigma = []
        S = []
        for i, (a, b) in enumerate(pairs):
            if bits & (1 << i):
                S.append(a)
                sigma.append(+1)
            else:
                S.append(b)
                sigma.append(-1)
        S = sorted(S)
        A = build_adj(p, S)
        H = count_ham_paths(A, p)
        results[tuple(sigma)] = H
    return results


def walsh_decomposition(H_dict, m):
    n = 1 << m
    h_hat = {}
    for bits in range(n):
        S = [i for i in range(m) if bits & (1 << i)]
        coeff = 0
        for sigma_bits in range(n):
            sigma = [(1 if sigma_bits & (1 << i) else -1) for i in range(m)]
            H = H_dict[tuple(sigma)]
            prod = 1
            for i in S:
                prod *= sigma[i]
            coeff += H * prod
        h_hat[tuple(S)] = coeff / n
    return h_hat


def evaluate_walsh_at_point(h_hat, sigma, m):
    """Evaluate the Walsh expansion at a specific orientation sigma.
    Returns contributions by degree."""
    by_deg = defaultdict(float)
    terms_by_deg = defaultdict(list)

    for S, coeff in h_hat.items():
        deg = len(S)
        prod = 1
        for i in S:
            prod *= sigma[i]
        term_val = coeff * prod
        by_deg[deg] += term_val
        terms_by_deg[deg].append((S, coeff, prod, term_val))

    return by_deg, terms_by_deg


def alignment_score(terms_by_deg):
    """Compute 'alignment score': fraction of terms that are positive.
    Perfect alignment = 1.0 (all positive)."""
    total = 0
    pos = 0
    for deg, items in terms_by_deg.items():
        if deg == 0:
            continue  # skip constant
        for S, coeff, prod, val in items:
            if abs(coeff) < 1e-10:
                continue
            total += 1
            if val > 0:
                pos += 1
    return pos / total if total > 0 else 1.0


def main():
    print("=" * 70)
    print("PALEY PERFECT ALIGNMENT THEOREM")
    print("=" * 70)

    for p in [7, 11, 13]:
        m = (p - 1) // 2
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        NQR = set(range(1, p)) - QR

        print(f"\n{'='*70}")
        print(f"p={p}, m={m}, p mod 4 = {p%4}")
        print(f"{'='*70}")
        print(f"  QR = {sorted(QR)}")
        print(f"  NQR = {sorted(NQR)}")

        # Paley orientation: sigma_i = +1 iff gap (i+1) is QR
        sigma_paley = tuple(1 if (i+1) in QR else -1 for i in range(m))
        print(f"  Paley sigma = {sigma_paley}")

        # Count NQR gaps among {1,...,m}
        nqr_in_first_m = [i for i in range(m) if (i+1) in NQR]
        print(f"  NQR chords in {{0,...,{m-1}}}: {nqr_in_first_m}")
        print(f"  Number of NQR chords: {len(nqr_in_first_m)}")

        # Compute H for all orientations
        H_dict = all_circulant_H(p)
        h_hat = walsh_decomposition(H_dict, m)

        # Evaluate at Paley point
        print(f"\n  --- Paley point analysis ---")
        paley_contrib, paley_terms = evaluate_walsh_at_point(h_hat, sigma_paley, m)
        H_paley = sum(paley_contrib.values())
        print(f"  H(Paley) = {H_paley:.0f}")
        for deg in sorted(paley_contrib):
            val = paley_contrib[deg]
            n_pos = sum(1 for _,_,_,v in paley_terms[deg] if v > 1e-10 and abs(h_hat.get(_, 0)) > 1e-10)
            # Actually recount properly
            items = [(S, c, pr, v) for S, c, pr, v in paley_terms[deg] if abs(c) > 1e-10]
            n_p = sum(1 for _,_,_,v in items if v > 0)
            n_n = sum(1 for _,_,_,v in items if v < 0)
            n_z = sum(1 for _,_,_,v in items if abs(v) < 1e-10)
            print(f"    deg {deg}: {val:>+14.2f}  ({n_p}+ {n_n}- {n_z}zero)")

        aln_paley = alignment_score(paley_terms)
        print(f"  Alignment score: {aln_paley:.4f}")
        if aln_paley == 1.0:
            print(f"  *** PERFECT ALIGNMENT! ALL non-zero Walsh terms positive ***")

        # Evaluate at Interval point
        sigma_int = tuple([1] * m)
        print(f"\n  --- Interval point analysis ---")
        int_contrib, int_terms = evaluate_walsh_at_point(h_hat, sigma_int, m)
        H_int = sum(int_contrib.values())
        print(f"  H(Interval) = {H_int:.0f}")
        for deg in sorted(int_contrib):
            val = int_contrib[deg]
            items = [(S, c, pr, v) for S, c, pr, v in int_terms[deg] if abs(c) > 1e-10]
            n_p = sum(1 for _,_,_,v in items if v > 0)
            n_n = sum(1 for _,_,_,v in items if v < 0)
            print(f"    deg {deg}: {val:>+14.2f}  ({n_p}+ {n_n}-)")

        aln_int = alignment_score(int_terms)
        print(f"  Alignment score: {aln_int:.4f}")

        # Alignment score for ALL orientations
        print(f"\n  --- Alignment scores for ALL {1 << m} orientations ---")
        scores = []
        for sigma, H in sorted(H_dict.items(), key=lambda x: -x[1]):
            _, terms = evaluate_walsh_at_point(h_hat, sigma, m)
            aln = alignment_score(terms)
            scores.append((H, sigma, aln))

        # Sort by alignment
        scores_by_aln = sorted(scores, key=lambda x: -x[2])
        print(f"\n    Top 5 by alignment:")
        for H, sigma, aln in scores_by_aln[:5]:
            S = sorted([(i+1) if s == 1 else (p-i-1) for i, s in enumerate(sigma)])
            label = ""
            if sigma == sigma_paley:
                label = " [PALEY]"
            elif sigma == sigma_int:
                label = " [INTERVAL]"
            # Check if sigma equals -sigma_paley (anti-Paley = same tournament)
            neg_paley = tuple(-s for s in sigma_paley)
            if sigma == neg_paley:
                label = " [ANTI-PALEY]"
            print(f"      H={H:>12}, alignment={aln:.4f}, sigma={sigma}{label}")

        print(f"\n    Bottom 3 by alignment:")
        for H, sigma, aln in scores_by_aln[-3:]:
            S = sorted([(i+1) if s == 1 else (p-i-1) for i, s in enumerate(sigma)])
            print(f"      H={H:>12}, alignment={aln:.4f}, sigma={sigma}")

        # Correlation between alignment and H
        H_vals = [s[0] for s in scores]
        aln_vals = [s[2] for s in scores]
        mean_H = sum(H_vals) / len(H_vals)
        mean_a = sum(aln_vals) / len(aln_vals)
        cov = sum((h - mean_H) * (a - mean_a) for h, a in zip(H_vals, aln_vals)) / len(H_vals)
        std_H = (sum((h - mean_H)**2 for h in H_vals) / len(H_vals)) ** 0.5
        std_a = (sum((a - mean_a)**2 for a in aln_vals) / len(aln_vals)) ** 0.5
        corr = cov / (std_H * std_a) if std_H > 0 and std_a > 0 else 0
        print(f"\n    Correlation(H, alignment) = {corr:.4f}")

    # ====== THE PROOF SKETCH ======
    print(f"\n{'='*70}")
    print("PALEY PERFECT ALIGNMENT: PROOF SKETCH")
    print("=" * 70)
    print("""
  CLAIM: At p=3 mod 4, the Paley orientation achieves PERFECT ALIGNMENT:
    all non-zero Walsh coefficients h_hat[S] * prod_{i in S} sigma_P[i] > 0.

  MECHANISM:
  1. The Paley orientation has sigma_i = chi(i+1) where chi is the
     Legendre symbol (QR -> +1, NQR -> -1).

  2. For a Walsh coefficient h_hat[S] with |S|=d, the sign at Paley is:
     sign = sign(h_hat[S]) * prod_{i in S} chi(i+1)
          = sign(h_hat[S]) * chi(prod_{i in S}(i+1))

  3. If h_hat[S] = f(|prod gaps|) * chi(prod gaps), then
     sign = sign(f) * chi(prod)^2 = sign(f) > 0  (assuming f > 0)

  This works iff the Walsh coefficient has the form:
     h_hat[S] = |c(S)| * chi(prod_{i in S}(i+1))

  i.e., the sign of h_hat is the Legendre symbol of the product of gaps.

  VERIFICATION:
""")

    for p in [7, 11]:
        m = (p - 1) // 2
        H_dict = all_circulant_H(p)
        h_hat = walsh_decomposition(H_dict, m)
        QR = set(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)

        print(f"  p={p}:")
        # For each non-zero Walsh coefficient, check if sign = chi(prod gaps)
        all_match = True
        for S, c in sorted(h_hat.items()):
            if len(S) < 2 or abs(c) < 0.01:
                continue
            prod_gaps = 1
            for i in S:
                prod_gaps = (prod_gaps * (i + 1)) % p
            chi_prod = 1 if prod_gaps in QR else -1
            sign_c = 1 if c > 0 else -1
            match = (sign_c == chi_prod)
            if not match:
                all_match = False
            label = "MATCH" if match else "MISMATCH"
            print(f"    S={set(S)}: h={c:>+10.2f}, "
                  f"prod_gaps={prod_gaps}, chi={chi_prod:>+2d}, "
                  f"sign(h)={sign_c:>+2d} [{label}]")

        if all_match:
            print(f"    *** ALL WALSH COEFFICIENTS MATCH: sign(h) = chi(prod gaps) ***")
        else:
            print(f"    Some mismatches found")
        print()


if __name__ == '__main__':
    main()
