#!/usr/bin/env python3
"""
hhat_magnitude_formula.py -- Investigate the formula for |h_hat[{a,b}]|
as a function of resonance level q.

KNOWN DATA:
  p=7:  |h_hat_2| = 3.5 = p/2 (all q=3)
  p=11: |h_hat_2| = 272.25 (q=3), 8.25 (q=5)
         Ratio = 33 = 3p
         272.25 = 9*p^2/4, 8.25 = 3*p/4

CONJECTURE: |h_hat[{a,b}]| depends only on q, and there's a formula
in terms of p and q.

From the resonance cascade: the D^{2n} contributions are nonzero starting
at n=(q+1)/2, and the sign is chi(q)*chi(ab). The magnitude of D^{2n}
grows geometrically with n.

Let's check at p=7 and p=11 whether h_hat values follow a simple formula.

Author: kind-pasteur-2026-03-12-S60
"""


def legendre(a, p):
    if a % p == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def classify_resonance(a, b, p):
    resonances = []
    for k in range(1, p):
        q = 2*k - 1
        if q >= p:
            break
        if (q*a - b) % p == 0:
            resonances.append((q, f"{q}a=b"))
        if (q*a + b) % p == 0:
            resonances.append((q, f"{q}a=-b"))
        if (a - q*b) % p == 0 and q != 1:
            resonances.append((q, f"a={q}b"))
        if (a + q*b) % p == 0 and q != 1:
            resonances.append((q, f"a=-{q}b"))
    return resonances


def count_ham_paths_fast(adj_list, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            cnt = dp[mask][v]
            if cnt == 0:
                continue
            for w in adj_list[v]:
                if not (mask & (1 << w)):
                    dp[mask | (1 << w)][w] += cnt
    full = (1 << n) - 1
    return sum(dp[full][v] for v in range(n))


def compute_all_H(p):
    m = (p - 1) // 2
    pairs = [(s, p - s) for s in range(1, m + 1)]
    H_vals = {}
    for bits in range(1 << m):
        sigma = [(1 if bits & (1 << i) else -1) for i in range(m)]
        S = sorted(pairs[i][0] if sigma[i] == 1 else pairs[i][1] for i in range(m))
        adj_list = [[] for _ in range(p)]
        for v in range(p):
            for s in S:
                adj_list[v].append((v + s) % p)
        H_vals[bits] = count_ham_paths_fast(adj_list, p)
    return H_vals, m


def compute_hhat_deg2(H_vals, m, a_idx, b_idx):
    """Compute h_hat[{a_idx, b_idx}] with standard normalization."""
    total = 0
    for bits in range(1 << m):
        sigma_a = 1 if bits & (1 << a_idx) else -1
        sigma_b = 1 if bits & (1 << b_idx) else -1
        total += H_vals[bits] * sigma_a * sigma_b
    return total / (1 << m)


def main():
    print("=" * 70)
    print("|h_hat| MAGNITUDE FORMULA INVESTIGATION")
    print("=" * 70)

    for p in [7, 11]:
        m = (p - 1) // 2
        print(f"\np={p}, m={m}")
        print(f"  H(Paley) = ?")

        H_vals, _ = compute_all_H(p)

        # Find Paley orientation
        QR = sorted(j for j in range(1, p) if pow(j, (p-1)//2, p) == 1)
        paley_bits = 0
        pairs = [(s, p-s) for s in range(1, m+1)]
        for i in range(m):
            if pairs[i][0] in QR:
                paley_bits |= (1 << i)
        print(f"  H(Paley) = {H_vals[paley_bits]}")
        print(f"  Mean H = {sum(H_vals.values())/(1<<m):.4f}")

        # Compute all degree-2 Walsh coefficients
        print(f"\n  Degree-2 h_hat:")
        from collections import defaultdict
        q_hhat = defaultdict(list)
        for a in range(m):
            for b in range(a+1, m):
                h = compute_hhat_deg2(H_vals, m, a, b)
                gap_a, gap_b = a+1, b+1
                res = classify_resonance(gap_a, gap_b, p)
                min_q = min(qq for qq, t in res) if res else 'inf'
                q_hhat[min_q].append(abs(h))
                chi_prod = legendre(gap_a * gap_b, p)
                print(f"    ({a},{b}): gaps=({gap_a},{gap_b}), q={min_q}, "
                      f"h_hat={h:>12.4f}, |h_hat|={abs(h):>10.4f}, "
                      f"chi(prod)={chi_prod:+d}")

        print(f"\n  |h_hat| by q-level:")
        for q in sorted(q_hhat.keys()):
            vals = q_hhat[q]
            val = vals[0]  # all should be same
            chi_q = legendre(q, p)
            print(f"    q={q}: |h_hat| = {val:.4f}")
            print(f"           = {val} = {val*4:.0f}/4")
            print(f"           val*4 = {val*4:.0f}")
            print(f"           val*4 / p = {val*4/p:.6f}")
            print(f"           val*4 / p^2 = {val*4/p**2:.6f}")

        # Check specific formulas
        print(f"\n  Formula check:")
        # D^{2n} magnitude at onset: |D^{2n}| = p^{n-1/2} * ...
        # For circulant: D^{2n} involves Gauss sums
        # Q_k(Paley) = (p+1)/4 for all k
        # The degree-2 Walsh of the CYCLE COUNT c_k (not H) is simpler
        # because c_k decomposes as a character sum

        # H(sigma) can be expanded: H = sum 2^j * alpha_j(Omega(sigma))
        # Each alpha_j is a product of indicator functions for independent sets
        # in the conflict graph Omega(sigma)

        # For a SIMPLER approach: consider the degree-2 Walsh of
        # c_3 = number of directed 3-cycles = (1/3) sum_k |S_hat(k)|^6 / |S_hat(k)|^2
        # No, that's wrong.

        # Actually, c_k = (1/k) sum_t (sum_{s in S} omega^{st})^k
        # For circulant: c_k = (1/k) sum_{t=0}^{p-1} S_hat(t)^k

        # Degree-2 Walsh of c_3:
        # c_3(sigma) = (p/3) + (1/3) sum_{t=1}^{p-1} S_hat(t;sigma)^3
        # S_hat(t;sigma) = sum_{i=0}^{m-1} sigma_i * (omega^{t(i+1)} - omega^{-t(i+1)})
        #                = -1/2 + i * D(t;sigma)  where D is the sine sum

        # This gets complicated. Let me just check empirical formulas.

    # PART 2: Ratio analysis
    print("\n\n" + "=" * 60)
    print("RATIO ANALYSIS")
    print("=" * 60)

    # p=7: only q=3, |h_hat|=3.5 = p/2
    # p=11: q=3: |h_hat|=272.25 = (3p)^2/4 = 9p^2/4
    #        q=5: |h_hat|=8.25 = 3p/4
    #        ratio = 3p

    # Let me check: does |h_hat(q)| = (D^{q+1} value) * correction?
    # D^4 at onset is proportional to ((p+1)/4)^2 = (p+1)^2/16
    # D^6 at onset is proportional to ((p+1)/4)^3 = (p+1)^3/64
    # Ratio D^4/D^6 = 64/16 * (4/(p+1)) = 4*4/(p+1) = 16/(p+1)
    # At p=11: 16/12 = 4/3. But actual ratio is 33. So D^4/D^6 doesn't explain it.

    print("\nFormula candidates for |h_hat_2(q)|:")
    for p in [7, 11]:
        m = (p-1)//2
        H_vals, _ = compute_all_H(p)
        for q in [3, 5]:
            if q >= p:
                continue
            # Find representative pair
            for a in range(1, m+1):
                for b in range(a+1, m+1):
                    res = classify_resonance(a, b, p)
                    if res and min(qq for qq, t in res) == q:
                        h = abs(compute_hhat_deg2(H_vals, m, a-1, b-1))
                        n_onset = (q+1)//2
                        print(f"  p={p}, q={q}: |h_hat|={h:.4f}")
                        print(f"    p^{n_onset-1} = {p**(n_onset-1)}")
                        print(f"    |h_hat| / p^{n_onset-1} = {h/p**(n_onset-1):.6f}")
                        print(f"    |h_hat| * 4 / p^{n_onset} = {h*4/p**n_onset:.6f}")
                        print(f"    H(Paley) = {H_vals[paley_bits] if p==11 else '?'}")
                        break
                else:
                    continue
                break

    # PART 3: Connection to H(Paley) - mean(H)
    print("\n\nPART 3: CONNECTION TO H VARIANCE")
    for p in [7, 11]:
        m = (p-1)//2
        H_vals, _ = compute_all_H(p)

        mean_H = sum(H_vals.values()) / len(H_vals)
        var_H = sum((H - mean_H)**2 for H in H_vals.values()) / len(H_vals)

        # Walsh: Var(H) = sum |h_hat[S]|^2 (Parseval)
        deg2_energy = 0
        for a in range(m):
            for b in range(a+1, m):
                h = compute_hhat_deg2(H_vals, m, a, b)
                deg2_energy += h**2

        print(f"\n  p={p}:")
        print(f"    mean(H) = {mean_H:.4f}")
        print(f"    Var(H)  = {var_H:.4f}")
        print(f"    deg-2 energy = {deg2_energy:.4f}")
        print(f"    deg-2 fraction of Var = {deg2_energy/var_H:.4f}")

        # Per-q breakdown
        from collections import defaultdict
        q_energy = defaultdict(float)
        q_count = defaultdict(int)
        for a in range(m):
            for b in range(a+1, m):
                h = compute_hhat_deg2(H_vals, m, a, b)
                gap_a, gap_b = a+1, b+1
                res = classify_resonance(gap_a, gap_b, p)
                min_q = min(qq for qq, t in res) if res else 'inf'
                q_energy[min_q] += h**2
                q_count[min_q] += 1

        for q in sorted(q_energy.keys()):
            print(f"    q={q}: energy={q_energy[q]:.4f}, "
                  f"count={q_count[q]}, "
                  f"per-pair={q_energy[q]/q_count[q]:.4f}, "
                  f"fraction={q_energy[q]/var_H:.4f}")


if __name__ == '__main__':
    main()
