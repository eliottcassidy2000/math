"""
schur_concavity_test.py — Test if H is Schur-concave on the spectral simplex

Schur-concavity: f(x) ≥ f(y) whenever y majorizes x (y is "more spread").
The center of the simplex is the MINIMUM in majorization order.
If H is Schur-concave → center maximizes H.

For p ≡ 3 mod 4: center IS the Paley tournament (spectrally flat).
Schur-concavity would PROVE Paley maximizes H among circulants.

For p ≡ 1 mod 4: center is not achievable as a tournament.
The maximum might occur at a boundary point.

Test at p=7, 11 (exhaustive) and check all pairs, not just center pairs.

Author: opus-2026-03-12-S60
"""
import sys
import time
import numpy as np
from collections import defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)


def hamiltonian_paths_dp(A, n):
    dp = defaultdict(lambda: defaultdict(int))
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        if not dp[mask]:
            continue
        for v in dp[mask]:
            if dp[mask][v] == 0:
                continue
            for w in range(n):
                if mask & (1 << w):
                    continue
                if A[v][w]:
                    dp[mask | (1 << w)][w] += dp[mask][v]
    return sum(dp[full][v] for v in range(n))


def circulant_adj(n, S):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for s in S:
            A[i][(i+s)%n] = 1
    return A


def circulant_eigenvalues(n, S):
    omega = np.exp(2j * np.pi / n)
    return [sum(omega ** (k * s) for s in S) for k in range(n)]


def circulant_tournaments(n):
    k = (n - 1) // 2
    pairs = [(d, n - d) for d in range(1, k + 1)]
    for bits in range(1 << k):
        S = set()
        for i in range(k):
            if bits & (1 << i):
                S.add(pairs[i][0])
            else:
                S.add(pairs[i][1])
        yield frozenset(S)


def is_paley(n, S):
    qr = set(pow(x, 2, n) for x in range(1, n))
    return set(S) == qr or set(S) == set(range(1, n)) - qr


def majorizes(x, y, tol=1e-6):
    """Check if x majorizes y: partial sums of sorted(x,↓) ≥ partial sums of sorted(y,↓).
    And total sum equal."""
    sx = sorted(x, reverse=True)
    sy = sorted(y, reverse=True)
    if abs(sum(sx) - sum(sy)) > tol:
        return False
    for k in range(len(sx)):
        if sum(sx[:k+1]) < sum(sy[:k+1]) - tol:
            return False
    return True


def main():
    print("SCHUR-CONCAVITY OF H ON THE SPECTRAL SIMPLEX")
    print("=" * 75)

    for p in [7, 11, 13]:
        print(f"\n{'='*75}")
        print(f"p = {p} (p mod 4 = {p % 4})")
        print(f"{'='*75}")

        m = (p - 1) // 2
        tournaments = list(circulant_tournaments(p))

        # Group by spectral class
        spectral_classes = {}
        for S in tournaments:
            evals = circulant_eigenvalues(p, S)
            y2_half = tuple(sorted(round(evals[k].imag**2, 8) for k in range(1, m + 1)))
            if y2_half not in spectral_classes:
                A = circulant_adj(p, S)
                H = hamiltonian_paths_dp(A, p)
                spectral_classes[y2_half] = {'y2': list(y2_half), 'H': H, 'S': S,
                                              'paley': is_paley(p, S)}

        classes = list(spectral_classes.values())
        nc = len(classes)
        print(f"  {nc} spectral classes")

        # Check all pairs for majorization and Schur-concavity
        comparable = 0
        concave_ok = 0
        concave_fail = 0
        failures = []

        for i in range(nc):
            for j in range(i+1, nc):
                xi = classes[i]['y2']
                xj = classes[j]['y2']
                Hi = classes[i]['H']
                Hj = classes[j]['H']

                if majorizes(xi, xj):
                    # xi ≻ xj (xi more spread)
                    comparable += 1
                    # Schur-concave: H(xi) ≤ H(xj) (more spread → less H)
                    if Hi <= Hj:
                        concave_ok += 1
                    else:
                        concave_fail += 1
                        failures.append((i, j, 'i≻j but H(i)>H(j)'))
                elif majorizes(xj, xi):
                    # xj ≻ xi (xj more spread)
                    comparable += 1
                    # Schur-concave: H(xj) ≤ H(xi)
                    if Hj <= Hi:
                        concave_ok += 1
                    else:
                        concave_fail += 1
                        failures.append((j, i, 'j≻i but H(j)>H(i)'))

        total_pairs = nc * (nc - 1) // 2
        incomparable = total_pairs - comparable
        print(f"  Total pairs: {total_pairs}")
        print(f"  Comparable (in majorization order): {comparable}")
        print(f"  Incomparable: {incomparable}")
        print(f"  Schur-concave OK: {concave_ok}/{comparable}")
        print(f"  Schur-concave FAIL: {concave_fail}/{comparable}")

        if failures:
            print(f"\n  FAILURES:")
            for i, j, msg in failures:
                ci = classes[i]
                cj = classes[j]
                print(f"    {msg}")
                print(f"      i: y²={ci['y2']}, H={ci['H']}")
                print(f"      j: y²={cj['y2']}, H={cj['H']}")

        # Summary
        if concave_fail == 0 and comparable > 0:
            print(f"\n  ✓ H is Schur-concave on all comparable pairs!")
            if any(c['paley'] for c in classes):
                print(f"  ✓ This PROVES Paley maximizes H among all circulants on Z_{p}")
        elif concave_fail > 0:
            print(f"\n  ✗ H is NOT Schur-concave on the spectral simplex")

        # Also check: is H a Schur-concave function of |λ|² values?
        # (alternative: use |λ|² instead of y²)
        print(f"\n  Classes sorted by H:")
        for c in sorted(classes, key=lambda c: -c['H']):
            label = " ★ PALEY" if c['paley'] else ""
            sigma2 = sum(y**2 for y in c['y2'])
            print(f"    H={c['H']:>10}, σ₂={sigma2:>10.2f}, y²={[round(y,3) for y in c['y2']]}{label}")


if __name__ == '__main__':
    main()
    print("\nDONE.")
