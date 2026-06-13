# ⚠️ WARNING: This script uses QR mod p for p ≡ 1 (mod 4), which does NOT
# produce a tournament (S ∩ (-S) ≠ ∅ gives bidirectional edges).
# Results for those primes are INVALID. See MISTAKE-011b.
# Valid Paley tournaments require p ≡ 3 (mod 4).

#!/usr/bin/env python3
"""
PROOF: Palindromic N(a,b,j) for circulant tournaments => M = (H/n)*I at odd n.

THEOREM: For any circulant tournament T on Z/nZ at odd n, M = (H/n)*I.

PROOF:

Step 1: Define f(d,j) = N(0, d, j) for d in {1,...,n-1}, j in {0,...,n-2}.
By translation symmetry of T: N(a,b,j) = f(b-a mod n, j).

Step 2: N is symmetric: N(a,b,j) = N(b,a,j).
So f(d,j) = N(0,d,j) = N(d,0,j) = N(0,-d,j) = f(n-d, j). [Translation by -d]

Step 3: Self-complementarity. The map sigma: i -> -i (mod n) sends T to T^op
(since T has generator S and T^op has generator {n-d : d in S} = -S mod n,
and sigma*S = -S).

Path reversal: Ham(T) <-> Ham(T^op) via (v_0,...,v_{n-1}) -> (v_{n-1},...,v_0).
Composition with sigma: phi(v_0,...,v_{n-1}) = (-v_{n-1},...,-v_0) is a bijection Ham(T) -> Ham(T).

For pair {a,b} at positions {j,j+1}: phi maps to pair {-a,-b} at {n-2-j, n-1-j}.
So N(a,b,j) = N(-a,-b, n-2-j).
In terms of f: f(d,j) = f(-d mod n, n-2-j) = f(n-d, n-2-j).

Step 4: Combine Steps 2 and 3:
f(d,j) = f(n-d,j)  ... (Step 2)
f(d,j) = f(n-d, n-2-j)  ... (Step 3)
=> f(n-d,j) = f(n-d, n-2-j)  [substitute (2) into LHS of (3)]
=> f(e,j) = f(e, n-2-j) for all e  [rename e = n-d]

This is palindromicity of f in j.

Step 5: At odd n, palindromic f(e,j) = f(e,n-2-j) with n-1 even implies:
sum_j (-1)^j f(e,j) = 0.
(Standard: substituting k = n-2-j gives sum_k (-1)^{n-2-k} f(e,k) = (-1)^n sum_k (-1)^k f(e,k).
At odd n: (-1)^n = -1, so alt_sum = -alt_sum => alt_sum = 0.)

Step 6: M[a,b] = sum_j (-1)^j N(a,b,j) = sum_j (-1)^j f(b-a,j) = 0 for all a != b.
Combined with M[a,a] = H/n (from position-uniformity, which all circulants have),
we get M = (H/n)*I. QED.

VERIFICATION at n=5, n=7, n=9, n=11, and n=13 (Paley).

kind-pasteur-2026-03-06-S25d
"""

import sys

def circulant_tournament(n, gen_set):
    T = {}
    for i in range(n):
        for j in range(n):
            if i == j: continue
            T[(i,j)] = 1 if (j - i) % n in gen_set else 0
    return T

def compute_f_dp(T, n):
    """Compute f(d,j) = N(0,d,j) using DP Hamiltonian paths.
    Returns f as dict, and H."""
    full = (1 << n) - 1

    # dp[mask][last] = # Ham paths visiting mask ending at last
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1

    for mask in range(1, 1 << n):
        for last in range(n):
            if not ((mask >> last) & 1):
                continue
            cnt = dp.get((mask, last), 0)
            if cnt == 0:
                continue
            for nxt in range(n):
                if (mask >> nxt) & 1:
                    continue
                if T.get((last, nxt), 0) == 0:
                    continue
                nkey = (mask | (1 << nxt), nxt)
                dp[nkey] = dp.get(nkey, 0) + cnt

    # dp_suf[mask][first] = # Ham paths visiting mask starting at first
    dp_suf = {}
    for v in range(n):
        dp_suf[(1 << v, v)] = 1
    for popcount in range(2, n + 1):
        for mask in range(1, 1 << n):
            if bin(mask).count('1') != popcount:
                continue
            for first in range(n):
                if not ((mask >> first) & 1):
                    continue
                total = 0
                for nxt in range(n):
                    if nxt == first or not ((mask >> nxt) & 1):
                        continue
                    if T.get((first, nxt), 0) == 0:
                        continue
                    total += dp_suf.get((mask & ~(1 << first), nxt), 0)
                if total > 0:
                    dp_suf[(mask, first)] = total

    H = sum(dp.get((full, v), 0) for v in range(n))

    # Compute f(d,j) = N(0,d,j) for d=1..n-1, j=0..n-2
    f = [[0]*(n-1) for _ in range(n)]  # f[d][j]

    for d in range(1, n):
        # Count paths where {0, d} appears at positions {j, j+1}
        # Two cases: 0->d at (j,j+1) if T[0,d]=1, and d->0 at (j,j+1) if T[d,0]=1
        for x, y in [(0, d), (d, 0)]:
            if T.get((x, y), 0) == 0:
                continue
            for j in range(n-1):
                pc = j + 1
                for prefix_mask in range(1, 1 << n):
                    if bin(prefix_mask).count('1') != pc:
                        continue
                    if not ((prefix_mask >> x) & 1):
                        continue
                    if (prefix_mask >> y) & 1:
                        continue
                    pcnt = dp.get((prefix_mask, x), 0)
                    if pcnt == 0:
                        continue
                    scnt = dp_suf.get((full ^ prefix_mask, y), 0)
                    if scnt == 0:
                        continue
                    f[d][j] += pcnt * scnt

    return f, H


def verify_theorem(n, gen_set, label):
    """Verify the palindromic N theorem for a circulant tournament."""
    T = circulant_tournament(n, gen_set)
    f, H = compute_f_dp(T, n)

    print(f"\n  {label}: n={n}, S={sorted(gen_set)}, H={H}")

    # Check 1: f(d,j) = f(n-d,j) (symmetry)
    sym_ok = True
    for d in range(1, n):
        for j in range(n-1):
            if f[d][j] != f[n-d][j]:
                sym_ok = False
                break
    print(f"    f(d,j) = f(n-d,j) [symmetry]: {sym_ok}")

    # Check 2: f(d,j) = f(d,n-2-j) (palindromic)
    pal_ok = True
    for d in range(1, n):
        for j in range(n-1):
            if f[d][j] != f[d][n-2-j]:
                pal_ok = False
                break
    print(f"    f(d,j) = f(d,n-2-j) [palindromic]: {pal_ok}")

    # Check 3: alternating sum = 0
    alt_ok = True
    for d in range(1, n):
        alt = sum((-1)**j * f[d][j] for j in range(n-1))
        if alt != 0:
            alt_ok = False
            print(f"    f({d}): alt_sum = {alt} != 0!")
    print(f"    sum (-1)^j f(d,j) = 0 for all d: {alt_ok}")

    # Show some f values
    for d in range(1, min(n, 4)):
        print(f"    f({d}) = {f[d]}")

    if pal_ok and alt_ok:
        print(f"    => M = ({H}/{n})*I = {H//n}*I")

    return pal_ok and alt_ok


# ============================================================
print("=" * 70)
print("PALINDROMIC N THEOREM: Verification for circulant tournaments")
print("=" * 70)

all_pass = True

# n=5 Paley: QR mod 5 = {1,4}
all_pass &= verify_theorem(5, {1, 4}, "Paley T_5")

# n=7 Paley: QR mod 7 = {1,2,4}
all_pass &= verify_theorem(7, {1, 2, 4}, "Paley T_7")

# n=7 non-Paley circulant
all_pass &= verify_theorem(7, {1, 2, 3}, "Circulant {1,2,3}")

# n=9 circulant
all_pass &= verify_theorem(9, {1, 2, 4, 6}, "Circulant {1,2,4,6}")

# n=11 Paley: QR mod 11 = {1,3,4,5,9}
all_pass &= verify_theorem(11, {1, 3, 4, 5, 9}, "Paley T_11")

print(f"\n{'='*70}")
if all_pass:
    print("ALL VERIFIED. Theorem holds for all tested circulant tournaments.")
else:
    print("SOME FAILURES!")

# ============================================================
# n=13 Paley: QR mod 13 = {1,3,4,9,10,12}
# ============================================================
print(f"\n{'='*70}")
print("n=13 Paley: Testing (may take a while)...")
print("=" * 70)
all_pass &= verify_theorem(13, {1, 3, 4, 9, 10, 12}, "Paley T_13")

print(f"\n{'='*70}")
print(f"FINAL RESULT: {'ALL PASS' if all_pass else 'FAILURES'}")
print("=" * 70)
