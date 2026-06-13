"""
Sequence speedups via CUT⊕CYCLE decomposition and H-weighting theorem.

Verified relationships:
  W(n) = Σ_{succ-free σ∈S_n} 2^{bp(σ)}
  ΣH(n) = W(n) · 2^{C(n-1,2)−n+1}        [verified n=2..7]
  |{labeled T: H=h}| = n!/h · |{tilings t: H=h}|  [H-weighting]
  E_tile[H] = W(n) / 2^{n−1}
"""

import sys
sys.path.insert(0, '/home/ubuntu/math/03-artifacts/code')
from tournament_lib import hamiltonian_path_count
from math import comb, factorial
from collections import Counter
from time import time

H = hamiltonian_path_count


def h_tiling(n, mask):
    """H for the tiling encoded by mask. Base path n-1 → n-2 → … → 0."""
    tiles = [(x, y) for y in range(n-1) for x in range(y+2, n)]
    A = [[0]*n for _ in range(n)]
    for k in range(1, n):
        A[k][k-1] = 1
    for bi, (x, y) in enumerate(tiles):
        if (mask >> bi) & 1:
            A[y][x] = 1
        else:
            A[x][y] = 1
    return H(A)


# ══════════════════════════════════════════════════════════════
# PART 1: W(n) via bitmask DP  (O(n² · 2^n))
# ══════════════════════════════════════════════════════════════
#
# W(n) = Σ_{σ∈S_n, succ-free} 2^{bp(σ)}
#   where "succ-free" means no position with σ[i+1]=σ[i]+1
#   and bp(σ) = #{i : σ[i]=σ[i+1]+1}  (base-path descents)
#
# Bitmask DP: dp[mask][j] = Σ 2^{bp} over succ-free partial perms
#   of elements in mask, ending at j.
# Transition: dp[mask][j] += dp[mask\{j}][k] * (2 if k==j+1 else 1)
#             skip if j==k+1 (succession)
#
# Total W(n) = Σ_j dp[(1<<n)-1][j].
#
print("=" * 65)
print("PART 1: W(n) via bitmask DP")
print("W(n) = Σ_{σ∈S_n succ-free} 2^{bp(σ)};  ΣH = W(n)·2^{C(n-1,2)-n+1}")
print("=" * 65)
print()

def compute_W_dp(n):
    """W(n) via bitmask DP over permutations of {0,…,n-1}."""
    dp = [[0] * n for _ in range(1 << n)]
    # Base: single-element permutations
    for j in range(n):
        dp[1 << j][j] = 1
    for mask in range(1, 1 << n):
        for j in range(n):
            if not (mask >> j & 1):
                continue
            if dp[mask][j] == 0:
                continue
            prev = mask ^ (1 << j)
            for k in range(n):
                if not (prev >> k & 1):
                    continue
                if j == k + 1:    # σ[i+1]=σ[i]+1 is a succession — skip
                    continue
                factor = 2 if k == j + 1 else 1  # k=j+1 means k↓j, base-path arc
                dp[mask][j] += dp[prev][k] * factor
    full = (1 << n) - 1
    return sum(dp[full][j] for j in range(n))

W = {}
print(f"{'n':>4} | {'W(n)':>24} | {'A000255(n)':>12} | {'2^C(n-1,2)':>16} | {'time':>8}")
print("-" * 75)

for n in range(1, 22):
    t0 = time()
    w = compute_W_dp(n)
    ms = (time() - t0) * 1000
    W[n] = w
    m = comb(n-1, 2)
    tile_space = 2**m if m <= 62 else float('inf')
    print(f"  {n:>2} | {w:>24,} | {'':>12} | {tile_space:>16,} | {ms:>7.1f}ms")
    if n >= 15 and ms > 30000:
        print(f"  (stopping at n={n})")
        break

print()
print(f"W sequence n=1..{max(W.keys())}:")
print(f"  {[W[n] for n in sorted(W.keys())]}")

# ══════════════════════════════════════════════════════════════
# PART 2: ΣH from W(n)
# ══════════════════════════════════════════════════════════════
print()
print("=" * 65)
print("PART 2: ΣH(n) = W(n) · 2^{C(n-1,2)−n+1}")
print("=" * 65)
print()
for n in sorted(W.keys()):
    f = comb(n-1, 2) - (n - 1)
    if f >= 0:
        sh = W[n] * (2**f)
        print(f"  n={n:>2}: W={W[n]:>24,}  ΣH = {sh:,}")
    else:
        print(f"  n={n:>2}: W={W[n]:>24,}  ΣH = {W[n]}·2^{{{f}}} = {W[n]/(2**(-f)):.3f}")

# ══════════════════════════════════════════════════════════════
# PART 3: H-distribution for n=3..8 via tiling enumeration
# ══════════════════════════════════════════════════════════════
print()
print("=" * 65)
print("PART 3: H-distribution — tiling enumeration speedup")
print("Formula: |{labeled T: H=h}| = n!/h · |{tilings: H=h}|")
print("=" * 65)

for n in range(3, 9):
    m = comb(n-1, 2)
    tiles = [(x, y) for y in range(n-1) for x in range(y+2, n)]
    assert len(tiles) == m

    t0 = time()
    tile_dist = Counter()
    for mask in range(1 << m):
        tile_dist[h_tiling(n, mask)] += 1
    ms = (time() - t0) * 1000

    nfact = factorial(n)
    full_dist = {}
    total_full = 0
    for h, cnt in sorted(tile_dist.items()):
        labeled = nfact * cnt // h
        assert nfact * cnt % h == 0, f"Non-integer: n={n}, h={h}"
        full_dist[h] = labeled
        total_full += labeled

    expected = 2**comb(n, 2)
    ok = "✓" if total_full == expected else f"✗ (got {total_full}, expected {expected})"
    forb = [h for h in [7, 21] if h <= max(tile_dist) and h not in tile_dist]

    print(f"\nn={n} ({ms:.0f}ms, {1<<m:,} tilings vs {2**comb(n,2):,} labeled, {ok}):")
    print(f"  Tiling dist:  {dict(sorted(tile_dist.items()))}")
    print(f"  Labeled T:    {full_dist}")
    if forb:
        print(f"  Forbidden:    {forb} ✓")
    # Sanity: ΣH from tiling dist matches W formula
    sh_check = sum(h * cnt for h, cnt in tile_dist.items())
    if n in W:
        f = comb(n-1, 2) - (n-1)
        w_sh = W[n] * (2**f) if f >= 0 else W[n] / (2**(-f))
        match = "✓" if abs(sh_check - w_sh) < 0.01 else f"✗ (W-formula gives {w_sh})"
        print(f"  ΣH check:     {sh_check:,} {match}")

# ══════════════════════════════════════════════════════════════
# PART 4: E[H] and CV² from W(n)
# ══════════════════════════════════════════════════════════════
print()
print("=" * 65)
print("PART 4: Statistical moments — E[H] and CV²[H]")
print("E_tile[H] = W(n)/2^{n−1};  E_full[H] = n!/2^{C(n,2)}")
print("W(n)/n! = 1 + CV²[H]")
print("=" * 65)
print()
print(f"{'n':>4} | {'E_tile[H]':>14} | {'E_full[H]':>14} | {'W(n)/n!':>12}")
print("-" * 50)
for n in sorted(W.keys()):
    if n > 20:
        break
    w = W[n]
    nfact = factorial(n)
    e_tile = w / (2**(n-1))
    e_full = nfact / (2**comb(n, 2)) if comb(n, 2) <= 62 else float('inf')
    ratio = w / nfact
    print(f"  {n:>2} | {e_tile:>14.6f} | {e_full:>14.8f} | {ratio:>12.8f}")

# ══════════════════════════════════════════════════════════════
# PART 5: W(n) recursion search
# ══════════════════════════════════════════════════════════════
print()
print("=" * 65)
print("PART 5: W(n) recursion search")
print("=" * 65)
wv = [W[n] for n in range(1, min(max(W.keys())+1, 20))]
ns = list(range(1, len(wv)+1))
print(f"W = {wv[:12]}")
print()

print("Checking W(n) = (an+b)·W(n-1) + (cn+d)·W(n-2):")
found_any = False
for i in range(2, min(10, len(wv))):
    n = ns[i]
    w0, w1, w2 = wv[i-2], wv[i-1], wv[i]
    for a in range(-3, 10):
        for b in range(-30, 40):
            A = a*n + b
            R = w2 - A*w1
            if w0 != 0 and R % w0 == 0:
                B = R // w0
                for c in range(-3, 6):
                    for d in range(-40, 40):
                        if c*n + d == B:
                            # Verify all remaining
                            valid = True
                            for j in range(2, min(len(wv), 15)):
                                nj = ns[j]
                                if (a*nj+b)*wv[j-1] + (c*nj+d)*wv[j-2] != wv[j]:
                                    valid = False; break
                            if valid:
                                print(f"  FOUND: W(n) = ({a}n+{b})·W(n-1) + ({c}n+{d})·W(n-2)")
                                found_any = True

if not found_any:
    print("  No 2-term linear-in-n recurrence found.")
    print()
    print("  W(n)/W(n-1) ratios:")
    for i in range(1, min(13, len(wv))):
        r = wv[i] / wv[i-1]
        print(f"    W({ns[i]})/W({ns[i-1]}) = {r:.6f}  ({ns[i]:.1f}× expected if factorial-like: {float(ns[i]):.1f})")

# ══════════════════════════════════════════════════════════════
# PART 6: OEIS summary
# ══════════════════════════════════════════════════════════════
print()
print("=" * 65)
print("PART 6: New sequences — OEIS candidates")
print("=" * 65)
print()
print("SEQUENCE A: W(n) = Σ_{σ∈S_n succ-free} 2^{bp(σ)}")
print("  = E_tile[H(n-vertex tournament)] · 2^{n-1}")
print("  = ΣH(n-vertex tournament) / 2^{C(n-1,2)−n+1}")
wlist = [W[n] for n in sorted(W.keys())]
print(f"  n=1..{max(W.keys())}: {wlist}")
print()
print("SEQUENCE B: ΣH(n) = W(n)·2^{C(n-1,2)−n+1}")
print("  = total Hamiltonian path count over all n-vertex tilings")
shlist = []
for n in sorted(W.keys()):
    f = comb(n-1, 2) - (n-1)
    if f >= 0:
        shlist.append(W[n] * (2**f))
    else:
        shlist.append(None)
print(f"  {[s for s in shlist if s is not None]}")
print()
print("SEQUENCE C: E_tile[H] = W(n)/2^{n-1} (rational)")
for n in sorted(W.keys()):
    if n > 12: break
    from fractions import Fraction
    e = Fraction(W[n], 2**(n-1))
    print(f"  n={n:>2}: W={W[n]:>12,} / 2^{n-1} = {e}")
