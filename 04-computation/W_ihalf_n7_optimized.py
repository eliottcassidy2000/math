#!/usr/bin/env python3
"""
Optimized W(i/2) analysis for n=7 tournaments.

Strategy: Use optimized bitmask DP with arrays instead of lists,
and only track what we need.

For n=7, W(i/2) = (2*N_5 - N_3)/8 where N_f = #{permutations with f forward edges}.

Actually, we only need H = N_6 and W*8 = 2*N_5 - N_3.
The full Nf distribution is expensive. Can we compute H and W*8 more directly?

H is computed by standard bitmask DP: O(2^n * n^2) = 128*49 ~ 6300 per tournament.
For W*8 we need N_3 and N_5. We need the full DP with f-tracking.

Alternative: use the polynomial evaluation approach.
2^6 * W(r) = P(2r) where P(u) = sum_P prod_edges ((u+1) or (u-1))
W(i/2) = P(i) / 64

P(i) = sum_P prod_edges ((i+1) or (i-1))
(i+1)^a * (i-1)^b where a+b=6, a = forward edges, b = backward edges.

Since i+1 = 1+i and i-1 = -(1-i):
(1+i)^a * (-(1-i))^b = (-1)^b * (1+i)^a * (1-i)^b

For a+b=6, (-1)^b = (-1)^{6-a} = (-1)^a (since (-1)^6 = 1).
= (-1)^a * (1+i)^a * (1-i)^{6-a}

(1+i)(1-i) = 2, so:
If a <= 3: (-1)^a * 2^a * (1-i)^{6-2a}
If a > 3:  (-1)^a * 2^{6-a} * (1+i)^{2a-6}

Powers of (1-i): (1-i)^0=1, (1-i)^2=-2i, (1-i)^4=-4, (1-i)^6=8i
Powers of (1+i): (1+i)^0=1, (1+i)^2=2i, (1+i)^4=-4, (1+i)^6=-8i

a=0: (-1)^0 * 1 * (1-i)^6 = 8i
a=1: (-1)^1 * 2 * (1-i)^4 = -2*(-4) = 8
a=2: (-1)^2 * 4 * (1-i)^2 = 4*(-2i) = -8i
a=3: (-1)^3 * 8 * (1-i)^0 = -8
a=4: (-1)^4 * 4 * (1+i)^2 = 4*(2i) = 8i
a=5: (-1)^5 * 2 * (1+i)^4 = -2*(-4) = 8
a=6: (-1)^6 * 1 * (1+i)^6 = -8i

Coefficient table (real, imag):
a=0: (0, 8)
a=1: (8, 0)
a=2: (0, -8)
a=3: (-8, 0)
a=4: (0, 8)
a=5: (8, 0)
a=6: (0, -8)

P(i) = sum_a N_a * c_a
Real part = 8*(N_1 - N_3 + N_5)  [matches our earlier derivation!]
Imag part = 8*(N_0 - N_2 + N_4 - N_6) = 0 by N_f = N_{6-f}

So P(i) is purely real = 8*(N_1 - N_3 + N_5)
W(i/2) = P(i)/64 = (N_1 - N_3 + N_5)/8

With N_1 = N_5: W(i/2) = (2*N_5 - N_3)/8. Confirmed.

OPTIMIZATION IDEA: We don't need full N_f. We need a WEIGHTED sum.
Define S = N_1 - N_3 + N_5 (= N_5 - N_3 + N_1).
Since N_f = N_{6-f}: S = N_5 - N_3 + N_5 = 2*N_5 - N_3.

We can compute S directly using a modified DP where instead of tracking
full f, we track (-1)^f weighted by whether f is odd or even, etc.

Actually, S = sum_P (-1)^{(a(P) mod 4 < 2 ? 0 : 1)} * (a(P) mod 2 == 1 ? 1 : 0)
Hmm, that's not simpler.

Better approach: use the POLYNOMIAL evaluation directly in the DP.
Instead of tracking N_f, track the contribution to P(i).

For each partial path ending at v through mask with current complex value z,
extending to u adds a factor of (i+1) if forward, (i-1) if backward.

z_new = z * (i+1) if T[v][u]=1
z_new = z * (i-1) if T[v][u]=0

Since P(i) is always real, we can track the real and imaginary parts.
At the end, P(i) = sum of real parts (imag cancels to 0).

DP state: dp[mask][v] = (real_sum, imag_sum)
Transitions:
  Forward (T[v][u]=1): multiply by (1+i): new_r = r - im, new_im = r + im
  Backward: multiply by (-1+i): new_r = -r - im, new_im = r - im

This is O(2^n * n^2) per tournament with small constant. Very fast!
"""

import sys
import time
from collections import defaultdict
import math

def main():
    n = 7
    m = n * (n - 1) // 2  # 21
    total = 1 << m  # 2,097,152

    print(f"W(i/2) vs H(T) for all {total} tournaments on n={n}")
    print(f"W(i/2) = P(i)/64 where P(i) = Re[sum_P prod((i+1) or (i-1))]")
    print("=" * 70)

    # Precompute: for each tournament (bits), compute H and P(i) using DP
    # T[a][b] is determined by bit position of pair (a,b) with a<b.

    # Bit ordering: k-th bit corresponds to pair (i,j) with i<j in lex order
    # k=0: (0,1), k=1: (0,2), k=2: (1,2), k=3: (0,3), ...
    # For i<j: k = i*n - i*(i+1)/2 + (j - i - 1) ... actually let's just precompute
    pair_to_bit = {}
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            pair_to_bit[(i, j)] = k
            k += 1

    # For the DP, we need T[v][u] for any v,u with v!=u.
    # T[v][u] = 1 if v<u and bit set, or v>u and bit NOT set.
    # T[i][j] = (bits >> pair_to_bit[(i,j)]) & 1  if i < j
    # T[j][i] = 1 - T[i][j]

    H_all = []
    Pi_all = []  # P(i) values (integer, always real)

    # Also compute c3 for each tournament (for correlation analysis)
    # c3 = C(n,3) - sum_{v} C(s_v, 2) where s_v is out-degree... wait that's for
    # undirected. For tournaments: c3 = C(n,3) - sum_v C(d+(v), 2) ... no.
    # Actually: number of 3-cycles = C(n,3) - sum_{v} C(d+(v), 2)
    # where d+(v) = out-degree of v. This is a known formula.
    # sum C(d+,2) = sum d+(d+-1)/2 = (sum d+^2 - sum d+)/2
    # sum d+ = C(n,2) = 21. sum d+^2 depends on tournament.
    # c3 = 35 - (sum d+^2 - 21)/2 = 35 - sum_d+^2/2 + 21/2 = 35 + 10.5 - sum_d+^2/2
    # = 45.5 - sum_d+^2/2
    # Hmm, c3 must be integer, so sum_d+^2 must be odd... let me recheck.
    # c3 = C(n,3) - sum_v C(outdeg(v), 2)
    # = 35 - sum_v outdeg(v)*(outdeg(v)-1)/2
    # sum outdeg(v) = 21, so sum outdeg^2 = sum (outdeg^2).
    # c3 = 35 - (sum_outdeg^2 - 21)/2

    c3_all = []

    t0 = time.time()
    step = total // 20

    for bits in range(total):
        if bits % step == 0 and bits > 0:
            elapsed = time.time() - t0
            rate = bits / elapsed
            eta = (total - bits) / rate
            print(f"  {bits}/{total} ({100*bits/total:.0f}%) [{elapsed:.0f}s, ETA {eta:.0f}s]", file=sys.stderr)

        # === Compute H(T) using standard DP ===
        # dp_H[mask][v] = number of paths through mask ending at v
        dp_H = [[0] * n for _ in range(1 << n)]
        for v in range(n):
            dp_H[1 << v][v] = 1

        full = (1 << n) - 1
        for mask in range(1, 1 << n):
            for v in range(n):
                cnt = dp_H[mask][v]
                if not (mask & (1 << v)) or cnt == 0:
                    continue
                for u in range(n):
                    if mask & (1 << u):
                        continue
                    # Check T[v][u]
                    if v < u:
                        fwd = (bits >> pair_to_bit[(v, u)]) & 1
                    else:
                        fwd = 1 - ((bits >> pair_to_bit[(u, v)]) & 1)
                    if fwd:
                        dp_H[mask | (1 << u)][u] += cnt

        H = sum(dp_H[full][v] for v in range(n))

        # === Compute P(i) using complex DP ===
        # dp[mask][v] = (real, imag) accumulator
        # Initialize: dp[{v}][v] = (1, 0)
        # Transition: multiply by (1+i) if forward, (-1+i) if backward

        dp_r = [[0] * n for _ in range(1 << n)]
        dp_i = [[0] * n for _ in range(1 << n)]
        for v in range(n):
            dp_r[1 << v][v] = 1

        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)):
                    continue
                r_val = dp_r[mask][v]
                i_val = dp_i[mask][v]
                if r_val == 0 and i_val == 0:
                    continue
                for u in range(n):
                    if mask & (1 << u):
                        continue
                    if v < u:
                        fwd = (bits >> pair_to_bit[(v, u)]) & 1
                    else:
                        fwd = 1 - ((bits >> pair_to_bit[(u, v)]) & 1)

                    new_mask = mask | (1 << u)
                    if fwd:  # multiply by (1+i)
                        dp_r[new_mask][u] += r_val - i_val
                        dp_i[new_mask][u] += r_val + i_val
                    else:  # multiply by (-1+i)
                        dp_r[new_mask][u] += -r_val - i_val
                        dp_i[new_mask][u] += r_val - i_val

        Pi_real = sum(dp_r[full][v] for v in range(n))
        Pi_imag = sum(dp_i[full][v] for v in range(n))

        # Pi_imag should be 0
        assert Pi_imag == 0, f"Pi_imag = {Pi_imag} at bits={bits}"

        # === Compute c3 via out-degree formula ===
        outdeg = [0] * n
        for v in range(n):
            for u in range(n):
                if v == u:
                    continue
                if v < u:
                    outdeg[v] += (bits >> pair_to_bit[(v, u)]) & 1
                else:
                    outdeg[v] += 1 - ((bits >> pair_to_bit[(u, v)]) & 1)
        sum_od2 = sum(d * d for d in outdeg)
        c3 = 35 - (sum_od2 - 21) // 2

        H_all.append(H)
        Pi_all.append(Pi_real)
        c3_all.append(c3)

    elapsed = time.time() - t0
    print(f"\nCompleted {total} tournaments in {elapsed:.1f}s", file=sys.stderr)

    N = total
    # W(i/2) = Pi / 64
    # W(i/2) * 8 = Pi / 8

    # =====================================================================
    print(f"\n{'='*70}")
    print("BASIC CHECKS")
    print(f"{'='*70}")

    # Pi divisibility by 8
    div8 = sum(1 for p in Pi_all if p % 8 == 0)
    print(f"  P(i) always divisible by 8: {div8 == N} ({div8}/{N})")

    # W = Pi/64
    W_all = [p // 8 for p in Pi_all]  # W * 8

    print(f"  W*8 range: [{min(W_all)}, {max(W_all)}]")
    print(f"  W*8 always even: {all(w % 2 == 0 for w in W_all)}")

    # =====================================================================
    print(f"\n{'='*70}")
    print("CORRELATIONS")
    print(f"{'='*70}")

    def corr(xs, ys):
        n = len(xs)
        mx = sum(xs) / n
        my = sum(ys) / n
        cov = sum((xs[i] - mx) * (ys[i] - my) for i in range(n)) / n
        vx = sum((xs[i] - mx)**2 for i in range(n)) / n
        vy = sum((ys[i] - my)**2 for i in range(n)) / n
        return cov / (vx**0.5 * vy**0.5) if vx > 0 and vy > 0 else 0

    absW = [abs(w) for w in W_all]

    print(f"  Corr(H, W*8)         = {corr(H_all, W_all):.6f}")
    print(f"  Corr(H, |W*8|)       = {corr(H_all, absW):.6f}")
    print(f"  Corr(H, c3)          = {corr(H_all, c3_all):.6f}")
    print(f"  Corr(c3, W*8)        = {corr(c3_all, W_all):.6f}")
    print(f"  Corr(c3, |W*8|)      = {corr(c3_all, absW):.6f}")

    # =====================================================================
    print(f"\n{'='*70}")
    print("JOINT DISTRIBUTION: H vs W(i/2)")
    print(f"{'='*70}")

    H_to_data = defaultdict(list)
    for i in range(N):
        H_to_data[H_all[i]].append((W_all[i], c3_all[i]))

    H_set = sorted(H_to_data.keys())

    print(f"\n  {'H':>5}  {'count':>7}  {'W*8 min':>8}  {'W*8 max':>8}  {'|W|*8 min':>9}  {'|W|*8 max':>9}  {'W/8 min':>8}  {'W/8 max':>8}  {'c3':>6}")
    for H in H_set:
        entries = H_to_data[H]
        Ws = [e[0] for e in entries]
        c3s = set(e[1] for e in entries)
        mn_W, mx_W = min(Ws), max(Ws)
        mn_abs = min(abs(w) for w in Ws)
        mx_abs = max(abs(w) for w in Ws)
        c3_str = str(sorted(c3s)) if len(c3s) <= 3 else f"{min(c3s)}-{max(c3s)}"

        print(f"  {H:>5}  {len(entries):>7}  {mn_W:>8}  {mx_W:>8}  {mn_abs:>9}  {mx_abs:>9}  {mn_W/8:>8.2f}  {mx_W/8:>8.2f}  {c3_str:>6}")

    # =====================================================================
    print(f"\n{'='*70}")
    print("MONOTONICITY OF MEAN |W| vs H")
    print(f"{'='*70}")

    prev_mean = float('inf')
    violations = 0
    for H in H_set:
        entries = H_to_data[H]
        absWs = [abs(e[0]) for e in entries]
        mean_w = sum(absWs) / len(absWs)
        marker = " ***" if mean_w > prev_mean else ""
        if mean_w > prev_mean:
            violations += 1
        print(f"  H={H:>3}: mean|W|*8={mean_w:>10.1f}, median={sorted(absWs)[len(absWs)//2]:>8}, n={len(entries):>6}{marker}")
        prev_mean = mean_w

    print(f"\n  Mean |W| monotonicity violations: {violations}/{len(H_set)-1}")

    # =====================================================================
    print(f"\n{'='*70}")
    print("W(i/2) = 0 ANALYSIS")
    print(f"{'='*70}")

    zero_idx = [i for i in range(N) if W_all[i] == 0]
    print(f"  W(i/2) = 0 count: {len(zero_idx)}")

    zero_H = defaultdict(int)
    for i in zero_idx:
        zero_H[H_all[i]] += 1
    print(f"  H distribution when W=0: {dict(sorted(zero_H.items()))}")

    # Max H
    max_H = max(H_all)
    max_H_idx = [i for i in range(N) if H_all[i] == max_H]
    print(f"\n  H-maximizers (H={max_H}): {len(max_H_idx)} tournaments")
    max_W = set(W_all[i] for i in max_H_idx)
    print(f"  W*8 values at max H: {sorted(max_W)}")
    print(f"  W(i/2) at max H: {sorted(w/8 for w in max_W)}")

    # =====================================================================
    print(f"\n{'='*70}")
    print("|W(i/2)|^2 ANALYSIS")
    print(f"{'='*70}")

    # W(i/2) = W*8 / 8
    # W^2 = (W*8)^2 / 64
    W2_scaled = sorted(set(w * w for w in W_all))
    print(f"  Distinct (W*8)^2 values: {len(W2_scaled)}")

    # Are they all divisible by 64?
    all_div64 = all(w2 % 64 == 0 for w2 in W2_scaled)
    print(f"  (W*8)^2 always divisible by 64: {all_div64}")

    if all_div64:
        W2_int = sorted(set(w2 // 64 for w2 in W2_scaled))
        print(f"  W^2 integer values: {W2_int}")
        # Perfect squares?
        for w2 in W2_int[:10]:
            sq = math.isqrt(w2)
            ps = sq*sq == w2
            print(f"    W^2={w2}: sqrt={sq if ps else math.sqrt(w2):.4f}, perfect_square={ps}")
    else:
        # Check what divides
        for d in [2, 4, 8, 16, 32]:
            cnt = sum(1 for w2 in W2_scaled if w2 % d == 0)
            print(f"  (W*8)^2 div by {d}: {cnt}/{len(W2_scaled)}")

    # =====================================================================
    print(f"\n{'='*70}")
    print("KEY FINDING SUMMARY")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
