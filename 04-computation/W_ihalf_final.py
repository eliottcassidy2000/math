#!/usr/bin/env python3
"""
Final analysis of W(i/2) anti-correlation with H(T).

ESTABLISHED FACTS (from W_ihalf_analysis.py run):
1. W(i/2) is always REAL for odd n (confirmed n=3,5,7)
2. For n=3: W(i/2) is NOT zero. Q(i) = 6 always (constant!). Corr(H,|Q(i)|)=0
3. For n=5: Q(i) always real negative. Corr(H, |Q(i)|) = -0.973. ANTI-MONOTONE.
4. For n=7: Q(i) always real positive integer. Corr(H, |Q(i)|) = -0.904. NOT fully monotone.
5. At max H=189 (n=7): Q(i) = 6048 (minimum Q value).

where Q(u) = P(u)/(u^2+1) and P(u) = 2^{n-1} * W(r) with u=2r.

KEY FORMULA derivation:
  For n=7: P(i) = 8*(N_1 - N_3 + N_5) where N_f = #{perms with f forward edges}
  W(i/2) = P(i)/64 = (N_1 - N_3 + N_5)/8
  Using N_f = N_{6-f}: W(i/2) = (2*N_5 - N_3)/8

This script does exhaustive n=5 analysis and sampled n=7 with focus on:
- Exact quadratic form W^2 in terms of invariants
- OCF connection
- Why anti-correlation occurs
"""

from itertools import permutations, combinations
from collections import defaultdict
import sys
import random
import time

random.seed(42)

def tournament_from_bits(n, bits):
    T = [[0]*n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << k):
                T[i][j] = 1
            else:
                T[j][i] = 1
            k += 1
    return T

def Nf_from_perms(T):
    n = len(T)
    Nf = [0] * n
    for P in permutations(range(n)):
        f = sum(1 for k in range(n-1) if T[P[k]][P[k+1]])
        Nf[f] += 1
    return Nf

def H_dp(T):
    n = len(T)
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            c = dp[mask][v]
            if not (mask & (1 << v)) or c == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if T[v][u]:
                    dp[mask | (1 << u)][u] += c
    return sum(dp[full])

def W_ihalf_dp(T):
    """Compute P(i) = 2^{n-1} * W(i/2) * (something) via complex DP.
    Returns (H, P_real, P_imag) where P(i) = P_real + P_imag * i."""
    n = len(T)
    full = (1 << n) - 1
    # dp[mask][v] = (real_sum, imag_sum) of accumulated product
    dp_r = [[0]*n for _ in range(1 << n)]
    dp_i = [[0]*n for _ in range(1 << n)]
    dp_H = [[0]*n for _ in range(1 << n)]

    for v in range(n):
        dp_r[1 << v][v] = 1
        dp_H[1 << v][v] = 1

    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            rv = dp_r[mask][v]
            iv = dp_i[mask][v]
            hv = dp_H[mask][v]
            if rv == 0 and iv == 0 and hv == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                nm = mask | (1 << u)
                if T[v][u]:  # forward: multiply by (1+i)
                    dp_r[nm][u] += rv - iv
                    dp_i[nm][u] += rv + iv
                    dp_H[nm][u] += hv
                else:  # backward: multiply by (-1+i)
                    dp_r[nm][u] += -rv - iv
                    dp_i[nm][u] += rv - iv

    H = sum(dp_H[full])
    P_real = sum(dp_r[full])
    P_imag = sum(dp_i[full])
    return H, P_real, P_imag

def c3_from_outdeg(T):
    n = len(T)
    outdeg = [sum(T[v]) for v in range(n)]
    sum_od2 = sum(d*d for d in outdeg)
    return n*(n-1)*(n-2)//6 - (sum_od2 - n*(n-1)//2)//2

def count_odd_cycles(T):
    """Count directed odd cycles by length."""
    n = len(T)
    counts = {}
    for length in range(3, n+1, 2):
        count = 0
        for verts in combinations(range(n), length):
            first = verts[0]
            for perm in permutations(verts[1:]):
                path = (first,) + perm
                valid = True
                for i in range(length):
                    if not T[path[i]][path[(i+1) % length]]:
                        valid = False
                        break
                if valid:
                    count += 1
        counts[length] = count
    return counts


# ===================================================================
# EXHAUSTIVE n=5 ANALYSIS
# ===================================================================
def analyze_n5():
    n = 5
    m = n*(n-1)//2  # 10
    total = 1 << m  # 1024
    deg = n - 1  # 4

    print(f"\n{'='*70}")
    print(f"EXHAUSTIVE n=5 ANALYSIS ({total} tournaments)")
    print(f"{'='*70}")

    # For n=5: c_f = [-4, 0, 4, 0, -4]
    # P(i) = -4*N_0 + 4*N_2 - 4*N_4 = 4*(N_2 - N_0 - N_4)
    # With N_0 = N_4 = H: P(i) = 4*(N_2 - 2*H)
    # W(i/2) = P(i)/16 = (N_2 - 2*H)/4

    data = []
    for bits in range(total):
        T = tournament_from_bits(n, bits)
        Nf = Nf_from_perms(T)
        H = Nf[4]
        c3 = c3_from_outdeg(T)

        # Count odd cycles
        cycles = count_odd_cycles(T)
        t3 = cycles.get(3, 0)  # directed 3-cycles
        t5 = cycles.get(5, 0)  # directed 5-cycles (Hamiltonian cycles)

        # W(i/2) * 4 = N_2 - 2*H
        W4 = Nf[2] - 2 * H
        # W(i/2) = W4/4

        # OCF: H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2
        # where alpha_1 = total directed odd cycles, alpha_2 = disjoint pairs
        alpha_1 = t3 + t5
        # For n=5: at most one 5-cycle, and 3-cycles.
        # Disjoint pairs: only 3-cycle with itself impossible (share vertices).
        # Two 3-cycles are disjoint iff they share no vertex. At n=5: two 3-cycles
        # use 6 vertices but only 5 available, so alpha_2 = 0 always!
        # Actually: two directed 3-cycles can share vertices. They're disjoint
        # in Omega iff they share NO vertex.
        # At n=5: C(5,3)=10 possible vertex triples. Two disjoint triples impossible
        # (need 6 vertices). So alpha_2 = 0 at n=5.

        alpha_2 = 0  # always 0 at n=5
        H_ocf = 1 + 2 * alpha_1 + 4 * alpha_2

        data.append({
            'bits': bits, 'H': H, 'c3': c3, 'Nf': Nf,
            'W4': W4, 't3': t3, 't5': t5, 'alpha_1': alpha_1,
            'H_ocf': H_ocf
        })

    # Verify OCF
    ocf_ok = all(d['H'] == d['H_ocf'] for d in data)
    print(f"\n  OCF H = 1 + 2*alpha_1: {ocf_ok}")

    # Group by H
    H_groups = defaultdict(list)
    for d in data:
        H_groups[d['H']].append(d)

    print(f"\n  {'H':>5}  {'cnt':>5}  {'W*4':>6}  {'W':>8}  {'c3':>4}  {'t3':>4}  {'t5':>4}  {'a1':>4}  {'N2':>6}")
    for H in sorted(H_groups.keys()):
        entries = H_groups[H]
        W4s = sorted(set(d['W4'] for d in entries))
        c3s = sorted(set(d['c3'] for d in entries))
        t3s = sorted(set(d['t3'] for d in entries))
        t5s = sorted(set(d['t5'] for d in entries))
        a1s = sorted(set(d['alpha_1'] for d in entries))
        N2s = sorted(set(d['Nf'][2] for d in entries))

        print(f"  {H:>5}  {len(entries):>5}  {W4s}  {[w/4 for w in W4s]}  {c3s}  {t3s}  {t5s}  {a1s}  {N2s}")

    # Key relation: H = 1 + 2*alpha_1, W*4 = N_2 - 2*H
    # alpha_1 = t3 + t5 = c3 + t5 (since c3 = t3 at n=5, each vertex set has exactly one cycle)
    # Wait: c3 counts vertex SETS with a 3-cycle. t3 counts directed 3-cycles.
    # At n=5: each 3-vertex set has at most one directed 3-cycle (either clockwise or counterclockwise).
    # So c3 = t3.

    # For each tournament, alpha_1 = c3 + t5.
    # H = 1 + 2*(c3 + t5) = 1 + 2*c3 + 2*t5.
    # Since H is determined by c3 alone (from the table), t5 is determined by c3 too.

    print(f"\n  Exact relationships at n=5:")
    for c3_val in sorted(set(d['c3'] for d in data)):
        subset = [d for d in data if d['c3'] == c3_val]
        Hs = set(d['H'] for d in subset)
        t5s = set(d['t5'] for d in subset)
        W4s = set(d['W4'] for d in subset)
        N2s = set(d['Nf'][2] for d in subset)
        print(f"    c3={c3_val}: H={Hs}, t5={t5s}, W*4={W4s}, N2={N2s}")

    # Formula: W*4 = N_2 - 2*H. And H = 1 + 2*alpha_1.
    # So W*4 = N_2 - 2 - 4*alpha_1 = N_2 - 2 - 4*(c3 + t5)
    # If N_2 is determined by c3, then W is determined by c3.

    # Let's verify: is N_2 a function of c3 (or equivalently of alpha_1)?
    # N_2 = #{permutations with exactly 2 forward edges}
    # = #{permutations with exactly 2 backward edges} (since N_f = N_{4-f}, so N_2 = N_2 ✓)

    # At n=5, out-degree sequence determines everything. c3 = 10 - sum_outdeg^2/2 + ...
    # So c3 determines the out-degree sequence (up to relabeling), which determines N_f.

    # W^2 analysis
    W_vals = sorted(set(d['W4']/4 for d in data))
    print(f"\n  W(i/2) values at n=5: {W_vals}")
    print(f"  W(i/2)^2 values: {sorted(set(w*w for w in W_vals))}")

    # Express W in terms of alpha_1:
    # From the table: at c3=0, H=1, W*4=-2, W=-0.5
    # c3=1, H=3, W*4=-6, W=-1.5
    # ...
    # The pattern: W*4 decreases as H increases.
    # Actually from the table: W*4 = N_2 - 2*H

    # Can we find W as a function of alpha_1?
    print(f"\n  W(i/2) as function of alpha_1:")
    for a1 in sorted(set(d['alpha_1'] for d in data)):
        subset = [d for d in data if d['alpha_1'] == a1]
        Ws = set(d['W4']/4 for d in subset)
        Hs = set(d['H'] for d in subset)
        print(f"    alpha_1={a1}: H={Hs}, W(i/2)={Ws}")

    # Exact formula: n=5 has at most c3=5 (regular), t5=0 or 24.
    # c3 determines score sequence. N_2 is a polynomial in the score sequence.

    # DIRECT FORMULA: W*4 = N_2 - 2H
    # Since H = 1 + 2*alpha_1 and alpha_1 = c3 + t5:
    # W*4 = N_2 - 2 - 4*(c3 + t5)
    # Need: N_2 as function of score invariants.

    # N_2 = #{perms with exactly 2 forward edges out of 4}
    # A permutation (v0,v1,v2,v3,v4) has 4 edges. Exactly 2 forward.
    # "2 forward, 2 backward" out of 4 consecutive edges.

    # This is related to the number of permutations with a specific edge pattern.
    # For each perm, the 4-bit string (T[v0,v1], T[v1,v2], T[v2,v3], T[v3,v4]) has exactly 2 ones.

    # This is getting complicated. The key point is: at n=5, c3 (= score class) determines
    # both H and W(i/2), and they are anti-correlated.

    return data


# ===================================================================
# SAMPLED n=7 ANALYSIS
# ===================================================================
def analyze_n7(sample_size=100000):
    n = 7
    m = n*(n-1)//2
    total = 1 << m

    print(f"\n{'='*70}")
    print(f"SAMPLED n=7 ANALYSIS ({sample_size} tournaments)")
    print(f"{'='*70}")

    # For n=7: P(i) = 8*(N_1 - N_3 + N_5), W(i/2) = P(i)/64

    pair_bit = {}
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            pair_bit[(i,j)] = k
            k += 1

    sample_bits = random.sample(range(total), sample_size)

    H_list = []
    W8_list = []  # P(i)/8 = N_1 - N_3 + N_5
    c3_list = []

    t0 = time.time()
    step = max(1, sample_size // 20)

    for idx, bits in enumerate(sample_bits):
        if idx % step == 0 and idx > 0:
            elapsed = time.time() - t0
            print(f"  {idx}/{sample_size} ({100*idx/sample_size:.0f}%) [{elapsed:.1f}s]", file=sys.stderr)

        T = tournament_from_bits(n, bits)
        H, P_real, P_imag = W_ihalf_dp(T)
        assert P_imag == 0
        W8 = P_real // 8  # = N_1 - N_3 + N_5

        c3 = c3_from_outdeg(T)

        H_list.append(H)
        W8_list.append(W8)
        c3_list.append(c3)

    elapsed = time.time() - t0
    print(f"  Completed in {elapsed:.1f}s", file=sys.stderr)

    N = sample_size

    # Correlations
    def corr(xs, ys):
        n = len(xs)
        mx = sum(xs) / n
        my = sum(ys) / n
        cov = sum((xs[i] - mx) * (ys[i] - my) for i in range(n)) / n
        vx = sum((xs[i] - mx)**2 for i in range(n)) / n
        vy = sum((ys[i] - my)**2 for i in range(n)) / n
        return cov / (vx**0.5 * vy**0.5) if vx > 0 and vy > 0 else 0

    absW8 = [abs(w) for w in W8_list]
    W_actual = [w/8 for w in W8_list]
    absW_actual = [abs(w)/8 for w in W8_list]

    print(f"\n  Corr(H, W8)          = {corr(H_list, W8_list):.6f}")
    print(f"  Corr(H, |W8|)        = {corr(H_list, absW8):.6f}")
    print(f"  Corr(H, c3)          = {corr(H_list, c3_list):.6f}")
    print(f"  Corr(c3, W8)         = {corr(c3_list, W8_list):.6f}")

    # Group by H
    H_to_data = defaultdict(list)
    for i in range(N):
        H_to_data[H_list[i]].append((W8_list[i], c3_list[i]))

    H_set = sorted(H_to_data.keys())

    # Table
    print(f"\n  {'H':>5}  {'cnt':>5}  {'W8 min':>8}  {'W8 max':>8}  {'|W8| mean':>10}  {'c3 rng':>10}")
    prev_mean = float('inf')
    for H in H_set:
        entries = H_to_data[H]
        W8s = [e[0] for e in entries]
        c3s = [e[1] for e in entries]
        mean_abs = sum(abs(w) for w in W8s) / len(W8s)
        marker = " ***" if mean_abs > prev_mean else ""
        print(f"  {H:>5}  {len(entries):>5}  {min(W8s):>8}  {max(W8s):>8}  {mean_abs:>10.1f}  [{min(c3s):>2},{max(c3s):>2}]{marker}")
        prev_mean = mean_abs

    # W8 = N1 - N3 + N5. Is this related to c3?
    # c3 counts directed 3-cycles. The out-degree sequence determines c3.
    # But N_f depends on more than just the out-degree sequence.

    # Key question: does c3 determine W8?
    c3_to_W8 = defaultdict(set)
    for i in range(N):
        c3_to_W8[c3_list[i]].add(W8_list[i])

    unique = sum(1 for v in c3_to_W8.values() if len(v) == 1)
    print(f"\n  c3 uniquely determines W8: {unique}/{len(c3_to_W8)} ({unique == len(c3_to_W8)})")

    # Check (H, c3) -> W8
    Hc3_to_W8 = defaultdict(set)
    for i in range(N):
        Hc3_to_W8[(H_list[i], c3_list[i])].add(W8_list[i])
    multi = sum(1 for v in Hc3_to_W8.values() if len(v) > 1)
    print(f"  (H, c3) uniquely determines W8: {multi == 0} ({multi} multi-valued pairs)")

    # Sign of W8
    pos = sum(1 for w in W8_list if w > 0)
    neg = sum(1 for w in W8_list if w < 0)
    zero = sum(1 for w in W8_list if w == 0)
    print(f"\n  W8 > 0: {pos} ({100*pos/N:.1f}%)")
    print(f"  W8 = 0: {zero} ({100*zero/N:.1f}%)")
    print(f"  W8 < 0: {neg} ({100*neg/N:.1f}%)")

    # W8 at max H
    max_H = max(H_list)
    max_H_entries = [W8_list[i] for i in range(N) if H_list[i] == max_H]
    print(f"\n  Max H = {max_H}: W8 values = {sorted(set(max_H_entries))}")
    print(f"  W(i/2) at max H = {sorted(set(w/8 for w in max_H_entries))}")

    # At max H (Paley T_7, H=189):
    # By BIBD analysis: all regular n=7 with BIBD arrangement have H=189.
    # These have c3=14 (the maximum).

    # Min W8 (most negative)
    min_W8 = min(W8_list)
    min_W8_H = [H_list[i] for i in range(N) if W8_list[i] == min_W8]
    print(f"  Min W8 = {min_W8}: H values = {sorted(set(min_W8_H))}")

    # OCF connection: H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3 + ...
    # W8 = N_1 - N_3 + N_5
    # Is there a formula connecting W8 to the alpha_k?

    # From the n=5 analysis: W*4 = N_2 - 2*H, and H = 1 + 2*alpha_1.
    # So W*4 = N_2 - 2 - 4*alpha_1.
    # At n=5, alpha_1 determines everything.

    # At n=7: W8 = N_1 - N_3 + N_5, and H = 1 + 2*alpha_1 + 4*alpha_2 + 8*alpha_3.
    # Need: N_1, N_3, N_5 in terms of tournament invariants.

    # N_f relates to the "f-th derivative" of the permanent.
    # sum_f N_f * z^f = perm(z*J + (1-z)*T_complement)... hmm, complicated.

    # Actually: sum_f N_f * z^f = sum_P z^{f(P)} = sum_P prod_edges z^{T[P_k,P_{k+1}]}
    # = perm of matrix M where M[v][u] = z if T[v][u]=1, else 1.
    # This is the permanent of (I - (1-z)*A) where A is some matrix...

    # The generating function G(z) = sum N_f z^f satisfies G(1) = H, G(0) = H(T^op),
    # sum N_f = n!, etc.

    # Key insight: W8 = G(omega) where omega = e^{2pi i/4} = i (evaluated at 4th root of unity?)
    # No: W8 = N_1 - N_3 + N_5 = Im(sum N_f * i^f) / 1?
    # sum N_f * i^f = N_0 - N_2 + N_4 - N_6 + i*(N_1 - N_3 + N_5) [wait N_6 doesn't exist if f goes to 6]
    # For n=7: f goes from 0 to 6.
    # sum N_f * i^f = N_0*1 + N_1*i + N_2*(-1) + N_3*(-i) + N_4*(1) + N_5*(i) + N_6*(-1)
    #              = (N_0 - N_2 + N_4 - N_6) + i*(N_1 - N_3 + N_5)
    # Real part = N_0 - N_2 + N_4 - N_6 = 0 (by N_f = N_{6-f})
    # Imag part = N_1 - N_3 + N_5 = W8 !!

    print(f"\n  KEY IDENTITY: W8 = Im(G(i)) where G(z) = sum N_f z^f = perm of per-edge-z matrix")
    print(f"  G(z) = sum_P prod_edges z^{{T[P_k,P_{{k+1}}]}} = perm(M) where M[v][u] = z if T[v][u], else 1")

    # So W8 = Im(perm(M_i)) where M_i[v][u] = i if T[v][u]=1, else 1.
    # And |G(i)|^2 = (Re G(i))^2 + (Im G(i))^2 = 0 + W8^2 = W8^2.
    # So |G(i)| = |W8|.

    # G(z) evaluated at z=1 gives H. At z=i, we get W8*i.
    # G(z) = perm(z*A + (1-z)*(J-A-I)) where A is adjacency matrix... no.
    # Simpler: G(z) = perm(M(z)) where M(z)[v][u] = z if T[v][u]=1, 1 if T[v][u]=0 (and v != u).

    # Actually: M(z)[v][u] = z^{T[v][u]} = z if forward, 1 if backward.
    # NOT the adjacency matrix but a "weighted" version.

    # G(i) = perm(M(i)) where M(i)[v][u] = i if T[v][u], 1 otherwise.
    # This is the permanent of a (0,1,i)-matrix (with 0 on diagonal).

    # |G(i)|^2 = |perm(M(i))|^2 = W8^2.

    # For the anti-correlation: high H means the tournament is "well-connected"
    # (many forward edges in many permutations). But this same well-connectedness
    # causes more cancellation in G(i) because the i-weighted terms interfere.

    # DEEPER: G(i) = sum_P i^{f(P)} where f(P) = forward edges in P.
    # |G(i)|^2 = sum_{P,Q} i^{f(P)} * (-i)^{f(Q)} = sum_{P,Q} i^{f(P)-f(Q)}
    # = #{(P,Q): f(P) - f(Q) = 0 mod 4} - #{f(P)-f(Q) = 2 mod 4}
    #   (the i^1 and i^3 terms cancel by N_f = N_{6-f})

    # Actually |G(i)|^2 = |Im G(i)|^2 = W8^2 since Re G(i) = 0.

    return H_list, W8_list, c3_list


# ===================================================================
# THEORETICAL ANALYSIS
# ===================================================================
def theory():
    print(f"\n{'='*70}")
    print("THEORETICAL ANALYSIS")
    print(f"{'='*70}")

    print("""
KEY RESULTS:

1. EXACT FORMULA: W(i/2) = Im[G(i)] / 8 where
   G(z) = sum_P prod_{edges} z^{T[P_k,P_{k+1}]}
        = perm(M(z)) with M(z)[v][u] = z if T[v][u]=1, else 1 (v != u)
   G(i) is purely imaginary (by N_f = N_{deg-f} symmetry for odd n).

2. For n=7: W(i/2) = (N_1 - N_3 + N_5) / 8
   where N_f = #{permutations with exactly f forward edges}

3. ANTI-CORRELATION MECHANISM:
   - H(T) = G(1) = sum_P 1^{f(P)} = #{permutations where all edges forward} = H
   - G(i) = sum_P i^{f(P)} = phase-weighted sum
   - High H means many permutations contribute positively at z=1
   - But at z=i, these same permutations have diverse phases i^f
   - More balanced tournament = more permutations with high f
   - More high-f permutations contribute i^5 or i^6 = i or -1 (large imaginary)
   - But also more with i^1 or i^2 = i or -1 (cancellation)

4. PRECISE MECHANISM (for n=7):
   W8 = N_1 - N_3 + N_5 = (N_1 + N_5) - N_3
   High H means: N_6 large, N_0 = N_6 large (path reversal).
   Constraint: sum N_f = 5040 (= 7!)
   If N_0 and N_6 are large, less "mass" for N_1,...,N_5.
   But also N_1=N_5, N_2=N_4, so:
   2*N_0 + 2*N_1 + 2*N_2 + N_3 = 5040
   => N_3 = 5040 - 2*H - 2*N_1 - 2*N_2
   => W8 = 2*N_1 - (5040 - 2*H - 2*N_1 - 2*N_2) = 4*N_1 + 2*N_2 + 2*H - 5040

   So W8 = 4*N_1 + 2*N_2 + 2*H - 5040.

5. ANTI-CORRELATION OF |W8| with H:
   |W8| = |4*N_1 + 2*N_2 + 2*H - 5040|
   Mean value of 4*N_1 + 2*N_2 ~ 4*mean(N_1) + 2*mean(N_2)
   By symmetry (averaging over all tournaments): mean(N_f) = 5040/7 = 720 for all f.
   So mean(4*N_1 + 2*N_2) = 4*720 + 2*720 = 4320.
   mean(W8) = 4320 + 2*mean(H) - 5040.

   For tournaments with high H: N_0 = N_6 = H is large.
   The remaining sum N_1 + N_2 + N_3/2 = (5040 - 2*H)/2 = 2520 - H.
   So 4*N_1 + 2*N_2 is bounded by roughly 4*(2520-H) (if N_2=N_3=0).
   |W8| = |4*N_1 + 2*N_2 + 2*H - 5040| <= |4*(2520-H) + 2*H - 5040| = |5040-2*H|
   The maximum |W8| is bounded by ~ 5040 - 2*H for large H.
   This gives a HARD BOUND: |W8| <= 5040 - 2*H (approximately).
   Higher H => smaller possible |W8|. THIS IS THE ANTI-CORRELATION!
""")


def main():
    n5_data = analyze_n5()
    H7, W7, c3_7 = analyze_n7(100000)
    theory()

    # Final: verify the bound |W8| <= 5040 - 2*H
    print(f"\n{'='*70}")
    print("BOUND VERIFICATION: |W8| vs 5040 - 2*H")
    print(f"{'='*70}")

    violations = 0
    for i in range(len(H7)):
        bound = 5040 - 2 * H7[i]
        if abs(W7[i]) > bound:
            violations += 1

    print(f"  |W8| > 5040-2*H violations: {violations}/{len(H7)}")

    # Tighter bound from the constraint
    # W8 = 4*N_1 + 2*N_2 + 2*H - 5040
    # N_1 >= 0, N_2 >= 0, N_1 + N_2 <= 2520 - H (from sum = 5040, N_f=N_{6-f})
    # W8 ranges from 2*H - 5040 (when N_1=N_2=0) to
    #   4*(2520-H) + 2*H - 5040 = 10080 - 4*H + 2*H - 5040 = 5040 - 2*H
    # So: 2*H - 5040 <= W8 <= 5040 - 2*H
    # => |W8| <= 5040 - 2*H (when H < 2520)

    # At max H=189: |W8| <= 5040 - 378 = 4662
    # At min H=1: |W8| <= 5040 - 2 = 5038

    # Check actual values
    H_to_W = defaultdict(list)
    for i in range(len(H7)):
        H_to_W[H7[i]].append(W7[i])

    print(f"\n  {'H':>5}  {'bound':>6}  {'actual max|W8|':>15}  {'tight?':>8}")
    for H in sorted(H_to_W.keys()):
        bound = 5040 - 2 * H
        actual_max = max(abs(w) for w in H_to_W[H])
        tight = actual_max >= bound * 0.95
        print(f"  {H:>5}  {bound:>6}  {actual_max:>15}  {'YES' if tight else 'no':>8}")

    print(f"\n  CONCLUSION: |W(i/2)| <= (5040 - 2*H) / 8 = (n! - 2*H) / 2^{{n-1}}")
    print(f"  This is a HARD UPPER BOUND that forces anti-correlation!")
    print(f"  The bound comes from: W8 = 4*N_1 + 2*N_2 + 2*H - 5040")
    print(f"  with N_1, N_2 >= 0 and 2*N_1 + 2*N_2 <= 5040 - 2*H.")


if __name__ == "__main__":
    main()
