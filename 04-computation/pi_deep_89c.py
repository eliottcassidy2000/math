#!/usr/bin/env python3
"""
pi_deep_89c.py — Deep dives into the most promising π connections

opus-2026-03-14-S89c

Following up on Part 1 findings:
1. Paley H-ratio H(P_p)/E[H] appears to converge — to what?
2. Eigenvalue phases arg(λ_k)/π — is there a pattern?
3. The doubly-regular property of Paley (exact common outneighbors)
4. E[H]·E[1/H] sequence: 1.25, 1.6, 1.9, 2.0... — converges to what?
5. IP(P_n, 2) = (2^{n+2} - (-1)^n)/3 — a Jacobsthal-like sequence
6. The "tournament Basel problem": is Σ_T 1/H(T) expressible?
"""

import math
import random
import cmath
from collections import Counter

def count_H(adj, n):
    dp = [0] * ((1 << n) * n)
    for v in range(n):
        dp[(1 << v) * n + v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            val = dp[mask * n + v]
            if val == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[(mask | (1 << u)) * n + u] += val
    full = (1 << n) - 1
    return sum(dp[full * n + v] for v in range(n))

def is_qr(a, p):
    if a % p == 0:
        return False
    return pow(a, (p-1)//2, p) == 1

def paley_tournament(p):
    adj = [[0]*p for _ in range(p)]
    for i in range(p):
        for j in range(p):
            if i != j and is_qr(j - i, p):
                adj[i][j] = 1
    return adj

def random_tournament(n):
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return adj

# ═══════════════════════════════════════════════════════════════
print("=" * 70)
print("DEEP DIVE 1: Paley H-ratio convergence")
print("=" * 70)

# H(P_p) / E[H] where E[H] = p!/2^{p-1}
# Known values:
paley_H = {3: 3, 7: 189, 11: 95095, 19: 1172695746915}

print("\nH(P_p) / (p!/2^{p-1}):")
for p, H in sorted(paley_H.items()):
    mean_H = math.factorial(p) / 2**(p-1)
    ratio = H / mean_H
    log_ratio = math.log(ratio)
    print(f"  p={p:2d}: H(P_p)={H}, E[H]={mean_H:.2f}, ratio={ratio:.6f}, log(ratio)={log_ratio:.4f}")

# Also compute H for the cyclic tournament (i→j if j-i in {1,...,⌊n/2⌋})
print("\nH(C_n) / (n!/2^{n-1}) for cyclic tournaments:")
for n in [3, 5, 7, 9, 11, 13]:
    adj = [[0]*n for _ in range(n)]
    forward = set(range(1, n//2 + 1))
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in forward:
                adj[i][j] = 1
    H = count_H(adj, n)
    mean_H = math.factorial(n) / 2**(n-1)
    ratio = H / mean_H
    print(f"  n={n:2d}: H(C_n)={H}, ratio={ratio:.6f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("DEEP DIVE 2: Eigenvalue Phases and π")
print("=" * 70)

# For Paley P_p, the non-trivial eigenvalues of the adjacency matrix are
# λ_k = Σ_{j∈QR} ω^{jk} for k=1,...,p-1
# where ω = e^{2πi/p}
# These are Gauss-type sums. For p ≡ 3 mod 4:
# λ_k = (-1 + i·√p)/2 if k is QR, (-1 - i·√p)/2 if k is QNR

print("\nPhase structure arg(λ_k)/π:")
for p in [7, 11, 19, 23, 43]:
    if p % 4 != 3:
        continue
    omega = cmath.exp(2j * cmath.pi / p)
    qr_set = {a for a in range(1, p) if is_qr(a, p)}

    phases_qr = []
    phases_nqr = []
    for k in range(1, p):
        lam = sum(omega**(j*k) for j in qr_set)
        phase = cmath.phase(lam) / math.pi
        if k in qr_set:
            phases_qr.append(phase)
        else:
            phases_nqr.append(phase)

    # For p ≡ 3 mod 4: QR eigenvalue phase = arctan(√p)/π
    # QNR eigenvalue phase = -arctan(√p)/π (= π - arctan(√p))/π
    expected_qr = math.atan(math.sqrt(p)) / math.pi
    # Actually: λ = (-1 ± i√p)/2
    # phase = atan2(±√p, -1)/π = ±(π - atan(√p))/π
    expected_phase = (math.pi - math.atan(math.sqrt(p))) / math.pi

    print(f"\n  P_{p}: QR phases: [{phases_qr[0]:.6f}, ...], QNR phases: [{phases_nqr[0]:.6f}, ...]")
    print(f"    Expected QR phase: ±{expected_phase:.6f}")
    print(f"    All QR same?  {all(abs(ph - phases_qr[0]) < 1e-10 or abs(ph + phases_qr[0]) < 1e-10 for ph in phases_qr)}")
    print(f"    All QNR same? {all(abs(ph - phases_nqr[0]) < 1e-10 or abs(ph + phases_nqr[0]) < 1e-10 for ph in phases_nqr)}")
    print(f"    QR phase + QNR phase = {abs(phases_qr[0]) + abs(phases_nqr[0]):.6f} (should be 1.0)")

# Key observation: as p → ∞, the phase → arctan(√p)/π → 1/2
# So eigenvalue phases approach ±π/2 (purely imaginary!)
print("\nPhase convergence as p → ∞:")
print("  Phase = (π - arctan(√p))/π → 1/2 as p → ∞")
for p in [3, 7, 11, 19, 43, 67, 103, 163, 283, 523, 1087]:
    if p % 4 != 3:
        continue
    phase = (math.pi - math.atan(math.sqrt(p))) / math.pi
    print(f"  p={p:5d}: phase = {phase:.8f}, 0.5 + 1/(π√p) ≈ {0.5 + 1/(math.pi*math.sqrt(p)):.8f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("DEEP DIVE 3: Jacobsthal Numbers — IP(path, 2)")
print("=" * 70)

# IP(P_n, 2) = (2^{n+2} - (-1)^n) / 3
# This is the Jacobsthal number J(n+2)!
# Jacobsthal: J(0)=0, J(1)=1, J(n) = J(n-1) + 2J(n-2)
# J(n) = (2^n - (-1)^n) / 3

# So IP(path_n, 2) = J(n+2) = (2^{n+2} - (-1)^{n+2})/3 = (2^{n+2} - (-1)^n)/3 ✓

jacobsthal = [0, 1]
for i in range(2, 20):
    jacobsthal.append(jacobsthal[-1] + 2*jacobsthal[-2])
print("\nJacobsthal numbers:")
print(f"  J(n) = {jacobsthal[:15]}")
print(f"  IP(path_n, 2) = J(n+2) = {jacobsthal[2:14]}")
print(f"  Formula: J(n) = (2^n - (-1)^n)/3")

# Jacobsthal identity: J(n)² + J(n+1)² = J(2n+1)·... nah
# More interesting: J(n) mod p
print("\nJ(n) mod 7:")
print("  " + " ".join(str(j % 7) for j in jacobsthal[:20]))
# Period?
for per in range(1, 20):
    if all(jacobsthal[i] % 7 == jacobsthal[i+per] % 7 for i in range(min(8, len(jacobsthal)-per))):
        print(f"  Period mod 7 = {per}")
        break

# Connection to H: for the PATH GRAPH, IP(P_n, 2) = J(n+2)
# For general tournament conflict graphs, IP(G, 2) = H(T) which is much more complex
# But the Jacobsthal sequence is the "baseline" — the simplest possible conflict graph

print("\nJacobsthal ratio J(n+1)/J(n) → 2:")
for i in range(2, 15):
    print(f"  J({i+1})/J({i}) = {jacobsthal[i+1]/jacobsthal[i]:.6f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("DEEP DIVE 4: The Tournament Sum 1/H")
print("=" * 70)

# Σ_T 1/H(T) over all tournaments on n vertices
# This is a "harmonic tournament sum"
# We know Σ H = n!·2^{m-n+1}
# What about Σ 1/H?

import itertools

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(pairs)
    for bits in range(1 << m):
        adj = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(pairs):
            if bits & (1 << idx):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj

from fractions import Fraction

print("\nΣ_T 1/H(T) as exact fraction:")
for n in range(3, 7):
    total = Fraction(0)
    for adj in all_tournaments(n):
        h = count_H(adj, n)
        total += Fraction(1, h)

    N = 2**(n*(n-1)//2)
    mean_inv = total / N
    print(f"\n  n={n}: Σ 1/H = {total} = {float(total):.6f}")
    print(f"    N = {N}, E[1/H] = {total}/{N} = {float(mean_inv):.8f}")
    print(f"    Simplified: {mean_inv}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("DEEP DIVE 5: H Distribution — Exact Moments")
print("=" * 70)

# Compute exact moments of H
for n in range(3, 7):
    h_values = []
    for adj in all_tournaments(n):
        h_values.append(count_H(adj, n))

    N = len(h_values)
    m1 = sum(h_values) / N
    m2 = sum(h**2 for h in h_values) / N
    m3 = sum(h**3 for h in h_values) / N
    m4 = sum(h**4 for h in h_values) / N

    # As exact fractions
    m1_f = Fraction(sum(h_values), N)
    m2_f = Fraction(sum(h**2 for h in h_values), N)

    var = m2 - m1**2
    # Var(H) — can we express this?
    # E[H²] = E[(Σ_π 1_π)²] = Σ_{π,σ} P(π and σ both HP)
    # = Σ_{π,σ} 2^{m - |constraints(π∪σ)|}

    print(f"\n  n={n}: E[H]={m1_f}, E[H²]={m2_f}")
    print(f"    Var(H) = {float(m2_f - m1_f**2):.4f}")
    print(f"    E[H²]/E[H]² = {float(m2_f)/(float(m1_f)**2):.6f}")
    print(f"    E[H³]/E[H]³ = {m3/m1**3:.6f}")
    print(f"    E[H⁴]/E[H]⁴ = {m4/m1**4:.6f}")

    # Var(H)/E[H]² is the coefficient of variation squared
    cv2 = var / m1**2
    print(f"    CV² = Var/E² = {cv2:.6f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("DEEP DIVE 6: E[H²] — Pairs of Hamiltonian Paths")
print("=" * 70)

# E[H²] = Σ_{π,σ} P(π,σ both HP of T)
# For each pair of permutations (π,σ), the pair is HP iff
# all n-1 arcs of π and all n-1 arcs of σ go the right way.
# If π and σ share k arcs, they constrain 2(n-1)-k arcs total.
# P(both HP) = 2^{m - (2(n-1)-k)} = 2^{m-2(n-1)+k}
#
# So E[H²] = Σ_{π,σ} 2^{k(π,σ) - 2(n-1) + m} / 2^m
#           = 2^{-2(n-1)} Σ_{π,σ} 2^{k(π,σ)}
#
# where k(π,σ) = number of shared arcs between permutations π and σ

print("\nE[H²] via pair counting:")
for n in range(3, 8):
    m = n*(n-1)//2
    perms = list(itertools.permutations(range(n)))

    # For each pair of perms, count shared directed edges
    total = 0
    for pi in perms:
        for sigma in perms:
            # Arcs of π: π(0)→π(1), π(1)→π(2), ..., π(n-2)→π(n-1)
            pi_arcs = set((pi[i], pi[i+1]) for i in range(n-1))
            sigma_arcs = set((sigma[i], sigma[i+1]) for i in range(n-1))
            k = len(pi_arcs & sigma_arcs)
            total += 2**k

    E_H2 = total / 2**(2*(n-1))
    E_H = math.factorial(n) / 2**(n-1)
    ratio = E_H2 / E_H**2

    print(f"  n={n}: E[H²] = {E_H2:.2f}, E[H]² = {E_H**2:.2f}, E[H²]/E[H]² = {ratio:.6f}")

    if n <= 6:
        # Verify against direct computation
        h_values = [count_H(adj, n) for adj in all_tournaments(n)]
        direct_E_H2 = sum(h**2 for h in h_values) / len(h_values)
        print(f"    Direct verification: E[H²] = {direct_E_H2:.2f} ✓" if abs(E_H2 - direct_E_H2) < 0.01 else f"    MISMATCH: {direct_E_H2:.2f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("DEEP DIVE 7: Shared Arc Distribution")
print("=" * 70)

# The shared arc count k(π,σ) has a distribution.
# k can range from 0 to n-1.
# The distribution of k is related to the "overlap" of two random permutations.

for n in range(3, 8):
    perms = list(itertools.permutations(range(n)))
    k_counts = Counter()
    for pi in perms:
        for sigma in perms:
            pi_arcs = set((pi[i], pi[i+1]) for i in range(n-1))
            sigma_arcs = set((sigma[i], sigma[i+1]) for i in range(n-1))
            k = len(pi_arcs & sigma_arcs)
            k_counts[k] += 1

    total_pairs = len(perms)**2
    print(f"\n  n={n}: shared arc distribution k(π,σ):")
    E_k = 0
    for k in sorted(k_counts.keys()):
        prob = k_counts[k] / total_pairs
        E_k += k * prob
        print(f"    k={k}: P={prob:.6f} ({k_counts[k]}/{total_pairs})")
    print(f"    E[k] = {E_k:.6f}")
    # Expected shared arcs: each arc appears in (n-1)!/n! = 1/(n) fraction of perms
    # Wait: each DIRECTED arc (u,v) appears in exactly (n-1)! · 1/2... no.
    # Actually, the number of permutations where u immediately precedes v
    # is (n-1)! / 1 — there are (n-1)! permutations with u at any position
    # ... let me think. For a random perm π, P(π(i)=u, π(i+1)=v) = 1/(n(n-1))
    # P(directed arc (u,v) in π) = (n-1)·1/(n(n-1)) = 1/n
    # P(two independent arcs both in π and σ) ≈ 1/n² (but not independent)
    # E[k(π,σ)] = Σ_{arcs} P(arc in π)·P(arc in σ) = (n-1) · (1/n)² = (n-1)/n²
    # Wait, that doesn't match. Let me count differently.
    # Number of directed arcs in a permutation = n-1
    # Each directed arc appears in (n-2)! permutations (fix the pair, permute rest)
    # No: permutations of {0,...,n-1} where π(i)=u, π(i+1)=v for SOME i
    # = (n-1) positions × (n-2)! / ... actually just (n-1)!
    # Hmm. P(specific directed arc in π) = (n-1)! / n! · (n-1 choose for positions)
    # Actually: there are n-1 consecutive pairs in a permutation. Each ordered pair
    # (u,v) with u≠v appears in exactly (n-2)! × (n-1) / ??? No.
    # Just: P(u→v in π) = #{perms with u immediately before v} / n!
    # = (n-1)! / n! = 1/n (treating the pair uv as one unit, (n-1)! arrangements)
    # So E[k] = (n-1) · (1/n)² ... no, that's E[k] for independent perms.
    # E[k(π,σ)] = Σ_{(u,v)} P((u,v) in π AND (u,v) in σ) = Σ P(in π)² = n(n-1)·(1/n)² = (n-1)/n
    expected_E_k = (n - 1) / n  # Hmm wait, there are n(n-1) directed arcs
    # E[k] = Σ_{directed (u,v)} P(in π)·P(in σ) = n(n-1)·(1/n)·(1/n) = (n-1)/n ... wait
    # P((u,v) in π) = 1/n for each of the n(n-1) directed pairs
    # But each perm has exactly n-1 arcs
    # E[k] = Σ_{(u,v)} P(in π)·P(in σ) where π,σ independent
    # = n(n-1) · (1/n)² = (n-1)/n
    # Hmm, n(n-1)/n² = (n-1)/n. For n=3: 2/3 = 0.667, actual E_k = ?
    print(f"    Predicted E[k] = (n-1)/n = {(n-1)/n:.6f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("DEEP DIVE 8: E[H²] Closed Form")
print("=" * 70)

# E[H²] = 2^{-2(n-1)} · Σ_{π,σ} 2^{k(π,σ)}
# = 2^{-2(n-1)} · Σ_k #{pairs with k shared arcs} · 2^k
# = 2^{-2(n-1)} · n!² · E_pairs[2^k]

# Let's see if E[H²] has a nice form
print("\nE[H²] values:")
for n in range(3, 8):
    perms = list(itertools.permutations(range(n)))
    total = sum(
        2**sum(1 for i in range(n-1) if pi[i]==sigma[i] and pi[i+1]==sigma[i+1])
        for pi in perms for sigma in perms
    )
    # Wait that's wrong — need to count shared DIRECTED ARCS, not positional matches
    # Let me redo
    total = 0
    for pi in perms:
        pi_arcs = set((pi[i], pi[i+1]) for i in range(n-1))
        for sigma in perms:
            sigma_arcs = set((sigma[i], sigma[i+1]) for i in range(n-1))
            k = len(pi_arcs & sigma_arcs)
            total += 2**k
    E_H2 = total / 2**(2*(n-1))
    print(f"  n={n}: E[H²] = {E_H2}, Σ 2^k = {total}, n!² = {math.factorial(n)**2}")
    # Ratio E[H²]/n!² / 2^{-2(n-1)}
    print(f"    Σ 2^k / n!² = {total / math.factorial(n)**2:.6f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("DEEP DIVE 9: Paley Self-Complementarity and H")
print("=" * 70)

# Paley tournaments are self-complementary: P_p ≅ P_p^op
# This means H(P_p) = H(P_p^op) (trivially true for ALL T by H(T)=H(T^op))
# But more: the NUMBER of odd k-cycles is the same in P_p and P_p^op

# The number of 3-cycles in P_p:
# t_3(P_p) = p(p-1)(p-3)/24 (for p ≡ 3 mod 4, doubly regular)
# Wait: for ANY regular tournament on p vertices, t_3 = p(p-1)(p-3)/24
# Actually: score sequence is all (p-1)/2
# t_3 = (1/6) Σ_v s_v(s_v - 1) ... no
# For a regular tournament: every vertex has score (p-1)/2
# Number of 3-cycles: t_3 = C(p,3) - Σ_v C(s_v, 2)
# = C(p,3) - p·C((p-1)/2, 2)

for p in [3, 7, 11, 19, 23, 43]:
    if p % 4 != 3:
        continue
    t3 = math.comb(p, 3) - p * math.comb((p-1)//2, 2)
    # Also: E[t_3] = C(p,3)/4 for random tournament
    E_t3 = math.comb(p, 3) / 4
    print(f"  P_{p:2d}: t₃ = {t3}, E[t₃] = {E_t3:.1f}, t₃/E[t₃] = {t3/E_t3:.4f}")
    # For doubly regular: t₃/E[t₃] → ?
    # t₃ = C(p,3) - p·C((p-1)/2, 2)
    # = p(p-1)(p-2)/6 - p·(p-1)/2·((p-1)/2-1)/2
    # = p(p-1)/6 · [p-2 - 3·((p-3)/4)]
    # = p(p-1)/6 · [(4p-8-3p+9)/4]
    # = p(p-1)/6 · (p+1)/4
    # = p(p-1)(p+1)/24
    formula = p*(p-1)*(p+1)//24
    E_t3_int = p*(p-1)*(p-2)//24  # C(p,3)/4 · 4/1 ... no
    # C(p,3)/4 = p(p-1)(p-2)/24
    # t_3/E_t3 = p(p²-1)/24 / (p(p-1)(p-2)/24) = (p+1)/(p-2)
    theoretical_ratio = (p+1)/(p-2)
    print(f"    Formula: t₃ = p(p²-1)/24 = {formula}, ratio = (p+1)/(p-2) = {theoretical_ratio:.4f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("DEEP DIVE 10: Tournament Entropy and π")
print("=" * 70)

# Shannon entropy of the H distribution
# S = -Σ P(H=h) log P(H=h)
# For uniform distribution on {1,3,...,max_h}: S = log(distinct_values)
# How close is the actual entropy to this maximum?

for n in range(3, 7):
    h_values = [count_H(adj, n) for adj in all_tournaments(n)]
    N = len(h_values)
    h_counts = Counter(h_values)

    S = -sum((c/N) * math.log(c/N) for c in h_counts.values())
    S_max = math.log(len(h_counts))

    print(f"\n  n={n}: H entropy = {S:.4f}, max entropy = {S_max:.4f}, ratio = {S/S_max:.4f}")
    print(f"    #distinct H = {len(h_counts)}, equiv uniform = e^S = {math.exp(S):.2f}")
    # The ratio S/log(2^m) gives the fraction of "tournament information" in H
    m = n*(n-1)//2
    print(f"    H entropy / m·log2 = {S/(m*math.log(2)):.4f}")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("DEEP DIVE 11: The H-Spectrum of Paley P_7")
print("=" * 70)

# For P_7 we know H = 189 = IP(G(P_7), 2)
# Break this down: how many independent sets of each size?
# From THM-209: H = 1 + 2·(t_3 + t_5 + t_7) + 4·d_{33}
# We know t_3(P_7) = 7·8/24 = ... let me compute

p = 7
adj = paley_tournament(p)

# Count 3-cycles
t3 = 0
for combo in itertools.combinations(range(p), 3):
    a, b, c = combo
    if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
        t3 += 1

# Count 5-cycles
t5 = 0
for combo in itertools.combinations(range(p), 5):
    for perm in itertools.permutations(combo):
        is_cycle = all(adj[perm[i]][perm[(i+1)%5]] for i in range(5))
        if is_cycle:
            t5 += 1
t5 //= 5  # Each cycle counted 5 times

# Count 7-cycles
t7 = 0
for perm in itertools.permutations(range(p)):
    if all(adj[perm[i]][perm[(i+1)%p]] for i in range(p)):
        t7 += 1
t7 //= p

# Count d33: disjoint pairs of 3-cycles
all_3cycles = []
for combo in itertools.combinations(range(p), 3):
    a, b, c = combo
    if adj[a][b] and adj[b][c] and adj[c][a]:
        all_3cycles.append((a, b, c))
    if adj[a][c] and adj[c][b] and adj[b][a]:
        all_3cycles.append((a, c, b))

d33 = 0
for i in range(len(all_3cycles)):
    for j in range(i+1, len(all_3cycles)):
        s1 = set(all_3cycles[i])
        s2 = set(all_3cycles[j])
        if not s1 & s2:
            d33 += 1

H_formula = 1 + 2*(t3 + t5 + t7) + 4*d33
H_actual = count_H(adj, p)

print(f"\nP_7 decomposition:")
print(f"  t_3 = {t3}")
print(f"  t_5 = {t5}")
print(f"  t_7 = {t7}")
print(f"  d_33 = {d33}")
print(f"  H = 1 + 2·({t3}+{t5}+{t7}) + 4·{d33}")
print(f"    = 1 + 2·{t3+t5+t7} + {4*d33}")
print(f"    = {H_formula}")
print(f"  Direct: H = {H_actual}")
print(f"  Match: {H_formula == H_actual}")

# Independence numbers of G(P_7)
# Total odd cycles = t3 + t5 + t7
total_cycles = t3 + t5 + t7
print(f"\n  Total directed odd cycles = {total_cycles}")
print(f"  These are the VERTICES of the odd-cycle disjointness graph G(P_7)")

# The f-vector of the independence complex of G(P_7):
# f_0 = 1 (empty set)
# f_1 = total_cycles (singletons)
# f_2 = d33 (disjoint pairs) -- but also d35, d55, d37, ...

# Wait — d33 only counts disjoint PAIRS OF 3-CYCLES.
# For a full account at n=7, we need d35 and d55 too.
# Actually at n=7, a disjoint pair of 5-cycles would need 10 > 7 vertices.
# So d55 = 0 at n=7. And d35 needs 8 > 7 vertices. So d35 = 0.
# And d37 needs 10 > 7. So only d33 contributes at n=7.
# At level 3: d333 would need 9 > 7. So d333 = 0.
# Therefore H = 1 + 2·(t3+t5+t7) + 4·d33 is COMPLETE for n=7!

print(f"\n  At n=7: d35 = d55 = d37 = d333 = ... = 0 (vertex count exceeds 7)")
print(f"  H = 1 + 2·total_cycles + 4·d33 is EXACT for n=7")
print(f"  The independence complex of G(P_7) has f-vector: (1, {total_cycles}, {d33})")

# ═══════════════════════════════════════════════════════════════
print("\n" + "=" * 70)
print("DEEP DIVE 12: Angle Sums and π")
print("=" * 70)

# Place tournament vertices on a unit circle at angles 2πk/n
# For each Hamiltonian path π, compute the "total turning angle"
# = Σ_{i=0}^{n-2} angle(π(i) → π(i+1))
# where angle(j → k) = 2π(k-j)/n mod 2π

# For the cyclic tournament, paths that follow the circle have small angles
# For random tournaments, the angles should average out

for n in [5, 7]:
    adj_paley = None
    if n == 7:
        adj_paley = paley_tournament(n)

    adj_cyclic = [[0]*n for _ in range(n)]
    forward = set(range(1, n//2 + 1))
    for i in range(n):
        for j in range(n):
            if i != j and (j - i) % n in forward:
                adj_cyclic[i][j] = 1

    for label, adj in [("Cyclic", adj_cyclic)] + ([("Paley", adj_paley)] if adj_paley else []):
        # Find all Hamiltonian paths and compute their total angle
        total_angle_sum = 0
        path_count = 0
        for perm in itertools.permutations(range(n)):
            is_hp = all(adj[perm[i]][perm[i+1]] for i in range(n-1))
            if is_hp:
                # Total turning angle
                total_angle = 0
                for i in range(n-1):
                    angle = (2 * math.pi * (perm[i+1] - perm[i]) % n) / n
                    if angle < 0:
                        angle += 2 * math.pi / n * n
                    total_angle += angle
                total_angle_sum += total_angle
                path_count += 1

        if path_count > 0:
            mean_angle = total_angle_sum / path_count
            # Expected: each arc goes to one of the n-1 possible vertices
            # Mean arc angle for random = π (middle of [0, 2π])
            # Mean total angle for random HP = (n-1)·π
            print(f"\n  {label} T_{n}: H = {path_count}")
            print(f"    Mean total angle = {mean_angle:.4f}")
            print(f"    (n-1)·π = {(n-1)*math.pi:.4f}")
            print(f"    Mean total angle / ((n-1)·π) = {mean_angle/((n-1)*math.pi):.4f}")

print("\n" + "=" * 70)
print("END OF DEEP DIVES — opus-2026-03-14-S89c")
print("=" * 70)
