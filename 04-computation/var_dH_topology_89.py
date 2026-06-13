#!/usr/bin/env python3
"""
Var(dH) SEQUENCE AND TOPOLOGICAL CONNECTIONS
opus-2026-03-14-S89

The sequence Var(dH) = E[dH^2] = 2, 4, 15, 84 for n=3,4,5,6.
What IS this sequence? What does it count?

Also: exploring the H-level set graph, the "tournament flow", and
the connection between cone towers and simplicial homology.
"""

from itertools import permutations
from math import factorial, comb
from collections import Counter, defaultdict
from fractions import Fraction

def compute_H(n, adj):
    count = 0
    for p in permutations(range(n)):
        ok = True
        for k in range(n-1):
            if adj[(p[k], p[k+1])] != 1:
                ok = False
                break
        if ok:
            count += 1
    return count

print("=" * 70)
print("Var(dH) SEQUENCE AND TOPOLOGICAL CONNECTIONS")
print("opus-2026-03-14-S89")
print("=" * 70)

# ======================================================================
# PART 1: THE Var(dH) SEQUENCE
# ======================================================================
print("\n" + "=" * 70)
print("PART 1: WHAT IS Var(dH) = 2, 4, 15, 84?")
print("=" * 70)

# Var(dH) = E[dH^2] where the expectation is over all T and all edges e.
# From Walsh: E_T[dH_e^2] = 4 * sum_{S: e in S} H_hat[S]^2
# By symmetry over edges: E_{T,e}[dH^2] = 4/m * sum_S |S| * H_hat[S]^2

# But wait: sum_S H_hat[S]^2 = Var(H) + Mean(H)^2 = E[H^2]
# Actually: sum_S H_hat[S]^2 = E[H^2] / 2^m? No.
# Parseval: sum_S H_hat[S]^2 = E[H^2] (with our normalization)

# More carefully:
# H(T) = sum_S H_hat[S] chi_S(T)
# E[H^2] = sum_S H_hat[S]^2 (Parseval)
# E_T[dH_e^2] = 4 * sum_{S: e in S} H_hat[S]^2
# E_{T,e}[dH^2] = (1/m) * sum_e 4 * sum_{S: e in S} H_hat[S]^2
#                = 4/m * sum_S |S| * H_hat[S]^2

# Now: Mean(H) = H_hat[empty set] = n!/2^{n-1}
# Mean(H)^2 for n=3: 9/4, n=4: 9, n=5: 225/4, n=6: 2025/4

# E[H^2] = Mean^2 + Var = Mean^2 * (1 + Var/Mean^2)

# Let me just compute it from the grand formula.
# E[dH^2] = 4/m * sum_S |S| * H_hat[S]^2
# The level-k contribution: sum_{|S|=k} H_hat[S]^2 = E_k = (energy at level k)
# E[dH^2] = 4/m * sum_k k * E_k

# From the grand formula: E_{2k}/E_0 = 2*(n-2k)^k / P(n,2k)
# And E_0 = Mean^2 = (n!/2^{n-1})^2
# So: E[dH^2] = 4/m * sum_k (2k) * E_0 * E_{2k}/E_0
#             = 8*E_0/m * sum_k k * 2*(n-2k)^k / P(n,2k)

# Wait, I should include E_0 contribution: level 0 has |S|=0, contributes 0.
# And odd levels vanish. So:
# E[dH^2] = 4/m * sum_{k=1}^{K} (2k) * E_{2k}
# where E_{2k} = E_0 * 2*(n-2k)^k / P(n,2k)

print("\n  The Var(dH) sequence (empirical):")
var_dH = {3: 2, 4: 4, 5: 15, 6: 84}
for n, v in var_dH.items():
    print(f"    n={n}: Var(dH) = {v}")

print("\n  Can we derive these from the grand formula?")
for n in range(3, 9):
    m = n * (n - 1) // 2
    K = (n - 1) // 2
    mean_H = Fraction(factorial(n), 2**(n-1))
    E0 = mean_H ** 2

    # Var(dH) = 4/m * sum_{k=1}^K (2k) * E_0 * 2*(n-2k)^k / P(n,2k)
    total = Fraction(0)
    for k in range(1, K + 1):
        j = n - 2*k  # = n - 2k
        P_n_2k = Fraction(factorial(n), factorial(j))
        E_2k_over_E0 = Fraction(2 * j**k, 1) / P_n_2k
        total += Fraction(2*k) * E_2k_over_E0

    var_formula = Fraction(4, m) * E0 * total

    print(f"    n={n}: formula = 4/{m} * {E0} * {float(total):.6f} = {float(var_formula):.4f}")
    if n in var_dH:
        print(f"           empirical = {var_dH[n]}, match = {float(var_formula) == var_dH[n]}")

print("\n  WAIT — the formula gives Var(dH) as a function of n.")
print("  Let me compute the simplified formula.")

for n in range(3, 12):
    m = n * (n - 1) // 2
    K = (n - 1) // 2
    mean_sq = Fraction(factorial(n)**2, 4**(n-1))

    # sum_k (2k) * 2*(n-2k)^k / P(n,2k)
    s = Fraction(0)
    for k in range(1, K + 1):
        j = n - 2*k
        P = Fraction(factorial(n), factorial(j))
        s += Fraction(4*k * j**k, 1) / P

    var_dH_formula = mean_sq * s / Fraction(m, 1)

    print(f"    n={n}: Var(dH) = {float(var_dH_formula):.6f}", end="")
    # Check if it's an integer
    if var_dH_formula.denominator == 1:
        print(f" = {var_dH_formula.numerator} (integer!)", end="")
    else:
        print(f" = {var_dH_formula}", end="")
    print()

# ======================================================================
# PART 2: OEIS LOOKUP FOR 2, 4, 15, 84
# ======================================================================
print("\n" + "=" * 70)
print("PART 2: SEQUENCE IDENTIFICATION")
print("=" * 70)

# 2, 4, 15, 84
# 84 = 4 * 21 = C(9,2) = 36? No. 84 = C(9,2) = 36? No, C(9,2)=36.
# 84 = 7 * 12 = C(7,1)*C(4,2)?
# 84 = 4!/4 * 14 = ... hmm
# Let me check ratios: 4/2=2, 15/4=3.75, 84/15=5.6
# Not clean ratios.

# Try: 2 = 2, 4 = 4, 15 = 15, 84 = 84
# 2 = 1*2, 4 = 2*2, 15 = 3*5, 84 = 4*21 = 4*C(7,2)
# Or: 2 = C(3,2)-1, 4 = C(4,2)-2, 15 = C(5,2)+5, 84 = C(6,2)+69. No pattern.

# Check: n!/2^(n-2):
# n=3: 6/2 = 3 (not 2)
# n=4: 24/4 = 6 (not 4)
# Not that.

# Maybe: (n-1)! - something?
# n=3: 2! = 2 = 2. Match!
# n=4: 3! = 6 != 4.
# n=5: 4! = 24 != 15.
# No.

# (n-1)*(n-2)/2:
# n=3: 1, n=4: 3, n=5: 6, n=6: 10. No.

# OK let me think differently.
# From the formula: Var(dH) = 4*Mean^2/m * sum_k k * E_{2k}/E_0
# = 4*Mean^2/m * (1*E_2/E_0 + 2*E_4/E_0 + ...)
# = 4*Mean^2/m * (1*(2(n-2)/(n(n-1))) + ...)

# For n=3: m=3, Mean^2=9/4
# sum = 1 * 2*1/(3*2) = 1/3
# Var(dH) = 4*(9/4)/3 * 1/3 = (9/3) * 1/3 = 3/3 = 1??? That's not 2.

# Wait, I think I'm confusing normalization. Let me redo.
# The empirical Var(dH) averages over both T and e.
# E_{T,e}[dH^2] = 1/m * sum_e E_T[dH_e^2]
# E_T[dH_e^2] = 4 * sum_{S: e in S} (H_hat[S])^2
# where H_hat[S] = (1/2^m) sum_T H(T) chi_S(T)
# Parseval: sum_S (H_hat[S])^2 = (1/2^m) sum_T H(T)^2 = E[H^2]

# Actually let me re-derive. H = sum_S c_S chi_S where c_S = H_hat[S].
# <H, H> = sum_T H(T)^2 = 2^m * sum_S c_S^2
# So E[H^2] = (1/2^m) sum_T H^2 = sum_S c_S^2

# dH_e(T) = 2 * sum_{S: e in S} c_S chi_S(T)
# dH_e^2(T) = 4 * (sum_{S: e in S} c_S chi_S(T))^2
# E_T[dH_e^2] = 4 * sum_{S,S': e in S, e in S'} c_S c_{S'} E[chi_S chi_{S'}]
#             = 4 * sum_{S: e in S} c_S^2  (orthogonality)

# E_{T,e}[dH^2] = (1/m) sum_e 4 sum_{S: e in S} c_S^2
#               = 4/m * sum_S |S| c_S^2

# Now: sum_S |S| c_S^2 = sum_{k=0}^m k * sum_{|S|=k} c_S^2
#                       = sum_{k=1}^m k * E_k
# where E_k = level-k energy = sum_{|S|=k} c_S^2

# Only even k: E_k = E_0 * formula for k = 2,4,...
# E_0 = c_emptyset^2 = Mean(H)^2 (level 0 has weight 0, doesn't contribute)

# Hmm I got a different thing earlier. Let me just compute from scratch for n=3.

print("\n  Direct verification for n=3:")
# n=3, m=3
# H values: T=000:1, 001:1, 010:3, 011:1, 100:1, 101:3, 110:1, 111:1
H3 = [1, 1, 3, 1, 1, 3, 1, 1]
N = 8
m3 = 3

# Walsh transform
c = [0.0] * N
for S in range(N):
    for T in range(N):
        chi = (-1) ** bin(S & T).count('1')
        c[S] += H3[T] * chi / N
print(f"  Walsh coefficients c_S: {[round(x,4) for x in c]}")
print(f"  c_S^2: {[round(x**2, 4) for x in c]}")

# sum_S |S| c_S^2
total = 0
for S in range(N):
    wt = bin(S).count('1')
    total += wt * c[S]**2
print(f"  sum_S |S| * c_S^2 = {total:.4f}")
print(f"  4/m * this = {4/m3 * total:.4f}")
print(f"  Empirical Var(dH) at n=3 = 2")

# Hmm, 4/3 * 1.5 = 2. YES!
# sum_S |S| * c_S^2 = level 2 contribution = 2 * sum_{|S|=2} c_S^2
# H_hat = [1.5, 0, 0, -0.5, 0, 0.5, -0.5, 0]
# |S|=0: S=000, c=1.5, |S|*c^2 = 0
# |S|=1: S=001,010,100, c = 0,0,0, contribution = 0
# |S|=2: S=011,101,110, c = -0.5, 0.5, -0.5
#   contribution = 2*(0.25 + 0.25 + 0.25) = 2*0.75 = 1.5
# |S|=3: S=111, c = 0, contribution = 0
# Total: 1.5
# 4/3 * 1.5 = 2. Matches!

print("\n  Detailed calculation:")
for S in range(N):
    wt = bin(S).count('1')
    contrib = wt * c[S]**2
    if contrib != 0:
        print(f"    S={S:03b} (|S|={wt}): c_S={c[S]:.4f}, "
              f"|S|*c_S^2 = {contrib:.4f}")

# ======================================================================
# PART 3: EXACT FORMULA FOR Var(dH)
# ======================================================================
print("\n" + "=" * 70)
print("PART 3: EXACT FORMULA FOR Var(dH)")
print("=" * 70)

# Var(dH) = 4/m * sum_{k>=1} k * E_k
# Only even k contribute: k = 2, 4, 6, ...
# E_{2j} = C(m, 2j) * (Mean^2) * E_{2j}/E_0 / C(m, 2j)
# Wait, that's circular. Let me use:
# E_{2j} = level-2j energy = (E_{2j}/E_0) * E_0
# where E_0 = Mean^2 and E_{2j}/E_0 = 2*(n-2j)^j / P(n,2j)

# Actually: E_{2j} = sum_{|S|=2j} c_S^2
# And level energy per coefficient: E_{2j}/C(m,2j)
# The grand formula says: E_{2j}/E_0 = 2*(n-2j)^j / P(n,2j)
# where this is the TOTAL level energy ratio.

# So: sum_{k>=1} k * E_k = sum_{j>=1} (2j) * E_{2j}
#                         = sum_{j>=1} 2j * E_0 * 2*(n-2j)^j / P(n,2j)
#                         = 4*E_0 * sum_{j>=1} j*(n-2j)^j / P(n,2j)

# Therefore: Var(dH) = 4/m * 4*E_0 * sum_j j*(n-2j)^j / P(n,2j)
#                     = 16*Mean^2/m * sum_j j*(n-2j)^j / P(n,2j)

print("\n  Var(dH) = 16*Mean^2/m * sum_j j*(n-2j)^j / P(n,2j)")
print()

for n in range(3, 10):
    m = n * (n - 1) // 2
    K = (n - 1) // 2
    mean = Fraction(factorial(n), 2**(n-1))
    mean_sq = mean * mean

    s = Fraction(0)
    for j in range(1, K + 1):
        r = n - 2*j
        P = Fraction(factorial(n), factorial(r))
        s += Fraction(j * r**j, 1) / P

    var = Fraction(16, 1) * mean_sq * s / Fraction(m, 1)

    print(f"  n={n}: Mean={float(mean):.4f}, m={m}, "
          f"sum={float(s):.6f}, Var(dH) = {float(var):.6f}", end="")
    if var.denominator == 1:
        print(f" = {var.numerator}", end="")
    print()

# ======================================================================
# PART 4: SIMPLIFY THE FORMULA
# ======================================================================
print("\n" + "=" * 70)
print("PART 4: SIMPLIFIED Var(dH)")
print("=" * 70)

# Var(dH) = 16 * (n!/2^{n-1})^2 / (n(n-1)/2) * sum_j j*(n-2j)^j / P(n,2j)
# = 32 * (n!)^2 / (4^{n-1} * n * (n-1)) * sum_j j*(n-2j)^j / (n!/(n-2j)!)
# = 32 * n! / (4^{n-1} * n * (n-1)) * sum_j j * (n-2j)^j * (n-2j)!

# Let me define: S(n) = sum_{j=1}^{K} j * (n-2j)^j * (n-2j)!

print("\n  Define S(n) = sum_{j=1}^K j * (n-2j)^j * (n-2j)!")
print("  Then Var(dH) = 32 * n! * S(n) / (4^{n-1} * n * (n-1))")
print()

for n in range(3, 12):
    K = (n - 1) // 2
    S_n = 0
    terms = []
    for j in range(1, K + 1):
        r = n - 2*j
        term = j * r**j * factorial(r)
        S_n += term
        terms.append(f"{j}*{r}^{j}*{r}! = {term}")

    var_exact = Fraction(32 * factorial(n) * S_n, 4**(n-1) * n * (n-1))

    print(f"  n={n}: S({n}) = {S_n}")
    for t in terms:
        print(f"         {t}")
    print(f"         Var(dH) = 32*{factorial(n)}*{S_n} / (4^{n-1}*{n}*{n-1}) "
          f"= {float(var_exact):.6f}")
    if var_exact.denominator == 1:
        print(f"         = {var_exact.numerator} (integer)")

# ======================================================================
# PART 5: THE SEQUENCE S(n) AND ITS PROPERTIES
# ======================================================================
print("\n" + "=" * 70)
print("PART 5: S(n) SEQUENCE AND RATIOS")
print("=" * 70)

S_vals = []
for n in range(3, 15):
    K = (n - 1) // 2
    S_n = 0
    for j in range(1, K + 1):
        r = n - 2*j
        S_n += j * r**j * factorial(r)
    S_vals.append(S_n)
    print(f"  S({n:2d}) = {S_n}")

print("\n  S(n) ratios:")
for i in range(1, len(S_vals)):
    n = i + 3
    print(f"  S({n})/S({n-1}) = {S_vals[i]/S_vals[i-1]:.4f}")

# ======================================================================
# PART 6: THE TOPOLOGY — H AS A MORSE FUNCTION
# ======================================================================
print("\n" + "=" * 70)
print("PART 6: H AS A MORSE FUNCTION ON THE TOURNAMENT HYPERCUBE")
print("=" * 70)

print("""
  View H as a function on the tournament hypercube Q_m.
  Q_m has 2^m vertices connected by m edges at each vertex.

  A vertex T is a LOCAL MAXIMUM if H(T) >= H(T^e) for all edges e.
  A vertex T is a LOCAL MINIMUM if H(T) <= H(T^e) for all edges e.
  Otherwise it's a SADDLE.

  In Morse theory, the topology of the hypercube is determined by
  the critical points (local max, min, saddle) of any "Morse function."
""")

for n in range(3, 7):
    m = n * (n - 1) // 2
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))

    H_list = [0] * (2**m)
    for bits in range(2**m):
        adj = {}
        for idx, (i, j) in enumerate(edges):
            if (bits >> idx) & 1:
                adj[(i,j)] = 1
                adj[(j,i)] = 0
            else:
                adj[(i,j)] = 0
                adj[(j,i)] = 1
        H_list[bits] = compute_H(n, adj)

    # Find local maxima, minima, saddles
    local_max = []
    local_min = []
    saddle = []

    for bits in range(2**m):
        is_max = True
        is_min = True
        for e_idx in range(m):
            nbr = bits ^ (1 << e_idx)
            if H_list[nbr] > H_list[bits]:
                is_max = False
            if H_list[nbr] < H_list[bits]:
                is_min = False
        if is_max:
            local_max.append((bits, H_list[bits]))
        elif is_min:
            local_min.append((bits, H_list[bits]))
        else:
            saddle.append((bits, H_list[bits]))

    print(f"\n  n={n} (m={m}):")
    print(f"    Local maxima: {len(local_max)}")
    max_h_vals = Counter(h for _, h in local_max)
    for h, cnt in sorted(max_h_vals.items()):
        print(f"      H={h}: {cnt} tournaments")
    print(f"    Local minima: {len(local_min)}")
    min_h_vals = Counter(h for _, h in local_min)
    for h, cnt in sorted(min_h_vals.items()):
        print(f"      H={h}: {cnt} tournaments")
    print(f"    Saddles: {len(saddle)}")

    # Morse inequality: #max - #saddle + #min = chi(Q_m)
    # chi(Q_m) = 0 for m >= 1 (Euler characteristic of hypercube = 0)
    # But this is the weak Morse inequality: #crit >= sum b_i
    # For Q_m: b_i = C(m,i), sum = 2^m
    print(f"    #max + #min + #saddle = {len(local_max) + len(local_min) + len(saddle)} = 2^{m}")

    # The INDEX of a critical point = number of descending directions
    # For a local max: index = m (all descending)
    # For a local min: index = 0 (all ascending)
    # For a saddle: 0 < index < m

    # Compute index for each tournament
    index_dist = Counter()
    for bits in range(2**m):
        descending = 0
        for e_idx in range(m):
            nbr = bits ^ (1 << e_idx)
            if H_list[nbr] < H_list[bits]:
                descending += 1
        index_dist[descending] += 1

    print(f"    Morse index distribution:")
    for idx in sorted(index_dist.keys()):
        frac = index_dist[idx] / (2**m)
        print(f"      index {idx:2d}: {index_dist[idx]:6d} ({frac:.4f})")

# ======================================================================
# PART 7: THE GRADIENT FLOW
# ======================================================================
print("\n" + "=" * 70)
print("PART 7: GRADIENT FLOW ON THE TOURNAMENT HYPERCUBE")
print("=" * 70)

print("""
  Define the GRADIENT FLOW: from each tournament T, follow the
  steepest ascent direction (flip the edge that increases H the most).

  This gives a FLOW on the hypercube that terminates at local maxima.
  The basin of attraction of each local maximum is a "cell" of the
  tournament space.
""")

for n in range(3, 7):
    m = n * (n - 1) // 2
    edges = []
    for i in range(n):
        for j in range(i+1, n):
            edges.append((i, j))

    H_list = [0] * (2**m)
    for bits in range(2**m):
        adj = {}
        for idx, (i, j) in enumerate(edges):
            if (bits >> idx) & 1:
                adj[(i,j)] = 1
                adj[(j,i)] = 0
            else:
                adj[(i,j)] = 0
                adj[(j,i)] = 1
        H_list[bits] = compute_H(n, adj)

    # Steepest ascent flow
    basin = {}  # maps tournament -> local max it flows to
    for bits in range(2**m):
        current = bits
        path = [current]
        while True:
            best_nbr = current
            best_h = H_list[current]
            for e_idx in range(m):
                nbr = current ^ (1 << e_idx)
                if H_list[nbr] > best_h:
                    best_h = H_list[nbr]
                    best_nbr = nbr
            if best_nbr == current:
                break  # local max
            current = best_nbr
            if current in path:
                break  # cycle (shouldn't happen with strict increase)
            path.append(current)
        basin[bits] = (current, H_list[current])

    # Count basin sizes
    basin_sizes = Counter(basin.values())
    print(f"\n  n={n}: Gradient flow basins (steepest ascent):")
    for (max_t, max_h), size in sorted(basin_sizes.items(), key=lambda x: -x[1]):
        print(f"    Local max H={max_h}: basin size = {size} ({size/2**m:.4f})")

# ======================================================================
# PART 8: COMPLEMENT SYMMETRY IN THE MORSE PICTURE
# ======================================================================
print("\n" + "=" * 70)
print("PART 8: COMPLEMENT SYMMETRY AND MORSE STRUCTURE")
print("=" * 70)

print("""
  THM-203: H(T) = H(T_bar) where T_bar = complement (all arcs reversed).
  In the hypercube: T_bar = T XOR 111...1 (flip ALL bits).

  This means: if T is a local max, so is T_bar.
  The Morse function H has a Z/2Z symmetry (the complement involution).

  The basins of T and T_bar have the SAME size (by symmetry).
  Every local max comes in a PAIR (T, T_bar), unless T = T_bar.
  But T = T_bar only if m is even and T is "self-complementary."

  For n=3: T_bar = T XOR 111 = 7 - T.
    T=010 (C_3, H=3): T_bar=101 (also C_3, H=3). These are the two 3-cycles.
    T=000 (trans, H=1): T_bar=111 (trans, H=1). Complement pairs.

  The complement involution is an ISOMETRY of the hypercube
  (distance-preserving, since it flips all bits).
""")

for n in range(3, 7):
    m = n * (n - 1) // 2
    complement_mask = (1 << m) - 1

    # Check: is there a self-complementary tournament?
    self_comp = []
    for bits in range(2**m):
        if bits == bits ^ complement_mask:
            self_comp.append(bits)

    print(f"\n  n={n} (m={m}):")
    if self_comp:
        print(f"    Self-complementary tournaments: {len(self_comp)}")
    else:
        print(f"    No self-complementary tournaments (m={m} is {'even' if m%2==0 else 'odd'})")
        # Note: self-comp exists only if m is even? No, T != T_bar always for tournaments.
        # Actually for tournaments, self-complementary means T is isomorphic to T_bar,
        # not that T = T_bar bitwise. Bitwise T = T XOR mask means all arcs flip,
        # which reverses all directions. T = T_bar bitwise means the tournament
        # is its own complement, which is impossible since the identity permutation
        # maps T to T_bar only if all edges are "undirected" which they're not.
        # For labeled tournaments: T = T_bar (bitwise) requires each bit = complement bit,
        # which is impossible.
    print(f"    All tournaments pair up: {2**m} = 2 * {2**(m-1)} pairs")

print("\n" + "=" * 70)
print("DONE -- Var(dH) TOPOLOGY AND MORSE THEORY")
print("=" * 70)
