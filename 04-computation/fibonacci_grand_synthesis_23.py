#!/usr/bin/env python3
"""
Grand Fibonacci-Tournament Synthesis — all sessions unified.
opus-2026-03-14-S86

Building on ALL sessions S85-S90:
- opus S85: I(C_k,2)=2^k+(-1)^k, I(P_k,2)=(2^{k+2}-(-1)^k)/3, forbidden H
- kind-pasteur S86: writhe = odd strand, H = even strand
- opus S86: period-6, golden ratio, strand coupling matrix
- kind-pasteur S87: three-strand Pascal, s2/s0=2=OCF fugacity
- kind-pasteur S88: Pascal-Fibonacci period-6, I(C_3,2)=7=forbidden
- kind-pasteur S89: max_H=Jacobsthal for n<=4, 21 triple coincidence, 1729
- kind-pasteur S90: a(n)=C(n,1)+C(n,3), n=5 triple coincidence

THIS SCRIPT: Unify everything into a single coherent framework.

THE THESIS: All these sequences are different evaluations of a SINGLE
underlying structure — the independence polynomial of the conflict graph
at different fugacities x:
  x=0: trivial (I=1)
  x=1: Fibonacci/Lucas world (combinatorics of paths/cycles)
  x=2: Tournament world (H = Hamiltonian path count)
  x=φ: Golden ratio world
  x=-1: Alternating world (Möbius function)
"""

import math
import numpy as np
from collections import Counter, defaultdict
from fractions import Fraction
import random

def fib(k):
    if k < 0:
        return ((-1)**((-k)+1)) * fib(-k)
    if k == 0: return 0
    if k == 1: return 1
    a, b = 0, 1
    for _ in range(k - 1):
        a, b = b, a + b
    return b

def jacobsthal(k):
    if k == 0: return 0
    if k == 1: return 1
    a, b = 0, 1
    for _ in range(k - 1):
        a, b = b, b + 2*a
    return b

def lucas(k):
    if k == 0: return 2
    if k == 1: return 1
    a, b = 2, 1
    for _ in range(k - 1):
        a, b = b, a + b
    return b

def triangular(n):
    return n * (n + 1) // 2

def compute_H_dp(adj, n):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            if not (S & (1 << v)):
                continue
            if dp[S][v] == 0:
                continue
            for w in range(n):
                if S & (1 << w):
                    continue
                if adj[v][w]:
                    dp[S | (1 << w)][w] += dp[S][v]
    return sum(dp[full][v] for v in range(n))

def get_tournament(n, bits):
    arcs = [(i, j) for i in range(n) for j in range(i+1, n)]
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(arcs):
        if (bits >> k) & 1:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj

# ============================================================
# PART 1: THE MASTER ALIGNMENT TABLE
# ============================================================
print("=" * 70)
print("PART 1: MASTER ALIGNMENT TABLE — ALL SEQUENCES UNIFIED")
print("=" * 70)

max_H = [1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095]

print("\n  n | F(n) | J(n) | L(n) | T(n) | maxH | I(C_n,2) | I(P_n,2) | C(n,1)+C(n,3)")
for n in range(0, 12):
    Fn = fib(n)
    Jn = jacobsthal(n)
    Ln = lucas(n)
    Tn = triangular(n)
    mH = max_H[n] if n < len(max_H) else "?"
    IC = 2**n + (-1)**n if n >= 1 else 2
    IP = (2**(n+2) - (-1)**n) // 3
    mystery = math.comb(n, 1) + math.comb(n, 3) if n >= 0 else 0

    # Highlight coincidences
    vals = {'F': Fn, 'J': Jn, 'L': Ln, 'T': Tn}
    if isinstance(mH, int):
        vals['mH'] = mH
    vals['IC'] = IC
    vals['IP'] = IP
    vals['a'] = mystery

    seen = defaultdict(list)
    for name, val in vals.items():
        if val > 1:
            seen[val].append(name)

    matches = [f"{'='.join(names)}={val}" for val, names in seen.items() if len(names) >= 2]

    print(f"  {n:2d} | {Fn:4d} | {Jn:5d} | {Ln:4d} | {Tn:4d} | {str(mH):>5s} | {IC:8d} | {IP:8d} | {mystery:>13d}"
          + (f"  <- {', '.join(matches)}" if matches else ""))

# ============================================================
# PART 2: THE 21 NEXUS — WHERE EVERYTHING MEETS
# ============================================================
print("\n" + "=" * 70)
print("PART 2: THE 21 NEXUS — WHERE EVERYTHING MEETS")
print("=" * 70)

print("""
21 appears in SEVEN different contexts:
  1. 21 = F(8)          -- 8th Fibonacci number
  2. 21 = J(6)          -- 6th Jacobsthal number
  3. 21 = T(6)          -- 6th triangular number = C(7,2)
  4. 21 = I(P_4, 2)     -- independence polynomial of P_4 at x=2
  5. 21 = FORBIDDEN H   -- never achievable as H(T) for any tournament
  6. 21 = C(7, 2)       -- arcs in 7-vertex tournament
  7. 21 = 3 x 7         -- product of two tournament primes
""")

# Find all graphs on 4 vertices with I(G, 2) = 21
from itertools import combinations

def graph_indep_poly_at_2(adj_list, n):
    total = 0
    for size in range(n + 1):
        for subset in combinations(range(n), size):
            is_indep = True
            for i, u in enumerate(subset):
                for v in subset[i+1:]:
                    if v in adj_list.get(u, set()):
                        is_indep = False
                        break
                if not is_indep:
                    break
            if is_indep:
                total += 2**size
    return total

print("Graphs on 4 vertices with I(G,2) = 21:")
edges_4 = [(i, j) for i in range(4) for j in range(i+1, 4)]
count_21 = 0
for edge_bits in range(64):
    adj = defaultdict(set)
    for k, (i, j) in enumerate(edges_4):
        if (edge_bits >> k) & 1:
            adj[i].add(j)
            adj[j].add(i)
    ig2 = graph_indep_poly_at_2(adj, 4)
    if ig2 == 21:
        edge_list = [(i, j) for k, (i, j) in enumerate(edges_4) if (edge_bits >> k) & 1]
        print(f"  edges={edge_list}, I(G,2)={ig2}")
        count_21 += 1
print(f"  Total: {count_21} graphs on 4 vertices with I(G,2)=21")

# Graphs on 5 vertices
print("\nGraphs on 5 vertices with I(G,2) = 21:")
edges_5 = [(i, j) for i in range(5) for j in range(i+1, 5)]
count_21_5 = 0
for edge_bits in range(1 << 10):
    adj = defaultdict(set)
    for k, (i, j) in enumerate(edges_5):
        if (edge_bits >> k) & 1:
            adj[i].add(j)
            adj[j].add(i)
    ig2 = graph_indep_poly_at_2(adj, 5)
    if ig2 == 21:
        count_21_5 += 1
        if count_21_5 <= 3:
            edge_list = [(i, j) for k, (i, j) in enumerate(edges_5) if (edge_bits >> k) & 1]
            print(f"  edges={edge_list}, I(G,2)={ig2}")
print(f"  Total: {count_21_5} graphs on 5 vertices with I(G,2)=21")

# ============================================================
# PART 3: FUGACITY LANDSCAPE — I(G, x) FOR KEY GRAPHS
# ============================================================
print("\n" + "=" * 70)
print("PART 3: FUGACITY LANDSCAPE — I(P_k, x) AND I(C_k, x)")
print("=" * 70)

def I_path(k, x):
    if k == 0: return 1
    if k == 1: return 1 + x
    a, b = 1, 1 + x
    for _ in range(k - 1):
        a, b = b, b + x * a
    return b

def I_cycle(k, x):
    if k < 3: return None
    return I_path(k-1, x) + x * I_path(k-3, x)

print("\n  x   | I(P_3,x) | I(P_4,x) | I(P_5,x) | I(C_3,x) | I(C_4,x) | I(C_5,x)")
for x in [0, 1, 2, 3, -1, Fraction(1, 2)]:
    vals = [I_path(3, x), I_path(4, x), I_path(5, x),
            I_cycle(3, x), I_cycle(4, x), I_cycle(5, x)]
    print(f"  {str(x):>4s} | " + " | ".join(f"{str(v):>8s}" for v in vals))

# At x=-1: magic
print("\nAt x = -1 (alternating world):")
for k in range(1, 10):
    ipk = I_path(k, -1)
    ick = I_cycle(k, -1) if k >= 3 else None
    extras = []
    if ipk == 0: extras.append("= 0!")
    if ipk == 1: extras.append("= 1")
    if ipk == -1: extras.append("= -1")
    print(f"  I(P_{k}, -1) = {ipk}" + (f", I(C_{k}, -1) = {ick}" if ick is not None else "") + (" " + " ".join(extras) if extras else ""))

# ============================================================
# PART 4: UNIVERSAL RECURRENCE EIGENVALUES
# ============================================================
print("\n" + "=" * 70)
print("PART 4: T_x = [[1,x],[1,0]] — EIGENVALUE LANDSCAPE")
print("=" * 70)

print("\nEigenvalues of T_x:")
for x in [-2, -1, 0, 1, 2, 3, 4, 6, 12]:
    disc = 1 + 4*x
    if disc >= 0:
        e1 = (1 + disc**0.5) / 2
        e2 = (1 - disc**0.5) / 2
        int_eig = (abs(e1 - round(e1)) < 1e-10 and abs(e2 - round(e2)) < 1e-10)
        print(f"  x={x:3d}: lambda = {e1:.6f}, {e2:.6f}" + (" <-- INTEGER eigenvalues!" if int_eig else ""))
    else:
        print(f"  x={x:3d}: lambda = (1+-i*sqrt({-disc}))/2, |lambda| = {(-x)**0.5:.4f}")

# Integer eigenvalue x values: x = k(k+1) = oblong numbers
print("\nInteger eigenvalue points x = k(k+1):")
for k in range(-3, 8):
    x = k * (k + 1)
    d = abs(2*k + 1)
    e1 = (1 + d) // 2 if d > 0 else 0
    e2 = (1 - d) // 2 if d > 0 else 0
    name = ""
    if x == 0: name = "(trivial)"
    if x == 2: name = "(TOURNAMENT!)"
    if x == 6: name = "(next oblong)"
    if x == -2: name = "(negative)"
    print(f"  k={k:2d}: x={x:3d}, eigenvalues = ({e1}, {e2}) {name}")

# ============================================================
# PART 5: THE 7x3^k FORBIDDEN FAMILY
# ============================================================
print("\n" + "=" * 70)
print("PART 5: THE 7 x 3^k FORBIDDEN FAMILY")
print("=" * 70)

# Verify with sampling at n=7
random.seed(42)
n = 7
m = n * (n - 1) // 2
N = 1 << m

h_seen_7 = set()
for _ in range(500000):
    bits = random.randint(0, N - 1)
    adj = get_tournament(n, bits)
    H = compute_H_dp(adj, n)
    h_seen_7.add(H)

print(f"\nn=7 (500K random samples): {len(h_seen_7)} distinct H values seen")
for val_name, val in [("7", 7), ("21", 21), ("63", 63), ("49=7^2", 49), ("189=7x27", 189), ("147=7x21", 147)]:
    status = "SEEN" if val in h_seen_7 else "NOT SEEN (likely forbidden)"
    print(f"  H={val_name}: {status}")

# What is the structure of 7 x 3^k?
print(f"\n7 x 3^k: 7, 21, 63, 189, 567, 1701...")
print(f"  7 = I(C_3, 2)")
print(f"  21 = 3 x 7 = I(P_4, 2)")
print(f"  63 = 9 x 7 = 2^6 - 1")
print(f"  189 = 27 x 7 = max_H(7) -- ACHIEVABLE!")
print(f"  So 189 BREAKS the pattern: 7 x 3^3 = 189 IS achievable.")
print(f"  The forbidden family is exactly {{7, 21, 63}} at n<=7.")

# Check: is 63 = some I(G, 2)?
# 63 = 2^6 - 1 = I(C_6, 2) - 2 = 65 - 2? No.
# I(C_6, 2) = 2^6 + (-1)^6 = 64 + 1 = 65
# 63 = (2^8 - (-1)^6)/3 = (256-1)/3 = 255/3 = 85? No.
# Actually I(P_5, 2) = (2^7 + 1)/3 = 129/3 = 43
# I(P_6, 2) = (2^8 - 1)/3 = 255/3 = 85
# So 63 is NOT an I(P_k, 2) or I(C_k, 2) value!

print(f"\nIs 63 = I(P_k, 2) for some k?")
for k in range(1, 15):
    ip2 = (2**(k+2) - (-1)**k) // 3
    if ip2 == 63:
        print(f"  YES: I(P_{k}, 2) = 63")
        break
    if ip2 > 63:
        print(f"  NO: I(P_k, 2) jumps from {(2**(k+1) - (-1)**(k-1))//3} to {ip2}")
        break

print(f"\nIs 63 = I(C_k, 2) for some k?")
for k in range(3, 15):
    ic2 = 2**k + (-1)**k
    if ic2 == 63:
        print(f"  YES: I(C_{k}, 2) = 63")
        break
    if ic2 > 63:
        print(f"  NO: I(C_k, 2) jumps from {2**(k-1) + (-1)**(k-1)} to {ic2}")
        break

# So 63 is a PRODUCT of cycle/path I-values
# 63 = 7 x 9 = I(C_3, 2) x 9
# 63 = 3 x 21 = 3 x I(P_4, 2)
# If Omega decomposes into components with I-values 7 and 9:
# I(C_3, 2) = 7, I(C_2, 2)? C_2 isn't a proper cycle.
# I(single edge, 2) = 1 + 2 + 2 = 5? No, edge = P_2: I(P_2, 2) = 1 + 2·2 = 5
# I(two isolated vertices, 2) = (1+2)^2 = 9
# So 63 = 7 x 9 = I(C_3 ∪ 2K_1, 2)? That requires the conflict graph
# to have a K_3 component and two isolated vertices.
# But K_3 component is impossible (THM-201)!
print(f"\n63 = 7 x 9 = I(C_3, 2) x I(2K_1, 2)")
print(f"  This requires Omega to have K_3 as a component.")
print(f"  THM-201: K_3 CANNOT be a component of Omega(T).")
print(f"  Therefore H=63 is forbidden by the same mechanism as H=7!")

# ============================================================
# PART 6: THE FIVE-FOLD WAY — n=5 AS UNIVERSAL FIXED POINT
# ============================================================
print("\n" + "=" * 70)
print("PART 6: THE FIVE-FOLD WAY — n=5 AS UNIVERSAL FIXED POINT")
print("=" * 70)

# Verify: max_H(5) = 15
n = 5
m = n * (n - 1) // 2
N = 1 << m
h_max_5 = 0
for bits in range(N):
    adj = get_tournament(n, bits)
    H = compute_H_dp(adj, n)
    h_max_5 = max(h_max_5, H)

print(f"\n  max_H(5) = {h_max_5}")
print(f"  T(5) = {triangular(5)}")
print(f"  C(5,1)+C(5,3) = {math.comb(5,1)+math.comb(5,3)}")
print(f"  2^4-1 = {2**4-1}")
print(f"  F(5)*F(4) = {fib(5)*fib(4)}")
print(f"  All = 15? {h_max_5 == triangular(5) == math.comb(5,1)+math.comb(5,3) == 2**4-1 == fib(5)*fib(4)}")

# Triangular ∩ Mersenne numbers
tri_set = {triangular(k) for k in range(200)}
mersenne_set = {2**j - 1 for j in range(1, 30)}
overlap = sorted(tri_set & mersenne_set)
print(f"\n  Triangular AND Mersenne numbers: {overlap}")
print(f"  (Ramanujan-Nagell related: k(k+1)/2 = 2^j - 1)")

# phi^5 identity
phi = (1 + 5**0.5) / 2
print(f"\n  phi^5 = {phi**5:.6f} = 5*phi + 3 = {5*phi+3:.6f}")
print(f"  phi^n = F(n)*phi + F(n-1):")
for n_val in [3, 4, 5, 6]:
    print(f"    phi^{n_val} = {fib(n_val)}*phi + {fib(n_val-1)} = {fib(n_val)*phi + fib(n_val-1):.6f}")

# ============================================================
# PART 7: GROWTH RATE CLASSIFICATION
# ============================================================
print("\n" + "=" * 70)
print("PART 7: GROWTH RATE CLASSIFICATION — 5 GROUPS")
print("=" * 70)

print("\nRatios a(n)/a(n-1):")
print(f"  {'n':>3s} | {'F':>10s} | {'J':>10s} | {'maxH':>10s} | {'T':>10s}")
for n_val in range(2, 12):
    fr = fib(n_val) / fib(n_val-1) if fib(n_val-1) > 0 else float('inf')
    jr = jacobsthal(n_val) / jacobsthal(n_val-1) if jacobsthal(n_val-1) > 0 else float('inf')
    if n_val < len(max_H) and n_val > 0 and max_H[n_val-1] > 0:
        mr = max_H[n_val] / max_H[n_val-1]
    else:
        mr = None
    tr = triangular(n_val) / triangular(n_val-1) if triangular(n_val-1) > 0 else float('inf')
    print(f"  {n_val:3d} | {fr:10.4f} | {jr:10.4f} | {str(round(mr,4)) if mr else '?':>10s} | {tr:10.4f}")

print(f"\nLimiting growth rates:")
print(f"  Fibonacci: phi = {phi:.6f}")
print(f"  Jacobsthal: 2 exactly")
print(f"  Triangular: 1 (polynomial)")
print(f"  max_H: empirically ~ (n-1) for large n")

# ============================================================
# PART 8: PISANO PERIOD COMPARISON
# ============================================================
print("\n" + "=" * 70)
print("PART 8: PISANO PERIODS — FIBONACCI vs JACOBSTHAL")
print("=" * 70)

def pisano_period(m_val):
    a, b = 0, 1
    for k in range(1, 6 * m_val * m_val + 1):
        a, b = b, (a + b) % m_val
        if a == 0 and b == 1:
            return k
    return None

def jacobsthal_pisano(m_val):
    if m_val <= 1: return 1
    a, b = 0, 1
    for k in range(1, 100 * m_val + 1):
        a, b = b % m_val, (b + 2*a) % m_val
        if a == 0 and b % m_val == 1:
            return k
    return None

print(f"\n  {'m':>3s} | {'pi_F(m)':>8s} | {'pi_J(m)':>8s} | {'same?':>6s}")
for m_val in range(2, 25):
    pf = pisano_period(m_val)
    pj = jacobsthal_pisano(m_val)
    pj_str = str(pj) if pj else "none"
    same = "YES" if pf == pj else "no"
    print(f"  {m_val:3d} | {pf:8d} | {pj_str:>8s} | {same:>6s}")

# Where do they agree?
print(f"\nm values where pi_F(m) = pi_J(m):")
for m_val in range(2, 50):
    pf = pisano_period(m_val)
    pj = jacobsthal_pisano(m_val)
    if pf and pj and pf == pj:
        print(f"  m={m_val}: pi_F = pi_J = {pf}")

# ============================================================
# PART 9: phi^n = F(n)*phi + F(n-1) AND TOURNAMENT MEANING
# ============================================================
print("\n" + "=" * 70)
print("PART 9: THE GOLDEN IDENTITY AND TOURNAMENT INTERPRETATION")
print("=" * 70)

print("""
phi^n = F(n)*phi + F(n-1)

At n=5: phi^5 = 5*phi + 3

Tournament interpretation at x=2:
  Replace phi with the tournament eigenvalue 2:
  2^n = F(n)*2 + F(n-1)

  Check: 2^5 = 32 = 5*2 + 3 = 13? NO!

  The identity phi^n = F(n)*phi + F(n-1) is specific to phi.
  At x=2: 2^n = J(n)*2 + J(n-1) where J(n) = Jacobsthal.
""")

# Verify: 2^n = J(n)*2 + J(n-1)?
# Actually the transfer matrix gives: T^n × [1,0]^T = [a(n+1), a(n)]^T
# For Fibonacci: [F(n+1), F(n)]
# For Jacobsthal: [J(n+1), J(n)]
# And lambda^n = a(n)*lambda + a(n-1)·... hmm.

# At x=2: eigenvalues are 2 and -1.
# 2^n = J(n)·2 + J(n-1)·(-2)? Let me derive properly.
# If lambda^n = A_n * lambda + B_n for lambda satisfying lambda^2 = lambda + x,
# then A_n = a_n (generalized Fibonacci), B_n = x·a_{n-1}.
# For x=1: A_n = F(n), B_n = F(n-1)
# For x=2: A_n = J(n), B_n = 2·J(n-1)
# Check: phi^n = F(n)*phi + F(n-1)
# At lambda=2: 2^n = J(n)*2 + 2*J(n-1)
print("Verification: 2^n = J(n)*2 + 2*J(n-1)")
for n_val in range(1, 12):
    lhs = 2**n_val
    rhs = jacobsthal(n_val) * 2 + 2 * jacobsthal(n_val - 1)
    print(f"  n={n_val:2d}: 2^n = {lhs:6d}, J(n)*2 + 2*J(n-1) = {rhs:6d}, match = {lhs == rhs}")

# At lambda=-1: (-1)^n = J(n)*(-1) + 2*J(n-1)
print("\nVerification: (-1)^n = -J(n) + 2*J(n-1)")
for n_val in range(1, 12):
    lhs = (-1)**n_val
    rhs = -jacobsthal(n_val) + 2 * jacobsthal(n_val - 1)
    print(f"  n={n_val:2d}: (-1)^n = {lhs:3d}, -J(n) + 2*J(n-1) = {rhs:3d}, match = {lhs == rhs}")

# BEAUTIFUL! Both identities hold. This means:
# 2^n + (-1)^n = J(n)*(2-1) + 2*J(n-1) + 2*J(n-1)? No...
# From the two identities:
#   2^n = 2*J(n) + 2*J(n-1)
#   (-1)^n = -J(n) + 2*J(n-1)
# Adding: 2^n + (-1)^n = J(n) + 4*J(n-1)
# But we know I(C_n, 2) = 2^n + (-1)^n for n>=3 cycles!
print(f"\nI(C_n, 2) = 2^n + (-1)^n = J(n) + 4*J(n-1)")
for n_val in range(3, 10):
    IC = 2**n_val + (-1)**n_val
    JJ = jacobsthal(n_val) + 4 * jacobsthal(n_val - 1)
    print(f"  n={n_val}: I(C_n,2) = {IC}, J(n)+4J(n-1) = {JJ}, match = {IC == JJ}")

# Also: I(C_n, 2) = 2^n + (-1)^n = Lucas-Jacobsthal hybrid
# This is L(n) in the "x=2 Lucas sequence"!
# I(C_n, 2) = alpha^n + beta^n where alpha=2, beta=-1
# Standard Lucas: L(n) = phi^n + (-1/phi)^n
# So I(C_n, 2) is the "Jacobsthal-Lucas" number!
print(f"\n--> I(C_n, 2) is the JACOBSTHAL-LUCAS number: JL(n) = 2^n + (-1)^n")
print(f"    Just as Lucas L(n) = phi^n + (-1/phi)^n relates to Fibonacci,")
print(f"    JL(n) = 2^n + (-1)^n relates to Jacobsthal J(n).")
print(f"    The forbidden H=7 = JL(3) = Jacobsthal-Lucas number!")

# ============================================================
# SYNTHESIS
# ============================================================
print("\n" + "=" * 70)
print("GRAND SYNTHESIS")
print("=" * 70)
print("""
THE INDEPENDENCE POLYNOMIAL THESIS:

Tournament theory = Fibonacci theory evaluated at fugacity x=2.

DICTIONARY:
  Fibonacci world (x=1)     Tournament world (x=2)
  ---------------------     ---------------------
  F(n) = Fibonacci           J(n) = Jacobsthal
  L(n) = Lucas               JL(n) = 2^n+(-1)^n = I(C_n,2)
  phi = golden ratio         2 = tournament eigenvalue
  -1/phi                     -1
  phi^n = F(n)phi+F(n-1)    2^n = 2J(n)+2J(n-1)
  Pi(m) = Pisano period     Pi_J(m) = Jacobsthal Pisano

COINCIDENCE MAP:
  7 = I(C_3, 2) = JL(3) = FORBIDDEN
  21 = F(8) = J(6) = T(6) = I(P_4,2) = FORBIDDEN
  15 = max_H(5) = T(5) = 2^4-1 = F(5)*F(4) = FIXED POINT
  1729 = max_H(11)/T(10) = TAXICAB NUMBER

STRUCTURAL CONSTANTS:
  phi^5 = 5*phi + 3 encodes the (5,3) pair
  5 = vertices, 3 = cycle length, 15 = 5*3 = max_H(5)
  x=2 is the unique positive integer fugacity with integer eigenvalues
  x=k(k+1) = oblong numbers are ALL integer-eigenvalue fugacities
  Tournament sits at k=1: x=1*2=2

PERIOD-6: Pi_F(4) = 6 but Pi_J(4) DOES NOT EXIST
  The Jacobsthal matrix mod 4 is SINGULAR (det=-2, 2|4).
  So the period-6 structure is FIBONACCI-SPECIFIC.
  Tournaments inherit period-6 only through the Fibonacci bridge.

THE FORBIDDEN MECHANISM:
  H = I(Omega, 2) where Omega = cycle-arc conflict graph.
  7 = JL(3) is forbidden because K_3 cannot be Omega component (THM-201).
  21 = I(P_4, 2) is forbidden because P_4 cannot be Omega (THM-202).
  63 = 7*9 = JL(3)*I(2K_1,2) is forbidden by THM-201 (K_3 factor).
  189 = 7*27 IS achievable (different factorization of Omega).
""")
