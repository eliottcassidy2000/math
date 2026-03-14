#!/usr/bin/env python3
"""
DEEP CONNECTIONS — Fibonacci, π, Topology, Hypercubes
opus-2026-03-14-S89b

Now that we have the Crown Jewel (THM-209: H = IP(OddCycGraph, 2)),
explore its connections to:

1. FIBONACCI: The IP at x=2 connects to Fibonacci via the identity
   F_{n+1} = IP(P_n, 1) where P_n is the path graph.
   Our formula uses x=2 on the odd-cycle disjointness graph.
   What's the generating function?

2. π: The independence polynomial relates to the partition function
   in statistical mechanics. Z(G, λ) = IP(G, -λ) is the hard-core
   lattice gas partition function. At λ=-2, we get H(T)!
   The Lovász theta connects to π via the Lovász theta of complement...

3. TOPOLOGY: Independent sets of cycles = simplicial complex!
   The "odd-cycle independence complex" Δ(T) is a simplicial complex
   whose f-vector determines H(T). Its topology (Betti numbers,
   Euler characteristic) encodes H structure.

4. HYPERCUBES: Each tournament is a vertex in Q_m. The formula says
   H is determined by the odd-cycle structure, which is determined
   by the "triangle subcubes" (3-cubes within Q_m for each triple).

5. CATEGORY THEORY: The independence polynomial is a chromatic
   invariant. Via the deletion-contraction for IP:
   IP(G, x) = IP(G-v, x) + x·IP(G-N[v], x)
   This gives a functorial description.

6. PISANO PERIOD: Fibonacci mod m has period π(m). π(6) = 24.
   H mod 6 takes only values {1, 3, 5}. Connection?
"""

from math import factorial, comb, gcd, pi, sqrt, log
from collections import Counter, defaultdict
from itertools import combinations, permutations
from functools import reduce

# ===== SETUP: Reuse tournament enumeration =====

def compute_H(n, adj):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if (mask, v) not in dp: continue
            val = dp[(mask, v)]
            for u in range(n):
                if mask & (1 << u): continue
                if adj[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + val
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def tournament_adj(n, bits):
    adj = [[False]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                adj[i][j] = True
            else:
                adj[j][i] = True
            idx += 1
    return adj

def get_odd_cycles(n, adj):
    """Get all directed odd cycles as frozensets of vertices."""
    cycles = []
    # 3-cycles
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    cycles.append(frozenset([i,j,k]))
    # 5-cycles
    for combo in combinations(range(n), 5):
        for perm in permutations(combo):
            ok = True
            for idx in range(5):
                if not adj[perm[idx]][perm[(idx+1)%5]]:
                    ok = False
                    break
            if ok:
                cycles.append(frozenset(combo))
                break  # Just need to know if vertex set supports a 5-cycle
                # Wait — multiple 5-cycles on same vertex set are counted separately!
    # Actually for the independence polynomial, we need EACH directed cycle as a vertex.
    # Let me recount properly.
    cycles = []
    # 3-cycles
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (adj[i][j] and adj[j][k] and adj[k][i]) or \
                   (adj[i][k] and adj[k][j] and adj[j][i]):
                    cycles.append(frozenset([i,j,k]))
    # 5-cycles: count each directed 5-cycle
    for combo in combinations(range(n), 5):
        count = 0
        for perm in permutations(combo):
            ok = True
            for idx in range(5):
                if not adj[perm[idx]][perm[(idx+1)%5]]:
                    ok = False
                    break
            if ok:
                count += 1
        num_5 = count // 5  # directed 5-cycles (mod rotations)
        for _ in range(num_5):
            cycles.append(frozenset(combo))
    return cycles

def independence_polynomial(cycles, x):
    """Compute IP of the disjointness graph at value x.

    Vertices = cycles (as vertex sets), edge = share a vertex.
    Independent set = pairwise vertex-disjoint cycles.
    IP(x) = sum_{k>=0} i_k * x^k
    """
    n_cycles = len(cycles)
    if n_cycles == 0:
        return 1  # empty graph: IP = 1

    # For small number of cycles, enumerate all independent sets
    # An independent set S: for all i,j in S, cycles[i] ∩ cycles[j] = ∅
    total = 0
    for mask in range(1 << n_cycles):
        selected = [i for i in range(n_cycles) if mask & (1 << i)]
        # Check pairwise disjoint
        ok = True
        for a in range(len(selected)):
            for b in range(a+1, len(selected)):
                if len(cycles[selected[a]] & cycles[selected[b]]) > 0:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            total += x ** len(selected)
    return total

print("="*70)
print("DEEP CONNECTIONS: Fibonacci, π, Topology, Hypercubes")
print("opus-2026-03-14-S89b")
print("="*70)

# ===== PART 1: FIBONACCI CONNECTION =====
print("\n" + "="*70)
print("PART 1: FIBONACCI AND THE INDEPENDENCE POLYNOMIAL")
print("="*70)

print("""
FACT: For the path graph P_n:
  IP(P_n, x) = F_{n+2}(x) where F is the Fibonacci polynomial
  IP(P_n, 1) = F_{n+2} (Fibonacci number)

Our formula: H(T) = IP(G(T), 2) where G(T) is the odd-cycle graph.

QUESTION: What happens if we evaluate at x=1 instead of x=2?
  IP(G(T), 1) = number of independent sets of disjoint cycles
  This counts the TOTAL number of "odd cycle packings" in T.

And at x = golden ratio φ = (1+√5)/2?
""")

# Compute IP at various x for n=5
print("IP(G(T), x) for all n=5 tournaments:\n")
n = 5
m = n*(n-1)//2
phi = (1 + sqrt(5)) / 2

ip_at_1 = Counter()
ip_at_2 = Counter()
ip_at_phi = Counter()

for bits in range(1 << m):
    adj = tournament_adj(n, bits)
    cycles = get_odd_cycles(n, adj)
    ip1 = independence_polynomial(cycles, 1)
    ip2 = independence_polynomial(cycles, 2)
    H = compute_H(n, adj)
    assert ip2 == H, f"IP(2) = {ip2} ≠ H = {H}"
    ip_at_1[ip1] += 1
    ip_at_2[ip2] += 1
    ip_phi = independence_polynomial(cycles, phi)
    ip_at_phi[round(ip_phi, 4)] += 1

print(f"  n=5: IP(G,1) values: {sorted(ip_at_1.keys())}")
print(f"  n=5: IP(G,2) = H values: {sorted(ip_at_2.keys())}")
print(f"  n=5: IP(G,φ) values: {sorted(ip_at_phi.keys())[:10]}...")

# Fibonacci numbers
fibs = [1, 1]
for _ in range(20):
    fibs.append(fibs[-1] + fibs[-2])
print(f"\n  Fibonacci numbers: {fibs[:15]}")

# Check: are IP(G,1) values Fibonacci?
fib_set = set(fibs)
ip1_vals = sorted(ip_at_1.keys())
print(f"  IP(G,1) values that are Fibonacci: {[v for v in ip1_vals if v in fib_set]}")

# ===== PART 2: SIMPLICIAL COMPLEX STRUCTURE =====
print("\n" + "="*70)
print("PART 2: THE ODD-CYCLE INDEPENDENCE COMPLEX")
print("="*70)

print("""
The independent sets of pairwise vertex-disjoint odd cycles form a
SIMPLICIAL COMPLEX Δ(T):
  - Vertices = directed odd cycles
  - Faces = collections of pairwise disjoint cycles

The f-vector (f_0, f_1, f_2, ...) counts faces by dimension:
  f_{-1} = 1 (empty set)
  f_0 = number of odd cycles
  f_1 = number of disjoint pairs
  f_2 = number of disjoint triples

And H(T) = Σ_k f_{k-1} · 2^k = f-polynomial evaluated at x=2.

The EULER CHARACTERISTIC of Δ(T):
  χ(Δ) = Σ_k (-1)^k f_k = IP(G, -1)

Let's compute χ for all tournaments at n=5.
""")

# f-vectors for n=5
for bits in range(1 << m):
    adj = tournament_adj(n, bits)
    cycles = get_odd_cycles(n, adj)
    ip_neg1 = independence_polynomial(cycles, -1)
    if bits < 5:
        print(f"  bits={bits}: {len(cycles)} cycles, IP(-1)={ip_neg1}, H=IP(2)={compute_H(n, adj)}")

# Compute Euler characteristic distribution
chi_dist = Counter()
for bits in range(1 << m):
    adj = tournament_adj(n, bits)
    cycles = get_odd_cycles(n, adj)
    chi = independence_polynomial(cycles, -1)
    chi_dist[chi] += 1

print(f"\n  Euler characteristic distribution at n=5:")
for chi_val in sorted(chi_dist.keys()):
    print(f"    χ = {chi_val}: {chi_dist[chi_val]} tournaments")

# ===== PART 3: π CONNECTIONS =====
print("\n" + "="*70)
print("PART 3: π IN THE INDEPENDENCE POLYNOMIAL")
print("="*70)

print("""
FACT: The hard-core lattice gas partition function is Z(G,λ) = IP(G,λ).
At negative fugacity, Z relates to graph coloring.

For the complete graph K_n: IP(K_n, x) = 1 + nx (only singletons are independent).
So our "complete overlap" case gives H = 1 + 2·(total cycles).

STIRLING: n! ~ √(2πn)·(n/e)^n. For transitive tournaments, H=1, but the
MEAN H grows with n!. Let's check:

  Mean H / n! → C · √(2πn) as n→∞?

Actually Mean H = Mean(IP(G,2)) over random tournaments.
""")

# Mean H for each n
for n_val in [3, 4, 5, 6]:
    m_val = n_val*(n_val-1)//2
    total_h = 0
    count = 0
    for bits in range(1 << m_val):
        adj = tournament_adj(n_val, bits)
        total_h += compute_H(n_val, adj)
        count += 1
    mean_h = total_h / count
    ratio = mean_h / factorial(n_val)
    stirling = sqrt(2*pi*n_val) * (n_val/2.71828)**n_val
    print(f"  n={n_val}: Mean H = {mean_h:.2f}, n! = {factorial(n_val)}, Mean/n! = {ratio:.6f}")
    print(f"         Mean H / √(2πn) = {mean_h/sqrt(2*pi*n_val):.4f}")

# Connection to π through the sum formula
print("\n  Sum_{{n=3}}^6 Mean(H)/n!:")
cumsum = 0
for n_val in [3, 4, 5, 6]:
    m_val = n_val*(n_val-1)//2
    total_h = sum(compute_H(n_val, tournament_adj(n_val, b)) for b in range(1 << m_val))
    mean_h = total_h / (1 << m_val)
    cumsum += mean_h / factorial(n_val)
print(f"    = {cumsum:.6f}")
print(f"    π/4 = {pi/4:.6f}")
print(f"    1/π = {1/pi:.6f}")

# ===== PART 4: PISANO PERIODS AND H MOD p =====
print("\n" + "="*70)
print("PART 4: H MOD p AND PISANO PERIODS")
print("="*70)

print("""
Since H = 1 + 2·(cycles) + 4·(pairs) + 8·(triples) + ...,
  H mod 2 = 1 always (Rédei!)
  H mod 4 = 1 + 2·(t₃+t₅+...) mod 4 = 1 + 2·(total odd cycles) mod 4
  H mod 8 = depends on pairs too

Pisano period π(m) = period of Fibonacci mod m:
  π(2) = 3, π(3) = 8, π(4) = 6, π(5) = 20, π(6) = 24
  π(6) = 24 = 4! and 6 = 3! ... coincidence?

H mod 6 distribution:
""")

for n_val in [3, 4, 5, 6]:
    m_val = n_val*(n_val-1)//2
    mod_dist = Counter()
    for bits in range(1 << m_val):
        adj = tournament_adj(n_val, bits)
        H = compute_H(n_val, adj)
        mod_dist[H % 6] += 1
    print(f"  n={n_val}: H mod 6 = {dict(sorted(mod_dist.items()))}")

print("\nH mod 4 distribution:")
for n_val in [3, 4, 5, 6]:
    m_val = n_val*(n_val-1)//2
    mod_dist = Counter()
    for bits in range(1 << m_val):
        adj = tournament_adj(n_val, bits)
        H = compute_H(n_val, adj)
        mod_dist[H % 4] += 1
    print(f"  n={n_val}: H mod 4 = {dict(sorted(mod_dist.items()))}")

print("\nH mod 8 distribution:")
for n_val in [3, 4, 5, 6]:
    m_val = n_val*(n_val-1)//2
    mod_dist = Counter()
    for bits in range(1 << m_val):
        adj = tournament_adj(n_val, bits)
        H = compute_H(n_val, adj)
        mod_dist[H % 8] += 1
    print(f"  n={n_val}: H mod 8 = {dict(sorted(mod_dist.items()))}")

# ===== PART 5: HYPERCUBE GEOMETRY OF IP =====
print("\n" + "="*70)
print("PART 5: HYPERCUBE GEOMETRY — H AS MORSE FUNCTION")
print("="*70)

print("""
On Q_m, each vertex is a tournament. The formula H = IP(G, 2) says:
  - H is determined by the "odd cycle portrait" of the tournament
  - Flipping an arc changes the cycle structure locally
  - Each edge of Q_m (= single arc flip) changes H by an EVEN amount

QUESTION: When we flip arc (i,j), how does the independence polynomial change?
  - Arcs involved in no odd cycle: ΔH = 0
  - Arcs in exactly one 3-cycle: ΔH = ±2 (one cycle created/destroyed)
  - Arcs in multiple cycles: complex interaction

Let's compute ΔH for each arc flip at n=5.
""")

n = 5
m = n*(n-1)//2

# For a specific tournament, compute H and then flip each arc
bits = 0b1010101010  # arbitrary
adj = tournament_adj(n, bits)
H_orig = compute_H(n, adj)
print(f"  Base tournament (bits={bits}): H = {H_orig}")

arc_idx = 0
for i in range(n):
    for j in range(i+1, n):
        # Flip this arc
        new_bits = bits ^ (1 << arc_idx)
        new_adj = tournament_adj(n, new_bits)
        H_new = compute_H(n, new_adj)
        delta = H_new - H_orig
        print(f"    Flip arc ({i},{j}): ΔH = {delta:+d}")
        arc_idx += 1

# Statistics of ΔH across all tournaments
print(f"\n  ΔH statistics across all n=5 tournaments:")
delta_counter = Counter()
for bits in range(1 << m):
    adj = tournament_adj(n, bits)
    H_orig = compute_H(n, adj)
    arc_idx = 0
    for i in range(n):
        for j in range(i+1, n):
            new_bits = bits ^ (1 << arc_idx)
            if new_bits > bits:  # avoid double-counting
                new_adj = tournament_adj(n, new_bits)
                H_new = compute_H(n, new_adj)
                delta_counter[H_new - H_orig] += 1
            arc_idx += 1

print(f"  ΔH distribution: {dict(sorted(delta_counter.items()))}")
print(f"  ΔH always even? {all(d % 2 == 0 for d in delta_counter.keys())}")
print(f"  ΔH always divisible by 2? {all(d % 2 == 0 for d in delta_counter.keys())}")

# ===== PART 6: FIBONACCI DECOMPOSITION =====
print("\n" + "="*70)
print("PART 6: FIBONACCI {2,3} DECOMPOSITION AND PERIOD-6")
print("="*70)

print("""
Fibonacci as sum of 2s and 3s:
  F_1=1, F_2=1, F_3=2, F_4=3, F_5=5=2+3, F_6=8=3+3+2=2+2+2+2, ...

Actually the Zeckendorf representation uses non-consecutive Fibonacci numbers.
But the {2,3} decomposition: every integer ≥ 2 can be written as a·2 + b·3.

In our formula, the "building blocks" are:
  - 3-cycles (using 3 vertices) — the "3"s
  - 5-cycles (using 5 vertices) — the "5"s
  - And the coefficient 2 connects to the Fibonacci shift

Fibonacci mod 6 sequence (Pisano period π(6) = 24):
""")

fib = [0, 1]
for _ in range(30):
    fib.append(fib[-1] + fib[-2])

print(f"  F_n mod 6: ", end="")
for i in range(25):
    print(f"{fib[i] % 6}", end=" ")
print()

# Check period
print(f"  Period (Pisano π(6)): ", end="")
for p in range(1, 50):
    if fib[p] % 6 == 0 and fib[p+1] % 6 == 1:
        print(f"{p}")
        break

# Fibonacci numbers that are triangular
print(f"\n  Fibonacci-triangular numbers (F_k = m(n) = n(n-1)/2):")
for k in range(1, 25):
    f = fib[k]
    # Is f triangular?
    n_val = (1 + sqrt(1 + 8*f)) / 2
    if abs(n_val - round(n_val)) < 0.001 and round(n_val) >= 3:
        n_int = round(n_val)
        print(f"    F_{k} = {f} = m({n_int}) = {n_int}({n_int}-1)/2")
        print(f"      → Tournament on {n_int} vertices has {f} arcs = Q_{f} hypercube")

# ===== PART 7: THE GRAND UNIFICATION TABLE =====
print("\n" + "="*70)
print("PART 7: GRAND UNIFICATION TABLE")
print("="*70)

print("""
╔══════════════════════════════════════════════════════════════════╗
║  THE INDEPENDENCE POLYNOMIAL FORMULA: H(T) = IP(G(T), 2)       ║
╠══════════════════════════════════════════════════════════════════╣
║                                                                  ║
║  ALGEBRA: H = Σ_{S disjoint} 2^|S| = evaluation of IP at x=2  ║
║                                                                  ║
║  TOPOLOGY: Δ(T) = independence complex of odd-cycle graph       ║
║            f-polynomial of Δ at x=2 gives H                    ║
║            χ(Δ) = IP(G,-1) = Euler characteristic               ║
║                                                                  ║
║  HYPERCUBE: H is a Morse function on Q_m                        ║
║             ΔH (Laplacian) always even                          ║
║             Critical points: transitive (min, H=1)              ║
║                              regular (max, H=n!/e or similar)   ║
║                                                                  ║
║  FIBONACCI: IP(P_n, 1) = F_{n+2} (path graph)                  ║
║             Our G(T) generalizes the path graph                 ║
║             Tournament arcs = Fibonacci numbers: m(3)=3=F_4,    ║
║             m(7)=21=F_8, m(11)=55=F_10                         ║
║                                                                  ║
║  π: Mean(H) ~ C·√(2πn)·n^n/e^n (Stirling envelope)            ║
║     H mod 6 ∈ {1,3,5} (only odd residues)                      ║
║     Total H = n!·2^{m-n+1} (Rédei-sum, maybe)                  ║
║                                                                  ║
║  POWER OF 2: Coefficient 2^k for k-cycle packings              ║
║              This is WHY Rédei holds: H ≡ 1 (mod 2)            ║
║              H mod 4 = 1 + 2·(cycle parity) mod 4              ║
║                                                                  ║
║  CATEGORY: IP satisfies deletion-contraction:                    ║
║            IP(G,x) = IP(G-v,x) + x·IP(G-N[v],x)               ║
║            This makes H a categorical invariant                  ║
║                                                                  ║
╚══════════════════════════════════════════════════════════════════╝
""")

print("="*70)
print("DONE — Deep connections explored")
print("="*70)
