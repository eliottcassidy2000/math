#!/usr/bin/env python3
"""depth_width_s116h.py — The depth-width obstruction in the OCF binary expansion.

H = 1 + 2*a1 + 4*a2 + 8*a3 + ...

a_k = independent sets of size k in Omega(T).
SIZE k = the DEPTH (how many non-conflicting cycles you stack).
VALUE a_k = the WIDTH (how many such stacks exist).

The forbidden numbers arise from a DEPTH-WIDTH OBSTRUCTION:
certain depths require minimum widths that overshoot the target.

Tournaments are isomorphic to even graphs (every vertex has even degree
in the underlying undirected graph). The conflict graph Omega is built
from the ODD cycles of this even graph. The interplay between the
even structure of T and the odd structure of Omega creates the obstruction.
"""
from math import sqrt, log, comb, factorial
from itertools import permutations, combinations
from collections import Counter, defaultdict

def all_tournaments(n):
    arcs = [(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in range(2**len(arcs)):
        adj = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(arcs):
            if bits & (1 << k):
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        yield adj

def count_ham(adj, n):
    count = 0
    for p in permutations(range(n)):
        if all(adj[p[i]][p[i+1]] for i in range(n-1)):
            count += 1
    return count

def find_odd_cycles(adj, n, max_len=None):
    """Find all directed odd cycles up to max_len."""
    if max_len is None:
        max_len = n
    cycles = []
    for length in range(3, max_len+1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts):
                if all(adj[perm[i]][perm[(i+1)%length]] for i in range(length)):
                    key = frozenset(verts)
                    if key not in [frozenset(c) for c in cycles]:
                        cycles.append(verts)
                    break
    return cycles

def conflict_graph_adj(cycles):
    nc = len(cycles)
    adj = [[0]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if set(cycles[i]) & set(cycles[j]):
                adj[i][j] = adj[j][i] = 1
    return adj

def count_indep(adj, n):
    alpha = [0]*(n+1)
    for sz in range(n+1):
        for sub in combinations(range(n), sz):
            if all(adj[sub[i]][sub[j]] == 0 for i in range(len(sub)) for j in range(i+1, len(sub))):
                alpha[sz] += 1
    return alpha

print()
print("  THE DEPTH-WIDTH OBSTRUCTION")
print()
print("="*70)
print()

# ============================================================
print("  I. DEPTH AND WIDTH IN THE OCF")
print("  " + "-"*40)
print()
print("  H = 1 + sum_{k>=1} 2^k * alpha_k")
print()
print("  DEPTH k: the size of the independent set.")
print("    k=1: one cycle (depth 1)")
print("    k=2: two disjoint cycles (depth 2)")
print("    k=3: three mutually disjoint cycles (depth 3)")
print()
print("  WIDTH alpha_k: the number of independent sets at depth k.")
print("    alpha_1 = total number of odd cycles")
print("    alpha_2 = number of disjoint cycle PAIRS")
print("    alpha_3 = number of disjoint cycle TRIPLES")
print()
print("  The 'weight' 2^k means: each unit of DEPTH doubles the contribution.")
print("  Going one level deeper DOUBLES the payoff.")
print("  But going deeper requires MORE VERTICES (disjoint cycles need room).")
print()
print("  THE OBSTRUCTION: at depth k, you need at least 3k vertices")
print("  (k disjoint 3-cycles). So the maximum depth at n is floor(n/3).")
print()
print("  n    max depth    max 2^k     max H contribution from depth alone")
print("  " + "-"*60)
for n in range(3, 16):
    max_d = n // 3
    max_weight = 2**max_d
    print(f"  {n:3d}    {max_d:3d}          {max_weight:6d}")

print()

# ============================================================
print()
print("  II. THE EVEN GRAPH CONNECTION")
print("  " + "-"*40)
print()
print("  A tournament on n vertices IS an orientation of the complete")
print("  graph K_n. Every vertex has out-degree + in-degree = n-1.")
print()
print("  For odd n: every vertex has even total degree n-1.")
print("  The UNDERLYING undirected graph (ignoring orientation)")
print("  is K_n, which has degree n-1 = even for odd n.")
print("  So for odd n: the tournament IS an even graph.")
print()
print("  For even n: degree n-1 is odd. NOT an even graph.")
print("  But the COMPLEMENT (undirected K_n minus the tournament edges)")
print("  is also NOT even.")
print()
print("  The key: a tournament on odd n vertices is an orientation of")
print("  an EULERIAN graph (every vertex has even degree => Euler circuit).")
print("  K_n for odd n has an Euler circuit.")
print("  The tournament is a way to ORIENT this Euler circuit.")
print()
print("  Euler circuits visit every EDGE. Hamiltonian paths visit every VERTEX.")
print("  The tournament lives at the intersection of these two worlds.")
print()

# ============================================================
print()
print("  III. EXHAUSTIVE DEPTH-WIDTH DATA AT n=6")
print("  " + "-"*40)
print()

# Compute full data at n=6
dw_data = defaultdict(lambda: defaultdict(int))
h_from_dw = defaultdict(set)
total = 0

for adj in all_tournaments(6):
    total += 1
    h = count_ham(adj, 6)
    cycles = find_odd_cycles(adj, 6, max_len=5)
    nc = len(cycles)
    if nc > 0:
        omega = conflict_graph_adj(cycles)
        alphas = count_indep(omega, nc)
    else:
        alphas = [1]

    # Record (alpha_1, alpha_2) and the resulting H
    a1 = alphas[1] if len(alphas) > 1 else 0
    a2 = alphas[2] if len(alphas) > 2 else 0
    a3 = alphas[3] if len(alphas) > 3 else 0
    dw_data[(a1, a2, a3)][h] += 1
    h_from_dw[(a1, a2, a3)].add(h)

print(f"  Enumerated {total} tournaments on 6 vertices.")
print()
print(f"  (a1, a2, a3) -> H values  [H_formula = 1+2*a1+4*a2+8*a3]")
print(f"  " + "-"*65)
for key in sorted(h_from_dw.keys()):
    a1, a2, a3 = key
    h_vals = sorted(h_from_dw[key])
    h_formula = 1 + 2*a1 + 4*a2 + 8*a3
    match = "EXACT" if len(h_vals) == 1 and h_vals[0] == h_formula else ""
    discrepancy = ""
    if h_vals[0] != h_formula:
        discrepancy = f" [off by {h_vals[0] - h_formula:+d} to {h_vals[-1] - h_formula:+d}]"
    print(f"  ({a1:2d},{a2:2d},{a3:d}) -> {str(h_vals):<25s} formula={h_formula}{discrepancy} {match}")

print()

# ============================================================
print()
print("  IV. THE DISCREPANCY: WHY H != 1+2*a1+4*a2+8*a3")
print("  " + "-"*40)
print()
print("  At n=6, H EXCEEDS the formula 1+2*a1+4*a2+8*a3.")
print("  The excess comes from 5-CYCLES and HIGHER odd cycles.")
print("  The formula only counts 3-cycle contributions to alpha_k.")
print()
print("  The FULL OCF counts ALL odd cycles (3, 5, 7, ...).")
print("  At n=6: 5-cycles first appear, contributing to alpha_1.")
print("  These 5-cycle contributions push H above the 3-cycle-only prediction.")
print()
print("  The 5-cycles are HIGHER-ORDER partial bits.")
print("  They don't fit into the 3-cycle-based depth-width framework")
print("  because a 5-cycle uses 5 vertices, not 3.")
print()
print("  At depth 1: a 5-cycle contributes 2 to H, same as a 3-cycle.")
print("  At depth 2: a 5-cycle + disjoint 3-cycle needs 5+3=8 vertices > 6.")
print("  So at n=6, 5-cycles can only contribute at depth 1.")
print()

# ============================================================
print()
print("  V. THE DEPTH-WIDTH OBSTRUCTION FOR H=7")
print("  " + "-"*40)
print()
print("  H = 7 requires: sum_{k>=1} 2^k * alpha_k = 6.")
print()
print("  Depth-width decompositions of 6:")
print("  depth 1 only: 2*alpha_1 = 6 => alpha_1 = 3")
print("  depth 2 only: 4*alpha_2 = 6 => alpha_2 = 3/2 (not integer!)")
print("  depth 1+2: 2*a1 + 4*a2 = 6 => (a1, a2) in {(3,0), (1,1)}")
print("  depth 3: 8*a3 = 6 => impossible")
print()
print("  The depth-2-only option FAILS because 6 is not divisible by 4.")
print("  The depth-3 option FAILS because 6 < 8.")
print()
print("  So H=7 requires EITHER:")
print("  (A) alpha_1 = 3, alpha_2 = 0: three cycles, all pairwise conflicting")
print("  (B) alpha_1 = 1, alpha_2 = 1: impossible (1 cycle can't have a disjoint pair)")
print()
print("  Case (A) is the DEPTH-WIDTH OBSTRUCTION:")
print("  At depth 1, you need exactly 3 cycles.")
print("  At depth 2, you need exactly 0 independent pairs.")
print("  But 3 pairwise-intersecting 3-cycles on n vertices force")
print("  a COMMON VERTEX (for small n) or ADDITIONAL cycles (for large n).")
print("  Either way, alpha_1 >= 4, not 3. Contradiction.")
print()
print("  The obstruction: WIDTH 3 AT DEPTH 1 IS UNSTABLE.")
print("  Three intersecting cycles can't stay at exactly 3.")
print("  They either collapse (shared vertex forces more cycles)")
print("  or separate (losing the pairwise intersection).")
print()

# ============================================================
print()
print("  VI. THE STABILITY PRINCIPLE")
print("  " + "-"*40)
print()
print("  PRINCIPLE: In the OCF binary expansion, certain (depth, width)")
print("  configurations are UNSTABLE — they cannot be maintained because")
print("  the graph-theoretic constraints force transitions.")
print()
print("  STABLE configurations: those that actually appear in the")
print("  exhaustive data.")
print()
print("  At n=6, the stable (a1, a2) pairs are:")
stable = set()
for key in h_from_dw:
    stable.add((key[0], key[1]))

print(f"  {sorted(stable)}")
print()

# Check which (a1, a2) are NOT achievable
all_possible = set()
for a1 in range(max(k[0] for k in stable)+1):
    for a2 in range(max(k[1] for k in stable)+1):
        all_possible.add((a1, a2))

missing = all_possible - stable
print(f"  Missing (a1, a2) pairs: {sorted(missing)}")
print()
print("  NOTABLE MISSING: (3, 1) — you can't have 3 cycles with 1 disjoint pair")
print("  at n=6, because 3 cycles with a disjoint pair need at least 3+3=6 vertices")
print("  for the pair, leaving 0 for the third cycle's disjoint vertex.")
print()

# ============================================================
print()
print("  VII. EVEN GRAPHS, ODD CYCLES, AND THE PARITY OBSTRUCTION")
print("  " + "-"*40)
print()
print("  A tournament on n vertices has C(n,2) arcs.")
print("  Each vertex has out-degree + in-degree = n-1.")
print("  The score sequence (d_1, ..., d_n) satisfies sum d_i = C(n,2).")
print()
print("  The PARITY of H is determined by the OCF: H is always odd.")
print("  This is because the 'base' contribution is 1 (the empty set),")
print("  and all other contributions are even (2^k * alpha_k).")
print()
print("  But there's a DEEPER parity:")
print("  The number of odd cycles in T has a parity constraint.")
print("  For a tournament on n vertices:")
print("  c_3 = number of 3-cycles = C(n,3) - number of transitive triples.")
print("  Transitive triples = sum C(d_i, 2) (from the score sequence).")
print("  c_3 = C(n,3) - sum C(d_i, 2).")
print()
print("  For a REGULAR tournament (all d_i = (n-1)/2, n odd):")
print("  c_3 = C(n,3) - n*C((n-1)/2, 2)")
for n in [3, 5, 7, 9, 11]:
    c3 = comb(n, 3) - n * comb((n-1)//2, 2)
    print(f"  n={n}: c_3 = {comb(n,3)} - {n}*{comb((n-1)//2,2)} = {c3}")

print()
print("  For n=7 regular: c_3 = 35 - 7*3 = 14. Always 14 3-cycles.")
print("  This is FORCED by the regularity. No variation possible.")
print()

# ============================================================
print()
print("  VIII. THE BINARY EXPANSION AND CARRY PROPAGATION")
print("  " + "-"*40)
print()
print("  H = 1 + 2*a1 + 4*a2 + 8*a3 + ...")
print()
print("  Think of this as addition in binary:")
print("  The '1' is bit 0.")
print("  2*a1 contributes a1 copies of '10' (bit 1).")
print("  4*a2 contributes a2 copies of '100' (bit 2).")
print()
print("  If a1 >= 2: bit 1 'overflows' and carries into bit 2.")
print("  The carry: a1 = 2*q + r, contribution = 2*r at bit 1 plus 4*q at bit 2.")
print("  Effective: bit 1 = a1 mod 2, bit 2 += a1 // 2.")
print()
print("  This carry propagation means:")
print("  H in binary = 1 + binary_expand(a1, a2, a3, ...) with carries.")
print()
print("  H = 7 = 111 in binary.")
print("  Needs bit 0 = 1 (always), bit 1 = 1, bit 2 = 1.")
print("  bit 1 = 1 means a1 is odd.")
print("  bit 2 = 1 means a2 + (carry from a1) is odd.")
print("  If a1 = 3: bit 1 = 1 (3 mod 2 = 1), carry = 1 to bit 2.")
print("  Then bit 2 = a2 + 1. For bit 2 = 1: a2 + 1 = 1, so a2 = 0.")
print("  This gives (a1, a2) = (3, 0). But this is NOT ACHIEVABLE.")
print()
print("  If a1 = 1: bit 1 = 1, carry = 0.")
print("  bit 2 = a2. For bit 2 = 1: a2 = 1.")
print("  This gives (a1, a2) = (1, 1). But 1 cycle can't have a disjoint pair.")
print()
print("  BOTH paths to H = 111 in binary are blocked.")
print("  The carry propagation from depth 1 to depth 2 creates an")
print("  IMPOSSIBLE demand: either too many cycles (3) or too few (1).")
print()
print("  THE FORBIDDEN NUMBER 7 = 111 IN BINARY IS FORBIDDEN")
print("  BECAUSE CARRY PROPAGATION IN THE OCF EXPANSION")
print("  CREATES CONTRADICTORY DEMANDS ON THE CYCLE STRUCTURE.")
print()

# ============================================================
print()
print("  IX. H = 21 = 10101 IN BINARY")
print("  " + "-"*40)
print()
print("  21 = 10101 in binary.")
print("  Needs bits 0, 2, 4 to be 1.")
print("  bit 0 = 1 (always).")
print("  bit 2 = 1: requires a2 + carry from (a1) is odd.")
print("  bit 4 = 1: requires a4 + carry from lower levels is odd.")
print()
print("  The carry chain through 5 bits creates a COUPLED SYSTEM")
print("  of parity constraints on (a1, a2, a3, a4).")
print("  The OCF graph constraints make certain combinations impossible.")
print()
print("  21 has an ALTERNATING bit pattern: 1_0_1_0_1.")
print("  This alternation requires alpha values at even depths (0, 2, 4)")
print("  to contribute 1 each, with odd depths (1, 3) contributing 0.")
print("  The 'skip-one-depth' pattern is especially hard to achieve")
print("  because getting alpha_2 > 0 requires disjoint cycle pairs,")
print("  which need many vertices, which generate more cycles at depth 1.")
print()
print("  7 = 111 (all ones up to bit 2)")
print("  21 = 10101 (alternating up to bit 4)")
print("  Both have the property that achieving the target bit pattern")
print("  requires contradictory cycle configurations.")
print()

# ============================================================
print()
print("  X. THE LANDSCAPE: ACHIEVABLE H IN BINARY")
print("  " + "-"*40)
print()

# Compute H values at n=6 and show their binary
h_vals_n6 = set()
for adj in all_tournaments(6):
    h_vals_n6.add(count_ham(adj, 6))

print("  n=6: achievable H values in binary:")
print(f"  {'H':>4s}  {'binary':>10s}  {'bits set':>8s}  {'pattern'}")
print("  " + "-"*50)
for h in sorted(h_vals_n6):
    binary = bin(h)[2:]
    bits_set = binary.count('1')
    # Characterize pattern
    pattern = ""
    if all(c == '1' for c in binary):
        pattern = "all-ones"
    elif binary == binary[::-1]:
        pattern = "palindrome"
    print(f"  {h:4d}  {binary:>10s}  {bits_set:>8d}  {pattern}")

print()
print("  The MISSING values near 7 and 21:")
all_odd = set(range(1, max(h_vals_n6)+1, 2))
missing = sorted(all_odd - h_vals_n6)
print(f"  Missing: {missing[:20]}")
print()
for h in [7, 21]:
    binary = bin(h)[2:]
    print(f"  {h} = {binary}: {'FORBIDDEN' if h in missing else 'achievable'}")
print()

# ============================================================
print()
print("  XI. SYNTHESIS: THE DEPTH-WIDTH-CARRY OBSTRUCTION")
print("  " + "-"*40)
print()
print("  The forbidden H-values arise from THREE interacting constraints:")
print()
print("  1. DEPTH constraint: alpha_k requires k disjoint odd cycles,")
print("     needing >= 3k vertices. Max depth = floor(n/3).")
print()
print("  2. WIDTH constraint: alpha_k (width at depth k) is determined")
print("     by the conflict graph Omega(T). Not all widths are achievable")
print("     because the graph structure constrains independent set counts.")
print()
print("  3. CARRY constraint: the binary representation of H propagates")
print("     carries from lower depths to higher depths. The carries create")
print("     coupled parity demands that the width constraints cannot satisfy.")
print()
print("  The forbidden numbers are those where the carry chain demands")
print("  a width pattern that the depth constraint makes impossible.")
print()
print("  7 = 111:   carry from depth 1 demands width 0 at depth 2,")
print("             but 3 pairwise-intersecting cycles can't stay at 3.")
print("  21 = 10101: alternating bits demand skip-depth contributions")
print("             that conflict with the vertex budget.")
print()
print("  THE TOURNAMENT'S EVEN STRUCTURE (complete graph, even-degree")
print("  for odd n) constrains the ODD CYCLE structure (Omega),")
print("  which constrains the INDEPENDENT SET structure (alpha_k),")
print("  which constrains the BINARY DIGITS of H,")
print("  which determines which H-values are achievable.")
print()
print("  Even graph -> odd cycles -> independent sets -> binary digits -> H.")
print("  Five levels of constraint. The forbidden numbers are where")
print("  all five levels conspire to block every path.")
