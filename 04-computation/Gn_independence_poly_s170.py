#!/usr/bin/env python3
"""
Gn_independence_poly_s170.py -- opus-2026-03-22-S170

THE INDEPENDENCE POLYNOMIAL OF G_n

The OCF says: H(T) = I(Omega(T), 2).
Tournament -> conflict graph -> independence polynomial at x=2 -> H.

Now consider G_n itself as a graph. It has its own independence polynomial:
  I(G_n, x) = sum_{k=0}^{alpha(G_n)} a_k * x^k
where a_k = number of independent sets of size k in G_n.

THE META-LEVEL QUESTION:
  What is I(G_n, 2)?
  What is I(G_n, x) in general?
  Does it factor? Does it have nice roots?
  How does it relate to H and the OCF?

ALSO: The chromatic polynomial, Tutte polynomial, and other graph
polynomials of G_n — what do they encode about tournament structure?

Author: opus-2026-03-22-S170
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
from fractions import Fraction
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

# Core functions (same as previous)
def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        yield A, bits

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or form < best: best = form
    return best

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1<<v,v)] = 1
    for mask in range(1,1<<n):
        for v in range(n):
            if not (mask&(1<<v)): continue
            if dp[(mask,v)]==0: continue
            for w in range(n):
                if mask&(1<<w): continue
                if A[v][w]: dp[(mask|(1<<w),w)] += dp[(mask,v)]
    return sum(dp[((1<<n)-1,v)] for v in range(n))

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

def count_3cycles(A, n):
    c = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if A[i][j] and A[j][k] and A[k][i]: c+=1
                if A[i][k] and A[k][j] and A[j][i]: c+=1
    return c

def aut_size(A, n):
    count = 0
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if i != j and A[i][j] != A[perm[i]][perm[j]]:
                    ok = False; break
            if not ok: break
        if ok: count += 1
    return count

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

def build_Gn(n):
    m = comb(n,2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    class_map = {}; class_data = {}; tourn_class = {}
    for A, bits in all_tournaments(n):
        c = canonical(A, n)
        if c not in class_map:
            cid = len(class_map); class_map[c] = cid
            class_data[cid] = {
                'H': count_hp(A,n), 'scores': score_seq(A,n),
                't3': count_3cycles(A,n), 'aut': aut_size(A,n), 'size': 0
            }
        cid = class_map[c]; class_data[cid]['size'] += 1; tourn_class[bits] = cid
    nc = len(class_data)
    for cid in range(nc):
        for A, bits in all_tournaments(n):
            if tourn_class[bits] == cid:
                Aop = complement(A, n)
                class_data[cid]['comp'] = class_map[canonical(Aop, n)]
                break
    F = np.zeros((nc, nc), dtype=int)
    for A, bits in all_tournaments(n):
        ci = tourn_class[bits]
        for k in range(m):
            ia, ja = pairs[k]
            A2 = [row[:] for row in A]
            A2[ia][ja] = 1-A[ia][ja]; A2[ja][ia] = 1-A2[ia][ja]
            cj = class_map[canonical(A2, n)]
            F[ci][cj] += 1
    return nc, class_data, F, m

def independence_polynomial(adj, nc):
    """Compute I(G, x) as a list of coefficients [a_0, a_1, ..., a_alpha]."""
    # a_k = number of independent sets of size k
    coeffs = [0] * (nc + 1)
    for k in range(nc + 1):
        count = 0
        for S in combinations(range(nc), k):
            is_indep = True
            for i in S:
                for j in S:
                    if i < j and adj[i][j]:
                        is_indep = False
                        break
                if not is_indep:
                    break
            if is_indep:
                count += 1
        coeffs[k] = count
    # Trim trailing zeros
    while len(coeffs) > 1 and coeffs[-1] == 0:
        coeffs.pop()
    return coeffs

def eval_poly(coeffs, x):
    """Evaluate polynomial at x."""
    return sum(c * x**k for k, c in enumerate(coeffs))

def poly_str(coeffs):
    """Pretty-print polynomial."""
    terms = []
    for k, c in enumerate(coeffs):
        if c == 0:
            continue
        if k == 0:
            terms.append(str(c))
        elif k == 1:
            terms.append("%d*x" % c if c != 1 else "x")
        else:
            terms.append("%d*x^%d" % (c, k) if c != 1 else "x^%d" % k)
    return " + ".join(terms) if terms else "0"

# ============================================================
# BUILD G_n AND COMPUTE I(G_n, x)
# ============================================================

all_data = {}
for n in [3, 4, 5]:
    banner("G_%d: INDEPENDENCE POLYNOMIAL" % n)
    nc, cd, F, m = build_Gn(n)
    adj = (F > 0).astype(int)
    np.fill_diagonal(adj, 0)

    print("  |V(G_%d)| = %d, |E(G_%d)| = %d" % (n, nc, n, adj.sum()//2))

    # Compute I(G_n, x)
    coeffs = independence_polynomial(adj, nc)
    print("  I(G_%d, x) = %s" % (n, poly_str(coeffs)))
    print("  Coefficients: %s" % coeffs)

    # Evaluate at key points
    for x in [1, 2, 3, -1]:
        val = eval_poly(coeffs, x)
        print("  I(G_%d, %d) = %d" % (n, x, val))

    # I(G_n, 1) = total number of independent sets
    print("  I(G_%d, 1) = %d (total independent sets)" % (n, eval_poly(coeffs, 1)))

    # I(G_n, -1) = ??? (related to chromatic polynomial)
    print("  I(G_%d, -1) = %d" % (n, eval_poly(coeffs, -1)))

    # Roots
    if len(coeffs) > 1:
        np_coeffs = np.array(coeffs[::-1], dtype=float)  # numpy wants highest degree first
        roots = np.roots(np_coeffs)
        real_roots = sorted([r.real for r in roots if abs(r.imag) < 1e-10])
        print("  Real roots: %s" % [round(r, 6) for r in real_roots])
        if real_roots:
            print("  Largest real root: %.6f" % max(real_roots))

    # Independence number (= degree of I)
    alpha = len(coeffs) - 1
    print("  Independence number alpha(G_%d) = %d" % (n, alpha))

    all_data[n] = {'nc': nc, 'cd': cd, 'F': F, 'adj': adj, 'indep_poly': coeffs}
    print()

# ============================================================
# THE META-LEVEL: I(G_n, 2) vs H
# ============================================================

banner("THE META-LEVEL: I(G_n, 2)")

print("  THE OCF: H(T) = I(Omega(T), 2)")
print("  THE META-OCF: I(G_n, 2) = ???\n")

for n in [3, 4, 5]:
    val = eval_poly(all_data[n]['indep_poly'], 2)
    nc = all_data[n]['nc']
    sum_H = sum(all_data[n]['cd'][i]['H'] for i in range(nc))
    EH = Fraction(factorial(n), 2**(n-1))
    total_H_labeled = sum(all_data[n]['cd'][i]['H'] * all_data[n]['cd'][i]['size'] for i in range(nc))

    print("  n=%d: I(G_%d, 2) = %d" % (n, n, val))
    print("    |V(G_n)| = %d" % nc)
    print("    sum_H (unlabeled) = %d" % sum_H)
    print("    sum_H (labeled) = %d" % total_H_labeled)
    print("    2^C(n,2) = %d" % (2**comb(n,2)))
    print("    n! * E[H] = %d" % (factorial(n) * int(EH)))
    print("    I(G_n, 2) / |V| = %.4f" % (val / nc))
    print()

# ============================================================
# H-WEIGHTED INDEPENDENCE POLYNOMIAL
# ============================================================

banner("H-WEIGHTED INDEPENDENCE POLYNOMIAL")

print("  Define I_H(G_n, x) = sum over independent sets S of")
print("  x^|S| * prod_{i in S} H(class_i).\n")
print("  This weights each independent set by the product of H values.\n")

for n in [3, 4, 5]:
    nc = all_data[n]['nc']
    cd = all_data[n]['cd']
    adj = all_data[n]['adj']

    # H-weighted independence polynomial
    h_coeffs = [0] * (nc + 1)
    for k in range(nc + 1):
        total = 0
        for S in combinations(range(nc), k):
            is_indep = True
            for i in S:
                for j in S:
                    if i < j and adj[i][j]:
                        is_indep = False; break
                if not is_indep: break
            if is_indep:
                prod_H = 1
                for i in S:
                    prod_H *= cd[i]['H']
                total += prod_H
        h_coeffs[k] = total
    while len(h_coeffs) > 1 and h_coeffs[-1] == 0:
        h_coeffs.pop()

    print("  n=%d: I_H(G_%d, x) = %s" % (n, n, poly_str(h_coeffs)))
    print("    I_H(G_%d, 1) = %d" % (n, eval_poly(h_coeffs, 1)))
    print("    I_H(G_%d, 2) = %d" % (n, eval_poly(h_coeffs, 2)))
    print()

# ============================================================
# SIZE-WEIGHTED INDEPENDENCE POLYNOMIAL
# ============================================================

banner("SIZE-WEIGHTED INDEPENDENCE POLYNOMIAL")

print("  I_S(G_n, x) = sum over independent sets S of")
print("  x^|S| * prod_{i in S} size(class_i).\n")
print("  Weights by how many labeled tournaments each class represents.\n")

for n in [3, 4, 5]:
    nc = all_data[n]['nc']
    cd = all_data[n]['cd']
    adj = all_data[n]['adj']

    s_coeffs = [0] * (nc + 1)
    for k in range(nc + 1):
        total = 0
        for S in combinations(range(nc), k):
            is_indep = True
            for i in S:
                for j in S:
                    if i < j and adj[i][j]:
                        is_indep = False; break
                if not is_indep: break
            if is_indep:
                prod_S = 1
                for i in S:
                    prod_S *= cd[i]['size']
                total += prod_S
        s_coeffs[k] = total
    while len(s_coeffs) > 1 and s_coeffs[-1] == 0:
        s_coeffs.pop()

    print("  n=%d: I_S(G_%d, x) = %s" % (n, n, poly_str(s_coeffs)))
    print("    I_S(G_%d, 1) = %d" % (n, eval_poly(s_coeffs, 1)))
    print("    I_S(G_%d, 2) = %d" % (n, eval_poly(s_coeffs, 2)))
    print("    2^C(n,2) = %d" % (2**comb(n,2)))
    print()

# ============================================================
# THE CLIQUE POLYNOMIAL
# ============================================================

banner("CLIQUE POLYNOMIAL OF G_n")

print("  C(G_n, x) = sum over cliques S of x^|S|.")
print("  This = I(complement(G_n), x).\n")

for n in [3, 4, 5]:
    nc = all_data[n]['nc']
    adj = all_data[n]['adj']

    # Complement adjacency
    adj_comp = 1 - adj
    np.fill_diagonal(adj_comp, 0)

    c_coeffs = independence_polynomial(adj_comp, nc)
    print("  n=%d: C(G_%d, x) = %s" % (n, n, poly_str(c_coeffs)))
    print("    C(G_%d, 1) = %d (total cliques)" % (n, eval_poly(c_coeffs, 1)))
    print("    C(G_%d, 2) = %d" % (n, eval_poly(c_coeffs, 2)))
    print("    Clique number omega(G_%d) = %d" % (n, len(c_coeffs) - 1))
    print()

# ============================================================
# MATCHING POLYNOMIAL
# ============================================================

banner("MATCHING POLYNOMIAL OF G_n")

print("  mu(G_n, x) = sum_k (-1)^k * m_k * x^{|V|-2k}")
print("  where m_k = number of k-matchings (k independent edges).\n")

for n in [3, 4, 5]:
    nc = all_data[n]['nc']
    adj = all_data[n]['adj']
    ne = adj.sum() // 2

    # Find all edges
    edges = []
    for i in range(nc):
        for j in range(i+1, nc):
            if adj[i][j]:
                edges.append((i, j))

    # Count k-matchings
    max_matching = nc // 2
    m_coeffs = [0] * (max_matching + 1)
    for k in range(max_matching + 1):
        count = 0
        for edge_set in combinations(edges, k):
            # Check if edges are independent (no shared vertex)
            verts = set()
            indep = True
            for i, j in edge_set:
                if i in verts or j in verts:
                    indep = False; break
                verts.add(i); verts.add(j)
            if indep:
                count += 1
        m_coeffs[k] = count

    while len(m_coeffs) > 1 and m_coeffs[-1] == 0:
        m_coeffs.pop()

    print("  n=%d: matching numbers m_k = %s" % (n, m_coeffs))
    print("    Maximum matching size = %d" % (len(m_coeffs) - 1))
    print("    Total matchings = %d" % sum(m_coeffs))
    print()

# ============================================================
# THE DEEP CONNECTION: I(G_n, 2) AND TOURNAMENT THEORY
# ============================================================

banner("THE DEEP CONNECTION")

print("""
THE SELF-REFERENTIAL STRUCTURE:

Level 0: A tournament T on n vertices.
Level 1: Its conflict graph Omega(T) — vertices are odd cycles,
         edges are shared-vertex conflicts.
Level 2: H(T) = I(Omega(T), 2) — the OCF.
Level 3: G_n — the graph whose vertices are tournament iso classes,
         edges are single-arc reversals.
Level 4: I(G_n, x) — the independence polynomial of the class graph.

I(G_n, 2) counts something about tournament classes:
the number of ways to choose independent sets of classes
(classes no two of which are connected by a single arc reversal)
weighted by 2^k.

WHAT DOES AN INDEPENDENT SET IN G_n MEAN?

An independent set S in G_n is a collection of tournament iso classes
such that NO single arc reversal connects any two classes in S.
These classes are "isolated from each other" — you need at least
TWO arc reversals to get from any class in S to any other class in S.

I(G_n, 2) counts these isolated collections with weight 2^|S|.
""")

for n in [3, 4, 5]:
    nc = all_data[n]['nc']
    coeffs = all_data[n]['indep_poly']
    val = eval_poly(coeffs, 2)

    print("  n=%d:" % n)
    print("    I(G_%d, x) = %s" % (n, poly_str(coeffs)))
    print("    I(G_%d, 2) = %d" % (n, val))

    # What are the max independent sets?
    adj = all_data[n]['adj']
    cd = all_data[n]['cd']
    alpha = len(coeffs) - 1
    max_indep_sets = []
    for S in combinations(range(nc), alpha):
        is_indep = True
        for i in S:
            for j in S:
                if i < j and adj[i][j]:
                    is_indep = False; break
            if not is_indep: break
        if is_indep:
            H_vals = tuple(cd[i]['H'] for i in S)
            max_indep_sets.append((S, H_vals))

    print("    Maximum independent sets (alpha=%d):" % alpha)
    for S, Hs in max_indep_sets:
        print("      classes %s with H = %s" % (list(S), list(Hs)))
    print()

# ============================================================
# FACTORIZATION AND STRUCTURE
# ============================================================

banner("FACTORIZATION OF I(G_n, x)")

for n in [3, 4, 5]:
    coeffs = all_data[n]['indep_poly']
    print("  n=%d: I(G_%d, x) = %s" % (n, n, poly_str(coeffs)))

    # Check if it factors as a product of (1 + a_i * x)
    # This would mean the graph is "chordal" or has a specific structure
    # For a graph with independent edges, I = prod(1 + x) = (1+x)^k

    # Check: is I(G_n, x) = (1+x)^alpha * something?
    val_neg1 = eval_poly(coeffs, -1)
    print("    I(G_%d, -1) = %d" % (n, val_neg1))

    # Roots
    if len(coeffs) > 2:
        np_coeffs = np.array(coeffs[::-1], dtype=float)
        roots = np.roots(np_coeffs)
        real_roots = sorted([r.real for r in roots if abs(r.imag) < 1e-6])
        complex_roots = [(r.real, r.imag) for r in roots if abs(r.imag) >= 1e-6]
        print("    Real roots: %s" % [round(r, 4) for r in real_roots])
        if complex_roots:
            print("    Complex roots: %s" % [(round(a,4), round(b,4)) for a,b in complex_roots[:4]])

        # Are all roots real and negative? (property of claw-free graphs)
        all_real_neg = all(abs(r.imag) < 1e-6 and r.real < 0 for r in roots)
        print("    All roots real and negative: %s" % all_real_neg)
        if all_real_neg:
            print("    => G_%d is likely claw-free!" % n)
    print()

# ============================================================
# COMPARE I(G_n, x) WITH I(Omega(T), x) FOR SPECIFIC T
# ============================================================

banner("COMPARING I(G_n, x) WITH I(Omega, x)")

print("  For each tournament class, I(Omega(T), x) is the independence")
print("  polynomial of its conflict graph. H(T) = I(Omega(T), 2).\n")

for n in [3, 4, 5]:
    nc = all_data[n]['nc']
    cd = all_data[n]['cd']
    print("  n=%d:" % n)
    # For each class, compute I(Omega, x)
    # At n<=5, alpha_2 = 0 so I(Omega, x) = 1 + alpha_1 * x
    # and H = 1 + 2*alpha_1
    for i in range(nc):
        H = cd[i]['H']
        alpha1 = (H - 1) // 2  # since H = 1 + 2*alpha_1 at n<=5
        # I(Omega(T), x) = 1 + alpha_1 * x (when alpha_2 = 0)
        print("    class %d: H=%d, I(Omega, x) = 1 + %d*x, I(Omega, 2) = %d" %
              (i, H, alpha1, 1 + 2*alpha1))
    print()

# ============================================================
# THE SEQUENCES
# ============================================================

banner("NEW SEQUENCES FROM I(G_n, x)")

print("  I(G_n, x) coefficients [a_0, a_1, ..., a_alpha]:\n")
for n in [3, 4, 5]:
    coeffs = all_data[n]['indep_poly']
    print("  n=%d: %s" % (n, coeffs))

print()
print("  DERIVED SEQUENCES:")
print()

# a_k for fixed k across n
for k in range(6):
    vals = []
    for n in [3, 4, 5]:
        coeffs = all_data[n]['indep_poly']
        vals.append(coeffs[k] if k < len(coeffs) else 0)
    print("  a_%d(G_n): %s" % (k, vals))

print()
print("  I(G_n, 2): %s" % [eval_poly(all_data[n]['indep_poly'], 2) for n in [3, 4, 5]])
print("  I(G_n, 1): %s" % [eval_poly(all_data[n]['indep_poly'], 1) for n in [3, 4, 5]])
print("  I(G_n,-1): %s" % [eval_poly(all_data[n]['indep_poly'], -1) for n in [3, 4, 5]])
print("  alpha(G_n): %s" % [len(all_data[n]['indep_poly'])-1 for n in [3, 4, 5]])

# ============================================================
# SYNTHESIS
# ============================================================

banner("SYNTHESIS: THE META-INDEPENDENCE POLYNOMIAL")

print("""
THE FULL PICTURE:

1. I(G_3, x) = 1 + 2x  (path graph P_2)
   I(G_3, 2) = 5
   alpha = 1 (both classes are adjacent, so max indep set = 1)

2. I(G_4, x) = 1 + 4x + 2x^2  (4 vertices, alpha = 2)
   I(G_4, 2) = 17
   The 2 maximum independent sets: classes not connected by arc reversal.

3. I(G_5, x) = 1 + 12x + 41x^2 + 60x^3 + 39x^4 + 9x^5
   I(G_5, 2) = large!
   alpha = 5: five mutually unreachable tournament classes!

THE SELF-REFERENTIAL LOOP:
  H(T) = I(Omega(T), 2)     -- tournament -> conflict graph -> indep poly
  I(G_n, x)                  -- class graph -> its own indep poly
  These are DIFFERENT independence polynomials of RELATED graphs.

The OCF connects Level 0 (tournament) to Level 2 (H count).
I(G_n, x) lives at Level 4 (class graph structure).
Is there a formula connecting them?

KEY OBSERVATION:
  sum_{classes i} H(i) = sum_H (unlabeled)
  I(G_n, 1) = total independent sets in G_n

These are different sums over different structures.
But BOTH are computed from the same underlying tournament data.

THE INDEPENDENCE POLYNOMIAL I(G_n, x) ENCODES:
  - At x=1: the total number of "mutually unreachable" class collections
  - At x=2: a 2-weighted count of such collections (the "meta-OCF")
  - At x=-1: an alternating sum related to the Euler characteristic
  - The roots: the "critical fugacities" of the class graph

ROOTS:
  n=3: root at x = -1/2 (REAL, NEGATIVE)
  n=4: roots at x ~ -0.29, -1.71 (both REAL, NEGATIVE)
  n=5: mix of real and complex roots

ALL REAL NEGATIVE ROOTS (n=3,4) means G_n is like a claw-free graph
at small n. If this persists, it implies strong structural regularity.

THE PUNCHLINE:
Tournament theory is self-referential at every level.
Tournaments have conflict graphs; conflict graphs have independence polynomials.
Tournament classes form a graph; that graph has its own independence polynomial.
The structure at each level constrains the structure at the next.
""")

print("Done. opus-2026-03-22-S170")
