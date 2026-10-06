#!/usr/bin/env python3
"""
trees5_dynkin_kervaire_milnor_20261005.py

Owner's prompt: Poincaré conjecture, exotic spheres, 992 and 16256, Milnor's S^7 (28 total) against the
trees-on-5-vertices work.  Exact facts frozen here:

 A. The three unlabeled trees on 5 vertices are Dynkin-type diagrams: path = A_5, fork S(2,1,1) = D_5,
    star K_{1,4} = affine D_4 (D4~).  Plumbing D^4-bundles over S^4 with Euler number 2 along a tree T
    gives an 8-manifold whose intersection form is the Cartan matrix C(T) = 2I - Adj(T); its boundary
    is a homotopy 7-sphere iff det C(T) = +-1.  det A_n = n+1, det D_n = 4, det E_6/E_7/E_8 = 3/2/1,
    det(affine) = 0.  So among trees on <= 8 vertices only E_8 plumbs to a homotopy sphere (Milnor's
    generator of Theta_7 = Z/28); the three 5-vertex trees plumb to boundaries with H_3 = Z/6, Z/4 and
    an infinite group.
 B. Kervaire--Milnor: |bP_{4k}| = 2^(2k-2) (2^(2k-1) - 1) numerator(4 B_k / k) with B_k = |B_{2k}| the
    Bernoulli numbers (B_1 = 1/6, B_2 = 1/30, ...).  Values 28, 992, 8128, 261632, 1448448, 67100672 for
    k = 2..7; Theta_7 = bP_8 = 28, Theta_11 = bP_12 = 992, Theta_15 = bP_16 + Z/2 = 16256.
    The factor 2^(2k-2)(2^(2k-1)-1) is Euclid's form 2^(p-1)(2^p - 1) with p = 2k-1, perfect exactly when
    2^p - 1 is a Mersenne prime: p = 3, 5, 7 (28, 496, 8128), fails at p = 9, 11, returns at p = 13
    (33550336, k = 7).  That is the whole content of "992 = 2*496, 16256 = 2*8128, 28 perfect".
 C. The numeral 28 in our work: C(8,2) = 28 arcs of the Sumner host (8-tournaments), the E_8 tree has 8
    vertices; no map carries |Theta_7| to either (NUMEROLOGY).  The one exact bridge is A: trees are
    plumbing graphs, and only the unimodular tree E_8 bounds an exotic sphere.
Reproduce: python3 trees5_dynkin_kervaire_milnor_20261005.py
"""
import sys, math
from fractions import Fraction
import networkx as nx

CHECKS = 0
def check(c, msg):
    global CHECKS
    CHECKS += 1
    if not c:
        print("CHECK FAILED:", msg); sys.exit(1)

def det_int(M):
    # Bareiss algorithm on integer matrix
    n = len(M); A = [row[:] for row in M]; sign = 1; prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]

def cartan(T):
    nodes = list(T.nodes); idx = {v: i for i, v in enumerate(nodes)}
    n = len(nodes); C = [[0] * n for _ in range(n)]
    for i in range(n): C[i][i] = 2
    for u, v in T.edges: C[idx[u]][idx[v]] = -1; C[idx[v]][idx[u]] = -1
    return C

def tree_name(T):
    degs = sorted((d for _, d in T.degree()), reverse=True)
    if degs[0] == T.number_of_nodes() - 1: return "star"
    if degs[0] == 2: return "path"
    return "fork" if T.number_of_nodes() == 5 else "other"

print("=" * 78)
print("A. Cartan determinants of trees (plumbing intersection forms with Euler number 2)")
print("=" * 78)
# named diagrams
def path(n): return nx.path_graph(n)
def Dn(n):   # path on n-1 vertices with an extra leaf at the second vertex
    T = nx.path_graph(n - 1); T.add_edge(1, n - 1); return T
def En(n):   # path on n-1 vertices with an extra leaf at the third vertex
    T = nx.path_graph(n - 1); T.add_edge(2, n - 1); return T
named = {"A_5": path(5), "D_5": Dn(5), "D4~ (star K_{1,4})": nx.star_graph(4), "E_6": En(6), "E_7": En(7), "E_8": En(8),
         "D_4": Dn(4), "D_8": Dn(8), "A_8": path(8), "E7~ (8 vertices)": None, "D7~ (8 vertices)": None}
E7t = nx.path_graph(7); E7t.add_edge(3, 7)      # affine E_7: path on 7 with a leaf at the middle
D7t = nx.path_graph(6); D7t.add_edge(1, 6); D7t.add_edge(4, 7)   # affine D_7: 8 vertices, two forks
named["E7~ (8 vertices)"] = E7t; named["D7~ (8 vertices)"] = D7t
for nm, T in named.items():
    d = det_int(cartan(T))
    print(f"  {nm:22s} vertices={T.number_of_nodes()}  det C = {d:3d}   boundary of the plumbing: " +
          ("homotopy sphere" if abs(d) == 1 else (f"rational homology sphere, |H_3| = {abs(d)}" if d != 0 else "H_3 infinite (affine)")))
check(det_int(cartan(path(5))) == 6 and det_int(cartan(Dn(5))) == 4 and det_int(cartan(nx.star_graph(4))) == 0, "A5/D5/D4~ determinants")
check(det_int(cartan(En(8))) == 1 and det_int(cartan(En(7))) == 2 and det_int(cartan(En(6))) == 3, "E6/E7/E8")
trees5 = {tree_name(T): T for T in nx.nonisomorphic_trees(5)}
check(nx.is_isomorphic(trees5["path"], path(5)) and nx.is_isomorphic(trees5["fork"], Dn(5)) and nx.is_isomorphic(trees5["star"], nx.star_graph(4)), "the three 5-trees are A_5, D_5, D4~")
print("  All 23 trees on 8 vertices: determinants of their Cartan matrices (positive definite iff Dynkin: A_8, D_8, E_8)")
cnt = {}
for T in nx.nonisomorphic_trees(8):
    d = det_int(cartan(T)); cnt[d] = cnt.get(d, 0) + 1
print("   det -> count:", dict(sorted(cnt.items())))
check(cnt.get(1, 0) == 1, "exactly one unimodular tree on 8 vertices (E_8)")
# which 8-vertex trees are positive definite (Dynkin)? check eigenvalues via Sylvester (leading minors)
def pos_def(C):
    n = len(C)
    return all(det_int([row[:k] for row in C[:k]]) > 0 for k in range(1, n + 1))
dynkin8 = [T for T in nx.nonisomorphic_trees(8) if pos_def(cartan(T))]
print(f"   positive definite (finite Dynkin) trees on 8 vertices: {len(dynkin8)}  (A_8, D_8, E_8)")
check(len(dynkin8) == 3, "A_8, D_8, E_8")

print("\n" + "=" * 78)
print("B. Kervaire--Milnor orders |bP_4k| = 2^(2k-2)(2^(2k-1)-1) num(4 B_k/k), perfect-number factor")
print("=" * 78)
def bernoulli(nmax):
    # B_0..B_nmax with B_1 = -1/2 convention; returns list of Fractions
    B = [Fraction(0)] * (nmax + 1); B[0] = Fraction(1)
    for m in range(1, nmax + 1):
        B[m] = -sum(Fraction(math.comb(m + 1, j)) * B[j] for j in range(m)) / (m + 1)
    return B
B = bernoulli(16)
known = {2: 28, 3: 992, 4: 8128, 5: 261632, 6: 1448424448, 7: 67100672}   # A001676: Theta_19 = 523264 = 2*261632, Theta_23 = 69524373504 = 48*1448424448, Theta_27 contains bP_28
def is_prime(n):
    if n < 2: return False
    i = 2
    while i * i <= n:
        if n % i == 0: return False
        i += 1
    return True
for k in range(2, 8):
    Bk = abs(B[2 * k])
    euclid = 2 ** (2 * k - 2) * (2 ** (2 * k - 1) - 1)
    numer = (4 * Bk / k).numerator
    order = euclid * numer
    p = 2 * k - 1
    perfect = is_prime(2 ** p - 1)
    print(f"  k={k}: dim 4k-1 = {4*k-1:2d}  B_k = {Bk}  Euclid factor 2^{2*k-2}(2^{p}-1) = {euclid:9d} ({'perfect: 2^'+str(p)+'-1 prime' if perfect else 'not perfect: 2^'+str(p)+'-1 = '+str(2**p-1)+' composite'})  num(4B_k/k) = {numer}  |bP_{4*k}| = {order}")
    check(order == known[k], f"|bP_{4*k}|")
print("  Theta_7 = bP_8 = 28;  Theta_11 = bP_12 = 992;  Theta_15 = bP_16 + Z/2 = 2*8128 = 16256.")
print("  28 = 2^2*7, 992 = 2*2^4*31, 16256 = 2*2^6*127: the perfect numbers enter through Euclid's factor,")
print("  and only because 7, 31, 127 are Mersenne primes; at k = 5, 6 (511, 2047 composite) the factor is not perfect.")

print("\n" + "=" * 78)
print("C. The numeral 28 in this session's objects")
print("=" * 78)
print(f"  C(8,2) = {math.comb(8,2)} arcs of an 8-tournament (the Sumner host at n = 5); E_8 has 8 vertices;")
print("  2^28 = labeled 8-tournaments; 6880 classes.  No map carries |Theta_7| = 28 to either: NUMEROLOGY.")
print(f"\nALL {CHECKS} CHECKS PASSED")
