#!/usr/bin/env python3
"""
TEST: Does even-r-powers require the COMPLETE GRAPH?

If we set some edge weights to zero (removing edges), does M[a,b] still
have only even r-powers? If NO, completeness is essential for the proof.

Also: test whether the property holds for WEIGHTED complete bipartite graphs,
or other graph families.

kind-pasteur-2026-03-06-S24
"""

from itertools import permutations
from sympy import symbols, expand, Symbol, Poly, Rational
from collections import defaultdict

def make_variables(n):
    r = Symbol('r')
    s = {}
    for i in range(n):
        for j in range(n):
            if i < j:
                s[(i,j)] = Symbol(f's{i}{j}')
                s[(j,i)] = -s[(i,j)]
    return r, s

def compute_M_general(a, b, n, r, edge_weights):
    """
    Compute M[a,b] with GENERAL edge weights (not necessarily c-tournament).
    edge_weights[(i,j)] = weight of edge i->j.
    """
    U = [v for v in range(n) if v != a and v != b]
    total = 0
    for mask in range(2**len(U)):
        S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
        R = [u for u in U if u not in S]
        sign = (-1)**len(S)
        Sa = set(S + [a])
        Rb = set(R + [b])

        ea = 0
        for p in permutations(sorted(Sa)):
            if p[-1] != a: continue
            w = 1
            for k in range(len(p)-1):
                w *= edge_weights.get((p[k], p[k+1]), 0)
            ea += w
        if len(Sa) == 1: ea = 1

        bb = 0
        for p in permutations(sorted(Rb)):
            if p[0] != b: continue
            w = 1
            for k in range(len(p)-1):
                w *= edge_weights.get((p[k], p[k+1]), 0)
            bb += w
        if len(Rb) == 1: bb = 1

        total += sign * ea * bb
    return expand(total)

# ==============================================================
# Test 1: Remove one edge from a c-tournament
# ==============================================================
print("=" * 70)
print("Test 1: Remove edges from c-tournament")
print("=" * 70)

n = 4
r, s = make_variables(n)

# Full c-tournament: t(i,j) = r + s_{ij}, t(j,i) = r - s_{ij}
full_weights = {}
for i in range(n):
    for j in range(n):
        if i != j:
            full_weights[(i,j)] = r + s[(i,j)] if i < j else r - s[(j,i)]

M_full = compute_M_general(0, 1, n, r, full_weights)
print(f"\nn=4 full c-tournament: M[0,1] = {M_full}")
p_full = Poly(M_full, r)
odd_powers = [k for k in range(p_full.degree()+1) if k % 2 == 1 and expand(p_full.nth(k)) != 0]
print(f"  Odd r-powers: {odd_powers if odd_powers else 'NONE'}")

# Remove edge (2,3): set t(2,3) = 0 and t(3,2) = 0
# This breaks the c-tournament property since t(2,3)+t(3,2) = 0 != c
missing_weights = dict(full_weights)
missing_weights[(2,3)] = 0
missing_weights[(3,2)] = 0

M_missing = compute_M_general(0, 1, n, r, missing_weights)
print(f"\nn=4 with edge (2,3) removed: M[0,1] = {M_missing}")
if M_missing != 0:
    p_miss = Poly(M_missing, r)
    odd_powers = [k for k in range(p_miss.degree()+1) if k % 2 == 1 and expand(p_miss.nth(k)) != 0]
    print(f"  Odd r-powers: {odd_powers if odd_powers else 'NONE'}")
    for k in range(p_miss.degree()+1):
        print(f"  [r^{k}] = {expand(p_miss.nth(k))}")

# Remove edge (0,2): endpoint-adjacent edge
missing_weights2 = dict(full_weights)
missing_weights2[(0,2)] = 0
missing_weights2[(2,0)] = 0

M_missing2 = compute_M_general(0, 1, n, r, missing_weights2)
print(f"\nn=4 with edge (0,2) removed: M[0,1] = {M_missing2}")
if M_missing2 != 0:
    p_miss2 = Poly(M_missing2, r)
    odd_powers = [k for k in range(p_miss2.degree()+1) if k % 2 == 1 and expand(p_miss2.nth(k)) != 0]
    print(f"  Odd r-powers: {odd_powers if odd_powers else 'NONE'}")

# ==============================================================
# Test 2: c-tournament but with NON-SKEW-SYMMETRIC s
# ==============================================================
print("\n" + "=" * 70)
print("Test 2: Keep c-tournament constraint, break skew-symmetry of s")
print("=" * 70)

# What if t(i,j) + t(j,i) = c but s_{ij} != -s_{ji}?
# Then t(i,j) = r + alpha_{ij}, t(j,i) = r + beta_{ij} with alpha + beta = 0.
# So beta_{ij} = -alpha_{ij}. This IS skew-symmetric.
# So the c-tournament constraint FORCES skew-symmetry of s.
print("  c-tournament constraint t(i,j)+t(j,i)=2r forces s skew-symmetric.")
print("  No separate test needed.")

# ==============================================================
# Test 3: Directed weights that DON'T satisfy c-tournament
# ==============================================================
print("\n" + "=" * 70)
print("Test 3: Non-c-tournament (arbitrary directed weights)")
print("=" * 70)

n = 4
r, s = make_variables(n)

# Arbitrary weights: t(i,j) = r + s_{ij} for i < j,
# but t(j,i) = r + s_{ji} where s_{ji} is INDEPENDENT of s_{ij}
# (not forced to be -s_{ij})

# For this test: t(0,1)=r+alpha, t(1,0)=r+beta, where alpha+beta != 0
alpha = Symbol('alpha')
beta = Symbol('beta')
arb_weights = {}
for i in range(n):
    for j in range(n):
        if i != j:
            arb_weights[(i,j)] = r + s[(min(i,j), max(i,j))] * (1 if i < j else -1)

# This is still a c-tournament. Let me make it non-c:
# Set t(2,3) = r + s23 + delta (asymmetric perturbation)
delta = Symbol('delta')
pert_weights = dict(full_weights)
pert_weights[(2,3)] = r + s[(2,3)] + delta
pert_weights[(3,2)] = r - s[(2,3)] - delta  # Still c-tournament: sum = 2r

M_pert = compute_M_general(0, 1, n, r, pert_weights)
print(f"\nn=4 with perturbation delta on (2,3): M[0,1] = {M_pert}")
p_pert = Poly(expand(M_pert), r)
odd_powers = [k for k in range(p_pert.degree()+1) if k % 2 == 1 and expand(p_pert.nth(k)) != 0]
print(f"  Odd r-powers: {odd_powers if odd_powers else 'NONE'}")
# Note: if delta just shifts s23, this is still a c-tournament
# so even r-powers should hold.

# Now TRULY break c-tournament: t(2,3) = r + s23 + delta, t(3,2) = r - s23 (no delta)
# So t(2,3) + t(3,2) = 2r + delta != 2r
pert_weights2 = dict(full_weights)
pert_weights2[(2,3)] = r + s[(2,3)] + delta
# Keep t(3,2) = r - s23 (unchanged)

M_pert2 = compute_M_general(0, 1, n, r, pert_weights2)
print(f"\nn=4 with t(2,3)+t(3,2)=2r+delta: M[0,1] = {M_pert2}")
p_pert2 = Poly(expand(M_pert2), r)
for k in range(p_pert2.degree()+1):
    ck = expand(p_pert2.nth(k))
    # Check if there's a delta-dependent odd-r term
    if k % 2 == 1:
        print(f"  [r^{k}] = {ck}")
        delta_coeff = expand(ck.coeff(delta))
        if delta_coeff != 0:
            print(f"    delta coefficient: {delta_coeff} != 0 --> ODD r-power from c-breaking!")

# ==============================================================
# Test 4: Does the property hold for t(i,j) = r + s(i,j), s skew-symmetric,
# but on a NON-COMPLETE graph?
# ==============================================================
print("\n" + "=" * 70)
print("Test 4: Skew-symmetric weights on non-complete graph")
print("=" * 70)

n = 5
r, s = make_variables(n)

# Path graph: 0-2-3-4-1 (edges: 02, 20, 23, 32, 34, 43, 41, 14)
# Plus endpoint edges: 01, 10 (keep these)
# Missing edges: 03, 04, 12, 13, 24

path_weights = {}
for i in range(n):
    for j in range(n):
        if i != j:
            path_weights[(i,j)] = 0  # start with no edges

# Add path graph edges with c-tournament weights
path_edges = [(0,2), (2,3), (3,4), (4,1), (0,1)]
for u, v in path_edges:
    path_weights[(u,v)] = r + s[(min(u,v), max(u,v))] * (1 if u < v else -1)
    path_weights[(v,u)] = r - s[(min(u,v), max(u,v))] * (1 if u < v else -1)

M_path = compute_M_general(0, 1, n, r, path_weights)
print(f"\nn=5 path graph (0-2-3-4-1): M[0,1] = {M_path}")
if M_path != 0:
    p_path = Poly(M_path, r)
    for k in range(p_path.degree()+1):
        ck = expand(p_path.nth(k))
        if ck != 0:
            print(f"  [r^{k}] = {ck}")
    odd_powers = [k for k in range(p_path.degree()+1) if k % 2 == 1 and expand(p_path.nth(k)) != 0]
    print(f"  Odd r-powers: {odd_powers if odd_powers else 'NONE'}")
else:
    print("  M = 0 (no Hamiltonian paths through this sparse graph)")

# Complete bipartite graph: K_{2,3} between {0,1} and {2,3,4}
bip_weights = {}
for i in range(n):
    for j in range(n):
        if i != j:
            bip_weights[(i,j)] = 0

for i in [0, 1]:
    for j in [2, 3, 4]:
        bip_weights[(i,j)] = r + s[(min(i,j), max(i,j))] * (1 if i < j else -1)
        bip_weights[(j,i)] = r - s[(min(i,j), max(i,j))] * (1 if i < j else -1)

M_bip = compute_M_general(0, 1, n, r, bip_weights)
print(f"\nn=5 complete bipartite K_{{2,3}}: M[0,1] = {M_bip}")
if M_bip != 0:
    p_bip = Poly(M_bip, r)
    for k in range(p_bip.degree()+1):
        ck = expand(p_bip.nth(k))
        if ck != 0:
            print(f"  [r^{k}] = {ck}")
    odd_powers = [k for k in range(p_bip.degree()+1) if k % 2 == 1 and expand(p_bip.nth(k)) != 0]
    print(f"  Odd r-powers: {odd_powers if odd_powers else 'NONE'}")
else:
    print("  M = 0")

# ==============================================================
# Test 5: Cycle graph
# ==============================================================
# Complete graph minus one edge at n=4
n = 4
r, s = make_variables(n)
full_w = {}
for i in range(n):
    for j in range(n):
        if i != j:
            full_w[(i,j)] = r + s[(min(i,j), max(i,j))] * (1 if i < j else -1)

# Remove edge (2,3)
almost_w = dict(full_w)
almost_w[(2,3)] = 0
almost_w[(3,2)] = 0

M_almost = compute_M_general(0, 1, n, r, almost_w)
print(f"\nn=4 K4 minus edge (2,3): M[0,1] = {M_almost}")
if M_almost != 0:
    p_alm = Poly(M_almost, r)
    for k in range(p_alm.degree()+1):
        ck = expand(p_alm.nth(k))
        if ck != 0:
            parity = "EVEN" if k % 2 == 0 else "ODD"
            print(f"  [r^{k}] = {ck}  ({parity})")
    odd_powers = [k for k in range(p_alm.degree()+1) if k % 2 == 1 and expand(p_alm.nth(k)) != 0]
    print(f"  Odd r-powers: {odd_powers if odd_powers else 'NONE'}")

# Also check M[0,1] vs M[1,0]
M_almost_10 = compute_M_general(1, 0, n, r, almost_w)
print(f"  M[1,0] = {M_almost_10}")
print(f"  M[0,1] = M[1,0]? {expand(M_almost - M_almost_10) == 0}")

print("\n" + "=" * 70)
print("CONCLUSION")
print("=" * 70)
print("""
If removing edges creates odd r-powers, then COMPLETENESS is essential.
If breaking c-tournament creates odd r-powers, then the c-constraint is essential.
Both completeness and the c-tournament constraint together force even r-powers.
""")
