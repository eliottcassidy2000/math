#!/usr/bin/env python3
"""
WRONSKIAN / EVEN-ODD DECOMPOSITION APPROACH to M[a,b] = M[b,a]

KEY REFORMULATION:
  M[a,b] = (-1)^{n-2} * G(r)
  where G(r) = sum_S E_a(S+a; r) * E_b(R+b; -r)

  M[a,b] = M[b,a]  iff  G(r) = G(-r)  (G is even in r)

  Write f_S(r) = E_a(S+a; r), g_S(r) = E_b(R+b; r).
  Then G(r) = sum_S f_S(r) * g_S(-r).
  G(-r) = sum_S f_S(-r) * g_S(r).

  G(r) - G(-r) = sum_S [f_S(r)*g_S(-r) - f_S(-r)*g_S(r)]

  Decompose f_S = f_S^+ + f_S^- (even/odd parts in r).
  Similarly g_S = g_S^+ + g_S^-.

  Then: f(r)g(-r) - f(-r)g(r) = 2[f^- g^+ - f^+ g^-]

  So G(r) = G(-r) iff sum_S [f_S^- * g_S^+ - f_S^+ * g_S^-] = 0.

  This is a "Wronskian-like" cancellation between the even and odd
  parts of the two endpoint Hamiltonian path polynomials.

QUESTION: Does this cancel term-by-term (for each S) or only after summing?
If term-by-term: f_S^-/f_S^+ = g_S^-/g_S^+ (same "slope ratio" for all S).
If only after summing: there's a more subtle identity.

kind-pasteur-2026-03-06-S24
"""

from itertools import permutations
from sympy import symbols, expand, Symbol, Poly, Rational
from collections import defaultdict

def make_tournament(n):
    r = Symbol('r')
    s = {}
    for i in range(n):
        for j in range(n):
            if i < j:
                s[(i,j)] = Symbol(f's{i}{j}')
                s[(j,i)] = -s[(i,j)]
    def t(i, j):
        if i == j: return 0
        return r + s[(i,j)]
    return r, s, t

def E_v(vertex_set, target, r, s):
    """Ham paths through vertex_set ending at target. Weight = prod(r + s_e)."""
    vs = list(vertex_set)
    if len(vs) == 1:
        return 1
    result = 0
    for perm in permutations([v for v in vs if v != target]):
        path = list(perm) + [target]
        w = 1
        for k in range(len(path)-1):
            u, v = path[k], path[k+1]
            w *= (r + s[(u,v)])
        result += w
    return expand(result)

def even_odd_parts(poly_expr, r):
    """Decompose polynomial into even and odd parts in r."""
    p = Poly(poly_expr, r)
    even_part = 0
    odd_part = 0
    for k in range(p.degree() + 1):
        coeff = expand(p.nth(k))
        if coeff == 0: continue
        if k % 2 == 0:
            even_part += coeff * r**k
        else:
            odd_part += coeff * r**k
    return expand(even_part), expand(odd_part)

print("=" * 70)
print("WRONSKIAN DECOMPOSITION OF G(r) = sum_S f_S(r) g_S(-r)")
print("=" * 70)

for n in [4, 5]:
    r, s, t = make_tournament(n)
    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]

    print(f"\n{'='*60}")
    print(f"n = {n}, a={a}, b={b}, U={U}")
    print(f"{'='*60}")

    G_r = 0
    wronskian_total = 0

    for mask in range(2**len(U)):
        S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
        R = [u for u in U if u not in S]

        f_S = E_v(set(S + [a]), a, r, s)  # E_a(S+a; r)
        g_S = E_v(set(R + [b]), b, r, s)  # E_b(R+b; r)

        # G contribution: f_S(r) * g_S(-r)
        from sympy import sympify
        g_neg = expand(sympify(g_S).subs(r, -r))
        contrib = expand(f_S * g_neg)
        G_r = expand(G_r + contrib)

        # Wronskian contribution: f^- g^+ - f^+ g^-
        f_even, f_odd = even_odd_parts(f_S, r)
        g_even, g_odd = even_odd_parts(g_S, r)

        wronsk = expand(f_odd * g_even - f_even * g_odd)
        wronskian_total = expand(wronskian_total + wronsk)

        print(f"\n  S={set(S)}: |S|={len(S)}, |R|={len(R)}")
        print(f"    f_S = E_a(S+a) = {f_S}")
        print(f"    g_S = E_b(R+b) = {g_S}")
        print(f"    f^+ = {f_even}, f^- = {f_odd}")
        print(f"    g^+ = {g_even}, g^- = {g_odd}")
        print(f"    f^- g^+ - f^+ g^- = {wronsk}")

        if wronsk == 0:
            print(f"    >> Wronskian = 0 (individual cancellation!)")
        else:
            # Check: does f^-/f^+ = g^-/g^+ (proportionality)?
            if f_even != 0 and g_even != 0:
                ratio_check = expand(f_odd * g_even - f_even * g_odd)
                if ratio_check == 0:
                    print(f"    >> Proportional: f^-/f^+ = g^-/g^+")

    print(f"\n  Total G(r) = {G_r}")
    G_neg = expand(G_r.subs(r, -r))
    print(f"  G(-r) = {G_neg}")
    print(f"  G(r) - G(-r) = {expand(G_r - G_neg)}")
    print(f"  G is even? {expand(G_r - G_neg) == 0}")

    print(f"\n  Sum of Wronskians: {wronskian_total}")
    print(f"  (Should be 0: {wronskian_total == 0})")

    # DEEPER ANALYSIS: pair S with complement U\S
    print(f"\n  Complement pairing analysis:")
    for mask in range(2**(len(U)-1)):  # only half, avoid double-counting
        S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
        R = [u for u in U if u not in S]
        comp_mask = (2**len(U) - 1) ^ mask

        f_S = E_v(set(S + [a]), a, r, s)
        g_S = E_v(set(R + [b]), b, r, s)
        f_R = E_v(set(R + [a]), a, r, s)  # E_a(R+a) -- complement
        g_R = E_v(set(S + [b]), b, r, s)  # E_b(S+b) -- complement

        f_Se, f_So = even_odd_parts(f_S, r)
        g_Se, g_So = even_odd_parts(g_S, r)
        f_Re, f_Ro = even_odd_parts(f_R, r)
        g_Re, g_Ro = even_odd_parts(g_R, r)

        w_S = expand(f_So * g_Se - f_Se * g_So)
        w_R = expand(f_Ro * g_Re - f_Re * g_Ro)

        pair_sum = expand(w_S + w_R)
        print(f"  S={set(S)} + S'={set(R)}: W_S={w_S != 0}, W_R={w_R != 0}, sum={pair_sum}")

        if pair_sum == 0 and w_S != 0:
            print(f"    >> COMPLEMENT PAIRING CANCELS!")
            print(f"    W_S = {w_S}")
            print(f"    W_R = {w_R}")

# ================================================================
# TEST: Is the Wronskian identity related to a DETERMINANT?
# ================================================================
print("\n" + "=" * 70)
print("TEST: Wronskian as determinant structure")
print("=" * 70)

n = 4
r, s, t = make_tournament(n)
a, b = 0, 1
U = [2, 3]

# For n=4, |U|=2 subsets: {}, {2}, {3}, {2,3}
# f_{} = E_a({0}) = 1, g_{} = E_b({1,2,3}) = paths ending at b=1 through {1,2,3}
# f_{2} = E_a({0,2}) = t(2,0), g_{2} = E_b({1,3}) = t(3,1)
# f_{3} = E_a({0,3}) = t(3,0), g_{3} = E_b({1,2}) = t(2,1)
# f_{2,3} = E_a({0,2,3}) = t(2,3)t(3,0) + t(3,2)t(2,0), g_{2,3} = E_b({1}) = 1

print(f"\nn=4 detailed:")
for mask in range(4):
    S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
    R = [u for u in U if u not in S]
    f = E_v(set(S + [a]), a, r, s)
    g = E_v(set(R + [b]), b, r, s)
    print(f"  S={set(S)}: f_S = {f}, g_S = {g}")

# The G function: G(r) = sum_S f_S(r) g_S(-r)
# = 1 * g_{}(-r) + (r+s20)*g_{2}(-r) + (r+s30)*g_{3}(-r)
#   + [(r+s23)(r+s30) + (r-s23)(r+s20)] * 1
#
# = g_{}(-r) + (r+s20)*(-r+s31) + (r+s30)*(-r+s21)
#   + (r+s23)(r+s30) + (r-s23)(r+s20)

# Let me compute each piece:
# g_{}(-r) = E_b({1,2,3}; -r)
#   Paths ending at 1 through {1,2,3}:
#   2->3->1: (-r+s23)(-r+s31) and 3->2->1: (-r-s23)(-r+s21)
g_empty_neg = expand((-r+s[(2,3)])*(-r+s[(3,1)]) + (-r+s[(3,2)])*(-r+s[(2,1)]))
print(f"\n  g_{{}}(-r) = {g_empty_neg}")

# The full G:
# = g_empty_neg + (r+s20)(-r+s31) + (r+s30)(-r+s21) + (r+s23)(r+s30) + (r-s23)(r+s20)
piece1 = expand((r+s[(2,0)])*(-r+s[(3,1)]))
piece2 = expand((r+s[(3,0)])*(-r+s[(2,1)]))
piece3 = expand((r+s[(2,3)])*(r+s[(3,0)]))
piece4 = expand((r-s[(2,3)])*(r+s[(2,0)]))

G_manual = expand(g_empty_neg + piece1 + piece2 + piece3 + piece4)
print(f"  G(r) manual = {G_manual}")

# Note: each piece has BOTH even and odd r-powers.
# But the sum should have ONLY even r-powers.

p_G = Poly(G_manual, r)
print(f"  G r-expansion:")
for k in range(p_G.degree() + 1):
    print(f"    [r^{k}] = {expand(p_G.nth(k))}")

# Check: [r^1] and [r^3] should be 0
assert expand(p_G.nth(1)) == 0, f"[r^1] = {expand(p_G.nth(1))} != 0"
assert expand(p_G.nth(3)) == 0 or p_G.degree() < 3, f"[r^3] nonzero"
print(f"\n  Confirmed: G has only even r-powers at n=4")

# Now: What IS G in matrix terms?
# G(r) = sum_S f_S(r) g_S(-r)
# where f_S(r) uses weights (r + s_e) and g_S(-r) uses weights (-r + s_e) = (s_e - r)

# So: G(r) = sum_S prod_{e in Pa}(r+s_e) * prod_{e in Pb}(s_e-r) * (sum over paths)

# For each edge e:
# If e is in P_a: weight = (r + s_e)
# If e is in P_b: weight = (s_e - r)

# Product over all edges = prod_a(r+s_e) * prod_b(s_e-r)
# = prod_a(r+s_e) * (-1)^|Pb| * prod_b(r-s_e)

# For |Pb| = |R| edges:
# = (-1)^|R| * prod_a(r+s_e) * prod_b(r-s_e)

# This is (-1)^|R| times a product where each factor is
# (r + eps_e * s_e) with eps = +1 for Pa, -1 for Pb.

# The product prod_e (r + eps_e * s_e) is a polynomial in r
# that is EVEN when all eps_e are simultaneously flipped
# (since flipping all eps = negating all s = negating r gets (-1)^m * original = (-1)^{n-2} * original).
# Wait, that's not right. Flipping all eps: prod(r - eps_e * s_e) vs prod(r + eps_e * s_e).
# These are different unless m = n-2 is even.

# Actually: prod(r + eps_e s_e) evaluated at -r = prod(-r + eps_e s_e) = (-1)^m prod(r - eps_e s_e).
# So: prod_{+eps}(-r) = (-1)^m prod_{-eps}(r).
# And summing over all (S, paths) with the proper signs:
# G(-r) = sum_S f_S(-r) g_S(r) = sum_S [(-1)^|S| f_S^{flip}(r)] * [(-1)^|R| g_S^{flip}(r)]
# Hmm, this is getting circular again.

# Let me try: NUMERICAL verification at n=6 with random s-values
print("\n" + "=" * 70)
print("NUMERICAL VERIFICATION: G(r) = G(-r) for n=6,7,8")
print("=" * 70)

import random
random.seed(42)

def make_numeric_tournament(n):
    """Random skew-symmetric s-values."""
    s = {}
    for i in range(n):
        for j in range(n):
            if i < j:
                s[(i,j)] = random.uniform(-1, 1)
                s[(j,i)] = -s[(i,j)]
    return s

def E_v_numeric(vertex_set, target, r_val, s):
    vs = list(vertex_set)
    if len(vs) == 1: return 1.0
    total = 0.0
    for perm in permutations([v for v in vs if v != target]):
        path = list(perm) + [target]
        w = 1.0
        for k in range(len(path)-1):
            u, v = path[k], path[k+1]
            w *= (r_val + s[(u,v)])
        total += w
    return total

def G_numeric(n, a, b, r_val, s):
    U = [v for v in range(n) if v != a and v != b]
    total = 0.0
    for mask in range(2**len(U)):
        S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
        R = [u for u in U if u not in S]
        f = E_v_numeric(set(S + [a]), a, r_val, s)
        g = E_v_numeric(set(R + [b]), b, -r_val, s)
        total += f * g
    return total

for n in [5, 6, 7, 8]:
    s_num = make_numeric_tournament(n)
    for trial in range(5):
        r_val = random.uniform(0.1, 2.0)
        G_pos = G_numeric(n, 0, 1, r_val, s_num)
        G_neg = G_numeric(n, 0, 1, -r_val, s_num)
        diff = abs(G_pos - G_neg)
        print(f"  n={n}, r={r_val:.3f}: G(r)={G_pos:.8f}, G(-r)={G_neg:.8f}, |diff|={diff:.2e}")

# ================================================================
# KEY INSIGHT: Express G as sum of PAIRED products
# ================================================================
print("\n" + "=" * 70)
print("PAIRED PRODUCT STRUCTURE")
print("=" * 70)

# For each edge set, define:
# P(eps; r) = prod_e (r + eps_e * s_e) where eps is a sign vector
# Then G(r) involves sums of P(eps; r) with different eps vectors.

# KEY: P(eps; r) + P(eps; -r) = P(eps; r) + (-1)^m P(-eps; r)
# So for m even: P(eps; r) + P(eps; -r) = P(eps; r) + P(-eps; r)
#                                        = prod(r+eps*s) + prod(r-eps*s)
# This is manifestly even in r (it's 2 * even part of P(eps; r)).

# For m odd: P(eps; r) + P(eps; -r) = P(eps; r) - P(-eps; r)
#                                    = prod(r+eps*s) - prod(r-eps*s)
# This is 2 * odd part = manifestly ODD in r.

# But G involves a SUM over different eps vectors (from different path pairs).
# The question is whether the sum remains even.

# OBSERVATION at n=4 (m=2, even):
# Each pair contributes P(eps; r) = (r+s_1)(s_2-r) or similar.
# The sum over all pairs somehow gives an even polynomial.

# At n=5 (m=3, odd): each pair contributes prod of 3 factors, odd degree.
# But G still has only even powers!

# This means: the CROSS TERMS between even-m and odd-m contributions
# must cancel the wrong-parity parts.

# Wait, m = n-2 is fixed! All terms have the same m.
# At n=4: m=2 (even). Each single product P(eps;r) has degree 2.
# P(eps;r) + P(-eps;r) is even. So P(eps;r) = P_even + P_odd,
# and summing PAIRS (eps, -eps) gives 2*P_even.

# But G doesn't pair eps with -eps! G pairs (f_S * g_S(-r)),
# which uses eps = (+...+, -...-) where + for Pa edges, - for Pb edges.

# Actually, for a FIXED path pair (Pa, Pb):
# Product = prod_{e in Pa}(r+s_e) * prod_{e in Pb}(s_e - r)
# = prod_{e in Pa}(r+s_e) * (-1)^|Pb| * prod_{e in Pb}(r - s_e)
# = (-1)^|R| * prod over all edges (r + eps_e * s_e)
# where eps_e = +1 for Pa edges, -1 for Pb edges.

# Now: under r -> -r:
# prod(r + eps_e s_e) -> prod(-r + eps_e s_e) = (-1)^m prod(r - eps_e s_e)

# So the product at -r = (-1)^m * product with flipped eps.

# For the SUM: sum over (S, paths) of [product at r]
# equals sum over (S, paths) of [(-1)^m * product at -r with flipped eps].

# "Flipped eps" means: edges that were in Pa are now in Pb and vice versa.
# This corresponds to: swapping Pa and Pb roles!

# But swapping Pa and Pb means: what was E_a through S+a becomes
# something through the SAME S+a but as a different path type.

# WAIT: flipping eps changes (r+s) to (r-s) for Pa edges and
# (s-r) to (-s+r) = (r-s) for Pb edges.
# Hmm wait: Pa edge has eps=+1, so (r+s). Flipped: (r-s).
# Pb edge has eps=-1, so (r-s). Flipped: (r+s).

# So flipping eps swaps the role of Pa and Pb edges.
# The product with flipped eps = prod_{Pa edges}(r-s_e) * prod_{Pb edges}(r+s_e)
# = prod_{Pa edges}(r-s_e) * prod_{Pb edges}(r+s_e)

# This is the product for the REVERSED tournament (s -> -s)!
# For a c-tournament: s -> -s means T -> T^op.

# So: G(-r) = (-1)^m * G^{op}(r) where G^{op} uses -s.

# But we proved M(-r, s) = M(b,a)(r, s) = (-1)^{n-2} M(a,b)(r, -s)... hmm, we're going in circles.

# Let me just verify the key structural fact:

print("\nn=4: Detailed eps-vector analysis")
n = 4
r, s, t = make_tournament(n)
a, b = 0, 1
U = [2, 3]

for mask in range(4):
    S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
    R = [u for u in U if u not in S]
    Sa = set(S + [a])
    Rb = set(R + [b])

    # All path pairs
    paths_a = []
    for p in permutations(sorted(Sa)):
        if p[-1] == a:
            edges_a = [(p[k], p[k+1]) for k in range(len(p)-1)]
            paths_a.append(edges_a)
    if len(Sa) == 1: paths_a = [[]]

    paths_b = []
    for p in permutations(sorted(Rb)):
        if p[-1] == b:  # Note: E_b means ending at b
            edges_b = [(p[k], p[k+1]) for k in range(len(p)-1)]
            paths_b.append(edges_b)
    if len(Rb) == 1: paths_b = [[]]

    print(f"\n  S={set(S)}, Sa={Sa}, Rb={Rb}")
    for pa in paths_a:
        for pb in paths_b:
            eps_str = []
            product = 1
            for e in pa:
                eps_str.append(f"+s{min(e)}{max(e)}")
                product *= (r + s[(e[0], e[1])])
            for e in pb:
                eps_str.append(f"-s{min(e)}{max(e)}")
                product *= (s[(e[0], e[1])] - r)
            product = expand(product)
            print(f"    Pa={pa}, Pb={pb}")
            print(f"    eps: {', '.join(eps_str)}")
            print(f"    prod = {product}")

            # Even/odd in r
            p_poly = Poly(product, r)
            p_even = sum(expand(p_poly.nth(k)) * r**k for k in range(0, p_poly.degree()+1, 2))
            p_odd = sum(expand(p_poly.nth(k)) * r**k for k in range(1, p_poly.degree()+1, 2))
            print(f"    even part: {expand(p_even)}, odd part: {expand(p_odd)}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
