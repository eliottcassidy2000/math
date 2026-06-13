#!/usr/bin/env python3
"""
Half-symmetry test: Is F_u^{[a,b]}(r) = F_u^{[b,a]}(r)?

If yes, then M[a,b] = F_u(r) + G_u(r) where G_u(r) = F_u(-r),
making M manifestly even in r!

Proof sketch (if the half-symmetry holds):
  From path reversal: F_u(-r) = G_u^{[b,a]}(r) and G_u(-r) = F_u^{[b,a]}(r)
  If F_u^{[a,b]} = F_u^{[b,a]}, then G_u(-r) = F_u(r), so G_u(r) = F_u(-r).
  Then M[a,b] = F_u(r) + F_u(-r) = manifestly even in r.
  Also M[a,b] = M[b,a] automatically.

Instance: opus-2026-03-06-S22
"""

from itertools import permutations
from sympy import symbols, expand, Symbol

def make_symbols(n):
    r = Symbol('r')
    s = {}
    for i in range(n):
        for j in range(n):
            if i < j:
                s[(i,j)] = Symbol(f's{i}{j}')
                s[(j,i)] = -s[(i,j)]
    return r, s

def edge_weight(i, j, r, s):
    return r + s[(i,j)]

def ham_paths_ending_at(vertex_set, target, r, s):
    vs = list(vertex_set)
    if len(vs) == 1:
        return 1
    total = 0
    others = [v for v in vs if v != target]
    for perm in permutations(others):
        path = list(perm) + [target]
        w = 1
        for k in range(len(path)-1):
            w *= edge_weight(path[k], path[k+1], r, s)
        total += w
    return expand(total)

def ham_paths_beginning_at(vertex_set, source, r, s):
    vs = list(vertex_set)
    if len(vs) == 1:
        return 1
    total = 0
    others = [v for v in vs if v != source]
    for perm in permutations(others):
        path = [source] + list(perm)
        w = 1
        for k in range(len(path)-1):
            w *= edge_weight(path[k], path[k+1], r, s)
        total += w
    return expand(total)


def compute_halves(a, b, u, n, r, s):
    """Compute F_u and G_u for M[a,b]."""
    U = [v for v in range(n) if v != a and v != b]
    U_minus_u = [v for v in U if v != u]

    F_u = 0  # sum over S with u ∈ S
    G_u = 0  # sum over S with u ∉ S

    for mask in range(2**len(U_minus_u)):
        Sp = [U_minus_u[i] for i in range(len(U_minus_u)) if (mask >> i) & 1]
        Rp = [v for v in U_minus_u if v not in Sp]

        # u in S case: S = Sp ∪ {u}
        S_full = Sp + [u]
        R_full = Rp
        sign = (-1)**len(S_full)
        Ea = ham_paths_ending_at(set(S_full + [a]), a, r, s)
        Bb = ham_paths_beginning_at(set(R_full + [b]), b, r, s)
        F_u += sign * Ea * Bb

        # u not in S case: S = Sp
        S_full2 = Sp
        R_full2 = Rp + [u]
        sign2 = (-1)**len(S_full2)
        Ea2 = ham_paths_ending_at(set(S_full2 + [a]), a, r, s)
        Bb2 = ham_paths_beginning_at(set(R_full2 + [b]), b, r, s)
        G_u += sign2 * Ea2 * Bb2

    return expand(F_u), expand(G_u)


def test_half_symmetry(n):
    """Test if F_u^{[a,b]} = F_u^{[b,a]} for all valid (a,b,u)."""
    print(f"\n{'='*60}")
    print(f"HALF-SYMMETRY TEST: n={n}")
    print(f"{'='*60}")

    r, s = make_symbols(n)

    # Test multiple (a,b) pairs
    pairs = [(0, 1), (0, 2), (1, 2)]
    if n >= 5:
        pairs = [(0, 1)]  # restrict for speed

    for a, b in pairs:
        U = [v for v in range(n) if v != a and v != b]
        for u in U:
            F_u_ab, G_u_ab = compute_halves(a, b, u, n, r, s)
            F_u_ba, G_u_ba = compute_halves(b, a, u, n, r, s)

            # Test 1: F_u^{[a,b]} = F_u^{[b,a]}?
            diff_F = expand(F_u_ab - F_u_ba)
            half_sym = (diff_F == 0)

            # Test 2: G_u(r) = F_u(-r)?
            F_u_neg = expand(F_u_ab.subs(r, -r))
            diff_GF = expand(G_u_ab - F_u_neg)
            gu_is_fu_neg = (diff_GF == 0)

            # Test 3: M[a,b] = F_u(r) + F_u(-r)?
            M_ab = expand(F_u_ab + G_u_ab)
            M_even = expand(F_u_ab + F_u_neg)
            m_is_even = expand(M_ab - M_even) == 0

            print(f"\n  a={a}, b={b}, u={u}:")
            print(f"    F_u^[a,b] = F_u^[b,a]? {half_sym}")
            if not half_sym:
                print(f"    F_u^[a,b] - F_u^[b,a] = {diff_F}")
            print(f"    G_u(r) = F_u(-r)?       {gu_is_fu_neg}")
            if not gu_is_fu_neg:
                print(f"    G_u - F_u(-r) = {diff_GF}")
            print(f"    M = F_u(r) + F_u(-r)?   {m_is_even}")

            # If half-symmetry fails, check weaker: F_u^[a,b](r) = G_u^[b,a](-r)?
            if not half_sym:
                G_u_ba_neg = expand(G_u_ba.subs(r, -r))
                weaker = expand(F_u_ab - G_u_ba_neg) == 0
                print(f"    F_u^[a,b](r) = G_u^[b,a](-r)? {weaker} (from path reversal)")

    # Also compute and display F_u and G_u explicitly for small n
    if n == 4:
        a, b, u = 0, 1, 2
        r_sym, s_sym = make_symbols(n)
        F, G = compute_halves(a, b, u, n, r_sym, s_sym)
        print(f"\n  EXPLICIT (n=4, a=0, b=1, u=2):")
        print(f"    F_u = {F}")
        print(f"    G_u = {G}")
        print(f"    F_u(-r) = {expand(F.subs(r_sym, -r_sym))}")


def main():
    for n in [3, 4, 5]:
        test_half_symmetry(n)

    print(f"\n\n{'='*60}")
    print("THEORETICAL IMPLICATIONS")
    print(f"{'='*60}")
    print("""
If F_u^{[a,b]}(r) = F_u^{[b,a]}(r) for all u ∈ U:

PROOF OF EVEN R-POWERS:
  1. Path reversal gives: G_u(-r) = F_u^{[b,a]}(r) = F_u^{[a,b]}(r) = F_u(r)
  2. So G_u(r) = F_u(-r)
  3. Therefore M[a,b](r) = F_u(r) + G_u(r) = F_u(r) + F_u(-r)
  4. This is manifestly EVEN in r (symmetric under r → -r)

PROOF OF SYMMETRY:
  5. M[b,a](r) = F_u^{[b,a]}(r) + G_u^{[b,a]}(r) = F_u(r) + F_u(-r) = M[a,b](r)

So BOTH conjectures follow from the single claim:
  "The u∈S half of M[a,b] equals the u∈S half of M[b,a]."
""")


if __name__ == '__main__':
    main()
