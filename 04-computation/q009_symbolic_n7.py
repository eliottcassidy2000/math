#!/usr/bin/env python3
"""
Symbolic verification of alpha_w^H = alpha_w^I at n=7.

At n=7: W = {w, u1, u2, u3, u4, u5}, B_w has 5 vertices.
5! = 120 permutations, each with 6 insertion positions = 720 terms.
Plus cycle derivatives from chains of length 3 and 5.

Instance: opus-2026-03-05-S5
"""

from sympy import symbols, expand, simplify, S as Sym, Rational
from itertools import permutations, combinations

def verify_n7():
    print("=== Symbolic verification at n=7 ===\n")

    # B_w vertices: a, b, c, d  (at n=7: |W|=5, |B_w|=4)
    B_w = ['a', 'b', 'c', 'd']
    m_minus_1 = 4  # |B_w|

    # Arc variables
    # From w to B_w: w_a, w_b, w_c, w_d
    w_a, w_b, w_c, w_d = symbols('w_a w_b w_c w_d')
    arc_from_w = {'a': w_a, 'b': w_b, 'c': w_c, 'd': w_d}

    # Internal arcs among B_w: C(4,2) = 6 arcs
    t_ab, t_ac, t_ad = symbols('t_ab t_ac t_ad')
    t_bc, t_bd = symbols('t_bc t_bd')
    t_cd = symbols('t_cd')

    internal = {
        ('a','b'): t_ab, ('b','a'): 1-t_ab,
        ('a','c'): t_ac, ('c','a'): 1-t_ac,
        ('a','d'): t_ad, ('d','a'): 1-t_ad,
        ('b','c'): t_bc, ('c','b'): 1-t_bc,
        ('b','d'): t_bd, ('d','b'): 1-t_bd,
        ('c','d'): t_cd, ('d','c'): 1-t_cd,
    }

    # Interface arcs: p_a, p_b, ... (q = 1-p)
    p_a, p_b, p_c, p_d, p_w = symbols('p_a p_b p_c p_d p_w')
    p_vals = {'a': p_a, 'b': p_b, 'c': p_c, 'd': p_d, 'w': p_w}
    q_vals = {'a': 1-p_a, 'b': 1-p_b, 'c': 1-p_c, 'd': 1-p_d, 'w': 1-p_w}

    def T(s, t):
        if s == t: return Sym(0)
        if s == 'w':
            return arc_from_w[t]
        if t == 'w':
            return 1 - arc_from_w[s]
        return internal.get((s,t), 1 - internal.get((t,s), Sym(0)))

    # --- Compute alpha_w^H ---
    print("Computing alpha_w^H (this may take a while)...")
    alpha_H = Sym(0)
    count = 0
    for pi in permutations(B_w):
        B_wt = Sym(1)
        for k in range(len(pi) - 1):
            B_wt *= T(pi[k], pi[k+1])

        for pos in range(m_minus_1 + 1):
            if pos == 0:
                val = -(T('w', pi[0]) + q_vals[pi[0]]) * B_wt
            elif pos == m_minus_1:
                val = -(T(pi[-1], 'w') + p_vals[pi[-1]]) * B_wt
            else:
                u_k = pi[pos-1]
                u_k1 = pi[pos]
                numer = T(u_k, 'w') * q_vals[u_k1] + p_vals[u_k] * T('w', u_k1)
                t_arc = T(u_k, u_k1)
                val = -numer * B_wt / t_arc
            alpha_H += val

        count += 1
        if count % 30 == 0:
            print(f"  ... {count}/120 permutations processed")

    print("Expanding alpha_H...")
    alpha_H = expand(alpha_H)
    print(f"  alpha_H has {len(alpha_H.as_ordered_terms())} terms")

    # --- Compute alpha_w^I ---
    print("Computing alpha_w^I...")

    # -2*H(B_w)
    H_Bw = Sym(0)
    for pi in permutations(B_w):
        wt = Sym(1)
        for k in range(len(pi)-1):
            wt *= T(pi[k], pi[k+1])
        H_Bw += wt
    alpha_I = -2 * H_Bw

    # 5-cycle derivatives: chains of length 3 with w at endpoint
    W_verts = ['w'] + B_w
    for combo in combinations(W_verts, 3):
        if 'w' not in combo:
            continue
        for perm in permutations(combo):
            chain = perm
            int_wt = Sym(1)
            for k in range(len(chain)-1):
                int_wt *= T(chain[k], chain[k+1])
            if chain[0] == 'w':
                alpha_I += 2 * (-p_vals[chain[-1]]) * int_wt
            elif chain[-1] == 'w':
                alpha_I += 2 * (-q_vals[chain[0]]) * int_wt

    # 7-cycle derivatives: chains of length 5 with w at endpoint
    for combo in combinations(W_verts, 5):
        if 'w' not in combo:
            continue
        for perm in permutations(combo):
            chain = perm
            int_wt = Sym(1)
            for k in range(len(chain)-1):
                int_wt *= T(chain[k], chain[k+1])
            if chain[0] == 'w':
                alpha_I += 2 * (-p_vals[chain[-1]]) * int_wt
            elif chain[-1] == 'w':
                alpha_I += 2 * (-q_vals[chain[0]]) * int_wt

    print("Expanding alpha_I...")
    alpha_I = expand(alpha_I)
    print(f"  alpha_I has {len(alpha_I.as_ordered_terms())} terms")

    # --- Check difference ---
    print("Computing difference...")
    diff = expand(alpha_H - alpha_I)
    n_terms = len(diff.as_ordered_terms()) if diff != 0 else 0
    print(f"  diff has {n_terms} terms")

    if diff == 0:
        print("\n*** PROVED: alpha_w^H = alpha_w^I at n=7 ***\n")
        return True
    else:
        print(f"  Trying simplify...")
        diff_s = simplify(diff)
        if diff_s == 0:
            print("\n*** PROVED: alpha_w^H = alpha_w^I at n=7 ***\n")
            return True
        else:
            # Numerical spot check
            test_subs = {
                w_a: Rational(3,10), w_b: Rational(7,10), w_c: Rational(1,2),
                w_d: Rational(2,5),
                t_ab: Rational(2,5), t_ac: Rational(3,5), t_ad: Rational(1,3),
                t_bc: Rational(4,5), t_bd: Rational(1,4),
                t_cd: Rational(1,5),
                p_a: Rational(1,3), p_b: Rational(2,3), p_c: Rational(1,4),
                p_d: Rational(3,4), p_w: Rational(2,5)
            }
            val = diff.subs(test_subs)
            print(f"  Numerical check: {val}")
            if val == 0:
                print("  (Numerically zero; symbolic simplification may need more work)")
            return False


if __name__ == "__main__":
    verify_n7()
