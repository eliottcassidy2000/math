#!/usr/bin/env python3
"""
Symbolic verification of alpha_w^H = alpha_w^I at n=8.

At n=8: |W|=6, |B_w|=5. 5!=120 perms, 6 insertion positions = 720 terms.
Cycles of odd length through {i,j}: 3,5,7 (chains of length 1,3,5).
Max chain length = 5 (using 5 of 6 W-vertices).

Instance: opus-2026-03-05-S5
"""

from sympy import symbols, expand, simplify, S as Sym, Rational
from itertools import permutations, combinations

def verify_n8():
    print("=== Symbolic verification at n=8 ===\n")

    # B_w has 5 vertices at n=8
    B_w = ['a', 'b', 'c', 'd', 'e']
    m_minus_1 = 5

    # Arc variables from w
    w_a, w_b, w_c, w_d, w_e = symbols('w_a w_b w_c w_d w_e')
    arc_from_w = {'a': w_a, 'b': w_b, 'c': w_c, 'd': w_d, 'e': w_e}

    # Internal arcs: C(5,2) = 10
    t_ab, t_ac, t_ad, t_ae = symbols('t_ab t_ac t_ad t_ae')
    t_bc, t_bd, t_be = symbols('t_bc t_bd t_be')
    t_cd, t_ce = symbols('t_cd t_ce')
    t_de = symbols('t_de')

    internal = {}
    pairs = [('a','b',t_ab), ('a','c',t_ac), ('a','d',t_ad), ('a','e',t_ae),
             ('b','c',t_bc), ('b','d',t_bd), ('b','e',t_be),
             ('c','d',t_cd), ('c','e',t_ce), ('d','e',t_de)]
    for x, y, t in pairs:
        internal[(x,y)] = t
        internal[(y,x)] = 1-t

    # Interface arcs
    p_a, p_b, p_c, p_d, p_e, p_w = symbols('p_a p_b p_c p_d p_e p_w')
    p_vals = {'a': p_a, 'b': p_b, 'c': p_c, 'd': p_d, 'e': p_e, 'w': p_w}
    q_vals = {k: 1-v for k,v in p_vals.items()}

    def T(s, t):
        if s == t: return Sym(0)
        if s == 'w': return arc_from_w[t]
        if t == 'w': return 1 - arc_from_w[s]
        return internal.get((s,t), 1 - internal.get((t,s), Sym(0)))

    # --- Compute alpha_w^H ---
    print("Computing alpha_w^H...")
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
        if count % 40 == 0:
            print(f"  ... {count}/120 permutations")

    print("Expanding alpha_H...")
    alpha_H = expand(alpha_H)
    n_terms_H = len(alpha_H.as_ordered_terms())
    print(f"  {n_terms_H} terms")

    # --- Compute alpha_w^I ---
    print("Computing alpha_w^I...")
    W_verts = ['w'] + B_w

    # -2*H(B_w)
    H_Bw = Sym(0)
    for pi in permutations(B_w):
        wt = Sym(1)
        for k in range(len(pi)-1):
            wt *= T(pi[k], pi[k+1])
        H_Bw += wt
    alpha_I = -2 * H_Bw

    # FULL inductive formula: chains of odd length >= 3 with w at endpoint,
    # weighted by H(comp) where comp = W \ {chain vertices}
    for chain_len in [3, 5]:  # max chain = 5 (uses 5 of 6 W-vertices)
        for combo in combinations(W_verts, chain_len):
            if 'w' not in combo:
                continue
            comp = [x for x in W_verts if x not in combo]
            # H(comp) = sum of Ham path weights in comp
            H_comp = Sym(0)
            if len(comp) == 0:
                H_comp = Sym(1)
            elif len(comp) == 1:
                H_comp = Sym(1)
            else:
                for pi in permutations(comp):
                    wt = Sym(1)
                    for k in range(len(pi)-1):
                        wt *= T(pi[k], pi[k+1])
                    H_comp += wt

            for perm in permutations(combo):
                int_wt = Sym(1)
                for k in range(len(perm)-1):
                    int_wt *= T(perm[k], perm[k+1])
                if perm[0] == 'w':
                    alpha_I += 2 * (-p_vals[perm[-1]]) * int_wt * H_comp
                elif perm[-1] == 'w':
                    alpha_I += 2 * (-q_vals[perm[0]]) * int_wt * H_comp

    print("Expanding alpha_I...")
    alpha_I = expand(alpha_I)
    n_terms_I = len(alpha_I.as_ordered_terms())
    print(f"  {n_terms_I} terms")

    # --- Check ---
    print("Computing difference...")
    diff = expand(alpha_H - alpha_I)
    n_terms = len(diff.as_ordered_terms()) if diff != 0 else 0
    print(f"  diff has {n_terms} terms")

    if diff == 0:
        print("\n*** PROVED: alpha_w^H = alpha_w^I at n=8 ***\n")
        return True
    else:
        print("  Trying simplify...")
        diff_s = simplify(diff)
        if diff_s == 0:
            print("\n*** PROVED (after simplify): alpha_w^H = alpha_w^I at n=8 ***\n")
            return True
        else:
            # Numerical spot check
            test_subs = {
                w_a: Rational(3,10), w_b: Rational(7,10), w_c: Rational(1,2),
                w_d: Rational(2,5), w_e: Rational(3,5),
                t_ab: Rational(2,5), t_ac: Rational(3,5), t_ad: Rational(1,3),
                t_ae: Rational(2,3), t_bc: Rational(4,5), t_bd: Rational(1,4),
                t_be: Rational(3,4), t_cd: Rational(1,5), t_ce: Rational(7,10),
                t_de: Rational(9,10),
                p_a: Rational(1,3), p_b: Rational(2,3), p_c: Rational(1,4),
                p_d: Rational(3,4), p_e: Rational(1,2), p_w: Rational(2,5)
            }
            val = diff.subs(test_subs)
            print(f"  Numerical: {val}")
            return val == 0


if __name__ == "__main__":
    verify_n8()
