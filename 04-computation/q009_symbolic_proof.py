#!/usr/bin/env python3
"""
Symbolic verification that alpha_w^H = alpha_w^I at n=5 (and attempt for general n).

At n=5: W = {w, u, v}, B_w = {u, v}.

alpha_w^H = -sum_pi B_wt * S(pi)  where S involves inserting w into perms of B_w.
alpha_w^I = -2*H(B_w) + cycle_derivs(w)

We need to show these are equal symbolically.

At s=0: p_x + q_x = 1 for all x in W, T(I,x) = T(J,x) = q_x, T(x,I) = T(x,J) = p_x.

Instance: opus-2026-03-05-S5
"""

try:
    from sympy import symbols, simplify, expand, factor, collect, Rational, S as Sym
    HAS_SYMPY = True
except ImportError:
    HAS_SYMPY = False
    print("sympy not available, using numerical verification only")

from itertools import permutations
import random


def verify_n5_symbolic():
    """Verify alpha_w^H = alpha_w^I symbolically at n=5."""
    if not HAS_SYMPY:
        return

    print("=== Symbolic verification at n=5 ===\n")

    # Variables: tournament arcs
    # T(u,v) = t_uv, T(v,u) = 1-t_uv
    # T(w,u) = a, T(w,v) = b  (internal arcs involving w)
    # T(u,v) = t  (internal arc not involving w)
    # p_u, p_v, p_w = T(u,I), T(v,I), T(w,I) (interface arcs)
    # q_u = 1-p_u, q_v = 1-p_v, q_w = 1-p_w (at s=0)

    a, b, t = symbols('a b t', positive=True)
    p_u, p_v, p_w = symbols('p_u p_v p_w', positive=True)

    # Derived
    q_u = 1 - p_u
    q_v = 1 - p_v
    q_w = 1 - p_w

    # alpha_w^H from insertion decomposition
    # Boundary start (pos 0): -(T(w,u_1) + q_{u_1}) * B_wt
    # Boundary end (pos m-1): -(T(u_{m-1},w) + p_{u_{m-1}}) * B_wt
    # Interior (pos 1): -(T(u_1,w)*q_{u_2} + p_{u_1}*T(w,u_2)) (already divided by T(u_1,u_2) and multiplied by B_wt)

    # pi = (u, v): B_wt = t
    # pos 0: -(a + q_u) * t
    # pos 2: -((1-a) + p_v) * t (T(v,w)=1-b → wait, T(u_{m-1},w) with u_{m-1}=v, so T(v,w)=1-b)
    #         Actually: -(T(v,w) + p_v) * t = -(1-b + p_v) * t
    # pos 1 (interior): -(T(u,w)*q_v + p_u*T(w,v)) = -((1-a)*q_v + p_u*b)

    alpha_H_uv = -(a + q_u) * t - (1 - b + p_v) * t - ((1 - a) * q_v + p_u * b)

    # pi = (v, u): B_wt = 1-t
    # pos 0: -(b + q_v) * (1-t)
    # pos 2: -((1-a) + p_u) * (1-t) ... wait, T(u_{m-1},w) with u_{m-1}=u, so T(u,w)=1-a
    #         -(1-a + p_u) * (1-t)
    # pos 1 (interior): -(T(v,w)*q_u + p_v*T(w,u)) = -((1-b)*q_u + p_v*a)

    alpha_H_vu = -(b + q_v) * (1 - t) - (1 - a + p_u) * (1 - t) - ((1 - b) * q_u + p_v * a)

    alpha_H = expand(alpha_H_uv + alpha_H_vu)
    print(f"alpha_H = {alpha_H}")

    # alpha_w^I
    # -2*H(B_w) = -2*(t + (1-t)) = -2
    alpha_I_3 = -2

    # Cycle derivatives for 5-cycles:
    # Chains of length 3 through W = {w,u,v} with w at endpoint:
    # (w,u,v): 2*(-p_v)*T(w,u)*T(u,v) = -2*p_v*a*t
    # (w,v,u): 2*(-p_u)*T(w,v)*T(v,u) = -2*p_u*b*(1-t)
    # (u,v,w): 2*(-q_u)*T(u,v)*T(v,w) = -2*q_u*t*(1-b)
    # (v,u,w): 2*(-q_v)*T(v,u)*T(u,w) = -2*q_v*(1-t)*(1-a)

    alpha_I_5 = -2 * p_v * a * t - 2 * p_u * b * (1 - t) - 2 * q_u * t * (1 - b) - 2 * q_v * (1 - t) * (1 - a)

    alpha_I = expand(alpha_I_3 + alpha_I_5)
    print(f"alpha_I = {alpha_I}")

    # Check difference
    diff = expand(alpha_H - alpha_I)
    print(f"\nalpha_H - alpha_I = {diff}")
    print(f"Simplified: {simplify(diff)}")

    if simplify(diff) == 0:
        print("\n*** PROVED: alpha_w^H = alpha_w^I at n=5 ***\n")
    else:
        print(f"\nNon-zero difference: {diff}")
        # Try collecting terms
        print(f"Collected by a: {collect(diff, a)}")
        print(f"Collected by t: {collect(diff, t)}")


def verify_n6_symbolic():
    """Verify alpha_w^H = alpha_w^I symbolically at n=6."""
    if not HAS_SYMPY:
        return

    print("=== Symbolic verification at n=6 ===\n")

    # W = {w, u, v, x}, B_w = {u, v, x}
    # Internal arcs: T(w,u)=a, T(w,v)=b, T(w,x)=c (arcs from w)
    # T(u,v)=t_uv, T(u,x)=t_ux, T(v,x)=t_vx (arcs among B_w)
    # Interface: p_u, p_v, p_x, p_w (and q = 1-p)

    a, b, c = symbols('a b c', positive=True)
    t_uv, t_ux, t_vx = symbols('t_uv t_ux t_vx', positive=True)
    p_u, p_v, p_x, p_w = symbols('p_u p_v p_x p_w', positive=True)

    q_u, q_v, q_x, q_w = 1-p_u, 1-p_v, 1-p_x, 1-p_w

    # T function
    arcs = {
        ('w','u'): a, ('u','w'): 1-a,
        ('w','v'): b, ('v','w'): 1-b,
        ('w','x'): c, ('x','w'): 1-c,
        ('u','v'): t_uv, ('v','u'): 1-t_uv,
        ('u','x'): t_ux, ('x','u'): 1-t_ux,
        ('v','x'): t_vx, ('x','v'): 1-t_vx,
    }
    p_vals = {'u': p_u, 'v': p_v, 'x': p_x, 'w': p_w}
    q_vals = {'u': q_u, 'v': q_v, 'x': q_x, 'w': q_w}

    def T(s, t):
        if s == t: return 0
        return arcs.get((s,t), 1 - arcs.get((t,s), Sym(0)))

    B_w = ['u', 'v', 'x']

    # Compute alpha_w^H via insertion decomposition
    alpha_H = Sym(0)
    for pi in permutations(B_w):
        # B_wt = product of arcs in pi
        B_wt = Sym(1)
        for k in range(len(pi) - 1):
            B_wt *= T(pi[k], pi[k+1])

        # Insert w at each position
        m_minus_1 = len(pi)  # 3
        for pos in range(m_minus_1 + 1):  # 0,1,2,3
            if pos == 0:
                # -(T(w,pi[0]) + q_{pi[0]}) * B_wt
                factor_val = -(T('w', pi[0]) + q_vals[pi[0]]) * B_wt
            elif pos == m_minus_1:
                # -(T(pi[-1],w) + p_{pi[-1]}) * B_wt
                factor_val = -(T(pi[-1], 'w') + p_vals[pi[-1]]) * B_wt
            else:
                # Interior: -(T(pi[pos-1],w)*q_{pi[pos]} + p_{pi[pos-1]}*T(w,pi[pos]))
                #   * B_wt / T(pi[pos-1], pi[pos]) ... but B_wt includes T(pi[pos-1],pi[pos])
                # Wait: the formula is:
                # -(T(u_k,w)*q_{u_{k+1}} + p_{u_k}*T(w,u_{k+1})) * L_wt * R_wt
                # where L_wt * R_wt = B_wt / T(u_k, u_{k+1})
                u_k = pi[pos-1]
                u_k1 = pi[pos]
                numer = T(u_k, 'w') * q_vals[u_k1] + p_vals[u_k] * T('w', u_k1)
                t_arc = T(u_k, u_k1)
                # L_wt * R_wt = B_wt / t_arc, but we need to be careful with division
                factor_val = -numer * B_wt / t_arc

            alpha_H += factor_val

    alpha_H = expand(alpha_H)

    # Compute alpha_w^I
    # -2*H(B_w)
    H_Bw = Sym(0)
    for pi in permutations(B_w):
        wt = Sym(1)
        for k in range(len(pi)-1):
            wt *= T(pi[k], pi[k+1])
        H_Bw += wt

    alpha_I = -2 * H_Bw

    # 5-cycle derivatives: chains of length 3 from W={w,u,v,x} with w at endpoint
    for combo in [('w','u','v'), ('w','u','x'), ('w','v','u'), ('w','v','x'),
                  ('w','x','u'), ('w','x','v'),
                  ('u','v','w'), ('u','x','w'), ('v','u','w'), ('v','x','w'),
                  ('x','u','w'), ('x','v','w')]:
        chain = combo
        internal = Sym(1)
        for k in range(len(chain)-1):
            internal *= T(chain[k], chain[k+1])
        if chain[0] == 'w':
            alpha_I += 2 * (-p_vals[chain[-1]]) * internal
        elif chain[-1] == 'w':
            alpha_I += 2 * (-q_vals[chain[0]]) * internal

    alpha_I = expand(alpha_I)

    diff = expand(alpha_H - alpha_I)
    diff_simplified = simplify(diff)

    print(f"alpha_H - alpha_I simplified: {diff_simplified}")

    if diff_simplified == 0:
        print("\n*** PROVED: alpha_w^H = alpha_w^I at n=6 ***\n")
    else:
        print(f"\nNon-zero difference, trying harder...")
        # Try substituting specific values to check
        test_subs = {a: Rational(3,10), b: Rational(7,10), c: Rational(1,2),
                     t_uv: Rational(2,5), t_ux: Rational(3,5), t_vx: Rational(4,5),
                     p_u: Rational(1,3), p_v: Rational(2,3), p_x: Rational(1,4), p_w: Rational(1,2)}
        val = diff.subs(test_subs)
        print(f"  Numerical check: {val}")


def verify_n4_symbolic():
    """Quick check at n=4."""
    if not HAS_SYMPY:
        return

    print("=== Symbolic verification at n=4 ===\n")

    # W = {w, u}, B_w = {u}
    # T(w,u) = a, p_u = T(u,I), q_u = 1-p_u
    a = symbols('a', positive=True)
    p_u = symbols('p_u', positive=True)
    q_u = 1 - p_u

    # alpha_H: pi = (u,), B_wt = 1
    # S(pi) = T(w,u)+q_u + T(u,w)+p_u = a+(1-p_u)+(1-a)+p_u = 2
    # alpha_H = -1 * 2 = -2

    alpha_H = -(a + q_u) - ((1-a) + p_u)
    alpha_H = expand(alpha_H)
    print(f"alpha_H = {alpha_H}")

    # alpha_I = -2*H(B_w) = -2*1 = -2 (no 5-cycles at n=4)
    alpha_I = -2
    print(f"alpha_I = {alpha_I}")

    diff = simplify(alpha_H - alpha_I)
    print(f"diff = {diff}")
    if diff == 0:
        print("*** PROVED at n=4 ***\n")


if __name__ == "__main__":
    if HAS_SYMPY:
        verify_n4_symbolic()
        verify_n5_symbolic()
        verify_n6_symbolic()
    else:
        print("Install sympy: pip install sympy")
