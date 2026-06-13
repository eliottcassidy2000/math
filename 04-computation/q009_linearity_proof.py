#!/usr/bin/env python3
"""
KEY INSIGHT: The identity P = interior + fE + chains = 0 is LINEAR in each
w-arc variable T(w,v). This means P = c_0 + sum_v c_v * T(w,v) where c_0
and c_v are polynomials in internal arcs and interface arcs only.

P = 0 iff c_0 = 0 AND c_v = 0 for each v in B_w.

This reduces the proof from O(m^2) variables to m separate polynomial
identities in O(m^2/2 + m) variables.

This script:
1. Verifies the linearity in w-arcs (no products T(w,u)*T(w,v))
2. Computes c_0 and c_v explicitly
3. Verifies c_0 = 0 and c_v = 0

Instance: opus-2026-03-05-S6
"""

from itertools import permutations, combinations
import random

try:
    from sympy import symbols, expand, S as Sym
    HAS_SYMPY = True
except ImportError:
    HAS_SYMPY = False


def make_T(n, arc_values):
    def T(a, b):
        if a == b: return 0.0
        return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))
    return T


def random_s0_tournament(n, seed=None):
    if seed is not None:
        random.seed(seed)
    I, J = 0, 1
    W = list(range(2, n))
    arc_values = {(I, J): 1.0, (J, I): 0.0}
    for w in W:
        q_w = random.uniform(0.1, 0.9)
        arc_values[(J, w)] = q_w
        arc_values[(w, J)] = 1 - q_w
        arc_values[(w, I)] = 1 - q_w
        arc_values[(I, w)] = q_w
    for a in W:
        for b in W:
            if a < b:
                v = random.uniform(0.1, 0.9)
                arc_values[(a, b)] = v
                arc_values[(b, a)] = 1 - v
    return arc_values


def H_total(T, verts):
    if len(verts) <= 1:
        return 1.0
    total = 0.0
    for p in permutations(verts):
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def compute_G_direct(T, B_w, u, v):
    total = 0.0
    for pi in permutations(B_w):
        for k in range(len(pi) - 1):
            if pi[k] == u and pi[k + 1] == v:
                wt = 1.0
                for j in range(len(pi) - 1):
                    if j != k:
                        wt *= T(pi[j], pi[j + 1])
                total += wt
    return total


def compute_cv(T, W, w_target, v_target):
    """Compute c_v = dP/dT(w,v) for the identity P = interior + fE + chains.

    c_v = sum_{u!=v} p_u*G(u,v) - sum_{v'!=v} q_{v'}*G(v,v') + E(v)
          + d(chains)/dT(w,v)
    """
    I, J = 0, 1
    p = {v: T(v, I) for v in W}
    q = {v: T(J, v) for v in W}
    B_w = [x for x in W if x != w_target]
    B_v = [x for x in B_w if x != v_target]

    # Term 1: sum_{u!=v} p_u * G(u,v)
    term1 = sum(p[u] * compute_G_direct(T, B_w, u, v_target) for u in B_v)

    # Term 2: -sum_{v'!=v} q_{v'} * G(v,v')
    term2 = -sum(q[vp] * compute_G_direct(T, B_w, v_target, vp) for vp in B_v)

    # Term 3: E(v)
    hs_val = 0.0
    he_val = 0.0
    for pi in permutations(B_w):
        wt = 1.0
        for k in range(len(pi) - 1):
            wt *= T(pi[k], pi[k + 1])
        if pi[0] == v_target:
            hs_val += wt
        if pi[-1] == v_target:
            he_val += wt
    E_v = hs_val - he_val

    # Term 4: d(chains)/dT(w,v)
    # Chains starting (w, v, ...): d/dT(w,v) gives -2*p_last * path_wt * H(comp)
    # Chains ending (..., v, w): d/dT(w,v) of T(v,w)=1-T(w,v) gives -1
    #   so contribution = 2*q_first * path_wt * H(comp)
    chain_deriv = 0.0

    for chain_len in range(3, len(W) + 1, 2):
        for combo in combinations(W, chain_len):
            if w_target not in combo or v_target not in combo:
                continue
            comp = [x for x in W if x not in combo]
            H_comp = H_total(T, comp)

            for perm in permutations(combo):
                # Only chains with T(w,v) dependence
                if perm[0] == w_target and perm[1] == v_target:
                    # Chain (w, v, ..., u_k): weight = T(w,v) * T(v,...) * ...
                    # d/dT(w,v) = T(v,perm[2])*...*T(perm[-2],perm[-1])
                    path_wt = 1.0
                    for k in range(1, chain_len - 1):
                        path_wt *= T(perm[k], perm[k + 1])
                    chain_deriv += 2 * (-p[perm[-1]]) * path_wt * H_comp

                elif perm[-1] == w_target and perm[-2] == v_target:
                    # Chain (..., v, w): weight = ... * T(v,w) = ... * (1-T(w,v))
                    # d/dT(w,v) = -1 * rest
                    path_wt = 1.0
                    for k in range(chain_len - 2):
                        path_wt *= T(perm[k], perm[k + 1])
                    chain_deriv += 2 * (-q[perm[0]]) * (-1) * path_wt * H_comp

    return term1 + term2 + E_v + chain_deriv


def compute_c0(T, W, w_target):
    """Compute c_0 = P evaluated at T(w,v)=0 for all v.

    At T(w,v)=0: T(v,w) = 1, f(v) = q_v, perm(u,v) = 1*q_v + p_u*0 = q_v.

    interior|_{w=0} = sum_{u!=v} q_v * G(u,v) = sum_v q_v * G_in(v)
    fE|_{w=0} = sum_v q_v * E(v)
    chains|_{w=0} = chains ending with w only (starting chains have T(w,*)=0)
    """
    I, J = 0, 1
    p = {v: T(v, I) for v in W}
    q = {v: T(J, v) for v in W}
    B_w = [x for x in W if x != w_target]

    # Interior: sum q_v * G_in(v)
    interior_0 = 0.0
    for v in B_w:
        G_in_v = sum(compute_G_direct(T, B_w, u, v) for u in B_w if u != v)
        interior_0 += q[v] * G_in_v

    # fE: sum q_v * E(v)
    fE_0 = 0.0
    for v in B_w:
        hs = 0.0
        he = 0.0
        for pi in permutations(B_w):
            wt = 1.0
            for k in range(len(pi) - 1):
                wt *= T(pi[k], pi[k + 1])
            if pi[0] == v:
                hs += wt
            if pi[-1] == v:
                he += wt
        fE_0 += q[v] * (hs - he)

    # Chains ending with w (T(v,w)=1, T(w,v)=0): chains starting with w vanish.
    chains_0 = 0.0
    for chain_len in range(3, len(W) + 1, 2):
        for combo in combinations(W, chain_len):
            if w_target not in combo:
                continue
            comp = [x for x in W if x not in combo]
            H_comp = H_total(T, comp)

            for perm in permutations(combo):
                if perm[-1] == w_target:
                    # Chain (..., w): weight uses T(perm[-2], w) = 1
                    # and NO T(w,*) arcs
                    path_wt = 1.0
                    for k in range(chain_len - 1):
                        if perm[k] == w_target:
                            continue  # T(w, next) = 0
                        if perm[k + 1] == w_target:
                            path_wt *= 1  # T(prev, w) = 1
                        else:
                            path_wt *= T(perm[k], perm[k + 1])
                    # Actually, let me be more careful.
                    # Chain perm = (u_1, ..., u_{k-1}, w)
                    # Weight = T(u_1,u_2)*...*T(u_{k-2},u_{k-1})*T(u_{k-1},w)
                    # T(u_{k-1},w) = 1-T(w,u_{k-1}) = 1 (since T(w,*)=0)
                    path_wt = 1.0
                    for k in range(chain_len - 2):  # edges between non-w vertices
                        path_wt *= T(perm[k], perm[k + 1])
                    # last edge T(perm[-2], w) = 1
                    chains_0 += 2 * (-q[perm[0]]) * path_wt * H_comp

                elif perm[0] == w_target:
                    # Chain (w, ...): weight includes T(w, perm[1]) = 0
                    pass  # vanishes

    return interior_0 + fE_0 + chains_0


def verify_coefficients(n, num_trials=5):
    """Verify c_0 = 0 and c_v = 0 for all v."""
    print(f"\n=== Coefficient verification at n={n} ===\n")

    W = list(range(2, n))
    max_err = 0.0

    for trial in range(min(num_trials, 3)):
        arc_values = random_s0_tournament(n, seed=42 + trial)
        T = make_T(n, arc_values)

        for w in W[:1]:
            B_w = [x for x in W if x != w]

            # Compute c_0
            c0 = compute_c0(T, W, w)
            max_err = max(max_err, abs(c0))

            # Compute c_v for each v
            cv_vals = {}
            for v in B_w:
                cv = compute_cv(T, W, w, v)
                cv_vals[v] = cv
                max_err = max(max_err, abs(cv))

            if trial < 2:
                print(f"  trial {trial}, w={w}:")
                print(f"    c_0 = {c0:.2e}")
                for v in B_w:
                    print(f"    c_{v} = {cv_vals[v]:.2e}")

    print(f"\n  Max |c| = {max_err:.2e}")
    return max_err < 1e-10


def verify_linearity_symbolic():
    """Symbolically verify that P is linear in each T(w,v)."""
    if not HAS_SYMPY:
        print("  (sympy not available)")
        return

    print("\n=== Symbolic linearity check at n=5 ===\n")

    # n=5: B_w = {u,v}, W = {w,u,v}
    a, b = symbols('a b')  # T(w,u) = a, T(w,v) = b
    t = symbols('t')  # T(u,v) = t
    p_u, p_v = symbols('p_u p_v')
    q_u, q_v = 1 - p_u, 1 - p_v

    # perm(u,v) = T(u,w)*q_v + p_u*T(w,v) = (1-a)*q_v + p_u*b
    perm_uv = (1 - a) * q_v + p_u * b
    perm_vu = (1 - b) * q_u + p_v * a

    f_u = a + q_u
    f_v = b + q_v

    # G values at m=2: G(u,v) = G(v,u) = 1
    G_uv = Sym(1)
    G_vu = Sym(1)

    # Interior
    interior = perm_uv * G_uv + perm_vu * G_vu

    # fE: E(u) = T(u,v) - T(v,u) = 2t-1, E(v) = 1-2t
    E_u = 2 * t - 1
    E_v = 1 - 2 * t
    fE = f_u * E_u + f_v * E_v

    # Chains
    # (w,u,v): 2*(-p_v)*a*t
    # (w,v,u): 2*(-p_u)*b*(1-t)
    # (u,v,w): 2*(-q_u)*t*(1-b)
    # (v,u,w): 2*(-q_v)*(1-t)*(1-a)
    chains = (2 * (-p_v) * a * t +
              2 * (-p_u) * b * (1 - t) +
              2 * (-q_u) * t * (1 - b) +
              2 * (-q_v) * (1 - t) * (1 - a))

    P = expand(interior + fE + chains)
    print(f"  P = {P}")

    # Check linearity in a and b
    from sympy import degree, Poly
    P_poly = Poly(P, a, b)
    print(f"  Degree in (a,b): {P_poly.total_degree()}")
    print(f"  Max degree in a: {degree(P, a)}")
    print(f"  Max degree in b: {degree(P, b)}")

    # Check if P contains a*b term
    P_ab = P.coeff(a * b)
    print(f"  Coefficient of a*b: {P_ab}")

    # Extract coefficients
    c_0 = P.subs(a, 0).subs(b, 0)
    c_a = P.coeff(a).subs(b, 0)  # but this may include b terms
    # Better: since P = c0 + c_a*a + c_b*b (linear):
    P_at_a0_b0 = expand(P.subs(a, 0).subs(b, 0))
    dP_da = expand(P.diff(a))
    dP_db = expand(P.diff(b))

    print(f"\n  c_0 = P(a=0,b=0) = {P_at_a0_b0}")
    print(f"  c_a = dP/da = {dP_da}")
    print(f"  c_b = dP/db = {dP_db}")

    # Verify c_0, c_a, c_b = 0
    if P_at_a0_b0 == 0 and dP_da == 0 and dP_db == 0:
        print("\n  *** All coefficients zero: P = 0 confirmed! ***")
    else:
        # Try simplifying with p+q=1
        from sympy import simplify
        c0_s = simplify(P_at_a0_b0.subs(p_u + q_u, 1).subs(q_u, 1 - p_u).subs(q_v, 1 - p_v))
        ca_s = simplify(dP_da.subs(q_u, 1 - p_u).subs(q_v, 1 - p_v))
        cb_s = simplify(dP_db.subs(q_u, 1 - p_u).subs(q_v, 1 - p_v))
        print(f"\n  After p+q=1 substitution:")
        print(f"  c_0 = {c0_s}")
        print(f"  c_a = {ca_s}")
        print(f"  c_b = {cb_s}")


def verify_cv_structure_n6():
    """At n=6, compute c_v symbolically for P = interior_raw + fE + chains.

    BUG FIX: interior_raw is ONLY the interior insertion terms (not boundary).
    """
    if not HAS_SYMPY:
        print("  (sympy not available)")
        return

    print("\n=== Symbolic c_v at n=6 (FIXED) ===\n")

    a, b, c_var = symbols('a b c')
    t_uv, t_ux, t_vx = symbols('t_uv t_ux t_vx')
    p_u, p_v, p_x = symbols('p_u p_v p_x')
    q_u, q_v, q_x = 1 - p_u, 1 - p_v, 1 - p_x

    arcs = {
        ('w', 'u'): a, ('u', 'w'): 1 - a,
        ('w', 'v'): b, ('v', 'w'): 1 - b,
        ('w', 'x'): c_var, ('x', 'w'): 1 - c_var,
        ('u', 'v'): t_uv, ('v', 'u'): 1 - t_uv,
        ('u', 'x'): t_ux, ('x', 'u'): 1 - t_ux,
        ('v', 'x'): t_vx, ('x', 'v'): 1 - t_vx,
    }
    p_vals = {'u': p_u, 'v': p_v, 'x': p_x}
    q_vals = {'u': q_u, 'v': q_v, 'x': q_x}

    def T(s, t):
        if s == t: return Sym(0)
        return arcs.get((s, t), 1 - arcs.get((t, s), Sym(0)))

    B_w = ['u', 'v', 'x']

    # Compute interior_raw (ONLY interior insertion terms)
    interior_raw = Sym(0)
    for pi in permutations(B_w):
        B_wt = Sym(1)
        for k in range(len(pi) - 1):
            B_wt *= T(pi[k], pi[k + 1])
        m1 = len(pi)
        for pos in range(1, m1):  # ONLY interior positions
            u_k = pi[pos - 1]
            u_k1 = pi[pos]
            numer = T(u_k, 'w') * q_vals[u_k1] + p_vals[u_k] * T('w', u_k1)
            val = numer * B_wt / T(u_k, u_k1)
            interior_raw += val
    interior_raw = expand(interior_raw)

    # fE = sum_v f(v) * E(v) where f(v) = T(w,v) + q_v, E(v) = h_start - h_end
    fE = Sym(0)
    for v in B_w:
        f_v = T('w', v) + q_vals[v]
        hs = Sym(0)
        he = Sym(0)
        for pi in permutations(B_w):
            wt = Sym(1)
            for k in range(len(pi) - 1):
                wt *= T(pi[k], pi[k + 1])
            if pi[0] == v:
                hs += wt
            if pi[-1] == v:
                he += wt
        fE += f_v * (hs - he)
    fE = expand(fE)

    # Chains (length 3 only at n=6)
    W_verts = ['w'] + B_w
    chains = Sym(0)
    for combo in combinations(W_verts, 3):
        if 'w' not in combo:
            continue
        comp = [x for x in W_verts if x not in combo]
        H_comp = Sym(0)
        if len(comp) <= 1:
            H_comp = Sym(1)
        else:
            for pi in permutations(comp):
                wt = Sym(1)
                for k in range(len(pi) - 1):
                    wt *= T(pi[k], pi[k + 1])
                H_comp += wt
        for perm in permutations(combo):
            internal = Sym(1)
            for k in range(2):
                internal *= T(perm[k], perm[k + 1])
            if perm[0] == 'w':
                chains += 2 * (-p_vals[perm[-1]]) * internal * H_comp
            elif perm[-1] == 'w':
                chains += 2 * (-q_vals[perm[0]]) * internal * H_comp
    chains = expand(chains)

    P = expand(interior_raw + fE + chains)
    print(f"  interior_raw terms: {len(interior_raw.as_ordered_terms())}")
    print(f"  fE terms: {len(fE.as_ordered_terms())}")
    print(f"  chains terms: {len(chains.as_ordered_terms())}")
    print(f"  P terms: {len(P.as_ordered_terms()) if P != 0 else 0}")
    print(f"  P = 0? {P == 0}")

    if P != 0:
        from sympy import simplify
        P_s = simplify(P.subs(q_u, 1 - p_u).subs(q_v, 1 - p_v).subs(q_x, 1 - p_x))
        print(f"  After p+q=1: P = {P_s}")

    if P == 0:
        print("\n  *** P = interior_raw + fE + chains = 0 PROVED at n=6 ***")
        # Compute derivatives
        from sympy import diff
        for var, name in [(a, 'T(w,u)'), (b, 'T(w,v)'), (c_var, 'T(w,x)')]:
            dP = expand(diff(P, var))
            print(f"  dP/d{name} = {dP}")
    else:
        # Check linearity
        from sympy import diff, degree
        print(f"\n  Degree in a: {degree(P, a)}")
        print(f"  Degree in b: {degree(P, b)}")
        print(f"  Degree in c: {degree(P, c_var)}")

        for var, name in [(a, 'T(w,u)'), (b, 'T(w,v)'), (c_var, 'T(w,x)')]:
            dP = expand(diff(P, var))
            dP_s = simplify(dP.subs(q_u, 1-p_u).subs(q_v, 1-p_v).subs(q_x, 1-p_x))
            print(f"  dP/d{name} = {dP_s}")

        P0 = expand(P.subs(a, 0).subs(b, 0).subs(c_var, 0))
        P0_s = simplify(P0.subs(q_u, 1-p_u).subs(q_v, 1-p_v).subs(q_x, 1-p_x))
        print(f"  P(0,0,0) = {P0_s}")


if __name__ == "__main__":
    # Numerical verification of coefficients
    for n in [5, 6, 7, 8]:
        verify_coefficients(n, 5 if n <= 7 else 3)

    # Symbolic verification
    verify_linearity_symbolic()

    # n=6 symbolic
    verify_cv_structure_n6()
