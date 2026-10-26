# O(q^4) evaluation of a seven-vertex Paley sum S_H by summing out a degree-3 vertex v in closed form.
# With N(v) = {a, b, c} and order orientation, the factors at v are s_v chi((t - x_a)(t - x_b)(t - x_c)),
# s_v = (-1)^(number of arcs v -> u, i.e. neighbours u > v).  Summed over t in F_q this is
# s_v chi(x_b - x_a) L(lambda), L(l) = sum_s chi(s(s-1)(s-l)), lambda = (x_c - x_a)/(x_b - x_a); the values t = x_z of
# the three non-neighbours z must be removed (injectivity).  Fixing x_a = 0, x_b = 1 (the summand stays AGL(1)-invariant)
#   S_H = s_v [ A - sum_z S'_z ],
# A = sum over the other four free variables of the remaining pair factors times L(x_c), and S'_z the six-vertex sum of
# H - v times chi((x_z)(x_z - 1)(x_z - x_c)).  Every term is a contraction of a complete K4 network: O(q^4).
# mac-mini-2026-10-25-pyth.
import itertools
import numpy as np
from sum7 import chi_table

def L_table(p, chi):
    t = np.arange(p)
    base = chi[(t * (t - 1)) % p]
    return np.array([float((base * chi[(t - l) % p]).sum()) for l in range(p)])

def contract4(p, free, u, Mf):
    """sum over injective-by-factor x_free (4 variables) of prod u_v(x_v) prod_{pairs} Mf[(v,w)][x_v, x_w]"""
    f1, f2, f3, f4 = free
    tot = 0.0
    M12, M13, M14 = Mf[(f1, f2)], Mf[(f1, f3)], Mf[(f1, f4)]
    M23, M24, M34 = Mf[(f2, f3)], Mf[(f2, f4)], Mf[(f3, f4)]
    for x1 in range(p):
        w1 = u[f1][x1]
        if w1 == 0: continue
        u2 = u[f2] * M12[x1]; u3 = u[f3] * M13[x1]; u4 = u[f4] * M14[x1]
        K = (M23 * u3[None, :]) @ M34                 # K[x2, x4] = sum_x3 M23[x2,x3] u3[x3] M34[x3,x4]
        tot += w1 * float(u2 @ (K * M24) @ u4)
    return tot

def S7_fast(p, edges, v=None, chi=None, L=None):
    chi = chi_table(p) if chi is None else chi
    L = L_table(p, chi) if L is None else L
    E = {tuple(sorted(e)) for e in edges}
    nb = lambda x: sorted(y for y in range(7) if tuple(sorted((x, y))) in E)
    if v is None:
        v = next(x for x in range(7) if len(nb(x)) == 3)
    a, b, c = nb(v)
    sv = (-1) ** sum(1 for u_ in (a, b, c) if u_ > v)
    others = [x for x in range(7) if x not in (v, a, b, c)]     # the three non-neighbours z
    idx = np.arange(p)
    val = {a: 0, b: 1}
    def build(extra_unary=None, toggles=None):
        """network on the six vertices != v with x_a = 0, x_b = 1; toggles: set of pairs whose arc status flips
        (with the sign bookkeeping done by the caller); extra_unary: dict vertex -> weight vector"""
        Ez = set(e for e in E if v not in e)
        if toggles:
            for e in toggles: Ez ^= {e}
        const = 1.0
        if (min(a, b), max(a, b)) in Ez:
            const = chi[(val[max(a, b)] - val[min(a, b)]) % p]
        free = [x for x in range(7) if x not in (v, a, b)]
        u = {}
        for x in free:
            w = np.ones(p)
            for f in (a, b):
                cst = val[f]; arc = (min(x, f), max(x, f)) in Ez
                if x < f: w = w * (chi[(cst - idx) % p] if arc else (idx != cst))
                else:     w = w * (chi[(idx - cst) % p] if arc else (idx != cst))
            if extra_unary and x in extra_unary: w = w * extra_unary[x]
            u[x] = w
        Mf = {}
        D = (idx[None, :] - idx[:, None]) % p
        for x, y in itertools.combinations(free, 2):
            F = chi[D] if (x, y) in Ez else 1.0 - np.eye(p)
            Mf[(x, y)] = F; Mf[(y, x)] = F.T
        return const * contract4(p, free, u, Mf)
    # main term: weight L(x_c) on vertex c (lambda = x_c since x_a = 0, x_b = 1); chi(x_b - x_a) = chi(1) = 1
    A = build(extra_unary={c: L})
    corr = 0.0
    for z in others:
        # multiply by chi(x_z - x_a) chi(x_z - x_b) chi(x_z - x_c): toggle {z,u} for u in a, b, c, with sign
        # chi(x_z - x_u) = (order factor of pair {z,u}) if u < z, and -(order factor) if z < u
        sign = (-1) ** sum(1 for u_ in (a, b, c) if z < u_)
        corr += sign * build(toggles={(min(z, u_), max(z, u_)) for u_ in (a, b, c)})
    return int(round(sv * (A - corr)))

if __name__ == '__main__':
    import json, time
    from sum7 import S7
    G = json.load(open('genuine37.json'))
    for p in (11, 19, 23, 31):
        chi = chi_table(p); L = L_table(p, chi)
        bad = [i for i, H in enumerate(G) if any(len([y for y in range(7) if tuple(sorted((x, y))) in {tuple(sorted(e)) for e in H['edges']}]) == 3 for x in range(7))
               and S7_fast(p, H['edges'], chi=chi, L=L) != S7(p, H['edges'])]
        print(f"p = {p}: S7_fast = S7 on every survivor with a degree-3 vertex: {not bad} {bad}")
    chi = chi_table(401); L = L_table(401, chi)
    t0 = time.time(); S7_fast(401, G[2]['edges'], chi=chi, L=L); print(f"timing at q = 401: {time.time() - t0:.1f}s per sum")
