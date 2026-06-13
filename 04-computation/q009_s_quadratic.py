#!/usr/bin/env python3
"""
Analyze the quadratic s-structure of B(Li,Rj).

Hypothesis: B = B_0(d,t) + sum_{w<u} c_{wu}(d,t) * s_w * s_u
(no diagonal s_w^2 terms, no linear s_w terms)

If true, what are the c_{wu} coefficients? Do they have a clean formula?

At n=4: B = 1/2 + (da-db)(1-2t)/2 + da*db/2 - sa*sb/2
So c_{ab} = -1/2 (constant! independent of d and t!)

Question: is c_{wu} = -1/2 for all w,u at all n?

Instance: opus-2026-03-05-S4
"""

from itertools import permutations, combinations
import random


def h_end(T, verts, v):
    if len(verts) == 1:
        return 1.0 if v == verts[0] else 0.0
    total = 0.0
    for p in permutations(verts):
        if p[-1] != v:
            continue
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def h_start(T, verts, v):
    if len(verts) == 1:
        return 1.0 if v == verts[0] else 0.0
    total = 0.0
    for p in permutations(verts):
        if p[0] != v:
            continue
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def compute_B(n, arc_values):
    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)
    def T(a, b):
        if a == b: return 0
        return arc_values.get((a, b), 1 - arc_values.get((b, a), 0))
    total = 0.0
    for smask in range(1 << m):
        S = [W[bit] for bit in range(m) if smask & (1 << bit)]
        R = [W[bit] for bit in range(m) if not (smask & (1 << bit))]
        sign = (-1) ** len(S)
        Li = h_end(T, S + [I], I)
        Rj = h_start(T, [J] + R, J)
        total += sign * Li * Rj
    return total


def extract_s_coefficients(n, d_vals, internal_vals, eps=1e-5):
    """Extract all s-polynomial coefficients of B at given (d, internal arcs)."""
    I, J = 0, 1
    W = list(range(2, n))
    internal_arcs = [(a, b) for a in W for b in W if a < b]

    def make_av(s_vals):
        av = {(I, J): 1.0, (J, I): 0.0}
        for w in W:
            pw = (1 - s_vals.get(w, 0) + d_vals[w]) / 2
            qw = (1 - s_vals.get(w, 0) - d_vals[w]) / 2
            av[(w, I)] = pw
            av[(I, w)] = 1 - pw
            av[(J, w)] = qw
            av[(w, J)] = 1 - qw
        for (a, b) in internal_arcs:
            av[(a, b)] = internal_vals.get((a, b), 0.5)
            av[(b, a)] = 1 - av[(a, b)]
        return av

    s0 = {w: 0 for w in W}
    B0 = compute_B(n, make_av(s0))

    # Linear coefficients
    linear = {}
    for w in W:
        sp = dict(s0); sp[w] = eps
        sm = dict(s0); sm[w] = -eps
        linear[w] = (compute_B(n, make_av(sp)) - compute_B(n, make_av(sm))) / (2 * eps)

    # Quadratic coefficients (cross and diagonal)
    quadratic = {}
    for w in W:
        sp = dict(s0); sp[w] = eps
        sm = dict(s0); sm[w] = -eps
        quadratic[(w, w)] = (compute_B(n, make_av(sp)) + compute_B(n, make_av(sm)) - 2 * B0) / eps**2

    for w1, w2 in combinations(W, 2):
        spp = dict(s0); spp[w1] = eps; spp[w2] = eps
        spm = dict(s0); spm[w1] = eps; spm[w2] = -eps
        smp = dict(s0); smp[w1] = -eps; smp[w2] = eps
        smm = dict(s0); smm[w1] = -eps; smm[w2] = -eps
        quadratic[(w1, w2)] = (
            compute_B(n, make_av(spp)) - compute_B(n, make_av(spm))
            - compute_B(n, make_av(smp)) + compute_B(n, make_av(smm))
        ) / (4 * eps**2)

    return B0, linear, quadratic


def test_constant_coefficient(n, num_trials=100):
    """Test whether c_{wu} = -1/2 for all pairs at n."""
    print(f"=== Quadratic s-coefficients at n={n} ===\n")

    W = list(range(2, n))
    internal_arcs = [(a, b) for a in W for b in W if a < b]

    random.seed(42)
    max_linear = 0.0
    max_diag = 0.0
    cross_vals = {pair: [] for pair in combinations(W, 2)}

    for trial in range(num_trials):
        d_vals = {w: random.uniform(-1, 1) for w in W}
        internal_vals = {(a, b): random.uniform(-1, 2) for (a, b) in internal_arcs}

        B0, linear, quadratic = extract_s_coefficients(n, d_vals, internal_vals)

        for w in W:
            max_linear = max(max_linear, abs(linear[w]))
            max_diag = max(max_diag, abs(quadratic[(w, w)]))

        for pair in combinations(W, 2):
            cross_vals[pair].append(quadratic[pair])

        if trial < 3:
            print(f"  trial {trial}: B0={B0:.4f}")
            print(f"    linear: {[f'{linear[w]:.2e}' for w in W]}")
            print(f"    diag:   {[f'{quadratic[(w,w)]:.2e}' for w in W]}")
            print(f"    cross:  {[(pair, f'{quadratic[pair]:.4f}') for pair in combinations(W,2)]}")

    print(f"\n  Max |linear|: {max_linear:.2e}")
    print(f"  Max |diagonal|: {max_diag:.2e}")

    # Check if cross coefficients are constant (-1/2)
    print(f"\n  Cross-term analysis (is c_{{wu}} = -1/2?):")
    for pair in combinations(W, 2):
        vals = cross_vals[pair]
        mean = sum(vals) / len(vals)
        var = sum((v - mean)**2 for v in vals) / len(vals)
        mn, mx = min(vals), max(vals)
        print(f"    c_{pair}: mean={mean:.4f}, std={var**0.5:.4f}, range=[{mn:.4f}, {mx:.4f}]")


def test_cross_formula(n, num_trials=100):
    """
    If c_{wu} != -1/2, find what it depends on.

    At n=4: c_{ab} = -1/2 (constant).
    At n=5: check if c_{wu} depends on internal arcs between w and u.
    """
    print(f"\n=== Cross-term formula exploration at n={n} ===\n")

    W = list(range(2, n))
    internal_arcs = [(a, b) for a in W for b in W if a < b]

    random.seed(42)

    # Test: fix all internal arcs except T[w1][w2], vary it
    for w1, w2 in combinations(W, 2):
        print(f"  Pair ({w1},{w2}):")
        for trial in range(3):
            d_vals = {w: random.uniform(-0.5, 0.5) for w in W}
            # Fix other internal arcs
            internal_vals = {(a, b): random.uniform(0, 1) for (a, b) in internal_arcs}

            c_at_t = []
            for t_val in [0.0, 0.25, 0.5, 0.75, 1.0]:
                internal_vals[(w1, w2)] = t_val
                _, _, quad = extract_s_coefficients(n, d_vals, internal_vals)
                c_at_t.append(quad[(w1, w2)])

            print(f"    trial {trial}: c at T[{w1}][{w2}]=0,.25,.5,.75,1: "
                  f"{[f'{v:.4f}' for v in c_at_t]}")

        # Also check: does c depend on d_w1, d_w2?
        internal_vals = {(a, b): random.uniform(0, 1) for (a, b) in internal_arcs}
        c_at_d = []
        for dw1 in [-0.5, 0, 0.5]:
            d_vals = {w: 0.0 for w in W}
            d_vals[w1] = dw1
            _, _, quad = extract_s_coefficients(n, d_vals, internal_vals)
            c_at_d.append(quad[(w1, w2)])
        print(f"    c at d_{w1}=-0.5,0,0.5 (d others=0): {[f'{v:.4f}' for v in c_at_d]}")
        print()


if __name__ == "__main__":
    test_constant_coefficient(4, 50)
    test_constant_coefficient(5, 50)
    test_constant_coefficient(6, 30)
    test_cross_formula(5, 10)
