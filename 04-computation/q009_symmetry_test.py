#!/usr/bin/env python3
"""
Test the sigma-symmetry of B(Li, Rj).

The substitution sigma: p_w -> 1-q_w, q_w -> 1-p_w (keeping internal arcs fixed)
negates s_w = 1-p_w-q_w and preserves d_w = p_w - q_w.

The main identity B(Li,Rj) = B(Lj,Ri) is equivalent to:
  B(Li,Rj)(T) = B(Li,Rj)(sigma(T))

i.e., B(Li,Rj) is invariant under sigma.

Reparametrize: p_w = (1-s_w+d_w)/2, q_w = (1-s_w-d_w)/2.
sigma: s -> -s, d -> d.

Test: is B(Li,Rj) an even function of s (with d and internal arcs fixed)?
If it depends on s only through even powers, the identity follows.

Instance: opus-2026-03-05-S4
"""

from itertools import permutations
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


def compute_B_Li_Rj(n, arc_values):
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


def test_sigma_invariance(n, num_trials=200):
    """Verify B(Li,Rj)(T) = B(Li,Rj)(sigma(T))."""
    print(f"=== Sigma Invariance at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))

    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]
    # Separate interface and internal arcs
    interface_arcs = [(w, I) for w in W] + [(J, w) for w in W]
    internal_arcs = [(a, b) for a in W for b in W if a < b]

    random.seed(42)
    max_diff = 0.0

    for trial in range(num_trials):
        # Random arc values
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.uniform(-1, 2)
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        B_orig = compute_B_Li_Rj(n, arc_values)

        # Apply sigma: p_w -> 1-q_w, q_w -> 1-p_w
        sigma_values = dict(arc_values)
        for w in W:
            old_pw = arc_values[(w, I)]  # T[w][i] = p_w
            old_qw = arc_values[(J, w)]  # T[j][w] = q_w
            new_pw = 1 - old_qw
            new_qw = 1 - old_pw
            sigma_values[(w, I)] = new_pw
            sigma_values[(I, w)] = 1 - new_pw
            sigma_values[(J, w)] = new_qw
            sigma_values[(w, J)] = 1 - new_qw

        B_sigma = compute_B_Li_Rj(n, sigma_values)

        diff = abs(B_orig - B_sigma)
        max_diff = max(max_diff, diff)

        if trial < 5:
            print(f"  trial {trial}: B_orig={B_orig:.6f}, B_sigma={B_sigma:.6f}, diff={diff:.2e}")

    print(f"\n  Max |B(T) - B(sigma(T))| = {max_diff:.2e}")
    if max_diff < 1e-8:
        print(f"  CONFIRMED: B(Li,Rj) is sigma-invariant (= polynomial identity)")
    else:
        print(f"  NOT sigma-invariant")


def test_s_parity(n, num_trials=100):
    """
    Fix d_w = p_w - q_w and internal arcs. Vary s_w.
    Check if B(Li,Rj) is an even function of s = (s_1,...,s_m).
    """
    print(f"\n=== S-Parity Test at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)
    internal_arcs = [(a, b) for a in W for b in W if a < b]

    random.seed(42)

    for trial in range(min(num_trials, 20)):
        # Fix d_w and internal arcs, vary s_w
        d = {w: random.uniform(-1, 1) for w in W}
        internal_vals = {(a, b): random.uniform(-1, 2) for (a, b) in internal_arcs}

        # With s-values as parameter
        def make_arc_values(s_vals):
            av = {(I, J): 1.0, (J, I): 0.0}
            for w in W:
                pw = (1 - s_vals[w] + d[w]) / 2
                qw = (1 - s_vals[w] - d[w]) / 2
                av[(w, I)] = pw
                av[(I, w)] = 1 - pw
                av[(J, w)] = qw
                av[(w, J)] = 1 - qw
            for (a, b) in internal_arcs:
                av[(a, b)] = internal_vals[(a, b)]
                av[(b, a)] = 1 - internal_vals[(a, b)]
            return av

        # Compute B at s and -s
        s_plus = {w: random.uniform(-1, 1) for w in W}
        s_minus = {w: -s_plus[w] for w in W}

        B_plus = compute_B_Li_Rj(n, make_arc_values(s_plus))
        B_minus = compute_B_Li_Rj(n, make_arc_values(s_minus))

        if trial < 10:
            print(f"  trial {trial}: B(+s)={B_plus:.6f}, B(-s)={B_minus:.6f}, "
                  f"diff={abs(B_plus - B_minus):.2e}, "
                  f"even={(B_plus+B_minus)/2:.4f}, odd={(B_plus-B_minus)/2:.4f}")

    print()
    # Check: B(+s) = B(-s) for all s? (i.e., B is even in s?)
    max_odd = 0.0
    for trial in range(num_trials):
        s_plus = {w: random.uniform(-2, 2) for w in W}
        s_minus = {w: -s_plus[w] for w in W}
        d = {w: random.uniform(-1, 1) for w in W}
        internal_vals = {(a, b): random.uniform(-1, 2) for (a, b) in internal_arcs}

        B_plus = compute_B_Li_Rj(n, make_arc_values(s_plus))
        B_minus = compute_B_Li_Rj(n, make_arc_values(s_minus))
        max_odd = max(max_odd, abs(B_plus - B_minus))

    print(f"  Max |B(+s) - B(-s)| = {max_odd:.2e}")
    if max_odd < 1e-8:
        print(f"  B is EVEN in s! The identity B(Li,Rj) = B(Lj,Ri) follows.")
    else:
        print(f"  B is NOT even in s. The odd part is non-zero.")
        print(f"  But sigma-invariance still holds (even though parity alone doesn't give it)")


def test_linear_s_coefficient(n, num_trials=50):
    """
    Decompose B(Li,Rj) as polynomial in s_w:
    B = B_0(d,t) + sum_w s_w * B_1w(d,t) + sum_{w<w'} s_w*s_w' * B_2ww'(d,t) + ...

    The sigma-invariance (s -> -s) means all ODD-degree terms vanish:
    B_1w = 0, B_3www' = 0, etc.

    Test: is the linear coefficient B_1w = dB/ds_w = 0?
    """
    print(f"\n=== Linear S-Coefficient Test at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    m = len(W)
    internal_arcs = [(a, b) for a in W for b in W if a < b]

    random.seed(42)

    for trial in range(min(num_trials, 10)):
        d = {w: random.uniform(-1, 1) for w in W}
        internal_vals = {(a, b): random.uniform(-1, 2) for (a, b) in internal_arcs}

        def make_arc_values(s_vals):
            av = {(I, J): 1.0, (J, I): 0.0}
            for w in W:
                pw = (1 - s_vals[w] + d[w]) / 2
                qw = (1 - s_vals[w] - d[w]) / 2
                av[(w, I)] = pw
                av[(I, w)] = 1 - pw
                av[(J, w)] = qw
                av[(w, J)] = 1 - qw
            for (a, b) in internal_arcs:
                av[(a, b)] = internal_vals[(a, b)]
                av[(b, a)] = 1 - internal_vals[(a, b)]
            return av

        # Compute dB/ds_w by finite difference
        eps = 1e-6
        s_base = {w: 0.0 for w in W}  # evaluate at s=0

        B_base = compute_B_Li_Rj(n, make_arc_values(s_base))

        derivs = {}
        for w in W:
            s_pert = dict(s_base)
            s_pert[w] = eps
            B_pert = compute_B_Li_Rj(n, make_arc_values(s_pert))
            derivs[w] = (B_pert - B_base) / eps

        if trial < 5:
            print(f"  trial {trial}: B(s=0)={B_base:.6f}, "
                  f"dB/ds = {[f'{derivs[w]:.4f}' for w in W]}")

    print()
    # Check: is dB/ds_w = 0 at ALL points (not just s=0)?
    max_deriv = 0.0
    for trial in range(num_trials):
        d = {w: random.uniform(-1, 1) for w in W}
        internal_vals = {(a, b): random.uniform(-1, 2) for (a, b) in internal_arcs}
        s_base = {w: random.uniform(-1, 1) for w in W}

        def make_arc_values_local(s_vals):
            av = {(I, J): 1.0, (J, I): 0.0}
            for w in W:
                pw = (1 - s_vals[w] + d[w]) / 2
                qw = (1 - s_vals[w] - d[w]) / 2
                av[(w, I)] = pw
                av[(I, w)] = 1 - pw
                av[(J, w)] = qw
                av[(w, J)] = 1 - qw
            for (a, b) in internal_arcs:
                av[(a, b)] = internal_vals[(a, b)]
                av[(b, a)] = 1 - internal_vals[(a, b)]
            return av

        B_base = compute_B_Li_Rj(n, make_arc_values_local(s_base))
        for w in W:
            s_pert = dict(s_base)
            s_pert[w] = s_base[w] + eps
            B_pert = compute_B_Li_Rj(n, make_arc_values_local(s_pert))
            deriv = (B_pert - B_base) / eps
            max_deriv = max(max_deriv, abs(deriv))

    print(f"  Max |dB/ds_w| at random points = {max_deriv:.4f}")
    if max_deriv < 1e-4:
        print(f"  B does NOT depend on s at all! B = B(d, internal arcs)")
    else:
        print(f"  B depends on s (gradient non-zero). "
              f"But the ALTERNATING parity might still hold.")


if __name__ == "__main__":
    test_sigma_invariance(4)
    test_sigma_invariance(5)

    test_s_parity(4)
    test_s_parity(5)

    test_linear_s_coefficient(4)
    test_linear_s_coefficient(5)
