#!/usr/bin/env python3
"""
Creative hypothesis testing around Claim (B) and OCF.

Claim (B): sum_{S⊆W\{v}} (-1)^|S| H(S) h_start(W\S, v) = (-1)^{m+1} h_end(W, v)

We test:
H1: Does h_end variant hold? (replace h_start with h_end on LHS)
H2: Does it hold for general digraphs (drop T[a][b]+T[b][a]=1)?
H3: Determinantal formulation — does the alternating sum relate to det of path matrix?
H4: Does a DUAL identity hold? (swap H and h_start/h_end roles)
H5: Does the identity hold "vertex by vertex"? (split H(S) = sum_u h_start(S,u))
H6: Is there a generating function F(x) = sum x^|S| H(S) h_start(R,v) with a nice form?
H7: Pfaffian connection — does H(T) relate to Pfaffian of skew-symmetric tournament matrix?

Instance: kind-pasteur-2026-03-05-S10
"""

from itertools import permutations, combinations
import random
import sys


def make_tournament(m, seed=None):
    """Random tournament on m vertices with real-valued arcs (T[a][b]+T[b][a]=1)."""
    rng = random.Random(seed)
    arcs = {}
    for a in range(m):
        for b in range(a + 1, m):
            val = rng.uniform(-1, 2)
            arcs[(a, b)] = val
            arcs[(b, a)] = 1 - val
    return arcs


def make_digraph(m, seed=None):
    """Random digraph (T[a][b]+T[b][a] NOT necessarily 1)."""
    rng = random.Random(seed)
    arcs = {}
    for a in range(m):
        for b in range(m):
            if a != b:
                arcs[(a, b)] = rng.uniform(0, 1)
    return arcs


def T_func(arcs):
    def T(a, b):
        if a == b:
            return 0
        return arcs.get((a, b), 0)
    return T


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


def h_start(T, verts, v):
    if len(verts) <= 1:
        return 1.0 if (len(verts) == 0 or verts[0] == v) else 0.0
    total = 0.0
    for p in permutations(verts):
        if p[0] != v:
            continue
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def h_end(T, verts, v):
    if len(verts) <= 1:
        return 1.0 if (len(verts) == 0 or verts[0] == v) else 0.0
    total = 0.0
    for p in permutations(verts):
        if p[-1] != v:
            continue
        w = 1.0
        for k in range(len(p) - 1):
            w *= T(p[k], p[k + 1])
        total += w
    return total


def claim_b_lhs(T, W, v):
    """Compute LHS of Claim (B): sum_{S⊆W\{v}} (-1)^|S| H(S) h_start(W\S, v)"""
    Wv = [w for w in W if w != v]
    total = 0.0
    for smask in range(1 << len(Wv)):
        S = [Wv[bit] for bit in range(len(Wv)) if smask & (1 << bit)]
        R = [w for w in W if w not in S]
        total += ((-1) ** len(S)) * H_total(T, S) * h_start(T, R, v)
    return total


def claim_b_rhs(T, W, v):
    """Compute RHS of Claim (B): (-1)^{m+1} h_end(W, v)"""
    m = len(W)
    return ((-1) ** (m + 1)) * h_end(T, W, v)


# ============================================================
# HYPOTHESIS 1: h_end variant
# sum_{S⊆W\{v}} (-1)^|S| H(S) h_end(W\S, v) = ???
# ============================================================
def test_h1(m_max=6, trials=50):
    print("=" * 60)
    print("H1: h_end variant — replace h_start with h_end on LHS")
    print("=" * 60)
    for m in range(2, m_max + 1):
        W = list(range(m))
        max_err_original = 0.0
        # Check if h_end variant equals something simple
        results = []
        for trial in range(trials):
            arcs = make_tournament(m, seed=trial * 100 + m)
            T = T_func(arcs)
            v = 0
            Wv = [w for w in W if w != v]

            # h_end variant LHS
            lhs_end = 0.0
            for smask in range(1 << len(Wv)):
                S = [Wv[bit] for bit in range(len(Wv)) if smask & (1 << bit)]
                R = [w for w in W if w not in S]
                lhs_end += ((-1) ** len(S)) * H_total(T, S) * h_end(T, R, v)

            # Compare to various RHS candidates
            rhs1 = ((-1) ** (m + 1)) * h_start(T, W, v)
            results.append((lhs_end, rhs1))
            max_err_original = max(max_err_original, abs(lhs_end - rhs1))

        if max_err_original < 1e-8:
            print(f"  m={m}: h_end variant = (-1)^(m+1) h_start(W,v)  [CONFIRMED, max_err={max_err_original:.2e}]")
        else:
            # Try other candidates
            print(f"  m={m}: h_end variant != (-1)^(m+1) h_start(W,v)  [max_err={max_err_original:.2e}]")
            # Check ratio
            ratios = [r[0] / r[1] if abs(r[1]) > 1e-10 else None for r in results[:5]]
            print(f"    Ratios LHS/candidate: {[f'{r:.4f}' if r else 'N/A' for r in ratios]}")


# ============================================================
# HYPOTHESIS 2: General digraphs (T[a][b]+T[b][a] != 1)
# ============================================================
def test_h2(m_max=5, trials=50):
    print("\n" + "=" * 60)
    print("H2: Does Claim (B) hold for general digraphs?")
    print("=" * 60)
    for m in range(2, m_max + 1):
        W = list(range(m))
        fails = 0
        max_err = 0.0
        for trial in range(trials):
            arcs = make_digraph(m, seed=trial * 100 + m)
            T = T_func(arcs)
            v = 0
            lhs = claim_b_lhs(T, W, v)
            rhs = claim_b_rhs(T, W, v)
            err = abs(lhs - rhs)
            max_err = max(max_err, err)
            if err > 1e-6:
                fails += 1
        print(f"  m={m}: {fails}/{trials} failures, max_err={max_err:.2e}")


# ============================================================
# HYPOTHESIS 3: Determinant connection
# For skew-symmetric A[i][j] = T[i][j] - 1/2, does the alternating sum
# relate to det(A) or Pfaffian(A)?
# ============================================================
def test_h3(m_max=6, trials=30):
    print("\n" + "=" * 60)
    print("H3: Determinant / Pfaffian connection")
    print("=" * 60)
    try:
        import numpy as np
    except ImportError:
        print("  numpy not available, skipping")
        return

    for m in range(2, m_max + 1):
        W = list(range(m))
        for trial in range(min(trials, 5)):
            arcs = make_tournament(m, seed=trial * 100 + m)
            T = T_func(arcs)

            # Build skew-symmetric matrix A[i][j] = T(i,j) - 1/2
            A = np.zeros((m, m))
            for i in range(m):
                for j in range(m):
                    if i != j:
                        A[i][j] = arcs[(i, j)] - 0.5

            det_A = np.linalg.det(A)
            H_W = H_total(T, W)

            # Test various relationships
            h_end_vals = [h_end(T, W, v) for v in W]
            h_start_vals = [h_start(T, W, v) for v in W]

            if trial == 0:
                print(f"\n  m={m}, trial {trial}:")
                print(f"    H(W) = {H_W:.4f}")
                print(f"    det(A) = {det_A:.6f}")
                print(f"    h_end = {[f'{x:.4f}' for x in h_end_vals]}")
                print(f"    h_start = {[f'{x:.4f}' for x in h_start_vals]}")
                print(f"    sum h_end = {sum(h_end_vals):.4f} (should = H(W) = {H_W:.4f})")
                print(f"    sum h_start = {sum(h_start_vals):.4f} (should = H(W) = {H_W:.4f})")

                # Check if Claim B LHS for each v gives something related to det
                for v in W[:3]:
                    lhs = claim_b_lhs(T, W, v)
                    rhs = claim_b_rhs(T, W, v)
                    print(f"    v={v}: Claim B LHS={lhs:.4f}, RHS={rhs:.4f}, err={abs(lhs-rhs):.2e}")

                # Cofactor expansion?
                if m <= 6:
                    # Test: does h_end(W,v) relate to cofactors of A?
                    for v in range(min(m, 3)):
                        # Minor: delete row v, col v
                        minor = np.delete(np.delete(A, v, 0), v, 1)
                        det_minor = np.linalg.det(minor)
                        ratio = h_end_vals[v] / det_minor if abs(det_minor) > 1e-10 else None
                        print(f"    h_end(v={v})/det(A_{{v,v}}) = {ratio if ratio else 'N/A'}")


# ============================================================
# HYPOTHESIS 5: Vertex-split — H(S) = sum_u h_start(S,u) = sum_u h_end(S,u)
# Does Claim B "localize" when we decompose H(S)?
# ============================================================
def test_h5(m_max=5, trials=30):
    print("\n" + "=" * 60)
    print("H5: Per-vertex decomposition of Claim B")
    print("H(S) = sum_u h_end(S,u). Substitute into Claim B:")
    print("sum_S (-1)^|S| [sum_u h_end(S,u)] h_start(R,v) = (-1)^{m+1} h_end(W,v)")
    print("Does this hold per-u? I.e., for each u:")
    print("  sum_{S containing u} (-1)^|S| h_end(S,u) h_start(R,v) = ???")
    print("=" * 60)
    for m in range(3, m_max + 1):
        W = list(range(m))
        v = 0
        Wv = [w for w in W if w != v]

        arcs = make_tournament(m, seed=42 + m)
        T = T_func(arcs)

        total_rhs = claim_b_rhs(T, W, v)
        print(f"\n  m={m}, v={v}: RHS = {total_rhs:.6f}")

        per_u = {}
        for u in W:
            val = 0.0
            for smask in range(1 << len(Wv)):
                S = [Wv[bit] for bit in range(len(Wv)) if smask & (1 << bit)]
                R = [w for w in W if w not in S]
                if u in S:
                    val += ((-1) ** len(S)) * h_end(T, S, u) * h_start(T, R, v)
                elif len(S) == 0 and u == v:
                    # S=empty, u=v: H(empty) = 1, h_start(W, v) — u=v is the v-term
                    pass
            per_u[u] = val

        # S=empty contributes H(empty)*h_start(W,v) = h_start(W,v)
        # This has no "u" in S, so it's a separate boundary term
        boundary = h_start(T, W, v)
        total_from_parts = boundary + sum(per_u.values())
        print(f"    Boundary (S=empty) = {boundary:.6f}")
        for u in W:
            print(f"    u={u}: per-u contribution = {per_u[u]:.6f}")
        print(f"    Total = {total_from_parts:.6f}, check = {abs(total_from_parts - total_rhs):.2e}")


# ============================================================
# HYPOTHESIS 6: Generating function F(x) = sum x^|S| H(S) h_start(R,v)
# What does F look like? Is D(x) = F(x) - RHS divisible by (1+x)?
# ============================================================
def test_h6(m_max=6, trials=5):
    print("\n" + "=" * 60)
    print("H6: Generating function F(x) = sum x^k [sum_{|S|=k} H(S) h_start(R,v)]")
    print("=" * 60)
    for m in range(2, m_max + 1):
        W = list(range(m))
        v = 0
        Wv = [w for w in W if w != v]

        arcs = make_tournament(m, seed=42 + m)
        T = T_func(arcs)

        coeffs = [0.0] * m  # coeffs[k] = sum_{|S|=k} H(S) h_start(R,v)
        for smask in range(1 << len(Wv)):
            S = [Wv[bit] for bit in range(len(Wv)) if smask & (1 << bit)]
            R = [w for w in W if w not in S]
            k = len(S)
            coeffs[k] += H_total(T, S) * h_start(T, R, v)

        h_end_v = h_end(T, W, v)
        print(f"\n  m={m}, v={v}:")
        for k in range(m):
            print(f"    F[{k}] = {coeffs[k]:.6f}")
        print(f"    h_end(W,v) = {h_end_v:.6f}")
        print(f"    F(-1) = {sum((-1)**k * coeffs[k] for k in range(m)):.6f}")
        print(f"    Expected F(-1) = {(-1)**(m+1) * h_end_v:.6f}")

        # Check: is F(x) - (-1)^{m+1}*h_end(W,v)/(1+x)^? divisible by (1+x)?
        # Compute F(x) at a few x values
        for x in [0, 1, 2, -0.5]:
            Fx = sum(x**k * coeffs[k] for k in range(m))
            print(f"    F({x}) = {Fx:.6f}")


# ============================================================
# HYPOTHESIS 7: Pfaffian connection
# For even m, does H(T) relate to Pfaffian of some matrix?
# ============================================================
def test_h7(trials=20):
    print("\n" + "=" * 60)
    print("H7: Pfaffian connection for {0,1}-valued tournaments")
    print("=" * 60)
    try:
        import numpy as np
    except ImportError:
        print("  numpy not available, skipping")
        return

    def pfaffian_4x4(A):
        """Pfaffian of 4x4 skew-symmetric matrix."""
        return A[0,1]*A[2,3] - A[0,2]*A[1,3] + A[0,3]*A[1,2]

    for m in [3, 4, 5, 6]:
        W = list(range(m))
        print(f"\n  m={m}:")
        for trial in range(min(trials, 5)):
            # {0,1} tournament
            arcs = {}
            rng = random.Random(trial * 100 + m)
            for a in range(m):
                for b in range(a + 1, m):
                    if rng.random() < 0.5:
                        arcs[(a, b)] = 1
                        arcs[(b, a)] = 0
                    else:
                        arcs[(a, b)] = 0
                        arcs[(b, a)] = 1
            T = T_func(arcs)

            H_W = H_total(T, W)

            # Skew matrix A[i][j] = T(i,j) - T(j,i) = 2*T(i,j) - 1
            A = np.zeros((m, m))
            for i in range(m):
                for j in range(m):
                    if i != j:
                        A[i][j] = 2 * arcs.get((i, j), 0) - 1

            det_A = np.linalg.det(A)

            # Score sequence
            scores = [sum(arcs.get((i, j), 0) for j in range(m) if j != i) for i in range(m)]

            if trial < 3:
                print(f"    trial {trial}: H={H_W:.0f}, det(A)={det_A:.1f}, scores={scores}")


# ============================================================
# HYPOTHESIS 8 (NEW): Claim B via TRANSFER MATRIX on subsets
# Define matrix M(v) where M(v)[S,R] = (-1)^|S| H(S) h_start(R,v) for S∪R=W\{v}
# Claim B says the sum of all entries = (-1)^{m+1} h_end(W,v)
# What are the eigenvalues of related matrices?
# ============================================================
def test_h8(m_max=5, trials=5):
    print("\n" + "=" * 60)
    print("H8: Transfer matrix structure")
    print("Define f(S) = H(S), g(R) = h_start(R∪{v}, v)")
    print("Claim B: sum_S (-1)^|S| f(S) g(W\\S) = (-1)^{m+1} h_end(W,v)")
    print("This is a CONVOLUTION. Check: does (f * mu_hat)(R) simplify?")
    print("=" * 60)
    try:
        import numpy as np
    except ImportError:
        print("  numpy not available, skipping")
        return

    for m in range(3, m_max + 1):
        W = list(range(m))
        v = 0
        Wv = [w for w in W if w != v]
        n_sub = len(Wv)

        arcs = make_tournament(m, seed=42 + m)
        T = T_func(arcs)

        # Build vectors f and g indexed by subsets of Wv
        f = np.zeros(1 << n_sub)  # f[S] = H(S)
        g = np.zeros(1 << n_sub)  # g[S] = h_start(S ∪ {v}, v)

        for smask in range(1 << n_sub):
            S = [Wv[bit] for bit in range(n_sub) if smask & (1 << bit)]
            R_verts = [w for w in W if w not in S]  # complement includes v

            f[smask] = H_total(T, S)
            g[smask] = h_start(T, R_verts, v)

        # Claim B: sum_S (-1)^|S| f[S] g[comp(S)] = target
        target = ((-1) ** (m + 1)) * h_end(T, W, v)

        # Compute the convolution
        conv = 0.0
        for smask in range(1 << n_sub):
            comp = ((1 << n_sub) - 1) ^ smask
            k = bin(smask).count('1')
            conv += ((-1) ** k) * f[smask] * g[comp]

        print(f"\n  m={m}, v={v}:")
        print(f"    Convolution = {conv:.6f}, target = {target:.6f}, err = {abs(conv - target):.2e}")

        # What does the "Hadamard transform" of f and g look like?
        # Walsh-Hadamard: F_hat[y] = sum_S (-1)^{<S,y>} f[S]
        def walsh_hadamard(vec, n):
            """In-place Walsh-Hadamard transform."""
            result = vec.copy()
            h = 1
            while h < (1 << n):
                for i in range(0, 1 << n, 2 * h):
                    for j in range(i, i + h):
                        x = result[j]
                        y = result[j + h]
                        result[j] = x + y
                        result[j + h] = x - y
                h *= 2
            return result

        f_hat = walsh_hadamard(f, n_sub)
        g_hat = walsh_hadamard(g, n_sub)

        print(f"    f_hat (WHT of H(S)): {[f'{x:.4f}' for x in f_hat]}")
        print(f"    g_hat (WHT of h_start): {[f'{x:.4f}' for x in g_hat]}")
        print(f"    f_hat * g_hat (pointwise): {[f'{f_hat[i]*g_hat[i]:.4f}' for i in range(len(f_hat))]}")

        # The signed convolution at x=-1 corresponds to evaluating the WHT at the all-1 vector
        # sum_S (-1)^|S| f[S] g[comp(S)] = (1/2^n) sum_y f_hat[y] g_hat[comp(y)] * (-1)^{popcount(y)}
        # Actually, let me be more careful.
        # The standard convolution f*g [mask] = sum_{S} f[S] g[mask ^ S]
        # Our sum is f * g at the all-1s mask, with the (-1)^|S| sign.
        # This is sum_S (-1)^|S| f[S] g[all\S] = sum_S (-1)^{popcount(S)} f[S] g[all^S]
        # Since all^S = comp(S), this is just the "signed correlation" at the all-1s point.

        # Under WHT: (f * g)[mask] = (1/2^n) * (f_hat .* g_hat)[mask]
        # But our sum has the (-1)^|S| sign, which is the same as multiplying f by chi_{all}(S) = (-1)^|S|
        # So our sum = (f * chi_all * g)[all] = ... this gets messy.
        # Skip the spectral analysis for now.


# ============================================================
# HYPOTHESIS 9 (CREATIVE): "Path reversal duality"
# Claim B has h_start on LHS and h_end on RHS.
# In tournament T^op (all arcs reversed):
# h_start(T^op, R, v) = h_end(T, R, v) and vice versa
# H(T^op) = H(T) (reversing a path in T^op gives a path in T)
# Does this give a SYMMETRY of Claim B?
# ============================================================
def test_h9(m_max=5, trials=20):
    print("\n" + "=" * 60)
    print("H9: Path reversal duality under T -> T^op")
    print("In T^op: h_start <-> h_end, H unchanged")
    print("Claim B in T^op: sum (-1)^|S| H(S) h_end(R,v) = (-1)^{m+1} h_start(W,v)")
    print("Compare this with H1 (h_end variant)...")
    print("=" * 60)
    for m in range(2, m_max + 1):
        W = list(range(m))
        max_err = 0.0
        for trial in range(trials):
            arcs = make_tournament(m, seed=trial * 100 + m)
            T = T_func(arcs)
            v = 0
            Wv = [w for w in W if w != v]

            # Apply Claim B to T^op:
            # sum (-1)^|S| H_Top(S) h_start_Top(R, v)
            # = sum (-1)^|S| H_T(S) h_end_T(R, v)  [since H(T^op)=H(T), h_start(T^op)=h_end(T)]
            # = (-1)^{m+1} h_end_Top(W, v) = (-1)^{m+1} h_start_T(W, v)

            lhs = 0.0
            for smask in range(1 << len(Wv)):
                S = [Wv[bit] for bit in range(len(Wv)) if smask & (1 << bit)]
                R = [w for w in W if w not in S]
                lhs += ((-1) ** len(S)) * H_total(T, S) * h_end(T, R, v)

            rhs = ((-1) ** (m + 1)) * h_start(T, W, v)
            err = abs(lhs - rhs)
            max_err = max(max_err, err)

        status = "CONFIRMED" if max_err < 1e-8 else "FAILED"
        print(f"  m={m}: Dual Claim B (h_end variant) = (-1)^(m+1) h_start(W,v)  [{status}, max_err={max_err:.2e}]")


# ============================================================
# HYPOTHESIS 10 (CREATIVE): "Recursive structure"
# Can we express the m-vertex Claim B in terms of (m-1)-vertex instances?
# Fix vertex a != v. Condition on whether a is first or last in the path on S.
# ============================================================
def test_h10(m_max=6, trials=20):
    print("\n" + "=" * 60)
    print("H10: Inductive decomposition — condition on a specific vertex")
    print("Fix vertex a != v. For S containing a, split H(S) by position of a:")
    print("  H(S) = h_start(S,a) + h_end(S,a) + h_mid(S,a)")
    print("Can we relate the a-start and a-end contributions to (m-1) case?")
    print("=" * 60)
    for m in range(3, m_max + 1):
        W = list(range(m))
        v = 0
        a = 1  # vertex to condition on
        Wv = [w for w in W if w != v]
        Wva = [w for w in W if w != v and w != a]

        arcs = make_tournament(m, seed=42 + m)
        T = T_func(arcs)

        # Split LHS by whether a is in S or not
        # Part 1: a not in S (a is in R)
        part_a_in_R = 0.0
        for smask in range(1 << len(Wva)):
            S = [Wva[bit] for bit in range(len(Wva)) if smask & (1 << bit)]
            R = [w for w in W if w not in S]  # contains v and a
            part_a_in_R += ((-1) ** len(S)) * H_total(T, S) * h_start(T, R, v)

        # Part 2: a in S (a not in R)
        part_a_in_S = 0.0
        for smask in range(1 << len(Wva)):
            S_without_a = [Wva[bit] for bit in range(len(Wva)) if smask & (1 << bit)]
            S = S_without_a + [a]
            R = [w for w in W if w not in S]  # contains v but not a
            part_a_in_S += ((-1) ** len(S)) * H_total(T, S) * h_start(T, R, v)

        total = part_a_in_R + part_a_in_S
        expected = claim_b_rhs(T, W, v)

        # Now try: when a not in S, h_start(R, v) involves paths starting at v in a set containing a.
        # Can decompose: h_start(R, v) = h_start on R with a somewhere.
        # When a is first after v: sum_... T(v,a) * h_remaining
        # When a is elsewhere: ...

        # Claim B on W\{a} (m-1 vertices):
        W_minus_a = [w for w in W if w != a]
        cb_minus_a_lhs = claim_b_lhs(T, W_minus_a, v)
        cb_minus_a_rhs = claim_b_rhs(T, W_minus_a, v)

        print(f"\n  m={m}, v={v}, a={a}:")
        print(f"    Part (a in R) = {part_a_in_R:.6f}")
        print(f"    Part (a in S) = {part_a_in_S:.6f}")
        print(f"    Total = {total:.6f}, expected = {expected:.6f}")
        print(f"    Claim B on W\\{{a}} (m-1): LHS={cb_minus_a_lhs:.6f}, RHS={cb_minus_a_rhs:.6f}")
        print(f"    Ratio part_R / CB(m-1)_RHS = {part_a_in_R / cb_minus_a_rhs:.6f}" if abs(cb_minus_a_rhs) > 1e-10 else "")

        # CREATIVE: Does part_a_in_R relate to Claim B at m-1 in a nice way?
        # If we could show part_a_in_R = T(v,a) * CB_rhs(m-1) + correction terms...
        # that would give an inductive proof.

        # Try: part_a_in_R = sum over S ⊆ Wva of (-1)^|S| H(S) h_start(S^c ∪ {v,a}, v)
        # h_start({v,a,...}, v) = T(v,a) * h_start({a,...}, a) + sum of paths starting v -> non-a
        # So part_a_in_R = T(v,a) * [sum_S (-1)^|S| H(S) h_start_from_a(Wva\S ∪ {a}, a)]
        #                 + [paths starting v -> non-a vertex]

        # The first term is T(v,a) * Claim_B(W\{v}, a, restricted)... hmm, not quite.


if __name__ == "__main__":
    # Run all hypotheses
    test_h1()
    test_h2()
    test_h9()  # Dual (T^op) — this should confirm H1
    test_h6()
    test_h5()
    test_h3()
    test_h8()
    test_h10()
