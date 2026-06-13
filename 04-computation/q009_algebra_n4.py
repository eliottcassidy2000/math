#!/usr/bin/env python3
"""
Exact algebraic analysis of B(Li,Rj) at n=4.

At n=4, W = {a, b} with a=2, b=3. Internal arc t = T[a][b].
Interface: p_w = T[w][i], q_w = T[j][w].

We compute B(Li,Rj) and B(Lj,Ri) symbolically and verify equality.
Then express B in (s,d) coordinates and verify s-evenness.

Instance: opus-2026-03-05-S4
"""

import random
from itertools import permutations


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


def compute_B_Lj_Ri(n, arc_values):
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
        Lj = h_end(T, S + [J], J)
        Ri = h_start(T, [I] + R, I)
        total += sign * Lj * Ri
    return total


def verify_n4_algebra():
    """Verify the correct algebraic formula for B(Li,Rj) at n=4."""
    print("=== n=4 Algebraic Verification ===\n")
    I, J = 0, 1

    random.seed(42)
    max_err_LiRj = 0.0
    max_err_LjRi = 0.0
    max_identity = 0.0

    for trial in range(200):
        pa = random.uniform(-2, 3)
        pb = random.uniform(-2, 3)
        qa = random.uniform(-2, 3)
        qb = random.uniform(-2, 3)
        t = random.uniform(-2, 3)

        arc_values = {
            (I, J): 1.0, (J, I): 0.0,
            (2, I): pa, (I, 2): 1 - pa,
            (3, I): pb, (I, 3): 1 - pb,
            (J, 2): qa, (2, J): 1 - qa,
            (J, 3): qb, (3, J): 1 - qb,
            (2, 3): t, (3, 2): 1 - t,
        }

        B_num = compute_B_Li_Rj(4, arc_values)
        # Correct formula:
        # B(Li,Rj) = pa + qb + t*(qa - qb + pb - pa) - pa*qb - pb*qa
        B_alg = pa + qb + t * (qa - qb + pb - pa) - pa * qb - pb * qa
        # = pa(1-t) + qb(1-t) + t*qa + t*pb - pa*qb - pb*qa

        err = abs(B_num - B_alg)
        max_err_LiRj = max(max_err_LiRj, err)

        B2_num = compute_B_Lj_Ri(4, arc_values)
        # B(Lj,Ri): swap p_w <-> 1-q_w in the formula
        # Lj uses T[u][j] = 1-q_u instead of T[u][i] = p_u
        # Ri uses T[i][w] = 1-p_w instead of T[j][w] = q_w
        B2_alg = (1 - qa) + (1 - pb) + t * ((1 - pa) - (1 - pb) + (1 - qb) - (1 - qa)) \
                 - (1 - qa) * (1 - pb) - (1 - qb) * (1 - pa)
        # = (1-qa)(1-t) + (1-pb)(1-t) + t*(1-pa) + t*(1-qb) - (1-qa)*(1-pb) - (1-qb)*(1-pa)

        err2 = abs(B2_num - B2_alg)
        max_err_LjRi = max(max_err_LjRi, err2)

        identity_err = abs(B_num - B2_num)
        max_identity = max(max_identity, identity_err)

        if trial < 5:
            print(f"  trial {trial}: B(Li,Rj)={B_num:.6f}, alg={B_alg:.6f}, err={err:.2e}")
            print(f"           B(Lj,Ri)={B2_num:.6f}, alg={B2_alg:.6f}, err={err2:.2e}")
            print(f"           Identity: {identity_err:.2e}")

    print(f"\n  Max err B(Li,Rj): {max_err_LiRj:.2e}")
    print(f"  Max err B(Lj,Ri): {max_err_LjRi:.2e}")
    print(f"  Max |B(Li,Rj)-B(Lj,Ri)|: {max_identity:.2e}")


def verify_s_coords_n4():
    """Express B in (s,d) coords and verify s-evenness at n=4."""
    print("\n=== n=4 in (s,d) coordinates ===\n")

    random.seed(42)
    max_diff = 0.0

    for trial in range(200):
        da = random.uniform(-2, 2)
        db = random.uniform(-2, 2)
        sa = random.uniform(-2, 2)
        sb = random.uniform(-2, 2)
        t = random.uniform(-2, 3)

        pa = (1 - sa + da) / 2
        qa = (1 - sa - da) / 2
        pb = (1 - sb + db) / 2
        qb = (1 - sb - db) / 2

        B = pa * (1 - t) + qb * (1 - t) + t * qa + t * pb - pa * qb - pb * qa

        # In (s,d) coords, the formula is:
        # B = 1/2 + (da-db)(1-2t)/2 + da*db/2 - sa*sb/2
        B_sd = 0.5 + (da - db) * (1 - 2 * t) / 2 + da * db / 2 - sa * sb / 2

        err = abs(B - B_sd)
        if trial < 5:
            print(f"  trial {trial}: B={B:.6f}, B_sd={B_sd:.6f}, err={err:.2e}")

        # Verify s-evenness: B(+s) = B(-s)
        B_neg = 0.5 + (da - db) * (1 - 2 * t) / 2 + da * db / 2 - (-sa) * (-sb) / 2
        # = same thing since (-sa)*(-sb) = sa*sb
        diff = abs(B - B_neg)
        max_diff = max(max_diff, diff)

    print(f"\n  Max |B(+s)-B(-s)|: {max_diff:.2e}")
    print(f"\n  n=4 formula: B = 1/2 + (da-db)(1-2t)/2 + da*db/2 - sa*sb/2")
    print(f"  s-dependence is ONLY through sa*sb (even). NO linear s terms.")
    print(f"  Confirmed: B is even in s at n=4.")


def verify_identity_larger_n(n, num_trials=500):
    """Verify B(Li,Rj) = B(Lj,Ri) at given n with extended range reals."""
    print(f"\n=== B(Li,Rj) = B(Lj,Ri) at n={n} ===\n")
    I, J = 0, 1
    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)
    max_diff = 0.0

    for trial in range(num_trials):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.uniform(-2, 3)
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        B1 = compute_B_Li_Rj(n, arc_values)
        B2 = compute_B_Lj_Ri(n, arc_values)
        diff = abs(B1 - B2)
        max_diff = max(max_diff, diff)

        if trial < 5:
            print(f"  trial {trial}: B(Li,Rj)={B1:.6f}, B(Lj,Ri)={B2:.6f}, diff={diff:.2e}")

    print(f"\n  Max |B(Li,Rj) - B(Lj,Ri)| over {num_trials} trials: {max_diff:.2e}")
    if max_diff < 1e-8:
        print(f"  CONFIRMED: B(Li,Rj) = B(Lj,Ri) as polynomial identity at n={n}")
    else:
        print(f"  FAILED!")


def analyze_s_structure_n5():
    """
    At n=5, W = {a, b, c}. Analyze the s-dependence of B(Li,Rj).

    Each monomial in B involves at most one p_u (left boundary) and one q_w
    (right boundary). In (s,d) coords, p_u = (1-s_u+d_u)/2 and q_w = (1-s_w-d_w)/2.

    The product p_u * q_w = (1-s_u+d_u)(1-s_w-d_w)/4 has terms:
    - degree 0 in s: (1+d_u)(1-d_w)/4
    - degree 1 in s_u: -(1-d_w)/4 * s_u
    - degree 1 in s_w: -(1+d_u)/4 * s_w
    - degree 2 (s_u*s_w): s_u*s_w/4

    For u=w: p_u*q_u = [(1-s_u)^2 - d_u^2]/4. Degree 0: (1-d_u^2)/4. Degree 2: s_u^2/4.
    NO degree-1 terms when u=w! (Because p_u*q_u is quadratic in s_u.)

    Question: at n=5, what is the s-polynomial structure? What are the degree-2 terms?
    """
    print("\n=== s-polynomial structure at n=5 ===\n")

    I, J = 0, 1
    W = [2, 3, 4]

    random.seed(42)

    # Test: compute B at s=0 and at various s perturbations to extract coefficients
    for trial in range(3):
        da = random.uniform(-1, 1)
        db = random.uniform(-1, 1)
        dc = random.uniform(-1, 1)
        # Internal arcs
        tab = random.uniform(-1, 2)
        tac = random.uniform(-1, 2)
        tbc = random.uniform(-1, 2)

        d = {2: da, 3: db, 4: dc}

        def make_arc_values(s):
            av = {(I, J): 1.0, (J, I): 0.0}
            for w in W:
                pw = (1 - s[w] + d[w]) / 2
                qw = (1 - s[w] - d[w]) / 2
                av[(w, I)] = pw
                av[(I, w)] = 1 - pw
                av[(J, w)] = qw
                av[(w, J)] = 1 - qw
            av[(2, 3)] = tab
            av[(3, 2)] = 1 - tab
            av[(2, 4)] = tac
            av[(4, 2)] = 1 - tac
            av[(3, 4)] = tbc
            av[(4, 3)] = 1 - tbc
            return av

        s0 = {2: 0, 3: 0, 4: 0}
        B0 = compute_B_Li_Rj(5, make_arc_values(s0))

        # Extract linear coefficients by finite differences
        eps = 1e-5
        linear = {}
        for w in W:
            sp = dict(s0)
            sp[w] = eps
            Bp = compute_B_Li_Rj(5, make_arc_values(sp))
            linear[w] = (Bp - B0) / eps

        # Extract quadratic coefficients
        quadratic = {}
        for i_idx, w1 in enumerate(W):
            for w2 in W[i_idx:]:
                if w1 == w2:
                    sp = dict(s0)
                    sp[w1] = eps
                    Bp = compute_B_Li_Rj(5, make_arc_values(sp))
                    sm = dict(s0)
                    sm[w1] = -eps
                    Bm = compute_B_Li_Rj(5, make_arc_values(sm))
                    quadratic[(w1, w2)] = (Bp + Bm - 2 * B0) / eps ** 2
                else:
                    spp = dict(s0)
                    spp[w1] = eps
                    spp[w2] = eps
                    spm = dict(s0)
                    spm[w1] = eps
                    spm[w2] = -eps
                    smp = dict(s0)
                    smp[w1] = -eps
                    smp[w2] = eps
                    smm = dict(s0)
                    smm[w1] = -eps
                    smm[w2] = -eps
                    Bpp = compute_B_Li_Rj(5, make_arc_values(spp))
                    Bpm = compute_B_Li_Rj(5, make_arc_values(spm))
                    Bmp = compute_B_Li_Rj(5, make_arc_values(smp))
                    Bmm = compute_B_Li_Rj(5, make_arc_values(smm))
                    quadratic[(w1, w2)] = (Bpp - Bpm - Bmp + Bmm) / (4 * eps ** 2)

        print(f"  trial {trial}: B(s=0) = {B0:.6f}")
        print(f"    Linear coefficients (should be ~0):")
        for w in W:
            print(f"      dB/ds_{w} = {linear[w]:.6e}")
        print(f"    Quadratic coefficients:")
        for (w1, w2), val in quadratic.items():
            label = f"s_{w1}*s_{w2}" if w1 != w2 else f"s_{w1}^2"
            print(f"      d2B/d{label} = {val:.6f}")
        print()


def prove_identity_n3():
    """
    Clean proof at n=3.

    W = {w}. B(Li,Rj) = Li(empty)*Rj({w}) - Li({w})*Rj(empty)
    = 1*q_w - p_w*1 = q_w - p_w

    B(Lj,Ri) = Lj(empty)*Ri({w}) - Lj({w})*Ri(empty)
    = 1*(1-p_w) - (1-q_w)*1 = q_w - p_w

    Identical! The proof uses T[a][b] + T[b][a] = 1 (complementarity).
    """
    print("\n=== Proof at n=3 ===\n")
    print("  B(Li,Rj) = q_w - p_w = -d_w")
    print("  B(Lj,Ri) = (1-p_w) - (1-q_w) = q_w - p_w = -d_w")
    print("  Identity holds by T[a][b] + T[b][a] = 1.")
    print("  B depends only on d_w, NO s-dependence.")


def prove_identity_n4():
    """
    Clean proof at n=4.

    B(Li,Rj) = pa(1-t) + qb(1-t) + t*qa + t*pb - pa*qb - pb*qa

    Under sigma (p_w -> 1-q_w, q_w -> 1-p_w):
    B(Lj,Ri) = (1-qa)(1-t) + (1-pb)(1-t) + t*(1-pa) + t*(1-qb) - (1-qa)(1-pb) - (1-qb)(1-pa)

    Expand:
    = (1-t) - qa(1-t) + (1-t) - pb(1-t) + t - t*pa + t - t*qb
      - [1-pb-qa+qa*pb] - [1-pa-qb+qb*pa]
    = 2(1-t) + 2t - qa(1-t) - pb(1-t) - t*pa - t*qb - 2 + pa + pb + qa + qb - qa*pb - qb*pa
    = 2 - 2 + pa + qa + pb + qb - qa(1-t) - pb(1-t) - t*pa - t*qb - qa*pb - qb*pa
    = pa(1-t) + qa*t + pb*t + qb(1-t) - qa*pb - qb*pa

    = B(Li,Rj). QED.
    """
    print("\n=== Proof at n=4 ===\n")
    print("  B(Li,Rj) = pa(1-t) + qb(1-t) + t*qa + t*pb - pa*qb - pb*qa")
    print("  B(Lj,Ri) obtained by p_w -> 1-q_w, q_w -> 1-p_w")
    print("  = (1-qa)(1-t) + (1-pb)(1-t) + t(1-pa) + t(1-qb) - (1-qa)(1-pb) - (1-qb)(1-pa)")
    print("  After expansion: = pa(1-t) + qa*t + pb*t + qb(1-t) - qa*pb - qb*pa")
    print("  = B(Li,Rj). QED.")
    print()
    print("  In (s,d) coords: B = 1/2 + (da-db)(1-2t)/2 + da*db/2 - sa*sb/2")
    print("  s-dependence: only sa*sb (degree 2, even). Zero linear s-terms.")


if __name__ == "__main__":
    prove_identity_n3()
    prove_identity_n4()
    verify_n4_algebra()
    verify_s_coords_n4()
    verify_identity_larger_n(5, 500)
    verify_identity_larger_n(6, 200)
    verify_identity_larger_n(7, 100)
    analyze_s_structure_n5()
