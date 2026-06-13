#!/usr/bin/env python3
"""
Prove C_w + D_w = 0 for each w in W, where:
  C_w = dB/dp_w  (sensitivity of B(Li,Rj) to p_w = T[w][i])
  D_w = dB/dq_w  (sensitivity of B(Li,Rj) to q_w = T[j][w])

This is the KEY REDUCTION: if C_w + D_w = 0 for all w, then B(Li,Rj) is
even in s_w = 1 - p_w - q_w, which proves the signed adjacency identity
(= OCF).

Approach: compute C_w and D_w SYMBOLICALLY as polynomials in the remaining
arc variables, then verify cancellation.

Instance: opus-2026-03-05-S4
"""

from itertools import permutations
import random
from fractions import Fraction


def h_end(T, verts, v):
    """Hamiltonian paths on verts ending at v, weighted by T."""
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
    """Hamiltonian paths on verts starting at v, weighted by T."""
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
    """B(Li, Rj) = sum_S (-1)^|S| Li(S) Rj(W\S)."""
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


def compute_Cw_Dw_numerically(n, w_target, arc_values, eps=1e-7):
    """Compute C_w = dB/dp_w and D_w = dB/dq_w by finite differences."""
    I, J = 0, 1

    # Base value
    B0 = compute_B_Li_Rj(n, arc_values)

    # Perturb p_w = T[w][i]
    av_p = dict(arc_values)
    av_p[(w_target, I)] = arc_values[(w_target, I)] + eps
    av_p[(I, w_target)] = 1 - av_p[(w_target, I)]
    Bp = compute_B_Li_Rj(n, av_p)
    Cw = (Bp - B0) / eps

    # Perturb q_w = T[j][w]
    av_q = dict(arc_values)
    av_q[(J, w_target)] = arc_values[(J, w_target)] + eps
    av_q[(w_target, J)] = 1 - av_q[(J, w_target)]
    Bq = compute_B_Li_Rj(n, av_q)
    Dw = (Bq - B0) / eps

    return Cw, Dw


def test_Cw_Dw_cancellation(n, num_trials=200):
    """Test C_w + D_w = 0 for all w in W."""
    print(f"=== C_w + D_w = 0 Test at n={n} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)
    max_sum = 0.0

    for trial in range(num_trials):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.uniform(-2, 3)  # extended range for polynomial identity
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        for w in W:
            Cw, Dw = compute_Cw_Dw_numerically(n, w, arc_values)
            s = abs(Cw + Dw)
            max_sum = max(max_sum, s)

            if trial < 3:
                print(f"  trial {trial}, w={w}: C_w={Cw:.6f}, D_w={Dw:.6f}, "
                      f"C_w+D_w={Cw+Dw:.2e}")

    print(f"\n  Max |C_w + D_w| = {max_sum:.2e}")
    if max_sum < 1e-4:
        print(f"  CONFIRMED: C_w + D_w = 0 for all w (polynomial identity)")
    else:
        print(f"  FAILED: C_w + D_w != 0")


def analyze_Cw_structure(n, w_target, num_trials=50):
    """
    Understand what C_w and D_w depend on.

    C_w = dB/dp_w: how does B change when we increase T[w][i]?
    D_w = dB/dq_w: how does B change when we increase T[j][w]?

    Since B is multilinear in each arc variable, C_w and D_w are polynomials
    in the REMAINING variables (not involving p_w or q_w directly, but they
    may involve d_w = p_w - q_w through other terms).

    Actually: B is LINEAR in p_w (since p_w = T[w][i] appears at most once
    per monomial — w can be at most one boundary). Similarly LINEAR in q_w.
    So C_w = dB/dp_w is INDEPENDENT of p_w, and D_w = dB/dq_w is INDEPENDENT
    of q_w.

    But C_w may depend on q_w, and D_w may depend on p_w!

    Let's check: is C_w independent of q_w too? Is D_w independent of p_w?
    """
    print(f"\n=== C_w Structure Analysis at n={n}, w={w_target} ===\n")

    I, J = 0, 1
    W = list(range(2, n))
    free_arcs = [(a, b) for a in range(n) for b in range(a + 1, n) if (a, b) != (I, J)]

    random.seed(42)

    # Fix all arcs except q_w, vary q_w, check if C_w changes
    for trial in range(min(num_trials, 5)):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.uniform(-1, 2)
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        # Compute C_w at different q_w values
        cw_vals = []
        for qw in [0.0, 0.25, 0.5, 0.75, 1.0]:
            av = dict(arc_values)
            av[(J, w_target)] = qw
            av[(w_target, J)] = 1 - qw
            Cw, Dw = compute_Cw_Dw_numerically(n, w_target, av)
            cw_vals.append(Cw)

        cw_range = max(cw_vals) - min(cw_vals)
        print(f"  trial {trial}: C_w at q_w=0,.25,.5,.75,1: "
              f"{[f'{v:.4f}' for v in cw_vals]}, range={cw_range:.2e}")

    print()
    # Similarly: vary p_w, check if D_w changes
    for trial in range(min(num_trials, 5)):
        arc_values = {(I, J): 1.0, (J, I): 0.0}
        for (a, b) in free_arcs:
            v = random.uniform(-1, 2)
            arc_values[(a, b)] = v
            arc_values[(b, a)] = 1 - v

        dw_vals = []
        for pw in [0.0, 0.25, 0.5, 0.75, 1.0]:
            av = dict(arc_values)
            av[(w_target, I)] = pw
            av[(I, w_target)] = 1 - pw
            Cw, Dw = compute_Cw_Dw_numerically(n, w_target, av)
            dw_vals.append(Dw)

        dw_range = max(dw_vals) - min(dw_vals)
        print(f"  trial {trial}: D_w at p_w=0,.25,.5,.75,1: "
              f"{[f'{v:.4f}' for v in dw_vals]}, range={dw_range:.2e}")


def explicit_Cw_formula_n4(num_trials=100):
    """
    At n=4, W = {2, 3}. Compute C_w and D_w EXACTLY using algebra.

    B(Li, Rj) = sum_{S subset W} (-1)^|S| Li(S) Rj(W\S)

    For n=4, W = {a, b} where a=2, b=3.
    Subsets: empty, {a}, {b}, {a,b}

    Li(empty) = h_end({i}, i) = 1
    Li({a}) = h_end({a,i}, i) = T[a][i] = p_a
    Li({b}) = h_end({b,i}, i) = T[b][i] = p_b
    Li({a,b}) = h_end({a,b,i}, i) = T[a][b]*T[b][i] + T[b][a]*T[a][i]
              = T[a][b]*p_b + (1-T[a][b])*p_a = T[a][b]*(p_b-p_a) + p_a

    Rj(W) = h_start({j,a,b}, j) = T[j][a]*T[a][b] + T[j][b]*T[b][a]
           = (1-q_a)*T[a][b] + (1-q_b)*(1-T[a][b])
           = 1 - q_b + T[a][b]*(q_b - q_a)
    Rj({b}) = h_start({j,b}, j) = T[j][b] = 1-q_b
    Rj({a}) = h_start({j,a}, j) = T[j][a] = 1-q_a
    Rj(empty) = h_start({j}, j) = 1

    B = (+1)*Li(empty)*Rj({a,b}) + (-1)*Li({a})*Rj({b}) + (-1)*Li({b})*Rj({a})
        + (+1)*Li({a,b})*Rj(empty)

    B = [1 - q_b + t*(q_b-q_a)] - p_a*(1-q_b) - p_b*(1-q_a) + [t*(p_b-p_a)+p_a]
      where t = T[a][b] = T[2][3].

    Let me expand:
    = 1 - q_b + t*q_b - t*q_a - p_a + p_a*q_b - p_b + p_b*q_a + t*p_b - t*p_a + p_a
    = 1 - q_b + t*q_b - t*q_a - p_b + p_b*q_a + p_a*q_b + t*p_b - t*p_a
    = 1 - q_b - p_b + t*(q_b - q_a + p_b - p_a) + p_a*q_b + p_b*q_a

    dB/dp_a = t*(-1) + q_b = q_b - t
    dB/dq_a = t*(-1) + p_b = p_b - t
    dB/dp_b = -1 + t + q_a = q_a + t - 1
    dB/dq_b = -1 + t + p_a = p_a + t - 1

    C_a + D_a = (q_b - t) + (p_b - t) = p_b + q_b - 2t
    C_b + D_b = (q_a + t - 1) + (p_a + t - 1) = p_a + q_a + 2t - 2

    Hmm, these are NOT zero in general! Let me recheck...
    """
    print(f"\n=== Explicit n=4 Formula Check ===\n")

    I, J = 0, 1
    a, b = 2, 3

    random.seed(42)
    max_err = 0.0

    for trial in range(num_trials):
        pa = random.uniform(-1, 2)
        pb = random.uniform(-1, 2)
        qa = random.uniform(-1, 2)
        qb = random.uniform(-1, 2)
        t = random.uniform(-1, 2)  # T[a][b] = T[2][3]

        arc_values = {
            (I, J): 1.0, (J, I): 0.0,
            (a, I): pa, (I, a): 1 - pa,
            (b, I): pb, (I, b): 1 - pb,
            (J, a): qa, (a, J): 1 - qa,
            (J, b): qb, (b, J): 1 - qb,
            (a, b): t, (b, a): 1 - t,
        }

        # Numerical
        Ca_num, Da_num = compute_Cw_Dw_numerically(4, a, arc_values)
        Cb_num, Db_num = compute_Cw_Dw_numerically(4, b, arc_values)

        # Algebraic formulas
        Ca_alg = qb - t
        Da_alg = pb - t
        Cb_alg = qa + t - 1
        Db_alg = pa + t - 1

        err_Ca = abs(Ca_num - Ca_alg)
        err_Da = abs(Da_num - Da_alg)
        err_Cb = abs(Cb_num - Cb_alg)
        err_Db = abs(Db_num - Db_alg)

        sum_a = Ca_alg + Da_alg  # = pb + qb - 2t
        sum_b = Cb_alg + Db_alg  # = pa + qa + 2t - 2

        if trial < 5:
            print(f"  trial {trial}: C_a+D_a = {sum_a:.4f} (= pb+qb-2t), "
                  f"C_b+D_b = {sum_b:.4f} (= pa+qa+2t-2)")
            print(f"    formula errors: {err_Ca:.2e}, {err_Da:.2e}, {err_Cb:.2e}, {err_Db:.2e}")

        max_err = max(max_err, err_Ca, err_Da, err_Cb, err_Db)

    print(f"\n  Max formula error: {max_err:.2e}")
    print(f"\n  C_a + D_a = p_b + q_b - 2t")
    print(f"  C_b + D_b = p_a + q_a + 2t - 2")
    print(f"\n  These are NOT zero! So C_w + D_w = 0 is NOT the right condition.")
    print(f"  The s-evenness must work differently than per-vertex C_w + D_w = 0.")


def recheck_s_evenness_n4():
    """
    Let me recheck what the s-evenness condition actually says at n=4.

    B(Li,Rj) as a function of (s_a, s_b, d_a, d_b, t):
    p_w = (1-s_w+d_w)/2, q_w = (1-s_w-d_w)/2

    B = 1 - q_b - p_b + t*(q_b - q_a + p_b - p_a) + p_a*q_b + p_b*q_a

    Substituting:
    p_a = (1-sa+da)/2, q_a = (1-sa-da)/2
    p_b = (1-sb+db)/2, q_b = (1-sb-db)/2

    p_a + q_a = 1 - sa
    p_b + q_b = 1 - sb
    p_a - q_a = da
    p_b - q_b = db
    q_b - q_a = (sa-sb-da+db)/2 = ((sa-sb)-(da-db))/2
    p_b - p_a = ((sa-sb)+(db-da))/2 = ((sa-sb)+(db-da))/2

    q_b-q_a + p_b-p_a = sa-sb (the s-difference!)

    1 - q_b - p_b = 1 - (1-sb) = sb

    p_a*q_b = (1-sa+da)/2 * (1-sb-db)/2
    p_b*q_a = (1-sb+db)/2 * (1-sa-da)/2

    p_a*q_b + p_b*q_a = [(1-sa+da)(1-sb-db) + (1-sb+db)(1-sa-da)] / 4
    = [(1-sa)(1-sb) - da*db + da*(1-sb) - (1-sa)*db
       + (1-sa)(1-sb) - da*db - da*(1-sb) + (1-sa)*db] / 4
    Wait, let me be more careful.

    (1-sa+da)(1-sb-db) = (1-sa)(1-sb) - (1-sa)db + da(1-sb) - da*db
    (1-sb+db)(1-sa-da) = (1-sb)(1-sa) - (1-sb)da + db(1-sa) - db*da
    Sum = 2(1-sa)(1-sb) - 2*da*db
        + [-(1-sa)db + da(1-sb)] + [-(1-sb)da + db(1-sa)]
    = 2(1-sa)(1-sb) - 2*da*db + 0
    = 2(1-sa)(1-sb) - 2*da*db

    So p_a*q_b + p_b*q_a = [(1-sa)(1-sb) - da*db] / 2

    B = sb + t*(sa-sb) + [(1-sa)(1-sb) - da*db] / 2

    Let me expand:
    B = sb + t*sa - t*sb + (1 - sa - sb + sa*sb - da*db)/2
    = sb(1-t) + t*sa + 1/2 - sa/2 - sb/2 + sa*sb/2 - da*db/2

    Group by s-degree:
    - s-degree 0: 1/2 - da*db/2
    - s-degree 1 in sa: t*sa - sa/2 = sa*(t - 1/2)
    - s-degree 1 in sb: sb*(1-t) - sb/2 = sb*(1/2 - t)
    - s-degree 2: sa*sb/2

    So dB/dsa = t - 1/2 + sb/2
       dB/dsb = 1/2 - t + sa/2

    These are NOT zero! But wait: B(+s) = B(-s) was confirmed numerically.
    Let me check: B(sa,sb) vs B(-sa,-sb):

    B(sa,sb) = 1/2 - da*db/2 + sa*(t-1/2) + sb*(1/2-t) + sa*sb/2
    B(-sa,-sb) = 1/2 - da*db/2 - sa*(t-1/2) - sb*(1/2-t) + (-sa)(-sb)/2
              = 1/2 - da*db/2 - sa*(t-1/2) - sb*(1/2-t) + sa*sb/2

    B(sa,sb) - B(-sa,-sb) = 2*sa*(t-1/2) + 2*sb*(1/2-t) = 2*(t-1/2)*(sa-sb)

    THIS IS NOT ZERO when sa != sb and t != 1/2!

    So either my algebra is wrong, or the sigma-invariance test had a bug.
    Let me recheck numerically...
    """
    print(f"\n=== Recheck s-evenness at n=4 ===\n")

    I, J = 0, 1
    W = [2, 3]

    random.seed(123)

    for trial in range(10):
        da = random.uniform(-1, 1)
        db = random.uniform(-1, 1)
        t = random.uniform(0, 1)
        sa = random.uniform(-1, 1)
        sb = random.uniform(-1, 1)

        def make_arc_values(sa_val, sb_val):
            pa = (1 - sa_val + da) / 2
            qa = (1 - sa_val - da) / 2
            pb = (1 - sb_val + db) / 2
            qb = (1 - sb_val - db) / 2
            return {
                (I, J): 1.0, (J, I): 0.0,
                (2, I): pa, (I, 2): 1 - pa,
                (3, I): pb, (I, 3): 1 - pb,
                (J, 2): qa, (2, J): 1 - qa,
                (J, 3): qb, (3, J): 1 - qb,
                (2, 3): t, (3, 2): 1 - t,
            }

        B_plus = compute_B_Li_Rj(4, make_arc_values(sa, sb))
        B_minus = compute_B_Li_Rj(4, make_arc_values(-sa, -sb))

        # Algebraic
        B_alg_plus = 0.5 - da * db / 2 + sa * (t - 0.5) + sb * (0.5 - t) + sa * sb / 2
        B_alg_minus = 0.5 - da * db / 2 - sa * (t - 0.5) - sb * (0.5 - t) + sa * sb / 2

        diff = B_plus - B_minus
        alg_diff = 2 * (t - 0.5) * (sa - sb)

        print(f"  trial {trial}: B(+s)-B(-s) = {diff:.6f}, "
              f"2(t-1/2)(sa-sb) = {alg_diff:.6f}, "
              f"B_num={B_plus:.6f}, B_alg={B_alg_plus:.6f}")

    print(f"\n  CONCLUSION: B is NOT even in (sa, sb) at n=4!")
    print(f"  The previous numerical test must have had a bug.")
    print(f"  Let me check what sigma-invariance actually means...")


def recheck_sigma_definition():
    """
    Recheck: sigma maps p_w -> 1-q_w, q_w -> 1-p_w.
    In (s,d) coords: s_w = 1 - p_w - q_w, d_w = p_w - q_w.
    sigma: p_w' = 1 - q_w, q_w' = 1 - p_w
    s_w' = 1 - p_w' - q_w' = 1 - (1-q_w) - (1-p_w) = p_w + q_w - 1 = -s_w. Good.
    d_w' = p_w' - q_w' = (1-q_w) - (1-p_w) = p_w - q_w = d_w. Good.

    But wait -- sigma also changes T[i][w] = 1 - T[w][i] = 1 - p_w.
    Under sigma: T[i][w]' = 1 - p_w' = 1 - (1-q_w) = q_w = T[j][w].
    And T[w][j] = 1 - T[j][w] = 1 - q_w.
    Under sigma: T[w][j]' = 1 - q_w' = 1 - (1-p_w) = p_w = T[w][i].

    So sigma is: swap i<->j in the interface arcs, but DON'T swap i<->j vertices
    themselves. In other words:
    sigma: T[w][i] <-> T[w][j], T[i][w] <-> T[j][w], T[w][i] -> 1-T[j][w], etc.

    NO wait. Let me reread. sigma: p_w -> 1-q_w means T[w][i] -> 1 - T[j][w].

    For w=a at n=4:
    p_a = T[a][i]. Under sigma: p_a' = 1 - q_a = 1 - T[j][a] = T[a][j].
    q_a = T[j][a]. Under sigma: q_a' = 1 - p_a = 1 - T[a][i] = T[i][a].

    So sigma maps: T[a][i] -> T[a][j], T[j][a] -> T[i][a].
    This is exactly swapping i and j in the interface! And since T[a][b]+T[b][a]=1,
    T[a][j] = 1 - T[j][a] = 1 - q_a, and T[i][a] = 1 - T[a][i] = 1 - p_a.

    So sigma(T) is the tournament where i and j trade their interface arcs.

    B(Li,Rj)(sigma(T)): under sigma, i still has index 0, j index 1. But the
    arc weights have changed. In particular, the Li function uses T[w][i] = p_w,
    which under sigma becomes T[w][j] = 1-q_w. But Li still uses vertex i.

    Hmm, so under sigma, Li(S)(sigma(T)) counts paths on S+{i} ending at i,
    but with T[w][i] replaced by 1-T[j][w]. This is NOT the same as Lj.

    Let me just recompute: B(Li,Rj)(sigma(T)).

    At n=4: Li(empty)=1, Li({a})=p_a'=1-q_a, Li({b})=p_b'=1-q_b,
    Li({a,b}) = t*(p_b'-p_a') + p_a' = t*((1-q_b)-(1-q_a)) + (1-q_a) = t*(q_a-q_b) + 1-q_a

    Rj(W) at sigma: q_a'=1-p_a, q_b'=1-p_b.
    Rj({a,b}) = 1 - q_b' + t*(q_b'-q_a') = 1-(1-p_b) + t*((1-p_b)-(1-p_a))
              = p_b + t*(p_a-p_b)
    Rj({b}) = 1-q_b' = 1-(1-p_b) = p_b
    Rj({a}) = 1-q_a' = p_a
    Rj(empty) = 1

    B(sigma) = Rj({a,b}) - Li({a})*Rj({b}) - Li({b})*Rj({a}) + Li({a,b})
    = [p_b+t*(p_a-p_b)] - (1-q_a)*p_b - (1-q_b)*p_a + [t*(q_a-q_b)+1-q_a]

    = p_b + t*p_a - t*p_b - p_b + q_a*p_b - p_a + q_b*p_a + t*q_a - t*q_b + 1 - q_a
    = t*p_a - t*p_b + q_a*p_b + q_b*p_a - p_a + t*q_a - t*q_b + 1 - q_a

    Compare to B(original):
    B = 1 - q_b - p_b + t*(q_b - q_a + p_b - p_a) + p_a*q_b + p_b*q_a
    = 1 - q_b - p_b + t*q_b - t*q_a + t*p_b - t*p_a + p_a*q_b + p_b*q_a

    B(sigma):
    = 1 - q_a - p_a + t*p_a - t*p_b + t*q_a - t*q_b + p_a*q_b + p_b*q_a

    Compare: B = 1 - q_b - p_b + t*(q_b+p_b-q_a-p_a) + p_a*q_b + p_b*q_a
    B(sigma) = 1 - q_a - p_a + t*(p_a+q_a-p_b-q_b) + p_a*q_b + p_b*q_a

    Note: t*(q_b+p_b-q_a-p_a) vs t*(p_a+q_a-p_b-q_b) = -t*(q_b+p_b-q_a-p_a).

    B - B(sigma) = (1-q_b-p_b) - (1-q_a-p_a) + 2t*(q_b+p_b-q_a-p_a)
    = -(q_b+p_b) + (q_a+p_a) + 2t*(q_b+p_b-q_a-p_a)
    = (1-2t)*[(q_a+p_a) - (q_b+p_b)]
    = (1-2t)*(-s_a - (-s_b))   since s_w = 1-p_w-q_w, so p_w+q_w = 1-s_w
    Wait: p_w+q_w = 1-s_w, so q_a+p_a = 1-s_a, q_b+p_b = 1-s_b
    = (1-2t)*[(1-s_a)-(1-s_b)] = (1-2t)*(s_b-s_a)

    So B != B(sigma) unless s_a = s_b or t = 1/2!

    This contradicts the numerical tests. Let me verify numerically.
    """
    print(f"\n=== Recheck sigma definition at n=4 ===\n")

    I, J = 0, 1

    random.seed(99)
    for trial in range(10):
        pa = random.uniform(0, 1)
        qa = random.uniform(0, 1)
        pb = random.uniform(0, 1)
        qb = random.uniform(0, 1)
        t = random.uniform(0, 1)

        arc_orig = {
            (I, J): 1.0, (J, I): 0.0,
            (2, I): pa, (I, 2): 1 - pa,
            (3, I): pb, (I, 3): 1 - pb,
            (J, 2): qa, (2, J): 1 - qa,
            (J, 3): qb, (3, J): 1 - qb,
            (2, 3): t, (3, 2): 1 - t,
        }

        # Apply sigma: p_w -> 1-q_w, q_w -> 1-p_w
        arc_sigma = dict(arc_orig)
        for w in [2, 3]:
            old_pw = arc_orig[(w, I)]
            old_qw = arc_orig[(J, w)]
            new_pw = 1 - old_qw
            new_qw = 1 - old_pw
            arc_sigma[(w, I)] = new_pw
            arc_sigma[(I, w)] = 1 - new_pw
            arc_sigma[(J, w)] = new_qw
            arc_sigma[(w, J)] = 1 - new_qw

        B_orig = compute_B_Li_Rj(4, arc_orig)
        B_sigma = compute_B_Li_Rj(4, arc_sigma)

        sa = 1 - pa - qa
        sb = 1 - pb - qb
        predicted = (1 - 2 * t) * (sb - sa)

        print(f"  trial {trial}: B-B(sigma)={B_orig-B_sigma:.6f}, "
              f"(1-2t)(sb-sa)={predicted:.6f}, match={abs(B_orig-B_sigma-predicted)<1e-10}")


if __name__ == "__main__":
    # First: verify C_w + D_w at the correct level
    test_Cw_Dw_cancellation(4)
    test_Cw_Dw_cancellation(5)

    # Explicit n=4 formula
    explicit_Cw_formula_n4()

    # Recheck whether sigma-invariance and s-evenness really hold
    recheck_s_evenness_n4()
    recheck_sigma_definition()
