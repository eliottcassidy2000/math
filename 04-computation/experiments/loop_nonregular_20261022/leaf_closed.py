"""Closed forms of the leaf classes of the linear-forest loop expansion (non-regular tournaments).
Session opus-2026-10-22-S54.  Derived by closed_forms.py from the leaf engine's decomposition weights; checked here
against the engine (leaf_engine.py), which is checked against brute force (leaf_engine_selftest.out).

Notation: Shat = S/n, G = -Shat^2 (I + Shat^2)^-1, H = Shat^2 (I - Shat^2)^-1, O = Shat (I - Shat^2)^-1,
Q = Shat (I + Shat^2)^-1, d = diag(H - G), g = G 1, q = Q 1 = (I + Shat^2)^-1 delta, delta = Shat 1 (score deviations / n).
Entrywise products are written u*v; M^(2) is the entrywise square.

  X(-1,2) = -1^T G 1 = -delta^T (I + Shat^2)^-1 delta                         (path trees)
  X(0,1)  = -2 d^T G 1                                                        (lollipop; dumbbell -d^T G d with a leaf for a loop)
  X(-1,3) = -(2/3) sum_w (g_w^3 + 3 g_w q_w^2)                                 (stars)
  X(0,2)  = sum_w [q_w^2 (3G - H)_ww - 3 g_w^2 d_w]                            (loop + two leaves at one vertex)
          - 2 (g*g + q*q)^T G d - 4 (g*q)^T Q d                                (two leaves, bridge, loop)
          + 2 sum_(u,v) [g_u g_v (G^(2) - H^(2) + O^(2) + Q^(2))_uv + q_u q_v (G^(2) - Q^(2))_uv]   (cut theta)
The dumbbell family closes into a square:  X_dumbbell + X(0,1) + X(-1,2) = -(d + 1)^T G (d + 1).
Usage: python leaf_closed.py"""
import sys, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def mats(Sh):
    n = Sh.shape[0]; I = np.eye(n); S2 = Sh @ Sh
    G = -S2 @ np.linalg.inv(I + S2); H = S2 @ np.linalg.inv(I - S2)
    O = Sh @ np.linalg.inv(I - S2); Q = Sh @ np.linalg.inv(I + S2)
    one = np.ones(n)
    return G, H, O, Q, np.diag(H - G), G @ one, Q @ one, one


def X_tree2(Sh):
    G, H, O, Q, d, g, q, one = mats(Sh)
    return -one @ g


def X_lollipop(Sh):
    G, H, O, Q, d, g, q, one = mats(Sh)
    return -2 * d @ g


def X_star(Sh):
    G, H, O, Q, d, g, q, one = mats(Sh)
    return -(2 / 3) * np.sum(g ** 3 + 3 * g * q ** 2)


def X_02(Sh):
    G, H, O, Q, d, g, q, one = mats(Sh)
    a = np.sum(q ** 2 * np.diag(3 * G - H) - 3 * g ** 2 * d)
    b = -2 * (g * g + q * q) @ G @ d - 4 * (g * q) @ Q @ d
    c = 2 * (g @ (G * G - H * H + O * O + Q * Q) @ g + q @ (G * G - Q * Q) @ q)
    return a + b + c


def X_dumbbell(Sh):
    G, H, O, Q, d, g, q, one = mats(Sh)
    return -d @ G @ d


if __name__ == '__main__':
    import leaf_engine as LE
    from cluster_brute_classes import almost_regular
    rng = np.random.default_rng(21)
    worst = 0.0
    for n in (8, 10, 12, 16, 20):
        A = almost_regular(n, rng)
        Sh = (A - A.T).astype(float) / n
        checks = [((-1, 2), X_tree2(Sh)), ((0, 1), X_lollipop(Sh)), ((-1, 3), X_star(Sh)), ((0, 2), X_02(Sh))]
        line = []
        for c, v in checks:
            e = LE.X(Sh, c[0], c[1]).real
            worst = max(worst, abs(v - e))
            line.append(f'{c}: {v:+.12e} (engine diff {v - e:+.1e})')
        G, H, O, Q, d, g, q, one = mats(Sh)
        sq = -(d + one) @ G @ (d + one)
        fam = X_dumbbell(Sh) + X_lollipop(Sh) + X_tree2(Sh)
        line.append(f'square identity diff {sq - fam:+.1e}')
        print(f'n={n}: ' + '; '.join(line), flush=True)
    print(f'max |closed - engine| = {worst:.2e}')
