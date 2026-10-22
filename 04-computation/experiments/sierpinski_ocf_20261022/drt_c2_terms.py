# Term-by-term check of the doubly regular second order: n^2 times each term of Y_2 for Paley P_q, against the exact
# limits X_1^2/2 -> 49/8, X_1 E_0 a/n -> -21/2, delta/n -> -49, E_0 a^2/2 -> 21, E_0 b -> 29/6, X_2 -> 29, so that
# n^2 Y_2 -> 35/24 and c_2 = 5/3 + 35/24 = 25/8.
import sys, math
import numpy as np
sys.path.insert(0, '.')
import second_order as SO
import loop_engine_fast as LEF
for q in [int(a) for a in sys.argv[1:]]:
    sq = {(t * t) % q for t in range(1, q)}
    A = np.array([[1.0 if (j - i) % q in sq else 0.0 for j in range(q)] for i in range(q)])
    S = A - A.T; n = q; Sh = S / n
    x = np.linalg.eigvals(Sh).imag
    D0 = SO.cumulants(lambda z: SO.X0_of(x, z)); D1 = SO.cumulants(lambda z: SO.X1_of(Sh.astype(complex), z))
    k1, k2, k3, k4 = D0[1], D0[2], D0[3], D0[4]
    EL = k1; EL2 = k2 + k1 ** 2; EL3 = k3 + 3 * k2 * k1 + k1 ** 3; EL4 = k4 + 4 * k3 * k1 + 3 * k2 ** 2 + 6 * k2 * k1 ** 2 + k1 ** 4
    Ea = (EL2 - EL) / 2; Ea2 = (EL4 - 2 * EL3 + EL2) / 4; Eb = (2 * EL3 - 3 * EL2 + EL) / 12
    delta = (D1[2] + 2 * k1 * D1[1] - D1[1]) / 2
    X1 = D1[0]; X2 = float(np.real(LEF.XB(Sh, 2)))
    terms = {'X2': X2, 'X1^2/2': X1 ** 2 / 2, 'X1 Ea/n': X1 * Ea / n, 'delta/n': delta / n, 'Ea2/2/n^2': Ea2 / 2 / n ** 2, 'Eb/n^2': Eb / n ** 2}
    Y2 = sum(terms.values())
    print(f"q={q}: " + '  '.join(f"{k}={n * n * v:+.5f}" for k, v in terms.items()) + f"   n^2 Y_2 = {n * n * Y2:+.6f}", flush=True)
print('exact limits: X2 29, X1^2/2 6.125, X1 Ea/n -10.5, delta/n -49, Ea2/2 21, Eb 4.83333; n^2 Y_2 -> 35/24 =', 35 / 24)
