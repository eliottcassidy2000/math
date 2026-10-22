"""Audit check (S53): THM-4759 (1) for D5, D6, E6, E7, E8 by direct enumeration (independent of principal_cores.py).
Root lattice Q with Gram = Cartan matrix A; rho = sum c_i alpha_i with A c = 1, so x.rho = ht(x).
  shifted:  sum_{x in Q} q^{|x + rho/h|^2/2} = eta(tau)^l eta_pi(tau/h):  coefficient of q^{k/h} after the common
            factor q^{|rho/h|^2/2} is #{x : h x^2/2 + ht(x) = k}, against prod (1-p^{hn})^l prod_a (1-p^{an})^{m_a}, p = q^{1/h};
  twined:   sum_{x in Q} zeta_h^{ht x} q^{x^2/2} = eta^l / eta_pi."""
import cmath, math
import numpy as np
from fractions import Fraction
from audit_common import shell_vecs


def cartan(t, n):
    A = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        A[i, i] = 2
    if t == 'D':
        for i in range(n - 2):
            A[i, i + 1] = A[i + 1, i] = -1
        A[n - 3, n - 1] = A[n - 1, n - 3] = -1
    else:   # E_n, Bourbaki: 1-3-4-5-..., 2 attached to 4
        edges = [(0, 2), (2, 3), (3, 4)] + [(i, i + 1) for i in range(4, n - 1)] + [(1, 3)]
        for i, j in edges:
            A[i, j] = A[j, i] = -1
    return A


FR = {('D', 5): {8: 1, 2: 1, 4: -1, 1: -1}, ('D', 6): {10: 1, 2: 1, 5: -1, 1: -1},
      ('E', 6): {12: 1, 3: 1, 2: 1, 6: -1, 4: -1, 1: -1}, ('E', 7): {18: 1, 3: 1, 2: 1, 9: -1, 6: -1, 1: -1},
      ('E', 8): {30: 1, 5: 1, 3: 1, 2: 1, 15: -1, 10: -1, 6: -1, 1: -1}}
H = {('D', 5): 8, ('D', 6): 10, ('E', 6): 12, ('E', 7): 18, ('E', 8): 30}


def prod_series(factors, L):
    c = [0] * (L + 1); c[0] = 1
    for step0, e in factors:
        for nn in range(1, L // step0 + 1):
            step = step0 * nn
            for _ in range(abs(e)):
                if e > 0:
                    for k in range(L, step - 1, -1):
                        c[k] -= c[k - step]
                else:
                    for k in range(step, L + 1):
                        c[k] += c[k - step]
    return c


allok = True
for (t, n), pi in FR.items():
    A = cartan(t, n)
    h = H[(t, n)]
    l = n
    Ainv = [[Fraction(x).limit_denominator(1000) for x in row] for row in np.linalg.inv(A.astype(float))]
    c = [sum(Ainv[i][j] for j in range(n)) for i in range(n)]
    Lden = 1
    for x in c:
        Lden = Lden * x.denominator // math.gcd(Lden, x.denominator)
    Yd = np.array([int(x * Lden) for x in c], dtype=np.int64)     # centre rho/h = (Yd/Lden)/h
    rho2 = sum(c[i] * sum(A[i, j] * c[j] for j in range(n)) for i in range(n))
    assert rho2 == Fraction(l * h * (h + 1), 12)            # strange formula
    K = 3 * h                                              # compare coefficients of p^k, k < K
    R2 = float(rho2) / h ** 2 + 2.0 * K / h + 1e-6
    vals, vecs = shell_vecs(A.astype(np.float64), A, Yd, h * Lden, R2)
    cnt = [0] * K
    for a in vecs:
        a = np.asarray(a, dtype=np.int64)
        k = int(a @ A @ a) * h // 2 + int(a.sum())        # h x^2/2 + ht(x)
        if 0 <= k < K:
            cnt[k] += 1
    rhs = prod_series([(h, l)] + [(a_, m_) for a_, m_ in pi.items()], K - 1)
    ok1 = cnt == rhs
    # twined: sum zeta^{ht x} q^{x^2/2} vs eta^l / eta_pi = prod (1-q^n)^l / prod_a (1-q^{an})^{m_a} (q^0 normalised)
    M = 6
    vals2, vecs2 = shell_vecs(A.astype(np.float64), A, np.zeros(n, dtype=np.int64), 1, 2.0 * M + 1e-6)
    tw = [0j] * (M + 1)
    for a in vecs2:
        a = np.asarray(a, dtype=np.int64)
        tw[int(a @ A @ a) // 2] += cmath.exp(2j * math.pi * int(a.sum()) / h)
    rhs2 = prod_series([(1, l)] + [(a_, -m_) for a_, m_ in pi.items()], M)
    ok2 = all(abs(tw[k] - rhs2[k]) < 1e-8 for k in range(M + 1))
    allok &= ok1 and ok2
    print(f'{"ok  " if ok1 and ok2 else "FAIL"} {t}{n} (h={h}): shifted theta = eta^l eta_pi(tau/h) through p^{K - 1}: {ok1} '
          f'{cnt[:12]}...; twined theta = eta^l/eta_pi through q^{M}: {ok2} {[round(z.real) for z in tw]}')
print('ALL OK' if allok else 'FAIL')
