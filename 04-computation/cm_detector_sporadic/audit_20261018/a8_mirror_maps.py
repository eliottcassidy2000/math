"""Audit of THM-4708 'Reading': mirror maps of Zagier B and D from the Frobenius log solution.
For (n+1)^2 u_{n+1} = (a n^2 + a n + b) u_n - c n^2 u_{n-1}: u_n(eps) from the recurrence at n+eps, v_n = u_n'(0);
y0 = sum u_n t^n, y1 = y0 log t + sum v_n t^n,  q = t exp(sum v_n t^n / y0);  t(q) by reversion (exact rationals).
Claims: t_B(q) = -P_9(-q), P_9 = (eta(9tau)/eta(tau))^3;   t_D(q) = q prod (1-q^n)^{5 (n/5)} = (chi_0/chi_{-1/5})^5."""
from fractions import Fraction as Fr
import sys

NT = int(sys.argv[1]) if len(sys.argv) > 1 else 40
ZAG = {'B': (9, 3, 27), 'D': (11, 3, -1)}


def mul(a, b, N):
    c = [Fr(0)] * N
    for i, x in enumerate(a[:N]):
        if x:
            for j, y in enumerate(b[:N - i]):
                c[i + j] += x * y
    return c


def exp_series(f, N):  # f[0] = 0
    e = [Fr(0)] * N
    e[0] = Fr(1)
    # e' = f' e
    df = [k * f[k] for k in range(N)]
    for n in range(1, N):
        e[n] = sum(df[k] * e[n - k] for k in range(1, n + 1)) / n
    return e


def inverse(a, N):
    b = [Fr(0)] * N
    b[0] = 1 / a[0]
    for n in range(1, N):
        b[n] = -sum(a[k] * b[n - k] for k in range(1, n + 1)) / a[0]
    return b


def revert(q_of_t, N):
    """q = t * (1 + ...): return t(q) as series in q with t = q * (...)."""
    # Lagrange-free iterative: t = q / h(t), h = q/t
    h = q_of_t[1:] + [Fr(0)]           # q/t as series in t
    t = [Fr(0)] * N
    t[1] = Fr(1)
    for _ in range(N):
        # compute h(t(q))
        ht = [Fr(0)] * N
        tp = [Fr(1)] + [Fr(0)] * (N - 1)
        for k in range(N):
            if h[k]:
                for i in range(N):
                    ht[i] += h[k] * tp[i]
            tp = mul(tp, t, N)
        inv = inverse(ht, N)
        tn = [Fr(0)] + inv[:N - 1]
        if tn == t:
            break
        t = tn
    return t


def mirror(name, N):
    a, b, c = ZAG[name]
    u = [Fr(1)]
    v = [Fr(0)]
    um1 = vm1 = Fr(0)
    for n in range(0, N):
        # (n+eps+1)^2 U_{n+1} = (a(n+eps)^2 + a(n+eps) + b) U_n - c (n+eps)^2 U_{n-1}, mod eps^2
        A0, A1 = Fr((n + 1) ** 2), Fr(2 * (n + 1))
        B0, B1 = Fr(a * n * n + a * n + b), Fr(2 * a * n + a)
        C0, C1 = Fr(c * n * n), Fr(2 * c * n)
        un, vn = u[n], v[n]
        up, vp = (u[n - 1], v[n - 1]) if n >= 1 else (um1, vm1)
        R0 = B0 * un - C0 * up
        R1 = B0 * vn + B1 * un - C0 * vp - C1 * up
        u1 = R0 / A0
        v1 = (R1 - A1 * u1) / A0
        u.append(u1)
        v.append(v1)
    u, v = u[:N], v[:N]
    ratio = mul(v, inverse(u, N), N)
    E = exp_series(ratio, N)
    q_of_t = [Fr(0)] + E[:N - 1]
    return revert(q_of_t, N), u


def prod_series(exps, N):  # prod_n (1 - q^n)^{exps(n)}
    s = [Fr(0)] * N
    s[0] = Fr(1)
    for n in range(1, N):
        e = exps(n)
        if e == 0:
            continue
        # multiply by (1-q^n)^e via binomial series
        f = [Fr(0)] * N
        k = 0
        coef = Fr(1)
        while n * k < N:
            f[n * k] = coef * (-1) ** k
            coef = coef * (e - k) / (k + 1)
            k += 1
        s = mul(s, f, N)
    return s


def leg5(n):
    r = n % 5
    return 0 if r == 0 else (1 if r in (1, 4) else -1)


if __name__ == '__main__':
    N = NT
    tB, uB = mirror('B', N)
    # P9(q) = q prod (1-q^{9n})^3 / (1-q^n)^3 ; -P9(-q)
    base = prod_series(lambda n: (3 if n % 9 == 0 else 0) - 3, N)  # prod (1-q^{9n})^3 (1-q^n)^{-3}
    P9 = [Fr(0)] + base[:N - 1]
    claimB = [-(P9[k] * (-1) ** k) for k in range(N)]
    print('Zagier B: u =', [int(x) for x in uB[:8]])
    print('  t_B(q) =', [str(x) for x in tB[:8]])
    print('  -P_9(-q) =', [str(x) for x in claimB[:8]])
    print(f'  t_B(q) == -P_9(-q) to O(q^{N}):', tB == claimB)
    tD, uD = mirror('D', N)
    claimD = [Fr(0)] + prod_series(lambda n: 5 * leg5(n), N)[:N - 1]
    print('Zagier D: t_D(q) =', [str(x) for x in tD[:8]])
    print(f'  q prod (1-q^n)^(5(n/5)) to O(q^{N}):', tD == claimD)
