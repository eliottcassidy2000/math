"""Audit check 2 (THM-4708 part 2).
 (i)  m_5 = (eta1/eta5)^6 + 22 + 125 (eta5/eta1)^6 == T_5A + 16 and
      m_7 = (eta1/eta7)^4 + 13 + 49 (eta7/eta1)^4 == T_7A + 9, against the OEIS b-files A007251 / A007264
      (Monster McKay-Thompson series 5A, 7A, constant term 0), to O(q^NQ).
 (ii) sum_n s7(n) m_7^{-n} == (7 E2(q^7) - E2(q))/6 and sum_n C(2n,n) b_n m_5^{-n} == (5 E2(q^5) - E2(q))/4
      as exact q-series to O(q^NG).
 (iii) the stated three-term recurrences, for n < NR, against the binomial-sum definitions.
"""
import sys, os
from math import comb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_seqs import s7, a274786

NQ = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
NG = int(sys.argv[2]) if len(sys.argv) > 2 else 160
NR = int(sys.argv[3]) if len(sys.argv) > 3 else 700
HERE = os.path.dirname(os.path.abspath(__file__))


def mul(a, b, N):
    c = [0] * N
    for i, x in enumerate(a[:N]):
        if x:
            bi = b[:N - i]
            for j, y in enumerate(bi):
                c[i + j] += x * y
    return c


def phi_pow(step, e, N):  # prod_{n>=1} (1-q^{step n})^e
    s = [0] * N
    s[0] = 1
    n = 1
    while step * n < N:
        d = step * n
        if e >= 0:
            for _ in range(e):
                for k in range(N - 1, d - 1, -1):
                    s[k] -= s[k - d]
        else:
            for _ in range(-e):
                for k in range(d, N):
                    s[k] += s[k - d]
        n += 1
    return s


def hauptmodul(N_lev, e, const, coef, N):
    """(eta1/etaN)^e + const + coef (etaN/eta1)^e as list c[k] = coefficient of q^{k-1}, k = 0..N-1.
    Valid since e (N_lev - 1)/24 = 1 for (5,6) and (7,4)."""
    assert e * (N_lev - 1) == 24
    A = mul(phi_pow(1, e, N), phi_pow(N_lev, -e, N), N)        # (eta1/etaN)^e = q^{-1} A
    B = mul(phi_pow(N_lev, e, N), phi_pow(1, -e, N), N)        # (etaN/eta1)^e = q B
    out = [0] * N
    for k in range(N):
        out[k] += A[k]                     # q^{k-1}
    out[1] += const
    for k in range(N - 2):
        out[k + 2] += coef * B[k]          # q^{k+1}
    return out


def read_bfile(fn, N):
    d = {}
    with open(os.path.join(HERE, fn)) as f:
        for line in f:
            n, a = line.split()
            d[int(n)] = int(a)
    return [d[k - 1] for k in range(N)]


def series_inverse(a, N):  # 1/a for a[0] = 1, integer
    assert a[0] == 1
    b = [0] * N
    b[0] = 1
    for n in range(1, N):
        b[n] = -sum(a[k] * b[n - k] for k in range(1, min(n, len(a) - 1) + 1))
    return b


def sigma(n):
    return sum(d for d in range(1, n + 1) if n % d == 0)


def gf_check(N_lev, e, const, coef, u, NG, norm):
    m = hauptmodul(N_lev, e, const, coef, NG + 2)     # m = q^{-1} * (m[0] + m[1] q + ...)
    inv = series_inverse(m, NG + 1)                   # x = 1/m = q * inv
    x = [0] + inv[:NG]                                # x as a power series, length NG+1
    Z = [0] * (NG + 1)
    xp = [1] + [0] * NG
    for n in range(NG + 1):
        for k in range(n, NG + 1):
            Z[k] += u[n] * xp[k]
        xp = mul(xp, x, NG + 1)
    # RHS: (N E2(q^N) - E2(q))/(N-1), E2 = 1 - 24 sum sigma(n) q^n
    rhs = [0] * (NG + 1)
    rhs[0] = 1
    for n in range(1, NG + 1):
        v = 24 * sigma(n)                             # from -E2(q)
        if n % N_lev == 0:
            v -= N_lev * 24 * sigma(n // N_lev)
        assert v % (N_lev - 1) == 0
        rhs[n] = v // (N_lev - 1)
    return Z == rhs, Z[:8], rhs[:8]


if __name__ == '__main__':
    m5 = hauptmodul(5, 6, 22, 125, NQ)
    m7 = hauptmodul(7, 4, 13, 49, NQ)
    t5 = read_bfile('b007251.txt', NQ)
    t7 = read_bfile('b007264.txt', NQ)
    d5 = [a - b for a, b in zip(m5, t5)]
    d7 = [a - b for a, b in zip(m7, t7)]
    print(f'(i) m_5 - T_5A to O(q^{NQ - 1}): nonzero diffs at q^k:', {k - 1: v for k, v in enumerate(d5) if v})
    print(f'    m_7 - T_7A to O(q^{NQ - 1}): nonzero diffs at q^k:', {k - 1: v for k, v in enumerate(d7) if v})
    print('    m_5 =', ' + '.join(f'{c}q^{k - 1}' for k, c in enumerate(m5[:5])), '...')
    print('    m_7 =', ' + '.join(f'{c}q^{k - 1}' for k, c in enumerate(m7[:5])), '...')
    # Fricke-elliptic constant check: other constants do not give 3-term recurrences (not re-checked here)
    U = a274786(max(NG + 1, NR))
    S = s7(max(NG + 1, NR))
    ok7, z7, r7 = gf_check(7, 4, 13, 49, S, NG, 6)
    ok5, z5, r5 = gf_check(5, 6, 22, 125, U, NG, 4)
    print(f'(ii) sum s7(n) m_7^-n == (7E2(q^7)-E2(q))/6 to O(q^{NG + 1}): {ok7}   head {z7} vs {r7}')
    print(f'     sum C(2n,n)b_n m_5^-n == (5E2(q^5)-E2(q))/4 to O(q^{NG + 1}): {ok5}   head {z5} vs {r5}')
    # (iii) recurrences
    bad5 = [n for n in range(1, NR - 1) if (n + 1) ** 3 * U[n + 1] != 2 * (2 * n + 1) * (11 * n * n + 11 * n + 3) * U[n] + 4 * n * (4 * n * n - 1) * U[n - 1]]
    bad7 = [n for n in range(1, NR - 1) if (n + 1) ** 3 * S[n + 1] != (2 * n + 1) * (13 * n * n + 13 * n + 4) * S[n] + 3 * n * (9 * n * n - 1) * S[n - 1]]
    print(f'(iii) U_n = C(2n,n) b_n satisfies (n+1)^3U_(n+1) = 2(2n+1)(11n^2+11n+3)U_n + 4n(4n^2-1)U_(n-1) for 1<=n<{NR - 1}: failures {bad5[:5]}')
    print(f'      s7 satisfies (n+1)^3s_(n+1) = (2n+1)(13n^2+13n+4)s_n + 3n(9n^2-1)s_(n-1) for 1<=n<{NR - 1}: failures {bad7[:5]}')
    # n = 0 step: u_1 = (coefficient at n=0) * u_0
    print('      n=0 step: U_1 =', U[1], '= 2*1*3 =', 2 * 1 * 3, '; s7_1 =', S[1], '= 1*4 =', 4)
