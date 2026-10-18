"""Referee's own sequence generator (audit 2026-10-18).  Written from the definitions in the
audit brief only; nothing is imported from the author's scripts.

All sequences are returned as exact Python ints, u[0..N-1].
"""
from functools import lru_cache
from math import comb, factorial

ZAGIER = {  # (a, b, c):  (n+1)^2 u_{n+1} = (a n^2 + a n + b) u_n - c n^2 u_{n-1}
    'A': (7, 2, -8), 'B': (9, 3, 27), 'C': (10, 3, 9),
    'D': (11, 3, -1), 'E': (12, 4, 32), 'F': (17, 6, 72),
}


def zagier(name, N):
    a, b, c = ZAGIER[name]
    u = [1, b]
    for n in range(1, N - 1):
        num = (a * n * n + a * n + b) * u[n] - c * n * n * u[n - 1]
        q, r = divmod(num, (n + 1) ** 2)
        assert r == 0, (name, n)
        u.append(q)
    return u[:N]


def binsum(N, term):
    return [sum(term(n, k) for k in range(n + 1)) for n in range(N)]


def apery(N):
    return binsum(N, lambda n, k: comb(n, k) ** 2 * comb(n + k, k) ** 2)


def domb(N):
    return binsum(N, lambda n, k: comb(n, k) ** 2 * comb(2 * k, k) * comb(2 * n - 2 * k, n - k))


def s7(N):
    return binsum(N, lambda n, k: comb(n, k) ** 2 * comb(n + k, k) * comb(2 * k, n))


def s10(N):
    return binsum(N, lambda n, k: comb(n, k) ** 4)


def apery_zeta2(N):  # b_n = sum C(n,k)^2 C(n+k,k)
    return binsum(N, lambda n, k: comb(n, k) ** 2 * comb(n + k, k))


def a274786(N):  # U_n = C(2n,n) b_n  (product formula; recurrence checked separately)
    b = apery_zeta2(N)
    return [comb(2 * n, n) * b[n] for n in range(N)]


def hulek_verrill(N, r=5):
    """hv_n = sum over m_1+...+m_r = n of (n!/prod m_i!)^2, from the multinomial definition,
    organised as S_r(n) = sum_k C(n,k)^2 S_{r-1}(n-k), S_1 = 1 (peel off the last part)."""
    S = [1] * N
    for _ in range(r - 1):
        S = [sum(comb(n, k) ** 2 * S[n - k] for k in range(n + 1)) for n in range(N)]
    return S


def hv_bruteforce(n, r=5):
    """Literal multinomial sum (for small n only), to validate hulek_verrill."""
    tot = 0
    def rec(left, parts_left, denom):
        nonlocal tot
        if parts_left == 1:
            d = denom * factorial(left)
            tot += (factorial(n) // d) ** 2
            return
        for m in range(left + 1):
            rec(left - m, parts_left - 1, denom * factorial(m))
    rec(n, r, 1)
    return tot


def hyp333(N):  # ((3k)!/k!^3)^2
    return [(factorial(3 * k) // factorial(k) ** 3) ** 2 for k in range(N)]


def binom2k4(N):  # C(2k,k)^4
    return [comb(2 * k, k) ** 4 for k in range(N)]


def quintic(N):  # (5k)!/k!^5
    return [factorial(5 * k) // factorial(k) ** 5 for k in range(N)]


FAMILIES = {
    # name: (generator, weight)
    'A': (lambda N: zagier('A', N), 1), 'B': (lambda N: zagier('B', N), 1),
    'C': (lambda N: zagier('C', N), 1), 'D': (lambda N: zagier('D', N), 1),
    'E': (lambda N: zagier('E', N), 1), 'F': (lambda N: zagier('F', N), 1),
    'Apery': (apery, 2), 'Domb': (domb, 2), 's7': (s7, 2), 's10': (s10, 2), 'A274786': (a274786, 2),
    'HV': (hulek_verrill, 3), 'H333': (hyp333, 3), 'B2k4': (binom2k4, 3), 'Quintic': (quintic, 3),
}

_cache = {}


def seq(name, N):
    key = name
    if key in _cache and len(_cache[key]) >= N:
        return _cache[key][:N]
    u = FAMILIES[name][0](N)
    _cache[key] = u
    return u


def primes_in(lo, hi):
    out = []
    for n in range(max(2, lo), hi + 1):
        if all(n % d for d in range(2, int(n ** 0.5) + 1)):
            out.append(n)
    return out


if __name__ == '__main__':
    # sanity: closed forms of the Zagier families
    from math import comb as C
    N = 30
    franel = [sum(C(n, k) ** 3 for k in range(n + 1)) for n in range(N)]
    assert zagier('A', N) == franel
    assert zagier('C', N) == [sum(C(n, k) ** 2 * C(2 * k, k) for k in range(n + 1)) for n in range(N)]
    assert zagier('D', N) == apery_zeta2(N)
    assert zagier('E', N) == [sum(C(n, k) * C(2 * k, k) * C(2 * (n - k), n - k) for k in range(n + 1)) for n in range(N)]
    assert zagier('F', N) == [sum((-1) ** k * 8 ** (n - k) * C(n, k) * franel[k] for k in range(n + 1)) for n in range(N)]
    assert zagier('B', N) == [sum((-1) ** k * 3 ** (n - 3 * k) * C(n, 3 * k) * C(3 * k, k) * C(2 * k, k)
                                  for k in range(n // 3 + 1)) for n in range(N)]
    for n in range(9):
        assert hulek_verrill(12)[n] == hv_bruteforce(n), n
    print('Zagier A-F closed forms OK for n<30; HV recursion = brute-force multinomial for n<9')
    for nm in FAMILIES:
        print(nm, seq(nm, 8))
