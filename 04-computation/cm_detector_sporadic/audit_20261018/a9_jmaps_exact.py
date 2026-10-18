"""Audit: the level-5/7 j-maps of THM-4708 part 3 as exact q-series identities (integer arithmetic, to O(q^NN)).
 t = (eta1/eta7)^4:  t^7 j(tau) = (t^2+13t+49)(t^2+245t+2401)^3,   t j(7tau) = (t^2+13t+49)(t^2+5t+1)^3
 s = (eta1/eta5)^6:  s^5 j(tau) = (s^2+250s+3125)^3,                s j(5tau) = (s^2+10s+5)^3
Both sides are modular functions on X_0(N) whose only possible pole is at infinity, of order <= N+1, so agreement of
the q-expansions through q^0 already proves each identity (valence formula); we check far beyond that."""
import sys
NN = int(sys.argv[1]) if len(sys.argv) > 1 else 120


def mul(a, b, N):
    c = [0] * N
    for i, x in enumerate(a[:N]):
        if x:
            for j, y in enumerate(b[:N - i]):
                c[i + j] += x * y
    return c


def phi_pow(step, e, N):
    s = [0] * N
    s[0] = 1
    n = 1
    while step * n < N:
        d = step * n
        for _ in range(abs(e)):
            if e > 0:
                for k in range(N - 1, d - 1, -1):
                    s[k] -= s[k - d]
            else:
                for k in range(d, N):
                    s[k] += s[k - d]
        n += 1
    return s


def sigma3(n):
    return sum(d ** 3 for d in range(1, n + 1) if n % d == 0)


def jq(N, step=1):
    """q^{step} j(q^step) as power series (j = E4^3/Delta), length N."""
    E4 = [1] + [240 * sigma3(n) for n in range(1, N)]
    E43 = mul(mul(E4, E4, N), E4, N)
    invphi24 = phi_pow(1, -24, N)          # 1/prod(1-q^n)^24 ; Delta = q prod^24
    J = mul(E43, invphi24, N)              # q j(q)
    if step == 1:
        return J
    out = [0] * N
    for k in range(N):
        if k * step < N:
            out[k * step] = J[k]
    return out                              # q^step j(q^step)


def laurent_poly_eval(coeffs_t, tser, shift, N):
    """sum_c coeffs_t[c] t^c with t = q^{-1} * tser: returns (series, shift) where value = q^{-deg} * series."""
    D = len(coeffs_t) - 1
    out = [0] * N
    tp = [1] + [0] * (N - 1)
    pw = []
    for c in range(D + 1):
        pw.append(tp)
        tp = mul(tp, tser, N)
    for c, a in enumerate(coeffs_t):
        if a:
            # t^c = q^{-c} pw[c] = q^{-D} q^{D-c} pw[c]
            for k in range(N - (D - c)):
                out[k + D - c] += a * pw[c][k]
    return out  # value = q^{-D} * out


def polymul(a, b):
    c = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i + j] += x * y
    return c


def cube(a):
    return polymul(polymul(a, a), a)


N = NN
for (lev, e, F1, F2, pj1, pj7) in [
        (7, 4, polymul([49, 13, 1], cube([2401, 245, 1])), polymul([49, 13, 1], cube([1, 5, 1])), 7, 1),
        (5, 6, cube([3125, 250, 1]), cube([5, 10, 1]), 5, 1)]:
    tser = mul(phi_pow(1, e, N), phi_pow(lev, -e, N), N)   # t = q^{-1} tser
    # identity 1: t^{pj1} j(tau) = F1(t).  LHS = q^{-pj1} tser^{pj1} * q^{-1} (q j) ; degree of F1 = pj1 + 1
    tp = [1] + [0] * (N - 1)
    for _ in range(pj1):
        tp = mul(tp, tser, N)
    lhs1 = mul(tp, jq(N), N)                                  # value = q^{-(pj1+1)} lhs1
    rhs1 = laurent_poly_eval(F1, tser, 0, N)                  # value = q^{-deg F1} rhs1, deg F1 = pj1 + 1
    assert len(F1) - 1 == pj1 + 1
    ok1 = lhs1 == rhs1
    # identity 2: t^{pj7} j(N tau) = F2(t); LHS = q^{-1} tser * q^{-lev} (q^lev j(q^lev)); deg F2 = lev + 1
    lhs2 = mul(tser, jq(N, lev), N)                           # value = q^{-(1+lev)} lhs2
    rhs2 = laurent_poly_eval(F2, tser, 0, N)
    assert len(F2) - 1 == lev + 1
    ok2 = lhs2 == rhs2
    print(f'level {lev}: j(tau) identity exact to O(q^{N - (pj1 + 1)}): {ok1};  j({lev} tau) identity exact to O(q^{N - lev - 1}): {ok2}')
