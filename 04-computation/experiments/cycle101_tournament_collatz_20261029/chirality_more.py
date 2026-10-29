"""Chirality of the Paley run-observer determinant beta = det B_1 in Q(sqrt(-q)) (THM-4848 (iv), HYP-9446 (3)).

beta = a + b*eta with eta = sum_{j in QR} zeta^j = (-1 + sqrt(-q))/2 (q = 3 mod 4, every QR < q - 1, so a and b are
the constant and linear coefficients of det B_1 reduced mod Phi_q).  For each prime p of N(beta) = a^2 - ab + ((q+1)/4) b^2
that splits in K = Q(sqrt(-q)), compute v_P(beta) and v_Pbar(beta) for the two primes P = (p, eta - r), Pbar = (p, eta - r')
above p (r, r' the Hensel-lifted roots of X^2 + X + (q+1)/4).  "One-sided" means one of the two valuations is 0.
Pre-registered before computing (2026-10-29): beta is one-sided at every split prime p not dividing 2q
(the C5-C6 pattern "one cubic, zero of its dual" of CLOCK_MINIMAL_PAIRING for P7, k = 3).
This copy runs extra cases given on the command line as q:k (default 7:4), after the audit-time checkpoint.
"""
import sys, itertools, json
from math import gcd
sys.path.insert(0, '../growth_gf_runs_fermat_20261029')
from fast_runs import m2_fast
from runs import circulant, run_matrix, gp_batch
from sympy import factorint, legendre_symbol


def beta_coeffs(q, k):
    QR = sorted({(x * x) % q for x in range(1, q)})
    A = circulant(q, set(QR))
    if k == 2:
        T, _ = m2_fast(A)
        tuples = [(a, b) for a in range(q) for b in range(q) if a != b]
        M = {(s, t): int(T[s[0], s[1], t[0], t[1]]) for s in tuples for t in tuples}
    else:
        tuples, Ml = run_matrix(A, k, list(itertools.permutations(range(q), k)))
        M = {(s, t): Ml[i][j] for i, s in enumerate(tuples) for j, t in enumerate(tuples)}
    shift = lambda t, s: tuple((x + s) % q for x in t)
    reps, seen = [], set()
    for t in tuples:
        if t not in seen:
            seen.update(shift(t, s) for s in range(q))
            reps.append(t)
    def entry(r, rp):
        terms = ['%d*z^%d' % (M[(r, shift(rp, s))], s) for s in range(q) if M[(r, shift(rp, s))]]
        return '(' + ('+'.join(terms) if terms else '0') + ')'
    B1 = 'Mat([' + ';'.join(','.join(entry(r, rp) for rp in reps) for r in reps) + '])'
    out = gp_batch(['my(d=lift(Mod(matdet(%s), polcyclo(%d, z)))); [polcoeff(d, 0), polcoeff(d, 1), d - polcoeff(d,0) - polcoeff(d,1)*(%s)]'
                    % (B1, q, '+'.join('z^%d' % j for j in QR))], timeout=7200)[0]
    a, b, rest = out.strip('[]').split(',', 2)
    assert rest.strip() == '0', 'det B_1 is not a + b*eta'
    return int(a), int(b)


def val_P(a, b, p, r_fn, e_max):
    v = 0
    for j in range(1, e_max + 1):
        r = r_fn(j)
        if (a + b * r) % p ** j == 0:
            v = j
        else:
            break
    return v


def roots_lifted(p, c, j):
    """the two roots of X^2 + X + c modulo p^j, Hensel-lifted from the roots mod p (p odd or p = 2, disc unit)."""
    rs = [x for x in range(p) if (x * x + x + c) % p == 0]
    out = []
    for r in rs:
        m = p
        while m < p ** j:
            m2 = min(m * p, p ** j)
            # Newton step modulo m2
            f, fp = r * r + r + c, 2 * r + 1
            r = (r - f * pow(fp, -1, m2)) % m2
            m = m2
        out.append(r % p ** j)
    return out


report = {}
cases = [tuple(int(x) for x in a.split(':')) for a in sys.argv[1:]] or [(7, 4)]
for q, k in cases:
    a, b = beta_coeffs(q, k)
    c = (q + 1) // 4
    N = a * a - a * b + c * b * b
    g = gcd(a, b)
    fac = factorint(N, limit=10 ** 7)
    lines = []
    for p, e in sorted(fac.items()):
        if p > 10 ** 7 or p == q:
            continue
        split = (p == 2 and (-q) % 8 == 1) or (p != 2 and legendre_symbol((-q) % p, p) == 1)
        if not split:
            lines.append('p=%d inert, exponent %d in N(beta), v_p(gcd(a,b)) = %d' % (p, e, factorint(g).get(p, 0)))
            continue
        vs = [val_P(a, b, p, lambda j, idx=idx: roots_lifted(p, c, j)[idx], e) for idx in (0, 1)]
        lines.append('p=%d split, exponent %d: (v_P, v_Pbar) = %s -> %s' % (p, e, tuple(vs), 'one-sided' if min(vs) == 0 else 'TWO-SIDED'))
    big = [p for p in fac if p > 10 ** 7]
    print('P%d, k=%d: beta = %d + (%d) eta ; gcd(a,b) = %s ; unfactored part > 10^7: %s' % (q, k, a, b, factorint(g), [len(str(p)) for p in big]))
    for ln in lines:
        print('   ' + ln)
