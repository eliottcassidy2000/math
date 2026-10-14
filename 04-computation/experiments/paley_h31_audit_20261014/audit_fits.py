#!/usr/bin/env python3
"""AUDIT of the R(q) fits, holdouts, q = 27 'universality' control and the pre-registered R(43), R(47) intervals
(paley_R_asymptotic_fit_20261014.py), recomputed in 50-digit arithmetic, plus the effect of the leading arithmetic
(W_4 / Apery-form) term R/e -= 7 t_q / q^6 found by the referee (t_q = a_q for primes, alpha^3 + beta^3 = 260 for q = 27).
"""
from mpmath import mp, mpf, e as E, matrix, lu_solve, factorial
mp.dps = 50

H = {7: 189, 11: 95095, 19: 1172695746915, 23: 15760206976379349,
     27: 418425945688463432031, 31: 19876729650849384218015865}
R = {q: mpf(H[q]) * mpf(2) ** (q - 1) / factorial(q) for q in H}
T = {7: 24, 11: -44, 19: 44, 23: -56, 27: 260, 31: -160, 43: 52, 47: 528}   # Frobenius traces (eta(2z)^4 eta(4z)^4)
PR = [7, 11, 19, 23, 31]


def fit(points, Rd, c1=None):
    m = len(points)
    if c1 is None:
        A = matrix([[mpf(1) / mpf(q) ** j for j in range(1, m + 1)] for q in points])
        b = matrix([1 - Rd[q] / E for q in points])
        return [lu_solve(A, b)[i] for i in range(m)]
    A = matrix([[mpf(1) / mpf(q) ** j for j in range(2, m + 2)] for q in points])
    b = matrix([1 - Rd[q] / E - mpf(c1) / q for q in points])
    return [mpf(c1)] + [lu_solve(A, b)[i] for i in range(m)]


def pred(c, q):
    return E * (1 - sum(cj / mpf(q) ** (j + 1) for j, cj in enumerate(c)))


def report(Rd, label):
    print(f'\n===== {label}')
    for q in sorted(Rd):
        print(f'   q = {q}: R = {mp.nstr(Rd[q], 12)};  q^2 (1 - R/e - 3/(2q)) = {mp.nstr(q * q * (1 - Rd[q] / E - mpf(3) / (2 * q)), 7)}')
    for target in (23, 31):
        below = [q for q in PR if q < target]
        for m in (1, 2, 3):
            c = fit(below[-m:], Rd)
            print(f'   holdout target {target}, order {m} on {below[-m:]}: error {mp.nstr(pred(c, target) - Rd[target], 3)}')
    for m in (2, 3):
        pts = [19, 23, 31][-m:]
        c = fit(pts, Rd)
        print(f'   q = 27 control, order {m} on {pts}: predicted - exact = {mp.nstr(pred(c, 27) - Rd[27], 4)}')
    # c1 = 3/2 fixed: holdouts at 31 and 23, and at 27
    c = fit([23], Rd, c1=1.5)
    print(f'   c1 fixed, c2 from 23 -> 31 error {mp.nstr(pred(c, 31) - Rd[31], 3)}; -> 27 error {mp.nstr(pred(c, 27) - Rd[27], 3)}')
    c = fit([19, 23], Rd, c1=1.5)
    print(f'   c1 fixed, (c2,c3) from 19,23 -> 31 error {mp.nstr(pred(c, 31) - Rd[31], 3)}; -> 27 error {mp.nstr(pred(c, 27) - Rd[27], 3)}')
    c = fit([19], Rd, c1=1.5)
    print(f'   c1 fixed, c2 from 19 -> 23 error {mp.nstr(pred(c, 23) - Rd[23], 3)}')
    c = fit([19, 23, 31], Rd, c1=1.5)
    print(f'   c1 fixed, (c2,c3,c4) from 19,23,31: c = {[mp.nstr(x, 6) for x in c]}; -> 27 error {mp.nstr(pred(c, 27) - Rd[27], 3)}')
    # free fits for c1
    for pts in ([23, 31], [19, 23, 31], [11, 19, 23, 31]):
        c = fit(pts, Rd)
        print(f'   free fit on {pts}: c1 = {mp.nstr(c[0], 7)}, c2 = {mp.nstr(c[1], 6)}')


def author_intervals(Rd):
    # reproduce the author's construction exactly
    preds = {43: [], 47: []}
    for m in (1, 2, 3, 4):
        c = fit(PR[-m:], Rd)
        for q in (43, 47):
            preds[q].append(pred(c, q))
    h31 = [abs(pred(fit([19, 23][-m:] if m <= 2 else [11, 19, 23], Rd), 31) - Rd[31]) for m in (2, 3)]
    pad = max(h31)
    out = {}
    for q in (43, 47):
        lo = min(preds[q][1], preds[q][2]) - pad * 31 / q
        hi = max(preds[q][1], preds[q][2]) + pad * 31 / q
        out[q] = (lo, hi)
        print(f'   model-agnostic R({q}) interval (unrounded): [{mp.nstr(lo, 9)}, {mp.nstr(hi, 9)}]')
    c2a = (1 - Rd[31] / E - mpf(1.5) / 31) * 31 ** 2
    cb = fit([23, 31], Rd, c1=1.5)
    p31a = pred(fit([23], Rd, c1=1.5), 31)
    p31b = pred(fit([19, 23], Rd, c1=1.5), 31)
    for q in (43, 47):
        va = E * (1 - mpf(1.5) / q - c2a / mpf(q) ** 2)
        vb = pred(cb, q)
        pad = max(abs(p31a - Rd[31]), abs(p31b - Rd[31])) * (mpf(31) / q) ** 2
        lo, hi = min(va, vb) - pad, max(va, vb) + pad
        print(f'   c1-fixed R({q}) interval (unrounded): [{mp.nstr(lo, 9)}, {mp.nstr(hi, 9)}]  centres {mp.nstr(va, 9)}, {mp.nstr(vb, 9)}')
    return out


def main():
    report(R, 'exact R(q), no arithmetic correction')
    print('\n   author interval construction reproduced:')
    author_intervals(R)
    Rc = {q: R[q] + 7 * E * T[q] / mpf(q) ** 6 for q in R}
    print('\n   leading W_4 term -7 e t_q/q^6 in R at each q:', {q: mp.nstr(-7 * E * T[q] / mpf(q) ** 6, 3) for q in sorted(T)})
    report(Rc, 'R(q) with the leading W_4 term removed (R + 7 e t_q / q^6)')


if __name__ == '__main__':
    main()
