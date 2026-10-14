#!/usr/bin/env python3
"""The rate of R(q) = H(P_q) 2^(q-1)/q! -> e (THM-438) after the exact H(P_31): fits, holdouts, and a pre-registered
second-generation prediction for R(43), R(47).

Ansatz (THM-438 addendum; HYP-2371): 1 - R(q)/e = c_1/q + c_2/q^2 + ... + c_m/q^m.  For each model order m we fit
the m largest exact PRIME points (q = 7, 11, 19, 23, 31; q = 27 = 3^3 is kept out of every fit and used as a control),
predict the next point, and report holdout errors:
  * the June prediction (fit through q <= 23) against the exact R(31);
  * every order's holdout at q = 31 (fit on q <= 23) and at q = 23 (fit on q <= 19);
  * the q = 27 control (prime power, not used in any fit).
The prediction interval for R(43), R(47) is [min, max] over the stable orders, widened by the largest holdout error
seen at the last two holdout points (q = 23 and 31), scaled by (31/q)^(m) for the extrapolation distance -- stated
explicitly so that the interval is falsifiable rather than fitted after the fact.
Session opus-2026-10-14-S34.  Exact inputs from paley_H31_exact_20261014.out.
"""
import math
from fractions import Fraction

E = math.e
H = {7: 189, 11: 95095, 19: 1172695746915, 23: 15760206976379349,
     27: 418425945688463432031, 31: 19876729650849384218015865}
R = {q: H[q] * 2 ** (q - 1) / math.factorial(q) for q in H}
PRIMES = [7, 11, 19, 23, 31]


def solve(A, b):
    n = len(A)
    M = [list(map(Fraction, A[i])) + [Fraction(b[i])] for i in range(n)]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        for r in range(n):
            if r != c and M[r][c] != 0:
                f = M[r][c] / M[c][c]
                M[r] = [x - f * y for x, y in zip(M[r], M[c])]
    return [float(M[i][n] / M[i][i]) for i in range(n)]


def fit(points):
    """exact interpolation of 1 - R/e = sum_{j=1..m} c_j q^-j through m points"""
    m = len(points)
    A = [[Fraction(1, q ** j) for j in range(1, m + 1)] for q in points]
    b = [Fraction(1) - Fraction(R[q] / E).limit_denominator(10 ** 15) for q in points]
    return solve(A, b)


def predict(c, q):
    return E * (1 - sum(cj / q ** (j + 1) for j, cj in enumerate(c)))


def main():
    print('exact R(q) and q (1 - R/e):')
    for q in sorted(R):
        print(f'   q = {q:2d}: R = {R[q]:.9f}   q(1 - R/e) = {q * (1 - R[q] / E):.6f}' + ('   (prime power; control only)' if q == 27 else ''))
    print('\nholdouts (fit the m largest primes below the target, predict the target):')
    worst = {}
    for target in (23, 31):
        below = [q for q in PRIMES if q < target]
        for m in (1, 2, 3, 4):
            if m > len(below):
                continue
            pts = below[-m:]
            c = fit(pts)
            pr = predict(c, target)
            err = pr - R[target]
            worst.setdefault(target, []).append(abs(err))
            print(f'   target q = {target}, order m = {m}, fit on {pts}: predicted {pr:.6f}, exact {R[target]:.6f}, error {err:+.6f}; c = {[round(x, 4) for x in c]}')
    print('\ncontrol q = 27 (prime power, never fitted):')
    for m in (1, 2, 3):
        pts = [19, 23, 31][-m:] if m <= 3 else None
        c = fit(pts)
        print(f'   order {m} through {pts}: predicted R(27) = {predict(c, 27):.6f}, exact {R[27]:.6f}, error {predict(c, 27) - R[27]:+.6f}')
    print('\nfits on all primes up to 31 (largest points), and the second-generation predictions:')
    preds = {43: [], 47: []}
    for m in (1, 2, 3, 4):
        pts = PRIMES[-m:]
        c = fit(pts)
        line = f'   order {m}, fit on {pts}: c = {[round(x, 5) for x in c]}'
        for q in (43, 47):
            v = predict(c, q)
            preds[q].append(v)
            line += f';  R({q}) = {v:.6f}'
        print(line)
    # stable orders: 2 and 3 (order 1 is biased low-order; order 4 uses q = 11 and overreacts) -- declared before seeing
    # any data at 43 (there is none).  Interval = [min, max] over orders 2, 3, widened by the worst holdout error of
    # orders 2, 3 at q = 31 scaled by (31/q).
    h31 = []
    below = [q for q in PRIMES if q < 31]
    for m in (2, 3):
        c = fit(below[-m:])
        h31.append(abs(predict(c, 31) - R[31]))
    pad = max(h31)
    print(f'\n   worst order-2/3 holdout error at q = 31: {pad:.6f}')
    for q in (43, 47):
        lo = min(preds[q][1], preds[q][2]) - pad * 31 / q
        hi = max(preds[q][1], preds[q][2]) + pad * 31 / q
        mid = (preds[q][1] + preds[q][2]) / 2
        print(f'   PRE-REGISTERED: R({q}) in [{lo:.5f}, {hi:.5f}]  (central {mid:.5f});  '
              f'H(P_{q}) in [{lo * math.factorial(q) / 2 ** (q - 1):.6e}, {hi * math.factorial(q) / 2 ** (q - 1):.6e}]')
    print(f'\n   limit check: e = {E:.9f};  c_1 (order 2 on 23, 31) = {fit([23, 31])[0]:.5f}')
    # c_1 = 3/2 fixed (first-order cluster derivation, paley_rate_c1_ingredients_20261014.py); fit c_2 (and c_3)
    print('\nwith c_1 = 3/2 FIXED (derived): 1 - R/e - 3/(2q) = c_2/q^2 + c_3/q^3 + ...')
    def resid(q):
        return 1 - R[q] / E - 1.5 / q
    for q in sorted(R):
        print(f'   q = {q:2d}: q^2 (1 - R/e - 3/(2q)) = {q * q * resid(q):+.5f}' + ('  (prime power)' if q == 27 else ''))
    c2a = resid(31) * 31 ** 2
    A = [[Fraction(1, 23 ** 2), Fraction(1, 23 ** 3)], [Fraction(1, 31 ** 2), Fraction(1, 31 ** 3)]]
    c2b, c3b = solve(A, [Fraction(resid(23)).limit_denominator(10 ** 15), Fraction(resid(31)).limit_denominator(10 ** 15)])
    print(f'   c_2 (from q = 31 alone) = {c2a:.4f};  (c_2, c_3) from q = 23, 31 = ({c2b:.4f}, {c3b:.4f})')
    c2h = resid(23) * 23 ** 2
    p31a = E * (1 - 1.5 / 31 - c2h / 31 ** 2)
    A = [[Fraction(1, 19 ** 2), Fraction(1, 19 ** 3)], [Fraction(1, 23 ** 2), Fraction(1, 23 ** 3)]]
    c2hb, c3hb = solve(A, [Fraction(resid(19)).limit_denominator(10 ** 15), Fraction(resid(23)).limit_denominator(10 ** 15)])
    p31b = E * (1 - 1.5 / 31 - c2hb / 31 ** 2 - c3hb / 31 ** 3)
    print(f'   holdout at 31 with c_1 = 3/2: (c_2 from 23) {p31a:.6f} (error {p31a - R[31]:+.6f});  (c_2, c_3 from 19, 23) {p31b:.6f} (error {p31b - R[31]:+.6f})')
    for q in (43, 47):
        va = E * (1 - 1.5 / q - c2a / q ** 2)
        vb = E * (1 - 1.5 / q - c2b / q ** 2 - c3b / q ** 3)
        pad = max(abs(p31a - R[31]), abs(p31b - R[31])) * (31 / q) ** 2
        lo, hi = min(va, vb) - pad, max(va, vb) + pad
        print(f'   c_1 = 3/2 model: R({q}) = {va:.6f} (c_2 only) / {vb:.6f} (c_2, c_3);  PRE-REGISTERED (c_1-fixed): R({q}) in [{lo:.5f}, {hi:.5f}]')


if __name__ == '__main__':
    main()
