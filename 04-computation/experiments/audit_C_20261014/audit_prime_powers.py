#!/usr/bin/env python3
"""Audit C: THM-4653 and HYP-9276 over prime-power fields with this audit's own GF(p^r) (audit_common.GF:
log/antilog tables from a generator of the unit group; defining polynomial chosen differently from the audited code).

Part 1 (THM-4653 internals at prime powers; the audited proof-step script covered primes 5..101 only):
  for q in {9, 25, 27, 49, 81, 121, 125} check
    t(0) = -eps, t(1) = -1;  unconstrained sum = -eps M;  TermA - main = 3 + 3 eps;  TermB = TermC = 2 + 5 eps - eps q;
    W = -eps M + 7 + 13 eps - 2 eps q;  h(x) = sum_y phi(y(1-y)) t(xy) (x != 0);
    t = eps q 2F1(phi,phi;e|x), h = q^2 3F2(phi,phi,phi;e,e|x) with Greene's definitions as exact rationals (x != 0);
    and (q <= 125) W from an affine-normalised brute force over injective tuples (audit_V_bruteforce.brute).
Part 2 (HYP-9276 Frobenius form): M_q (exact, O(q^2)) against -(alpha^(2r) + (alpha beta)^r + beta^(2r)) - 2q - 1, with
  alpha + beta = b(p) from the eta/theta construction, alpha beta = chi_-4(p) p^2, for 22 prime powers up to 2401.
Usage: python3 -u audit_prime_powers.py > audit_prime_powers.out
"""
import sys
import json
import time
from fractions import Fraction as Fr
import numpy as np
sys.path.insert(0, '/Users/e/Documents/GitHub/math-wt-wilson-20261014/04-computation/experiments/audit_C_20261014')
from audit_common import GF, M_field
from audit_V_bruteforce import brute, field_tables

D = '/Users/e/Documents/GitHub/math-wt-wilson-20261014/04-computation/experiments/audit_C_20261014/'
bp = {int(k): tuple(v) for k, v in json.load(open(D + 'audit_newform_bp.json')).items()}


def internals(F, do_brute):
    q, c, eps = F.q, F.chi, F.eps
    R = np.arange(q)
    sub = np.empty((q, q), dtype=np.int64)
    for a in range(q):
        sub[:, a] = F.sub_vec(R, a)                    # sub[x, a] = x - a
    one = 1
    g = c * c[sub[:, one]]                              # chi(x) chi(x-1)
    base = g[None, :] * c[sub.T]                        # base[a, x] = chi(x)chi(x-1)chi(x-a) = f4(x; a)
    t = base.sum(axis=1)
    F5 = base @ c[sub]                                  # F5[a, b] = sum_x f4(x;a) chi(x - b)
    K = c[:, None] * base
    M = -int(t @ K @ t)
    A = R[:, None]
    B = R[None, :]
    w = c[A] * c[B] * c[sub[B, one]]
    main = int((w * F5 * t[:, None]).sum())
    mask = (A != 0) & (A != one) & (B != 0) & (B != one) & (A != B)
    TA = int((mask * w * F5 * t[:, None]).sum())
    TB = -int((mask * w * F5 * base).sum())             # base[a, b] = f4(b; a)
    negb = np.array([F.neg(int(b)) for b in R])
    one_minus_b = np.array([F.sub_scalar(one, int(b)) for b in R])
    # TermC: chi(-b) + chi(1-b) + chi(a-b);  sub[a, b] = a - b
    TC = int((mask * w * (c[negb][None, :] + c[one_minus_b][None, :] + c[sub[R[:, None], R[None, :]]])).sum())
    W = TA + TB + TC
    ok = {
        't(0) = -eps, t(1) = -1': (t[0] == -eps and t[one] == -1),
        'main = -eps M': main == -eps * M,
        'TermA - main = 3 + 3eps': TA - main == 3 + 3 * eps,
        'TermB = 2 + 5eps - eps q': TB == 2 + 5 * eps - eps * q,
        'TermC = 2 + 5eps - eps q': TC == 2 + 5 * eps - eps * q,
        'W = -eps M + 7 + 13 eps - 2 eps q': W == -eps * M + 7 + 13 * eps - 2 * eps * q,
    }
    # h by convolution and by the multiplicative form
    h = np.array([int(np.dot(c * t, c[sub[x, :]])) for x in R])      # h(x) = sum_a chi(a) t(a) chi(x - a)
    mul = lambda x, y: F.mul(int(x), int(y))
    hk = [sum(int(c[y]) * int(c[one_minus_b[y]]) * int(t[mul(x, y)]) for y in R) for x in R]
    ok['h(x) = sum_y phi(y(1-y)) t(xy), x != 0'] = all(h[x] == hk[x] for x in range(1, q))
    ok['M = -sum chi(x(x-1)) t h'] = (M == -int(np.sum(g * t * h)))
    # Greene 2F1, 3F2 as exact rationals
    phim1 = eps
    F21 = [Fr(0)] * q
    for x in range(1, q):
        s = 0
        for y in range(q):
            s += int(c[y]) * int(c[one_minus_b[y]]) * int(c[F.sub_scalar(one, mul(x, y))])
        F21[x] = Fr(phim1, q) * s
    F32 = [Fr(0)] * q
    for x in range(1, q):
        F32[x] = Fr(phim1, q) * sum(F21[mul(x, y)] * int(c[y]) * int(c[one_minus_b[y]]) for y in range(q))
    ok['t = eps q 2F1 (x != 0)'] = all(t[x] == eps * q * F21[x] for x in range(1, q))
    ok['h = q^2 3F2 (x != 0)'] = all(h[x] == q * q * F32[x] for x in range(1, q))
    ok['M = -q^3 sum phi(x)phi(1-x) 2F1 3F2'] = (M == -q ** 3 * sum(int(c[x]) * int(c[one_minus_b[x]]) * F21[x] * F32[x]
                                                                  for x in range(q)))
    if do_brute:
        subT, chi = field_tables(F)
        Vb = brute(q, subT, chi, 'affine')
        ok['brute V = q(q-1) W'] = (Vb == q * (q - 1) * W)
    return M, W, ok


def predicted_M(p, r):
    a, b = bp[p]
    s1 = (a, b)
    nb = (1 if p % 4 == 1 else -1) * p * p
    s = [(2, 0), s1]
    for k in range(2, 2 * r + 1):
        u = (s1[0] * s[k - 1][0] - s1[1] * s[k - 1][1], s1[0] * s[k - 1][1] + s1[1] * s[k - 1][0])
        s.append((u[0] - nb * s[k - 2][0], u[1] - nb * s[k - 2][1]))
    assert s[2 * r][1] == 0
    S = s[2 * r][0] + nb ** r
    return -S - 2 * p ** r - 1, S


def main():
    t0 = time.time()
    print("Part 1: THM-4653 internals over prime-power fields")
    for p, r in ((3, 2), (5, 2), (3, 3), (7, 2), (3, 4), (11, 2), (5, 3)):
        F = GF(p, r)
        M, W, ok = internals(F, do_brute=(p ** r <= 125))
        bad = [k for k, v in ok.items() if not v]
        print(f"  q = {p}^{r} = {F.q:4d} (defining poly tail {F.poly_tail}), eps = {F.eps:+d}: M = {M}, W = {W}; "
              f"{len(ok)} identities, failures {bad}  [{time.time() - t0:.1f}s]", flush=True)
    print("Part 2: HYP-9276 Frobenius form at prime powers (own field arithmetic, own b(p))")
    cases = [(3, 2), (3, 3), (5, 2), (7, 2), (3, 4), (11, 2), (5, 3), (13, 2), (3, 5), (17, 2), (7, 3), (19, 2),
             (23, 2), (5, 4), (3, 6), (29, 2), (31, 2), (11, 3), (37, 2), (3, 7), (13, 3), (7, 4)]
    nok = 0
    for p, r in cases:
        F = GF(p, r)
        M, _ = M_field(F)
        Mp, S = predicted_M(p, r)
        W = -F.eps * M + 7 + 13 * F.eps - 2 * F.eps * F.q
        Wp = F.eps * (S + 14) + 7
        good = (M == Mp) and (W == Wp)
        nok += good
        print(f"  q = {p}^{r} = {F.q:5d}  eps = {F.eps:+d}  M_q = {M}  predicted {Mp}  S_q = {S}  W_q = {W}  "
              f"{'ok' if good else 'FAIL'}  [{time.time() - t0:.1f}s]", flush=True)
    print(f"HYP-9276 prime-power form: {nok}/{len(cases)} prime powers agree")


if __name__ == '__main__':
    main()
