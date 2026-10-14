#!/usr/bin/env python3
"""AUDIT of HYP-9274 (D): which run configurations carry the wheel W_4 at the leading arithmetic order?

For every multiset of even run lengths with total 8 edges -- (2,6), (4,4), (2,2,4), (2,2,2,2) -- and the single run
(8) (control: polynomial), compute the injective joint Sigma*(union of paths) = P(q) + C(q) Sigma(W_4) with the
referee's symbolic engine (audit_A10_cores.py).  Each r-run configuration contributes to R(q)/e at relative order
  (placements ~ q^r / prod(mult!)) * C * Sigma(W_4) / q^(8 + r),   Sigma(W_4) = -q(q-1)(a_q + q),
i.e. -(C / prod(mult!)) a_q q^(-6): all r = 2, 3, 4 enter at the SAME order a_q/q^6.
Also: single runs of length 10 (A_10) and 12 (A_12, if feasible) for the single-run channel.
"""
import sys, time
from math import factorial
from collections import Counter
sys.path.insert(0, '.')
from audit_A10_cores import injective_sum_paths
from audit_run_polys_symbolic import pstr
from fractions import Fraction


def main():
    total = Fraction(0)
    for lens in ([8], [2, 6], [4, 4], [2, 2, 4], [2, 2, 2, 2]):
        t0 = time.time()
        res, cnt = injective_sum_paths(lens)
        cores = {k: v for k, v in res.items() if k != 'P'}
        mult = 1
        for v in Counter(lens).values():
            mult *= factorial(v)
        line = f'runs {lens}: {cnt} partitions ({time.time() - t0:.1f} s); polynomial part degree {max(res.get("P", {0: 0}).keys())}; '
        if not cores:
            line += 'no non-SP core'
        for k, v in cores.items():
            line += f'core {k}: C(q) = {pstr(v)}; '
            if len(lens) > 1 and set(v.keys()) == {0}:
                total += Fraction(v[0], mult)
        print(line, flush=True)
    print(f'\nleading W_4 coefficient in R/e:  -(sum over configurations of C/prod(mult!)) * a_q / q^6 = {-total} * a_q / q^6'
          f'   (Sigma(W_4) ~ -q^2 a_q in the canonical orientation)')


if __name__ == '__main__':
    main()
