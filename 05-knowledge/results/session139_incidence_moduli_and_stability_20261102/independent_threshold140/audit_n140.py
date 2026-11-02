#!/usr/bin/env python3
"""Independent exact controls for the tensor-reserve n >= 140 proof.

Explicit checks remain active under -O. Negative controls: --threshold 139
and --mutate-perfect-indicator must both fail. Uniform proof is in AUDIT.md.
"""
import argparse
from fractions import Fraction as F
import json
from math import factorial
import sys


class CheckFailure(Exception):
    pass


def run(n, mutate_perfect_indicator=False):
    evidence = {}
    bootstrap = []

    def check(name, condition, value=None):
        if not condition:
            raise CheckFailure(f"{name}: {value}")
        evidence[name] = True if value is None else str(value)

    R = n // 2
    eps = F(1, 10000)
    check('order_domain', n >= 8 and R > 20)
    check('small_variance_excluded', F(49, 200) > F(4, R), F(49, 200)-F(4, R))
    check('exp4_partial_lower', sum((F(4)**j/factorial(j) for j in range(8)), F(0)) > 50)
    check('variance_above4', F(49*R, 800) > 4, F(49*R, 800)-4)
    c0 = F(249999, 1000000)
    check('initial_nonlinear_margin', (1-F(1, 50**4))/4 > c0,
          (1-F(1, 50**4))/4-c0)
    check('perfect_penalty_envelope', F(1, 4**R) < eps)
    # Test the correct perfect-matching indicator at representative cardinalities.
    # The universal product estimate itself follows from every w <= 1/4.
    for order in (2*R, 2*R+1):
        for r in (1, R-1, R, 0):
            product_max = F(1, 4**r)
            correction = product_max if (mutate_perfect_indicator or order == 2*r) else F(0)
            check(f'perfect_penalty_case_n{order}_r{r}', correction <= F(1, 4**R))
    constants = [c0, F(7,25), F(57,100), F(99,100)]
    for i, (c, target, degree) in enumerate(zip(constants, constants[1:], (2,3,4))):
        m = c*R/4 - 4 - F(2,R)
        check(f'bootstrap_mass_positive_{i}', m > 0, m)
        check(f'bootstrap_mass_uniform_slope_{i}', c/4 > 0)
        T = sum((m**j/factorial(j) for j in range(degree+1)), F(0))
        remaining = 1-target-eps
        check(f'bootstrap_remaining_positive_{i}', remaining > 0, remaining)
        coarse_slack = remaining*T-1
        exact_slack = (1-target-F(1,4**R))*T-1
        check(f'bootstrap_step_{i}', coarse_slack > 0, coarse_slack)
        check(f'bootstrap_exact_stronger_{i}', exact_slack > coarse_slack)
        bootstrap.append({'from':str(c), 'to':str(target), 'mass':str(m),
                          'degree':degree, 'coarse_slack':str(coarse_slack)})
    c, d = constants[-1], F(19,20)
    t, lam = d/8, d/64
    h = d/16
    heavy = c*(n-1)/2 - F(16*n,n-1)-d*n/4
    check('heavy_count_exceeds19', heavy > 19, heavy-19)
    check('heavy_count_uniform_slope', c/2-d/4 > 0, c/2-d/4)
    Q = F(R)+F(33,2)+F(8,R)
    check('variance_even_envelope', R+16 < Q)
    check('variance_odd_envelope', F(2*R+1,2)+2*(2*R+1)*F(4,R) == Q)
    check('tail_interval_nonempty', R-20 >= 1)
    rho = F(3,4)
    check('left_tail_endpoint', Q*rho**(R-1) < c-d, c-d-Q*rho**(R-1))
    check('right_tail_endpoint', Q*rho**20/(R-20) < c-d, c-d-Q*rho**20/(R-20))
    check('left_tail_uniform_ratio', rho*(1+F(1,R)) < 1, 1-rho*(1+F(1,R)))
    check('right_tail_decreasing_form', Q/(R-20) == 1+(F(73,2)+F(8,R))/(R-20))
    check('light_deficit_positive', d-4*t > 0, d-4*t)
    left={(0,2):F(1,4),(1,1):F(1,2),(0,1):-t,(0,0):lam*d,(1,0):-4*lam}
    right={(0,2):F(1,4),(1,1):F(1,2),(0,1):-t,(0,0):t*t,(1,0):-t/2}
    check('dual_polynomial_identity', left == right)
    check('dual_remainder_domain', 0 < t <= 2*t)
    check('h_identity', h == t-4*lam)
    check('final_margin', h > F(4,R), h-F(4,R))
    check('no_light_closure', t > F(4,R), t-F(4,R))
    for order in (n,n+1):
        r=order//2
        kappa=F(4,r)*(1-F(3,4)**r-(F(1,4)**r if order%2==0 else 0))
        check(f'parity_coefficient_{order}', 0 < kappa <= F(4,r))
    return {'status':'PASS', 'threshold':n, 'R_min':R,
            'coarse_penalty_bound':str(eps), 'bootstrap':bootstrap,
            'active_checks':evidence}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--threshold',type=int,default=140)
    p.add_argument('--mutate-perfect-indicator',action='store_true')
    args=p.parse_args()
    try:
        result=run(args.threshold,args.mutate_perfect_indicator)
    except CheckFailure as exc:
        print(json.dumps({'status':'REJECTED','reason':str(exc)},sort_keys=True))
        return 1
    print(json.dumps(result,sort_keys=True,indent=2))
    return 0


if __name__ == '__main__':
    sys.exit(main())
