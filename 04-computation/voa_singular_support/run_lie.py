#!/usr/bin/env python3
"""Driver for lie_voa_ss: python run_lie.py MODEL params Dmax [prime]
MODEL in {vir p q, sl2 kp kq (k=kp/kq), ns p q}."""
import sys, time
from fractions import Fraction
from lie_voa_ss import *

def main():
    model = sys.argv[1]
    P = int(sys.argv[-1]) if len(sys.argv) > 5 else 2147483647
    if model == 'vir':
        p, q, Dmax = map(int, sys.argv[2:5])
        c = 1 - Fraction(6 * (p - q) ** 2, p * q)
        alg = Virasoro(c, P)
        D0 = (p - 1) * (q - 1)
        sv = find_singular(alg, D0, ())
        print(f"Vir({p},{q}) c={c}: singular space dim {len(sv)} at weight {D0}")
        ss = SingularSupport(alg, Dmax, sv)
        ss.run(gen_modes=[(1, 0), (2, 0)])
    elif model == 'sl2':
        kp, kq, Dmax = map(int, sys.argv[2:5])
        k = Fraction(kp, kq)
        alg = AffineSL2(k, P)
        # admissible k = -2 + a/b, gcd(a,b)=1, a>=2 : singular vector weight (a-1)b, charge 2(a-1)
        t = k + 2
        a, b = t.numerator, t.denominator
        D0, ch0 = (a - 1) * b, (2 * (a - 1),)
        sv = find_singular(alg, D0, ch0)
        print(f"L_k(sl2) k={k} (a/b={a}/{b}): singular space dim {len(sv)} at weight {D0} charge {ch0}")
        ss = SingularSupport(alg, Dmax, sv)
        ss.run(gen_modes=[(1, 0), (1, 1), (1, 2)], extra_ops=[(0, 2), (0, 0)])
    elif model == 'ns':
        p, q, Dmax = map(int, sys.argv[2:5])
        c = Fraction(3, 2) * (1 - Fraction(2 * (p - q) ** 2, p * q))
        alg = NeveuSchwarz(c, P)
        D0 = Fraction((p - 1) * (q - 1), 2)
        sv = find_singular(alg, D0, ())
        print(f"NS({p},{q}) c={c}: singular space dim {len(sv)} at weight {D0}")
        ss = SingularSupport(alg, Dmax, sv)
        ss.run(gen_modes=[(Fraction(1, 2), 1), (Fraction(3, 2), 1)])
    print("generators (D, Li, charge, #):", ss.gens_found)

if __name__ == '__main__':
    main()
