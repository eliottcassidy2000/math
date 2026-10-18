import sys
from fractions import Fraction
from w3 import W3
from lie_voa_ss import SingularSupport, find_singular
from vir_ss import ratrec
P = 2147483647
def show(vec):
    items = sorted(vec.items(), key=lambda kv: kv[0], reverse=True)
    lead = items[0][1]; inv = pow(lead, -1, P)
    out = []
    for mo, co in items:
        nx = sum(1 for (w, t) in mo if t == 0); ny = sum(1 for (w, t) in mo if t == 1)
        out.append(f"({ratrec(co*inv % P, P)}) x^{nx} y^{ny}")
    return " + ".join(out)
for q_ in map(int, sys.argv[1:]):
    c = 2 * (1 - Fraction(12 * (3 - q_) ** 2, 3 * q_))
    alg = W3(c, P)
    sv = find_singular(alg, q_ - 2, ())
    ss = SingularSupport(alg, q_ + 1, sv, verbose=False)
    ss.run(gen_modes=[(1, 0), (2, 0), (1, 1), (2, 1), (3, 1)])
    print(f"W3(3,{q_}) c={c}:")
    for D in range(q_ - 2, q_ + 2):
        g = ss.grI.get((D, 0, ()), [])
        for v in g:
            print(f"   weight {D}: {show(v)}")
