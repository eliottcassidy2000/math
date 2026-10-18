#!/usr/bin/env python3
"""Li-degree-0 relations (Zhu C_2-ideal) of boundary W_3 minimal models W_3(3,q), with coefficients reconstructed by
CRT from TWO primes (audit 2 / MISTAKE-633: a single-prime rational reconstruction produced a wrong coefficient for
q = 13). A coefficient is printed only if the CRT reconstruction is consistent with both residues; otherwise it is
flagged.  Usage: python w3_c2_relations.py q [q ...]"""
import sys
from fractions import Fraction
from w3 import W3
from lie_voa_ss import SingularSupport, find_singular
from vir_ss import ratrec

PRIMES = (2147483647, 2147483629)

def relations_mod(q_, P):
    c = 2 * (1 - Fraction(12 * (3 - q_) ** 2, 3 * q_))
    alg = W3(c, P)
    sv = find_singular(alg, q_ - 2, ())
    ss = SingularSupport(alg, q_ + 1, sv, verbose=False)
    ss.run(gen_modes=[(1, 0), (2, 0), (1, 1), (2, 1), (3, 1)])
    out = {}
    for D in range(q_ - 2, q_ + 2):
        rows = []
        for v in ss.grI.get((D, 0, ()), []):
            items = sorted(v.items(), key=lambda kv: kv[0], reverse=True)
            lead = items[0][1]; inv = pow(lead, -1, P)
            rows.append({tuple(sorted(t for (w, t) in mo)): co * inv % P for mo, co in items})
        out[D] = rows
    return c, out

def crt(a1, a2, P1, P2):
    M = P1 * P2
    x = (a1 * P2 * pow(P2, -1, P1) + a2 * P1 * pow(P1, -1, P2)) % M
    r = ratrec(x, M)
    ok = r is not None and all(
        (r.numerator - a * r.denominator) % P == 0 for a, P in ((a1, P1), (a2, P2)))
    return r, ok

def show_monomial(key):
    nx = sum(1 for t in key if t == 0); ny = sum(1 for t in key if t == 1)
    return f"x^{nx} y^{ny}"

if __name__ == '__main__':
    for q_ in map(int, sys.argv[1:]):
        res = [relations_mod(q_, P) for P in PRIMES]
        c = res[0][0]
        print(f"W3(3,{q_}) c={c}:")
        for D in sorted(res[0][1]):
            for r1, r2 in zip(res[0][1][D], res[1][1][D]):
                terms = []
                for k in r1:
                    val, ok = crt(r1[k], r2.get(k, 0), *PRIMES)
                    terms.append(f"({val}{'' if ok else ' ??UNCONFIRMED'}) {show_monomial(k)}")
                print(f"   weight {D}: " + " + ".join(terms))
