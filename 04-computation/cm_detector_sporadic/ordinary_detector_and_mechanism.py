#!/usr/bin/env python3
"""P6 and P7 of PREREGISTRATION.md (third pre-registration, pushed in 4a60d0195a before this script was written).
kind-pasteur-2026-10-18-S2.

P6: ordinary-restricted Lucas-block detector (B1'_ord): a prime p counts for m only if p does not divide m and
    F_p(1/m) != 0 (mod p); scan 0 < |m| <= M, filter p in 7..43, confirm p in 47..200.
P7: over F_p, for 7 <= p <= 199: rad(F_p) | V and every root of F_p multiple (weight 2); F_p squarefree and
    gcd(F_p, V) = 1 (weight 1).  F_p = sum_{k<p} u_k x^k, V = sum_{j<p} ((u_{p+j} - u_1 u_j)/p) x^j, both mod p.
Usage: python ordinary_detector_and_mechanism.py [M]
"""
import sys, json
import numpy as np
from sympy import Poly, symbols, gcd as pgcd, div
from sporadic_cm_detector import FAMILIES, primes_between

X = symbols('x')
WEIGHT2 = {'Apery zeta3', 'Domb', 'Cooper s7', 'Cooper s10'}


def block_test(u, ms, p, ordinary_only):
    pp = p * p
    A = [u[k] % pp for k in range(2 * p)]
    r = np.mod(ms, pp)
    invtab = np.zeros(pp, dtype=np.int64)
    for x in range(pp):
        if x % p:
            invtab[x] = pow(x, -1, pp)
    w = invtab[r]
    S0 = np.zeros_like(ms); S1 = np.zeros_like(ms); pw = np.ones_like(ms)
    for k in range(2 * p):
        if k < p:
            S0 = (S0 + A[k] * pw) % pp
        else:
            S1 = (S1 + A[k] * pw) % pp
        pw = (pw * w) % pp
    ok = (S1 - (A[1] * w % pp) * S0) % pp == 0
    skip = (r % p == 0)
    if ordinary_only:
        skip = skip | (S0 % p == 0)
    return ok | skip


def scan(u, M, ordinary_only):
    ms = np.array([m for m in range(-M, M + 1) if m != 0], dtype=np.int64)
    for p in primes_between(7, 43):
        ms = ms[block_test(u, ms, p, ordinary_only)]
    for p in primes_between(47, 200):
        ms = ms[block_test(u, ms, p, ordinary_only)]
    return sorted(int(x) for x in ms)


def mechanism(u, p, weight2):
    F = Poly([u[k] % p for k in range(p)][::-1], X, modulus=p)
    V = Poly([((u[p + j] - u[1] * u[j]) // p) % p for j in range(p)][::-1], X, modulus=p)
    dF = F.diff(X)
    g = pgcd(F, dF)
    rad = div(F, g)[0]
    rad_divides_V = div(V, rad)[1].is_zero
    all_multiple = div(g, rad)[1].is_zero if g.degree() > 0 else (F.degree() == 0)
    squarefree = g.degree() == 0
    coprime = pgcd(F, V).degree() == 0
    return {'deg F': F.degree(), 'deg rad': rad.degree(), 'rad|V': rad_divides_V, 'all roots multiple': all_multiple,
            'squarefree': squarefree, 'gcd(F,V)=1': coprime}


def main():
    M = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    out = {}
    for name, gen in FAMILIES.items():
        u = gen()
        w2 = name in WEIGHT2
        surv_ord = scan(u, M, True)
        surv = scan(u, M, False)
        mech = {}
        for p in primes_between(7, 199):
            if any((u[p + j] - u[1] * u[j]) % p for j in range(p)):
                mech[p] = 'Lucas mod p fails'
                continue
            mech[p] = mechanism(u, p, w2)
        if w2:
            p7 = all(isinstance(v, dict) and v['rad|V'] and v['all roots multiple'] for v in mech.values())
        else:
            p7 = all(isinstance(v, dict) and v['squarefree'] and v['gcd(F,V)=1'] for v in mech.values())
        bad = {p: v for p, v in mech.items() if not (isinstance(v, dict) and (
            (w2 and v['rad|V'] and v['all roots multiple']) or ((not w2) and v['squarefree'] and v['gcd(F,V)=1'])))}
        out[name] = {'weight': 2 if w2 else 1, 'survivors_B1prime': surv, 'survivors_B1prime_ord': surv_ord,
                     'P7_holds_all_p': p7, 'P7_exceptions': {str(k): v for k, v in bad.items()}}
        print(f"{name} (weight {2 if w2 else 1}): B1' {surv} | B1'_ord {surv_ord} | P7 {'HOLDS' if p7 else 'FAILS'} "
              f"exceptions at p = {sorted(bad)[:12]}", flush=True)
        if bad:
            for p in sorted(bad)[:4]:
                print('     ', p, bad[p])
    json.dump(out, open(f'ordinary_detector_and_mechanism_M{M}.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
