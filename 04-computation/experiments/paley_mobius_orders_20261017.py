#!/usr/bin/env python3
"""Moebius orders of the Paley tournament: rank the vertices of P_p by the least residue of phi(x), phi in PGL_2(F_p)
(the pole of phi first). Which excesses occur? (S40, THM-4666 statement 1b)

Checks first the marked-coefficient lemma used in the proof of the inverse-order identity,
    N(t) = sum_(a=1)^(p-1) chi(a) [<a t> > a] = h(-p) (1 + chi(t - 1) - chi(t))   (t != 0, 1),
then, for every phi = (a x + b)/(c x + d) up to scalars (or a random sample for larger p), computes the excess of the order
"pole first, then x sorted by <phi(x)>" and records the set of values of excess / h(-p).
Usage: python paley_mobius_orders_20261017.py PMAX [SAMPLE]
Session opus-2026-10-17-S40.
"""
import sys, os, json
import numpy as np
from numba import njit
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from paley_class_number_ranking_20261017 import kronecker_table, class_number, is_prime


@njit(cache=True)
def excess_of(order, chi, p):
    s = 0
    n = order.shape[0]
    for i in range(n):
        u = order[i]
        for j in range(i + 1, n):
            s += chi[(order[j] - u) % p]
    return s


def mobius_order(a, b, c, d, p):
    if c % p == 0:
        key = {x: (a * x + b) * pow(d, -1, p) % p for x in range(p)}
        return sorted(range(p), key=lambda x: key[x])
    pole = (-d * pow(c, -1, p)) % p
    rest = [x for x in range(p) if x != pole]
    key = {x: (a * x + b) * pow((c * x + d) % p, -1, p) % p for x in rest}
    return [pole] + sorted(rest, key=lambda x: key[x])


def main():
    PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 31
    SAMPLE = int(sys.argv[2]) if len(sys.argv) > 2 else 400
    rng = np.random.default_rng(20261017)
    out = {}
    for p in [q for q in range(7, PMAX + 1, 4) if is_prime(q)]:
        chi = kronecker_table(-p)
        chil = [int(x) for x in chi]
        h = class_number(-p)
        # marked-coefficient lemma
        lemma = all(sum(chil[a] for a in range(1, p) if (a * t) % p > a) == h * (1 + chil[(t - 1) % p] - chil[t])
                    for t in range(2, p))
        vals = {}
        mats = []
        if p <= 23:
            for a in range(p):
                for b in range(p):
                    for c in range(p):
                        for d in range(p):
                            if (a * d - b * c) % p:
                                mats.append((a, b, c, d))
        else:
            while len(mats) < SAMPLE:
                a, b, c, d = [int(x) for x in rng.integers(0, p, 4)]
                if (a * d - b * c) % p:
                    mats.append((a, b, c, d))
        for (a, b, c, d) in mats:
            e = int(excess_of(np.array(mobius_order(a, b, c, d, p), dtype=np.int64), chi, p))
            key = 'affine' if c % p == 0 else 'proper'
            vals.setdefault(key, {})
            r = e / h
            vals[key][r] = vals[key].get(r, 0) + 1
        out[p] = dict(h=h, lemma=lemma, values={k: {str(r): n for r, n in sorted(v.items())} for k, v in vals.items()})
        print('p = %d (h = %d, p mod 8 = %d): marked-coefficient lemma %s; excess/h for affine phi: %s; for non-affine phi: %s'
              % (p, h, p % 8, lemma, sorted(vals.get('affine', {}).keys()),
                 sorted(vals.get('proper', {}).keys())[:12] + (['...'] if len(vals.get('proper', {})) > 12 else [])), flush=True)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'paley_mobius_orders_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
