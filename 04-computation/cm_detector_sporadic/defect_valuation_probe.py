#!/usr/bin/env python3
"""EXPLORATORY (post hoc).  kind-pasteur-2026-10-18-S2.

p-adic valuation of the Lucas-block defect  Delta_p(m) = sum_{k=p}^{2p-1} u_k m^{-k} - u_1 m^{-1} sum_{k<p} u_k m^{-k}
at rational CM points (and rational cusps / Fricke-elliptic points), split by the type of p: 'ord' (F_p(1/m) != 0 mod p)
vs 'ss' (F_p(1/m) == 0 mod p).  Computed exactly mod p^5.  S43 (HYP-9311(O)) found, in weight 2 (hypergeometric):
v >= 3 at split primes, v = 2 at inert primes (up to sporadic coincidences).  Question: weight 1?
"""
import sys, json
from collections import Counter
from sporadic_cm_detector import FAMILIES, primes_between
from level5_leeyang_family import level5_U

POINTS = {
    'Zagier A (Franel)': [-4, 2, -1, 8], 'Zagier B': [3, 9], 'Zagier C': [-3, 3, 1, 9], 'Zagier E': [4, 8],
    'Zagier F': [6, 12, 8, 9],
    'Apery zeta3': [1, -1], 'Domb': [-32, -8, -2, 1, 4, 8, 16, 64], 'Cooper s7': [-10648, -64, -8, -1, 2, 8, 24, 27, 125],
    'Cooper s10': [-324, -64, -20, -9, -4, -2, 1, 4, 12, 16, 36, 196, 1296, 5776],
    'level-5 A274786': [-15228, -828, -28, -3, 4, 18, 36, 72, 147],
}
WEIGHT = {k: (1 if k.startswith('Zagier') else 2) for k in POINTS}
CONTROL = 5   # m = 3 is a CM point of Zagier B and C, so 5 is used as the generic control


def val(u, m, p, K=5):
    mod = p ** K
    w = pow(m % mod, -1, mod)
    S0 = 0; S1 = 0; pw = 1
    for k in range(2 * p):
        if k < p:
            S0 = (S0 + u[k] * pw) % mod
        else:
            S1 = (S1 + u[k] * pw) % mod
        pw = pw * w % mod
    d = (S1 - u[1] * w * S0) % mod
    v = 0
    while v < K and d % p == 0:
        d //= p; v += 1
    typ = 'ss' if S0 % p == 0 else 'ord'
    return typ, v


def main():
    seqs = {name: gen() for name, gen in FAMILIES.items()}
    seqs['level-5 A274786'] = level5_U()
    out = {}
    for name, pts in POINTS.items():
        u = seqs[name]
        for m in pts + [CONTROL]:
            c = {'ord': Counter(), 'ss': Counter()}
            for p in primes_between(7, 199):
                if m % p == 0 or u[1] % p == 0:
                    continue
                typ, v = val(u, m, p)
                c[typ][v] += 1
            out[f'{name} m={m}'] = {t: dict(sorted(cc.items())) for t, cc in c.items()}
            tag = 'control' if m == CONTROL else ''
            print(f"w{WEIGHT[name]} {name:18s} m={m:7d} {tag:7s} ordinary v-histogram {dict(sorted(c['ord'].items()))}   "
                  f"supersingular v-histogram {dict(sorted(c['ss'].items()))}", flush=True)
    json.dump(out, open('defect_valuation_probe.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
