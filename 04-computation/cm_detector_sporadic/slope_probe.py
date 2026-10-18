#!/usr/bin/env python3
"""EXPLORATORY (post hoc, not pre-registered).  kind-pasteur-2026-10-18-S2.

Rewriting of the block test (B1').  With x = 1/m, F_p(x) = sum_{k<p} u_k x^k, Lucas defect
v_j = (u_{p+j} - u_1 u_j)/p and V(x) = sum_{j<p} v_j x^j, and using x^p == x (1 - p q_p(m)) (mod p^2) with
q_p(m) = (m^{p-1}-1)/p the Fermat quotient:

    (B1')  <=>  V(x) == q_p(m) u_1 F_p(x)   (mod p).

Define the SLOPE sigma_p(m) = V(x) / (u_1 q_p(m) F_p(x)) in F_p (ordinary p, q_p(m) != 0). The weight-2 detector says
sigma_p = 1 at CM points.  Question: at weight-1 CM points, is sigma_p a fixed rational (e.g. 1/2) for all p?
At supersingular p (F_p(x) == 0) the test reduces to V(x) == 0; we record that separately.
"""
import json, sys
from fractions import Fraction
from sporadic_cm_detector import FAMILIES, primes_between
import numpy as np

PR = primes_between(7, 257)


def lucas_ok(u, p):
    """Lucas mod p at the first carry: u_{p+j} == u_1 u_j (mod p) for all j < p (needed for v_j integrality)."""
    return all((u[p + j] - u[1] * u[j]) % p == 0 for j in range(p))


def slope_data(u, m):
    rows = []
    for p in PR:
        if m % p == 0 or u[1] % p == 0:
            continue
        x = pow(m % p, -1, p)
        Fp = 0; V = 0; w = 1
        for j in range(p):
            Fp = (Fp + u[j] * w) % p
            V = (V + ((u[p + j] - u[1] * u[j]) // p) * w) % p
            w = w * x % p
        qp = ((pow(m, p - 1, p * p) - 1) // p) % p
        if Fp == 0:
            rows.append((p, 'ss', V))
        elif qp == 0:
            rows.append((p, 'wieferich', None))
        else:
            sig = V * pow(u[1] * qp * Fp % p, -1, p) % p
            rows.append((p, 'ord', sig))
    return rows


def as_small_rational(res_by_p):
    """Find a/b with |a|,|b| <= 60 matching sigma_p at every ordinary p (CRT-free: test candidates)."""
    cands = []
    for b in range(1, 61):
        for a in range(-60, 61):
            if np.gcd(a, b) != 1:
                continue
            if all((a - s * b) % p == 0 for p, s in res_by_p):
                cands.append(Fraction(a, b))
    return cands


def main():
    targets = json.loads(sys.argv[1]) if len(sys.argv) > 1 else None
    seqs = {name: gen() for name, gen in FAMILIES.items()}
    for name, u in seqs.items():
        bad = [p for p in PR if not lucas_ok(u, p)]
        if bad:
            print(f"{name}: Lucas mod p FAILS at {bad[:5]} -> v_j not integral; skipped")
            continue
        ms = targets.get(name, []) if targets else []
        for m in ms:
            rows = slope_data(u, m)
            ordp = [(p, s) for p, t, s in rows if t == 'ord']
            ss = [(p, s) for p, t, s in rows if t == 'ss']
            fixed = as_small_rational(ordp)
            ss_zero = sum(1 for p, s in ss if s == 0)
            print(f"{name:24s} m={m:7d}: ordinary {len(ordp):3d} primes, slope fixed rational = {[str(f) for f in fixed]}; "
                  f"slopes p<=43: {[s for p, s in ordp if p <= 43]}; supersingular {len(ss)}, V==0 at {ss_zero}", flush=True)


if __name__ == '__main__':
    main()
