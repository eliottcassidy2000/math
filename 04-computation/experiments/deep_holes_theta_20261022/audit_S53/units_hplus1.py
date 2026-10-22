"""Audit check (S53): Galois/unit multiples of the (h+1)-points c''.

The principal-part proof of THM-4759 (10) needs Gamma_0(e)-invariance of Z_{c''}, i.e. rationality of
Theta_Lambda^[c'']. A necessary condition (THM-4760 (2)) is theta(Lambda + u c'') = theta(Lambda + c'') for every
unit u mod e. Here we compute, independently of the session's certificate, the norm profile below and at norm 2 of
Lambda + u c'' for every unit u (1 <= u <= e/2; -u gives the same profile), for all 23 types.
Usage: python units_hplus1.py [types...]"""
import sys, time, math
from audit_common import ROWS, h_plus_one_point, coset_profile

ONLY = sys.argv[1:]
bad = []
for name, h, n0, z2 in ROWS:
    if ONLY and name not in ONLY:
        continue
    t0 = time.time()
    Ym, m = h_plus_one_point(z2, h)
    assert m == h + 1
    e = m
    profs = {}
    for u in range(1, e // 2 + 1):
        if math.gcd(u, e) != 1:
            continue
        profs[u] = coset_profile(u * Ym, e, 2.0)
    base = profs[1]
    same = all(p == base for p in profs.values())
    fr = e * e * 2 - 2 * e          # e^2 (2 - 2/e)
    ok1 = min(base) == fr and base[fr] == 1 and base.get(2 * e * e, 0) == 24
    if not (same and ok1):
        bad.append(name)
    desc = {u: {k / (e * e): c for k, c in sorted(p.items())} for u, p in profs.items()}
    print(f'{"ok  " if (same and ok1) else "DIFF"} {name:9s} e={e:2d} units<=e/2: {sorted(profs)}; '
          f'profile(u=1) = {dict(sorted((round(k / (e * e), 6), c) for k, c in base.items()))}; all units equal: {same} '
          f'[{time.time() - t0:.1f}s]', flush=True)
    if not same:
        for u, p in profs.items():
            print(f'       u={u}: {dict(sorted((round(k / (e * e), 6), c) for k, c in p.items()))}')
print('ALL UNIT MULTIPLES MATCH' if not bad else f'MISMATCH: {bad}')
