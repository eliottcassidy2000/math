#!/usr/bin/env python3
"""zsig_wheels.py: Zsigmondy wheels.  For a period N the wheel W(N) is the set of primes p dividing 2^N - 1, i.e. the
primes with ord_p(2) | N (one residue class modulo ord_p(2) each).  A covering of Z/N by one class per prime is a
covering set of period N for numbers k*2^n +- 1 (Sierpinski/Riesel type).  For each N <= NMAX: the wheel, the orders,
the density sum 1/ord, whether a one-class-per-prime cover of Z/N exists (exact bitmask search with pruning), the number
of covering patterns up to translation (when the search is small enough), and whether a cover exists using only the
PRIMITIVE prime divisors (one per d | N, choosing one prime when several).  usage: zsig_wheels.py NMAX"""
import sys
from sympy import factorint, n_order
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 60
def covers_exist(N, mods, count_limit=2_000_000):
    """exact search: for primes with moduli mods (each modulus divides N), choose residues r_i mod d_i so the classes
    cover Z/N.  Returns (exists, number_of_covers_or_None).  Classes as bitmasks over Z/N."""
    full = (1 << N) - 1
    classes = []
    for d in mods:
        cl = []
        for r in range(d):
            m = 0
            for x in range(r, N, d): m |= 1 << x
            cl.append(m)
        classes.append(cl)
    # order primes by decreasing modulus coverage? put small moduli (big classes) first
    order = sorted(range(len(mods)), key=lambda i: mods[i])
    classes = [classes[i] for i in order]; mods2 = [mods[i] for i in order]
    # remaining coverage bound: sum of class sizes of remaining primes
    rem = [0] * (len(mods2) + 1)
    for i in range(len(mods2) - 1, -1, -1): rem[i] = rem[i + 1] + N // mods2[i]
    count = [0]; found = [False]; nodes = [0]
    def dfs(i, covered):
        nodes[0] += 1
        if covered == full: count[0] += 1; found[0] = True; return
        if i == len(mods2): return
        missing = N - bin(covered).count('1')
        if missing > rem[i]: return
        if count[0] >= count_limit or nodes[0] > 50_000_000: return
        for m in classes[i]:
            dfs(i + 1, covered | m)
    dfs(0, 0)
    exact = nodes[0] <= 50_000_000 and count[0] < count_limit
    return found[0], (count[0] if exact else None)
print("N | wheel primes (order) | sum 1/ord | one-class-per-prime cover of Z/N? | #covers (translations counted) | primitive-only cover?")
for N in range(1, NMAX + 1):
    f = factorint(2**N - 1)
    primes = sorted(f)
    ords = {p: n_order(2, p) for p in primes}
    dens = sum(1 / ords[p] for p in primes)
    prim = {}
    for p in primes:
        prim.setdefault(ords[p], []).append(p)
    if dens < 1 - 1e-12:
        print(f"{N} | {[(p, ords[p]) for p in primes]} | {dens:.3f} | no (density < 1) | 0 | no")
        continue
    ex, cnt = covers_exist(N, [ords[p] for p in primes])
    # primitive-only: one prime per order present
    mods_prim = sorted(prim)
    exp, cntp = covers_exist(N, mods_prim) if ex else (False, 0)
    print(f"{N} | {[(p, ords[p]) for p in primes]} | {dens:.3f} | {'YES' if ex else 'no'} | {cnt if cnt is not None else '>limit'} | {'YES' if exp else 'no'}{' (primes per order ' + str({d: len(v) for d, v in prim.items() if len(v) > 1}) + ')' if any(len(v) > 1 for v in prim.values()) else ''}", flush=True)
