#!/usr/bin/env python3
"""wheel_theorems_check.py NMAX: numerical checks of the base-2 wheel characterization and the self-dual reductions.
(1) characterization: W(N) covers Z/N iff (a >= 6) or (odd wheel covers Z/m); compared with the exact search for N <= NMAX.
(2) prime-power odd part: odd wheel covers Z/q^e iff sum of q^(-j) over its classes >= 1 (Kraft); checked by exact search.
(3) self-dual: for N = 2m (m odd) a self-dual cover exists iff the flipper wheel {primes of order exactly 2m', m' | m, m' > 1}
    covers Z/m by classes modulo m' (reduction); for 4 | N the parity-sheet bound  sum_keepers 1/d + 2 sum_flippers(d>2) 1/d >= 1
    is necessary.  Both compared with the exact two-sheet search where it finishes."""
import sys
from sympy import factorint, n_order, divisors
NMAX = int(sys.argv[1])
def wheel(N):
    out = []
    for d in divisors(N):
        if d == 1: continue
        for p in factorint(2**d - 1):
            if n_order(2, p) == d: out.append(d)
    return sorted(out)
def cover(N, mods, limit=3_000_000):
    """one class per modulus in mods (each dividing N) covering Z/N? exact; None on timeout"""
    full = (1 << N) - 1
    cm = {}
    for d in set(mods):
        for r in range(d):
            m = 0
            for x in range(r, N, d): m |= 1 << x
            cm[(d, r)] = m
    used = [False] * len(mods); nodes = [0]
    def rec(c):
        nodes[0] += 1
        if c == full: return True
        if nodes[0] > limit: raise TimeoutError
        x = (~c & full); x = (x & -x).bit_length() - 1; tried = set()
        for i, d in enumerate(mods):
            if used[i] or d in tried: continue
            tried.add(d); used[i] = True
            if rec(c | cm[(d, x % d)]): return True
            used[i] = False
        return False
    if sum(1 / d for d in mods) < 1 - 1e-12: return False
    try: return rec(0)
    except TimeoutError: return None
def two_sheet(N, mods, limit=3_000_000):
    full = (1 << N) - 1
    def cls(d, r):
        m = 0
        for x in range(r % d, N, d): m |= 1 << x
        return m
    used = [False] * len(mods); nodes = [0]
    def rec(cm_, cp):
        nodes[0] += 1
        if cm_ == full and cp == full: return True
        if nodes[0] > limit: raise TimeoutError
        if cm_ != full: x = (~cm_ & full); x = (x & -x).bit_length() - 1; sheet = 0
        else: x = (~cp & full); x = (x & -x).bit_length() - 1; sheet = 1
        tried = set()
        for i, d in enumerate(mods):
            if used[i] or d in tried: continue
            tried.add(d); used[i] = True
            r = x if sheet == 0 else x - d // 2
            if rec(cm_ | cls(d, r), cp | cls(d, r + d // 2)): return True
            used[i] = False
        return False
    if sum(1 / d for d in mods) < 1 - 1e-12: return False
    try: return rec(0, 0)
    except TimeoutError: return None
def oddpart(d):
    while d % 2 == 0: d //= 2
    return d
fmt = lambda v: 'YES' if v else ('?' if v is None else 'no')
mism1 = []; mism2 = []; mism3 = []
print("N | exact | characterization (a>=6 or odd wheel) | Kraft(prime power) | self-dual exact | reduction/bound | notes")
for N in range(2, NMAX + 1):
    mods = wheel(N); a = 0; m = N
    while m % 2 == 0: m //= 2; a += 1
    odd = [oddpart(d) for d in mods if oddpart(d) > 1]
    ex = cover(N, mods)
    oddcov = cover(m, odd) if m > 1 else False
    char = True if a >= 6 else oddcov
    # Kraft for prime-power odd part
    kraft = ''
    f = factorint(m) if m > 1 else {}
    if len(f) == 1:
        q = list(f)[0]; kraft = fmt(sum(1 / d for d in odd) >= 1 - 1e-12)
    # self-dual
    even = [d for d in mods if d % 2 == 0]
    sd = two_sheet(N, even) if sum(1 / d for d in even) >= 1 - 1e-12 else False
    if a == 1:
        flip = [d // 2 for d in even if d > 2]          # flipper odd parts
        red = cover(m, flip) if sum(1 / x for x in flip) >= 1 - 1e-12 else False
        note = f"flipper wheel on Z/{m}: {fmt(red)}"
        if sd is not None and red is not None and bool(sd) != bool(red): mism3.append(N)
    else:
        keepers = [d for d in even if (d // 2) % 2 == 0]; flippers = [d for d in even if d > 2 and (d // 2) % 2 == 1]
        bound = sum(1 / d for d in keepers) + 2 * sum(1 / d for d in flippers)
        note = f"parity-sheet bound {bound:.3f}" + (" (<1: impossible)" if bound < 1 - 1e-12 else " (>=1: undecided by bound)")
        if sd and bound < 1 - 1e-12: mism3.append(N)
    if ex is not None and bool(ex) != bool(char): mism1.append(N)
    if kraft and ex is not None and (kraft == 'YES') != bool(ex): mism2.append(N)
    if ex or char or sum(1 / d for d in mods) >= 1:
        print(f"{N} | {fmt(ex)} | {fmt(char)} | {kraft or '-'} | {fmt(sd)} | {note}", flush=True)
print("mismatches characterization:", mism1, "; mismatches Kraft:", mism2, "; self-dual mismatches:", mism3)
