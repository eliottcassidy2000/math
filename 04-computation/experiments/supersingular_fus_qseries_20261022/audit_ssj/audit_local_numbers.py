#!/usr/bin/env python3
"""Independent enumeration of the local numbers s_X(O) of section 6.3 (audit_ssj), over a COMPLETE
residue system: s(floor) depends only on p mod N, s(surface) only on c = (1+p)/4 mod N, i.e. on p mod 4N,
so enumerating every admissible class proves the table entries for all primes p (p >= 5, p not dividing N).
Own arithmetic in O/NO with O = Z[pi] (pi^2 = -p) or Z[w] (w^2 = w - c, pi = 2w - 1); no matrices from the note.
  X_0(N): number of pi-stable cyclic subgroups of order N of O/NO;
  X_1(N): number of +-classes of elements of exact order N with pi x = +-x.
Run: python3 -B audit_local_numbers.py > audit_local_numbers.out
"""
from math import gcd

bad = 0


def fail(m):
    global bad
    bad += 1
    print("FAIL:", m)


def ring(N, kind, r):
    """elements (a, b) of O/NO and multiplication by pi, for p == r (mod N) [floor] or c == r (mod N) [surface]"""
    if kind == "floor":
        mp = lambda a, b: ((-r * b) % N, a % N)               # pi (a + b pi) = -p b + a pi
    else:
        c = r
        mw = lambda a, b: ((-c * b) % N, (a + b) % N)         # w (a + b w) = a w + b (w - c)
        mp = lambda a, b: ((2 * mw(a, b)[0] - a) % N, (2 * mw(a, b)[1] - b) % N)   # pi = 2w - 1
    return mp


def order(a, b, N):
    return N // gcd(gcd(a, b), N)


def s_local(N, kind, r, typ):
    mp = ring(N, kind, r)
    elts = [(a, b) for a in range(N) for b in range(N) if order(a, b, N) == N]
    if typ == 1:   # X_1
        cnt = sum(1 for (a, b) in elts if mp(a, b) in ((a, b), ((-a) % N, (-b) % N)))
        return cnt // 2
    # X_0: stable cyclic subgroups <x>: pi x in <x>
    subs = set()
    for (a, b) in elts:
        S = frozenset(((k * a) % N, (k * b) % N) for k in range(N))
        x = mp(a, b)
        if x in S:
            subs.add(S)
    return len(subs)


def legendre(a, l):
    a %= l
    if a == 0:
        return 0
    return 1 if pow(a, (l - 1) // 2, l) == 1 else -1


results = {}
for N, typ, name in ((2, 0, "X_0(2)"), (4, 0, "X_0(4)"), (8, 0, "X_0(8)"), (8, 1, "X_1(8)"), (5, 0, "X_0(5)"), (5, 1, "X_1(5)"),
                     (7, 0, "X_0(7)"), (7, 1, "X_1(7)"), (9, 0, "X_0(9)"), (9, 1, "X_1(9)"), (3, 0, "X_0(3)"), (16, 1, "X_1(16)")):
    table = {}
    for p in range(1, 4 * N * 4):          # residues p mod 16N (contains a complete system mod 4N)
        if p % 2 == 0 or gcd(p, N) != 1:
            continue
        sf = s_local(N, "floor", p % N, typ)
        ss = s_local(N, "surface", ((1 + p) // 4) % N, typ) if p % 4 == 3 else None
        table[p % (4 * N)] = (sf, ss)
    results[name] = table
    # compare with the note's table (section 6.3)
    for r, (sf, ss) in table.items():
        if name == "X_0(2)":
            want = (1, 3 if r % 4 == 3 else None)
        elif name == "X_0(4)":
            want = (2 if r % 4 == 3 else 0, (4 if r % 8 == 7 else 0) if r % 4 == 3 else None)
        elif name == "X_0(8)":
            want = (4 if r % 8 == 7 else 0, (4 if r % 8 == 7 else 0) if r % 4 == 3 else None)
        elif name == "X_1(8)":
            want = (4 if r % 8 == 7 else 0, (8 if r % 16 == 15 else 0) if r % 4 == 3 else None)
        elif name in ("X_0(5)", "X_0(7)", "X_0(3)"):
            l = int(name[4])
            v = 1 + legendre(-r, l)
            want = (v, v if r % 4 == 3 else None)
        elif name == "X_0(9)":
            v = 1 + legendre(-r, 3)
            want = (v, v if r % 4 == 3 else None)
        elif name in ("X_1(5)", "X_1(7)", "X_1(9)"):
            Nn = int(name[4])
            v = (Nn - 1) if r % Nn == Nn - 1 else 0     # 4 for N = 5, 6 for N = 7, 9 (two eigenlines)
            if Nn == 9:
                v = 6 if r % 9 == 8 else 0
            want = (v, v if r % 4 == 3 else None)
        else:
            want = (sf, ss)
        if (sf, ss) != want:
            fail(f"{name} p = {r} mod {4 * N}: got {(sf, ss)}, note/prediction {want}")
for name, table in results.items():
    vals = sorted(set(table.values()), key=str)
    print(f"{name}: [s(floor), s(surface)] over all classes p mod {max(table) + 1 if table else 0}: distinct values {vals}")
print("X_1(8) surface by p mod 32:", {r: results['X_1(8)'][r][1] for r in sorted(results['X_1(8)']) if r % 4 == 3})
print("X_1(16) (beyond the note, for curiosity): floor/surface by p mod 64 (classes with nonzero entries):",
      {r: v for r, v in sorted(results['X_1(16)'].items()) if v[0] or (v[1] or 0)})
print("\nfailures:", bad)
print("ALL AUDIT CHECKS PASSED" if bad == 0 else "AUDIT FOUND FAILURES")
