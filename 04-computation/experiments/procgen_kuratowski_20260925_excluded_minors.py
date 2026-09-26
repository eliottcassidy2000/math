#!/usr/bin/env python3
"""procgen_kuratowski_20260925, part M: an excluded-minor calculus for Collatz-type maps and proof methods.

Two-branch maps T_{q,r}(x) = x/2 (x even), (qx+r)/2 (x odd), q, r odd, on N (and on -N, which is T_{q,-r}
on N after conjugation by nu).  The rcwa minor order: (q, r') <= (q, r) iff r' | r (restriction to the
invariant set r'' N, r'' = r/r', followed by x -> x/r''), plus nu.
Sections (every check raises on failure):
  M1  minor relation: T_{q,r}(d x) = d T_{q,r/d}(x) for odd d | r (window check); scaled cycles of (q, r/d)
      are cycles of (q, r); the multiplier q is invariant, so different q are incomparable.
  M2  census q in {1,3,5,...,15}, r in {1,...,25}: cycles met from starts 1..3000 on N and on -N, escapes
      past 1e40; P_fin evidence; Conn status (single grand orbit on N).
  M3  Conn: for primes p not dividing 2q, pN and N \\ pN are both invariant, so (q,p) is never Conn;
      each q-slice therefore has at least one excluded minor and the excluded-minor set is infinite
      (an antichain across q) -- PROVED; in the slice q = 3 it is finite iff Collatz fails.
  M4  the three controls and their independence: C = (3,1), sigma C = (3,-1), delta C = (5,1), and
      epsilon C = THM-4470's pair-flip maps (a single flip creates a cycle; density-zero flips create a
      divergent orbit).  Invariants: drift sign, sheet sign b*side, rcwa (finite description).
Run: python3 04-computation/experiments/procgen_kuratowski_20260925_excluded_minors.py  (~1 min, < 300 MB)
"""
import math
from collections import Counter


def check(cond, msg):
    if not cond:
        raise SystemExit("CHECK FAILED: " + msg)


def Tqr(q, r):
    return lambda x: (q * x + r) // 2 if x % 2 else x // 2


def census(q, r, sign, X=3000, cap=10 ** 40, maxsteps=100000):
    F = Tqr(q, r)
    known = {}
    cycles = {}
    esc = 0
    for s0 in range(1, X + 1):
        x = sign * s0
        path = []
        pos = {}
        steps = 0
        res = None
        while True:
            if x in known:
                res = known[x]
                if res == "ESC":
                    esc += 1
                break
            if x in pos:
                cyc = path[pos[x]:]
                res = min(cyc, key=lambda z: (abs(z), z))
                cycles[res] = tuple(cyc)
                break
            if abs(x) > cap or steps > maxsteps:
                esc += 1
                res = "ESC"
                break
            pos[x] = len(path)
            path.append(x)
            x = F(x)
            steps += 1
        for y in path:
            known[y] = res
    return cycles, esc


print("=" * 100)
print("M1  the rcwa minor relation on two-branch maps")
print("=" * 100)
for q in (1, 3, 5, 7):
    for r in (3, 5, 9, 15, 21, 25):
        for d in (dd for dd in range(3, r + 1, 2) if r % dd == 0):
            F, G = Tqr(q, r), Tqr(q, r // d)
            check(all(F(d * x) == d * G(x) for x in range(-3000, 3001)), "T_{q,r}(dx) = d T_{q,r/d}(x)")
print("  T_{q,r}(d x) = d T_{q,r/d}(x) for every odd d | r (q = 1,3,5,7; r = 3,5,9,15,21,25; |x| <= 3000):")
print("  (q, r/d) is the restriction of (q, r) to the invariant set dZ, conjugated by x -> x/d")
print("  the multiplier q is a minor invariant (restriction and affine conjugation keep the affine slopes q/2, 1/2),")
print("  so the q-slices are pairwise incomparable: {(q, 1) : q odd} is an infinite ANTICHAIN (no wqo, no")
print("  Robertson-Seymour finiteness for this order).")

print()
print("=" * 100)
print("M2  census of two-branch maps: cycles from starts 1..3000 on N and on -N, escapes past 1e40")
print("=" * 100)
CEN = {}
for q in (1, 3, 5, 7, 9, 11, 13, 15):
    for r in range(1, 26, 2):
        cp, ep = census(q, r, 1)
        cn, en = census(q, r, -1)
        CEN[(q, r)] = (cp, ep, cn, en)
print("  q | r=1..25: (#cycles on N, #escapes) ; (#cycles on -N, #escapes)")
for q in (1, 3, 5, 7, 9, 11, 13, 15):
    row = []
    for r in range(1, 26, 2):
        cp, ep, cn, en = CEN[(q, r)]
        row.append(f"{len(cp)}/{ep}")
    print(f"  q={q:2d} N : " + " ".join(f"{s:>7s}" for s in row))
    row = []
    for r in range(1, 26, 2):
        cp, ep, cn, en = CEN[(q, r)]
        row.append(f"{len(cn)}/{en}")
    print(f"  q={q:2d} -N: " + " ".join(f"{s:>7s}" for s in row))
for r in range(1, 26, 2):
    check(CEN[(1, r)][1] == 0 and CEN[(1, r)][3] == 0, "q = 1: no escapes (contracting)")
    check(CEN[(3, r)][1] == 0 and CEN[(3, r)][3] == 0, "q = 3: no escapes in the census")
for q in (5, 7, 9, 11, 13, 15):
    check(all(CEN[(q, r)][1] > 0 for r in range(1, 26, 2)), "q >= 5: escapes on N")
print("  q = 1: every orbit enters a cycle (PROVED: (x+r)/2 < x for x > r); q = 3: no escape (Lagarias's 3x+k")
print("  conjectures, OPEN); q >= 5: most starts escape past 1e40 (divergence expected, never PROVED for any q >= 5)")
# minor-closure consistency: scaled cycles
cnt = 0
for q in (1, 3, 5):
    for r in range(3, 26, 2):
        F = Tqr(q, r)
        for d in (dd for dd in range(3, r + 1, 2) if r % dd == 0):
            for key, cyc in CEN[(q, r // d)][0].items():
                z = [d * t for t in cyc]
                check(all(F(z[i]) == z[(i + 1) % len(z)] for i in range(len(z))), "scaled cycle is a cycle")
                cnt += 1
print(f"  {cnt} scaled cycles d * (cycle of (q, r/d)) verified as cycles of (q, r): cycle sets grow along the order")
print("  (3,1) on N: cycles", sorted(CEN[(3, 1)][0]), "; on -N (= (3,-1) on N after nu):", sorted(CEN[(3, 1)][2]))
print("  (5,1) on N: cycles", sorted(CEN[(5, 1)][0]))
check(sorted(CEN[(3, 1)][2]) == [-17, -5, -1] and sorted(CEN[(5, 1)][0]) == [1, 13, 17], "control cycles")

print()
print("=" * 100)
print("M3  Conn (single grand orbit on N): excluded minors")
print("=" * 100)


def primes_upto(n):
    return [p for p in range(2, n + 1) if all(p % k for k in range(2, int(p ** 0.5) + 1))]


for q in (1, 3, 5, 7):
    for p in primes_upto(25):
        if p == 2 or q % p == 0:
            continue
        F = Tqr(q, p)
        for x in range(1, 20001):
            check((F(x) % p == 0) == (x % p == 0), "pN and its complement invariant")
print("  codes: X2 = two cycles found (NOT Conn), Xp = invariant pN with invariant complement (NOT Conn),")
print("         ok = one cycle and no escape (consistent with Conn), ?e = escapes (status open)")
print("  for primes p not dividing 2q, x in pN <=> T_{q,p}(x) in pN (checked q = 1,3,5,7, p <= 25, x <= 20000):")
print("  both pN and N \\ pN are unions of grand orbits, so (q, p) is never Conn (PROVED for all such q, p).")
status = {}
for q in (1, 3, 5, 7):
    for r in range(1, 26, 2):
        cp, ep, cn, en = CEN[(q, r)]
        if len(cp) >= 2:
            st = "NOT Conn (>= 2 cycles)"
        elif any(r % p == 0 for p in primes_upto(r) if p != 2 and q % p != 0 and p > 1):
            st = "NOT Conn (invariant pN)"
        elif ep > 0:
            st = "open (escapes)"
        elif len(cp) == 1:
            st = "consistent with Conn"
        else:
            st = "open"
        status[(q, r)] = st
exc = {}
for q in (1, 3, 5, 7):
    mins = []
    for r in range(1, 26, 2):
        if status[(q, r)].startswith("NOT") and not any(status[(q, d)].startswith("NOT")
                                                        for d in range(1, r, 2) if r % d == 0):
            mins.append(r)
    exc[q] = mins
    code = {"NOT Conn (>= 2 cycles)": "X2", "NOT Conn (invariant pN)": "Xp", "open (escapes)": "?e",
            "consistent with Conn": "ok", "open": "?"}
    print(f"  q = {q}: r=1..25 ->", " ".join(f"{r}:{code[status[(q, r)]]}" for r in range(1, 26, 2)))
    print(f"         minimal non-Conn r (excluded minors of the slice, r <= 25): {mins}")
print("  (q = 1: r = 1 is PROVED Conn; q = 3: r = 1 is Collatz; q = 7: r = 1 is open, so its list is conditional)")
check(exc[3][:6] == [5, 7, 11, 13, 17, 19], "q = 3 excluded minors are the primes >= 5 (given Collatz)")
check(exc[5] == [1], "(5,1) is the unique excluded minor of the q = 5 slice")
print("  q = 3: if Collatz holds, (3,1) is Conn and the excluded minors of the slice are {(3,p): p prime >= 5},")
print("  an infinite antichain; if Collatz fails, (3,1) itself is the unique excluded minor of the slice.")
print("  So 'the q = 3 slice has a finite Kuratowski set' <=> 'Collatz is false' -- a restatement, not a reduction.")
print("  The SHEET control (3,-1) on N and the DRIFT control (5,1) are ATOMS (r = 1) and non-Conn: genuine excluded")
print("  minors (PROVED by their extra cycles); they are incomparable (different sheet / different multiplier).")

print()
print("=" * 100)
print("M4  the three controls and their independence")
print("=" * 100)


def flip_map(flipped):
    """THM-4470 pairing family: pair i = {2i-1, 2i}; unflipped: T; flipped: 2i-1 -> i-1, 2i -> 3i."""
    def F(n):
        i = (n + 1) // 2
        if i in flipped:
            return i - 1 if n % 2 else 3 * i
        return (3 * n + 1) // 2 if n % 2 else n // 2
    return F


# single flip (pair 4) creates a new cycle
F4 = flip_map({4})
cyc4 = set()
for s in range(1, 3000):
    x = s
    seen = []
    while x not in seen and x > 0 and len(seen) < 10000:
        seen.append(x)
        x = F4(x)
    if x in seen:
        cyc4.add(tuple(sorted(seen[seen.index(x):])))
print("  epsilon (single flip of pair 4 = {7,8}): cycles met from starts < 3000:", sorted(cyc4))
check(len(cyc4) >= 2, "single flip creates a second cycle")
# density-zero flips: the divergent chain
c = [3]
for t in range(300):
    c.append(c[-1] + (c[-1] + 1) // 2)
flips = {(ct + 1) // 2 for ct in c if ct % 2 == 0}
Fe = flip_map(flips)
check(all(Fe(c[t]) == c[t + 1] for t in range(300)), "the modified map sends c_t -> c_{t+1}")
Xd = 10 ** 30
dens = sum(1 for i in flips if 2 * i <= Xd)
print(f"  epsilon (density-zero flips, THM-4470(4)): c_0 = 3, c_(t+1) = c_t + ceil(c_t/2), flip the pair of every even")
print(f"  c_t: the orbit of 3 is c_0 < c_1 < ... (300 steps verified, c_300 ~ 10^{len(str(c[300])) - 1}); flipped pairs")
print(f"  below 10^30: {dens} (density zero)")
print("  epsilon C is not rcwa: an rcwa map agreeing with T on a density-one set agrees with T on every residue")
print("  class where both are affine, i.e. everywhere; epsilon C differs from T on the flipped pairs.")
rows = [
    ("C = (3,+1) on N", math.log(3) - 2 * math.log(2), +1, True, "Collatz (OPEN)"),
    ("sigma C = (3,-1) on N", math.log(3) - 2 * math.log(2), -1, True, "NOT Conn: 3 cycles"),
    ("delta C = (5,+1) on N", math.log(5) - 2 * math.log(2), +1, True, "NOT Conn: 3 cycles"),
    ("epsilon C = flips of T", math.log(3) - 2 * math.log(2), +1, False, "NOT Conn: new cycle / divergent orbit"),
]
print()
print(f"  {'system':26s} {'drift sign':>10s} {'b*side':>7s} {'rcwa':>5s}   status")
for name, dr, bs, rc, st in rows:
    print(f"  {name:26s} {('-' if dr < 0 else '+'):>10s} {bs:+7d} {str(rc):>5s}   {st}")
sep = {"drift sign": 1, "b*side": 2, "rcwa": 3}
base = rows[0]
table = {}
for k, idx in sep.items():
    table[k] = [(r[idx] if idx != 1 else (r[1] < 0)) != (base[idx] if idx != 1 else (base[1] < 0)) for r in rows[1:]]
print()
print("  which invariant separates C from which control (True = separates):")
for k in sep:
    print(f"    {k:10s}: sigma {table[k][0]!s:5s}  delta {table[k][1]!s:5s}  epsilon {table[k][2]!s:5s}")
check(table["drift sign"] == [False, True, False] and table["b*side"] == [True, False, False]
      and table["rcwa"] == [False, False, True], "diagonal independence table")
print("  => each invariant separates exactly one control: the controls are pairwise independent (separating one")
print("     never forces separating another), and a sound method proving Collatz must see all three kinds of data:")
print("     the archimedean order (sheet), size with the multiplier (drift), and pointwise/non-density data (defect).")
print()
print("ALL CHECKS PASSED (excluded_minors)")
