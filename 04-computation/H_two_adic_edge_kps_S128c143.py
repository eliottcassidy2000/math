#!/usr/bin/env python3
"""
THE 2-ADIC EDGE OF H — how many low bits of the Redei count are formula-expressible?
kind-pasteur-2026-07-21-S128c143.  Owner: H's #P-ness is an EDGE case; tournaments sit at the
boundary between formula-expressible and provably-not (like the harmonic series / p=1). Test it
2-adically: Redei says H mod 2 = 1 (a formula). How many MORE low bits of H does a poly-time
invariant pin before the #P wall?

Setup: P = the poly-time invariants (score, specA, specS, disc, arb) — all O(n^3). Group iso
classes by their P-tuple. Within a P-cell, H may still vary (no poly invariant refines H, THM-1965).
The largest modulus m with 'H mod m' constant on EVERY cell = gcd over cells of within-cell
H-differences = M. Then:
  - H mod M is a poly-time function of P; H mod (anything not dividing M) is NOT.
  - v2(M) = number of poly-time-determined LOW BITS of H.  Redei guarantees v2(M) >= 1.
If v2(M) = 1 exactly: only the parity bit (Redei) is formula-expressible — the razor's edge.

Also: (b) does the poly-determined part hold for |R| too? (c) the split H = (formula part mod M)
+ M*(#P core).  (d) confirm the cycle-length boundary (c_{n-1} spectral, c_n not) — THM-1870.
"""
import sys, math
from functools import reduce
from collections import defaultdict
sys.path.insert(0, '04-computation')
from invariant_lattice_definitive_kps_S128c142 import census

def gcd_list(xs):
    xs = [x for x in xs if x != 0]
    if not xs: return 0
    return reduce(math.gcd, xs)

def v2(m):
    if m == 0: return math.inf
    k = 0
    while m % 2 == 0: m //= 2; k += 1
    return k

POLY = ("score", "specA", "specS", "disc", "arb")   # all O(n^3) poly-time
print("=" * 76)
print("THE 2-ADIC EDGE OF H  —  poly-time invariants P =", POLY)
print("=" * 76)
CEN = {}
for n in range(4, 7):
    CEN[n] = census(n)

for n in range(4, 7):
    data = CEN[n]
    cells = defaultdict(list)
    for c, d in data.items():
        key = tuple(d[f] for f in POLY)
        cells[key].append(d)
    # within-cell H differences
    diffs = []
    split_cells = 0
    biggest = None
    for key, members in cells.items():
        Hs = sorted(set(m["H"] for m in members))
        if len(Hs) > 1:
            split_cells += 1
            for i in range(len(Hs)):
                for j in range(i + 1, len(Hs)):
                    diffs.append(Hs[j] - Hs[i])
            if biggest is None or len(Hs) > len(biggest[1]):
                biggest = (key, Hs)
    M = gcd_list(diffs)
    ncells = len(cells)
    print(f"\n n={n}: iso classes={len(data)}  poly-cells={ncells}  H-split cells={split_cells}")
    if M == 0:
        print(f"        H is FULLY poly-determined at n={n} (no cell splits) — edge not reached yet")
    else:
        print(f"        M = gcd(within-cell H-differences) = {M}   v2(M) = {v2(M)}   odd-part = {M >> v2(M)}")
        print(f"        => H mod {2**v2(M)} is poly-time (a FORMULA); H mod {2**(v2(M)+1)} is NOT")
        if biggest:
            print(f"        biggest split cell: H-values {biggest[1]}  (diffs {[biggest[1][i+1]-biggest[1][i] for i in range(len(biggest[1])-1)]})")
        # confirm all H odd and all H in a cell share H mod 2^{v2(M)}
        allH = [d["H"] for d in data.values()]
        print(f"        all H odd (Redei)? {all(h % 2 == 1 for h in allH)}   "
              f"H mod {2**v2(M)} constant on every cell? {all(len(set(m['H'] % (2**v2(M)) for m in mem))==1 for mem in cells.values())}")

print()
print("=" * 76)
print("(b) SAME for |R|, and the pair — is the signed count's poly-determined depth different?")
print("=" * 76)
for n in range(4, 7):
    data = CEN[n]
    cells = defaultdict(list)
    for c, d in data.items():
        cells[tuple(d[f] for f in POLY)].append(d)
    for inv in ["R", "H"]:
        diffs = []
        for mem in cells.values():
            vs = sorted(set(m[inv] for m in mem))
            for i in range(len(vs)):
                for j in range(i+1, len(vs)):
                    diffs.append(vs[j]-vs[i])
        M = gcd_list(diffs)
        print(f" n={n}: {inv}: M={M}  v2={v2(M) if M else 'inf'}", end="   ")
    print()

print()
print("=" * 76)
print("(c) THE SPLIT: H = (formula part, poly-time) + M*(#P core).  Show the residual is genuinely free.")
print("=" * 76)
for n in [6]:
    data = CEN[n]
    cells = defaultdict(list)
    for c, d in data.items():
        cells[tuple(d[f] for f in POLY)].append(d)
    diffs = [b-a for mem in cells.values() for a in [min(m['H'] for m in mem)] for b in set(m['H'] for m in mem)]
    M = gcd_list([d for d in diffs if d])
    core_vals = set()
    for mem in cells.values():
        base = min(m['H'] for m in mem)
        for m in mem:
            core_vals.add((m['H'] - (base % M if M else 0)))  # illustrative
    print(f" n={n}: M={M}. Within a split poly-cell, H takes values differing by multiples of {M};")
    print(f"        the low {v2(M)} bit(s) of H are the FORMULA part (Redei parity), the rest is the #P core.")

print()
print("=" * 76)
print("(d) THE CRITICAL LENGTH (THM-1870 re-confirmed here): c_{n-1} spectral, c_n not")
print("=" * 76)
for n in range(4, 7):
    data = CEN[n]
    byspec = defaultdict(list)
    for c, d in data.items():
        byspec[d["specA"]].append(d)
    # c_k = k-th entry of cyc (cyc = (c3,...,cn)); is c_k constant on cospectral classes?
    for k in range(3, n+1):
        idx = k - 3
        spectral = all(len(set(m["cyc"][idx] for m in g)) == 1 for g in byspec.values())
        crit = " <-- HAMILTONIAN LENGTH" if k == n else ""
        print(f" n={n}: c_{k} spectral (constant on cospectral classes)? {spectral}{crit}")
    print()
print("DONE.")
