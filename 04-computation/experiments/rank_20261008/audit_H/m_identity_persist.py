#!/usr/bin/env python3
"""audit H: THM-4611 proof of 2, escape word, bullet "identity runs end, because (M, e) != (1, 0) loses one digit of agreement
with (1, 0) per identity step".  Counterexample to the stated reason (not to the conclusion): on Z_5 (1,2,3,7,1) with the
standard constants (r_0 = 0), a state (M, 0) with M = 1 mod 5, M != 1, is FIXED by the digit 0 (identity coupling forever
along the word 000...), so the agreement with (1, 0) never drops.  Also: such states are reachable from integer offsets."""
from fractions import Fraction as Fr
from hcore import Map, residue
mp = Map(5, [1, 2, 3, 7, 1])
def agreement(M, e, cap=40):
    n = 0
    while n < cap and residue(M, 5 ** (n + 1)) == 1 and residue(e, 5 ** (n + 1)) == 0: n += 1
    return n
M, e = Fr(6), Fr(0)
print("start (M, e) = (6, 0): coupling (M mod 5, e mod 5) =", (residue(M, 5), residue(e, 5)), " agreement", agreement(M, e))
for t in range(6):
    M2, e2, i = mp.step(M, e, 0)
    print(f"  digit 0: u-digit {i}, new state ({M2}, {e2}), coupling {(residue(M2, 5), residue(e2, 5))}, agreement {agreement(M2, e2)}")
    M, e = M2, e2
# reachability of (M, 0) with M = 1 mod 5, M != 1 from integer starts (1, e0), words of length <= 4
found = []
for e0 in range(1, 200):
    level = [((), Fr(1), Fr(e0))]
    for depth in range(4):
        nxt = []
        for w, M, e in level:
            for j in range(5):
                M2, e2, i = mp.step(M, e, j)
                nxt.append((w + (j,), M2, e2))
                if e2 == 0 and M2 != 1 and residue(M2, 5) == 1:
                    found.append((e0, w + (j,), M2))
        level = nxt
    if found: break
print("reachable example (e0, digit word, M with e = 0, M = 1 mod 5, M != 1):", found[:3])
