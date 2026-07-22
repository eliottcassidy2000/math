#!/usr/bin/env python3
"""Exact controls for the corrected HYP-8840 constant-term thread.

The diagonal energy CT[P^m Pbar^m] sees only augmentation-zero relations.  The
full mixed table CT[P^r Pbar^s], or equivalently an observer at speed zero,
restores the augmentation coordinate needed to compare with LRC.  This script
checks the surviving identities and the decisive hostile examples; it proves
no GMC-to-LRC implication.
"""

from fractions import Fraction as F
from itertools import combinations, product
from math import comb, gcd


def section(title: str) -> None:
    print("\n" + "=" * 72)
    print(title)
    print("=" * 72)


def pmul(a: dict[int, int], b: dict[int, int]) -> dict[int, int]:
    out: dict[int, int] = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            out[ea + eb] = out.get(ea + eb, 0) + ca * cb
    return out


def ppow(a: dict[int, int], m: int) -> dict[int, int]:
    out = {0: 1}
    for _ in range(m):
        out = pmul(out, a)
    return out


def poly(speeds: tuple[int, ...]) -> dict[int, int]:
    return {v: 1 for v in speeds}


def bar(p: dict[int, int]) -> dict[int, int]:
    return {-e: c for e, c in p.items()}


def mixed_moment(speeds: tuple[int, ...], r: int, s: int) -> int:
    p = poly(speeds)
    return pmul(ppow(p, r), ppow(bar(p), s)).get(0, 0)


def direct_energy(speeds: tuple[int, ...], m: int) -> int:
    counts: dict[int, int] = {}
    for row in product(speeds, repeat=m):
        total = sum(row)
        counts[total] = counts.get(total, 0) + 1
    return sum(value * value for value in counts.values())


def schur_count(speeds: tuple[int, ...]) -> int:
    values = set(speeds)
    return sum(a + b in values for a in speeds for b in speeds)


def norm_residue(a: int, q: int) -> int:
    r = a % q
    return min(r, q - r)


def exact_lrc_max(speeds: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    candidates = {
        F(a, v + w) % 1
        for v in speeds
        for w in speeds
        for a in range(v + w)
    }
    values = {
        t: min(F(norm_residue(t.numerator * v, t.denominator), t.denominator) for v in speeds)
        for t in candidates
    }
    maximum = max(values.values())
    return maximum, tuple(sorted(t for t, value in values.items() if value == maximum))


section("P1  Universal constant-term identity")
for speeds, m in [((1, 2, 3, 4), 2), ((1, 2, 4, 8), 3), ((0, 1, 3, 7), 3)]:
    ct = mixed_moment(speeds, m, m)
    direct = direct_energy(speeds, m)
    assert ct == direct
    print(f"S={speeds}, m={m}: CT[P^m Pbar^m]={ct}=direct additive energy")
print("Survivor: the diagonal CT identity is exact; its relation support is affine (augmentation zero).")


section("P2  Finite AP-energy scout, with its actual scope")
for n, k in ((8, 5), (9, 5), (10, 6)):
    for m in (2, 3, 4):
        rows = list(combinations(range(n), k))
        scores = [(mixed_moment(row, m, m), row) for row in rows]
        best = max(score for score, _ in scores)
        winners = [row for score, row in scores if score == best]
        assert all(len({row[i + 1] - row[i] for i in range(k - 1)}) == 1 for row in winners)
        print(f"N={n}, k={k}, m={m}: best={best}; AP maximizers={len(winners)}/{len(rows)} rows")
print("Scope: nine finite censuses, not an all-k or all-order AP theorem; THM-730 counts a+b=c instead.")


section("P3  No-carry tensorization and the AP carry obstruction")
s1, s2, m, scale = (0, 1, 2), (0, 1, 2, 3), 2, 5
joined = tuple(a + scale * b for b in s2 for a in s1)
left = mixed_moment(joined, m, m)
right = mixed_moment(s1, m, m) * mixed_moment(s2, m, m)
assert left == right
print(f"no-carry scale {scale}: E2(S1+scale*S2)={left}=E2(S1)E2(S2)")
ap12 = tuple(range(12))
carry_left = mixed_moment(ap12, 2, 2)
carry_right = mixed_moment(s1, 2, 2) * mixed_moment(s2, 2, 2)
assert (carry_left, carry_right) == (1156, 836)
print(f"target AP at carry scale 3: E2([0..11])={carry_left}!={carry_right}=19*44")
print("Thus separated-digit factorization is conditional and removes the carry relations driving AP energy.")


section("P4  The missing augmentation coordinate")
a = (1, 2, 3, 4, 5)
b = (2, 3, 4, 5, 6)
c = (1, 3, 5, 7, 9)
for m in (2, 3, 4):
    energies = tuple(mixed_moment(row, m, m) for row in (a, b, c))
    assert energies[0] == energies[1] == energies[2]
    print(f"m={m}: unanchored energies A,B,C all equal {energies[0]}")
maxima = tuple(exact_lrc_max(row)[0] for row in (a, b, c))
assert maxima == (F(1, 6), F(1, 4), F(1, 2))
print(f"but exact LRC maxima are {maxima}")
anchored = tuple(mixed_moment((0,) + row, 2, 2) for row in (a, b, c))
schur = tuple(schur_count(row) for row in (a, b, c))
assert anchored == (146, 130, 106)
assert schur == (10, 6, 0)
assert tuple(mixed_moment(row, 2, 1) for row in (a, b, c)) == schur
print(f"anchored E2 values={anchored}; Schur counts={schur}=CT[P^2 Pbar]")
print("Repair: retain M_(r,s)=CT[P^r Pbar^s]; r-s is the augmentation sidecar.")


section("P5  Exact volume ceiling")
small = (1, 2, 3)
maximum, tops = exact_lrc_max(small)
assert maximum == F(1, 4)
assert tops == (F(1, 4), F(3, 4))
print(f"S={small}: M={maximum}; top fiber={tops}; Haar measure=0; Euler components={len(tops)}")
print("Positive safe measure certifies a strict exit; Euler characteristic also detects isolated equality.")


section("SUMMARY")
print("Exact survivors: diagonal and mixed CT identities, the finite AP scout, conditional no-carry")
print("tensorization, and the volume/Euler boundary distinction.  Refuted: diagonal energy as the full")
print("LRC lattice, THM-730 as an E2 theorem, arbitrary rank reduction, and a GMC radial-kernel swap.")
print("Live repair: the mixed (r,s) table or an observer at zero restores augmentation; no implication yet.")
