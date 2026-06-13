#!/usr/bin/env python3
"""
lrc_gabor_zero_column_s542b.py

Finite Gabor trienerments for the Lonely Runner sector model.

The upstream S539/S541 Gabor notes posed a trienerment on atoms
(sector, harmonic).  This script makes one exact finite version.

Input state:
    c = (c_0, ..., c_{n-1}),  sum c_i = n-1,
    the sector occupancy vector from the fixed n-sector circle.

Two-sector finite Gabor window:
    G_c(a,m) = c_a zeta^(ma) + c_{a+1} zeta^(m(a+1)),
    where a is a section-boundary window and m is a harmonic.

The observer danger window is a=n-1, covering sectors n-1 and 0.
Thus the open LRC sector condition c_0=c_{n-1}=0 is exactly:

    G_c(n-1,m) = 0 for every harmonic m.

This converts the runner-space local trienerment condition
"observer has no near ties" into a Gabor-space marked zero column.  The
Gabor trienerment below puts vertices on atoms (a,m).  Edges are directed by
phase half-turn order; ties mean phase-unresolved atoms, either because one
coefficient vanishes or because two phases lie inside one atom-resolution
cell.

The goal is not to prove LRC.  It is to test whether this Gabor lift is a
meaningfully restricted isomorphism-class target rather than another copy of
the full sector-vector simplex.
"""

from __future__ import annotations

import cmath
import math
from collections import Counter
from functools import lru_cache, reduce
from itertools import combinations
from math import gcd
from statistics import mean

EPS = 1e-10


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for k in range(total + 1):
        for rest in compositions(total - k, parts - 1):
            yield (k,) + rest


def dihedral_canon_vector(c: tuple[int, ...]) -> tuple[int, ...]:
    n = len(c)
    images = []
    for r in range(n):
        images.append(tuple(c[(i + r) % n] for i in range(n)))
        images.append(tuple(c[(r - i) % n] for i in range(n)))
    return min(images)


def frac_num_den(num: int, den: int):
    return (num % den, den)


def frac_float(num: int, den: int) -> float:
    return (num % den) / den


def circular_distance_float(a: float, b: float) -> float:
    d = abs(a - b) % 1.0
    return min(d, 1.0 - d)


def sector_vector_at_fraction(speeds: tuple[int, ...], n: int, num: int, den: int) -> tuple[int, ...]:
    counts = [0] * n
    for v in speeds:
        # position = frac(v * num / den), sector = floor(n * position).
        rem = (v * num) % den
        sector = (n * rem) // den
        # Midpoints never lie on walls in this script, but clamp defensively.
        if sector == n:
            sector = n - 1
        counts[sector] += 1
    return tuple(counts)


def chamber_midpoints(speeds: tuple[int, ...], n: int) -> list[tuple[int, int]]:
    walls = {(0, 1), (1, 1)}
    for v in speeds:
        den = n * v
        for k in range(den + 1):
            g = gcd(k, den)
            walls.add((k // g, den // g))
    ordered = sorted(walls, key=lambda x: x[0] / x[1])
    mids = []
    for (a_num, a_den), (b_num, b_den) in zip(ordered, ordered[1:]):
        num = a_num * b_den + b_num * a_den
        den = 2 * a_den * b_den
        g = gcd(num, den)
        mids.append((num // g, den // g))
    return mids


def closed_lonely_at_fraction(speeds: tuple[int, ...], n: int, num: int, den: int) -> bool:
    for v in speeds:
        rem = (v * num) % den
        d = min(rem, den - rem)
        if n * d < den:
            return False
    return True


def wall_times(speeds: tuple[int, ...], n: int) -> list[tuple[int, int]]:
    walls = {(0, 1), (1, 1)}
    for v in speeds:
        den = n * v
        for k in range(den + 1):
            g = gcd(k, den)
            walls.add((k // g, den // g))
    return sorted(walls, key=lambda x: x[0] / x[1])


def sector_menu(speeds: tuple[int, ...], n: int):
    out = []
    for num, den in chamber_midpoints(speeds, n):
        out.append((sector_vector_at_fraction(speeds, n, num, den), num, den))
    return out


def is_good_sector(c: tuple[int, ...]) -> bool:
    return c[0] == 0 and c[-1] == 0


def finite_gabor_coeffs(c: tuple[int, ...]) -> tuple[complex, ...]:
    n = len(c)
    zeta = cmath.exp(2j * math.pi / n)
    vals = []
    for a in range(n):
        b = (a + 1) % n
        for m in range(n):
            vals.append(c[a] * (zeta ** (m * a)) + c[b] * (zeta ** (m * b)))
    return tuple(vals)


def observer_zero_atoms(c: tuple[int, ...]) -> int:
    n = len(c)
    vals = finite_gabor_coeffs(c)
    start = (n - 1) * n
    return sum(1 for z in vals[start : start + n] if abs(z) < EPS)


def gabor_phase_trienerment(c: tuple[int, ...]) -> tuple[int, ...]:
    """Return labeled trienerment over atom order (a,m), a,m=0..n-1.

    State over i<j:
      0  tie/unresolved;
      1  i -> j;
     -1  j -> i.
    """
    n = len(c)
    vals = finite_gabor_coeffs(c)
    phases = []
    for z in vals:
        if abs(z) < EPS:
            phases.append(None)
        else:
            phases.append((math.atan2(z.imag, z.real) / (2 * math.pi)) % 1.0)

    # One angular cell in the n x n atom torus.
    phase_tie_threshold = 1.0 / (n * n)
    states = []
    for i in range(n * n):
        pi = phases[i]
        for j in range(i + 1, n * n):
            pj = phases[j]
            if pi is None or pj is None:
                states.append(0)
                continue
            delta = (pj - pi) % 1.0
            if min(delta, 1.0 - delta) < phase_tie_threshold:
                states.append(0)
            elif delta < 0.5:
                states.append(1)
            else:
                states.append(-1)
    return tuple(states)


def states_to_matrix(states: tuple[int, ...], N: int) -> list[list[int]]:
    M = [[0] * N for _ in range(N)]
    it = iter(states)
    for i in range(N):
        for j in range(i + 1, N):
            s = next(it)
            M[i][j] = s
            M[j][i] = -s if s else 0
    return M


def scaffold_permutations(n: int):
    """Dihedral permutations of sector-window coordinate, with harmonic sign on reflection."""
    def idx(a: int, m: int) -> int:
        return a * n + m

    for r in range(n):
        yield [idx((a + r) % n, m) for a in range(n) for m in range(n)]
    for r in range(n):
        yield [idx((r - a - 1) % n, (-m) % n) for a in range(n) for m in range(n)]


@lru_cache(maxsize=None)
def gabor_class(c: tuple[int, ...]) -> tuple[int, ...]:
    n = len(c)
    N = n * n
    M = states_to_matrix(gabor_phase_trienerment(c), N)
    best = None
    for perm in scaffold_permutations(n):
        seq = tuple(M[perm[i]][perm[j]] for i in range(N) for j in range(i + 1, N))
        if best is None or seq < best:
            best = seq
    return best


def trienerment_fingerprint(states: tuple[int, ...], N: int) -> dict[str, object]:
    M = states_to_matrix(states, N)
    ties = sum(1 for s in states if s == 0)
    out_hist = Counter()
    tie_hist = Counter()
    for i in range(N):
        out = sum(1 for j in range(N) if M[i][j] == 1)
        tie = sum(1 for j in range(N) if j != i and M[i][j] == 0)
        out_hist[out] += 1
        tie_hist[tie] += 1

    directed_triangles = 0
    for i, j, k in combinations(range(N), 3):
        if M[i][j] == 0 or M[i][k] == 0 or M[j][k] == 0:
            continue
        outs = [
            (M[i][j] == 1) + (M[i][k] == 1),
            (M[j][i] == 1) + (M[j][k] == 1),
            (M[k][i] == 1) + (M[k][j] == 1),
        ]
        if sorted(outs) == [1, 1, 1]:
            directed_triangles += 1

    return {
        "ties": ties,
        "out_hist": tuple(sorted(out_hist.items())),
        "tie_hist": tuple(sorted(tie_hist.items())),
        "directed_triangles": directed_triangles,
    }


def runner_trienerment_counts(speeds: tuple[int, ...], n: int, num: int, den: int) -> tuple[int, int]:
    pts = [0.0] + [frac_float(v * num, den) for v in speeds]
    ties = 0
    observer = 0
    for i, j in combinations(range(n), 2):
        near = circular_distance_float(pts[i], pts[j]) < 1.0 / n - 1e-12
        if near:
            ties += 1
            if i == 0 or j == 0:
                observer += 1
    return ties, observer


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2:
        return None
    mx, my = mean(xs), mean(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx <= EPS or vy <= EPS:
        return None
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)


def analyze_global(n: int):
    vectors = list(compositions(n - 1, n))
    raw_orbits = {dihedral_canon_vector(c) for c in vectors}
    classes = {gabor_class(c) for c in vectors}
    good_vectors = [c for c in vectors if is_good_sector(c)]
    good_classes = {gabor_class(c) for c in good_vectors}

    tie_counts = []
    good_tie_counts = []
    zero_obs_hist = Counter()
    for c in vectors:
        states = gabor_phase_trienerment(c)
        fp = trienerment_fingerprint(states, n * n)
        tie_counts.append(fp["ties"])
        zero_obs_hist[observer_zero_atoms(c)] += 1
        if is_good_sector(c):
            good_tie_counts.append(fp["ties"])

    return {
        "raw": len(vectors),
        "raw_orbits": len(raw_orbits),
        "gabor_classes": len(classes),
        "good_raw": len(good_vectors),
        "good_gabor_classes": len(good_classes),
        "tie_range": (min(tie_counts), max(tie_counts)),
        "good_tie_range": (min(good_tie_counts), max(good_tie_counts)) if good_tie_counts else None,
        "observer_zero_hist": tuple(sorted(zero_obs_hist.items())),
    }


def primitive_speed_sets(n: int, B: int):
    for speeds in combinations(range(1, B + 1), n - 1):
        if reduce(gcd, speeds) == 1:
            yield speeds


def analyze_clocks(n: int, B: int):
    global_classes = {gabor_class(c) for c in compositions(n - 1, n)}
    global_good_classes = {gabor_class(c) for c in compositions(n - 1, n) if is_good_sector(c)}
    seen_classes = set()
    seen_good_classes = set()
    menu_sizes = []
    open_good = 0
    wall_only = 0
    missing = 0
    runner_tie_samples = []
    gabor_tie_samples = []
    lonely_gabor_ties = []
    nonlonely_gabor_ties = []
    ap_record = None

    sets = list(primitive_speed_sets(n, B))
    for speeds in sets:
        menu = sector_menu(speeds, n)
        local_classes = set()
        has_open = False
        has_closed = False
        for c, num, den in menu:
            cls = gabor_class(c)
            local_classes.add(cls)
            seen_classes.add(cls)
            if is_good_sector(c):
                has_open = True
                seen_good_classes.add(cls)

            rties, oties = runner_trienerment_counts(speeds, n, num, den)
            gtie = trienerment_fingerprint(gabor_phase_trienerment(c), n * n)["ties"]
            runner_tie_samples.append(rties)
            gabor_tie_samples.append(gtie)
            if oties == 0:
                lonely_gabor_ties.append(gtie)
            else:
                nonlonely_gabor_ties.append(gtie)

        for num, den in wall_times(speeds, n):
            if closed_lonely_at_fraction(speeds, n, num, den):
                has_closed = True
                break

        if has_open:
            open_good += 1
        elif has_closed:
            wall_only += 1
        else:
            missing += 1
        menu_sizes.append(len(local_classes))

        if speeds == tuple(range(1, n)):
            ap_record = {
                "menu_classes": len(local_classes),
                "open_good": has_open,
                "closed": has_closed,
            }

    corr = pearson([float(x) for x in runner_tie_samples], [float(y) for y in gabor_tie_samples])
    return {
        "sets": len(sets),
        "classes_seen": len(seen_classes),
        "classes_global": len(global_classes),
        "good_classes_seen": len(seen_good_classes),
        "good_classes_global": len(global_good_classes),
        "open_good": open_good,
        "wall_only": wall_only,
        "missing": missing,
        "menu_min_avg_max": (min(menu_sizes), mean(menu_sizes), max(menu_sizes)) if menu_sizes else None,
        "runner_gabor_tie_corr": corr,
        "lonely_gabor_ties_avg": mean(lonely_gabor_ties) if lonely_gabor_ties else None,
        "nonlonely_gabor_ties_avg": mean(nonlonely_gabor_ties) if nonlonely_gabor_ties else None,
        "ap_record": ap_record,
    }


def main():
    print("Finite Gabor trienerments for LRC -- codex S542")
    print("=" * 78)
    print("Atoms are (two-sector window a, harmonic m).")
    print("Observer lonely in sectors iff the observer Gabor column is all zero.")
    print("Edges: phase half-turn order; ties: zero/unresolved or same atom-angle cell.")
    print()

    print("GLOBAL SECTOR-VECTOR IMAGE")
    print("-" * 78)
    for n in range(3, 7):
        g = analyze_global(n)
        print(f"n={n}: raw={g['raw']} raw_dihedral_orbits={g['raw_orbits']} "
              f"Gabor scaffold-classes={g['gabor_classes']}")
        print(f"     good_raw={g['good_raw']} good_Gabor_classes={g['good_gabor_classes']} "
              f"observer_zero_atoms_hist={g['observer_zero_hist']}")
        print(f"     Gabor tie range all={g['tie_range']} good={g['good_tie_range']}")
    print()

    print("FIXED-CLOCK MENUS")
    print("-" * 78)
    for n, B in [(4, 10), (5, 8), (6, 7)]:
        a = analyze_clocks(n, B)
        mn, av, mx = a["menu_min_avg_max"]
        print(f"n={n}, B<={B}: speed_sets={a['sets']} "
              f"open_good={a['open_good']} wall_only={a['wall_only']} missing={a['missing']}")
        print(f"     Gabor classes seen={a['classes_seen']}/{a['classes_global']} "
              f"good seen={a['good_classes_seen']}/{a['good_classes_global']}")
        print(f"     menu classes min/avg/max={mn}/{av:.2f}/{mx}; "
              f"runner-tie vs Gabor-tie corr={a['runner_gabor_tie_corr']:.3f}")
        print(f"     avg Gabor ties: observer-tie-free={a['lonely_gabor_ties_avg']:.1f}, "
              f"observer-tied={a['nonlonely_gabor_ties_avg']:.1f}")
        if a["ap_record"]:
            print(f"     AP row: {a['ap_record']}")
    print()

    print("INTERPRETATION")
    print("-" * 78)
    print("1. HYP-2028 survives the lift: every sector-vector gives a finite")
    print("   Gabor trienerment, so global existence is still not the obstruction.")
    print("2. The LRC target becomes sharper: not a favorite sector vector, but a")
    print("   marked zero column in the two-sector Gabor atom grid.")
    print("3. Runner trienerments and Gabor trienerments use ties oppositely:")
    print("   observer loneliness deletes real near-ties, while it creates a large")
    print("   unresolved zero block among observer-window Gabor atoms.")
    print("4. The bounded fixed-clock menus see only part of the global Gabor-class")
    print("   image, so the Gabor angle is useful only as a fixed-clock menu/fiber")
    print("   restriction, not as global realizability.")


if __name__ == "__main__":
    main()
