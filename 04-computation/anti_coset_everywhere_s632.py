#!/usr/bin/env python3
"""
anti_coset_everywhere_s632.py

Look for the anti-coset pattern across the current research carriers.

Common schema:

    Aut(X)  = Iso(X, X)
    Anti(X) = Iso(X, JX)

where J is an opposite/converse/conjugation/transport operation.  Whenever
Anti(X) is nonempty it is a transporter, hence a coset/torsor over Aut(X).

This scout checks three finite faces:

1. self-converse tournaments, using the S630 exact generator through n=7;
2. LRC pair-sum shells under the <2,-1> action on Z/(2n-1);
3. Eisenstein/triangular point prefixes under the D6 lattice action.

It then records a proof-lens Tournament Analysis whose vertices are carrier
choices rather than tournament vertices.
"""

from __future__ import annotations

from collections import Counter, deque
import contextlib
import importlib.util
import io
import math
import os
from itertools import combinations


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULT_PATH = os.path.join(ROOT, "05-knowledge", "results", "anti_coset_everywhere_s632.out")


def load_s630():
    path = os.path.join(ROOT, "04-computation", "sc_complement_perspective_s630.py")
    spec = importlib.util.spec_from_file_location("sc_complement_perspective_s630", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def compose_perm(p: tuple[int, ...], q: tuple[int, ...]) -> tuple[int, ...]:
    """Composition for the old-for-new permutation convention used in S630."""
    return tuple(q[p[i]] for i in range(len(p)))


def is_perm_coset(
    auts: list[tuple[int, ...]], antis: list[tuple[int, ...]]
) -> tuple[bool, bool]:
    """Return right/left coset checks for a nonempty anti set."""
    if not antis:
        return False, False
    tau = antis[0]
    anti_set = set(antis)
    right = {compose_perm(a, tau) for a in auts}
    left = {compose_perm(tau, a) for a in auts}
    return right == anti_set, left == anti_set


def sc_tournament_atlas():
    s630 = load_s630()
    rows = []
    global_cycle_hist: Counter[tuple[int, ...]] = Counter()
    failures = []
    h_gap_hits = []

    for n in range(1, 8):
        classes = s630.sc_classes(n)
        aut_sizes = Counter()
        anti_sizes = Counter()
        perspective_orbit_swaps = Counter()
        for bits in classes:
            auts = s630.automorphisms(bits, n)
            antis = s630.anti_automorphisms(bits, n)
            aut_sizes[len(auts)] += 1
            anti_sizes[len(antis)] += 1
            global_cycle_hist.update(s630.cycle_type(a) for a in antis)
            right, left = is_perm_coset(auts, antis)
            if len(auts) != len(antis) or not right or not left:
                failures.append((n, bits, len(auts), len(antis), right, left))

            anti_images = s630.anti_image_sets(antis, n)
            perspective_orbit_swaps[len(set(anti_images))] += 1

            H = s630.hamiltonian_path_count(bits, n)
            if H in {7, 21}:
                h_gap_hits.append((n, bits, H))

        rows.append(
            {
                "n": n,
                "classes": len(classes),
                "aut_sizes": dict(sorted(aut_sizes.items())),
                "anti_sizes": dict(sorted(anti_sizes.items())),
                "orbit_swap_counts": dict(sorted(perspective_orbit_swaps.items())),
            }
        )

    return rows, dict(global_cycle_hist), failures, h_gap_hits


def fold_residue(x: int, mod: int) -> int:
    x %= mod
    if x == 0:
        raise ValueError("zero residue has no folded shell")
    return min(x, mod - x)


def generated_mult_group(mod: int, gens: list[int]) -> list[int]:
    seen = {1}
    q = deque([1])
    while q:
        a = q.popleft()
        for g in gens:
            b = (a * g) % mod
            if b not in seen:
                seen.add(b)
                q.append(b)
    return sorted(seen)


def lrc_shell_atlas():
    rows = []
    detail_n14 = []
    for n in (6, 8, 10, 12, 14, 16, 18, 20):
        mod = 2 * n - 1
        group = generated_mult_group(mod, [2, mod - 1])
        shells = list(range(1, n))

        shell_orbits = []
        seen = set()
        for shell in shells:
            if shell in seen:
                continue
            orbit = sorted({fold_residue(g * shell, mod) for g in group})
            seen.update(orbit)
            shell_orbits.append(orbit)

        transporter_ok = True
        transporter_sizes = Counter()
        reflection_collapses = 0
        for shell in shells:
            aut = [g for g in group if fold_residue(g * shell, mod) == shell]
            doubled = fold_residue(2 * shell, mod)
            anti2 = [g for g in group if fold_residue(g * shell, mod) == doubled]
            reflected = fold_residue(-shell, mod)
            anti_reflect = [g for g in group if fold_residue(g * shell, mod) == reflected]
            transporter_sizes[(len(aut), len(anti2))] += 1
            if len(aut) != len(anti2) or not anti2:
                transporter_ok = False
            if reflected == shell and set(anti_reflect) == set(aut):
                reflection_collapses += 1

        rows.append(
            {
                "n": n,
                "mod": mod,
                "group_size": len(group),
                "raw_shells": len(shells),
                "orbit_count": len(shell_orbits),
                "orbit_gcds": [math.gcd(orb[0], mod) for orb in shell_orbits],
                "transporter_sizes": dict(sorted(transporter_sizes.items())),
                "transporter_ok": transporter_ok,
                "reflection_collapses": reflection_collapses,
            }
        )

        if n == 14:
            for orbit in shell_orbits:
                detail_n14.append(
                    {
                        "orbit": orbit,
                        "gcd": math.gcd(orbit[0], mod),
                        "stabilizer": [
                            g for g in group if fold_residue(g * orbit[0], mod) == orbit[0]
                        ],
                        "to_double": [
                            g
                            for g in group
                            if fold_residue(g * orbit[0], mod)
                            == fold_residue(2 * orbit[0], mod)
                        ],
                    }
                )

    return rows, detail_n14


def hex_norm(p: tuple[int, int]) -> int:
    q, r = p
    return max(abs(q), abs(r), abs(q + r))


def axial_to_xy(p: tuple[int, int]) -> tuple[float, float]:
    q, r = p
    return (q + 0.5 * r, (math.sqrt(3) / 2) * r)


def eisenstein_order(limit_shell: int = 3) -> list[tuple[int, int]]:
    out = [(0, 0)]
    for radius in range(1, limit_shell + 1):
        shell = [
            (q, r)
            for q in range(-radius, radius + 1)
            for r in range(-radius, radius + 1)
            if hex_norm((q, r)) == radius
        ]
        shell.sort(key=lambda p: (math.atan2(axial_to_xy(p)[1], axial_to_xy(p)[0]), p))
        out.extend(shell)
    return out


def rot(i: int, p: tuple[int, int]) -> tuple[int, int]:
    q, r = p
    rots = (
        (q, r),
        (-r, q + r),
        (-q - r, q),
        (-q, -r),
        (r, -q - r),
        (q + r, -q),
    )
    return rots[i % 6]


def reflect(p: tuple[int, int]) -> tuple[int, int]:
    q, r = p
    return (r, q)


def d6_transforms():
    transforms = []
    for i in range(6):
        transforms.append((f"r{i}", +1, lambda p, i=i: rot(i, p)))
    for i in range(6):
        transforms.append((f"s r{i}", -1, lambda p, i=i: reflect(rot(i, p))))
    return transforms


def transform_key(fn) -> tuple[tuple[int, int], tuple[int, int]]:
    return fn((1, 0)), fn((0, 1))


def compose_transform(f, g):
    return lambda p: f(g(p))


def unit_edges(points: set[tuple[int, int]]) -> int:
    count = 0
    for a, b in combinations(points, 2):
        if hex_norm((a[0] - b[0], a[1] - b[1])) == 1:
            count += 1
    return count


def eisenstein_prefix_atlas():
    order = eisenstein_order(3)
    transforms = d6_transforms()
    rows = []
    for n in range(1, 22):
        points = set(order[:n])
        aut_plus = [(name, fn) for name, sign, fn in transforms if sign == 1 and {fn(p) for p in points} == points]
        anti = [(name, fn) for name, sign, fn in transforms if sign == -1 and {fn(p) for p in points} == points]

        torsor_ok = False
        if anti:
            anti_keys = {transform_key(fn) for _, fn in anti}
            tau = anti[0][1]
            coset_keys = {
                transform_key(compose_transform(tau, aut_fn)) for _, aut_fn in aut_plus
            }
            torsor_ok = anti_keys == coset_keys and len(anti) == len(aut_plus)

        rows.append(
            {
                "n": n,
                "edges": unit_edges(points),
                "aut_plus": len(aut_plus),
                "anti": len(anti),
                "torsor_ok": torsor_ok if anti else None,
                "full_shell": n in {1, 7, 19},
                "anti_names": [name for name, _ in anti],
            }
        )
    return rows


def directed_three_cycles(adj: list[list[int]]) -> int:
    n = len(adj)
    total = 0
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            total += 1
        elif adj[a][c] and adj[c][b] and adj[b][a]:
            total += 1
    return total


def ham_paths(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            if not dp[mask][v]:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[-1])


def proof_lens_tournament():
    vertices = [
        "SC Anti(T)=Iso(T,T^op)",
        "LRC shell transporter",
        "Eisenstein D6 conjugation",
        "coimage sign-reversing involution",
        "raw complement merge",
        "raw scalar 7/21 match",
    ]
    gauges = [
        [
            "SC Anti(T)=Iso(T,T^op)",
            "LRC shell transporter",
            "Eisenstein D6 conjugation",
            "coimage sign-reversing involution",
            "raw complement merge",
            "raw scalar 7/21 match",
        ],
        [
            "LRC shell transporter",
            "coimage sign-reversing involution",
            "SC Anti(T)=Iso(T,T^op)",
            "Eisenstein D6 conjugation",
            "raw complement merge",
            "raw scalar 7/21 match",
        ],
        [
            "Eisenstein D6 conjugation",
            "SC Anti(T)=Iso(T,T^op)",
            "LRC shell transporter",
            "coimage sign-reversing involution",
            "raw scalar 7/21 match",
            "raw complement merge",
        ],
    ]
    ranks = [{name: len(vertices) - i for i, name in enumerate(g)} for g in gauges]
    adj = [[0] * len(vertices) for _ in vertices]
    for i, v in enumerate(vertices):
        for j, w in enumerate(vertices):
            if i == j:
                continue
            if sum(r[v] > r[w] for r in ranks) >= 2:
                adj[i][j] = 1
    scores = Counter(sum(row) for row in adj)
    ordered = sorted(((sum(adj[i]), v) for i, v in enumerate(vertices)), reverse=True)
    return vertices, dict(sorted(scores.items())), directed_three_cycles(adj), ham_paths(adj), ordered


def main():
    print("S632 Anti-Coset Everywhere Atlas")
    print("================================")
    print()
    print("Definition")
    print("----------")
    print("For a carrier X and a dual/opposite/conjugation operation J,")
    print("Anti(X)=Iso(X,JX). If nonempty, Anti(X) is a transporter/coset")
    print("over Aut(X)=Iso(X,X). The recurring mistake is to choose one")
    print("anti-map and forget the coset, or to quotient until the anti-map")
    print("looks like an ordinary stabilizer.")
    print()

    print("1. Self-converse tournament face")
    print("--------------------------------")
    sc_rows, cycle_hist, failures, h_gap_hits = sc_tournament_atlas()
    print("n  SC classes  |Aut| histogram       |Anti| histogram      orbit-image counts")
    for row in sc_rows:
        print(
            f"{row['n']:1d}  {row['classes']:10d}  {row['aut_sizes']!s:20s} "
            f"{row['anti_sizes']!s:20s} {row['orbit_swap_counts']}"
        )
    print(f"coset failures: {failures}")
    print(f"H=7/21 hits inside checked SC classes: {h_gap_hits}")
    print(f"global anti-cycle-type histogram: {dict(sorted(cycle_hist.items()))}")
    print()

    print("2. LRC shell transporter face")
    print("-----------------------------")
    lrc_rows, detail_n14 = lrc_shell_atlas()
    print("n  C=2n-1  |<2,-1>|  shells->orbits  gcd reps       transporter sizes")
    for row in lrc_rows:
        print(
            f"{row['n']:2d} {row['mod']:7d} {row['group_size']:9d} "
            f"{row['raw_shells']:3d}->{row['orbit_count']:<3d} "
            f"{row['orbit_gcds']!s:14s} {row['transporter_sizes']} "
            f"ok={row['transporter_ok']} reflection-folds={row['reflection_collapses']}"
        )
    print("n=14 C=27 detail: folded shells collapse to gcd strata.")
    for item in detail_n14:
        print(
            f"  gcd={item['gcd']}: orbit={item['orbit']} "
            f"Stab(rep)={item['stabilizer']} Anti_2(rep)={item['to_double']}"
        )
    print("Note: after folding j~-j, reflection anti-data becomes a stabilizer;")
    print("the quotient has destroyed the orientation side channel.")
    print()

    print("3. Eisenstein/unit-distance D6 face")
    print("-----------------------------------")
    e_rows = eisenstein_prefix_atlas()
    print("n  unit_edges  Aut+  Anti  torsor  full_shell  anti maps")
    for row in e_rows:
        marker = "*" if row["full_shell"] else " "
        print(
            f"{row['n']:2d}{marker} {row['edges']:10d} {row['aut_plus']:5d} "
            f"{row['anti']:5d} {str(row['torsor_ok']):7s} "
            f"{str(row['full_shell']):10s} {row['anti_names']}"
        )
    print("Full centered hex shells (n=1,7,19) keep the full D6 anti-coset.")
    print("Partial shell prefixes usually retain only a small or empty anti-coset,")
    print("which is the finite carrier version of the 'Eisenstein pattern ends'")
    print("warning beyond the first complete shell.")
    print()

    print("4. Tournament Analysis over anti-coset lenses")
    print("---------------------------------------------")
    print("Vertices are carrier lenses. Pairwise observable: which lens keeps")
    print("the transporter, side-channel, and quotient-loss data most explicitly.")
    print("Switches/gauges: SC torsor theorem, LRC shell action, Eisenstein D6")
    print("conjugation. Tie HP: majority tournament over the three gauges.")
    vertices, score_hist, tri, hp, ordered = proof_lens_tournament()
    print(f"score_hist={score_hist} directed_3cycles={tri} H={hp}")
    for score, vertex in ordered:
        print(f"  score={score}: {vertex}")
    print()

    print("Assumption Challenge")
    print("--------------------")
    print("Alternate vertices considered: tournament vertices, anti-maps,")
    print("rooted perspectives, residue shells, folded shell orbits, unit-distance")
    print("points, D6 symmetries, coimage cancellations, and proof obligations.")
    print("Chosen vertices: carrier lenses. They preserve the predicate")
    print("'Anti(X) is a transporter to JX and quotienting can hide it.'")
    print("They destroy exact embeddings, complete LRC certificates, and raw")
    print("Hamiltonian-path spectra beyond the checked SC range.")
    print("Challenged assumption: an anti-coset is special to SC tournaments.")
    print("It is the generic transporter shape behind complements, converses,")
    print("reflections, conjugations, and sign-reversing cancellations.")
    print()

    print("Interpretation")
    print("--------------")
    print("- SC tournaments give the clean theorem: Anti(T) is a coset over Aut(T).")
    print("- LRC shell work already uses the same object as transporters under")
    print("  doubling/reflection on Z/(2n-1); folding can erase reflection into Aut.")
    print("- Eisenstein/unit-distance finite point sets keep Anti as the orientation")
    print("  reversing D6 coset exactly at full centered shells and lose pieces on")
    print("  partial-shell prefixes.")
    print("- The practical next test is to replace every raw complement merge in the")
    print("  repo by an explicit transporter: name X, name J, compute Aut(X), compute")
    print("  Anti(X), and record what the quotient forgot.")


def run_and_store():
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        main()
    text = buf.getvalue()
    print(text, end="")
    with open(RESULT_PATH, "w", encoding="utf-8") as f:
        f.write(text)


if __name__ == "__main__":
    run_and_store()
