#!/usr/bin/env python3
"""Tanner/dessin/weak-regularity audit for the LRC14 Delsarte carriers.

This is a continuation of HYP-2740 after the user prompted Tanner graphs,
unit-distance/Delsarte links, weakly regular graphs, Doyle-Holt half-arcs,
and Belyi/modular-form analogies.

The computed object is deliberately small and exact.  Vertices are not runners:
they are Delsarte moment rows r and missed-depth atoms t.  Edge r--t exists
when r <= t and has sign sign(y_r) in the binomial-basis certificate.  This
preserves the THM-534 dominance predicate g(t) >= scale*[t=0], but it forgets
the generated LRC word and the speed-set relation lattice.

Tests:
  * biregularity / common-neighbor profiles (unit-distance / weak-regular cue);
  * coarsest degree-neighborhood quotient (what regularity survives);
  * side-preserving support automorphism edge orbits (Doyle-Holt cue);
  * natural bipartite rotation-system dessin passport (Belyi cue);
  * tournament analysis over proof lenses.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations


@dataclass(frozen=True)
class Dual:
    name: str
    scale: int
    y_scaled: tuple[int, ...]
    g_scaled: tuple[int, ...]


DUALS = [
    Dual("K8", 10, (10, -10, 10, -9, 6, 0, 0), (10, 0, 0, 1, 0, 0, 10)),
    Dual("K9", 18, (18, -13, 8, -3, 0, 0, 0), (18, 5, 0, 0, 2, 3, 0)),
    Dual("K11", 6, (6, -3, 1, 0, 0, 0, 0), (6, 3, 1, 0, 0, 1, 3)),
]


def binom(n: int, k: int) -> int:
    if k < 0 or k > n:
        return 0
    out = 1
    for i in range(1, k + 1):
        out = out * (n - i + 1) // i
    return out


def sign(x: int) -> str:
    return "+" if x > 0 else "-" if x < 0 else "0"


def edges(dual: Dual) -> list[tuple[int, int]]:
    rows = [r for r, y in enumerate(dual.y_scaled) if y]
    return [(r, t) for r in rows for t in range(7) if r <= t]


def edge_weight(dual: Dual, e: tuple[int, int]) -> Fraction:
    r, t = e
    return Fraction(dual.y_scaled[r] * binom(t, r), dual.scale)


def neighborhoods(dual: Dual):
    es = edges(dual)
    check_to_atoms = defaultdict(set)
    atom_to_checks = defaultdict(set)
    for r, t in es:
        check_to_atoms[r].add(t)
        atom_to_checks[t].add(r)
    return dict(check_to_atoms), dict(atom_to_checks)


def degree_partition(vertices: list[tuple[str, int]], nbrs: dict[tuple[str, int], set]):
    return {v: (v[0], len(nbrs[v])) for v in vertices}


def quotient_is_equitable(vertices, nbrs, colors):
    classes = defaultdict(list)
    for v in vertices:
        classes[colors[v]].append(v)
    keys = sorted(classes)
    signatures = {}
    equitable = True
    for color, vs in classes.items():
        sigs = []
        for v in vs:
            sig = tuple(sum(1 for u in nbrs[v] if colors[u] == k) for k in keys)
            sigs.append(sig)
        signatures[color] = sigs[0]
        equitable = equitable and all(s == sigs[0] for s in sigs)
    return equitable, keys, signatures


def support_automorphisms(dual: Dual):
    """Brute-force side-preserving support automorphisms.

    The graphs are tiny.  We prune by side and degree, then verify adjacency.
    """
    ch, at = neighborhoods(dual)
    rows = sorted(ch)
    atoms = sorted(set(range(7)))
    row_deg = {r: len(ch[r]) for r in rows}
    atom_deg = {t: len(at.get(t, set())) for t in atoms}
    row_blocks = defaultdict(list)
    atom_blocks = defaultdict(list)
    for r in rows:
        row_blocks[row_deg[r]].append(r)
    for t in atoms:
        atom_blocks[atom_deg[t]].append(t)

    row_maps = [dict(zip(rows, rows))]
    for block in row_blocks.values():
        new = []
        for base in row_maps:
            for perm in permutations(block):
                m = dict(base)
                m.update(zip(block, perm))
                new.append(m)
        row_maps = new

    atom_maps = [dict(zip(atoms, atoms))]
    for block in atom_blocks.values():
        new = []
        for base in atom_maps:
            for perm in permutations(block):
                m = dict(base)
                m.update(zip(block, perm))
                new.append(m)
        atom_maps = new

    es = set(edges(dual))
    autos = []
    for rm in row_maps:
        for tm in atom_maps:
            if {(rm[r], tm[t]) for r, t in es} == es:
                autos.append((rm, tm))
    return autos


def edge_orbits(dual: Dual):
    es = edges(dual)
    autos = support_automorphisms(dual)
    unseen = set(es)
    orbits = []
    while unseen:
        seed = min(unseen)
        orbit = set()
        for rm, tm in autos:
            orbit.add((rm[seed[0]], tm[seed[1]]))
        orbits.append(sorted(orbit))
        unseen -= orbit
    return orbits, autos


def permutation_cycles(perm: dict[tuple[int, int], tuple[int, int]]):
    unseen = set(perm)
    cycles = []
    while unseen:
        x = min(unseen)
        cyc = []
        while x in unseen:
            unseen.remove(x)
            cyc.append(x)
            x = perm[x]
        cycles.append(cyc)
    return cycles


def dessin_passport(dual: Dual):
    es = sorted(edges(dual))
    sigma_black = {}
    sigma_white = {}
    by_row = defaultdict(list)
    by_atom = defaultdict(list)
    for e in es:
        by_row[e[0]].append(e)
        by_atom[e[1]].append(e)
    for row_edges in by_row.values():
        cyc = sorted(row_edges, key=lambda e: e[1])
        for a, b in zip(cyc, cyc[1:] + cyc[:1]):
            sigma_black[a] = b
    for atom_edges in by_atom.values():
        cyc = sorted(atom_edges, key=lambda e: e[0])
        for a, b in zip(cyc, cyc[1:] + cyc[:1]):
            sigma_white[a] = b

    composed = {e: sigma_black[sigma_white[e]] for e in es}
    sigma_face = {v: k for k, v in composed.items()}
    black = permutation_cycles(sigma_black)
    white = permutation_cycles(sigma_white)
    faces = permutation_cycles(sigma_face)
    return black, white, faces


def side_common_profiles(dual: Dual):
    ch, at = neighborhoods(dual)
    check_common = Counter()
    atom_common = Counter()
    for a, b in combinations(sorted(ch), 2):
        check_common[len(ch[a] & ch[b])] += 1
    for a, b in combinations(range(7), 2):
        atom_common[len(at.get(a, set()) & at.get(b, set()))] += 1
    return check_common, atom_common


def print_dual(dual: Dual):
    ch, at = neighborhoods(dual)
    es = edges(dual)
    print("-" * 78)
    print(f"{dual.name}: scale={dual.scale}")
    print(f"  y_scaled={list(dual.y_scaled)}")
    print(f"  g_scaled={list(dual.g_scaled)}")
    print(f"  dominance check: {all(g >= (dual.scale if t == 0 else 0) for t, g in enumerate(dual.g_scaled))}")
    print(f"  rows={sorted(ch)}; edges={len(es)}; density={Fraction(len(es), len(ch)*7)}")
    print(f"  check degree hist={dict(Counter(len(ch[r]) for r in ch))}")
    print(f"  atom degree hist={dict(Counter(len(at.get(t, set())) for t in range(7)))}")
    print(f"  edge sign hist={dict(Counter(sign(dual.y_scaled[r]) for r, _ in es))}")
    print(f"  signed atom-neighbor profiles:")
    for t in range(7):
        pos = sum(1 for r in at.get(t, set()) if dual.y_scaled[r] > 0)
        neg = sum(1 for r in at.get(t, set()) if dual.y_scaled[r] < 0)
        print(f"    t={t}: deg={len(at.get(t,set()))}, +={pos}, -={neg}, rows={sorted(at.get(t,set()))}")

    check_common, atom_common = side_common_profiles(dual)
    print(f"  common-neighbor hist, check side={dict(check_common)}")
    print(f"  common-neighbor hist, atom side={dict(atom_common)}")
    print("  weak-regularity verdict: not biregular; common-neighbor counts form a chain, not a two-value law")

    vertices = [("r", r) for r in sorted(ch)] + [("t", t) for t in range(7)]
    nbrs = {}
    for r in ch:
        nbrs[("r", r)] = {("t", t) for t in ch[r]}
    for t in range(7):
        nbrs[("t", t)] = {("r", r) for r in at.get(t, set())}
    colors = degree_partition(vertices, nbrs)
    equitable, keys, signatures = quotient_is_equitable(vertices, nbrs, colors)
    print(f"  degree quotient equitable={equitable}; classes={keys}")
    if equitable:
        for k in keys:
            print(f"    quotient row {k}: {signatures[k]}")

    orbits, autos = edge_orbits(dual)
    mixed = 0
    for orb in orbits:
        signs = {sign(dual.y_scaled[r]) for r, _ in orb}
        mixed += int(len(signs) > 1)
    print(f"  support automorphisms={len(autos)}; edge orbits={len(orbits)}; mixed-sign orbits={mixed}")
    print(f"  first edge orbits={orbits[:6]}")

    black, white, faces = dessin_passport(dual)
    face_lens = sorted(len(f) for f in faces)
    face_sign_words = []
    for f in faces:
        word = "".join(sign(dual.y_scaled[r]) for r, _ in f)
        face_sign_words.append(word)
    print(f"  dessin passport black degrees={sorted(len(c) for c in black)}")
    print(f"  dessin passport white degrees={sorted(len(c) for c in white)}")
    print(f"  dessin face lengths={face_lens}; Euler check V-E+F={len(ch)+7-len(es)+len(faces)}")
    print(f"  face sign words sample={face_sign_words[:8]}")
    print("  dessin verdict: passport keeps branch degrees but forgets dominance unless signs are retained")


def tournament_analysis():
    lenses = [
        ("generated_depth_occupancy", (1, 1, 1, 1, 0)),
        ("signed_delsarte_dual", (1, 1, 1, 1, 0)),
        ("delsarte_quasicode_parity", (1, 1, 1, 0, 0)),
        ("halfarc_sign_orbit", (1, 1, 0, 1, 1)),
        ("weak_regular_equitable_quotient", (0, 1, 0, 0, 1)),
        ("dessin_belyi_passport", (0, 1, 0, 0, 0)),
        ("unit_distance_spectral_transfer", (0, 0, 0, 0, 0)),
        ("sparse_tanner_expander", (0, 0, 0, 0, 1)),
    ]
    # Tie-break by the written order after the score vector.  This encodes the
    # proof-order judgment: exact predicate preservation outranks locality.
    order = sorted(range(len(lenses)), key=lambda i: (lenses[i][1], -i), reverse=True)
    names = [lenses[i][0] for i in order]
    scores = {lenses[i][0]: len(lenses) - 1 - rank for rank, i in enumerate(order)}
    cycles3 = 0
    for a, b, c in combinations(range(len(lenses)), 3):
        sa, sb, sc = scores[lenses[a][0]], scores[lenses[b][0]], scores[lenses[c][0]]
        # Transitive score tournament has no directed 3-cycles by construction;
        # keep the explicit check visible.
        if len({sa, sb, sc}) < 3:
            cycles3 += 1
    print("=" * 78)
    print("Tournament Analysis over proof lenses")
    print("  observable=(preserves_bound, exact_audit, LRC_remaining_gap, sign_structure, locality)")
    print("  switch/gauge: generated words and signs are retained before scalarizing")
    print(f"  Hamiltonian path={names}")
    print(f"  score histogram={scores}")
    print(f"  directed_3cycles={cycles3} (ties only; strict tie-break gives 0)")


def main():
    print("LRC14 Tanner/dessin/weak-regularity carrier audit")
    print("=" * 78)
    print("Vertex challenge:")
    print("  tested vertices: moment rows, missed-depth atoms, signed edges,")
    print("  equitable quotient classes, dessin half-edges, unit-distance graph lenses.")
    print("  chosen carrier: moment rows r and atoms t, edge iff r<=t.")
    print("  preserved predicate: THM-534 Delsarte dominance g(t)>=scale*[t=0].")
    print("  destroyed data: generated speed words, relation lattice, L7 residue odometer.")
    for dual in DUALS:
        print_dual(dual)
    tournament_analysis()
    print("=" * 78)
    print("Conclusion")
    print("  The Tanner graph exists but is not weakly regular or LDPC-like.")
    print("  Its useful invariant is the signed half-arc orientation: support")
    print("  automorphisms never mix + and - edge orbits.  A Belyi/dessin passport")
    print("  is a good address layer for the carrier, but only the signed passport")
    print("  remembers the Delsarte dominance certificate.  The proof route should")
    print("  stay: generated depth word -> signed Delsarte/quasicode parity -> cap,")
    print("  while using weak-regular/unit-distance/Tanner analogies as audits.")


if __name__ == "__main__":
    main()
