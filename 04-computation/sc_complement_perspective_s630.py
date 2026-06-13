#!/usr/bin/env python3
"""
sc_complement_perspective_s630.py

Scout the mechanism behind self-complementary tournaments as node-perspective
swaps under global edge reversal, and compare the resulting carrier with the
unit-distance n=7/n=21 echoes.

The computation deliberately separates:

1. the anti-automorphism coset Anti(T) = Iso(T, T^op);
2. vertex anti-image sets {sigma(v): sigma in Anti(T)};
3. anti-automorphism cycle types from the Burnside self-converse formula;
4. raw scalar echoes 7 and 21 on the unit-distance side.

This is a small exact scout, not a large enumerator: SC classes are generated
through n=7 from one representative complementing permutation per allowed
cycle type, then canonicalized by brute force.
"""

from __future__ import annotations

from collections import Counter, deque
import contextlib
import io
import os
from itertools import permutations


KNOWN_SC = {1: 1, 2: 1, 3: 2, 4: 2, 5: 8, 6: 12, 7: 88}


def odd_partitions(n: int, max_part: int | None = None):
    if max_part is None:
        max_part = n
    if n == 0:
        yield ()
        return
    top = min(n, max_part)
    if top % 2 == 0:
        top -= 1
    for first in range(top, 0, -2):
        for rest in odd_partitions(n - first, first):
            yield (first,) + rest


def pairs_of(n: int):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def edge(bits: int, n: int, i: int, j: int) -> int:
    if i == j:
        raise ValueError("loop")
    if i > j:
        return 1 - edge(bits, n, j, i)
    idx = i * n - i * (i + 1) // 2 + (j - i - 1)
    return (bits >> idx) & 1


def bits_to_matrix(bits: int, n: int):
    adj = [[0] * n for _ in range(n)]
    for i, j in pairs_of(n):
        if edge(bits, n, i, j):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj


def permute_bits(bits: int, n: int, perm: tuple[int, ...]) -> int:
    out = 0
    for idx, (i, j) in enumerate(pairs_of(n)):
        val = edge(bits, n, perm[i], perm[j])
        if val:
            out |= 1 << idx
    return out


def complement_bits(bits: int, n: int) -> int:
    m = n * (n - 1) // 2
    return ((1 << m) - 1) ^ bits


def canonical_bits(bits: int, n: int) -> int:
    """Canonical form, restricting permutations to equal-score blocks."""
    scores = [sum(edge(bits, n, v, u) for u in range(n) if u != v) for v in range(n)]
    blocks = []
    for score in sorted(set(scores)):
        blocks.append([v for v in range(n) if scores[v] == score])

    best = None

    def extend(idx: int, prefix: list[int]):
        nonlocal best
        if idx == len(blocks):
            val = permute_bits(bits, n, tuple(prefix))
            if best is None or val < best:
                best = val
            return
        for p in permutations(blocks[idx]):
            extend(idx + 1, prefix + list(p))

    extend(0, [])
    assert best is not None
    return best


def is_automorphism(bits: int, n: int, perm: tuple[int, ...]) -> bool:
    for i, j in pairs_of(n):
        if edge(bits, n, perm[i], perm[j]) != edge(bits, n, i, j):
            return False
    return True


def is_anti_automorphism(bits: int, n: int, perm: tuple[int, ...]) -> bool:
    for i, j in pairs_of(n):
        if edge(bits, n, perm[i], perm[j]) != 1 - edge(bits, n, i, j):
            return False
    return True


def automorphisms(bits: int, n: int):
    return [p for p in permutations(range(n)) if is_automorphism(bits, n, p)]


def anti_automorphisms(bits: int, n: int):
    return [p for p in permutations(range(n)) if is_anti_automorphism(bits, n, p)]


def cycle_type(perm: tuple[int, ...]) -> tuple[int, ...]:
    used = set()
    out = []
    for i in range(len(perm)):
        if i in used:
            continue
        v = i
        length = 0
        while v not in used:
            used.add(v)
            length += 1
            v = perm[v]
        out.append(length)
    return tuple(sorted(out, reverse=True))


def anti_cycle_type_reps(n: int):
    """One complementing-permutation representative for each Burnside type."""
    if n <= 1:
        return [tuple(range(n))]
    reps = []
    m = n // 2
    for part in odd_partitions(m):
        lengths = [2 * x for x in part]
        if n % 2:
            lengths.append(1)
        perm = list(range(n))
        cursor = 0
        for length in lengths:
            block = list(range(cursor, cursor + length))
            if length == 1:
                perm[block[0]] = block[0]
            else:
                for a, b in zip(block, block[1:] + block[:1]):
                    perm[a] = b
            cursor += length
        reps.append(tuple(perm))
    return reps


def fixed_tournaments_for_anti_perm(n: int, sigma: tuple[int, ...]):
    """Enumerate tournaments T satisfying sigma: T -> T^op."""
    pairs = pairs_of(n)
    pair_index = {p: k for k, p in enumerate(pairs)}

    def transform(idx: int):
        i, j = pairs[idx]
        a, b = sigma[i], sigma[j]
        if a < b:
            jdx = pair_index[(a, b)]
            xor = 1
        else:
            jdx = pair_index[(b, a)]
            xor = 0
        return jdx, xor

    visited = set()
    orbits = []
    for start in range(len(pairs)):
        if start in visited:
            continue
        cur = start
        parity = 0
        seen = {}
        orbit = []
        ok = True
        while cur not in seen:
            seen[cur] = parity
            visited.add(cur)
            orbit.append((cur, parity))
            nxt, x = transform(cur)
            cur = nxt
            parity ^= x
        if parity != seen[cur]:
            ok = False
        if not ok:
            return []
        orbits.append(orbit)

    out = []
    for assignment in range(1 << len(orbits)):
        bits = 0
        for orbit_idx, orbit in enumerate(orbits):
            base = (assignment >> orbit_idx) & 1
            for idx, parity in orbit:
                val = base ^ parity
                if val:
                    bits |= 1 << idx
        out.append(bits)
    return out


def sc_classes(n: int):
    reps = {}
    for sigma in anti_cycle_type_reps(n):
        for bits in fixed_tournaments_for_anti_perm(n, sigma):
            can = canonical_bits(bits, n)
            reps.setdefault(can, bits)
    return sorted(reps)


def hamiltonian_path_count(bits: int, n: int) -> int:
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for v in range(n):
            val = dp[mask][v]
            if not val:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if edge(bits, n, v, u):
                    dp[mask | (1 << u)][u] += val
    return sum(dp[(1 << n) - 1])


def score_sequence(bits: int, n: int) -> tuple[int, ...]:
    return tuple(sorted(sum(edge(bits, n, i, j) for j in range(n) if j != i) for i in range(n)))


def c3_count(bits: int, n: int) -> int:
    total = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if (edge(bits, n, i, j) and edge(bits, n, j, k) and edge(bits, n, k, i)) or (
                    edge(bits, n, i, k) and edge(bits, n, k, j) and edge(bits, n, j, i)
                ):
                    total += 1
    return total


def is_strong(bits: int, n: int) -> bool:
    def reach(start: int, reverse: bool = False):
        seen = {start}
        q = deque([start])
        while q:
            v = q.popleft()
            for u in range(n):
                if u == v or u in seen:
                    continue
                ok = edge(bits, n, u, v) if reverse else edge(bits, n, v, u)
                if ok:
                    seen.add(u)
                    q.append(u)
        return len(seen) == n

    return reach(0, False) and reach(0, True)


def vertex_profile(bits: int, n: int, v: int):
    out = [u for u in range(n) if u != v and edge(bits, n, v, u)]
    inn = [u for u in range(n) if u != v and edge(bits, n, u, v)]
    scores = [sum(edge(bits, n, x, y) for y in range(n) if y != x) for x in range(n)]
    return {
        "score": len(out),
        "out_scores": tuple(sorted(scores[u] for u in out)),
        "in_scores": tuple(sorted(scores[u] for u in inn)),
    }


def orbit_partition(group, n: int):
    seen = set()
    parts = []
    for v in range(n):
        if v in seen:
            continue
        orb = {g[v] for g in group}
        changed = True
        while changed:
            changed = False
            for g in group:
                new = {g[x] for x in orb}
                if not new <= orb:
                    orb |= new
                    changed = True
        seen |= orb
        parts.append(tuple(sorted(orb)))
    return tuple(sorted(parts, key=lambda x: (len(x), x)))


def anti_image_sets(antis, n: int):
    return tuple(tuple(sorted({a[v] for a in antis})) for v in range(n))


def route_tournament_fingerprint():
    vertices = [
        "anti-coset torsor",
        "vertex perspective orbits",
        "Burnside cycle types",
        "Eisenstein Phi3 norm",
        "unit spine carrier",
        "tile/bulk side channel",
        "OCF H-gap guardrail",
        "raw scalar 7/21",
    ]
    gauges = [
        [
            "anti-coset torsor",
            "vertex perspective orbits",
            "Burnside cycle types",
            "OCF H-gap guardrail",
            "Eisenstein Phi3 norm",
            "unit spine carrier",
            "tile/bulk side channel",
            "raw scalar 7/21",
        ],
        [
            "unit spine carrier",
            "tile/bulk side channel",
            "Eisenstein Phi3 norm",
            "anti-coset torsor",
            "vertex perspective orbits",
            "OCF H-gap guardrail",
            "Burnside cycle types",
            "raw scalar 7/21",
        ],
        [
            "OCF H-gap guardrail",
            "anti-coset torsor",
            "Eisenstein Phi3 norm",
            "Burnside cycle types",
            "vertex perspective orbits",
            "unit spine carrier",
            "tile/bulk side channel",
            "raw scalar 7/21",
        ],
    ]
    rank = [{name: len(vertices) - i for i, name in enumerate(g)} for g in gauges]
    adj = [[0] * len(vertices) for _ in vertices]
    for i, a in enumerate(vertices):
        for j, b in enumerate(vertices):
            if i == j:
                continue
            wins = sum(r[a] > r[b] for r in rank)
            if wins >= 2:
                adj[i][j] = 1
    scores = Counter(sum(adj[i]) for i in range(len(vertices)))
    tri = 0
    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            for k in range(j + 1, len(vertices)):
                cyc = (
                    adj[i][j] and adj[j][k] and adj[k][i]
                ) or (
                    adj[i][k] and adj[k][j] and adj[j][i]
                )
                tri += int(cyc)
    hp = 0
    for p in permutations(range(len(vertices))):
        if all(adj[p[i]][p[i + 1]] for i in range(len(vertices) - 1)):
            hp += 1
    return vertices, scores, tri, hp


def centered_hex(k: int) -> int:
    return 3 * k * (k + 1) + 1


def unit_distance_carrier_rows():
    rows = [
        ("exact", 5, 7, "Phi3(2) edge echo; legal as 4 spine + 3 bulk"),
        ("exact", 7, 12, "first centered hexagon; canonical unit-flip order flops"),
        ("Harborth", 11, 21, "Phi3(2)*3 lattice echo; exact planar row is already 23"),
        ("exact", 11, 23, "small planar optimum after the lattice 21 echo"),
        ("triangular", 21, 47, "Eisenstein lattice carrier from S625"),
        ("Moser/exact", 21, 57, "THM-408 P_2^- row; not an H=21 scalar"),
        ("Moser lane", 22, 60, "known dense lane; exact 60/61 frontier remains open"),
    ]
    out = []
    for label, n, edges, note in rows:
        spine = n - 1
        bulk = edges - spine
        shell = None
        for k in range(1, 8):
            if bulk == centered_hex(k):
                shell = f"C_hex({k})"
        out.append((label, n, edges, spine, bulk, shell or "-", note))
    return out


def main():
    print("S630 SC COMPLEMENT PERSPECTIVE AND UNIT-DISTANCE CARRIER SCOUT")
    print("=" * 78)
    print()
    print("Core object:")
    print("  Anti(T) = {sigma in Sym(V): i->j in T iff sigma(j)->sigma(i) in T}.")
    print("  If T is self-complementary and tau in Anti(T), then Anti(T)=Aut(T)*tau.")
    print("  Programmatically, this is a complementing coset of the automorphism group.")
    print()

    print("Cyclotomic guardrail")
    print("--------------------")
    phi3_2 = 2 * 2 + 2 + 1
    print(f"Phi3(2)=2^2+2+1={phi3_2}")
    print(f"3*Phi3(2)={3 * phi3_2}")
    print(f"centered_hex(3)=3*3*(3+1)+1={centered_hex(3)}")
    print("Reading: Phi3 supplies a visible carrier scalar, not permission")
    print("to identify unit-distance edge counts with forbidden H evaluations.")
    print()

    print("SC class generation sanity check")
    print("-------------------------------")
    sc_by_n = {}
    for n in range(1, 8):
        classes = sc_classes(n)
        sc_by_n[n] = classes
        print(f"n={n}: generated SC classes={len(classes)} known={KNOWN_SC[n]}")
    print()

    print("H-gap check inside exact SC classes through n=7")
    print("----------------------------------------------")
    for n in range(1, 8):
        hs = [hamiltonian_path_count(bits, n) for bits in sc_by_n[n]]
        strong_hs = [hamiltonian_path_count(bits, n) for bits in sc_by_n[n] if is_strong(bits, n)]
        print(
            f"n={n}: H values={sorted(set(hs))[:14]}"
            f"{'...' if len(set(hs)) > 14 else ''} "
            f"contains7={7 in hs} contains21={21 in hs} "
            f"strong_minmax={((min(strong_hs), max(strong_hs)) if strong_hs else '-')}"
        )
    print()

    print("Complementing coset diagnostics")
    print("-------------------------------")
    for n in range(3, 8):
        type_incidence = Counter()
        aut_sizes = Counter()
        anti_sizes = Counter()
        rigid = 0
        score_reversal_fail = 0
        coset_size_fail = 0
        anti_image_equals_dual_aut_orbit = 0
        total_vertices = 0
        for bits in sc_by_n[n]:
            autos = automorphisms(bits, n)
            antis = anti_automorphisms(bits, n)
            aut_sizes[len(autos)] += 1
            anti_sizes[len(antis)] += 1
            rigid += int(len(autos) == 1)
            if len(autos) != len(antis):
                coset_size_fail += 1
            type_incidence.update(set(cycle_type(a) for a in antis))
            scores = [sum(edge(bits, n, v, u) for u in range(n) if u != v) for v in range(n)]
            aut_orbits = orbit_partition(autos, n)
            anti_sets = anti_image_sets(antis, n)
            for v in range(n):
                total_vertices += 1
                if any(scores[a[v]] != n - 1 - scores[v] for a in antis):
                    score_reversal_fail += 1
                dual = tuple(sorted(anti_sets[v]))
                if any(dual == orb for orb in aut_orbits):
                    anti_image_equals_dual_aut_orbit += 1
        print(f"n={n}: classes={len(sc_by_n[n])} rigid={rigid}")
        print(f"  |Aut| distribution:  {dict(sorted(aut_sizes.items()))}")
        print(f"  |Anti| distribution: {dict(sorted(anti_sizes.items()))}")
        print(f"  anti cycle-type incidence: {dict(sorted(type_incidence.items()))}")
        print(f"  coset size failures: {coset_size_fail}")
        print(f"  score-reversal vertex failures: {score_reversal_fail}")
        print(
            "  anti-image set is an Aut-orbit for "
            f"{anti_image_equals_dual_aut_orbit}/{total_vertices} vertices"
        )
    print()

    print("Two n=7 perspective examples")
    print("----------------------------")
    n = 7
    regular = 0
    for idx, (i, j) in enumerate(pairs_of(n)):
        if (j - i) % n in {1, 2, 3}:
            regular |= 1 << idx
    examples = [("regular cyclic R_7", canonical_bits(regular, n))]
    for bits in sc_by_n[n]:
        autos = automorphisms(bits, n)
        if len(autos) == 1:
            examples.append(("rigid SC example", bits))
            break
    for name, bits in examples:
        autos = automorphisms(bits, n)
        antis = anti_automorphisms(bits, n)
        print(f"{name}: H={hamiltonian_path_count(bits, n)} c3={c3_count(bits, n)}")
        print(f"  score={score_sequence(bits, n)} |Aut|={len(autos)} |Anti|={len(antis)}")
        print(f"  Aut orbits: {orbit_partition(autos, n)}")
        print(f"  anti cycle types: {dict(Counter(cycle_type(a) for a in antis))}")
        print(f"  anti-image sets: {anti_image_sets(antis, n)}")
        p0 = vertex_profile(bits, n, 0)
        a0 = antis[0][0]
        p1 = vertex_profile(bits, n, a0)
        print(f"  v=0 profile {p0}")
        print(f"  first anti sends 0->{a0}; image profile {p1}")
    print()

    print("Unit-distance carrier split table")
    print("---------------------------------")
    print("label       n   edges  spine  bulk  shell       note")
    for label, n, edges, spine, bulk, shell, note in unit_distance_carrier_rows():
        print(f"{label:10s} {n:2d} {edges:6d} {spine:6d} {bulk:5d} {shell:10s} {note}")
    print()

    print("Tournament Analysis over proof lenses")
    print("-------------------------------------")
    print("Vertices are proof/carrier lenses, not tournament vertices.")
    print("Pairwise observable: which lens preserves complement-perspective data,")
    print("unit-distance side channels, and the H-gap guardrail with less collapse.")
    print("Switches/gauges: SC-perspective, unit-distance carrier, H-gap guardrail.")
    print("Ties: insertion Hamiltonian path over the majority tournament.")
    vertices, score_hist, tri, hp = route_tournament_fingerprint()
    print(f"score_hist={dict(sorted(score_hist.items()))} directed_3cycles={tri} H={hp}")
    print("tie Hamiltonian path by majority score:")
    majority_scores = []
    # Recompute majority scores from the tournament by sorting the score histogram data
    # through the same helper's pairwise relation.
    gauges = [
        ["anti-coset torsor", "vertex perspective orbits", "Burnside cycle types",
         "OCF H-gap guardrail", "Eisenstein Phi3 norm", "unit spine carrier",
         "tile/bulk side channel", "raw scalar 7/21"],
        ["unit spine carrier", "tile/bulk side channel", "Eisenstein Phi3 norm",
         "anti-coset torsor", "vertex perspective orbits", "OCF H-gap guardrail",
         "Burnside cycle types", "raw scalar 7/21"],
        ["OCF H-gap guardrail", "anti-coset torsor", "Eisenstein Phi3 norm",
         "Burnside cycle types", "vertex perspective orbits", "unit spine carrier",
         "tile/bulk side channel", "raw scalar 7/21"],
    ]
    ranks = [{name: len(vertices) - i for i, name in enumerate(g)} for g in gauges]
    for v in vertices:
        score = 0
        for w in vertices:
            if v == w:
                continue
            if sum(r[v] > r[w] for r in ranks) >= 2:
                score += 1
        majority_scores.append((score, v))
    for score, v in sorted(majority_scores, reverse=True):
        print(f"  score={score}: {v}")
    print()

    print("Assumption Challenge")
    print("--------------------")
    print("Alternate vertices considered: raw tournament vertices, SC classes,")
    print("anti-automorphisms, automorphism orbits, anti-cycle types, unit-distance")
    print("points, unit edges, Hamiltonian paths, shell/bulk packets, and proof")
    print("obligations.  The chosen tournament vertices are proof/carrier lenses.")
    print("Preserved predicate: complementing-perspective duality plus the H-gap")
    print("guardrail under unit-distance carrier decompositions.")
    print("Destroyed data: full embeddings, all n>=8 SC classes, exact graph6")
    print("representatives for n=21, and continuous distance inequalities.")
    print("Challenged assumption: self-complementarity gives one canonical node swap;")
    print("the intrinsic object is the anti-coset, with a map chosen only after")
    print("coordinates are fixed.")
    print()

    print("Interpretation")
    print("--------------")
    print("- The SC node-swap object is not a single permutation; it is Anti(T),")
    print("  a coset/torsor for Aut(T). A chosen anti-map is a coordinate choice.")
    print("- Every tested anti-map reverses scores: score(sigma(v)) = n-1-score(v).")
    print("- The images of a vertex under all anti-maps are exactly the dual")
    print("  automorphism orbit, so node perspectives are swapped orbit-to-orbit.")
    print("- H=7 and H=21 remain absent in exact SC classes through n=7; unit-distance")
    print("  sees 7/21 only as decomposed edge-carrier echoes, not as tournament H.")
    print("- The n=21 Moser/exact row splits 57 = 20 spine + 37 bulk, and 37 is a")
    print("  centered-hex shell number. That is a recursive carrier signal, not an")
    print("  H=21 realization.")


def run_and_store():
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        main()
    text = buf.getvalue()
    print(text, end="")
    outpath = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..",
        "05-knowledge",
        "results",
        "sc_complement_perspective_s630.out",
    )
    with open(outpath, "w", encoding="utf-8") as f:
        f.write(text)


if __name__ == "__main__":
    run_and_store()
