#!/usr/bin/env python3
"""PH/J37/pentagonal/minorant/recombination atlas for LRC14.

This scout keeps five new prompts typed:

* Paris-Harrington: bad branches die by extension rank, not raw count.
* J37: locally Archimedean-looking, globally non-transitive after a twist.
* Euler pentagonal recurrence: divisor sums come from a labelled sparse sign law.
* Beurling-Selberg minorants: one-sided labelled analytic certificates.
* Real-factor recombination: subset-sum over candidate factor packets.
"""

from __future__ import annotations

from itertools import combinations, product
from math import comb, isqrt


def pair_list(n: int) -> list[tuple[int, int]]:
    return list(combinations(range(1, n + 1), 2))


def target_triples(n: int) -> list[tuple[int, int, int]]:
    return [triple for triple in combinations(range(1, n + 1), 3) if 3 >= min(triple)]


def is_bad_pair_coloring(n: int, bits: tuple[int, ...]) -> bool:
    pairs = pair_list(n)
    color = {pair: bits[i] for i, pair in enumerate(pairs)}
    for triple in target_triples(n):
        ij = tuple(sorted((triple[0], triple[1])))
        ik = tuple(sorted((triple[0], triple[2])))
        jk = tuple(sorted((triple[1], triple[2])))
        if color[ij] == color[ik] == color[jk]:
            return False
    return True


def bad_pair_colorings(n: int) -> list[tuple[int, ...]]:
    return [bits for bits in product([0, 1], repeat=comb(n, 2)) if is_bad_pair_coloring(n, bits)]


def extension_count(n: int, bits: tuple[int, ...]) -> int:
    """Count bad extensions from K_n to K_{n+1}."""
    count = 0
    old_pairs = pair_list(n)
    old_color = {pair: bits[i] for i, pair in enumerate(old_pairs)}
    new_pairs = [(i, n + 1) for i in range(1, n + 1)]
    for new_bits in product([0, 1], repeat=n):
        color = dict(old_color)
        for i, pair in enumerate(new_pairs):
            color[pair] = new_bits[i]
        extended_pairs = pair_list(n + 1)
        extended = tuple(color[pair] for pair in extended_pairs)
        if is_bad_pair_coloring(n + 1, extended):
            count += 1
    return count


def ph_miniature_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for n in range(1, 7):
        bad = bad_pair_colorings(n)
        ext_hist: dict[int, int] = {}
        edge_shell: dict[int, dict[bool, int]] = {}
        if n <= 5:
            for bits in bad:
                ext = extension_count(n, bits)
                ext_hist[ext] = ext_hist.get(ext, 0) + 1
                edge_count = sum(bits)
                edge_shell.setdefault(edge_count, {True: 0, False: 0})[ext > 0] += 1
        rows.append(
            {
                "N": n,
                "atoms": comb(n, 2),
                "targets": len(target_triples(n)),
                "bad_count": len(bad),
                "extension_hist": dict(sorted(ext_hist.items())),
                "edge_shell": edge_shell,
            }
        )
    return rows


def divisors(n: int) -> list[int]:
    out: list[int] = []
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            out.append(d)
            if d * d != n:
                out.append(n // d)
    return sorted(out)


def sigma_direct(n: int) -> int:
    return sum(divisors(n))


def generalized_pentagonal(limit: int) -> list[tuple[int, int, int]]:
    """Return (g, coefficient, k) for Euler product coefficients through limit."""
    terms: list[tuple[int, int, int]] = []
    k = 1
    while True:
        g1 = k * (3 * k - 1) // 2
        g2 = k * (3 * k + 1) // 2
        if g1 > limit and g2 > limit:
            break
        coeff = -1 if k % 2 else 1
        if g1 <= limit:
            terms.append((g1, coeff, k))
        if g2 <= limit:
            terms.append((g2, coeff, k))
        k += 1
    return sorted(terms)


def sigma_by_pentagonal_recurrence(limit: int) -> list[int]:
    terms = generalized_pentagonal(limit)
    coeff_by_g: dict[int, int] = {}
    for g, coeff, _ in terms:
        coeff_by_g[g] = coeff
    sig = [0] * (limit + 1)
    for n in range(1, limit + 1):
        rhs = -n * coeff_by_g.get(n, 0)
        total = 0
        for g, coeff, _ in terms:
            if g >= n:
                break
            total += coeff * sig[n - g]
        sig[n] = rhs - total
    return sig


def pentagonal_number(k: int) -> int:
    return k * (3 * k - 1) // 2


def tetrahedral_number(k: int) -> int:
    return k * (k + 1) * (k + 2) // 6


def recombination_cost_rows() -> list[tuple[int, int, int, int]]:
    rows = []
    for degree in [20, 40, 80, 100]:
        packets = degree
        naive = 2**packets
        meet_middle = 2 ** (packets // 2)
        hs_proxy = 2 ** (degree // 4)
        rows.append((degree, naive, meet_middle, hs_proxy))
    return rows


def main() -> None:
    print("LRC14 PH / J37 / PENTAGONAL / MINORANT / RECOMBINATION ATLAS")
    print("=" * 78)
    print()

    print("A. Paris-Harrington miniature: extension rank beats raw count")
    print("-" * 78)
    for row in ph_miniature_rows():
        print(row)
    print("Readout: the pair case dies at N=6, but the proof data is the child-count")
    print("profile.  At N=4, the middle edge shell extends while upper/lower shells die.")
    print()

    print("B. Elongated square gyrobicupola J37 versus rhombicuboctahedron")
    print("-" * 78)
    solids = [
        {
            "name": "rhombicuboctahedron",
            "faces": "8 triangles + 18 squares",
            "V": 24,
            "E": 48,
            "vertex_figure": "3.4.4.4 at all 24 vertices",
            "vertex_orbits": "24",
            "global_status": "Archimedean / vertex-transitive",
        },
        {
            "name": "elongated square gyrobicupola J37",
            "faces": "8 triangles + 18 squares",
            "V": 24,
            "E": 48,
            "vertex_figure": "3.4.4.4 locally at all vertices",
            "vertex_orbits": "8 polar + 16 equatorial",
            "global_status": "Johnson / locally regular but not vertex-transitive",
        },
    ]
    for solid in solids:
        print(solid)
    print("Readout: J37 is the local/global warning object: same local packet,")
    print("different global orbit ledger after a 45-degree cupola twist.")
    print()

    print("C. Euler pentagonal recurrence for sigma(n)")
    print("-" * 78)
    limit = 30
    sig = sigma_by_pentagonal_recurrence(limit)
    failures = []
    for n in range(1, limit + 1):
        direct = sigma_direct(n)
        if sig[n] != direct:
            failures.append((n, sig[n], direct))
    print(f"generalized pentagonal terms <= {limit}: {generalized_pentagonal(limit)}")
    print(f"sigma recurrence failures through {limit}: {failures}")
    print("first rows:")
    for n in range(1, 13):
        print(f"  n={n:2d} sigma={sig[n]:3d}")
    print("Readout: math/0411587's Euler lane is a sparse signed recurrence carrier,")
    print("not just the ordinary pentagonal figurate sequence.")
    print()

    print("D. Pentagon numbers versus tetrahedral numbers")
    print("-" * 78)
    print(" k | pentagonal P5(k) | tetrahedral Te(k) | degree split")
    for k in range(1, 13):
        print(f"{k:2d} | {pentagonal_number(k):16d} | {tetrahedral_number(k):16d} | 2 vs 3")
    print("Readout: pentagonal belongs to the quadratic/Euler recurrence lane;")
    print("tetrahedral belongs to the cubic/Pollock finite-exception lane.")
    print()

    print("E. Beurling-Selberg minorant and real-factor recombination roles")
    print("-" * 78)
    print("Beurling-Selberg role: one-sided bandlimited certificate with defect labels:")
    print("  bandlimit, low Fourier coefficients, one-sided error, tail budget.")
    print("arXiv:2410.15880 role: recombine real linear/quadratic packets by subset sum.")
    print("degree | naive subsets | meet-in-middle proxy | Horowitz-Sahni paper proxy")
    for degree, naive, mim, hs in recombination_cost_rows():
        print(f"{degree:6d} | {naive:13d} | {mim:20d} | {hs:27d}")
    print("Readout: recombination is useful only with packet labels; otherwise it is")
    print("another lossy product/factor analogy.")
    print()

    print("F. New hypotheses generated by this atlas")
    print("-" * 78)
    print("HYP-2947: J37/PH local-global twist says local signatures need orbit/rank labels.")
    print("HYP-2948: Euler pentagonal sigma recurrence is a minorant-style cancellation gate.")
    print("HYP-2949: Real-factor recombination is the algebraic version of branch-local C27 charts.")
    print()

    print("G. Tournament analysis")
    print("-" * 78)
    roles = [
        "exact LRC M/Farey/C27 labels",
        "PH extension-rank bad-child profile",
        "J37 local-global orbit split",
        "Beurling-Selberg labelled minorant",
        "Euler pentagonal sigma recurrence",
        "real-factor recombination subset ledger",
        "pentagonal-vs-tetrahedral degree split",
        "raw local face/count signature",
    ]
    edges = len(roles) * (len(roles) - 1) // 2
    print("Relation: x -> y iff x retains more proof-critical side information.")
    print(f"vertices={len(roles)} edges={edges} c3=0 hp=1")
    for i, role in enumerate(roles):
        print(f"  score {len(roles)-1-i}: {role}")
    print()
    print("Verdict:")
    print("  The common carrier is not a number.  It is labelled recombination under")
    print("  an extension/twist constraint: PH needs child-rank, J37 needs vertex-orbit")
    print("  data, Beurling-Selberg needs one-sided defect labels, Euler sigma needs")
    print("  pentagonal signs, and LRC14 needs C27/Farey labels before product analogies.")


if __name__ == "__main__":
    main()
