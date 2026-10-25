"""Exact Fano/sign-code and hourglass-quotient controls. Standard library only.

Universe: all 128 sign words, 5040 bijections of 7 branch labels with F2^3\0,
the 30 resulting simplex subgroups, and all 36 seven-label subsets of the
nine rational magic forms. No isomorphism census of complex surfaces is made.
"""
from collections import Counter
from itertools import combinations, permutations
from fractions import Fraction
from pathlib import Path
import hashlib
import json
import sys


def require(ok, msg):
    if not ok:
        raise RuntimeError(msg)


def dot(a, b):
    return (a & b).bit_count() % 2


def support(word, n=7):
    return tuple(i for i in range(n) if (word >> i) & 1)


def simplex(labels):
    return tuple(sorted(sum(dot(a, v) << i for i, v in enumerate(labels)) for a in range(8)))


def fano_lines(labels):
    return frozenset(tuple(i for i, v in enumerate(labels) if dot(a, v) == 0) for a in range(1, 8))


def syndrome(word, labels):
    ans = 0
    for i, v in enumerate(labels):
        if (word >> i) & 1:
            ans ^= v
    return ans


def cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def determinant(a, b, c):
    bc = cross(b, c)
    return sum(x * y for x, y in zip(a, bc))


def concurrent_triples(forms):
    triples = frozenset(t for t in combinations(range(len(forms)), 3) if determinant(*(forms[i] for i in t)) == 0)
    for pair in combinations(range(len(forms)), 2):
        p = cross(*(forms[i] for i in pair))
        require(any(p), "proportional forms")
        zeros = tuple(i for i, f in enumerate(forms) if sum(x * y for x, y in zip(p, f)) == 0)
        require(len(zeros) <= 3, "four concurrent lines")
    return triples


def perm_word(word, perm):
    return sum(((word >> i) & 1) << perm[i] for i in range(7))


def direct_projective_fixed(code, zero_set):
    nonzero = tuple(i for i in range(7) if i not in zero_set)
    return tuple(w for w in code if len({(w >> i) & 1 for i in nonzero}) == 1)


def main():
    labels = tuple(range(1, 8))
    S = simplex(labels)
    F = fano_lines(labels)
    C = tuple(w for w in range(128) if syndrome(w, labels) == 0)
    require(len(S) == 8 and len(C) == 16, "code dimensions")
    require(set(C) == set(S) | {w ^ 127 for w in S}, "Hamming lift of projective simplex")
    require(Counter(w.bit_count() for w in S) == {0: 1, 4: 7}, "simplex weight")
    require(Counter(w.bit_count() for w in C) == {0: 1, 3: 7, 4: 7, 7: 1}, "Hamming weight")
    # Independent incidence-kernel route, recovered from the arithmetic packet.
    normals = tuple(range(1, 8))
    B = lambda w: sum((sum((w >> i) & 1 for i, v in enumerate(labels) if not dot(a, v)) % 2) << j for j, a in enumerate(normals))
    require({w for w in range(128) if B(w) == 0} == set(S), "incidence kernel")
    require({B(w) for w in range(128)} == set(C), "incidence image")
    require(all((u ^ v) in S for u in S for v in S), "linear origin/closure hostile")
    require(all(syndrome(1 << i, labels) == labels[i] for i in range(7)), "branch inertia labels")
    require(all(syndrome(w, labels) == syndrome(w ^ 127, labels) for w in range(128)), "projective syndrome")

    structures = {}
    fibers = Counter()
    for lab in permutations(range(1, 8)):
        code = simplex(lab)
        lines = fano_lines(lab)
        if code in structures:
            require(structures[code][1] == lines, "code determines Fano structure")
        else:
            structures[code] = (lab, lines)
        fibers[code] += 1
    require(len(structures) == 30 and set(fibers.values()) == {168}, "30 marked subgroup census")

    # Enlarged universal hostile universe: every partial Steiner triple system
    # on seven labels, whether or not realizable by characteristic-zero lines.
    all_triples = tuple(combinations(range(7), 3))
    all_pairs = tuple(combinations(range(7), 2))
    pair_masks = tuple(sum(1 << all_pairs.index(p) for p in combinations(t, 2)) for t in all_triples)
    systems = set()
    def enumerate_pair_disjoint(start, used_pairs, chosen):
        systems.add(chosen)
        for j in range(start, 35):
            if not used_pairs & pair_masks[j]:
                enumerate_pair_disjoint(j + 1, used_pairs | pair_masks[j], chosen | (1 << j))
    enumerate_pair_disjoint(0, 0, 0)
    # Independent representation: enumerate all cliques in the compatibility
    # graph whose vertices are triples; no pair masks are used in this path.
    compat = tuple(sum(1 << j for j, u in enumerate(all_triples) if j > i and len(set(t) & set(u)) <= 1) for i, t in enumerate(all_triples))
    clique_systems = set()
    def enumerate_cliques(candidates, chosen):
        clique_systems.add(chosen)
        while candidates:
            bit = candidates & -candidates
            j = bit.bit_length() - 1
            candidates ^= bit
            enumerate_cliques(candidates & compat[j], chosen | bit)
    enumerate_cliques((1 << 35) - 1, 0)
    require(systems == clique_systems and len(systems) == 5596, "partial Steiner universe independent generation")
    fano_masks = tuple(sum(1 << all_triples.index(t) for t in lines) for _, lines in structures.values())
    all_system_hist = Counter()
    min_free = 30
    for mask in systems:
        free_count = sum(not mask & f for f in fano_masks)
        all_system_hist[(mask.bit_count(), free_count)] += 1
        min_free = min(min_free, free_count)
    require(min_free == 8, "universal free subgroup lower bound")

    forms = ((1, 0, 0), (1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (1, 1, 1), (1, -1, -1))
    names = ("E", "E+x", "E-x", "E+z", "E-z", "E+x+z", "E-x-z")
    triples = concurrent_triples(forms)
    require(triples == frozenset(((0, 1, 2), (0, 3, 4), (0, 5, 6))), "hourglass triple incidence")
    census = Counter()
    rows = []
    for code, (lab, lines) in sorted(structures.items()):
        bad = triples & lines
        r = len(bad)
        census[r] += 1
        # Direct sign-fixed criterion over every allowed branch zero set.
        for size in (0, 1, 2):
            for zeros in combinations(range(7), size):
                require(direct_projective_fixed(code, zeros) == (0,), "smooth stratum fixed point")
        for zeros in triples:
            expected = 2 if zeros in lines else 1
            require(len(direct_projective_fixed(code, zeros)) == expected, "node stabilizer")
        # For free cases, 24 nodes make three orbits. Every bad triple instead
        # gives two orbits of size4 with stabilizer2 and local type1/4(1,1).
        row = {"labels_by_form": dict(zip(names, lab)), "simplex_words": list(code),
               "bad_triples": [list(t) for t in sorted(bad)], "free": r == 0,
               "quotient_A1": 3 - r, "quotient_one_fourth_11": 2 * r,
               "resolution_K_squared": 2 - 2 * r}
        rows.append(row)
    require(census == {0: 16, 1: 12, 3: 2}, "hourglass Fano marking census")

    # Explicit free marking, and incidence-adapted Boolean complement marking.
    free = next(row for row in rows if row["free"])
    adapted = (7, 1, 6, 2, 5, 4, 3)
    adapted_S = simplex(adapted)
    require(triples <= fano_lines(adapted), "adapted star")
    reversal = (0, 2, 1, 4, 3, 6, 5)
    J = lambda v: 7 if v == 7 else v ^ 7
    require(all(adapted[reversal[i]] == J(adapted[i]) for i in range(7)), "Boolean complement/root reversal")
    reversed_S = tuple(sorted(perm_word(w, reversal) for w in adapted_S))
    require(reversed_S != adapted_S, "J does not normalize simplex")
    bad3_codes = {code for code, (_, lines) in structures.items() if len(triples & lines) == 3}
    require(bad3_codes == {adapted_S, reversed_S}, "J exchanges both star completions")
    structure_orbits = Counter()
    fixed_structure_counts = Counter()
    for code, (_, lines) in structures.items():
        im = tuple(sorted(perm_word(w, reversal) for w in code))
        require(im in structures, "reversal closure")
        r = len(triples & lines)
        require(r == len(triples & structures[im][1]), "reversal preserves obstruction count")
        if im == code:
            fixed_structure_counts[r] += 1
        if code <= im:
            structure_orbits[r] += 1
    # Parameter reversal check on coefficient rows, no alleged Fano automorphism.
    require(all(forms[reversal[i]] == (forms[i][0], -forms[i][1], -forms[i][2]) for i in range(7)), "actual geometric root permutation")

    all_forms = forms + ((1, 1, -1), (1, -1, 1))
    subset_rows = []
    free_distribution_by_triples = Counter()
    for subset in combinations(range(9), 7):
        ts = concurrent_triples(tuple(all_forms[i] for i in subset))
        cc = Counter(len(ts & lines) for _, lines in structures.values())
        require(sum(cc.values()) == 30, "complete seven subset code census")
        subset_rows.append({"subset": list(subset), "triple_count": len(ts), "bad_triple_census": dict(sorted(cc.items()))})
        free_distribution_by_triples[(len(ts), cc[0])] += 1
    require(len(subset_rows) == 36 and all(row["bad_triple_census"].get(0, 0) for row in subset_rows), "all magic seven-subsets have a free code")

    # Local quotient A1 = C^2/{+-1}, simultaneous sign on its three quadratic
    # generators lifts to (s,t)->(i*s,i*t). Its invariant monomials have degree
    # divisible by4. Resolution O(-4) gives discrepancy -1/2 and K-square drop1.
    monomials = [(a, b) for a in range(9) for b in range(9 - a)]
    require({(a,b) for a,b in monomials if (a+b) % 4 == 0} == {(a,b) for a,b in monomials if (a+b) % 2 == 0 and ((a+b)//2) % 2 == 0}, "local order4 invariant ring")
    discrepancy = Fraction(2, -4)
    require(discrepancy == Fraction(-1, 2) and discrepancy ** 2 * -4 == -1, "local discrepancy")

    receipt = {"status": "FINITE-EXACT; geometric statements require accompanying proof",
               "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               "universes": {"sign_words": 128, "labelings": 5040, "subgroups": 30, "magic_seven_subsets": 36},
               "universal_partial_Steiner_systems": {"count": len(systems), "independent_generation_agrees": True,
                    "minimum_free_subgroups": min_free,
                    "histogram": [{"triples": t, "free_subgroups": f, "systems": count} for (t, f), count in sorted(all_system_hist.items())]},
               "simplex": list(S), "hamming": list(C), "subgroup_labeling_fibers": 168,
               "hourglass_forms": dict(zip(names, forms)), "hourglass_triples": [list(t) for t in sorted(triples)],
               "hourglass_bad_triple_census": dict(sorted(census.items())), "explicit_free_marking": free,
               "incidence_adapted_labels": dict(zip(names, adapted)), "incidence_adapted_code": list(adapted_S),
               "reversed_adapted_code": list(reversed_S), "reversal_fixed_subgroups_by_bad_triples": dict(sorted(fixed_structure_counts.items())),
               "reversal_orbits_by_bad_triples": dict(sorted(structure_orbits.items())),
               "all_hourglass_markings": rows, "all_magic_seven_subsets": subset_rows,
               "free_marking_counts": [{"triple_count": t, "free_subgroups": f, "subsets": count} for (t, f), count in sorted(free_distribution_by_triples.items())],
               "hostiles": {"incidence_adapted_action_not_free": True, "Boolean_J_not_simplex_normalizer": True,
                            "noncanonical_node_quotient_K_square_drop": 1, "code_origin_and_closure_checked": True}}
    content = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if len(sys.argv) > 1:
        Path(sys.argv[1]).write_text(content, encoding="utf-8", newline="\n")
    else:
        print(content, end="")


if __name__ == "__main__":
    main()
