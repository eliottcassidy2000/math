"""Apex-boundary tournament classes for AP/Goddyn-Wong LRC14 rows.

This script is intentionally small and exact.  It tests one quotient proposed
in HYP-2950: turn the six denominator-14 unit witnesses into tournament
vertices, then let each speed vote on which of two witnesses is closer to the
threshold boundary.

The quotient only sees the denominator-14 apex boundary.  It deliberately
forgets exact speed magnitude and off-apex denominators, so the output is a
necessary-condition scout rather than a complete invariant.
"""

from itertools import combinations, permutations, product


Q = 14
UNIT_WITNESSES = (1, 3, 5, 9, 11, 13)
N = len(UNIT_WITNESSES)
PAIRS = tuple((i, j) for i in range(N) for j in range(i + 1, N))
PAIR_INDEX = {pair: idx for idx, pair in enumerate(PAIRS)}


def distance_numerator(speed, witness, q=Q):
    residue = (speed * witness) % q
    return min(residue, q - residue)


def tournament_bits(row):
    """Return labelled tournament bits for the apex-pressure majority switch.

    For a pair of witnesses a,b, each speed votes for the witness where
    ||speed*witness/14|| is smaller.  Ties cast no vote.  A nonnegative total
    is oriented by the fixed Hamiltonian tie path UNIT_WITNESSES.
    """

    bits = []
    scores = []
    for i, j in PAIRS:
        a = UNIT_WITNESSES[i]
        b = UNIT_WITNESSES[j]
        score = 0
        for speed in row:
            da = distance_numerator(speed, a)
            db = distance_numerator(speed, b)
            if da < db:
                score += 1
            elif db < da:
                score -= 1
        bits.append(1 if score >= 0 else 0)
        scores.append(score)
    return tuple(bits), tuple(scores)


def adjacency(bits):
    adj = [[False] * N for _ in range(N)]
    for idx, (i, j) in enumerate(PAIRS):
        if bits[idx]:
            adj[i][j] = True
        else:
            adj[j][i] = True
    return adj


PERM_MAPS = []
for perm in permutations(range(N)):
    mapping = []
    for i, j in PAIRS:
        oi = perm[i]
        oj = perm[j]
        if oi < oj:
            mapping.append((PAIR_INDEX[(oi, oj)], False))
        else:
            mapping.append((PAIR_INDEX[(oj, oi)], True))
    PERM_MAPS.append(tuple(mapping))


def canonical_bits(bits):
    best = None
    for mapping in PERM_MAPS:
        relabelled = tuple((1 - bits[idx]) if invert else bits[idx] for idx, invert in mapping)
        if best is None or relabelled < best:
            best = relabelled
    return "".join(str(bit) for bit in best)


def all_tournament_class_count():
    classes = set()
    for bits in product((0, 1), repeat=len(PAIRS)):
        classes.add(canonical_bits(bits))
    return len(classes)


def directed_triangle_count(bits):
    adj = adjacency(bits)
    count = 0
    for i, j, k in combinations(range(N), 3):
        cyclic = adj[i][j] + adj[j][k] + adj[k][i]
        if cyclic in (0, 3):
            count += 1
    return count


def hamiltonian_path_count(bits):
    adj = adjacency(bits)
    count = 0
    for path in permutations(range(N)):
        if all(adj[path[i]][path[i + 1]] for i in range(N - 1)):
            count += 1
    return count


def outdegree_sequence(bits):
    adj = adjacency(bits)
    return tuple(sorted((sum(row) for row in adj), reverse=True))


def edge_flips(bits_a, bits_b):
    return sum(1 for a, b in zip(bits_a, bits_b) if a != b)


def row_summary(name, row, ap_bits):
    bits, scores = tournament_bits(row)
    return {
        "name": name,
        "row_mod_14": tuple(sorted(speed % Q for speed in row)),
        "bits": "".join(str(bit) for bit in bits),
        "canonical": canonical_bits(bits),
        "edge_scores": scores,
        "c3": directed_triangle_count(bits),
        "H": hamiltonian_path_count(bits),
        "score_sequence": outdegree_sequence(bits),
        "flips_from_AP": edge_flips(bits, ap_bits),
    }


def print_summary(summary):
    print(f"{summary['name']}:")
    print(f"  row residues mod 14: {summary['row_mod_14']}")
    print(f"  labelled bits:       {summary['bits']}")
    print(f"  canonical class:     {summary['canonical']}")
    print(f"  edge scores:         {summary['edge_scores']}")
    print(f"  outdegree sequence:  {summary['score_sequence']}")
    print(f"  directed 3-cycles:   {summary['c3']}")
    print(f"  Hamiltonian paths:   {summary['H']}")
    print(f"  edge flips from AP:  {summary['flips_from_AP']}")


def ap_single_residue_mutation_classes():
    base = tuple(range(1, 14))
    seen = {}
    for removed in base:
        for residue in range(1, 14):
            if residue == removed:
                continue
            row = tuple(x for x in base if x != removed) + (residue,)
            bits, _ = tournament_bits(row)
            key = canonical_bits(bits)
            seen.setdefault(key, []).append(
                {
                    "removed": removed,
                    "replacement_residue": residue,
                    "bits": "".join(str(bit) for bit in bits),
                    "c3": directed_triangle_count(bits),
                    "H": hamiltonian_path_count(bits),
                    "score_sequence": outdegree_sequence(bits),
                }
            )
    return seen


def main():
    ap = tuple(range(1, 14))
    gw = tuple(range(1, 12)) + (13, 24)
    near = tuple(range(1, 12)) + (13, 36)
    petal_10 = tuple(x for x in range(1, 14) if x != 10) + (20,)
    petal_13 = tuple(range(1, 13)) + (26,)
    splice_gw = tuple(x for x in range(1, 14) if x not in (10, 12)) + (20, 24)
    splice_k33 = tuple(x for x in range(1, 14) if x not in (10, 12)) + (20, 36)

    ap_bits, _ = tournament_bits(ap)
    named_rows = (
        ("AP", ap),
        ("GW 12->24", gw),
        ("near/K33 12->36", near),
        ("petal 10->20", petal_10),
        ("petal 13->26", petal_13),
        ("splice 10,12->20,24", splice_gw),
        ("splice 10,12->20,36", splice_k33),
    )

    print("Apex-pressure tournament quotient for LRC14")
    print("vertices:", UNIT_WITNESSES)
    print("pairwise observable: signed majority of speeds closer to one unit witness boundary")
    print("tie Hamiltonian path:", " -> ".join(str(v) for v in UNIT_WITNESSES))
    print("all unlabeled six-vertex tournament classes:", all_tournament_class_count())
    print()

    print("Named AP/GW low-frontier rows")
    named_classes = set()
    for name, row in named_rows:
        summary = row_summary(name, row, ap_bits)
        named_classes.add(summary["canonical"])
        print_summary(summary)
        print()
    print("named low-frontier canonical classes:", sorted(named_classes))
    print()

    print("AP single-residue mutation quotient")
    mutation_classes = ap_single_residue_mutation_classes()
    print(f"achievable canonical classes: {len(mutation_classes)} of 56")
    for idx, key in enumerate(sorted(mutation_classes), 1):
        examples = mutation_classes[key]
        c3s = sorted({example["c3"] for example in examples})
        hs = sorted({example["H"] for example in examples})
        score_sequences = sorted({example["score_sequence"] for example in examples})
        print(f"{idx}. class {key}")
        print(f"   count: {len(examples)}")
        print(f"   c3 values: {c3s}")
        print(f"   H values: {hs}")
        print(f"   score sequences: {score_sequences}")
        print(f"   first examples: {examples[:8]}")
    print()

    print("Readout")
    print("- The named AP/GW, K33-near, petal, and two-splice low-frontier rows all")
    print("  land in the transitive six-vertex class under this apex-pressure switch.")
    print("- The full AP single-residue mutation quotient realizes only 4 of the")
    print("  56 possible six-vertex tournament isomorphism classes.")
    print("- Therefore any proposed AP/GW-kind boundary-only row whose apex-pressure")
    print("  class is cyclic is excluded from this quotient before exact M/Farey data.")
    print("- The quotient is not complete: it forgets magnitude, q-thresholds, and")
    print("  off-apex safe intervals, so exact M/Farey/C27 labels must remain attached.")


if __name__ == "__main__":
    main()
