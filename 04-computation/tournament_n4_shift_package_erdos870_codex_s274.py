"""HYP-3146 scout: n=4 tournament shift packages and Erdos-870 filler/canary logic.

This script compares two presentations of the four unlabeled tournaments on
four vertices:

1. fixed Hamiltonian path / half-tiling presentation with three free chords;
2. finite-filler scaffold presentation with partial score sequence 0,1,1,2 and
   two remaining endpoint arcs.

The point is not the n=4 census itself.  The point is the controlled-forgetting
lesson: a raw fixed-path cube is a cover with a large S fiber, while a finite
scaffold can turn the same four class labels into a two-bit shift package.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations, product


VERTICES = tuple(range(4))
EDGES = tuple(combinations(VERTICES, 2))

CLASS_BY_SCORE = {
    (0, 1, 2, 3): "T",
    (1, 1, 2, 2): "S",
    (0, 2, 2, 2): "+",
    (1, 1, 1, 3): "-",
}


def edge_bit_to_arrow(edge: tuple[int, int], bit: int) -> str:
    i, j = edge
    return f"{i}->{j}" if bit == 0 else f"{j}->{i}"


def scores(bits: tuple[int, ...]) -> tuple[int, ...]:
    out = [0] * 4
    for bit, (i, j) in zip(bits, EDGES):
        if bit == 0:
            out[i] += 1
        else:
            out[j] += 1
    return tuple(sorted(out))


def label(bits: tuple[int, ...]) -> str:
    return CLASS_BY_SCORE[scores(bits)]


def bits_from_assignments(assignments: dict[tuple[int, int], int]) -> tuple[int, ...]:
    return tuple(assignments[e] for e in EDGES)


def flip(bits: tuple[int, ...], idxs: tuple[int, ...] | list[int]) -> tuple[int, ...]:
    out = list(bits)
    for idx in idxs:
        out[idx] ^= 1
    return tuple(out)


def fixed_path_model() -> dict[str, object]:
    """Fixed Hamiltonian path 0->1->2->3, free chords a,b,c.

    The default state is transitive.  The fixed path edges are
    (0,1),(1,2),(2,3).  The free chords are named so that single flips match
    the user's convention:

    a=(0,2) gives +, b=(1,3) gives -, c=(0,3) gives S.
    """

    base = (0, 0, 0, 0, 0, 0)
    free = {
        "a": EDGES.index((0, 2)),
        "b": EDGES.index((1, 3)),
        "c": EDGES.index((0, 3)),
    }

    states = {}
    for mask in range(8):
        names = tuple(name for bit, name in enumerate(("a", "b", "c")) if mask & (1 << bit))
        idxs = [free[name] for name in names]
        states[names or ("E",)] = flip(base, idxs)

    return {"base": base, "free": free, "states": states}


def restricted_generator_table(model: dict[str, object]) -> list[list[str]]:
    names = ("E", "a", "b", "c")
    free = model["free"]
    base = model["base"]
    rows = []
    for r in names:
        row = []
        for c in names:
            idxs = []
            if r != "E":
                idxs.append(free[r])
            if c != "E":
                idxs.append(free[c])
            row.append(label(flip(base, idxs)))
        rows.append(row)
    return rows


def full_fixed_path_fibers(model: dict[str, object]) -> dict[str, list[tuple[str, ...]]]:
    fibers: dict[str, list[tuple[str, ...]]] = defaultdict(list)
    for name_tuple, bits in model["states"].items():
        pretty = tuple() if name_tuple == ("E",) else name_tuple
        fibers[label(bits)].append(pretty)
    return dict(fibers)


def class_pgf_for_fixed_path(model: dict[str, object]) -> dict[str, Counter[int]]:
    pgf: dict[str, Counter[int]] = defaultdict(Counter)
    for name_tuple, bits in model["states"].items():
        weight = 0 if name_tuple == ("E",) else len(name_tuple)
        pgf[label(bits)][weight] += 1
    return dict(pgf)


def deletion_profile_for_s_fiber(model: dict[str, object]) -> list[tuple[tuple[str, ...], list[str], bool]]:
    free = model["free"]
    base = model["base"]
    out = []
    for reps in full_fixed_path_fibers(model)["S"]:
        if not reps:
            continue
        child_labels = []
        for remove in reps:
            kept = [name for name in reps if name != remove]
            bits = flip(base, [free[name] for name in kept])
            child_labels.append(label(bits))
        out.append((reps, child_labels, all(x == "S" for x in child_labels)))
    return out


def find_two_bit_scaffold() -> dict[str, object]:
    """Find a two-variable scaffold satisfying the prompt's partial-score condition."""

    target_orders = (("T", "-", "+", "S"), ("T", "+", "-", "S"))
    for var_idxs in combinations(range(len(EDGES)), 2):
        fixed = [i for i in range(len(EDGES)) if i not in var_idxs]
        for fixed_vals in product((0, 1), repeat=4):
            bits = [None] * len(EDGES)
            partial = [0] * 4
            for idx, val in zip(fixed, fixed_vals):
                bits[idx] = val
                i, j = EDGES[idx]
                if val == 0:
                    partial[i] += 1
                else:
                    partial[j] += 1
            if tuple(sorted(partial)) != (0, 1, 1, 2):
                continue
            labels = []
            for x_bit, y_bit in product((0, 1), repeat=2):
                trial = list(bits)
                trial[var_idxs[0]] = x_bit
                trial[var_idxs[1]] = y_bit
                labels.append(label(tuple(trial)))
            if tuple(labels) in target_orders:
                # Normalize variable names so x gives + and y gives -.
                if tuple(labels) == ("T", "-", "+", "S"):
                    x_idx, y_idx = var_idxs[0], var_idxs[1]
                else:
                    y_idx, x_idx = var_idxs[0], var_idxs[1]
                return {
                    "fixed_bits": tuple(bits),
                    "partial_scores": tuple(partial),
                    "x_idx": x_idx,
                    "y_idx": y_idx,
                }
    raise RuntimeError("no scaffold found")


def scaffold_state(scaffold: dict[str, object], x: int, y: int) -> tuple[int, ...]:
    bits = list(scaffold["fixed_bits"])
    bits[scaffold["x_idx"]] = x
    bits[scaffold["y_idx"]] = y
    return tuple(bits)


def scaffold_table(scaffold: dict[str, object], restricted: bool) -> list[tuple[str, list[str]]]:
    names = ("E", "x", "y") if restricted else ("E", "x", "y", "xy")
    coords = {"E": (0, 0), "x": (1, 0), "y": (0, 1), "xy": (1, 1)}
    rows = []
    for r in names:
        row = []
        for c in names:
            x = coords[r][0] ^ coords[c][0]
            y = coords[r][1] ^ coords[c][1]
            row.append(label(scaffold_state(scaffold, x, y)))
        rows.append((r, row))
    return rows


def monotone_or_compression() -> list[tuple[tuple[int, int, int], tuple[int, int], str, str]]:
    """Class-preserving compression from fixed-path bits (a,b,c) to scaffold bits.

    x = a OR c, y = b OR c.  The c chord is the clustered canary: activating it
    forces both endpoint variables live.
    """

    model = fixed_path_model()
    free = model["free"]
    base = model["base"]
    out = []
    for a, b, c in product((0, 1), repeat=3):
        bits = flip(base, [idx for flag, idx in zip((a, b, c), (free["a"], free["b"], free["c"])) if flag])
        x = int(bool(a or c))
        y = int(bool(b or c))
        scaffold_class = {  # T,+,-,S under x,y convention.
            (0, 0): "T",
            (1, 0): "+",
            (0, 1): "-",
            (1, 1): "S",
        }[(x, y)]
        out.append(((a, b, c), (x, y), label(bits), scaffold_class))
    return out


def tournament_analysis() -> dict[str, object]:
    carriers = [
        {
            "name": "finite_filler_scaffold_shift_package",
            "predicate": 5,
            "exact": 5,
            "quotient": 5,
            "fiber_memory": 3,
            "deletion": 2,
            "lrc_transfer": 5,
            "cost": 2,
        },
        {
            "name": "clustered_canary_S_fiber",
            "predicate": 4,
            "exact": 4,
            "quotient": 2,
            "fiber_memory": 5,
            "deletion": 5,
            "lrc_transfer": 5,
            "cost": 3,
        },
        {
            "name": "monotone_OR_shift_package",
            "predicate": 4,
            "exact": 4,
            "quotient": 4,
            "fiber_memory": 4,
            "deletion": 3,
            "lrc_transfer": 5,
            "cost": 2,
        },
        {
            "name": "fixed_path_half_tiling_cube",
            "predicate": 4,
            "exact": 5,
            "quotient": 1,
            "fiber_memory": 5,
            "deletion": 4,
            "lrc_transfer": 4,
            "cost": 4,
        },
        {
            "name": "edge_tip_tail_witness_packet",
            "predicate": 5,
            "exact": 4,
            "quotient": 4,
            "fiber_memory": 5,
            "deletion": 4,
            "lrc_transfer": 5,
            "cost": 4,
        },
        {
            "name": "fiber_PGF_moment_packet",
            "predicate": 5,
            "exact": 5,
            "quotient": 3,
            "fiber_memory": 4,
            "deletion": 2,
            "lrc_transfer": 5,
            "cost": 3,
        },
        {
            "name": "raw_score_sequence_class",
            "predicate": 2,
            "exact": 3,
            "quotient": 3,
            "fiber_memory": 1,
            "deletion": 1,
            "lrc_transfer": 1,
            "cost": 1,
        },
        {
            "name": "raw_fixed_path_class_count",
            "predicate": 1,
            "exact": 2,
            "quotient": 1,
            "fiber_memory": 1,
            "deletion": 1,
            "lrc_transfer": 1,
            "cost": 1,
        },
    ]
    axes = ("predicate", "exact", "quotient", "fiber_memory", "deletion", "lrc_transfer")
    path = [c["name"] for c in carriers]
    scores_out = Counter()
    edges = []
    flips_against_raw_count = 0
    for i, a in enumerate(carriers):
        wins = 0
        for j, b in enumerate(carriers):
            if i == j:
                continue
            score = sum((a[x] > b[x]) - (a[x] < b[x]) for x in axes)
            if score == 0:
                # Lower proof cost wins, then listed Hamiltonian path.
                if a["cost"] != b["cost"]:
                    a_wins = a["cost"] < b["cost"]
                else:
                    a_wins = path.index(a["name"]) < path.index(b["name"])
            else:
                a_wins = score > 0
            if a_wins:
                wins += 1
                edges.append((a["name"], b["name"]))
                if path.index(a["name"]) > path.index(b["name"]):
                    flips_against_raw_count += 1
        scores_out[wins] += 1

    directed_3cycles = 0
    edge_set = set(edges)
    for a, b, c in combinations([x["name"] for x in carriers], 3):
        if (
            ((a, b) in edge_set and (b, c) in edge_set and (c, a) in edge_set)
            or ((a, c) in edge_set and (c, b) in edge_set and (b, a) in edge_set)
        ):
            directed_3cycles += 1

    selected_path = sorted(carriers, key=lambda c: (-sum(c[x] for x in axes), c["cost"], path.index(c["name"])))
    return {
        "score_hist": dict(sorted(scores_out.items())),
        "directed_3cycles": directed_3cycles,
        "edge_flips_against_listed_path": flips_against_raw_count,
        "selected_path": [c["name"] for c in selected_path],
    }


def print_table(headers: tuple[str, ...], rows: list[tuple[str, list[str]]]) -> None:
    print("        " + " ".join(f"{h:>3}" for h in headers))
    for name, row in rows:
        print(f"{name:>5} | " + " ".join(f"{x:>3}" for x in row))


def main() -> None:
    print("HYP-3146 / codex-2026-06-27-S274")
    print("n=4 tournament shift-package scout inspired by Erdos-870 filler/canary logic")
    print()
    print("External trigger: davidturturean/erdos-870 negative answer to Erdos Problem #870")
    print("Transferable proof pattern: sparse core + finite filler gadget + shift package + canary cluster + deletion-stable nonminimality")
    print()

    model = fixed_path_model()
    print("FIXED-PATH / HALF-TILING MODEL")
    print("base path: 0->1->2->3; transitive base class T")
    print("free chords: a=(0,2), b=(1,3), c=(0,3)")
    print("restricted generator table matching the prompt:")
    print_table(("E", "a", "b", "c"), [(n, r) for n, r in zip(("E", "a", "b", "c"), restricted_generator_table(model))])
    print()
    print("full F2^3 fixed-path fibers:")
    fibers = full_fixed_path_fibers(model)
    for class_name in ("T", "+", "-", "S"):
        print(f"  {class_name}: {fibers[class_name]}")
    print(f"fiber sizes: {dict((k, len(v)) for k, v in fibers.items())}")
    print("class PGFs by flip weight:")
    for class_name in ("T", "+", "-", "S"):
        pgf = class_pgf_for_fixed_path(model)[class_name]
        terms = " + ".join(f"{coef}*z^{deg}" for deg, coef in sorted(pgf.items()))
        print(f"  {class_name}: {terms}")
    print("S-fiber deletion robustness:")
    for reps, child_labels, stable in deletion_profile_for_s_fiber(model):
        print(f"  {reps}: delete-one children={child_labels}; all_S={stable}")
    print()

    scaffold = find_two_bit_scaffold()
    print("FINITE-FILLER / TWO-BIT SHIFT-PACKAGE MODEL")
    print("fixed scaffold arrows:")
    for idx, bit in enumerate(scaffold["fixed_bits"]):
        if idx not in (scaffold["x_idx"], scaffold["y_idx"]):
            print(f"  {EDGES[idx]}: {edge_bit_to_arrow(EDGES[idx], bit)}")
    print(f"partial outscore vector={scaffold['partial_scores']}; sorted={tuple(sorted(scaffold['partial_scores']))}")
    print(f"x variable edge={EDGES[scaffold['x_idx']]}; y variable edge={EDGES[scaffold['y_idx']]}")
    print("restricted table from the prompt:")
    print_table(("E", "x", "y"), scaffold_table(scaffold, restricted=True))
    print("full Klein-four table after adding xy=S:")
    print_table(("E", "x", "y", "xy"), scaffold_table(scaffold, restricted=False))
    print()

    print("CLASS-PRESERVING COMPRESSION FROM MODEL A TO MODEL B")
    print("x = a OR c, y = b OR c.  The c chord is the clustered canary: it activates both endpoint variables.")
    for abc, xy, old, new in monotone_or_compression():
        print(f"  abc={abc} -> xy={xy}; fixed-path={old}; scaffold={new}; ok={old == new}")
    print()

    print("INTERPRETATION")
    print("1. Model A is a fixed-path cover: F2^3 -> {T,+,-,S} with S fiber size 5.")
    print("2. Model B is a finite-filler section: fixing one more arc makes {T,+,-,S} a two-bit shift package.")
    print("3. The user tables are not both the same algebra: Model A is a generator-local absorbing-ideal table, not a congruence; Model B becomes the Klein four group once xy=S is named.")
    print("4. Erdos-870 suggests a proof move for LRC packets: add finite filler/scaffold data until the quotient becomes a shift package, but keep canary-cluster redundancy when deletion stability is the predicate.")
    print("5. Generating-function signal: Model A has S_PGF=z+3z^2+z^3, while Model B has S_PGF=z^2.  The whole difference is hidden fiber mass.")
    print()

    ta = tournament_analysis()
    print("TOURNAMENT ANALYSIS OVER PROOF-CARRIER INTERPRETATIONS")
    print("vertices=proof carriers, not runners or raw arcs")
    print("pairwise observable=majority over predicate retention, exactness, quotient legality, fiber memory, deletion stability, and LRC transfer")
    print("switch=orient A->B when A wins more axes; ties prefer lower proof cost, then the listed guardrail path")
    print(f"score_hist={ta['score_hist']}")
    print(f"directed_3cycles={ta['directed_3cycles']}")
    print(f"edge_flips_against_listed_path={ta['edge_flips_against_listed_path']}")
    print("selected_hamiltonian_path=" + " -> ".join(ta["selected_path"]))


if __name__ == "__main__":
    main()
