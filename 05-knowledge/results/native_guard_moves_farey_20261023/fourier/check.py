"""Exact ternary Fourier compiler for authenticated Z/8 guard tables.

Run from any directory; no third-party dependencies.  Numerical roots of unity
are deliberately absent.  The native consumer is added to the receipt only
when the adjacent guard packet supplies its explicitly named fixture.
"""
from itertools import product
from pathlib import Path
from array import array
import csv
import gzip
import hashlib
import json
import sys

MOD = 8
ZERO, ONE = (0, 0), (1, 0)
ROOTS = (ONE, (0, 1), (7, 7))
CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


def add(x, y):
    return ((x[0] + y[0]) % MOD, (x[1] + y[1]) % MOD)


def mul(x, y):
    a, b = x
    c, d = y
    return ((a*c-b*d) % MOD, (a*d+b*c-b*d) % MOD)


def dot(x, y):
    return sum(a*b for a, b in zip(x, y)) % 3


def points(r):
    return list(product(range(3), repeat=r))


def transform(table, r):
    """Normalized DFT; table values may be ring elements or integer residues."""
    space = points(r)
    scale = (pow(3**r, -1, MOD), 0)
    result = {}
    for k in space:
        total = ZERO
        for a in space:
            value = table[a]
            if isinstance(value, int):
                value = (value % MOD, 0)
            total = add(total, mul(value, ROOTS[-dot(k, a) % 3]))
        result[k] = mul(scale, total)
    return result


def inverse(spectrum, r):
    return {
        a: ring_sum(mul(spectrum[k], ROOTS[dot(k, a)]) for k in points(r))
        for a in points(r)
    }


def ring_sum(values):
    result = ZERO
    for value in values:
        result = add(result, value)
    return result


def row_basis(rows, r):
    rows = [[x % 3 for x in row] for row in rows if any(x % 3 for x in row)]
    lead = 0
    for col in range(r):
        pivot = next((j for j in range(lead, len(rows)) if rows[j][col]), None)
        if pivot is None:
            continue
        rows[lead], rows[pivot] = rows[pivot], rows[lead]
        unit = pow(rows[lead][col], -1, 3)
        rows[lead] = [(unit*x) % 3 for x in rows[lead]]
        for j in range(len(rows)):
            if j != lead:
                c = rows[j][col]
                rows[j] = [(x-c*y) % 3 for x, y in zip(rows[j], rows[lead])]
        lead += 1
        if lead == len(rows):
            break
    return [tuple(row) for row in rows[:lead]]


def kernel_basis(equations, r):
    rows = row_basis(equations, r)
    pivots = [next(i for i, x in enumerate(row) if x) for row in rows]
    result = []
    for free in range(r):
        if free in pivots:
            continue
        v = [0]*r
        v[free] = 1
        for pivot, row in zip(pivots, rows):
            v[pivot] = -row[free] % 3
        result.append(tuple(v))
    return result


def shift(a, v):
    return tuple((x+y) % 3 for x, y in zip(a, v))


def weight(v):
    return sum(x != 0 for x in v)


def analyse(records, r, direct=True):
    """Each record is a tuple of mod-eight data; all original sidecars stay external."""
    space = points(r)
    require(set(records) == set(space), "guard table must include every source coordinate")
    arity = len(records[space[0]])
    require(all(len(records[a]) == arity for a in space), "record arity")
    require(all(isinstance(v, int) and 0 <= v < 8 for a in space for v in records[a]), "lossless residue components only")
    support = set()
    for j in range(arity):
        table = {a: records[a][j] for a in space}
        spectrum = transform(table, r)
        recovered = inverse(spectrum, r)
        require(all(recovered[a] == (table[a] % 8, 0) for a in space), "Fourier inversion")
        support.update(k for k in space if spectrum[k] != ZERO)
    equations = row_basis(sorted(support), r)
    basis = kernel_basis(equations, r)
    kernel = [v for v in space if all(dot(k, v) == 0 for k in equations)]
    if direct:
        literal = [v for v in space if all(records[shift(a, v)] == records[a] for a in space)]
        require(kernel == literal, "exact literal translation stabilizer")
    footprint = [i for i in range(r) if any(k[i] for k in support)]
    literal_footprint = [i for i in range(r) if any(
        records[a] != records[shift(a, tuple(int(j == i) for j in range(r)))] for a in space)]
    require(footprint == literal_footprint, "minimal coordinate footprint")
    quotient = {}
    for a in space:
        key = tuple(dot(k, a) for k in equations)
        require(key not in quotient or quotient[key] == records[a], "lossless guard quotient")
        quotient[key] = records[a]
    require(len(quotient) == 3**len(equations), "surjective linear observer")
    require(len(kernel) == 3**len(basis), "kernel lift multiplicity")
    budget = 0
    ranks = []
    for b in range(r+1):
        rank = len(row_basis([v for v in kernel if weight(v) <= b], r))
        ranks.append(rank)
        if rank == len(basis) and budget == 0 and basis:
            budget = b
    # Explicit Cayley connectivity is independent of the rank calculation.
    for b in range(r+1):
        moves = [v for v in kernel if weight(v) <= b]
        seen, todo = {tuple([0]*r)}, [tuple([0]*r)]
        while todo:
            a = todo.pop()
            for v in moves:
                c = shift(a, v)
                if c not in seen:
                    seen.add(c)
                    todo.append(c)
        require(len(seen) == 3**ranks[b], "low-budget component size")
    # A complete residual table (including zeros) is the exact suffix signature.
    residual_counts = []
    for depth in range(r+1):
        residuals = {tuple(records[p+s] for s in points(r-depth)) for p in points(depth)}
        residual_counts.append(len(residuals))
    return {"coordinates": r, "components": arity,
            "footprint": footprint, "support_rank": len(equations),
            "observer_rows": equations, "kernel_basis": basis,
            "kernel_lifts": len(kernel), "kernel_move_budget": budget,
            "low_budget_span_ranks": ranks, "residual_tables_by_depth": residual_counts}


def controls():
    require(mul((7, 1), (6, 7)) == (3, 0), "root differences have unit product three")
    require(mul(ROOTS[1], ROOTS[2]) == ONE, "cubic roots")
    # Exhaust every one-coordinate residue guard, and every binary two-coordinate guard.
    exhaustive = 0
    for r, alphabet in ((1, range(8)), (2, range(2))):
        space = points(r)
        for values in product(alphabet, repeat=len(space)):
            analyse({a: (v,) for a, v in zip(space, values)}, r)
            exhaustive += 1
    examples = {}
    space = points(3)
    examples["sum_residue"] = analyse({a: (sum(a) % 3,) for a in space}, 3)
    require(examples["sum_residue"]["kernel_move_budget"] == 2, "paired sum compensation")
    examples["two_differences"] = analyse({a: ((a[1]-a[0]) % 3, (a[2]-a[1]) % 3) for a in space}, 3)
    require(examples["two_differences"]["kernel_move_budget"] == 3, "common ternary phase needs all three coordinates")
    space = points(2)
    examples["zero_only"] = analyse({a: (int(a == (0, 0)),) for a in space}, 2)
    require(examples["zero_only"]["kernel_lifts"] == 1, "same-value fibres need not be kernel cosets")
    # Four cell phases, endpoints fixed to zero: valid differences are000 or111.
    binary_bridge = {}
    for a, b in space:
        phases = (0, a, b, 0)
        binary_bridge[(a, b)] = (int(all((y-x) % 3 in (0, 1) for x, y in zip(phases, phases[1:]))),)
    examples["pinned_binary_bridge"] = analyse(binary_bridge, 2)
    accepted = [a for a in space if binary_bridge[a] == (1,)]
    require(accepted == [(0, 0), (1, 2)], "000 to111 endpoint-pinned bridge")
    require(examples["pinned_binary_bridge"]["kernel_lifts"] == 1, "adaptive bridge is not a universal translation")
    # Arithmetic obstruction at modulus9: constant3 has zero unnormalized DFT.
    roots9 = ((1, 0), (0, 1), (8, 8))
    annihilated = [tuple(sum(3*roots9[-k*a % 3][j] for a in range(3)) % 9 for j in range(2)) for k in range(3)]
    require(annihilated == [(0, 0)]*3, "mod9 transform loses a nonzero constant")
    return {"exhaustive_guard_tables": exhaustive, "examples": examples,
            "mod9_nonzero_kernel_function": [3, 3, 3], "mod9_transform": annihilated}


def native_controls():
    packet = Path(__file__).resolve().parents[1]
    results = packet.parent
    guard_bytes = (packet / "guards/receipt.json").read_bytes()
    data = json.loads(guard_bytes)
    result = {"guard_receipt_sha256": hashlib.sha256(guard_bytes).hexdigest()}
    # Recover the original owner gauge and frozen potential independently.
    folder = results / "metagraph_order8_20261023"
    owner = array("H")
    owner.frombytes(gzip.decompress((folder / "owner_map.u16.gz").read_bytes()))
    if sys.byteorder != "little":
        owner.byteswap()
    with (folder / "merged_nodes.tsv").open() as handle:
        nodes = list(csv.DictReader(handle, delimiter="\t"))
    merged = {}
    for node in nodes:
        for key in ("class_a", "class_b"):
            merged[int(node[key])] = int(node["merged_id"])
    original = json.loads((results / "metagraph_order8_electrical_20261023/receipt.json").read_text())
    potential = original["order8"]["higher_layers"][-1]["lifts"][0]
    chords = [(i, j) for i in range(8) for j in range(i+2, 8)]
    full = (1 << len(chords))-1
    for name, fixture in (("fresh_character_n6", data["two_cyclic_cells"]),
                          ("original_character_pinned_n8", data["original_eight"]["pinned_padded_two_cells"])):
        A = fixture["adjacency"]
        records = {}
        target_records = {}
        for row in fixture["rows"]:
            a, b = row["a"], row["b"]
            path = row["path"]
            valid = int(all(A[u][v] for u, v in zip(path, path[1:])))
            require(valid == row["valid"] == int((b-a) % 3 in (0, 1)), "native junction replay")
            if len(A) == 8:
                require(path[0] == 0 and path[-1] == 7, "actual endpoint clock pins")
                if valid:
                    word = sum(A[path[i]][path[j]] << bit for bit, (i, j) in enumerate(chords))
                    source, target = merged[owner[word]], merged[owner[full ^ word]]
                    phase = (potential[source]-potential[target]) % 8
                    require((word, source, target, phase) == (row["normalized_forward_mask"], row["source_merged"], row["target_merged"], row["phase_mod8"]), "literal original phase and owner replay")
                    # Tournament converse reverses all arcs AND the observer order.
                    rev = path[::-1]
                    converse_word = sum(A[rev[j]][rev[i]] << bit for bit, (i, j) in enumerate(chords))
                    require(merged[owner[converse_word]] == source and merged[owner[full ^ converse_word]] == target, "converse preserves merged source/target, not a port-orientation reversal")
            records[a, b] = (valid, row["phase_mod8"] if valid else 0)
            target = row["target_merged"] if valid else 0
            require(0 <= target < 8**4, "exact target ID fits retained digit bank")
            digits = tuple((target // (8**j)) % 8 for j in range(4))
            require(sum(v*8**j for j, v in enumerate(digits)) == target, "target ID reconstructed exactly")
            target_records[a, b] = records[a, b]+digits
        summary = analyse(records, 2)
        full_summary = analyse(target_records, 2)
        require(summary["support_rank"] == full_summary["support_rank"] == 1, "one actual relative-phase observer")
        require(summary["kernel_basis"] == [(1, 1)] and summary["kernel_move_budget"] == 2, "native diagonal kernel and sharp coupled budget")
        values = [records[0, d] for d in range(3)]
        # Every indexed first-cell owner determines a; the remaining valid
        # second-cell choice is distinguished by the original phase.
        pinned = analyse({(b,): records[0, b] for b in range(3)}, 1)
        require(pinned["kernel_lifts"] == 1, "actual owner pin removes common-phase freedom")
        result[name] = {"actual_guard": summary, "guard_plus_exact_target_id": full_summary,
                        "relative_phase_records": values, "first_cell_pinned": pinned}
    return result


def main():
    result = {"scope": "exact finite guard tables, not a native realizability oracle", "controls": controls(), "native": native_controls()}
    result["checks"] = CHECKS
    output = Path(__file__).with_name("receipt.json")
    rendered = json.dumps(result, indent=2) + "\n"
    if "--write" in sys.argv:
        output.write_text(rendered, encoding="utf-8")
    else:
        require(json.loads(output.read_text(encoding="utf-8")) == json.loads(rendered), "frozen Fourier receipt")
    print(f"PASS {CHECKS} exact checks; receipt={output.name}")


if __name__ == "__main__":
    main()
