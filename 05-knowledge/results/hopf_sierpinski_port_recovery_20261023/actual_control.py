"""Reconstruct section obstructions from literal incoming order-31 block data."""
from collections import Counter
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import gzip
import json
import sys

HERE = Path(__file__).resolve().parent
CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def run():
    source = HERE.parent / "endpoint_overlap_ports_20261023/receipt.json.gz"
    raw = gzip.decompress(source.read_bytes()).replace(b"\r\n", b"\n")
    data = json.loads(raw)
    rows = data["matrix_outmasks"]
    n = len(rows)
    A = [[rows[u] >> v & 1 for v in range(n)] for u in range(n)]
    require(n == 31, "declared actual order")
    require(all(not A[u][u] for u in range(n)), "zero diagonal")
    require(all(A[u][v]+A[v][u] == 1 for u, v in combinations(range(n), 2)),
            "literal source is a tournament")
    profiles, total_loops, section_counts = Counter(), Counter(), Counter()
    witnesses = []
    for rec in data["cases"]:
        if "parallel_sections" not in rec:
            continue
        masks = rec.get("after_overlap_cells", rec["initial_cells"])
        cells = [[v for v in range(n) if mask >> v & 1] for mask in masks]
        require(sorted(v for C in cells for v in C) == list(range(n)), "actual partition")
        require(sorted(map(len, cells)) == [1]*7+[3]*8, "seven singletons eight triples")
        C = [c for c in cells if len(c) == 3]
        require(C[:3] == rec["parallel_sections"]["selected_boundary_cells"], "retained root and boundary")
        maps = {}
        for cell in C:
            require(all(sum(A[u][v] for v in cell) == 1 for u in cell), "directed triple")
        for i, j in combinations(range(8), 2):
            q = sum(A[C[i][0]][v] for v in C[j])
            require(q in (1, 2), "nonuniform regular block")
            require(all(sum(A[u][v] for v in C[j]) == q for u in C[i]), "common row sum")
            pi = {}
            for u in C[i]:
                targets = [v for v in C[j] if A[u][v] == (q == 1)]
                require(len(targets) == 1, "owned minority matching")
                pi[u] = targets[0]
            require(len(set(pi.values())) == 3, "selected matching bijection")
            maps[i, j] = pi
            maps[j, i] = {v: u for u, v in pi.items()}
        profile, common = Counter(), set(C[0])
        bad_cycle = None
        for i, j in combinations(range(1, 8), 2):
            hol = {x: maps[j, 0][maps[i, j][maps[0, i][x]]] for x in C[0]}
            fixed = {x for x in C[0] if hol[x] == x}
            common &= fixed
            profile[len(fixed)] += 1
            if not fixed and bad_cycle is None:
                bad_cycle = {"cells": [0, i, j], "root_permutation": sorted(hol.items())}
        direct = []
        for x in C[0]:
            assignment = [x]+[maps[0, j][x] for j in range(1, 8)]
            if all(maps[i, j][assignment[i]] == assignment[j] for i, j in combinations(range(8), 2)):
                direct.append(x)
        require(set(direct) == common, "direct all-edge and based-loop agreement")
        require(direct == rec["parallel_sections"]["root_fixed_vertices"], "incoming section record")
        require(not direct and bad_cycle is not None, "one actual triangle already obstructs sections")
        profiles[tuple(sorted(profile.items()))] += 1
        total_loops.update(profile)
        section_counts[len(direct)] += 1
        witnesses.append({"marks": rec["marks"], "class": rec["class"],
                          "root_star_profile": dict(profile), "one_obstructing_triangle": bad_cycle})
    require(profiles == Counter({((0, 9), (1, 12)): 168, ((0, 21),): 252}), "two exact profiles")
    require(section_counts == Counter({0: 420}), "420 empty native section sets")
    return {"status": "FINITE-EXACT", "input_receipt_sha256_lf": sha256(raw).hexdigest(),
            "scope": "420 declared actual singleton/triple partitions of the fixed order31 host",
            "not_claimed": "all order31 tournaments or a bounded-width decomposition of K8",
            "source_matrix_outmasks": rows, "section_counts": dict(section_counts),
            "based_loop_fixed_count_totals": dict(total_loops), "witnesses": witnesses,
            "checks": CHECKS}


if __name__ == "__main__":
    result = run()
    output = HERE / "actual_control.json"
    raw = json.dumps(result, indent=2, sort_keys=True)+"\n"
    if "--write" in sys.argv:
        output.write_text(raw, encoding="utf-8")
    else:
        require(json.loads(output.read_text(encoding="utf-8")) == json.loads(raw), "frozen actual control")
    print(json.dumps({"checks": CHECKS, "sections": result["section_counts"],
                      "loop_counts": result["based_loop_fixed_count_totals"]}))
