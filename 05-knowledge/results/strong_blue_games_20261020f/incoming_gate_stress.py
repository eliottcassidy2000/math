"""A sharp nontrivial pinned-stress witness from the incoming spindle gate.

Imports the already audited exact coordinates and complete contacts from one
fixed Git commit. Does not run or import its producer and does not repeat its
geometric census. Independently checks local unit rows, a full-rank pinned
minor, and every core-vertex deletion under the fixed boundary colouring.
All guards remain active under python -O.
"""
from fractions import Fraction as F
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
import subprocess


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
INCOMING = "411c1fd35a"
SOURCE = "05-knowledge/results/kernels_games_detours_20261020f/forcing_check.json"
ZERO = (F(0),) * 4
ONE = (F(1), F(0), F(0), F(0))


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT)


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def mul(a, b):
    # Basis index bits represent sqrt(3) and sqrt(11).
    out = [F(0)] * 4
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            common = i & j
            factor = (3 if common & 1 else 1) * (11 if common & 2 else 1)
            out[i ^ j] += factor * x * y
    return tuple(out)


def vector(points, v, w):
    return tuple(sub(points[v][i], points[w][i]) for i in (0, 1))


def determinant(u, v):
    return sub(mul(u[0], v[1]), mul(u[1], v[0]))


def valid_core_colours(vertices, edges, lists):
    for colours in product(*(lists[v] for v in vertices)):
        assignment = dict(zip(vertices, colours))
        if all(assignment[v] != assignment[w] for v, w in edges
               if v in assignment and w in assignment):
            yield assignment


def main():
    commit = git("rev-parse", INCOMING).decode().strip()
    raw = git("show", f"{commit}:{SOURCE}")
    data = json.loads(raw)
    need(data["status"] == "PASS", "incoming contact receipt status")
    need(data["field_basis"] == ["1", "sqrt3", "sqrt11", "sqrt33"], "field basis")
    points = [tuple(tuple(F(c) for c in coord) for coord in p) for p in data["points"]]
    need(len(points) == len(set(points)) == 19, "distinct imported points")
    edges = [tuple(e) for e in data["graph19_edges"]]
    need(len(edges) == len(set(edges)) == 34, "full contact count")
    core = tuple(range(7))
    boundary = {int(v): c for v, c in data["boundary_colouring"].items()}
    need(set(boundary) == set(range(7, 19)), "twelve boundary vertices")
    core_edges = [(v, w) for v, w in edges if w < 7]
    port_edges = [(v, w) for v, w in edges if v < 7 <= w]
    boundary_edges = [(v, w) for v, w in edges if v >= 7]
    need(core_edges == [tuple(e) for e in data["spindle_edges"]], "spindle contact agreement")
    need(len(core_edges) == 11 and len(port_edges) == 21 and len(boundary_edges) == 2,
         "edge partition")
    need(all(boundary[v] != boundary[w] for v, w in boundary_edges), "proper fixed boundary")
    lists = {}
    minor = []
    full_determinant = ONE
    a = t = 0
    for v in core:
        ports = [w for u, w in port_edges if u == v]
        forbidden = {boundary[w] for w in ports}
        need(len(ports) == 3 and forbidden == {3, 4, 5}, "three distinct forbidden port colours")
        lists[v] = sorted(set(range(6)) - forbidden)
        a += len(ports)
        t += len(forbidden)
        # Recheck only the 21 imported core-port contacts, not the full census.
        rows = {w: vector(points, v, w) for w in ports}
        need(all(add(mul(r[0], r[0]), mul(r[1], r[1])) == ONE for r in rows.values()),
             "imported pinned edge has unit length")
        selected = next(((w, u, determinant(rows[w], rows[u]))
                         for w, u in combinations(ports, 2)
                         if determinant(rows[w], rows[u]) != ZERO), None)
        need(selected is not None, "a nonparallel pair at each core vertex")
        w, u, det = selected
        full_determinant = mul(full_determinant, det)
        minor.append({"core_vertex": v, "ports": [w, u],
                      "two_by_two_determinant": [str(x) for x in det]})
    need(full_determinant != ZERO, "nonzero 14-row block-diagonal determinant")
    # All chosen rows meet exactly one core vertex, so their 14x14 minor is
    # block diagonal with the seven determinants above. Fourteen is also
    # the number of columns, certifying the exact pinned rank.
    rank = 14
    h = len(core)
    row_count = len(core_edges) + a
    stress = row_count - rank
    sharp_bound = h + (a + 1) // 2
    need(a == t == 21 and row_count == 32 and stress == sharp_bound == 18,
         "sharp six-colour stress bound")
    need(not list(valid_core_colours(core, core_edges, lists)), "blocked full core")
    deletion_checks = []
    for deleted in core:
        vertices = tuple(v for v in core if v != deleted)
        colourings = list(valid_core_colours(vertices, core_edges, lists))
        need(bool(colourings), "every core deletion restores extension")
        witness = colourings[0]
        total = dict(boundary)
        total.update(witness)
        need(all(total[v] != total[w] for v, w in edges if deleted not in (v, w)),
             "full fixed-boundary deletion witness")
        deletion_checks.append({"deleted_core_vertex": deleted,
                                "extension_count": len(colourings), "core_witness": witness})
    witness4 = {int(v): c for v, c in data["global_four_colouring"].items()}
    need(set(witness4) == set(range(19)) and len(set(witness4.values())) == 4,
         "complete four-colouring")
    need(all(witness4[v] != witness4[w] for v, w in edges), "proper whole-graph four-colouring")
    result = {
        "status": "PASS",
        "incoming_commit": commit,
        "imported_contact_file": SOURCE,
        "imported_contact_blob": git("rev-parse", f"{commit}:{SOURCE}").decode().strip(),
        "imported_contact_sha256": hashlib.sha256(raw).hexdigest(),
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "geometry_scope": "Imported audited complete contacts; independently rechecked 21 pinned unit rows.",
        "core_vertices": h, "core_edges": len(core_edges), "boundary_vertices": len(boundary),
        "a": a, "t": t, "pinned_rows": row_count, "pinned_rank": rank,
        "pinned_stress_dimension": stress, "six_colour_critical_bound": sharp_bound,
        "nonzero_minor_blocks": minor,
        "full_minor_determinant": [str(x) for x in full_determinant],
        "blocked_core_extension_count": 0,
        "all_seven_core_deletions": deletion_checks,
        "ordinary_chromatic_number": 4,
        "consequence": "The incoming 19-point gate attains the critical pinned-stress bound with a vertex-critical four-chromatic interior."
    }
    target = HERE / "incoming_gate_stress.json"
    target.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8", newline="\n")
    print(json.dumps({k: result[k] for k in ("status", "incoming_commit", "core_vertices", "a", "t",
                                            "pinned_rows", "pinned_rank", "pinned_stress_dimension",
                                            "six_colour_critical_bound", "ordinary_chromatic_number")}, sort_keys=True))
    print("Core-deletion extension counts:", [x["extension_count"] for x in deletion_checks])


if __name__ == "__main__":
    main()
