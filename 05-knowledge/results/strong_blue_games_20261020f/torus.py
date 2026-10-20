"""Independent explicit Euler circuits for marked positive torus grids.

The p=1 or q=1 cases retain direction-labelled loops; p=2 or q=2 retains
the two distinct torus edges between the corresponding adjacent vertices.
These are cellulation multigraphs, not a simple undirected K_n game.
"""
from collections import Counter
from hashlib import sha256
from math import gcd
from pathlib import Path
import json


def check(test, message):
    if not test:
        raise RuntimeError(message)


def circuit(p, q):
    edges = []
    outgoing = [[] for _ in range(p * q)]
    for i in range(p):
        for j in range(q):
            source = i * q + j
            for direction, target in [("E", ((i + 1) % p) * q + j),
                                      ("N", i * q + (j + 1) % q)]:
                outgoing[source].append(len(edges))
                edges.append((source, target, direction))
    stack = [(0, None)]
    reverse = []
    while stack:
        vertex = stack[-1][0]
        if outgoing[vertex]:
            edge_id = outgoing[vertex].pop()
            stack.append((edges[edge_id][1], edge_id))
        else:
            _, incoming = stack.pop()
            if incoming is not None:
                reverse.append(incoming)
    tour = list(reversed(reverse))
    check(sorted(tour) == list(range(2 * p * q)), "each edge once")
    current = 0
    directions = Counter()
    visited = {current}
    for edge_id in tour:
        source, target, direction = edges[edge_id]
        check(source == current, "continuous Euler circuit")
        current = target
        visited.add(current)
        directions[direction] += 1
    check(current == 0 and len(visited) == p * q, "closed and vertex spanning")
    check(directions == {"E": p * q, "N": p * q}, "direction counts")
    check(directions["E"] % p == directions["N"] % q == 0, "integer winding")
    winding = [directions["E"] // p, directions["N"] // q]
    check(winding == [q, p], "marked winding")
    raw = [abs(q * q - p * p), 2 * p * q, q * q + p * p]
    check(raw[0] ** 2 + raw[1] ** 2 == raw[2] ** 2, "Euclid identity")
    content = gcd(gcd(*raw[:2]), raw[2])
    if gcd(p, q) == 1:
        check(content == (2 if p % 2 == q % 2 == 1 else 1), "primitive parity content")
    return {"p": p, "q": q, "vertices": p * q, "Euler_length": len(tour),
            "winding": winding, "primitive_winding": gcd(p, q) == 1,
            "raw_Euclid_triple": raw, "triple_content": content,
            "reduced_triple": [entry // content for entry in raw],
            "Euler_edge_ids": tour}


def main():
    records = [circuit(p, q) for p in range(1, 10) for q in range(1, 10)]
    examples = [circuit(*pair) for pair in [(1, 2), (3, 4), (3, 5), (13, 21), (5, 12), (6, 10)]]
    check(examples[3]["Euler_length"] == 546 and examples[3]["reduced_triple"] == [136, 273, 305],
          "golden 273 parity sidecar")
    check(examples[-2]["Euler_length"] == examples[-1]["Euler_length"] == 120,
          "equal length hostile")
    check(examples[-2]["primitive_winding"] and not examples[-1]["primitive_winding"],
          "length forgets primitive winding")
    receipt = {"status": "FINITE-EXACT explicit circuits; general identities proved in torus.md",
               "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
               "parameter_universe": "all 81 ordered pairs 1<=p,q<=9, plus named examples",
               "tested_circuits": len(records) + len(examples),
               "simple_undirected_cellulation_boundary": "p,q>=3; smaller parameters retain loops/parallel torus edges",
               "examples": examples,
               "scope": "Circuits carry their torus embedding and E/N arc directions. No knot-isotopy or game-winning computation."}
    Path(__file__).with_suffix(".json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in receipt.items() if key != "examples"}, indent=2))
    print("example summaries:")
    print(json.dumps([{key: value for key, value in item.items() if key != "Euler_edge_ids"} for item in examples], indent=2))


if __name__ == "__main__":
    main()
