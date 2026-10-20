"""Independent exact audit of the claw-strip reduction.

Universe: all 27 input/output triples on three labelled colours; all powers
through 12.  Two constructions of the one-module transfer matrix are compared.
No random tests or snark classification library is used.
Run: python3 04-computation/boundary_support_20261020_independent.py
"""
from collections import Counter
from itertools import permutations, product
import json

import networkx as nx
import sympy as sp


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


words = list(product(range(3), repeat=3))
index = {word: i for i, word in enumerate(words)}
T = sp.zeros(27)
# Directly test the two boundary triples: the missing colour at each leaf
# is the spoke colour, and the three spoke colours must be distinct.
for x in words:
    for y in words:
        if all(a != b for a, b in zip(x, y)):
            spokes = tuple(3 - a - b for a, b in zip(x, y))
            T[index[x], index[y]] = int(len(set(spokes)) == 3)

U = sp.zeros(27)
# Independent enumeration chooses the centre star first, then each leaf.
for spokes in permutations(range(3)):
    for x in words:
        if all(a != b for a, b in zip(x, spokes)):
            y = tuple(3 - a - b for a, b in zip(x, spokes))
            U[index[x], index[y]] += 1
require(T == U, "independent constructions disagree")

graph = nx.Graph()
graph.add_nodes_from(range(27))
graph.add_edges_from((i, j) for i in range(27) for j in range(i) if T[i, j])
component_types = []
for nodes in nx.connected_components(graph):
    part = graph.subgraph(nodes)
    component_types.append((len(nodes), part.number_of_edges()))
require(Counter(component_types) == Counter({(1, 0): 3, (3, 3): 2, (6, 6): 3}),
        "component classification")

square = T * T
square_graph = nx.Graph()
square_graph.add_nodes_from(range(27))
square_graph.add_edges_from((i, j) for i in range(27) for j in range(i)
                            if square[i, j])
blocks = [sorted(c) for c in nx.connected_components(square_graph) if len(c) > 1]
require(len(blocks) == 8 and all(len(b) == 3 for b in blocks), "eight blocks")
for b in blocks:
    require(square.extract(b, b) == sp.eye(3) + sp.ones(3), "I+J block")

powers = {0: sp.eye(27), 1: T}
for k in range(2, 13):
    powers[k] = powers[k - 1] * T
for r in range(1, 7):
    coeff = (4**r - 1) // 3
    for b in blocks:
        require(powers[2*r].extract(b, b) == sp.eye(3) + coeff*sp.ones(3),
                "closed power formula")
require(all(bool(powers[2][i, j]) == bool(powers[4][i, j])
            for i in range(27) for j in range(27)), "support reduction")
require(powers[2][0, 0] == 0 and powers[0][0, 0] == 1,
        "identity-wire hostile")

# Close the strip by transposing the last two strands.
P = sp.zeros(27)
for x in words:
    P[index[x], index[(x[0], x[2], x[1])]] = 1
closed_counts = {}
for k in range(1, 13):
    count = int(sp.trace(powers[k] * P))
    closed_counts[k] = count
    require(count == (0 if k % 2 else 6 * 2**k), f"twisted trace law k={k}: {count}")

z = sp.Symbol("z")
# The transfer denominator below is Q, not the user's unidentified D(z).
Q = (1 - 2*z)**5 * (1 + z)**10 * (1 - z)**6 * (1 + 2*z)**3
require(sp.expand(Q - (sp.eye(27) - z*T).det(method="domain-ge")) == 0,
        "characteristic denominator")

result = {
    "status": "FINITE-EXACT; closed formulas justified separately",
    "component_type_multiplicities":
        {str(k): v for k, v in sorted(Counter(component_types).items())},
    "even_power_blocks": [[list(words[i]) for i in b] for b in blocks],
    "support_sizes_powers_0_to_6":
        [sum(bool(v) for v in powers[k]) for k in range(7)],
    "twisted_closure_counts_k_1_to_12": closed_counts,
    "Q_transfer_degree": int(sp.degree(Q, z)),
    "Q_transfer_coefficients_low_to_high":
        [int(sp.expand(Q).coeff(z, i)) for i in range(25)],
    "support_equal_but_weights_not_proportional": {
        "T2_diagonal_offdiagonal": [2, 1],
        "T4_diagonal_offdiagonal": [6, 5],
    },
    "all_checks": "PASS",
}
print(json.dumps(result, indent=2))
