# Improved explicit distributions on the order-four native core

Status: independently audited exact certificate (PASS). This supplements the two-cube/six-square packet, not its lower bounds or structural claims. It uses no main theorem from arXiv:2611.05189v1.

Use precisely the labelled graph, core edge orbits A (48) and B (24), canonical capacities and48 relabelling/converse automorphisms defined and verified in session208/core-geometry. WITNESS.json gives two actual core spanning trees X,Y, each39 edges. Restore the24 fixed leaf channels for the full graph.

Their exact path-incidence sums, in order (stretch_A,stretch_B,load_A,load_B), are

    X: (186,36,122,100),
    Y: (140,96,168,68).

These trees were found by deterministic seeded edge-exchange searches optimizing weighted tree-response costs. Search optimality is neither asserted nor needed: the saved trees are standalone finite witnesses. Their validity and sums can be checked by exact path enumeration.

Uniformly relabel and optionally converse each seed, giving orbit distributions. For stretch mix X with weight26/83 and Y with57/83. Then each A-edge expectation is

    (26*186+57*140)/(83*48)=267/83,

and each B-edge expectation is

    (26*36+57*96)/(83*24)=267/83.

Thus S_*<=267/83, improving157/48. At most96 trees suffice.

For capacities, use X with16/55 and Y with39/55. Each A-edge expected canonical load is

    (16*122+39*168)/(55*48)=1063/330,

and each B-edge expected load is

    (16*100+39*68)/(55*24)=1063/330.

Thus C_*<=1063/330, improving13/4, again with at most96 trees. The smaller support48 of the earlier capacity distribution remains a separate advantage. As before leaves have stretch and load exactly1, and edge stretch extends to all pairs. Capacity domination and fractional routing implications carry over unchanged. The stretch and capacity mixtures intentionally have different weights.

These are improved feasible upper bounds, not optimal values, nor a general convergence claim for weighted edge exchange. The lower bounds19/8 and79/36 remain unchanged. The meaningful bridge to the paper is a small exact weighted-response search followed by rational mixture and symmetry averaging, with every final claim reduced to finite integer path data.
