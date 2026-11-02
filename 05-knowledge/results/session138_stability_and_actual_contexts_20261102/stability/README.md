# Session138: sharp and order-optimal weighted stability

This is the reviewed stability packet. It builds on Session137's global
weighted determinant inequality and signed insertion identities. All files
in `PUBLIC_FILES.txt` are intended for publication; exploratory searches,
executables, and pending order-eight work are intentionally excluded.

## Main results

- **All orders:** `MULTISTAR_STABILITY.md` proves
  `Delta >= E Delta(binary rounding) + gamma_(n,k) V`, where
  `gamma_(n,k)=4[n-k+2-(n+k+2)/2^k]/[n(n-1)]`, for every n>=3 and1<=k<=n.
  Optimizing k gives a coefficient asymptotic to4/n, improving the inherited
  order1/n^2 bound to the correct order1/n. Fair-pair examples give the
  asymptotic upper bound8/n, so the remaining universal constant gap is
  a factor two, not a solved optimality claim.
- **All-order extremal structure:** `H1_STABILITY.md` classifies H=1 as
  ordinal sums of weighted pairs and singletons and solves its exact
  variance ratio. The multistar note's appendix proves the same sharp
  coefficients on the larger class whose uncertainty graph is a matching,
  with an additional expected binary gap term. The generic graph-rounding
  framework lives in the companion `../rounding_graph/PROOF.md`.
- **Sharp strong coefficients:** orders3..7 are respectively
  `1, 3/4, 7/8, 3/4, 37/48`. The new order5,6 proofs are in
  `SHARP_LOW_ORDER_STABILITY.md`; order7 is in
  `SHARP_SEVEN_GRAPH_CERTIFICATE.md`. Equality with positive variance at
  orders4..7 is exactly a maximal number of fair pairs and at most one
  singleton, in ordinal sum.
- **Exact coefficient regions:** at order5, all real pairs(a,b) for which
  `Delta>=a tau+b V` holds are exactly
  `a<=1, b<=7/8, a/4+b<=17/16`. At order6 the exact region is the rectangle
  `a<=1, b<=3/4`. The order5 slanted facet has a continuous equality family,
  so its equality is not restricted to the two endpoint witnesses.
- **Sharp weak order6:** `SHARP_WEAK_SIX.md` proves
  `H-u>=tau+(37/48)V+(3/4)C4+C5+C6+3M33` through a small graph certificate.
- **Local bridge:** `TRIANGLE_VARIANCE.md` proves
  `tau>=sum_incident_edge_pairs q_e q_f>=2V^2/n-V/4`. It also preserves the
  failed order7 reduction that omitted the wedge-plus-disjoint-edge term.

Here `u=det(I+A)`, `v=(-1)^n det(A)`, `Delta=H-u-|v|`,
`V=sum_edges q_e`, `q_e=Aij Aji`, and `tau` is the total directed-triangle
weight. Cycles and disjoint cycle collections are counted as specified in
the proof notes; there are no silent labeled-cycle multiplicities.

## What is analytic versus finite-exact

The all-order theorems, H=1 structure, order5 `(3/4,7/8)` corner and order6
strong bound have analytic proofs. The order5 `(1,13/16)` corner uses a
complete tensor Bernstein certificate over all orientation orthants.
Order6 weak and order7 strong reduce to separately affine polynomials in
q, then use exact enumeration of all undirected graphs on6 and7 vertices.
These are continuous-cube proofs by polynomial interpolation, not weighted
matrix samples. Each load-bearing certificate has been independently
reconstructed by the audit lane with a different kernel.

## Reproduction

Python certificate scripts require NumPy; the graph checks require a C++17
compiler. No network, repository mutation, or external data is used.

```sh
python3 -B bernstein.py 5 7/8 > /tmp/bernstein5.replay.jsonl
python3 -B reserved_cycles.py 5 13/16 > /tmp/joint5.replay.jsonl
cmp bernstein5.jsonl /tmp/bernstein5.replay.jsonl
cmp triangle5_13over16.jsonl /tmp/joint5.replay.jsonl

g++ -O3 -std=c++17 graph6.cpp -o /tmp/graph6
/tmp/graph6 > /tmp/graph6.replay.json
cmp graph6.json /tmp/graph6.replay.json

g++ -O3 -std=c++17 sharp7_graph_certificate.cpp -o /tmp/sharp7
/tmp/sharp7 > /tmp/sharp7.replay.json
cmp sharp7_graph_certificate.json /tmp/sharp7.replay.json
```

The complete order7 universe is2,097,152 labeled simple graphs, no filters;
its sharp residual has minimum0 on421 graphs. The complete order6 graph
universe is32,768 graphs, minimum0 on76 graphs. The corroborative order6
Bernstein replay (`bernstein.py 6 3/4`) takes more memory and time but does
not create a large output file. Coefficient arrays are never stored.

## Open boundary

The exact all-order stability coefficient remains open beyond the proved
low-order scopes. Pending order8 work and further graph-rounding strategies
are not part of this frozen checkpoint. No statement here proves LRC(14)
or claims a novel literature priority.
