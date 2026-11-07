# Session211: score permutohedra and the literal gate boundary

The proof gives an exact bridge between the old score-zonotope research and
Session210's canonical tournament gate:

- Tight score inequalities recover the ordered strong-component face, a
  product of smaller permutohedra. They do not recover the internal arcs.
- Adding the gate translates the entire source/sink face strictly into the
  interior, with unnormalized facet margin1; the full degree slice is a
  larger, explicitly identified permutohedron.
- Reversing a directed triangle preserves the exact score and changes the
  gated cycle count from H(T) to H(T-v), erasing the gate for every mixed v.
- A directed four-cycle exposes an endpoint-refined two-path response.
  Explicit equal-score pairs have HC11 versus7, and HC37 versus21.
  All four tournaments are strong; the 7 and21 realizers have no gate.

These are precise correspondences and obstructions, not a solution of the
Chan–Pak Hamilton-path spectrum problem. In particular, no tournament here
has H=7 or H=21. The score projection loses the class-defining microscopic
arc information, while the full coherent orientation retains it.

Reproduce with Python's standard library:

    python3 verify.py > /tmp/session211-scores.json
    cmp RESULT.json /tmp/session211-scores.json
    python3 -O verify.py > /tmp/session211-scores-O.json
    cmp RESULT.json /tmp/session211-scores-O.json

PROOF.md states the all-order results and the exact finite verification
universes. SOURCES.json pins inherited results and historical passages.
