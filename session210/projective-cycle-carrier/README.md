# Projective cells behind (H-1)/2

A structured geometric carrier, not a proof that all values except 3 and 10 occur.

[The proof](PROOF.md) turns signed disjoint odd-cycle packings, modulo global sign reversal, into a finite regular CW complex with exactly `(H(T)-1)/2` cells. Its full cell grading remembers the packing polynomial. Its canonical mod-two cohomology class has height equal to the directed-triangle packing number.

For the inherited leaf grammar, r triangle SCCs give RP^(r-1); one S4 center plus r triangles gives RP^r wedge S^r (two points when r=0). A same-order H=9 example separates four isolated points from a circle, proving that even H together with tournament order loses the topology.

The exact replay checks every labelled tournament through order six: 19 H-values but 27 packing polynomials at order six. Run `python3 verify.py` or `python3 -O verify.py`; all checks survive optimization. [Results](RESULT.json) include explicit witnesses and quotient-cell boundary homology controls.

The odd-cycle formula and leaf classification are inherited. General octahedralization and projective-space cohomology are standard. No literature-priority claim is made, and arbitrary simplicial complexes have every possible cell total, so tournament realizability cannot be omitted.
