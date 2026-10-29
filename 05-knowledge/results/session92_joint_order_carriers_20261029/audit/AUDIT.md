# Session 92 independent audit

Source pin: `b07ec51b204e5c0e09e52eed206cb2bd3cc29ebc`.
Status: **FINITE-EXACT independent computations; scoped proof audit**.
No external priority or major-conjecture claim.

## Verified independently

### Folded W9 observer carrier

`check_folded_independent.py` constructs lexicographic C3[C3,C3,C3] directly,
with involution (i,j) ↦ (−i,−j), rather than importing producer code.
Four involution-pairs give 64 subset/current states and 96 transitions.
Backward recursive signatures preserving level, terminal flag and outgoing
multiplicities give 16 cells, distributed by levels as 3,7,4,2. Adding
`current vertex == 3` gives 35. Backward induction on the acyclic levels proves
these signatures give the coarsest allowed partition, not just a stable one.
Ordinary full Hamilton DP gives 3159; folded paths give 33.

Independent SymPy rational elimination confirms the pullback span of all base
colour indicators has rank 12. **Products of degree at most two of a basis of
that space already span rank 16.** This is a finite exact nonlinear-closure
fact, not a statement that 16 is the minimum dimension of any linear algorithm.

### Native-complement seam family

`check_seams_independent.py` uses edge sets and all 120 permutations for each
five-vertex partner, independently of the producer's subset DP. It verifies
H=11 and15. For shared-start padding it checks exact edge equality to the
universal-source join through order 13, which explains H preservation for all
sizes. The original padded tournament is strong because its old strong core
reaches the new vertex and is reached through the old shared start.

**Scope:** this substitute does not claim identical individual cursor replies,
nor does it replay the user's 2483993/2491877 or 2899675/2906091 examples.
Their exact witness data has not been recovered at this audit checkpoint.

### Arithmetic source multiplicities under regrouping

`check_divisor_lca_independent.py` verifies all 256 ordered four-factor rows
from {2,3,4,6}, under all five binary tree shapes, with 17925 fibre/determinant
checks. Conditional split weights C_L(a)C_R(b)/C_v(ab) telescope to uniform
leaf-tuple probability 1/C_root(d). They therefore preserve the native divisor
measure after multiplying by d. No coprimality assumption is used.
For A=B=p the corrected 2×2 value-weight array has determinant3p²/4, obstructing
a rank-one local tensor representation. This does not obstruct arbitrary
algorithms or imply a general computational lower bound.

### Signed boundary flux

`check_flux_independent.py` independently constructs the incidence boundary
for every graph K of orders1..4 (75 graphs). Elimination over F101 supplies
a lower bound on rational incidence rank; one explicit integer kernel vector
per independent subset supplies the matching upper bound. Thus these finite
cases are exact over Q, not merely modular evidence. All give3^m−I(K).

The general proof checks out: fixed occupied subsetS produces an exchange
cube; elementary-square flux is−1 exactly on edges of K[S]. A negative square
annihilates any kernel amplitude, propagated across the connected cube. IfS
is independent, the alternating-parity vector spans its kernel. Boundary
hypercube coefficients can be assigned symmetrically because each incidence
edge has a unique forbidden endpoint. Over R/C, rank(LL*)=rank(L). This last
identity is **not** exported to finite-field Gram matrices. These are
state-dependent couplings, not unchanged on-site physical drives.

## Reproduction

Run each of the following eight scripts normally and with `python3 -O`:

1. `check_folded_independent.py`: W9 64-state recurrence, 16/35 quotients,
   rational observer ranks and all 256 named-vertex menus; **requires SymPy**.
2. `check_algebra_controls_independent.py`: six-/seven-state closure controls;
   **requires SymPy**.
3. `check_seams_independent.py`: 240 base permutation checks and source-join
   edge equality through order 13; standard library.
4. `check_parity_independent.py`: all 1098 labelled tournaments orders2..5;
   standard library.
5. `check_divisor_lca_independent.py`: 256 factor rows, five bracketings,
   17925 fibre/determinant checks; standard library.
6. `check_flux_independent.py`: all 75 graphs orders1..4; standard library.
7. `check_hyperflux_independent.py`: 30 hypergraph families; standard library,
   importing the audit flux checker and consequently repeating its 75 graphs.
8. `check_fano_independent.py`: all 128 active subsets and the simplex-code
   distinction; standard library.

All eight pass both modes. Checks use explicit exceptions, not assertions.
The current `normal.jsonl` and `optimized.jsonl` streams agree byte-for-byte.
The scopes above use different units and are not summed into one check count.

## Outstanding provenance boundary

The cited `joint_port_germs_20261028/README.md` concerns movable selector
registers and LRC germs; it does not contain the quoted large Hamilton counts.
A missing witness is an unverified inheritance item, not a refutation.

## Additional closure audits

- Named W9 linear pullback space has rational rank21 and34 distinct evaluation
  signatures, still short of the35 equitable cells. States(mask8,vertex7) and
  (mask4,vertex8) agree on every linear future-colour count. If f indicates
  terminal/reject/not-vertex3, the nonlinear observable A((A²f)²) is3 versus9.
  Independently verified by exact elimination and direct transitions.
- `check_parity_independent.py` checks all1098 labelled tournaments orders2..5
  by permutation constraints, not producer subset DP. The full spanning-path
  restoration coefficient is (1+(-1)^(n−1)) times the product of native path
  orientation signs. Thus odd n has the stated parity-ensemble expectation
  gap; even n cancels. The richer operation menu is explicit.
- The all-edge polynomial exact degree argument audits as n−1 for odd n≥3,
  n−2 for even n≥4; n2 has constant H=1 and degree0 separately. The even
  lower-bound proof uses a path on n−1 vertices and the remaining singleton:
  four traversal/placement bridge indicators sum2. The n2 boundary must not
  reuse this four-distinct-permutations argument.
- `check_hyperflux_independent.py` adds30 hypergraph cases: all families of
  hyperedges of size≥2 on2/3 vertices, every singleton hyperedge on4 vertices,
  and the complete4-vertex family. It uses the independent rank routine and
  exact integer kernel bounds. The general hypergraph extension proof audits:
  highest distinguished-index squarefree polynomial has a nonzero discrete
  derivative whenever any hyperedge survives on the occupied subset. This
  supplies a negative square, hence zero kernel in that component. Independent
  subsets have the same parity kernel as before. A full-set hyperedge changes
  one boundary coefficient and is invisible in proper occupied restrictions.

These additional controls are included in the eight-script reproduction
inventory above. The hyperflux import repeats 75 graph cases; they are not
additional independent graph instances.
- Independent backward signatures now also replay all 256 named-vertex menus;
  the six inclusion-minimal full64 menus are exactly126,190,222,237,245,249.
  Both the stated submodularity and supermodularity counterexamples pass.

## Final bounded publication checks

The parent `algebra/JOINT_CLOSURE_ALGEBRA.md` theorem is sound over Q:
finite rational hyperplane avoidance supplies a separating observable in any
unital pointwise subalgebra; Lagrange interpolation supplies every cell
indicator. Transition invariance then exactly expresses weighted row-sum
equitability. Inclusion of algebras reverses partition coarsening as required.
The constructive closure terminates by dimension, without a bit-complexity
claim. Parent basis insertion and snapshot-closure controls were inspected.
An additional independent SymPy checker verifies the six-state linear/joint
dimensions4/6, correlation1/0 and seven-state linear blindness with variance
separator4/2; normal and optimized modes pass.

`check_fano_independent.py` independently enumerates Fano lines by XOR in
F2³\{0}. There are7 minimal forbidden triples and64 allowed active subsets,
with size counts1,7,21,28,7. The audited general hypergraph theorem therefore
gives boundary/real-Gram rank3^7−64=2123. This is not a direct2187-dimensional
matrix elimination. The seven four-element faces are complements of Fano
lines; with zero they are the eight words of the binary simplex[7,3,4] code.
The full64-face family is not closed under symmetric difference and is not
a linear code. No equivalence with the folded64-state carrier is inferred.
Normal and optimized enumeration agree.

The current `runs/SEAM_CARRIER.md` now explicitly restricts the four-traversal
even lower-bound argument to n≥4 and treats n2 separately as H=1, degree0.
The previously flagged boundary issue is repaired.
