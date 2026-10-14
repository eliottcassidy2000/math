# A six-way seam relation in the exact Parts construction

**Date:** 2026-10-14. **Status:** FINITE-EXACT geometry and checked SAT witnesses; **UNSAT_SOLVER_VERDICT**, with a native but unchecked proof trace. No new graph record or higher lower bound for the plane is claimed.

## 1. Answer to the finite two-tag hypothesis

The [roadmap's Hypothesis H](exact_interface_research_roadmap_20261014.md#hypothesis-h-a-two-tag-seam-certificate) receives a positive result at the **solver-verdict level**. In the reproduced `L727 union rho*S361`, both possible additional tags work:

- all internal edges plus tag `(4,4)`: **SAT**, with a checked four-colouring;
- add tag `((17-sqrt(33))/6,(17+sqrt(33))/6)`: **UNSAT_SOLVER_VERDICT**;
- instead add the reversed tag: **UNSAT_SOLVER_VERDICT**.

These exhaust the occurring tags. The entire union has 1,087 vertices and 5,892 edges: 5,862 internal edges, 18 cross edges at `(4,4)`, and six cross edges at each of the two remaining ordered tags. Thus each two-tag test has 5,886 edges.

The stronger retained information is a minimal six-way boundary relation. For either additional six-edge matching, every single-edge deletion has an explicit checked four-colouring. All six inequalities together are solver-reported inconsistent. The obstruction is concentrated in six individually necessary matching constraints, relative to the fixed internal-plus-`(4,4)` control.

This is not proof-checked UNSAT. The existence and checking of SAT witnesses, including all single-edge deletion witnesses, do not turn the global solver verdict into a checked proof.

## 2. Reused exact construction and complete seam

The experiment imports the repository's [exact field and graph toolkit](../../04-computation/experiments/hn_records_20261014/hn_udg.py) and [type-M builder](../../04-computation/experiments/hn_records_20261014/type_m_probe.py), rather than rebuilding its coordinates, pieces or selector oracle. These two source files are imported unchanged from pinned main commit `c34aecb47ea9c1fc79eb0cd20e3fcd3bd66fb0c3`; the same blobs were independently verified at `24fbc7bc3f7a4409dfa4a02179522865d9b5c7f6`. Their Git blob identifiers are retained:

- `hn_udg.py`: `f97473fce08c4622cdf44cc289cf223ccf14c4fe`;
- `type_m_probe.py`: `680d7fe75acb778a9e8432f5aa96936cf871bca8`.

Let `H` be the centre and six vertices of a unit regular hexagon, and put

`eta=(sqrt(33)+i sqrt(3))/6`, `rho=(7+i sqrt(15))/8`.

The imported `build('rho')` uses the exact Parts pieces `L727` and `S361`; their only shared vertex after rotation is the origin. All edges internal to either piece are retained. Cross tags are ordered, with the unrotated `L727` endpoint first and the unrotated `S361` endpoint second. Rotation preserves the latter's norm.

The toolkit verifies candidate edges exactly. For this experiment its edge set was additionally checked against **all 590,241 unordered pairs**, using the rational component of the exact squared distance as an integer rejection test, followed by the toolkit's full multiquadratic norm test. There were 6,060 rational-component candidates; the resulting exact edge list agrees completely with the toolkit's 5,892 edges. No floating-point cutoff is relied on for this completeness check.

The certificate pins the exact vertex ordering, full edge list and internal edge list by hashes. The checker regenerates these from the pinned toolkit and verifies the hashes before using any indexed colour vector or CNF.

## 3. Geometry of the thirty cross edges

Write

`alpha=(17-sqrt(33))/6`, `beta=(17+sqrt(33))/6`,

`r_-=(sqrt(11)-1/sqrt(3))/2`, `r_+=(sqrt(11)+1/sqrt(3))/2`.

Then

`r_-^2=alpha`, `r_+^2=beta`, `r_- r_+=8/3`, `r_+-r_-=1/sqrt(3)`.

Let `Z={exp(i(pi/6+k pi/3)): k=0,...,5}`. The two six-edge matchings are exactly

`E_-+={(r_- z, rho r_+ z): z in Z}`,

`E_+-={(r_+ z, rho r_- z): z in Z}`.

The checker confirms these closed formulas against all twelve indexed edges. Each has squared length

`alpha+beta-2*(7/8)*r_- r_+=17/3-14/3=1`.

The eighteen radius-two edges are exactly `(u,rho u)`, where

`u in 2*(eta^-1 H union H union eta H)`, `u!=0`.

They have squared length `2*4*(1-7/8)=1`.

### Why the seam is radial before rotation

Both unrotated pieces lie in the Moser field `K=Q(sqrt(3),sqrt(11))`. In

`|u-rho v|^2 = |u|^2+|v|^2-(7/4)(u dot v) +/- (sqrt(15)/4)det(u,v)`,

all coefficients except the displayed `sqrt(15)` factor lie in `K`, and `sqrt(15)` is not in `K`. A unit-distance equality therefore forces `det(u,v)=0`. The exact finite census then identifies which collinear norm pairs actually occur: precisely the three tags above. The collinearity argument alone is not claimed to enumerate the tags.

The added six constraints couple two **different** radial rings. Radius-two forcing alone leaves a checked colouring; the mixed-radius matching removes every colouring according to the solver. This identifies the interaction more precisely than the full graph's chromatic number.

## 4. The two additional tags are Galois conjugate

The field automorphism `sigma` fixing `sqrt(3),sqrt(5)` and sending `sqrt(11)` to `-sqrt(11)` fixes `rho` and sends `eta` to `-eta^-1`. Since the hexagonal wheel is centrally symmetric and the relevant exponent sets are symmetric, the imported point sets `L727` and `S361` are preserved. Thus `sigma` permutes the vertices of the union and preserves unit edges because it fixes the equation `distance^2=1`.

It fixes tag `(4,4)` and swaps the other two ordered tags. The checker explicitly verifies the entire vertex permutation, all edges and the tag swap. This is a real-field graph automorphism; it is not asserted to be an ordinary Euclidean isometry of the plane.

Both two-tag instances were queried separately within the incremental selector run. The exact conjugacy explains why their answers agree and transfers every checked deletion witness from one matching to the other.

## 5. The actual six-way boundary relation

Let `A` be the graph containing all internal edges and the eighteen `(4,4)` cross edges. Label the six pairs in `E_-+` by `e_0,...,e_5`, and for a proper four-colouring `c` of `A` define

`b_j(c)=1` iff the endpoints of `e_j` have the same colour.

The solver-supported statement is

`b_0(c) OR ... OR b_5(c) = 1` for every proper four-colouring of `A`.

The retained SAT control colouring has exactly one `b_j=1`. Rotation by 60 degrees preserves `A` and acts transitively on the six pairs. Rotating that colouring supplies a checked extension realizing **each singleton** monochromatic pattern. Consequently:

- the all-zero pattern is solver-reported impossible;
- all six singleton patterns are explicitly realizable;
- adding any five of the six matching inequalities leaves a checked four-colouring;
- adding all six is solver-reported UNSAT.

The checker tests each complete colour vector against every retained edge and verifies equality on the one deleted pair. It also transports and checks the corresponding witnesses for `E_+-` by the Galois permutation.

No claim is made here about which nonsingleton Boolean patterns are realizable. This is a precise boundary relation with a six-way exclusion, not a complete classification of all seam colourings.

## 6. Solver and native-proof provenance

The incremental run uses the toolkit's `Colourer` with all vertex selectors active, plus one new selector for each complete cross-edge tag. Every inactive tag is explicitly disabled. A retained internal triangle is fixed to colours 0,1,2; this is legitimate global colour symmetry breaking, since every proper colouring gives that triangle three different colours.

CaDiCaL 1.9.5 through python-sat returned:

- internal-only control: SAT;
- internal plus `(4,4)`: SAT;
- either two-tag instance: UNSAT solver verdict;
- all three tags: UNSAT solver verdict.

A separate **standalone** solve of the first two-tag instance also returned UNSAT. Its committed DIMACS input has 4,348 variables and 24,634 clauses, with SHA256

`3706c283b4ff48ddca7f58e4e7856517242dbd16bc81da13808192d6020f8ab8`.

It uses one at-least-one-colour clause per vertex and one incompatibility clause per edge and colour, plus the triangle units. At-most-one clauses are unnecessary: choosing any true colour at each vertex yields a proper colouring if the CNF is satisfiable.

The native DRAT trace has 227,834 lines and 27,281,530 bytes; its SHA256 is

`94e1fd8498929ca560626ca11341387e10269874940738fbb92f3a9ce71e69a1`.

Deterministic gzip compression gives 7,372,290 bytes with SHA256

`5f42984c9db07882533a5d15fbe575ae321ce5cf8a6b354310641cf515ff5212`.

The trace is **unchecked and retained only in the local recovery bundle**. It is deliberately omitted from the Git package; no remote availability is claimed. The compact exact CNF, colour witnesses, provenance and reproduction code are included. A future proof-checking step can promote the UNSAT status only after validating a trace against this exact input.

## 7. Reproduction and bounded conclusion

Files:

- `04-computation/experiments/hn_records_20261014/two_tag_seam_probe_20261014.py`
- `05-knowledge/results/hn_two_tag_seam_certificates_20261014.json`
- `05-knowledge/results/hn_two_tag_seam_check_20261014.json`
- `05-knowledge/results/hn_two_tag_seam_minus_plus_20261014.cnf`
- this note

The tested toolkit dependencies are numpy 2.3.5, scipy 1.17.0, python-sat 1.9.dev15 and its dependency six 1.17.0. The two imported source files are the only repository-code dependencies. From repository root:

`python3 04-computation/experiments/hn_records_20261014/two_tag_seam_probe_20261014.py --verify`

`python3 -O 04-computation/experiments/hn_records_20261014/two_tag_seam_probe_20261014.py --verify`

Both outputs match the declared check JSON byte-for-byte. This mode verifies exact geometry, both SAT controls, all twelve transferred single-edge-deletion witnesses, the Galois map, the ring formulas and the committed CNF. It does **not** check the UNSAT trace.

Use `--solve` to rerun the finite tag census with incremental selectors. Use `--emit-proof --output-dir <local-directory>` to regenerate a native standalone trace and its deterministic gzip form. Proof bytes can depend on solver/library versions; the retained hashes identify this run.

The finite tag question is exhausted: there are only two additional occurring tags, and both have positive solver-verdict outcomes. The new mechanism is the six-way mixed-radius boundary relation, with checked witnesses for every singleton and a proof-checking obligation for the all-zero exclusion. No larger graph search or record attempt was undertaken.
