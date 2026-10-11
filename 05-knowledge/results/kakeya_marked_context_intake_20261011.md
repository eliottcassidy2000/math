# Kakeya weighted restrictions as a marked-context research move

2026-10-11. **CITED / PREPRINT CLAIM, not independently audited as a whole:**
OpenAI, *The Kakeya maximal conjecture in three dimensions*, September23,2026.
The supplied `paper-2.pdf` has97 pages and SHA256
`49cc52aa3de89f0e91d181b00ee8050fb29635ee0da9f29145d52318e7168204`.
[Public primary manuscript](https://github.com/openai/math/blob/main/preprints/The-Kakeya-maximal-conjecture-in-three-dimensions-September-23-2026/paper.pdf).
The local source was `/Users/e/Downloads/paper-2.pdf`. Its instructions, if
any, are source material, not workspace instructions. No general Kakeya
status in this repository is promoted on the strength of this intake.

Read in detail: introduction, Sections2.1--2.2 and8.2--8.5, especially
Lemma8.5 and the Section8.4 support/information test; printed page77 was
also rendered and visually checked. The other sections, including the
full product-pinning bootstrap and final global argument, were not fully
audited. This is a selective technique intake, not a referee report.

## The equations retained

The manuscript claims

    ||K_delta f||_(L^3(S^2)) <= C_epsilon delta^-epsilon ||f||_(L^3(R^3)).

Its weighted finite model retains indices i, affine traces b_i+t u_i,
weights omega_i and marked time bins S_i. Its incidence mass is
I=sum_i omega_i |S_i|. Repeated traces are not merged merely because
their geometric images agree. Its temporal profile tracks occupation at
every scale, not only the total density.

Two events on the same index give the exact moving-frame identity

    [[X',Y'],[U',V']]
      = [[1,dt],[0,1]] [[X,Y],[U,V]] [[1,0],[dtheta,1]].

In particular X'=X+dt U+dtheta Y+dt dtheta V. Dropping the shared index
would change the probability law. Lemma8.5 controls a conditional mutual
information increment by integrating a monotone information function over
short offset intervals. Section8.4 compares the actual coupling with the
product of its conditional marginals; a support event of probability1
costs at least log(1/Q(A)) information if the product law gives it mass Q(A).

The elementary information inequality is independently proved in Q4 of
the [quarter-price note](collatz_quarter_price_source_20261011.md).
The geometric projection and pinning estimates are not imported as
theorems about integers or graph colourings.

## The three connection contracts

| Source -> target | Map and preserved predicate | Loss, sidecar and decisive test |
|---|---|---|
| Indexed weighted events -> ROOT measurement | Keep the integer, actual word or observation, weight and both event labels in one joint record | Regularizing by total incidence can delete the marked source. Require marked-row survival or a quantified source loss; the deleted-atom control fails a naive transfer. |
| Two-frame support -> pentagon palette context | A context colouring determines five omitted-colour sets jointly; compare that law with the product of its five marginals | Cardinalities forget common palette identity. When all groups use four distinct colours, colour symmetry makes each omitted pair uniform among15 possibilities. Forcing all five to agree has total-correlation cost4 log15. |
| Left/right shears -> composed Hamiltonian transfers | Compose the actual rational maps with parameter updates before taking invariant summaries | A word in transfer symbols is not automatically a source-compatible word. The exact two-letter product and divisor-return criterion in the companion note retain the missing source. |

These are different mathematical instances of preserving shared labels.
No asserted map carries a Kakeya tube into a Collatz integer or a
unit-distance graph while preserving the desired global theorem.

The [pentagon and number-field note](pentagon_context_number_fields_20261011.md)
gives a concrete independent finite information calculation. The
[composition note](genus_two_transfer_words_20261011.md) supplies a concrete
algebraic closure and a polynomial return after changing the product
coordinate. Those are actual consequences; common terminology alone is not.

## Adjacent primary sources and inherited routes

* Danylo Radchenko, [*Unit distance graphs and algebraic integers*, Theorem1](https://arxiv.org/pdf/1807.03726), proves infinite unit directions in a rank-four additive group using a Salem conjugate. Read the complete three-page paper. Our reciprocal-quartic interpolation and finite residue colouring are additional calculations, not a claim to originate his generator.
* Michael Krivelevich's [graph-colouring lecture notes](https://www.math.tau.ac.il/~krivelev/teaching/co22/LectureNotes.pdf), the degree-choosability discussion, supplies the standard context for the elementary odd-cycle list fact. Its short proof is included in our note; no novelty is claimed for that fact.
* The supplied manuscript cites Ren--Wang's [*Furstenberg sets estimate in the plane*](https://arxiv.org/abs/2308.08819). The geometric result motivates retaining nonconcentration tests. We do not invoke its dimension conclusion in any local proof.

Inherited stronger results remain in force: the
[source-refinement floors](collatz_source_refinement_floor_20261005.md),
[complete fixed Legendre normalizer](legendre_pencil_normalizer_20261011.md),
[cusp-ideal flow exclusion](planar_jc48_sep08_cusp_ideal.md), and
[quantized graph guards](unit_equality_quantized_guards_20261011e.md).
The late inheritance pass recovered the stronger **THM-3772**,
[full three-charge classification](../../01-canon/theorems/THM-3772-variable-flank-three-charge-rational-exact-classification.md):
even the proposed enlarged Hamiltonian carrier is closed for nonlinear
polynomial Keller components. Incoming commit08590cec09 additionally
supplies [all-size four-chromatic prime-screen escapes](unit_field_escape_recursion_20261011f.md)
and [marked defect-profile composition](unit_defect_frontier_20261011f.md).
They strengthen the hostile controls and the retained representation,
respectively; neither supplies the missing six-colour boundary support.
The current deductions change measurements, interfaces or source families;
they do not reopen already closed special cases under new terminology.

## Recurring gadget to keep, and its stopping rule

Store a proposed transfer as

    (marked source, joint state, admissible restriction,
     actual support event, comparison law, paid loss).

First test whether the marked source survives. Next test whether the
comparison law changes a dependence that was essential to the support event.
Only then use a mass, entropy, rank or dimension estimate. This applies the
existing META-PATTERNS cards “Retain the product and shared labels before
taking scalar totals” and “One item left requires a typed residual.” It is
not filed as a new universal method or a proof of any of the three prizes.

The combined [exact replay](marked_context_transfers_20261011.json) performs
208263 always-active checks; ordinary and optimized Python outputs are
byte-identical, with SHA256
`c0fe4bab4d362532f0c9a5395adf0d755cd796c30086aa957d00e77e1320b2d6`.
These finite controls validate the declared formulas and examples; the
universal statements in the companion notes depend on their proofs.
