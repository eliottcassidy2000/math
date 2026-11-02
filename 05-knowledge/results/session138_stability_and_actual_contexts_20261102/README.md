# Session138: exact stability and actual contexts

**The optimal weighted determinant stability coefficient is now exact for every order `n>=1500`, with a complete positive-variance equality classification.** Together with the finite-order proofs, exact coefficients are established for `n=3,...,9` and `n>=1500`; orders **10 through 1499 remain unresolved**.

**Follow-up, 2026-11-02:** the [new matching-budget proof](../stability_priority_moduli_20261102/STABILITY_N276.md), incorporating incoming Session139, lowers the sufficient threshold to **276**, with the same equality classification. The current open interval is **10 through 275**. The ranges below describe this packet's independently audited n>=1500 checkpoint; its frozen proofs and audit outputs are preserved.

[Results and scope](SYNTHESIS.md) · [Verification and checkpoint history](VERIFICATION.md) · [Final replay: PASS](REPLAY_CHECKPOINT4.json)

## 1. Eventual exact stability, not only an asymptotic constant

For weighted tournament matrices, let `Delta=H-det(I+A)-|det A|`, `V=sum_(i<j) A_ij A_ji`, and `K_n=inf_(V>0) Delta/V`. Write `r=floor(n/2)` and

```text
kappa_n=(4/r)[1-(3/4)^r-1_(n even)*(1/4)^r].
```

**PROVED, independently audited:** for every integer `n>=1500`,

```text
K_n=kappa_n.
```

Equality with `V>0` consists exactly of an ordinal sum of the maximum number of fair pairs, with one singleton at odd order. The proof first forces the uncertainty graph to be a matching and the expected binary gap to vanish; the accepted H=1 minimization then supplies the exact equality locus. At these orders, V=0 equality consists exactly of deterministic transitive tournaments by the inherited binary theorem.

This is an all-order argument above a sufficient threshold, not a census at order 1500 or an extrapolation from small orders. The threshold is not asserted optimal. The proof uses already accepted stability and matching theorems directly.

[Exact theorem and equality](threshold_refinement/PROOF_N1500.md) · [Independent audit](independent/threshold1500_audit/AUDIT.md)

## 2. Established finite orders and the remaining interval

| Order n | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Sharp K_n | 1 | 3/4 | 7/8 | 3/4 | 37/48 | 87/128 | 175/256 |

At orders four through nine, positive-variance equality has the same maximal fair-pair form. Order three is exceptional: `Delta=V` for every weighted tournament matrix. **Exact optimality for 10<=n<=1499 is not claimed.**

Parity is structural: an odd singleton makes the full determinant vanish, whereas an even full pair cover contributes the extra subtraction `(1/4)^r`.

The final order-nine continuation retains a positive triangle reserve:

```text
Delta >= (175/256)V+(105/256)tau+C4/16+(5/8)C5+C6/4+C7+C8.
```

Its independent complete graph certificate covers 3,160,576 nine-vertex augmentations, and the reviewed separately affine reduction proves the real-weight inequality. Equality uses the positive triangle reserve—not a false claim that the graph polynomial has only isolated zero endpoints.

[Sharp low orders](stability/SHARP_LOW_ORDER_STABILITY.md) · [Sharp seven](stability/SHARP_SEVEN_GRAPH_CERTIFICATE.md) · [Sharp eight](stability/SHARP_EIGHT_GRAPH_CERTIFICATE.md) · [Sharp nine](stability/SHARP_NINE_GRAPH_CERTIFICATE.md)

## 3. Uniform bounds and quantitative mechanisms

For every `n>=4`, the earlier continuation proves, with `c0=(1-exp(-4))/4`,

```text
kappa_n/[1+n*kappa_n^2/(2*c0)] <= K_n <= kappa_n,
K_n=8/n+O(n^(-2)).
```

These uniform bounds still apply in the unresolved intermediate range. The retained binary baseline, multiplicative Pfaffian margins, and greedy matching extraction also give quantitative near-candidate probability and variance-mass controls. The new exact theorem strengthens those conclusions at sufficiently large orders; it does not retroactively turn the earlier quantitative argument into exact rigidity at every order.

The first checkpoint retains its separate multistar coefficient asymptotic to `4/n`, arbitrary binary cross-orientation matching-support theorem, uncertainty-graph certificates, fixed-parameter evaluator, and H=1 classification.

[Uniform asymptotic theorem](asymptotic_stability/PROOF.md) · [Independent resolvent mechanism](stability/RESOLVENT_ROBUSTNESS_CROSSCHECK.md) · [Matching support](stability/MULTISTAR_STABILITY.md) · [Uncertainty graph](rounding_graph/PROOF.md)

## 4. Actual contexts and precise no-zero families

Literal one-edge reflectors realize the strong tournament counts `6E(U)-H(U)` and `6E(U)+H(U)` for every nonempty U. Full-profile operators, a source-certified inverse, minimal exterior budget, and unbounded fixed-gap families are proved. Bounded actual height/interior source banks fail even profinite universality; this is not an unrestricted spectral omission.

The ordinal-profile theorem proves internal star-edit zero-rigidity in `C3[TT_k[C3],A,B]` when outside attachments stay fixed. A separate theorem classifies three explicit actual strong-Q4 mixed-profile seeds for every positive side multiplicity: two are positive on the thin boundary and negative throughout the thick quadrant, never zero; the third is positive for arbitrary nonempty side tournaments, with sharp lower gap `135618 H(A)H(B)`.

[Reflector](exact_resources/PROOF.md) · [Resource bounds](exact_resources/RESOURCE_BUDGET.md) · [Restricted profile rigidity](zero_rigidity/PROFILE_DICHOTOMY.md) · [All-parameter sign chambers](zero_balance_later/NO_ZERO_THEOREM.md)

The later reusable framework gives finite **sufficient** certificates for universal nonnegative cyclic response and for a negative polynomial-multiplier chamber. A genuine strong same-card pair of order seven shows the row-crossing condition cannot be omitted: its response changes from −36 with singleton sides to +828 with transitive two-vertex sides. Complete independent coverage through module order six proves this hostile's module order is minimal. Failed certification remains inconclusive; this is not a counterexample to scalar zero-rigidity.

[General criteria and minimal actual hostile](zero_framework_later/GENERAL_CRITERION.md) · [Independent framework audit](independent_zero_framework_later/AUDIT.md)

## Boundaries and checkpoint provenance

**OPEN:** exact stability coefficients at orders 10 through 1499; the unrestricted positive-odd H-spectrum; global boundedness of exact spectral width; arbitrary strong-quotient scalar zero-rigidity; and LRC(14).

The [first checkpoint](CHECKPOINT1.md)'s through-seven results and the [second checkpoint](CHECKPOINT2.md)'s through-eight/asymptotic results retain their original lane bytes, selectors, and audit scopes. The [third checkpoint](CHECKPOINT3.md) adds the exact large-order theorem and general cyclic criteria; the final order-nine continuation has its own separately curated proof and graph audit. No unaccepted finite-order graph candidate or historical threshold candidate is required for this release. A separately audited [Pfaffian hostile](stability/PFAFFIAN_UNIT_COEFFICIENT_HOSTILE.md) refutes an unused termwise shortcut; it does not invalidate any accepted theorem.
