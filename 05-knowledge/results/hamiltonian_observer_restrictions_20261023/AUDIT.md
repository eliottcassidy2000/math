# Validity and scope audit

Status: **AUDITED** on the declared carriers. The exact controls support the proofs below; their finite ranges are not substituted for an all-size argument. No theorem namespace or reserved statement is used as a new dependency.

## Proof reviews

- **Native two-component incidence theorem:** root and the coupling author developed the degree-based descent together. Independent review is recorded in [compiler/AUDIT_NATIVE_DESCENT.md](compiler/AUDIT_NATIVE_DESCENT.md). The crucial missing premise in the initial attempt was the outgoing Hamiltonian boundary arc. The proved repair selects a nonneighbor of the degree-two-or-higher center; matchings, spanning stars and complete graphs are handled separately. Exact controls replay every nonzero labelled compressed graph through native order six, not merely representative counts.
- **Restriction and soft weights:** root reviewed the conditional-variance identity, weighted PSD inequality and bounded-density comparison. The order-four native kernel has centered eigenvector `(4,-3,-3)` and eigenvalue `8/15`, so removing transitivity still fails the proposed `8/9` bound. Unchanged-space coarsening and hard conditioning are explicitly separated.
- **Conditional compiler:** [history/AUDIT_COMPILER.md](history/AUDIT_COMPILER.md) independently checks shared-prefix constraints, right source-kernel action, target-prefix realization, targeted edge-cover optimality and the visible/invisible guard product. An audit corrected the selector scope: a retained distinguished reference path already breaks symmetry; the equivariant-selector obstruction applies to an unpointed fibre and its observations.
- **Trace and clocks:** [history/AUDIT_TRACE.md](history/AUDIT_TRACE.md) reviews the weighted harmonic extension, PSD proof, correct stationary conditional clock and actual owner-order lift. Independent checks use a nonuniform lazy-refresh closed form and all 64 labelled four-vertex tournaments.
- **Cycle/gain history:** [coupling/AUDIT_HISTORY.md](coupling/AUDIT_HISTORY.md) audits the original cycle normalization and gain lift. Its recorded check count predates the additional integer-clock corollary, which root separately reviewed: each turn divides by 13 and multiplies by 4, so a positive integer source permits exactly its 13-adic valuation in complete turns. Signed rational reconstruction additionally needs the sign. The final checker includes this extension.

## Exact universes and controls

| Lane | Universe / positive control | Hostile / retained boundary |
|---|---|---|
| Incidence coupling | Every fixed-path word at native orders 3 through 6; full root-marked cubes through 5 | Transitive singleton; nontransitive order-four gap `7/15`; paired gap `8/15` |
| Native descent | All `1+7+63+1023` nonempty labelled compressed inputs through native order 6; every intermediate Hamiltonian junction | Four-prefix exit guard, matching preparation without falsely charging a descent, complete/star exceptions |
| Compiler | Actual R9, mixed 27-vertex cells and regular 63-vertex source; all 32,768 target-slot sets | Nine flags required at order 63; disconnected guard; reference-path selector scope |
| Guard footprint | All 512 acceptance sets of a two-coordinate actual R9 owner observation | Visible disconnected factor remains disconnected after adding every invisible lift |
| Trace | All 255 nonempty subsets of eight words, seven unequal-weight subsets, 46 nonempty actual owner subsets at order four | Rejection gap `1/4` versus return gap one for the same ambient refresh and acceptance set |
| Gain lift | Exact valuation chains, rational labels, mod-8/mod-9 unit records, 597 positive integer sources | Changed centralizer normalization; joint-register hostile 73; all-permutation mod-9 cancellation fails on a Hamiltonian subset |

All checker failures use explicit exceptions; optimized Python does not remove the checks. Algebraic matrix tests use exact rational arithmetic, never a floating-point eigenvalue tolerance. Source-class and target-class identities in the native descent are independently checked using full permutation orbits.

## Type and quotient audit

1. A Hamiltonian path supplies an order, not a canonical invariant choice. Full owner lists and integer positions are retained wherever an actual port is required.
2. The auxiliary tournament has one fewer vertex than the native tournament. Their class observations are never identified, although their free-bit domain is the same.
3. The native incidence descent may change the native source between moves when the compressed graph is repositioned. It does not preserve fixed-source ports or electrical guards throughout.
4. Automorphism-frame moves preserve the complete ordered arc word. They change owner lifts, not the tiling word or target owner class. Their coordinate parameterization does not make the group law abelian.
5. A guard factorization must be authenticated before invisible coordinates can be removed from its search. No electrical mod-8 fingerprint is replaced by a residue of a chart count or a path clock.
6. The trace theorem preserves accepted endpoints; actual excursions may leave the source/port fibre. Return time and ambient time have different units. The duration matrix and harmonic extension remain available.
7. The valuation quotient has five base states under one normalization and two under another. The new gain decoration retains the endpoint ratio, but is not the historical arithmetic map or an automatic covering map.
8. The log branch cut changes an argument chart. It cannot cancel the real logarithmic defect. A gain graph does not determine a spherical/Euclidean/hyperbolic rotation system.

## Reproduction

From the repository root:

```text
python -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/reproduce.py
python -O -B 05-knowledge/results/hamiltonian_observer_restrictions_20261023/reproduce.py --optimized
python -B agents/check_docs.py
```

The reproduction driver verifies inherited file hashes and checks every frozen receipt after execution. The current-frontier edit only adds the new route, retaining the explicit open status of guarded wiring and Collatz. The two inherited CRLF-only TSV status entries are unrelated and are not staged.

Final replay: all five checkers pass in normal and optimized modes with byte-identical frozen receipts. The documentation check passes its bounded-startup, status, correction and router checks; all local Markdown links in this packet resolve.
