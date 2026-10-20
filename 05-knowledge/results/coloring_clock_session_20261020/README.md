# Coloring, finite clocks, and information lost by projection

Research session: 2026-10-20. This map joins the two supplied blueprints, binary LFSRs, the repository's tournament source calculus, an underused arithmetic Hankel sequence, and the proposed sphaleron–supernatural-number connection. Proofs, source attribution, exact checks, and unresolved gates live in the linked packets. No major open conjecture is claimed resolved. “Proved” refers to the precisely scoped statements in those packets, not a historical novelty claim.

## Main results and where to read them

1. **An exact finite-plane strengthening.** In the norm-one graph on F_17[sqrt(3)], every binary coloring has at least four monochromatic unordered triangles and at most 1730 cut edges. An independently replayed finite proof certificate strengthens the base's lower bound of two. Every distinct three-point template has at least two monochromatic ordered placements among its 5202 translated, norm-rotated placements. The exact minimum remains open. This finite-field graph cannot be transported edge-preservingly into the real unit-distance graph, so this is not a Hadwiger–Nelson improvement.
   - [Plane packet](../plane_binary_transfer_20261020/README.md)
   - [Four-defect addendum](../plane_binary_transfer_20261020/addendum-four-defects/README.md)

2. **Cycle graphs can erase primary structure.** For a binary companion update, one incoming edge per state holds exactly when the characteristic polynomial has constant term one. Fixed-point counts follow gcd(f,x^m−1), and a complete dyadically sampled primary profile classifies the cycle graph. Degree 18 is the sharp first cyclic-companion ambiguity in primary exponent multiplicities after grouping irreducibles by order: p²q⁴ and p³q³ have the same directed cycle graph but different primary structures, for the two irreducible cubics p,q over F_2. Hasse jets retain the missing thickness.
   - [LFSR and tournament packet](../binary_lfsr_tournament_clocks_20261020/README.md)

3. **Tournament source adjunction builds a genuine 2-adic clock.** On completed formal class coordinates, U=I+S has every nonzero orbit closure a pointed Z_2 torsor; the empty-class orbit reads binary digits through binomial coefficients. The correct inverse system uses dual formal coordinates, since raw statistic truncation need not commute with source pullback.
   - [Source-adjunction odometer](../binary_lfsr_tournament_clocks_20261020/TOURNAMENT_SOURCE_ADJUNCTION_ODOMETER.md)
   - [New functional-graph and regularity synthesis](FUNCTIONAL_GRAPHS_AND_REGULARITY.md): an explicit Borel proper two-coloring, continuous away from zero, but no proper discrete-valued coloring continuous at zero.

4. **An old arithmetic sequence hides a precision-dependent register.** The normalized THM-3484 sequence has rational recurrence degree 22, binary recurrence degree 6, and exact minimal recurrence orders 6,12,14,18,20,22 as 2-adic precision increases through specified thresholds. Exact Smith data and prefix-image equalities prove the full staircase. Its orbit clock is Z/3 × Z_2, with period 3·2^k modulo 2^k; its full state lattice has rank 22. These are different invariants.
   - [Observable Hankel clocks](../observable_hankel_clocks_20261020/OBSERVABLE_HANKEL_CLOCKS.md)

5. **Pair data end at three roots; even stars recover the rest.** Rooted minimum monochromatic-edge cost composes by min-plus elimination. Pair correlations determine every globally sign-symmetric cost on at most three roots; a four-leaf star distinguishes parity laws with identical pair marginals. All even-star probes, together with constants, form a triangular basis recovering the entire symmetrized boundary law. A finite minimax argument supplies an explicit positive-slack patching bound rather than assuming exact one-copy gluing.
   - [Rooted graph costs and tomography](../rooted_graph_costs_20261020/README.md)

6. **The valid sphaleron bridge is a completion of winding.** A degree-one loop in the specified SU(2) gauge–Higgs configuration model supplies a cyclic winding subgroup. Its profinite completion is Z-hat, isomorphic to Gal(F_2-bar/F_2), and has greatest supernatural order product_p p^infinity. The physical path, its energy, and its midpoint are not identified with this totally disconnected group. Primitive field multiplication sees only odd clock orders; Frobenius in normal bases supplies the missing even orders.
   - [Typed sphaleron–Steinitz bridge](../sphaleron_supernatural_bridge_20261020/SPHALERON_STEINITZ_BRIDGE.md)

## Further themes developed from the attachments

- **Binary transfer with no weighted loss:** independent cell thresholds turn a bounded measurable field into signs while preserving the selected expected triangle objective. This needs actual bounded-field realizability; positivity of a correlation kernel alone is weaker.
- **An arity boundary:** triangle monochromaticity is quadratic in signs. At four vertices, a fourth-order correlation becomes essential. The same threshold appears independently in rooted graph costs.
- **Circle compatibility:** simultaneous measurable sign antiperiodicity is equivalent to a common odd-frequency condition. With rotation pi/3, a rational second angle is compatible precisely in the stated odd/odd case. Irrational incompatibility can coexist with arbitrarily small defect, so qualitative impossibility does not supply a positive uniform gap.
- **PSD versus Boolean cuts:** an explicit radial positive-definite kernel on the q=17 norm plane predicts cut value 1734, violating the Boolean bound 1730. Spectral feasibility therefore cannot replace realizability by a binary law.
- **Finite norm spectra:** Kloosterman multipliers preserve the shared dual norm. Exact interval arithmetic gives the q=23 positivity bound; the classical Weil estimate gives the stated large-q regime.
- **Functional graph versus algebraic linearity:** every functional graph, even for nonlinear dynamics, is a pseudoforest after loops and duplicate edges are removed. Its finite bipartization cost is at most one third of its vertices. This easy special case distinguishes graph-linear hypotheses from linear algebra.

## The unifying proposal, with a testable scope

The most productive common structure is not a coincidence among small integers. It is the tension between a finite observable and the extra data needed to lift it:

    observable cycles → primary jets;
    pair moments → higher even moments;
    finite precision → invariant lattice filtration;
    finite winding quotients → profinite completion;
    finite colorings → regularity of a limiting choice.

A precise research program is to identify the smallest enriched observable that makes each proposed lifting map faithful. This session supplies successful enrichments (Hasse jets, even-star tomography, prefix-Smith observability), plus sharp counterexamples to several cheaper replacements. This is a research direction, not a theorem asserting that all the major conjectures are equivalent.

## Next substantial proof gates

1. Determine whether the q=17 lower bound four is attainable, or derive a structural inequality improving it. The existing certificate settles only exclusion through two and parity.
2. Find graph embeddings or gadget reductions that genuinely preserve the relevant geometric constraints before transporting finite-field forcing to a real-plane problem. The naive full norm-graph lift is impossible.
3. Classify root-gadget families whose entire boundary cost has degree at most two. The “each interior component touches at most three roots” condition is sufficient; necessity requires a more refined classification.
4. Extend the prefix-Smith observability criterion to other recovered global recurrences, separating sequence orbit period, state-module size, and minimal monic recurrence degree.
5. For Collatz, an affine-residue clock still needs a mechanism retaining carries and controlling integer descent. For planar Jacobian, finite clock/trace analogies still need the geometric properness or infinity obstruction. None of the present identifications supplies those missing steps.

## Reproducibility and source boundaries

The supplied blueprints were read in full and used as attributed sources; their claimed formal developments were not rebuilt here. Full attachment texts are not duplicated in this packet. See [source ledger](SOURCE_LEDGER.md). Each result packet retains its scripts, exact certificates or receipts, and independent review. Certificate verification is finite evidence for the stated finite theorem; global theorems are accompanied by arguments, not inferred from a fitted finite prefix.

All session changes are intended as ordinary commits on main. Publication status and final commit identity are recorded separately by the publication ledger; links here describe the completed research artifacts and must not be read as evidence of a successful push before that ledger is verified.
