# Session 244: prescribed even gaps, exact anchors, and root menus

**Session COMPLETE, 10 November 2026. All-target Chan–Pak coverage remains OPEN.**

The requested goal is a proof that every positive odd integer except 7 and 21 is a tournament Hamilton-path count. This packet does not claim that proof. It develops stronger physical gap constructions, tests the missing anchor condition, and separates exact theorems from bounded computational evidence.

## Strongest proved results

1. **Efficient endpoint realization.** Every positive N is the endpoint count at a root of a strong tournament with H≤3N, response U≤4N, trivial root stabilizer, and at most ceil(log2 N)+3 vertices. The complete four-coordinate profile is retained. Endpoint existence was inherited; the quantitative profile is the improvement. [THM-5113](../../../01-canon/theorems/THM-5113-linear-cost-endpoint-universality-from-compact-gap-roots.md), [proof](selector/ROUND3.md).
2. **Pure-count gap synchronization.** Two seven-letter physical words give strong same-order counts 17g(N)+40N and 17g(N)+42N, with identical endpoint, response and deletion coordinates. Their gap is 2N and the upper count is at most 59N. Every common O/Z tail preserves exactly that gap and all three equal sidecars; any common word containing J loses this synchronization. [THM-5115](../../../01-canon/theorems/THM-5115-synchronized-endpoint-gap-pairs-and-exact-repair-rigidity.md), [proof](composition/ROUND3.md).
3. **A sharp structural boundary.** Freely rerooted sole-exit births, sink births and ordinal joins generate precisely the outdegree-1-degenerate tournaments. They are characterized by avoiding six explicit induced cores of orders 5 and 6. More generally an inclusion-minimal minimum-outdegree-k tournament has at most (k+1)(k+2)/2 vertices, and this is sharp for every k. [THM-5116](../../../01-canon/theorems/THM-5116-sharp-outdegree-core-bound-and-six-core-forest-construction-class.md), [proof](logic/ROUND3.md).
4. **Why gap existence alone cannot prove coverage.** Explicit abstract multiplicatively closed odd sets can have every even gap, compact upper bounds, arbitrarily long odd intervals, relative density one among odd integers and constant-error anchors while omitting arbitrarily large primes. These are logical countermodels, not tournaments. [THM-5114](../../../01-canon/theorems/THM-5114-universal-even-gaps-and-constant-error-anchors-do-not-force-odd-coverage.md), [proof](logic/PROOF.md).

A shorter channel additionally gives every gap 2N in strong pairs with response difference exactly one, equal deletion count coprime to three, and upper count at most 24N−13. The endpoint difference is 2N−1, so it is not a pure-count pair. Its exact higher interface exposes a rooted three-path response that the older four-tuple does not retain. [Proof and composition laws](composition/INTERFACES.md).

The pure-count pair's all-tail anchors now have an exact two-index inverse. Even allowing all 32,768 gap inputs N≤32768 and unlimited O/Z tails omits at least **3/94 of odd targets asymptotically**. The full unbounded family remains open. [Formula, prime-input filter and density proof](selector/ROUND4.md).

## Exact target repair and the remaining obligation

Free rerooting repairs the target 173 that is absent from the inherited unit-start rooted grammar. The explicit eight-vertex witness has outgoing masks

    [218,252,57,240,224,193,132,4]

and start-count vector

    [23,72,5,27,17,14,10,5].

The sum is 173; no root has start count one. Thus this graph is outside Session197's carrier, not a contradiction of its complete bounded-image algorithm. Old output calling 173 “NEVER (forbidden)” was bounded at order seven and must not be read globally. [Construction and exact inverse bounds](selector/ROUND2.md), [independent permutation certificate](audit/FINAL_INDEPENDENT.json).

A one-vertex induction needs an actual smaller graph T and eligible root r satisfying

    U(T,r)-E(T,r) = target-H(T).

Universal even differences do not fix H(T); endpoint universality does not fix U(T,r). Matching a desired scalar difference is therefore not enough. This precise root-menu intersection remains open. The inherited Session199 four-block balanced-coverage lemma is another sufficient all-target route and remains open.

## Complete finite forest census and a late counterexample

All 479,001,600 chronological histories on twelve vertices were enumerated, not 479 million isomorphism classes. The graph-class image has 53,549 distinct counts, maximum 115,781, and first allowed missing count 91,081. This is a new exact census of the specified construction class, not an improvement of the known global coverage range.

Depth-two scalar equality survives every order through eleven but **fails at order twelve**: count **97,309** alone requires depth three at that order. A separate normalized shallow enumeration of 23,670,073 histories reproduces all 53,548 shallow counts and independently confirms the single exception. All 138,488 stored positive/minimum-depth witnesses across the twelve orders were literally recounted.

The scope matters: 97,309=31·43·73 has an explicit **depth-one, twenty-vertex** ordinal-join witness. The failure concerns preserving the number of vertices, not an all-order exclusion. A proved adjacent-swap normalization preserves every graph and never increases its reversal depth; it reduces the full forest-history search to (n−1)! Hamilton-path-ordered histories. [Census, proofs, certificates and reproduction](logic/forest_census/release/README.md).

## Further stress tests

- [Binary-ray invariant](dyadic/PROOF.md): a graph can initially support both 2H−1 and 2H+1 extensions but lose the opposite branch after either extension, even after rerooting. A usable induction needs closure of the whole root menu.
- [Forest normal-form obstruction](logic/forest_census/release/NORMAL_FORM.md): a concrete H=217 graph cannot be flattened to depth two at the same order. More strongly, its transition to 373 cannot be supplied by any shallow order-nine graph with H=217. This rules out a scalar-only fixed-order induction independently of the later twelve-vertex scalar counterexample.
- [Saturated-face audit](face/PROOF.md): the outdegree/response defect identity and classification were inherited and independently recovered. Their scoped corollaries are not new universal spectrum results.

## Provenance and verification

All-even-gap existence is inherited from Sessions199/221; the compact strong upper bound 5N is inherited from Session243. Session199 also supplies the zero-J density obstruction, the chronological caterpillar interpretation, joint modular surjectivity, and the still-open four-block lemma. Session231 supplies the nonsink zero-defect classification. New costs, interfaces and failures are distinguished from those inherited statements throughout this packet.

See [independent audit](audit/REPORT.md), [research-round synthesis](SYNTHESIS.md), and [aggregate verification receipt](AGGREGATE.json). Finite images are not promoted to eventual or all-order coverage.

From this packet directory, run:

    python3 -I -B verify_all.py > /tmp/session244-replay.json
    cmp AGGREGATE.json /tmp/session244-replay.json
    sha256sum -c SHA256SUMS

The default aggregate replays twenty-one Python checks in ordinary and optimized modes, strict native controls, every stored forest witness, and the small normalized-search controls. Full exhaustive regeneration is available with `--full-census`; it additionally repeats all 12! histories and the independent normalized shallow census. The full original and independent census receipts are preserved in the release. No unbounded coverage claim is based on those finite runs.
