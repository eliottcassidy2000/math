# Independent audit of the forest-cone dichotomy

**Verdict: PASS, with the stated boundaries; not an unrestricted circuit lower-bound advance.**

Reviewed `circuits/FOREST_CONE_DICHOTOMY.md`. The checker here imports no producer code. This audit separates an analytic proof review from bounded verification; enumeration does not prove the unbounded theorem.

## Proof review

1. **Leaf identity and matching.** For a leaf x adjacent to y, the identity
   `Q_F = y(x + sum_{u in N(y)-x}u) + Q_{F-{x,y}}`
   is exact in the Boolean polynomial ring. A maximum matching contains xy after exchanging y's matched edge, if necessary. If y is unmatched, x is also unmatched and adding xy contradicts maximal cardinality. Hence the residual matching number is exactly nu-1. Isolated residual vertices must be ignored but retained as dimensions. All pivot-neighbor degrees sum to e, since each edge disappears once.

2. **Original-coordinate cost.** Computing each parenthesized sum costs deg(y)-1 XOR gates, including zero for an isolated edge. There are e-nu transform XORs, nu ANDs, and nu+delta-1 combining XORs. All transform operands are original surviving boundary wires, so this is a literal acyclic circuit; no inverse transform or free affine gate is hidden. External users of boundary wires are unaffected. Deleting old products and XORs is valid only under the explicit exclusivity hypothesis; root fanout can be redirected as a whole. The bound is an upper bound for a general forest, not asserted optimal there.

3. **Owned M/D4 realization.** Process pivot chains in their construction order. At a current pivot, the leaf coordinate has exactly its first chain XOR as child; it has no outside equations or consumers by ownership. M moves that gate into an equation. Exchanging this coordinate with its output makes the old coordinate occur only in that equation, so D4 removes it. The fresh coordinate is the input to the next chain XOR. Earlier chains have already disappeared, so they leave no computational uses of later leaf coordinates. Their information persists only in the affine embedding, where it is not computational charge. The new equation is a 1-by-1 invertible block; the *whole ambient* auxiliary matrix need not be 1-by-1, but its old invertible block is unaffected because all forest coordinates are absent from old equations. Each step is an invertible shear and removes exactly one charged equation. This justifies the additional e-nu saving without treating GL transforms as free in ordinary circuits.

4. **Usage accounting.** Initially all v graph parameters occur in products, and all optional z usage remains. Finally exactly 2nu graph coordinates occur in the nu disjoint canonical products. Thus Delta I=-(v-2nu), provided the stated ownership holds. Keeping the unused coordinates in the embedding preserves m, q and a=m-q; no dimension deletion has been smuggled into I. An optional outside use of z does not change this accounting because z is unchanged.

5. **Affine optimality.** Leaf shears split off nu hyperbolic planes, giving polar rank 2nu. A k-coordinate representation has rank at most k. If delta=1, a radical direction has derivative one; in a representation with exactly 2nu essential coordinates the induced polar form there is nonsingular, forcing every radical vector to project to zero, contradiction. Therefore at least 2nu+delta coordinates remain essential in every invertible affine chart. Any fan-in-two single-output DAG depending on k distinct input coordinates needs at least k-1 computational gates, even with sharing and arbitrary binary gates. The canonical construction attains this. The assertion is about computational circuits after changing coordinates, not a lower bound on arbitrary charged equation/gate presentations where coordinate dependencies require separate accounting.

6. **Degenerate cases.** The theorem explicitly assumes e>0 and no isolated vertices; therefore nu>=1 and all displayed gate counts are nonnegative. A matching gives e=nu and v=2nu: both savings are zero. A single edge has one gate (delta=0) or two (delta=1). Empty forests are outside scope: constant and lone-literal terminals have zero cost, and the negative constant-case formula must not be extended to them. Isolated dimensions can be appended without changing gate counts, but cannot be included in the initially-used v count.

7. **Domains and limitations.** All forest parameters are free/unpaired. Arbitrary shears involving marked-pair coordinates need not preserve the prescribed product domain. Shared product outputs invalidate the deletion inventory. Dependent boundary signals invalidate the independent-input optimality claim but not the fixed-coordinate algebraic upper bound. Minimality of an initial circuit cannot be assumed for later restricted states. Nothing here establishes universal exposure of forest cones in elimination sequences or improves the general B2 record.

### Minor wording clarification

The boundary says all product outputs have no outside consumers, and also allows arbitrary root consumers. For the one-edge, delta=0 case the root itself is a product. Read the root exception as overriding the first clause (or explicitly say “all nonroot product outputs”). This does not affect any proof or count.

## Independent bounded evidence

`check_independent.py` uses component-count acyclicity, brute-force edge-subset matchings, binary linear-span cardinality, literal DAG evaluation with packed truth tables, forward and inverse shears, and endpoint essential-input counts. Its pivot construction uses the highest leaf, unlike the producer's lowest leaf. It also tests the matching recurrence for **every** legal first leaf in every enumerated forest.

- 1,819 labeled simple forests without isolates on 2–6 vertices.
- 330,948 input assignments across delta=0 and delta=1, tested through packed exact truth tables.
- 5,716 legal first-leaf choices checked for matching recurrence.
- All 20,160 GL(4,2) charts for each of four examples: four-vertex path, four-vertex star, two disjoint edges, and a three-vertex path plus an independent radical-linear input.
- Minimum essential-coordinate counts respectively 4, 2, 4, 3, agreeing with 2nu+delta. Affine translations need not be separately enumerated: translating a function never changes which coordinates are essential.
- Normal and `python -O` receipts agree byte-for-byte; guards use explicit exceptions, not assertions.

This checker does not enumerate all arbitrary B2 circuits, mechanically validate the external elimination framework, or certify a new complexity record. The matching and gate-count claims for unbounded forests rest on the reviewed argument above.
