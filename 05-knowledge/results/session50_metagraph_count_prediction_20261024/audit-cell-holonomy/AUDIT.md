# Independent audit: equitable-cell errors and corrected order-eight prediction

24 October 2026. **PASS**, subject to the provenance correction and exact scope
below. The new mathematical consequence is a structural correction to a
refinement upper trace. The signed affine-cell carrier and its automorphism
formula were already proved in Session39.

## 1. Corrected inheritance

The initial producer draft credited the Session40 two-prime matching family
but missed a closer source. Independent review identified Session40
`marked-ports/PATH_TRANSPORT_AND_FOREST_RULE.md`, which already proves coprime
crossing rigidity, uses odd stable cells in its small-core argument, and
explicitly invokes the triple-cell forest and joint affine-cycle carrier.
The producer then recovered the non-materialized Session39
`marked-cells/AFFINE_CONNECTIONS.md` from its Git object. It explicitly already
contains:

- the affine F3 crossing action, not merely its two-cell special case;
- the balanced/unbalanced signed kernel and rotational automorphism order;
- the fixed-sign phase-orbit count;
- the 216 three-cell presentations giving sixteen ordinary classes, twelve
  with automorphism group C3 and four rigid.

These mechanisms and finite family counts must be attributed as inherited.
The all-r complete-support formula is a specialization and sum of the
inherited fixed-sign orbit formula, not a discovery of affine holonomy.
The new extracted consequences audited here are the refinement-error weight,
its ordered-SCC transfer, strong-order-eight vanishing, and the resulting
child-refinement-only order-eight prediction and higher-order lower bounds.
No public literature-priority claim is reviewed or established.

Historical source bytes and their supplied commit/path pins are copied under
`source/` and recorded in `historical-sources.txt`. This audit does not modify
or run Git in the producer or publication worktree.

## 2. Independent mathematical proof review

### Stable colors, not anonymous blocks

Starting from one color and retaining the old color at every refinement step
makes each final color class individually invariant under every automorphism.
An arbitrary equivariant unlabelled block system would not suffice. A final
cell induces a regular tournament by equitability, hence has odd size.
For cells of sizes a,b, counting the same directed cross arcs gives ar=bs;
if gcd(a,b)=1 this forces r=0 or b. The crossing is uniform.

### Signed compatibility gives the exact error

A stable three-cell is a directed C3, whose orientation-preserving
self-isomorphisms are its three rotations. Nonuniform equitable 3-by-3 blocks
are matchings or their complements. Write the normalized matching as
sigma(x)=epsilon*x+c over F3. A rotation k_i at one end and k_j at the other
preserves the block exactly when k_j=epsilon*k_i. Translation c cancels.

The automorphism kernel therefore has one free F3 coordinate on each balanced
signed component and none on an unbalanced one. Stable colors prevent any
additional permutation of cells. A triangle in a balanced component is a
single orbit; all three vertices in an unbalanced component are individually
fixed. Thus the refinement excess

    e(T) = sum_vertex_orbits O (size(stable_cell(O)) - size(O))

is exactly six times the number of triangle cells in unbalanced components.
This is an extraction from the inherited automorphism theorem, not a claim
that signs classify the tournament. Phases can distinguish isomorphism types
even when the full automorphism action and error coincide.

### All possible cells through order eight are accounted for

A regular five-vertex tournament is unique and vertex-transitive: choose v,
order its two outneighbors a->b and its two inneighbors c->d, then regularity
forces b->c,d, a->c and d->a. This is the cyclic regular five-tournament.

With no seven-cell, an order-at-most-eight stable partition consists of
singletons and at most two triangles, or of one five-cell and residual cells
of size one or three. The triangle support has no cycle. The five-cell is
uniformly attached to all other cells by coprimality, so all its internal
rotations extend independently. Consequently all these cells are exact
vertex orbits. At order seven the only remaining case is a regular
seven-cell. At order eight a seven-cell leaves a singleton with uniform
crossing, which prevents strong connectivity. This proves zero error through
six, for nonregular order seven, and for every strong order-eight tournament
without any order-eight census.

### Ordered SCC transfer is exact

If the SCCs in dominance order are C_i and L_i is the number of vertices in
later components, each score in C_i lies in [L_i,L_i+|C_i|-1]. These integer
intervals are pairwise disjoint. Degree refinement therefore separates SCCs
immediately; external counts are constant on a component, so subsequent
refinement restricted to it is exactly standalone refinement. Automorphisms
fix the ordered SCCs and are the direct product of their internal groups.
Hence e is additive under ordered SCC composition.

Writing A(z) for all ordinary tournament classes, including the empty one,
and E_s(z) for the error-weighted strong classes, marking the component
responsible for an additive contribution gives

    E(z) = A(z)^2 E_s(z).

There are no missing multiplicity or symmetry factors: SCC order is intrinsic,
and the left and right component sequences form arbitrary ordinary classes.
All series are formal ordinary generating series.

### The correction and later lower bounds

The exhaustive regular-seven census has three classes with automorphism
orders 7,3,21. Only the order-three group has multiple vertex orbits, of
sizes 3,3,1, and contributes 14. Thus E_s,7=14. The preceding proof gives
E_s,8=0 and E_s,n=0 below seven. Since a_0=a_1=1, E_8=28. Precisely the two
ordered sums singleton+hostile and hostile+singleton contribute.

For r intrinsically degree-distinguished C3 cells with all forward crossings
permutation matchings, each cell i has degree r+i. Ordinary isomorphism thus
reduces to the inherited rotation action. With m=binom(r,2), there are
2^(r-1) balanced signings, each with 3^(m-r+1) phase classes; all other
signings have 3^(m-r) phase classes. This verifies the stated all-r family
formula. Every such tournament is strong. Its rigid classes contribute 6r
apiece, yielding E_s,9>=72 and E_9>=72+3*14=114.

For Session40's two cyclic prime cells, the p-1 affine ordinary classes have
Aut=Cp and zero error. Every other matching class is rigid, has two stable
p-cells, and contributes 2p(p-1). The inherited family count N_p therefore
gives contribution 2p(p-1)(N_p-p+1), hence E_s,10>=160 and E_s,14>=8568.
The coefficient [z^3]A(z)^2 is 6, so E_10>=160+2*72+6*14=388. These are
lower bounds, not unrestricted higher-order error censuses.

## 3. Independent finite tests

`independent.py` does not import any producer code and uses adjacency matrices
rather than the producer's outgoing-bitset implementation. Every test uses an
explicit exception, so Python optimization does not remove it.

- Complete labelled universes through order six are partitioned by literal
  full-permutation orbits. This gives ordinary class counts 1,1,2,4,12,56,
  and tests every class's actual automorphism orbits, stable cells, SCC
  restriction and zero error.
- An independently implemented vertex-by-vertex degree-constrained generator
  finds all 24 regular-five labels and 2,640 regular-seven labels. Full
  permutation class removal and literal automorphism tests give the unique
  regular-five class and all three regular-seven groups and orbit data.
- All 216 forward-matching three-cell configurations are canonicalized using
  all 216 degree-preserving S3^3 vertex maps, not only rotations. Their full
  literal groups agree with the signed kernel; the quotient has sixteen
  classes split twelve/ four, with total error 72.
- All 46,656 four-cell affine configurations are partitioned by their actual
  translation-gauge orbits. Eight balanced signings each have 27 classes;
  56 unbalanced signings each have nine classes, giving 720 classes. Two
  order-twelve controls also check all 6^4 degree-preserving vertex maps.
- Exact double-rotation orbit enumeration of all permutations at p=3,5,7
  independently verifies the inherited prime-family class and stabilizer
  counts and the stated new error lower bounds.
- Ordered joins with the regular-seven hostile check the SCC mechanism at
  orders eight and nine.

Normal and optimized receipts are byte-identical: 6,830 explicit checks in
the independent suite. This finite evidence
supports, but is not a substitute for, the all-order proofs above.

## 4. Independent end-to-end order-eight prediction

`independent_predictor.py` constructs its own parent atlas, rather than
loading the producer's atlas. It starts from the full-permutation partition
of all 32,768 labelled order-six tournaments. Extending its 56 representatives
by every one of the 64 attachment masks and canonicalizing only at order
seven gives exactly 456 parents. Completeness follows because deleting any
vertex of any order-seven tournament leaves a represented order-six class.
The canonicalizer orders degree cells and exhausts every internal permutation.
It has an explicit runtime prohibition on order-eight input.

For each parent it computes its actual degree-pruned, literal automorphism
permutation action and partitions all 128 attachment subsets into orbits.
For one extension from each orbit it computes only stable color refinement
at order eight and sums the new root's cell size. It does not canonicalize
an order-eight child, delete a child vertex, or test rerooted isomorphisms.

The independently generated result is:

    Parent classes:                  456
    Labelled parent mass:            2,097,152 = 2^21
    Pointed attachment states:       54,256
    Sum of new-root color-cell size: 55,068
    Strong-seven error from parents: 14
    Structural order-eight error:    28
    Corrected trace:                 55,040
    Predicted ordinary classes:      55,040 / 8 = 6,880

The root-cell upper sum is not inferred from the target count. The known
value 6,880 occurs only in an outer post-computation comparison. Neither it
nor an order-eight class list is an input. Normal and optimized output and
receipts agree byte for byte.

The independent `compare_atlas.py` additionally authenticates all 533 frozen
producer records across orders zero through seven against complete independently
generated class sets. It checks every supplied action list equals the entire
literal degree-preserving automorphism group, not merely a valid subgroup.
Both modes agree. The frozen producer itself replays identically in normal and
optimized Python; its predictor receipt hash is
`d19fea022c093103d0ba745ffefa8fadd478e994d2accd940039b0b76145fa85`.
Its AST/input guards constrain the inspected program but are not presented as
a general-purpose proof system for arbitrary programs. The actual arithmetic
and call graph were reviewed directly.

## 5. Scope and stopping boundary

This is a real finite-order simplification over full root transport: order
eight requires a complete order-seven parent atlas, its actual automorphism
actions, attachment subset orbits, and ordinary child refinement, plus an
independently established strong-error theorem. It avoids child classification
and root-change isomorphism tests in this stage. It is not a recurrence from
prior A000568 scalar values or a count from n alone. Constructing the parent
atlas remains classification-strength work.

No measured complexity advantage, full order-eight independent child census,
exact unrestricted order-nine error, or bare blue/black metagraph closure is
claimed. At higher orders the strong-error coefficients are the missing data,
not a dispensable side effect. The native nine-vertex rigid classes and
larger-prime matching families demonstrate why the order-eight zero-error
argument cannot simply be extrapolated.

## 6. Reproduction and closure

Run `python verify.py --deep` and `python -O verify.py --deep` to hash-check
the closed packet and rerun its independent computations and frozen producer
controls from a different temporary working directory. `relocation_test.py`
copies only manifest-listed files before running both modes, then confirms
that extra files, changed payload bytes, and symlinks are rejected. No network,
repository checkout, external atlas, or third-party package is required.
