# Arbitrarily large bounded-edit divisibility barriers

**Status:** PROVED orbit-counting lemma; FINITE-EXACT zero-response isolation for the named substitution controls. This is not a global-spectrum obstruction or a proof about arbitrary zero-response walks.

This is a continuation of Theorem D in `RESULT.md` and the deletion-deck obstruction in `PHASE2.md`. It reuses independent triangle rotations, but retains the **boundary vertices of a Hamiltonian path**, not just one fixed exposed gap.

## 1. A marked-boundary theorem

Let Q be any tournament on k vertices. Replace each quotient vertex by a directed triangle and orient every pair of blocks according to Q:

    T=Q[C3,...,C3],                  |T|=3k.

Choose a set R of r vertices. Let T' be any tournament on the same labelled vertex set obtained by changing arcs only when at least one endpoint belongs to R. No strongness or preservation of modules in T' is assumed. Set

    e=max(k-3r,0).

### Theorem

    3^e divides H(T'),
    2*3^e divides H(T')-H(T).                            (1)

### Proof by cells of marked path boundaries

Partition the directed Hamiltonian paths of T' by the following finite data:

1. the position of every labelled vertex in R in the path word;
2. the identity of every old vertex immediately preceding or following one of those positions.

There are at most 2r such old boundary vertices, including possible coincidences. At most 3r original triangle blocks meet R or one of its boundary vertices. Consequently at least k-3r blocks avoid all marked vertices.

Rotate any of these unmarked triangle blocks independently. These rotations fix every vertex of R and every old boundary vertex. They preserve all old-old arcs, because `T'-R=T-R` and an untouched triangle block has the same original external orientation to every other old block. An adjacency of the path involving R is also preserved: its old endpoint, if any, is one of the fixed boundary vertices. Thus the product of the unmarked rotations acts within each path cell.

This action is free on linear path words. A nonidentity permutation of vertices cannot fix an ordering in which every vertex appears exactly once. Each cell cardinality is therefore divisible by `3^(k-3r)` when this exponent is positive. Summing cells proves the first assertion. The same rotation action on T itself gives `3^k | H(T)`, so `3^e` divides the difference. Both tournament path counts are odd by Rédei parity; coprimality of 2 and 3 proves the second assertion. QED.

**The cell partition matters.** The altered tournament T' may have trivial automorphism group. We do not assert that the rotations are global automorphisms of T': they only need to preserve the used adjacencies inside each path cell. Forgetting the marked neighbors would invalidate that step.

The argument generalizes verbatim if each block admits a chosen automorphism group of the same odd order g: replace `3^e` by `g^e`. The C3 formulation needs no group-classification input.

## 2. Consequences and the zero-walk boundary

If `k>=3r+1`, then no r-vertex replacement can change H by `+2` or `-2`. More generally, any target T' with `3` not dividing `H(T')` must satisfy

    r >= ceil(k/3)                                      (2)

for every vertex set covering all arcs on which it differs from T.

For q individual arc reversals there is a simpler, stronger bound. At most 2q blocks meet their endpoints; all rotations of the remaining blocks are actual automorphisms of T'. Hence

    3^max(k-2q,0) divides H(T'),
    2*3^max(k-2q,0) divides H(T')-H(T).                  (3)

Thus a nonmultiple-of-three target requires `q>=ceil(k/2)` arc edits.

When Q is strong, T is strong. Equations (1)--(3) therefore give arbitrarily large local barriers inside the starting strong-tournament universe. They do not assert that a desired adjacent H target exists at that same order.

Most importantly, the theorem does **not** prove that a whole `0,+/-2` move component is trapped. A zero-response edit might leave the symmetric family and later allow a unit change. That missing zero-response closure is exactly what the explicit H147 leaf in `PHASE2.md` verifies separately.

## 3. Finite zero-response closure for the smallest substitution families

The script exhausts every strongly completing cut of every deletion-card orbit of `Q[C3^k]`, for **all** strong quotient isomorphism classes at k=3,4,5. Independent rotation of each original C3 block makes its three deletion choices equivalent, so one deleted vertex per block is complete.

Across these eight quotient classes and 37 deletion-orbit rows, every parent has exactly one cut restoring the original H, namely the original attachment cut. No cut gives H+2 or H-2. Consequently all these particular labelled substitution tournaments are isolated under arbitrary one-vertex-replacement unit moves, even before passing to isomorphism classes.

The first controls are:

| quotient order | quotient code | H(Q) | order of substitution | H(T) | every deletion response gcd |
|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 3 | 9 | 3159 | 2 |
| 4 | 4 | 5 | 12 | 241137 | 6 |
| 5 | 8 | 9 | 15 | 15731091 | 18 |

The other five strong order-five quotient rows are retained in the output. At k=3 the gcd is saturated, so its isolation is a finite feasible-state phenomenon rather than an arithmetic consequence. At k>=4 the theorem excludes the +/-2 levels; the finite unique-zero check supplies the additional closure needed for actual isolation.

**OPEN:** unique zero-response reconstruction for every strong quotient Q. The eight finite rows do not prove that assertion. It is equally possible that a larger quotient has a nonisomorphic zero-response replacement and escapes the family. No all-order isolated-component claim is made.

### Extension through every strong quotient of order six

The independent C++ compiler `zero_escape.cpp` extends the finite unique-zero test to all **35** strong order-six quotient classes. Its 210 deletion-card orbits represent order-seventeen parents of order-eighteen children. Every one again has exactly the original cut at H(T).

This extension is exact without enumerating every labelled cut separately. An intact C3 block has four subset orbits, indexed by cardinality 0,1,2,3 with multiplicities 1,3,3,1; the two vertices of the broken block have four literal patterns. Thus each row needs `4^6=4096` profiles representing `2^17=131072` cuts. Totals are

    35 quotient classes,
    210 deletion-orbit rows,
    860160 orbit profiles,
    27525120 represented labelled cuts,
    zero nonoriginal H-restoring cuts.

The compiler constructs the endpoint tables and exposed-gap capacities from scratch, directly counts each full order-eighteen original tournament, and checks that every deletion row reproduces that original count. No unsigned overflow is possible in this k<=6 range: all path and exposed-word counts, and their nonnegative convolution partial sums, are bounded by the corresponding permutation universes, at most 18!, which is below 2^63. The quadratic evaluator's linear coefficient total is also bounded by the complete exposed-word sum. The source is intended only for the documented k<=7 audit range, not arbitrary k.

The separate k=7 run has now completed with its own frozen `zero_escape7.out`:

    353 strong quotient classes,
    2471 deletion-card orbits,
    40484864 orbit profiles,
    2591031296 represented labelled cuts,
    zero nonoriginal H-restoring cuts.

This is producer FINITE-EXACT, not an independently replicated full census. A selected order-twenty parent is additionally checked by counting all212 literal order-twenty-one children needed to recover every quadratic coefficient; the independent parent-appointed auditor also checked one complete order-seventeen parent by all131072 literal cuts, without orbit compression. These are cross-route controls, not a relabelling of the full k=7 run as independently replicated.

The universal divisor need not equal the exact gcd. At k=6 all210 rows have gcd54. At k=7,2470 rows have gcd162 and exactly one has gcd486: quotientcode755846 with H(Q)=99, rootblock3. Its quotient masks are `(12,9,10,112,71,23,39)`, with homogeneous triples `{0,1,2}` and `{4,5,6}`; it is `C3[C3,K1,C3]`. Thus the exceptional deletion parent has the shape `C3[U9,TT2,U9]`, where `U9=C3[C3,C3,C3]`. This refutes equality in the generic divisor bound within the checked family. No general explanation of that extra factor3 is asserted.

All strong quotients through order seven therefore pass the scalar zero-rigidity test. The theorem for all quotient orders remains OPEN. `TRIANGLE_CARD_RIGIDITY.md` removes a separate ambiguity: any nonoriginal zero cut would necessarily leave the entire C3-substitution family, not merely change its presentation.
 For its arithmetic range, a universal bound `H(T)<=n!/2^floor(n/2)` follows by requiring only the disjoint consecutive-position arcs `(1,2),(3,4),...` in a permutation. At n=21 this is below 2^63. At parent order20, every Q entry is at most19!, and the full exposed sum is at most `19*H(B)+20!`, again below2^63, bounding every capacity and quadratic intermediate. Thus the documented k=7 run also fits signed64-bit arithmetic.

### Why the coarse invariants do not prove zero closure

Strongness is genuinely relevant: for the nonstrong quotient of order two, `C3=>C3` has H=9 and a deletion parent with masks `(30,28,8,16,4)`. Four distinct cuts `14,22,26,29` restore H=9. This is an actual tournament countercontrol to dropping the strong-quotient hypothesis, not to the strong candidate.

There is a stronger **formal-network** countercontrol to a tempting proof route. Take the actual order-eleven parent obtained by deleting vertex0 from the order-twelve control. Its masks are

    (1794,1792,11,19,7,95,159,63,764,1276,508).

It has H=24651, and its actual ear at cut1793 has H=241137. Starting with its actual capacities, change only these symmetric edge orbits in local parent labels:

    c_(0,j): 7182 -> 6684 for j=2,3,4,
    c_(0,j): 5616 -> 6114 for j=5,6,7.

The changes are -498 and +498 on three edges each. They preserve the C3 symmetries, every capacity's positivity, every actual pair-deletion lower bound, exact ear gcd6, total capacity582822, **every cardinality-layer ear mean**, parent H, and the original child value at cut1793. Yet the formal cut response has four states at241137:

    480,736,1248,1793,

where the actual network has only1793. In particular, the actual value at480 is239643; the formal perturbation raises it by1494.

The global maximum of the actual pair-deletion lower bounds is4194, below both altered capacities, so positivity and the deletion floor have not been silently violated. The compensating change preserves the total capacity and hence all fixed-cardinality means; the original cut uses neither changed directed-edge family.

This is **not** a tournament counterexample. It deliberately ceases to be the path-capacity kernel of the fixed parent. It shows exactly why symmetry, positive capacities, deletion floors, congruence, and even all layer means cannot alone prove zero closure. A successful proof must use an additional realizability identity from the actual exposed-path / rooted-pair tensor, or a different argument that distinguishes actual kernels from these formal perturbations.

The current exact reformulation is a quadratic equation on the box

    {0,1}^2 x {0,1,2,3}^{k-1},

with coefficients supplied by a genuine deleted-triangle path kernel. The observed original point is unique in the audited strong quotient universe. Generic positive, symmetric, congruence-compatible coefficients do not have that uniqueness property.

## 4. Replay and provenance

From this directory:

    python3 -B symmetry_edit.py > symmetry_edit.out
    python3 -B -O symmetry_edit.py > symmetry_edit_optimized.out
    cmp symmetry_edit.out symmetry_edit_optimized.out
    g++ -O3 -std=c++17 -Wall -Wextra -Werror zero_escape.cpp -o zero_escape
    ./zero_escape 6 quotients6.txt > zero_escape6.out
    ./zero_escape 7 quotients7.txt > zero_escape7.out
    python3 -B zero_closure_controls.py > zero_closure_controls.out

The all-order proof needs only independent block automorphisms, a free action on linear orders, and Rédei parity. Finite tests use the inherited orientation-supported ear polynomial, with the local independent endpoint-Q compiler and direct child DP controls. They do not import inherited research code. `phase2_order5.json` supplies the complete six-class strong order-five universe already generated and audited in this packet.
