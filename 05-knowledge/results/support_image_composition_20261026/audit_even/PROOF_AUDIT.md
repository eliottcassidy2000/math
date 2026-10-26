# Independent audit of signed root cards and the odd-component defect

2026-10-26. **PASS, with no mathematical correction requested.** This review
covers E1–E4 in [the producer note](../even_graphs/RESULT.md), its stated
failure boundaries and the distinction between an unrestricted graph-state
augmentation formula and a recurrence in the preceding scalar A000568
values. It does not assert literature priority or a tournament/graph bijection.

The reviewed proof-candidate SHA256 is
`ffb4c3e17563d4c7c51941eafb8bced11cd63700781675c1813df5fa73b0661c`.
Reviewed source and artifact hashes are bound in [MANIFEST.json](MANIFEST.json).
A subsequent status/link-only promotion is covered if its mathematical body
is unchanged; that comparison and the promoted hash will be recorded there.

## Mathematical review

**E1 passes in both directions.** The marked neighbourhood condition is
exactly what extends a deletion-card isomorphism to a child automorphism.
Deleting the distinguished root is its inverse. Hence the card fibres,
over all reachable roots, partition the entire child automorphism group.
Dividing its character sum by the root stabilizer order gives the orbit
size for an even child and zero for an odd child. Parent-neighbourhood
orbits classify rooted child types, so the final sum ranges once over each
vertex orbit of a fixed unrooted child. The orbit sizes sum to the child
order, explaining the factor `1/N`. No automorphism-free or transitive
action assumption is needed. Empty parents and orders one and two are valid.

The root is last in the producer's reference ordering, so its stabilizer
reverses no root-incident edge. Its character therefore restricts to the
parent character. This validates the relative-even filter; the filter does
not substitute relative evenness for global evenness.

**E2 passes, including its stronger transporter statements.** Nonempty
card fibres are cosets of the root stabilizer. If that stabilizer has an
odd element, each fibre cancels internally. Otherwise its coset sign is
well-defined. For a globally odd graph the index-two sign kernel has two
equal orbits on the root orbit, because the root stabilizer lies in the
kernel. Thus exactly half the reachable roots have each sign. This argument
would fail without the relative-even hypothesis, which the producer retains.

For any three vertices of such an orbit, composing two transporters gives
the third transporter sign. The resulting matrix is `s*s^T`; changing the
base root only changes the common sign of `s`. This proves rank one,
`B^2=|O|B`, and the zero row sum in the odd case. None of these matrices is
being identified with a tournament orientation.

**E3 passes, including repeated components.** Orient corresponding copies
of a connected component consistently. A permutation of entire copies then
has positive edge sign; internal component signs multiply. Any odd symmetry
of a component away from the root remains an odd root-fixing symmetry of
the full graph. Conversely, when the off-root components and the root
component's stabilizer are even, every root-fixing symmetry is even. This
is an iff, not merely a necessary test.

If the graph is odd but the root is relatively even, exactly one odd
component exists and it contains the root. A second isomorphic copy would
be another odd component away from the root and is excluded. Consequently
the root cannot move to a different component, and all off-root factors
have positive sign. Both the orbit equality and sign localization follow.

**E4 passes as an ordinary generating-function identity.** Taking the root
component and the remaining even graph is a unique decomposition of
unlabelled rooted types. Its inverse is disjoint union. The even remainder
cannot contain a component isomorphic to the odd rooted component, so
there is no omitted exchange factor or factorial. The Cauchy product
`D(z)=O(z)E(z)` therefore counts exactly the declared types, including the
zero conventions at order zero.

## Dependencies and scope

The primary [Royle–Praeger–Glasby–Freedman–Devillers v2](https://arxiv.org/pdf/2204.01947v2)
was read directly for Theorem 1.1 (printed page 3), Theorem 3.3 (page 6),
Lemma 4.1 and Corollary 4.2 (page 7). They support unlabelled equinumerosity,
the classical Burnside sum and the edge-sign character. Their main theorem
is stated for `n>2`; the producer supplies the remaining tiny cases directly.
The reference does not provide a natural local attachment bijection.

I also read the inherited signed-root augmentation proof and its root review.
Its earlier safe recurrence already admits all parents and corrects root
forgetting. The new card formula concretely supplies its missing root-moving
sign information. It is correctly presented as a refinement. The connected
defect factorization is proved here rather than inferred from count equality.

The classical odd-cycle partition formula has the correct denominator and
exponent: an odd vertex cycle contributes `(length-1)/2` pair orbits, and
two cycles contribute their greatest common divisor. The pointed identity
uses vertex **orbits**, not `N` times an unrooted count.

The child has already been explicitly constructed from its parent and
neighbourhood. Its marked cards can be compared using only order-`n`
isomorphism searches, but this effectively recovers child automorphisms.
Thus the statement “no order-`N` automorphism oracle is needed” is correct
as an implementation claim; it is not a complexity improvement. Keeping
all graph parents is essential, and the note explicitly preserves this cost.

## Independent exact checks

[check.py](check.py) imports no repository or producer code. It starts from
the complete labelled graph universe, removes relabelling orbits, and then
works **root first** on each child. It constructs induced deletion-card
maps on original vertex names and applies the neighbourhood test separately.
This reverses the producer's parent-neighbourhood augmentation traversal.
It then compares the resulting maps against independently enumerated full
child automorphisms and checks that the rooted types exhaust exactly the
marked parent-orbit keys. Graphs use explicit edge sets; the edge character
is cross-checked by an exact rational product of integer differences.

The universe includes all 209 graph types through order six, all 663 rooted
child types of orders one through six, all 6,169 card transport maps for those
root representatives, and an independent enumeration of all 77 tournament
types. It checks the charge of every rooted child, every relatively-even
transport matrix, the component iff, actual restricted transporter signs,
the defect convolution and the pointed identity. It reproduces

```
even graphs:              1,1,1,2,4,12,56
connected rooted defect:  0,0,1,0,4,9,58
all rooted defect:        0,0,1,1,5,15,75
pointed tournaments:      0,1,2,4,12,48,296.
```

The even-parent-only order-three count is exactly `4/3`, confirming the
claimed hostile rather than silently imposing the erroneous parent filter.

Normal and optimized runs are byte-identical: **25,468 explicit checks,
PASS**. A separate producer replay matches its frozen 17,990-check receipt.

```bash
python3 -B check.py > receipt.json
python3 -B -O check.py > receipt_optimized.json
cmp receipt.json receipt_optimized.json
```

The all-order verdict rests on the proofs above; these exhaustive finite
checks corroborate the implementation and the equality/failure boundaries.
