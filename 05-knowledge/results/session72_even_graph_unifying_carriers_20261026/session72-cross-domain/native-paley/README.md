# A native Paley7 three-port obstruction to pairwise alignment

**Status: PROVED elementary finite theorem; FINITE-EXACT independent controls.** This supplement leaves the frozen cross-domain core unchanged. The parent proposed the explicit pair below; this lane supplies direct proof, a native distinguishing guard, and exhaustive adjacency verification. Session69 already identified the full group and the two weight-three Fano families. The contribution here is their use as a native relative-alignment counterexample, not discovery of Paley7 or those counts.

## The theorem

Let T be the Paley7 tournament on Z/7 with

    i -> j iff j-i belongs to {1,2,4}.

Consider the two ordered port triples

    x=(0,1,3),     y=(0,1,5).

**P7-ALIGN.**

1. For every pair of positions i<j, the full ambient two-pointed tournaments `(T;x_i,x_j)` and `(T;y_i,y_j)` are isomorphic. This is stronger than merely having identical induced two-vertex arcs.
2. The three-pointed tournaments `(T;x)` and `(T;y)` are not isomorphic.
3. An actual native three-port guard distinguishes them:

       #{v: all three marked ports point to v} = 0 at x, 1 at y.

   The unique common out-neighbour of y is2. Conversely x has common in-neighbour6, while y has none.

Thus separate complete pair-frame types cannot always be glued into one coherent triple frame. This occurs in a genuine tournament with an odd-order automorphism group, not only in the A4 toy action.

## Proof from the native adjacency

All affine maps `j -> a*j+b`, with `a in {1,2,4}` and `b in Z/7`, preserve the difference predicate. They form a group of order21.

There are no further automorphisms. The out-neighbours of0 are `{1,2,4}` and induce a directed triangle. An automorphism fixing0 has at most three possible restrictions there. If it fixes those three vertices individually, the remaining in-neighbours `{3,5,6}` are distinguished by which one of `{1,2,4}` each dominates: respectively4,2,1. It therefore fixes every vertex. Hence the point stabilizer has size at most3 and the full group at most7*3=21; the displayed affine group attains this bound.

The affine group is transitive, in fact regular, on the21 directed arcs: the image of ordered arc(0,1) uniquely determines a and b. It likewise has one orbit on the21 reverse-arc ordered pairs. Both x and y have the same cyclic pattern `port1 -> port2 -> port3 -> port1`; each pair of corresponding ports consequently lies in one common ordered-pair orbit. Explicit affine witnesses are available in the receipt.

But an affine automorphism fixing0 and1 has b=0 and a=1, hence is identity. It cannot send3 to5. This proves the triple obstruction. Direct adjacency gives

    N+(0) intersect N+(1) intersect N+(3) = empty,
    N+(0) intersect N+(1) intersect N+(5) = {2},
    N-(0) intersect N-(1) intersect N-(3) = {6},
    N-(0) intersect N-(1) intersect N-(5) = empty.

The guard counts are invariant under a simultaneous ambient automorphism, proving the same obstruction through an actual extension predicate. QED.

## Exact relation to the cyclic alignment theorem

CA1 and CA2 concern one **cyclic** group acting by simultaneous translation on its cyclic quotient factors. They remain true. Here the permitted alignment group is the larger noncyclic group `C7 semidirect C3` of order21. In F8 exponent coordinates it is `GammaL1(8)`, including both scalar multiplication and Frobenius; it is not the pure cyclic Frobenius action of CA1. Changing the acting group changes the phase problem.

The general carrier is still the diagonal quotient of the three pointed-state orbits, but its projections to every pair quotient lose information. The loss is exactly visible in a native ternary adjacency conjunction. This is deterministic alignment loss, distinct from equality of pairwise probability marginals on the Boolean laws in Sessions70/71.

## The inherited Fano families, recovered as the obstruction

P7 has14 unordered cyclic triangles. The two common-out-neighbour counts0 and1 split them into two automorphism orbits of size7. Each orbit is a Steiner triple system on the same seven vertices: it has7 triples, and every unordered vertex pair occurs exactly once. Indeed the group is transitive on the21 unordered pairs, while7 triples provide21 pair incidences, so their uniform positive incidence is one. These are the two Fano families already present in Session69; the new purpose is that their ordered triples have indistinguishable complete pair-frame types.

For orientation-sensitive finite bookkeeping, there are210 ordered distinct triples and10 full-group orbits, each of size21. The six transitive pair-orientation patterns each determine one triple orbit; each of the two cyclic patterns permits two. The independent checker establishes these bounded counts by full permutation enumeration, not by a presumed group formula.

## The exact F8 dependence detector

The parent supplied a further concrete interpretation, verified here directly. Use `F8=F2[alpha]/(alpha^3+alpha+1)` and label vertex j by `alpha^j`. Its seven polynomial-basis bit masks are `(1,2,4,3,6,7,5)`. The actual tournament predicate is

    i -> j iff Tr(alpha^(j-i))=0, i != j.

Thus this is the Session69 field dictionary, not an arbitrary bijection. Under it,

    (0,1,3) -> (1,alpha,1+alpha),
    (0,1,5) -> (1,alpha,1+alpha+alpha^2).

The first triple has F2-rank2 and sum0; the second has rank3 and sum alpha^2. Every pair in either triple has rank2. The full automorphism group `GammaL1(8)` acts F2-linearly by `z -> beta*z^(2^k)`, so rank is a genuine full-group invariant. We do **not** replace this group by GL3(F2), whose extra transformations need not preserve the tournament.

Within the cyclic-triangle locus, the common-out count is exactly F2-rank minus2: the seven additive-line triangles have no common out-neighbour; the other seven cyclic triangles have one. This statement is scoped to cyclic triangles; rank3 by itself does not force a triple to be cyclic. Hence the lost joint coordinate is a literal ternary linear dependence, observed natively by a three-port adjacency conjunction. No field-addition law for arbitrary tournament composition is asserted.

## Transfer contract and reproduction

- **Source:** P7 with three named distinct vertex ports.
- **Operation/symmetry:** one simultaneous ambient tournament automorphism.
- **Coarse map:** the three separate full ambient ordered-pair isomorphism classes.
- **Preserved:** each named two-port type and its orientation.
- **Lost:** coherent three-port alignment; concretely the orbit of the cyclic triangle.
- **Sidecar/terminal:** common-out-neighbour guard, or the full ordered triple orbit.
- **Cheapest decisive test:** the explicit x/y pair above.

Run `python3 -B check_native_paley.py` and `python3 -O -B check_native_paley.py`; compare to the paired receipts. Both modes enumerate all5,040 vertex permutations, reconstruct all21 automorphisms and all210 triple states, check the14 cyclic triples and every pair incidence in each family, and verify the native guard directly. `--inject-failure` must fail in both modes. Source pin verification is mandatory. This supplement imports no producer code and does not mutate the earlier frozen manifest.
