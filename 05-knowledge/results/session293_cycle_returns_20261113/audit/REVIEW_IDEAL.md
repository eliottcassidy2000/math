# Independent review: exact uniform-module arithmetic ideal

**PASS.** Reviewed control/MODULE_IDEAL.md independently, including the
extra excursion factorial, physical probe basis, arithmetic examples and
prime-power cutoff. This report proves neither positive-semigroup saturation
nor any new Chan–Pak coverage claim.

## Free internal excursions

The context substitution coefficient a_k counts exterior orders together
with k occupied C-slots, not C-path choices. The k−1 internal exterior
excursions between slots are all nonempty, since the C-runs are maximal.
Each begins at an exterior vertex beaten by C and ends at a vertex beating C.
Permuting the entire disjoint labelled excursions retains every traversed
arc: no direct edge between two different excursions is used. Prefix and
suffix stay fixed. Distinct permutations yield distinct skeletons, so the
action is free and `(k−1)! | a_k`. This works with arbitrary exterior arcs;
no transitivity hypothesis has entered the universal divisibility step.

## Unitriangular physical probes

In the alternating transitive exterior of size2m, a maximal m+1-slot pattern
uses every exterior vertex in m internal excursions and has empty prefix and
suffix. Every excursion then has two vertices and must be xi→yj with i<=j.
The matching condition forces j=i by summing inequalities. The m excursions
can occur in all m! orders, establishing top coefficient m! on F_(m+1),
not one. Scaling by the independently established excursion factorial makes
the diagonal exactly one. All lower coefficients are integral.

For |C|=n, F_k(C)=0 beyond k=n. Integer unitriangular inversion of the first
n physical probes therefore recovers each scaled cover count D_k exactly.
Every context is a nonnegative integer combination of these D_k, while
probes are themselves contexts. Both inclusions of the claimed generated
ideals follow. Equality of probe vectors likewise recovers the complete
cover profile; it does not assert that H alone suffices.

## Arithmetic scope

The empty context is valid and forces the gcd to divide odd H(C). Each
module automorphism extends to any uniform context and acts freely on its
Hamilton paths, so its order divides every context value. The two explicit
examples show that this is not the whole obstruction.

The prime-power precision criterion follows directly from the exact D_k
formula. Factorial valuations are monotone, so the finite cutoff is valid.
At prime precision the relevant k are1,...,p−1. These are necessary/sufficient
conditions for the ideal divisor, not for attainable positive counts.

A fixed module with nontrivial ideal cannot be retained uniformly in both
members of a +2 repair, even if other modules are varied. Any proposed
simultaneous compensation must replace/break the obstructing module or
otherwise leave that common fixed-module hypothesis. The author was notified
to keep this qualification explicit in the design discussion.

## Independent certificates

`verify_ideal_independent.py` imports no construction-lane code. It enumerates
actual exterior permutations and occupied slots, checks all1,099 labelled
colored contexts through four exterior vertices, and derives probe coefficient
rows independently through m=4:

```
scaled b rows:
[1]
[3,1]
[9,10,1]
[31,83,25,1]
[117,684,454,56,1]
```

The corresponding unscaled top coefficients are1,1,2,6,24, as predicted.
Direct subset Hamilton-path DP checks the physical probes for both five-vertex
examples, all three-vertex exterior contexts around them, and the triangular
recovery. Automorphism counts use all120 vertex permutations, independently.

- mask8: H=9, Aut=1, F=(9,66,186,240,120), ideal3.
  Probe values:9,93,1113,16497,298605.
- mask76: H=15, Aut=5, F=(15,90,210,240,120), ideal15.
  Probe values:15,135,1455,19875,337515.

The independent order-eight strengthening witness has H=513, Aut=1,
raw cover gcd 9 and scaled context gcd 27. The automorphism count checks all
40,320 vertex permutations. No smallest-witness claim is made.

All5,716 explicit checks pass with byte-identical normal and optimized output.
The all-context theorem rests on the free-action and unitriangular proofs;
the finite certificate is an independent control, not an extrapolation.
