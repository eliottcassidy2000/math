---
id: THM-2703
title: "C3 boundary-tree arm determinant and standard-plane gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Let a finite
  weighted tree carry an order-three automorphism and
  let M be its nonsingular invariant symmetric integral tree matrix.  The
  fixed vertices form a subtree and every complementary component occurs in
  a free triple.  If C runs through representatives of those triples, then
  the standard F2[C3]-part of (Z^V/MZ^V)[2] is exactly the direct sum of
  W tensor ker(M_C mod 2).  Thus its standard-plane multiplicity is the sum
  of the arm nullities, and a standard plane occurs iff some moving-arm
  determinant is even.  For a plumbing chain this determinant is the
  negative-continued-fraction numerator.  For the stated full-rank independent
  rational-surface boundary, localization also kills nonconstant units, so
  all-odd moving arms directly exclude the quartic V4 Kummer carrier.  D4 is
  sharp; an odd-arm star fails.  Non-tree and non-full-rank boundaries,
  positive realization, reflection completion, and all arbitrary
  A4/S4/JC2/DC2 conclusions remain open.
source: a4-resolvent-next-gate-scout-2026-07-28
depends_on:
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
related:
  - THM-2056-kelvin-polar-farey-defect-certificate
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar
  - THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate
  - THM-2686-coprime-galois-invariant-kummer-vanishing-and-a4-standard-plane-gate
  - THM-2700-danielewski-s3-resolvent-standard-plane-exclusion
script: 04-computation/c3_boundary_tree_arm_determinant_thm2703.py
output: 05-knowledge/results/c3_boundary_tree_arm_determinant_thm2703.out
script_sha256: 1db2ae94330e0e14760bf1c12910d2a9587378d499dcdb14516dacc039535eaa
output_sha256: 8e6207f3c43ba26beaaab1a9791a07f26d7ea067bb9afa2e01834da366f17fa7
hash_basis: LF-normalized bytes
---

# THM-2703 -- C3 boundary-tree arm determinant and standard-plane gate

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The quartic gate in THM-2655 is a two-dimensional `F2` representation on
which `C3` cycles the three nonzero vectors.  This theorem locates that plane
exactly in a `C3`-equivariant boundary tree.  The answer is local but not
vertexwise: it is the parity of the determinant of an entire moving arm.

The mechanism isolates why the primes two and three recur together:

```text
two:   the required V4-character packet is 2-torsion;
three: moving boundary components occur in C3-orbits;
F2[C3] = F2 direct_sum F4 separates fixed data from the standard plane.
                                                                    (1)
```

This is an algebraic discriminant-module theorem.  It does not turn a coarse
nerve, a virtual homology class, or a Farey parity into a physical LRC
carrier.

## 1. The weighted-tree decomposition

Let `Gamma` be a finite tree with vertex set `V`, and let `sigma` be an
automorphism of exact order three.  Let

```text
M: Z^V -> Z^V                                                (2)
```

be a symmetric integral matrix such that

```text
M_(uv)=0 when u!=v are not adjacent in Gamma,
M_(sigma u,sigma v)=M_(uv),
det(M)!=0.                                                   (3)
```

No definiteness assumption is needed.  Write `Gamma^sigma` for the fixed
vertex subgraph.  It is a nonempty connected subtree.  Indeed, a finite tree
has a fixed vertex or a fixed edge under every finite automorphism group.  An
odd-order automorphism cannot exchange the ends of an edge, so `sigma` fixes
a vertex.  The unique path between two fixed vertices is pointwise fixed.

Every component of `Gamma\Gamma^sigma` meets the fixed subtree in exactly one
edge and occurs in a free orbit of three components.  A stabilized component
would have a uniquely distinguished attachment vertex, hence a fixed vertex
outside `Gamma^sigma`, a contradiction.  Choose one representative `C` from
each component orbit and let `M_C` be the principal submatrix on its vertices.

Let

```text
W=ker(1+sigma+sigma^2 : F2[C3] -> F2[C3]).                  (4)
```

Then `W` is the unique nontrivial irreducible `F2[C3]`-module.  It has
dimension two, and `C3` cycles its three nonzero vectors.

### The arm-block theorem

There is an exact `F2[C3]`-module identity

```text
ker(M mod 2)_std
   = direct_sum_[C] W tensor_F2 ker(M_C mod 2),              (5)
```

where the sum is over the chosen moving-component representatives and
`std` means the `W`-isotypic part.  Consequently

```text
multiplicity_W ker(M mod 2)
   = sum_[C] nullity_F2(M_C mod 2).                          (6)
```

In particular,

```text
ker(M mod 2) contains W
 iff det(M_C) is even for at least one moving component C.   (7)
```

### Proof

Because `2` does not divide `3`,

```text
F2[C3]=F2 direct_sum W                                      (8)
```

is semisimple.  The node module decomposes as the trivial module on the
fixed subtree together with, for each component orbit,

```text
F2[C3] tensor F2^C
 = F2^C direct_sum (W tensor F2^C).                          (9)
```

The matrix `M mod 2` commutes with `sigma`.  There are no edges between
different complementary components.  The three attachment coefficients in
one component orbit are equal by invariance.  The map from the three root
coordinates to the fixed attachment vertex is therefore their sum, and the
map in the other direction is the diagonal vector.  Both maps factor through
the trivial summand of `(9)` and vanish on the standard summand.  Thus on
`W tensor F2^C`, the matrix is exactly

```text
1_W tensor (M_C mod 2).                                    (10)
```

Taking kernels over all component orbits proves `(5)` and `(6)`.  A square
integer matrix is singular modulo two exactly when its determinant is even,
which proves `(7)`.  Notice that `det(M)!=0` also rules out a rational kernel
on a moving arm, because `(10)` holds over characteristic zero as well.

## 2. The discriminant group, not merely the mod-two matrix

Put

```text
A_M=Z^V/MZ^V.                                               (11)
```

It is finite by `(3)`.  Tensoring the free resolution

```text
0 -> Z^V --M--> Z^V -> A_M -> 0                            (12)
```

with `F2` gives the equivariant identification

```text
A_M[2]=Tor_1^Z(A_M,F2)=ker(M mod 2).                        (13)
```

Combining `(5)` and `(13)` gives the promised exact formula

```text
(A_M[2])_std
 = direct_sum_[C] W tensor ker(M_C mod 2).                  (14)
```

Thus an even determinant is not a heuristic parity shadow.  It is exactly
the condition that the corresponding three-arm orbit contributes at least
one standard plane to the two-torsion discriminant module.

The formula also records multiplicity.  If an arm has mod-two nullity `e`,
its orbit contributes `e` independent standard planes, not just one Boolean
pass/fail bit.

## 3. Chains, continuants, and the precise Farey bridge

Suppose a moving component is a plumbing chain with self-intersections

```text
-a_1,...,-a_r,       a_i>=2,                               (15)
```

and adjacent intersections one.  Define

```text
D_0=1,
D_1=a_1,
D_j=a_j D_(j-1)-D_(j-2).                                  (16)
```

Then

```text
abs(det M_C)=D_r,                                          (17)
```

and `D_r` is the numerator of the negative continued fraction
`[a_r,...,a_1]^-` (or of the reversed convention after reversing the chain).
Therefore the chain contributes a standard plane exactly when

```text
D_r=0 mod 2.                                               (18)
```

For the all-minus-two chain,

```text
D_r=r+1,                                                   (19)
```

so precisely the odd-length arms contribute.

The exact modular/Farey connection is the transfer identity

```text
[ D_j      * ]
[ D_(j-1)  * ]
 = product_(i=j down to 1) [ a_i  -1 ]
                              [  1    0 ],                  (20)
```

Each factor is `T^(a_i)S` in `SL_2(Z)`, where projectively `S` is the order-two
generator and `TS` is the order-three generator of

```text
PSL_2(Z)=C2 * C3.                                          (21)
```

So the binary continued-fraction/Farey tree and the ternary symmetry really
do meet on the same transfer word.  What transfers to the present theorem is
only the continuant numerator and its parity.  THM-2056's LRC predicate also
depends on the signed Gram matrix, hull owner, and cone.  THM-2596 proves that
those data are not preserved by an unlabelled Farey move.  Hence `(20)` is a
faithful arithmetic bridge, not a new LRC safety certificate.

## 4. Rational-surface boundary gate

Let `U` be a smooth complex surface with an order-three action and suppose it
has a smooth projective rational `C3`-equivariant completion

```text
U subset Xbar,       D=Xbar\U.                              (22)
```

Assume:

1. `D` is a simple-normal-crossings divisor whose dual graph is a tree;
2. its component classes are independent and span a finite-index,
   `C3`-stable lattice `L` in `Pic(Xbar)`;
3. `M` is their intersection matrix.

The rational surface lattice `P=Pic(Xbar)` is free and unimodular.  The
divisor localization sequence and the full-rank hypothesis give

```text
Pic(U)=P/L,                                                 (23)
L subset P=P^* subset L^*,
Pic(U) injects C3-equivariantly into L^*/L=A_M.             (24)
```

The unit row is not an additional escape under these same hypotheses.  The
localization sequence begins

```text
C^* -> O(U)^* -> Z^D -> Pic(Xbar).                         (25)
```

Independence of the boundary component classes makes the last displayed map
injective.  Hence `O(U)^*/C^*=0`; and over `C`, every constant is a square.
Therefore

```text
O(U)^*/O(U)^(*)2=0.                                       (26)
```

It follows from `(14)` that

```text
Pic(U)[2] contains W
 ==> some moving boundary component C has det(M_C) even.    (27)
```

Equivalently, if every moving arm determinant is odd, the class-group branch
cannot supply a standard plane.

Now place this in the quartic setting of THM-2655.  Let `U` be the regular
locus of the full `C3`- or `S3`-Galois cubic-resolvent normalization of a
hypothetical degree-four Keller map, and suppose its `C3` restriction admits
a completion satisfying `(22)`--`(24)`.  The Kummer row is

```text
0 -> O(U)^*/O(U)^(*)2 -> H^1_et(U,mu_2)
  -> Pic(U)[2] -> 0.                                       (28)
```

Restriction to `C3` is semisimple, and the required `V4` character module is
`W`.  Together with `(26)`, this gives the direct exclusion

```text
all moving arms odd
 ==> this completion cannot be a quartic Keller resolvent.  (29)
```

An even arm remains only a necessary gate: it is not an existence theorem for
a Keller map.

## 5. Sharp controls

### 5.1 `D4` triality passes exactly

Take the negative `D4` Cartan tree: one fixed central vertex of weight `-2`
and three weight-`-2` leaves cycled by `C3`.  Each moving arm has determinant
two.  The full matrix has

```text
det(M)=4,
Smith(M)=(1,1,2,2),
A_M[2]=W.                                                   (30)
```

This is the smallest sharp model.  Its three nonzero discriminant classes
are cycled by triality, exactly as the three nonzero `V4` characters are
cycled in `S4/V4`.  The `D4` surface control in THM-2700 shows that this is
not merely an abstract module coincidence.

### 5.2 Odd moving arms fail despite the same visible symmetry

Keep the same three-arm star and central weight `-2`, but give every leaf
weight `-3`.  Then

```text
abs(det M)=27,
A_M[2]=0.                                                   (31)
```

The graph still has literal triality, but there is no two-torsion on which
it can act.  Thus a `C3`-symmetric boundary graph alone is insufficient; the
even arm determinant is load-bearing.

### 5.3 Failure boundaries

The proof does not cover:

- a boundary dual graph with cycles;
- boundary components which do not form a full-rank independent lattice
  (where nonconstant global units may reappear);
- the reflection needed to upgrade the `C3` plane to a specified full `S3`
  action;
- non-tree incidence complexes, virtual nerves, or coarse carry-forgetting
  homology;
- realization of an even arm by an actual quartic resolvent;
- arbitrary `A4`, `S4`, planar Jacobian, or Dixmier candidates.

The first two failures destroy the direct-sum arm decomposition, the embedding
`(24)`, or the unit vanishing `(26)`.  The reflection item is why an even arm
does not classify the whole `S3` action.

## 6. Exact companion

The exact script

```text
python 04-computation/c3_boundary_tree_arm_determinant_thm2703.py
python -O 04-computation/c3_boundary_tree_arm_determinant_thm2703.py
```

checks:

- `6,752` weighted trees with a fixed subtree and one or two moving arm
  orbits, directly comparing the full standard kernel with `(6)`;
- `19,530` chains of lengths one through six and weights two through six,
  comparing determinants, transfer numerators, and mod-two singularity;
- the all-minus-two recurrence through length twelve;
- the `D4` Smith form and standard plane;
- the determinant-`27` odd-arm hostile.

Normal and optimized transcripts byte-match the stored output.  The finite
audit checks the theorem's consequence object; the general proof is
Sections 1--4.

## 7. Exact connection contract

```text
source object:
  a C3-equivariant weighted boundary tree;

target object:
  the standard part of its two-torsion discriminant module;

map:
  remove the fixed subtree, choose one component per free triple, and reduce
  each principal arm matrix modulo two;

preserved predicate:
  multiplicity of the nontrivial irreducible F2[C3] plane;

destroyed information:
  the reflection action, non-full-rank units, non-tree gluing, realizability,
  LRC owner/phase/height, and physical endpoint transport;

needed sidecars:
  an S3 reflection completion and a geometric completion theorem for the
  quartic lane; Gram-owner and physical-current data for any LRC transfer;

cheapest decisive next test:
  compute an equivariant SNC boundary lattice for each surviving quartic
  resolvent normal form and inspect the continuants of its moving arms.
```

That test is strictly sharper than comparing discriminant polynomials: it
asks whether the boundary can physically carry the `V4` character plane.

An independent hostile audit rederived the fixed-subtree/free-triple
decomposition, the exact standard-block identity, the discriminant-group and
Picard-lattice steps, continuant parity, both sharp controls, and every scope
boundary.  Normal and optimized runs byte-match the stored transcript, the LF
hashes match, and the companion contains no optimized-away assertions.

A subsequent independent localization audit sharpened the rational-surface
corollary: the already-assumed independence of the boundary classes makes
`Z^D -> Pic(Xbar)` injective, so the unit Kummer term vanishes as in
`(25)`--`(26)`.  This changes the full-rank conclusion from a class-group-only
gate to the direct all-odd exclusion `(29)`; it does not change the tree
module formula or the exact companion.

QED.
