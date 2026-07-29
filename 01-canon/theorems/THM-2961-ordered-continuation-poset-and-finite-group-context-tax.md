---
id: THM-2961
title: "Ordered continuation poset and finite-group context tax"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT. For a poset-valued observable on any right monoid, pointwise
  comparison under every continuation is the largest right-compatible
  preorder inside the original scalar order. Its symmetric kernel is the
  continuation equivalence, and for a totally ordered observable its
  intrinsic crossing graph is a co-comparability graph and is perfect on
  every finite packet. The T(2,7) mirror pair is an exact contextual crossing
  with a minimal two-query certificate. For the discrete response on every
  nontrivial finite group G, the crossing graph is K_|G| and every
  order-complete physical context dictionary has exactly |G| entries.
  Consequently C2, C3, and V4 have exact context taxes 2, 3, and 4. The V4
  response has affine S4 symmetry, normal translation V4, quotient S3, no
  invariant tournament orientation, and no partial-cube conclusion.
  Pointwise order on the full two-ended metric kernel is equality, while a
  non-group metric hostile proves that contextual crossing itself need not
  survive translation.
source: codex-knot-contextual-order-v4-2026-07-29
depends_on:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2242-tournament-complement-transport-and-knot-kernel-green-rigidity
  - THM-2281-common-optimal-context-for-finite-catalytic-families
related:
  - THM-2348-prime-type-rectangularity-and-target-token-conditioning
  - THM-2646-braid-three-modular-central-pullback-and-full-twist-knot-fibre
script: 04-computation/knot_ordered_continuation_poset_thm2961.py
output: 05-knowledge/results/knot_ordered_continuation_poset_thm2961.out
script_sha256: 31859066c2d3008f971fe9c5f65d62d4637cf03b44ee97e16144754d62fe5bf5
output_sha256: 7c25a38fd500650354e13916c863b784d9d1f197d7cc5183dc4ef0986318a29c
hash_basis: working-tree bytes (LF)
---

# THM-2961 -- ordered continuation is a poset, not a tournament

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

THM-2176 proves that equality of complete continuation profiles is the
largest operation congruence inside equality of an observable. This theorem
adds the ordered layer. It identifies the intrinsic binary relation, proves
its finite graph class, and isolates exactly where the `C_2`, `C_3`, and
affine `V_4` pictures are faithful.

## 1. Universal ordered continuation

Let `M` be a monoid, not necessarily commutative. Write its operation
additively and use **right** continuations. Let `(Y,<=)` be a poset and let

```text
f:M->Y.
```

Define

```text
x <=_f^ctx y
  iff f(x+z)<=f(y+z) for every z in M.                 (1)
```

> **Ordered continuation theorem.** The relation `<=_f^ctx` is the largest
> right-translation-compatible preorder contained in the pullback of `<=`
> along `f`.

### Proof

Reflexivity and transitivity hold pointwise. If (1) holds, then for every
`a,z in M`,

```text
f((x+a)+z)
 =f(x+(a+z))
 <=f(y+(a+z))
 =f((y+a)+z).                                         (2)
```

Thus right translation preserves the preorder. Taking `z=0` in (1) gives

```text
f(x)<=f(y).                                           (3)
```

Conversely, let `R` be a right-translation-compatible preorder such that

```text
x R y  implies  f(x)<=f(y).
```

If `x R y`, then `x+z R y+z` for every `z`, and hence
`f(x+z)<=f(y+z)`. Therefore `R` is contained in (1), proving maximality.

Its symmetric kernel is

```text
x <=_f^ctx y and y <=_f^ctx x
 iff f(x+z)=f(y+z) for every z
 iff Pi_f(x)=Pi_f(y),                                 (4)
```

where `Pi_f(x)(z)=f(x+z)`. This is exactly THM-2176's continuation
equivalence in the commutative setting and its right-congruence analogue in
general.

Let

```text
P_f=M/equiv_f.
```

The preorder descends to a partial order on `P_f`, and

```text
[x] -> Pi_f(x) in Y^M                                (5)
```

is an order embedding for the pointwise product order. QED.

## 2. Contextual crossing and the finite graph class

Assume now that `Y` is totally ordered. Define

```text
x cross_f y
```

when there are contexts `z_+,z_-` such that

```text
f(x+z_+)>f(y+z_+),
f(x+z_-)<f(y+z_-).                                   (6)
```

On `P_f`, this is exactly incomparability. It is symmetric and has no
intrinsic orientation.

For any finite `A subset P_f`, form the graph whose vertices are `A` and
whose edges are contextual crossings. Its cliques are precisely antichains
of the induced poset, while a proper color class is precisely a chain.
Dilworth's theorem therefore gives

```text
chi(cross_f|A)
 =minimum number of chains covering A
 =maximum size of an antichain in A
 =omega(cross_f|A).                                   (7)
```

The same argument applies to every induced subgraph. Hence every finite
crossing graph is a co-comparability graph and is perfect.

The preorder, not the crossing graph, is the operation-ready object. If

```text
S_f(x,y)
 ={sgn(f(x+z)-f(y+z)):z in M},
```

then

```text
S_f(x+a,y+a)
 ={sgn(f(x+t)-f(y+t)):t in a+M}
 subset S_f(x,y).                                    (8)
```

A noninvertible right translation can therefore resolve a crossing into
dominance or a tie. If `M` is a group, `a+M=M`, so equality holds in (8)
and every pair type is translation invariant.

## 3. Physical context dictionaries

For a finite labelled packet `A subset M`, call `D subset M`
**order-complete for `A`** when

```text
x <=_f^ctx y
 iff f(x+z)<=f(y+z) for every z in D,
                                      for all x,y in A. (9)
```

The least possible `|D|`, when finite, is the physical context-query number
of the packet. It is not the abstract order dimension: its coordinate chains
must be supplied by actual continuation evaluations.

Failure of dominance has a one-context certificate. Crossing has a
two-context certificate. Dominance itself retains the universal quantifier
in (1); neither a favorable single context nor a finite common-translate
statement proves it.

## 4. The exact knot crossing

Specialize to the commutative monoid of oriented knots under connected sum
and put `f=u`, the unknotting number. Let

```text
K=T(2,7),               Kbar=mirror(K).
```

THM-2176 gives

```text
u(K)=u(Kbar)=3,
u(K#K)=u(Kbar#Kbar)=6,
u(K#Kbar)<=5.                                           (10)
```

The diagonal equalities follow from the signature lower bound and separate
unknotting. The off-diagonal entry is only the certified upper bound from
the Brittenham--Hermiller crossing-change sequence; its exact value is not
used.

Using contexts `K` and `Kbar` gives the opposite strict comparisons

```text
u(K#K)>u(Kbar#K),
u(K#Kbar)<u(Kbar#Kbar).                                (11)
```

Therefore

```text
[K] cross_u [Kbar].                                    (12)
```

This sharpens “their continuation profiles differ”: the scalar tie
`u(K)=u(Kbar)` splits into a contextual antichain. One scalar context cannot
order-reflect a two-element antichain, while the two displayed contexts do.
Thus this packet has physical context-query number exactly two.

Mirroring exchanges the two knots and the two witness contexts while
preserving unknotting number. The symmetric crossing edge is therefore
mirror-gauge invariant. An orientation of that edge is reversed by the
gauge, so no mirror-invariant tournament relation exists on the packet.

THM-2281 remains complementary. It realizes every finite catalytic metric
packet in one common ordinary translated slice. Equation (1) instead
quantifies over the entire future cone. The finite-slice theorem neither
proves nor refutes contextual dominance.

## 5. The full metric kernel has only equality order

THM-2242 uses the two-ended continuation kernels

```text
P_x(a,b)=d(x+a,b).
```

Their pointwise order is rigid:

```text
P_x<=P_y pointwise
 => d(x,y)=P_x(0,y)<=P_y(0,y)=0
 => x=y.                                                (13)
```

Thus no nontrivial dominance survives on the full faithful kernel. The
contextual poset arises only after projecting to the root observable
`f(x)=d(x,0)` and retaining all its translates. This is controlled
forgetting: the root projection creates dominance and incomparability, while
the continuation family repairs the scalar ties.

## 6. The finite-group context tax

Let `G` be any nontrivial group with the discrete metric and root response

```text
f(e)=0,
f(g)=1 for g!=e.                                       (14)
```

Write the group operation multiplicatively in this section. For distinct
`x,y`, the unique right context satisfying

```text
f(xz)>f(yz)
```

is

```text
z=y^(-1).                                              (15)
```

Indeed, at `y^(-1)` the right side is zero and the left side is one. At
every other context the right side is one, so a strict inequality in that
direction is impossible.

It follows that:

```text
the continuation quotient is an |G|-element antichain;
the crossing graph is K_|G|;
if G is finite, every order-complete dictionary has at least |G| entries;
the full dictionary G is order-complete.                (16)
```

For the lower bound, fix `y`. Refuting `x<=_f^ctx y` for any `x!=y` requires
the unique context `y^(-1)`, so every inverse must occur in the dictionary.
Thus the physical context-query number is exactly `|G|`.

For an infinite group, no finite dictionary is order-complete for the whole
group. In particular the discrete response on

```text
C_2*C_3
```

has infinite context tax. Every finite quotient identifies distinct free
words and hence cannot preserve these discrete continuation profiles.
This is a precise stopping boundary: the binary and ternary factors do not
become one finite grammar merely because a finite quotient contains elements
of orders two and three.

The first finite cases are sharp:

```text
G=C_2:  crossing graph K_2, context tax 2;
G=C_3:  crossing graph K_3, context tax 3;
G=V_4:  crossing graph K_4, context tax 4.             (17)
```

Thus `C_2` gives the last complete graph which is a partial cube, and `C_3`
gives the first non-bipartite complete graph.

## 7. The affine `V_4` quotient and its tournament boundary

Identify

```text
V_4=F_2^2.
```

For the discrete response, its four continuation profiles are the four rows
of `J_4-I_4`: each row has one zero at its inverse context and ones
elsewhere. Simultaneously permuting object and context labels preserves this
response, so the full response automorphism group is

```text
S_4=AGL(2,2).
```

The affine translation subgroup is a normal regular `V_4`, and

```text
S_4/V_4=GL(2,2) isomorphic to S_3.                     (18)
```

The quotient permutes the three nonzero affine directions and forgets the
origin. But Section 6 proves that all four affine contexts are necessary for
an order-complete physical dictionary. The abstract order dimension of a
four-element antichain is only two, so neither two abstract linear orders nor
the three quotient directions replace the missing affine origin.

The unweighted crossing graph is `K_4`. Every transposition reverses the
edge between the two points it swaps, so no orientation of `K_4` is invariant
under full `S_4`. Moreover `K_4` contains a triangle and is not a partial
cube. The intrinsic object is the symmetric incomparability relation with
its labelled contexts, not a tournament.

## 8. A sharp non-group hostile

Crossing itself is not a right congruence. Let `M=N` under addition and use
the integer metric

```text
d(m,n)=0,                m=n;
       =1,                |m-n|=2;
       =2,                otherwise.                 (19)
```

This is the restriction of a translation-invariant metric on `Z`. If both
legs of a triangle are nonzero, their sum is at least two, which dominates
every direct distance; a zero leg gives equality. Translation invariance and
the triangle inequality also give joint nonexpansivity:

```text
d(a+c,b+e)<=d(a,b)+d(c,e).
```

Its root length is

```text
f(0)=0,
f(2)=1,
f(n)=2 for n>0, n!=2.                                (20)
```

The elements `1,2` cross:

```text
f(1)>f(2),
f(2)<f(3).                                           (21)
```

After translating both by one, however,

```text
2 <_f^ctx 3,                                         (22)
```

because the inequality is strict at the zero context and all later
responses are equal to two. This realizes strict inclusion in (8) even in a
jointly nonexpansive integer metric monoid.

## 9. Scope and information ledger

The proved mechanism is:

```text
source:
  a right continuation profile Pi_f;

target:
  its pointwise preorder, quotient poset, and symmetric incomparability
  graph;

preserved:
  every future scalar comparison and every continuation tie;

gauge:
  continuation equivalence; for the knot crossing, simultaneous mirroring
  of object and context;

ties:
  equality of complete profiles, not equality at the empty context;

lost by the crossing graph:
  response magnitudes, witness contexts, dominance direction, and the monoid
  operation;

sidecar:
  actual context coordinates, rather than abstract linear extensions or a
  quotient direction set.                              (23)
```

This theorem does **not** compute a new unknotting number, prove positive
translation catalysis, classify knot continuation profiles, make the knot
crossing graph a partial cube, or produce a tournament. The free-product tax
is for the discrete response only. Perfectness is a finite induced-graph
statement. The non-group hostile shows why only the preorder and its tie
kernel, not crossing edges, may be transported as operation relations.

## 10. Exact companion

The exact companion checks:

- the `C_2`, `C_3`, and `V_4` discrete response profiles, complete crossing
  graphs, unique directed witnesses, and exact context taxes `2,3,4`;
- all `24` affine permutations of `F_2^2`, normality of the four
  translations, and the quotient size six;
- all `64` orientations of `K_4` against all `24` affine permutations,
  finding no invariant tournament;
- the full-kernel pointwise-order rigidity on the `V_4` control;
- `9,261` metric triangles, `1,728` common translations, and `20,736` joint
  nonexpansivity inequalities for the non-group hostile;
- the crossing-to-dominance transition and the knot mirror upper-bound
  certificate.

Reproduce with

```bash
python 04-computation/knot_ordered_continuation_poset_thm2961.py
python -O 04-computation/knot_ordered_continuation_poset_thm2961.py
```

Both transcripts must match

```text
05-knowledge/results/knot_ordered_continuation_poset_thm2961.out
```

byte-for-byte after LF normalization.
