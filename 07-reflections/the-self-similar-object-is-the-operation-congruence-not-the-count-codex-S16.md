# The self-similar object is the operation congruence, not the count

**Session:** codex-2026-07-15-S16
**Scope:** H-drift, Farey/toothpick ladders, scale-three Hamming-six sheets,
Fano/`chi_7`, `j=4` floods, Kakeya needles, and the black self-line guardrail

## 1. The common diagnosis

The recurring mistake across these threads is to recognize the right finite
set but the wrong dynamics.  Equal counts, isomorphic incidence graphs, equal
charge, or a common tournament fingerprint do not justify identifying two
states.  A compression is proof-facing only when it intertwines the legal
operations and retains the terminal predicate.

Let `X` be a state space with partial labelled transitions `T_a`, let `P` be
the terminal predicate, and let `q:X->Y` be a proposed quotient.  The exact
test is:

```text
q(x)=q(x')
  => [a is legal at x iff a is legal at x']
     and q(T_a x)=q(T_a x') whenever a is legal,            (1)

P(x)=Pbar(q(x)).                                            (2)
```

Under (1)--(2), every finite recursive decision built from the transitions
factors through `q`, by induction on the remaining depth.  Conversely, a
single kernel pair violating legality, successor class, or terminal verdict
is a certificate that no deterministic recursion on `Y` alone can be both
sound and complete.  This elementary operation-congruence criterion is the
right meaning of self-similarity in this workspace.

It separates four phenomena that had repeatedly been conflated.

| phenomenon | exact example | what really recurs | what does not |
|---|---|---|---|
| Koopman observable conjugacy | THM-848 H-drift | future averaged `H` observables in the Mobius radial coordinate | the state transition law, diffusion, or an LRC witness |
| conditional fibre recursion | THM-862 toothpick code | affine matching constraints under sign pinning | literal component erosion and future rays |
| sourced self-affinity | THM-850/841 dyadic `chi_7` stalk | a three-step scaled stalk plus a birth term | the full Farey multiplier ladder |
| census coincidence | early toothpick totals and `2 selfK=SC` at `n=5,6,7` | nothing yet | an all-size bijection or quotient |

The missing object in the LRC branches is therefore not one more scalar.  It
is a fibred transition object: a symbolic operation base together with the
literal metric stalk needed to make (1)--(2) true.

## 2. H-drift is an observable positive control, not a state quotient

THM-848 gives

```text
A_T(x)=sum_r H_(2r)(1-x)^(2r)(1+x)^(m-2r),   m=n-1,
z=(1-x)/(1+x),
B_T(z)=sum_r H_(2r)z^(2r).                                (3)
```

The flip-sum generator on these observables becomes exactly

```text
-2z d/dz,                                                 (4)
```

so the Poissonized continuous flow is radial dilation
`B(z)->B(exp(-2t)z)`.  The discrete uniform-flip operator is instead
`P B=B-(2/M)zB'`, `M=binom(n,2)`, with modes `(1-4r/M)^t`; it is not substitution by one
scalar dilation once degrees two and four both occur.

More importantly, even the full polynomial is not a transition quotient.
At `n=5`, full-arc masks `8` and `10` have the same
`A=(9,30,42,30,9)`, `B=(3/2,0)`, `H=9`, `c3=3`, score sequence, and mean
drift `-3/5`, but their one-step `H` diffusions are `14` and `78/5` and their
target-`B` histograms differ.  The same labelled flip `{1,2}` already sends
their stalks to `(-3/2,-1)` and `(-9/2,0)`.  Thus radiality is a Koopman
conjugacy for averaged observables, not an operation-congruence of tournament states.  It
still explains why `(H,K)` eventually fails: those are only the first Krylov
coordinates, and at `n=8` a higher even layer becomes visible.

This sharpens a general rule.  A functional form may close a generator on
observables without closing the underlying transition kernel.  Matching a few
moments, future means, or small-`n` values is not strong lumpability.

## 3. Toothpicks occur in three different senses

### 3.1 The full Farey ladder is a renewal object

THM-841's one-kink and totient laws show that the full endpoint ladder is not
the classical toothpick sequence.  Its early dyadic agreements are count
coincidences; the total already differs at the next audited power
(`159 != 171` at `n=16`).  The exact defect is informative: every rescaling
creates new odd-reciprocal births.  Thus the honest equation is a renewal or
sourced recursion, not homogeneous self-similarity.

The `chi_7`-labelled dyadic subladder has the exact restricted law

```text
W_a^H(s)=zeta 1_(s<a)+zeta^2 1_(s<a/2)+zeta^4 1_(s<a/4)
         +W_a^(H-3)(8s).                                  (5)
```

The last term is a scaled copy; the first three are the source.  The full
Farey ladder adds the odd births omitted from (5).  This is toothpick
self-similarity with immigration.

At the local rewrite level there is a cleaner exact carrier.  If `(i,j)` are
the ordered denominators of a gap of `F_k`, then the next order replaces it by
`(i,i+j),(i+j,j)` exactly when `i+j=k+1`, and otherwise leaves it fixed.  The
coarse chain-length code already fails at `k=6`: pairs `(4,3)` and `(5,3)`
both map to `(1,2)`, but only the first splits at order seven.  Thus the
ordered Stern--Brocot address is operation-congruent for the interval core;
the toothpick count is not.

### 3.2 The scale-three sheet bank has affine toothpick codes

In THM-862, the common-sheet predicate really does factor through disjoint
signed pair equations.  For matching `M`,

```text
1_C(e)=2^(-|M|) product_(ij in M)(1+s_ij e_i e_j).         (6)
```

Pinning signs recursively halves matching/free components exactly.  This is
a genuine conditional self-similarity of the sheet fibre, and its Walsh
support is precisely the unions of matching edges.

It is not a metric self-similarity.  The same sheet orbit contains packets
with exact maxima `1/4`, `7/26`, and `14/65`.  After the first insertion,
`22,262` literal child geometries cache `146,912` logical lanes, but adjoining
the five actual future step-39 rays separates every lane.  The correct object
is a cache bundle, not a quotient tree.

### 3.3 The self-line law is a census guardrail

The proposed all-size black law is false.  The black quasi-fixed endpoint
counts at `n=5,6,7,8` are

```text
2 selfK(n)=8,12,88,404,
SC(n)      =8,12,88,176.                                  (7)
```

The `n=8` excess is entirely black; it is not caused by accidentally counting
blue self-lines.  What survives all size is the Klein-four orbit structure,
evenness/divisibility, and the path-relative holonomy law in THM-849/852/854.
The subsequent live-main freeness proof sharpens that replacement: on the
non-grid-symmetric quasi-fixed locus, `g` is excluded by definition, `kappa`
flips every tile bit, and `g kappa` cannot fix an anti-diagonal tile.  Thus the
Klein four acts freely for every size, so

```text
selfK(n)=2 orbitK(n),       orbitK(5..8)=2,3,22,101.
```

These are the merged-carrier orbit counts; they are not `SC(n)/2` at `n=8`.
This is the same lesson as the Farey ladder: an initial count identity can be
the shadow of a real group action without being the correct orbit count.  The
all-size theorem belongs to freeness of the action, not equality of the two
twisted diagonals.

## 4. The Fano carrier is an atlas, not the metric object

The `chi_7` operation charge is exact because it is a monoid homomorphism:

```text
q(A),q(B),q(C)=1,2,4 in F_7.                              (8)
```

It recovers the Boolean Mobius signs and produces the Paley tournament on
the seven nonempty face strata.  But one charge pencil has rank only seven on
the 49-point operation plane.  All eight projective directions recover an
arbitrary residue-plane function, and at depth eleven the residue pushforward
still loses thirty A/B/C carry contrasts.  The reflection-invariant six-pencil
repair works only after those carry channels are retained.

The flood bank has a parallel two-level structure.  Its 21 root edges are
the 21 Fano flags.  Point-star and Fano-line marginals have rank thirteen;
the eight missing dimensions are exactly `H_1` of the Heawood graph and are
recovered by oriented hexagon currents.  This reconstructs every scalar edge
field, including `m` and `V1`.

Yet a completed family has several root-edge presentations.  Re-rooting over
the 28 non-Fano triples coequalizes all 21 edge charts, so every edge-only
scalar descending to completed families is constant.  Heawood currents and
re-root equalities therefore move in opposite categorical directions:

```text
Heawood cycles: recover information inside the root atlas;
re-root triangles: identify charts of one completed family.                (9)
```

The only proved positive transport between flood roots is literal
containment.  The three independently closed bodies `(5,6)`, `(5,7)`, and
`(6,7)` form a triangle of anchors.  They shadow 18 of 21 five-small bases;
the other three require 260 exact one-speed sweeps.  Every six-small base is
anchored.  At four small labels the triangle shadows 22 of 35 bases.  On each
of the thirteen residual bases, THM-732 reduces the two-external-speed tail to
a finite bank; all 29,183 exact pair measures are positive.  Thus every
completion with at least four small labels is lonely, and the unresolved tail
has at most three small labels.  This is useful because containment and the
tail inequality commute with the actual extension operation.  It still does
not identify any of the eighteen remaining whole bodies.

The natural peel-role repair also fails the operation square.  Classify an
insertion as split-dominant, balanced, or kill-dominant (`A/B/C`) by component
Euler drift.  In the root `{1,2,8,...,14}`, histories `(3)` and `(4)` have the
same role `C` and the same unordered birth/death profile, but the common legal
next speed `26` produces `CC` and `CB`; their `chi_7` charges split to `1` and
`6`.  What the quotient forgot is not another character: the killed
components occur in cyclic run patterns `(8,20)` and `(12,16)`.  Position
along the component circle is the first concrete repair.

## 5. Kakeya needles need offsets, owners, and chronology

A direction-only Kakeya picture is too coarse in both the finite-operation
and metric-comb settings.  A line direction forgets its affine offset; a
charge fibre forgets position along the fibre; a residue address forgets its
carry; an endpoint needle forgets its opposite endpoint and owner; an
unordered collection of needles forgets wall chronology.

The proof-facing needle should instead be written as a decorated obligation

```text
(direction, offset, owner, open/closed side, sheet, birth scale,
 component incidence, legal future ray, chronology).       (10)
```

Different problems can forget different coordinates, but only after testing
(1)--(2).  The three-endpoint-needle obstruction in THM-850 is then correctly
scoped: three one-sided needles visit at most four Boolean masks and at most
two points of the negative Fano line.  A full runner has two cyclic ends and
can switch owners, so the obstruction is not a three-runner theorem.

This reframing also explains the seven-comb wall.  Pairwise projective ratios
and a maximum-spanning-tree rank word recover the second-order Hunter credit,
but the irreducible residual is a prefix--comb--comb third moment.  It is the
incidence and chronology part of (10), not another pairwise direction.

## 6. The underlying LRC state

The views above converge on the following minimal current state for the
AP-centred ramified Hamming-six recursion:

```text
symbolic base:
  common scale and complement-lcm fibre,
  missing labels and effective orders,
  signed provider/matching code,
  remaining labelled arithmetic rays;

metric stalk:
  literal strict-safe components with endpoint ownership,
  last inserted speed and numerical clock,
  sheet/carry data,
  shortcut/equality ancestry;

transition:
  intersect by one legal ray, update endpoints, then apply only
  operation-congruent certificates.                         (11)
```

The symbolic base is where Fano charge, toothpick codes, tournaments, and
Walsh projectors live.  The stalk is where loneliness is decided.  Geometry
caches may be shared across fibres, but the labelled future language must
remain attached.

This is better described as an **operation-intertwining orbit bundle** than
as a graph quotient.  It is recursive without pretending to be homogeneous;
it permits sourced births, ramification, and cache reuse; and every proposed
compression has a concrete kernel-pair test.

## 7. Ranked next tests

1. Shard complete retained-root blocks with the now-validated scale-three
   geometry-batched depth-three evaluator.  Its bounded prototype matches
   1,685,358 raw lanes at 1.599x CPU speedup; propagate every labelled future
   language and do not turn the cache into a quotient.
2. For `j=4`, move to the 22 unshadowed three-small completed-family bases now
   that all thirteen four-small two-external-speed banks are closed.  Seek a
   genuine three-speed tail inequality, or use `(3,4)` as the coverage-first
   whole-body job; do not infer it from the Fano orbit.
3. Replace the refuted Euler-drift role word by a cyclic component-incidence
   history retaining killed-block phase, then retest the operation-congruence
   square before using `chi_7` or Radon coordinates on flood histories.
4. Treat the full Farey ladder as a renewal equation: isolate its odd-birth
   source and ask whether the sourced operator, rather than its total count,
   has a contractive norm.
5. For every tournament viewpoint, name the vertices, pair observable,
   switch, tie path, preserved predicate, and a kernel pair.  A transitive
   fingerprint is a clock until proved otherwise.

The mathematical frontier has sharpened from “find the right invariant” to
“find the smallest state on which the legal operations are a congruence.”
That is a narrower and more falsifiable target.
