---
id: THM-765
title: Safe-component tooth decks force hereditary primitivity
status: PROVED (elementary circle geometry; the twelve-speed corollary cites LRC for eleven speeds)
source: codex-2026-07-14-S3 (n=12 component/deck audit)
renumber_note: first pushed as THM-762 after THM-762 was already occupied; moved to THM-765 by the first-pusher protocol
depends_on:
  - THM-592   # exact safe-component language and endpoint geometry
  - THM-619   # earlier component-tooth band criterion
  - LRCUpTo13 # only for M(P) >= 1/12 in the leave-one-out corollary
related:
  - THM-759   # one-component width bound / ratio bound
  - HYP-6775  # tight twelve-speed rigidity
  - HYP-6800  # the super-lonely-core (sporadic) branch
  - HYP-6820  # uniformity audit
  - MISTAKE-144 # corrected self-cusp support for the bounded census
---

# THM-765 — Safe-component tooth decks force hereditary primitivity

## Statement

For a finite nonempty set of positive integer speeds `P`, put

```text
Phi_P(t) = min_{p in P} ||p t||,       M(P) = max_t Phi_P(t),
E_L(P)   = {t in R/Z : Phi_P(t) > L}.
```

Let `0 < L < 1/2`, let `w` be another positive integer, and suppose

```text
M(P union {w}) <= L < M(P).                              (1)
```

Then the following hold.

### A. Exact safe-component / tooth containment

Every connected component `J` of `E_L(P)` has its closure in one unique
closed `w`-tooth.  Precisely, the teeth on `R/Z` are the images, indexed by
`k in Z/wZ`, of the real intervals

```text
T_k = [(k-L)/w,(k+L)/w].
```

After choosing compatible real lifts, this is the usual condition
`|wt-k|<=L`.  Write `J=(a,b)`, `c=(a+b)/2`, and `h=(b-a)/2`.  Tooth
containment is equivalent to the exact midpoint band

```text
||w c|| + w h <= L,                                    (2)
```

or, equivalently,

```text
w (b-a) <= 2L,
||w c|| <= L - w(b-a)/2.                               (3)
```

Thus the length tax and the midpoint residue band are not separate
heuristics: together they are exactly equivalent to coverage of the whole
component by one tooth.

### B. Quantitative deck obstruction

Let

```text
d = gcd(P),       g = gcd(d,w),       D = d/g.
```

If `D >= 2`, then

```text
M(P union {w}) >= min(M(P), (D-1)/(2D)).                (4)
```

Consequently, when `0 < L < 1/4`, (1) forces

```text
d | w.                                                  (5)
```

In words: a single runner can push a super-`L` core down to `L` only if it
respects every translation sheet of the core's gcd deck.

### C. Hereditary primitivity at twelve speeds

Let `A={a_1,...,a_12}` be a primitive twelve-speed set with

```text
M(A) = 1/13.
```

Then every leave-one-out core is primitive:

```text
gcd(A \ {a_i}) = 1                 for every i.          (6)
```

Equivalently, for every prime `ell`, at least two speeds of `A` are not
divisible by `ell`.  In particular, the max-peeled super-lonely core in the
putative `n=12` sporadic branch is necessarily primitive.  No imprimitive
core, at any height or winding number, can inhabit that branch.

### D. Infinite high-winding complete-residue family eliminated

Let `d>=2`, choose `j in {1,...,12}`, and choose an integer `h` such that
`w=dj+13h>0` and `gcd(d,13h)=1`.  Then

```text
A(d,j,h) = {d r : 1<=r<=12, r!=j} union {dj+13h}
```

is primitive and satisfies the uniform strict floor

```text
M(A(d,j,h)) >= 1/12 > 1/13.                            (7)
```

If additionally `d=1 (mod 13)`, this is a complete nonzero residue system
modulo `13`: it is exactly a one-coordinate defect of the dilation
`d*{1,...,12}`.  Thus no such complete-residue lift is tight, at any winding
height.

## Proof

### Proof of A

Assumption (1) says that every `t` with `Phi_P(t)>L` must satisfy
`||w t||<=L`; otherwise the same `t` would give clearance strictly larger
than `L` for `P union {w}`.  Hence

```text
E_L(P) subset {t : ||w t|| <= L}.                        (8)
```

For `L<1/2`, the set on the right is a disjoint union of closed teeth of
length `2L/w`, separated by positive gaps.  A connected component `J` of
the open set on the left therefore lies in one tooth.  Since the tooth is
closed, `closure(J)` lies in that same tooth, and separation makes the tooth
unique.

For real lifts, containment in the tooth centred at `k/w` is

```text
(k-L)/w <= c-h <= c+h <= (k+L)/w,
```

which is equivalent to `|wc-k|+wh<=L`.  Because `L<1/2`, such an integer
`k` is the unique nearest integer to `wc`; this gives (2), and (3) is just
the same inequality split into its nonnegative length and midpoint parts.
The converse follows by choosing that nearest integer, so the criterion is
an equivalence.  This proves A.

### Proof of B

Write every `p in P` as `p=d b_p`, and choose a maximizer `t_0` of the core,
so `Phi_P(t_0)=M(P)`.  For `j=0,...,d-1`, set

```text
t_j = t_0 + j/d  (mod 1).
```

Every core phase is unchanged on this orbit:

```text
p t_j = p t_0 + b_p j  (mod 1),
```

and hence `Phi_P(t_j)=M(P)` for all `j`.

The killer phases on the same orbit are

```text
w t_j = w t_0 + wj/d  (mod 1).
```

After repetitions are removed, these are a translate of the complete
`D`-grid `{0,1/D,...,(D-1)/D}`.  Any circular arc containing a complete
`D`-grid has length at least `(D-1)/D`: one may omit at most one of its
`D` equal gaps, and every gap has length `1/D`.  Therefore, for every
translate of the grid,

```text
max_j ||w t_j|| >= (D-1)/(2D).                           (9)
```

(The constant is sharp: rotate the grid so the shortest containing arc is
centred at zero.)  Choose a `j` realizing (9).  At `t_j`, the core clearance
is `M(P)` and the killer clearance is at least `(D-1)/(2D)`, proving (4).

If `D>=2`, then `(D-1)/(2D)>=1/4`.  Thus (4), `M(P)>L`, and `L<1/4`
would give `M(P union {w})>L`, contrary to (1).  Hence `D=1`, which is
exactly `d|w`.  This proves B.

### Proof of C

Fix `a_i in A` and let `P=A\{a_i}`.  This core has eleven speeds, so the
settled eleven-speed Lonely Runner Conjecture gives

```text
M(P) >= 1/12 > 1/13.
```

Apply B with `L=1/13` and `w=a_i`.  It gives `gcd(P)|a_i`.  But every member
of `P` is also divisible by `gcd(P)`, so `gcd(P)` divides all of `A`.
Primitivity of `A` forces `gcd(P)=1`.  This works for every `i`, proving
(6).  The prime-divisibility formulation is equivalent: a prime divides
some leave-one-out gcd exactly when it divides at least eleven of the twelve
speeds.  This proves C.  ∎

### Proof of D

Take

```text
P = d*({1,...,12}\{j}),       w=dj+13h.
```

Scaling does not change `M`, and the eleven-speed LRC gives `M(P)>=1/12`.
Moreover `gcd(P)=d`, `gcd(d,w)=gcd(d,13h)=1`, so the deck order in B is
`D=d`.  Formula (4) gives

```text
M(P union {w}) >= min(1/12,(d-1)/(2d)) = 1/12,
```

because `d>=2`.  The same gcd computation makes the full family primitive.
When `d=1 (mod 13)`, the unchanged coordinates have residues `r` and the
defect has residue `j`, so the residues are exactly `1,...,12`.  This proves
D.  ∎

## Relation to THM-759 and the sporadic branch

THM-759 uses one safe component near a core maximizer and only its guaranteed
width; it yields a ratio bound.  Part A retains the exact component length
*and* its tooth phase.  Part B then uses all translation lifts of that same
component simultaneously.  The resulting obstruction is scale-uniform:
increasing the common core scale creates more deck sheets, rather than a
harder unbounded-height case.

This does not yet eliminate primitive super-lonely eleven-cores, so it does
not by itself prove the full `n=12` rigidity.  It does remove, uniformly and
without a height cutoff, every imprimitive core and every proposed
high-winding construction whose max peel leaves a common scale.

Two guardrails show why the hypotheses and the minimum in (4) should not be
strengthened casually.

- The strict core gap `L<M(P)` is essential.  The primitive family
  `{15,26,52,...,312}` consists of the imprimitive tight core
  `26*{1,...,12}` plus a safe runner `15`; both the core and the full family
  have `M=1/13`, although `26` does not divide `15`.
- The `min` in (4) is sharp.  For `P=2*{1,...,11}` and `w=1`, one has
  `d=D=2`, `M(P)=1/12`, and
  `M({1,2,4,...,22})=1/12=min(1/12,1/4)`.  The deck forces one good killer
  phase, but adding the killer can never raise `M` above the core value.

## Tournament reading and assumption challenge

Runner vertices obscure the proof: the preserved predicate is simultaneous
coverage of a strict safe component and all of its deck translates.  Use the
`D` deck lifts as vertices instead.  The pairwise observable is the signed
circular displacement of their killer phases; reversing the circle is the
switch/gauge, and antipodal ties (when `D` is even) are resolved by the deck
order `0,1,...,D-1`, the tie Hamiltonian path.  For odd `D` the phase-order
tournament is the cyclic tournament (one SCC, equal score sequence, directed
cycles); for even `D` the antipodal tie pairs record the same wrap.  A single
tooth of length `<1/2` would linearize all vertices into one short arc, while
the cyclic fingerprint says that the deck wraps the circle.  Inequality (9)
is the metric strengthening of this tournament obstruction.

Alternate vertices considered were runners, tooth labels, safe components,
and proof obligations.  Deck lifts are the smallest choice preserving the
needed LRC predicate.  Quotienting them merely to runner residues destroys
the component-width tax in (2); quotienting to the gcd deck preserves the
coverage obstruction while discarding irrelevant absolute height.

## Computational audit

The companion exact census
`04-computation/lrc13_n12_tight_census_codex_S3.cpp` checks bounded speed boxes
using all rational heights `q<=2N`.  This range is exact because a peak of the
lower envelope occurs at a pair crossing (`q=v_i+v_j` or `|v_i-v_j|`) or at a
self-cusp (`q=2v_i`).  The self-cusp family is essential; the older S107/S108
local scripts omitted it from their denominator sets.  Re-running their
bounded domains with the corrected engine did not change their conclusions.

The census through `N=30` finds only `{1,...,12}` among `86,492,770` primitive
twelve-subsets.  This is exact finite evidence, not the uniform theorem: the
uniform contribution of THM-765 is the deck/divisibility reduction above.
