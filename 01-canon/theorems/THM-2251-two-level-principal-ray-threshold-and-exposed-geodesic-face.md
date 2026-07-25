---
id: THM-2251
title: "Two-level principal-ray threshold and exposed geodesic face"
status: >
  PROVED / HOSTILE-AUDITED (abstract threshold, sharpness model, and
  exposed-geodesic lemma) + CITED APPLICATION (the conditional
  10_6/trefoil ray). In an integer-valued nonexpansive commutative metric
  monoid, a
  translation-invariant floor only one unit below the raw length makes every
  principal-ray continuation response a two-level Heaviside sequence, with
  one threshold N in N_0 union {infinity}. The threshold can be any prescribed
  finite integer or infinity even for a cancellative finitely generated
  monoid, so the abstract kernel laws cannot force a one-copy catalyst.
  Every floor-attaining unit geodesic lies simultaneously on the exposed
  faces of all tight dominated pseudometrics and additive 1-Lipschitz
  potentials. CITED APPLICATION: for K=10_6 and the trefoil T in the
  Brittenham--Hermiller certificate, d_G(K#nT,nT) is 3 before one threshold
  and 2 after it; the cited root path does not determine that threshold.
source: codex-2026-07-25-two-level-principal-ray
depends_on:
  - THM-2176-gordian-continuation-profile-and-interaction-cocycle
  - THM-2191-catalytic-localization-of-the-gordian-metric
related:
  - THM-2220-fixed-context-stable-response-and-catalyst-complexity
  - THM-2242-tournament-complement-transport-and-knot-kernel-green-rigidity
  - THM-2248-higher-interaction-defect-complex-and-tropical-trace-spectrum
external:
  - "Mark Brittenham and Susan Hermiller, Unknotting number and connected sums: The knots 4_1 and 5_1, arXiv:2601.18757v1."
  - "Brendan Owens and Saso Strle, Immersed disks, slicing numbers and concordance unknotting numbers, Communications in Analysis and Geometry 24 (2016), 1107--1138, DOI 10.4310/CAG.2016.v24.n5.a8, arXiv:1311.6702v3."
---

# THM-2251 -- the conditional `10_6` ray has one cut

THM-2191 proves that a catalyst response is antitone under adding common
summands. On a general catalyst monoid this leaves a large upper ideal.
At the conditional `10_6` boundary, however, the Owens--Strle floor and the
raw unknotting bound are only one integer apart. This collapses every fixed
principal ray to one cut.

The result is exact but deliberately does not pretend to locate the cut.
A hostile model below proves that no formal metric-monoid argument can force
the cut to occur at the first copy, or at any finite copy.

## 1. Abstract two-level ray theorem

Let `(M,+,0)` be a commutative monoid and let `d` be an integer-valued metric
satisfying joint nonexpansiveness

```text
d(a+c,b+e) <= d(a,b)+d(c,e).                         (1)
```

Let `delta` be a translation-invariant pseudometric with

```text
delta(a,b) <= d(a,b),
delta(a+c,b+c)=delta(a,b).                           (2)
```

Fix `x,y in M` and an integer `m>=1`, and suppose

```text
delta(x,0)=m-1,
d(x,0)<=m.                                           (3)
```

For `n>=0`, define the principal-ray response

```text
r_n=d(x+n y,n y).                                    (4)
```

Translation invariance of `delta`, domination (2), and (1) give

```text
m-1=delta(x+n y,n y)
   <=r_n
   <=d(x,0)
   <=m.                                              (5)
```

The common translation by `y` gives

```text
r_(n+1)<=r_n.                                        (6)
```

Since `d` is integer-valued, (5) leaves only the two values `m-1` and `m`.
Put

```text
N=min{n>=0:r_n=m-1},                                 (7)
```

with `N=infinity` when the set is empty. Equation (6) proves the exact
normal form

```text
r_n = m                         for n<N,
r_n = m-1                       for n>=N.             (8)
```

Equivalently, in the continuation-kernel notation of THM-2176,

```text
P_x(n y,n y)=m - 1_[n>=N].                           (9)
```

Here the indicator is zero when `N=infinity`. The threshold orientation is
therefore fixed: large contexts can only move the response down. Moreover,

```text
d(x,0)=m-1  => N=0,
d(x,0)=m    => N>=1 or N=infinity.                  (10)
```

If `N<infinity` and a path at `n=N` realizes the lower value, adding
`(n-N)y` to every vertex gives a path of no greater length at every later
`n`: common translation may weakly shorten an individual edge. The floor
in (5) forces the translated path's total length back to exactly `m-1`, so
it is geodesic. Thus `N` is the unique primitive context size on this ray;
all later witnesses are automatic translates.

More generally, if the integer floor and ceiling in (3) are `a<=b`, the
same proof shows that `r_n` is a nonincreasing integer sequence in
`{a,...,b}` and has at most `b-a` strict drops. The one-unit case is special
because the full infinite ray has only one missing coordinate, the cut (7).

## 2. Equality exposes every geodesic

The floor supplies more than the endpoint lower bound.

> **Exposed-geodesic lemma.** Let `q<=d` be any pseudometric. Suppose
> `z_0,...,z_r` is a unit `d`-geodesic:
>
> ```text
> d(z_(i-1),z_i)=1,
> d(z_0,z_r)=r.
> ```
>
> If `q(z_0,z_r)=r`, then
>
> ```text
> q(z_(i-1),z_i)=1,
> q(z_0,z_i)=i,
> q(z_i,z_r)=r-i                                   (11)
> ```
>
> for every `i`.

Indeed,

```text
r=q(z_0,z_r)
 <=sum_i q(z_(i-1),z_i)
 <=sum_i d(z_(i-1),z_i)=r.                          (12)
```

Every summand in the middle is at most one, so equality forces every
summand to equal one. For a fixed `i`,

```text
r<=q(z_0,z_i)+q(z_i,z_r)<=i+(r-i)=r,                 (13)
```

which proves the prefix and suffix statements.

There is a simultaneous linear version. Let `phi:M->R` be additive and
1-Lipschitz for `d`. Suppose `r_n=r>0` and `|phi(x)|=r`. Write
`epsilon=sign(phi(x))`. On every unit geodesic

```text
z_0=x+n y, z_1,...,z_r=n y,
```

each edge must satisfy

```text
epsilon(phi(z_(i-1))-phi(z_i))=1,
phi(z_i)=phi(n y)+epsilon(r-i).                      (14)
```

To see this, sum the `r` quantities on the left. Their sum is
`epsilon phi(x)=r`, while each is at most one by the Lipschitz property.
The argument applies at once to every dominated pseudometric and every
additive potential that is tight on the same endpoint pair. Hence a
floor-attaining geodesic lies in the intersection of all their exposed
faces, not merely in one scalar level set.

This is the metric analogue of a tight Kantorovich potential. It is a
necessary carrier for a shortcut, not a claim that the exposed face is
geodesically convex or nonempty.

## 3. The cut can occur arbitrarily late

The normal form (8) is sharp. Fix any

```text
N in N_0 union {infinity}.
```

Take the cancellative finitely generated commutative monoid

```text
M=N_0^2
```

and make it an undirected graph with the following unit edges:

```text
(a,b) -- (a+1,b),                                    (H)
(a,b) -- (a,b+1),                                    (V)
(a+3,b) -- (a+1,b)       for b>=N,                  (S_N)
```

where the third family is omitted when `N=infinity`. Let `d_N` be the
graph metric. These are three translation-orbits of edges: `(S_N)` is the
orbit of

```text
(3,N)--(1,N)
```

when `N` is finite. Translating an edge by any monoid element remains an
edge. Consequently common translation is 1-Lipschitz, and a triangle
through the intermediate point `b+c` proves (1).

Define

```text
delta((a,b),(c,e))=ceil(|a-c|/2).                    (15)
```

This is a translation-invariant pseudometric. The triangle inequality
follows from

```text
ceil(|a-c|/2)
 <=ceil((|a-b|+|b-c|)/2)
 <=ceil(|a-b|/2)+ceil(|b-c|/2).
```

Each horizontal or shortcut edge has `delta`-length one and each vertical
edge has `delta`-length zero. Applying the pseudometric triangle inequality
along a graph path and then minimizing its length gives `delta<=d_N`.

Set

```text
x=(3,0), y=(0,1), m=3.                              (16)
```

Then `delta(x,0)=2`. Three horizontal edges show
`d_N(x+n y,n y)<=3` for every `n`. If `N` is finite and `n>=N`, the path

```text
(3,n)--(1,n)--(0,n)                                  (17)
```

has length two, and the floor makes its length exact.

If `n<N`, no two-edge path can change the first coordinate by three.
Such a path would need one shortcut, changing that coordinate by two, and
one horizontal edge, changing it by one. Both preserve the second
coordinate, but no shortcut exists at height `n<N`. A vertical edge cannot
replace either of them. Thus the three-horizontal path is geodesic. With no
shortcuts, the same is immediate for `N=infinity`. Therefore

```text
d_N(x+n y,n y)=3       for n<N,
d_N(x+n y,n y)=2       for n>=N.                    (18)
```

Every possible threshold, including infinity, occurs under the abstract
hypotheses and even with cancellativity and finitely many edge orbits. In
particular:

```text
antitonicity does not imply N=1;
finite generation does not imply N<infinity;
the full-kernel composition law does not locate N.  (19)
```

Any proof that the `10_6` cut is one must therefore use a knot-specific
coordinate not present in the bare metric-monoid axioms.

## 4. Exact `10_6` trefoil-ray reduction

Let

```text
K=10_6
```

and let `T` be the trefoil chirality used in the
Brittenham--Hermiller three-crossing certificate. Define

```text
r_n=d_G(K#nT,nT)=P_K(nT,nT).                        (20)
```

The source audit inherited from THM-2191 gives

```text
c*(K)=2,
u(K) in {2,3}.                                       (21)
```

Owens--Strle's immersed-concordance distance `d_*` is exactly
connected-sum-translation invariant and is bounded above by Gordian
distance. Hence

```text
2=d_*(K,U)
 =d_*(K#nT,nT)
 <=r_n
 <=u(K)
 <=3.                                                (22)
```

Apply Section 1 with `delta=d_*` and `m=3`. There is a unique

```text
N_T(K) in N_0 union {infinity}
```

such that

```text
r_n=3                         for n<N_T(K),
r_n=2                         for n>=N_T(K).          (23)
```

This completely classifies the values on the trefoil ray without claiming
to compute its cut:

```text
u(K)=2
  => N_T(K)=0 and r_n=2 for every n;

u(K)=3
  => N_T(K)>=1 or N_T(K)=infinity;

u(K)=3 and T catalyzes K
  <=> N_T(K)=1;

u(K)=3 and some trefoil power catalyzes K
  <=> N_T(K)<infinity.                               (24)
```

If the last line holds, (22) and THM-2191's localization give
`u_cat(K)=2`. Its amplification argument gives `u_hash(K)<=2`, while the
cited additive signature floor gives the reverse inequality. Hence

```text
u_cat(K)=u_hash(K)=2.                                (25)
```

If `N_T(K)=infinity`, only the trefoil ray is excluded. Another knot may
still catalyze `K`; (23) does not compute the infimum over all contexts.

## 5. Root bypass is not the catalytic diagonal

Brittenham--Hermiller exhibit a length-three crossing-change path from
`K#T` to the unknot. It proves exactly

```text
u(K#T)<=3,
```

which is the pinned root entry

```text
P_(K#T)(U,U)<=3.                                     (26)
```

The proposed one-trefoil catalytic shortcut is instead

```text
d_G(K#T,T)=P_K(T,T)=2.                               (27)
```

Appending one unknotting change for `T` proves

```text
(27) => (26),
```

but the reverse implication is invalid: the published three-step path has
the wrong target. In the possible branch `u(K)=3`, (26) makes the labelled
pair `{K,T}` a positive-defect nonface in the sense of THM-2248, because

```text
u(K)+u(T)-u(K#T)>=3+1-3=1.
```

Nevertheless,
(23) still permits either value two or three at the catalytic diagonal.
Thus the higher root-defect complex and the principal-ray threshold are
different traces of the continuation kernel.

This distinction also prevents the named intermediate knots on the
published root path from being silently promoted to a two-edge path ending
at `T`.

## 6. What a two-crossing witness must look like

Suppose `r_n=2`, and let

```text
K#nT = L_0 -- L_1 -- L_2 = nT                      (28)
```

be any two-crossing geodesic. Equality in (22) and the
exposed-geodesic lemma force

```text
d_*(L_0,L_1)=d_*(L_1,L_2)=1.                        (29)
```

The classical signature is additive and changes by at most two under one
crossing change. Since the cited data give

```text
|sigma(K)|/2=2,
```

the potential `sigma/2`, after choosing its sign, is also tight. Equation
(14) then forces both crossing changes to be signature-extremal and to have
the same signature direction. The intermediate signature is not free:

```text
sigma(L_1)
 =n sigma(T)+sigma(K)/2.                             (30)
```

Equivalently, a candidate intermediate must lie simultaneously in the two
unit Gordian spheres, the two unit `d_*` faces, and the signature midpoint
hyperplane. Every later witness obtained by adding trefoils remains on
these same translated faces.

These conditions are necessary, not sufficient. The equality `c*(K)=2`
records a two-double-point immersed-concordance class after common summands
are cancelled. Gordian distance additionally demands a level-preserving
movie of crossing changes with the specified endpoint knot. That endpoint
and movie coordinate is exactly what the translation-invariant floor
forgets. Likewise, (30) does not construct the intermediate knot.

## 7. Failure boundary and next decisive test

The theorem proves:

```text
infinite trefoil-ray kernel data -> one threshold N_T(K);
floor-attaining path            -> simultaneous exposed-face path;
published root shortcut         -> no value of N_T(K);
abstract monoid laws            -> no bound at all on N_T(K).   (31)
```

It does not decide `u(10_6)`, exhibit a positive knot catalyst, or prove
that `N_T(K)` is finite. The sharp knot-specific test is now:

```text
find or exclude a knot L with

d_G(10_6#T,L)=d_G(L,T)=1,
d_*(10_6#T,L)=d_*(L,T)=1,
sigma(L)=sigma(T)+sigma(10_6)/2.                    (32)
```

For higher trefoil powers replace `T` by `nT`; the first successful `n` is
exactly `N_T(10_6)`. A partial-cube coordinate, tournament orientation, or
root-defect label is useful only if it preserves this simultaneous
endpoint-sensitive face. None follows from the abstract continuation
kernel alone.

QED.
