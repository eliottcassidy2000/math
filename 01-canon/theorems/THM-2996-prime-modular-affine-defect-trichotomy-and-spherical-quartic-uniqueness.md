---
id: THM-2996
title: "Prime-modular affine-defect trichotomy and spherical quartic uniqueness"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For the displayed integral C2*C3 -> S3 module reduced modulo any
  prime p, translation-gauge classes of affine lifts form F_p, with defect
  lambda=2b2-b1 and cusp translation (st)^2=lambda e1.  A nonzero class has
  transitive image F_p^2 semidirect S3, order 6p^2, and cusp order 2p.
  The resulting triangle quotient is spherical only at p=2, Euclidean at
  p=3, and hyperbolic for p>=5; p=2 is uniquely characterized by the
  faithful transitive S3 action on nonzero module vectors and gives S4.
  At the bad prime p=3 the dual lattice has two, not one, cohomology
  coordinates and a separate four-stratum affine census.  The p=13
  consequence is an abstract 169-point affine shape only, not an LRC
  target/Bockstein or physical-current identification.
source: codex-prime-modular-affine-defect-2026-07-31
depends_on: []
related:
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2743-quartic-root-affine-plane-and-six-sheet-flag-rigidity
  - THM-2768-oriented-modular-triangle-quotients-and-keller-lrc-boundary
  - THM-2950-three-conjugate-pair-v-four-torsor-and-quartic-resolvent-frame
  - THM-2975-modular-six-sheet-schreier-orbifold-and-cusp-width-boundary
  - THM-2992-signed-quartic-edge-block-discriminant-parity-and-keller-owner-line-boundary
  - THM-2993-quartic-signed-edge-star-triangle-cube-and-derivative-square-reembedding
script: 04-computation/prime_modular_affine_defect_trichotomy_thm2996.py
output: 05-knowledge/results/prime_modular_affine_defect_trichotomy_thm2996.out
script_sha256: 0558163e210043a3e87e4dc818bd2356ff072cedff4fa053f9e644258be462e3
output_sha256: 4f15833270e6450dfde29de9839ca98d4cd06bebe2e8ec6562e34a584abc4588
hash_basis: LF-normalized bytes
---

# THM-2996 -- prime-modular affine defect and spherical quartic uniqueness

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**  This file is not a proved dependency until that audit promotes it.

## 1. Inheritance and statement

[THM-2595](THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go.md)
proves the binary `V4` case: the two free-factor restrictions are trivial,
but their common translation gauge may be obstructed, producing either a
point-stabilizing `S3` or transitive `S4`.  The closest hostile is already
visible at the next primes: the same compatibility class produces finite
affine quotients on `p^2` points without producing a quartic.  The least-used
sidecar is the integral root/weight-lattice distinction, whose determinant
three becomes singular exactly where the ternary factor resonates.

Let

```text
G = <s,t | s^2=t^3=1> = C2*C3,

        [0  1]              [0 -1]
S   =   [1  0],       T  =  [1 -1].                         (1)
```

Reduce `(1)` modulo a prime `p`, put `V_p=F_p^2`, and let `rho(s)=S`,
`rho(t)=T`.  Then `S^2=T^3=(ST)^2=1`, and the linear image is `S3` for
every prime.  An affine lift is

```text
s(v)=Sv+a,                    t(v)=Tv+b.                    (2)
```

Actions are on the left, and throughout `st=s∘t` (first apply
`t`, then `s`).  Reversing this convention conjugates the conclusions but
exchanges the displayed `e_1/e_2` cusp coordinate.

Up to translation conjugacy, all such lifts form

```text
H^1(G,V_p) = F_p.                                            (3)
```

After gauging `a=0`, the exact class coordinate is

```text
lambda(b)=2 b_2-b_1.                                       (4)
```

If `lambda=0`, the affine image is a linear `S3` complement and fixes a
point.  If `lambda!=0`, the image is

```text
V_p semidirect S3,
|image|=6p^2,
one orbit of size p^2,
ord(st)=2p.                                                 (5)
```

Moreover

```text
kappa=(st)^2=(1,lambda e_1).                               (6)
```

Thus the nonzero lift factors through `Delta(2,3,2p)`.  This triangle group
is spherical at `p=2`, Euclidean at `p=3`, and hyperbolic at every `p>=5`.
At `p=3`, (5) is a **finite affine quotient** of the Euclidean triangle
group `Delta(2,3,6)`; it is not the Euclidean group itself.

## 2. The cocycle and gauge calculation

Write a cocycle as `(a,b)=(c(s),c(t))`.  The free-product presentation gives
exactly

```text
(I+S)a=0,                    (I+T+T^2)b=0.                 (7)
```

The second matrix in `(7)` is identically zero over the integers.  The first
kernel is one-dimensional, so

```text
|Z^1(G,V_p)|=p^3.                                             (8)
```

Conjugating by `v -> v+w` adds

```text
((I-S)w,(I-T)w).                                             (9)
```

There is no common nonzero `S,T` invariant, hence `(9)` has `p^2` values.
This proves the count `(3)`.  More explicitly,

```text
ker(I+S)=im(I-S),                                            (10)
```

so one may set `a=0`.  The remaining gauge vectors lie on
`Fix(S)=<(1,1)>`, and

```text
(I-T)(u,u)=u(2,1).                                          (11)
```

The functional `(4)` annihilates exactly the line `(11)`, proving that it is
the gauge coordinate and that its kernel is the coboundary space.  A useful
uniform representative is `b=(-lambda,0)`.

For the two free factors separately,

```text
H^1(C2,V_p)=0                                  for every p,
H^1(C3,V_p)=0                                  for p!=3,
H^1(C3,V_3)=F_3.                                            (12)
```

Indeed the ternary norm is zero and `det(I-T)=3`.  At `p=3`, restriction of
the global root-module class `(3)` to the ternary factor is an isomorphism.
For `p!=3`, every global class is invisible on each factor separately: the
defect is pure binary/ternary **co-occurrence**, not local factor data.

## 3. Cusp, image, and fixed-locus proof

With `a=0`, direct multiplication gives `(6)`.  If `lambda!=0`, the normal
closure of its translation vector under the linear `S3` spans all of `V_p`.
Therefore the translation kernel is `V_p`, proving the image order and
transitivity in `(5)`.  Since `(st)^2` has order `p`, while the linear part
of `st` has order two, `st` has order `2p`.

If `lambda=0`, the cocycle is a coboundary, so a translation moves the action
to the linear complement and exposes a fixed origin.  Conversely a common
fixed point trivializes the cocycle.

There is a geometric local version.  The binary generator always fixes an
affine line of `p` points.  If `p!=3`, the ternary generator has one fixed
point, and

```text
lambda=0
  iff the ternary fixed point lies on the binary fixed line.               (13)
```

At `p=3`, the split ternary map has three fixed points, whereas every
nonzero root-module class makes it fixed-point-free.  For example `b=e_1`
has `lambda=-1`, no ternary fixed point, transitive image of order `54`, and
cusp order six.

The curvature trichotomy now follows from

```text
1/2+1/3+1/(2p)  >,=,<  1                                  (14)
```

at `p=2`, `p=3`, and `p>=5`, respectively.

## 4. Why the spherical case is uniquely quartic

For `p=2`, the following statements occur together:

1. `V_2` has four points and additive group `V4`;
2. `GL(2,F_2)=S3`;
3. its six elements act simply transitively on the six ordered bases of
   `V_2`;
4. the nonzero affine class has image
   `AGL(2,F_2)=V4 semidirect S3=S4` on four points; and
5. its marked presentation closes through the spherical group
   `Delta(2,3,4)=S4`.

This is not merely the first member of an identical family.  Suppose a
finite-dimensional vector space `V=F_p^d` over its prime field admits a
faithful `F_p`-linear `S3` action transitive on `V\{0}`.  Orbit-stabilizer
gives

```text
|V|-1 divides 6.                                             (15)
```

The prime-power possibilities for `|V|` are `2,3,4,7`.  The cases of size
`2,3,7` are one-dimensional, and their linear groups are cyclic; none
contains a faithful nonabelian `S3`.  Hence

```text
V=F_2^2                                                     (16)
```

is the unique possibility.  This proves the precise sense in which the
binary free factor, the three nonzero `V4` directions, the six modular
frames, and the quartic `S4/V4=S3` quotient coincide only in the spherical
case.

The `p=5` control is the first clean stopping boundary.  Taking `b=e_1`
gives a class invisible on both free factors, yet its image has order `150`,
is transitive on `25` points, and has cusp order `10`.  Co-occurrence and a
hyperbolic triangle quotient do not imply a quartic, a `V4` torsor, or `S4`.

## 5. The discriminant character, exactly scoped

On the `p^2` affine points, every translation is an even permutation.  A
linear reflection has `p` fixed points and `(p^2-p)/2` transpositions, while
an order-three linear element is even.  Consequently the affine permutation
sign factors through the quotient `S3` as

```text
sgn_affine = sgn_S3^epsilon,
epsilon = p(p-1)/2 mod 2.                                  (17)
```

It is the `S3` sign for `p=2` and for odd `p congruent 3 mod 4`; it is
trivial for `p congruent 1 mod 4`.

At `p=2`, `(17)` says that the sign character of the four-point `S4` action
is exactly the sign character of the induced `S3` action on the three
nonzero directions, equivalently the three perfect matchings.  Conditional
on an **actual connected generically separable four-sheet cover** over a
base of characteristic different from two whose monodromy realizes this
affine `S4` action and its matching `S3` quotient, this identifies the
quartic and matching-resolvent discriminant **square-class characters**,
because the matching quotient `S4 -> S3` has kernel `V4` contained in `A4`,
so the quartic sign is the pullback of the `S3` sign.
It does not by itself prove a literal polynomial identity.  In particular,
this abstract affine algebra constructs no cover, quartic polynomial, graph
quartic, Keller map, or resolvent morphism; those realizations and their
normalizations require the separate sidecars routed above.

## 6. The bad-prime root/dual-lattice boundary

Over the integers the contragredient module is

```text
S^vee=(S^-1)^T=S,                 T^vee=(T^-1)^T=[-1 -1; 1 0].   (18)
```

The Cartan matrix

```text
        [ 2 -1]
C   =   [-1  2],                 det C=3                      (19)
```

intertwines the root and dual modules:

```text
S^vee C=CS,                       T^vee C=CT.                 (20)
```

Thus they are isomorphic modulo every prime except three.  At `p=3`, the
root module has no nonzero global invariant and `dim H^1=1`, while the dual
module has common invariant line

```text
Fix(S^vee) intersect Fix(T^vee)=<(1,1)>                      (21)
```

and `dim H^1=2`.

In the dual module one may again set `a=0`, but now the residual gauge line
is the global invariant line `(21)` and does not move `b`.  Therefore
`b=(b_1,b_2)` itself labels the nine classes.  The correct two coordinates
are

```text
x=b_2                         (cusp coordinate),
y=b_1-b_2                     (C3-restriction coordinate),
(st)^2=(1,2x(1,1)).                                           (22)
```

The complete class table is:

| dual `b` | classes | image order | point orbits | cusp order | `t` fixed points |
|:--|--:|--:|:--|--:|--:|
| `00` | 1 | 6 | `1+1+1+6` | 2 | 3 |
| `10,20` | 2 | 6 | `3+3+3` | 2 | 0 |
| `11,22` | 2 | 18 | `3+6` | 6 | 3 |
| `01,02,12,21` | 4 | 18 | `9` | 6 | 0 |

Hence the dual action is transitive exactly when `xy!=0`.  Its two
transitive projective lines are `<(0,1)>` and `<(1,2)>`; the involution
`t -> t^2`, followed by the displayed module re-identification, exchanges
them.  The Cartan map `(19)` sends the root `H^1` line into the dual diagonal
`<(1,1)>`, which is the kernel of restriction to `C3` and is **not** the
transitive locus.  The root coordinate `lambda` must therefore never be
reused as a dual-lattice transitivity test.

This is the sharp weight-lattice hostile: the same integral hexagonal
picture has different cohomology, fixed loci, image orders, and orbit
geometry at its bad prime.

## 7. The exact `p=13` abstract control

For the root module at `p=13`, every nonzero class has

```text
169 points,
image order 6*169=1014,
cusp order 26,
one transitive affine orbit,                                  (23)
```

and both free-factor restrictions are trivial.  Thus `(23)` is an exact
abstract model of a 169-flag binary/ternary compatibility defect.

It is only an abstract shape.  No map in this theorem identifies `V_13`
with the LRC target/Bockstein plane, a common physical atom, an owner packet,
an ancestry address, a terminal Fourier current, or a common cospan.  In
particular `(23)` yields no LRC row exclusion or ledger decrement.  A lawful
application needs an explicit common module/cospan map that preserves the
physical carrier and the relevant phase.

## 8. Exact evidence and transfer contract

Run

```bash
python 04-computation/prime_modular_affine_defect_trichotomy_thm2996.py
python -O 04-computation/prime_modular_affine_defect_trichotomy_thm2996.py
```

The dependency-free companion uses only exact finite integer arithmetic and
explicit `require` gates.  For `p=2,3,5,7,11,13` it enumerates every cocycle,
coboundary, cohomology class, affine image, point orbit, fixed locus, cusp
order, translation kernel, and permutation sign, and computes the exact
ordered-basis count.  It also
checks the integral Cartan intertwiner, the complete nine-class dual `p=3`
table and restriction fibres, the two transitive projective lines and their
exchange, and the faithful-transitive uniqueness lemma.  Normal and
optimized transcripts byte-match the stored output; hashes above use
LF-normalized bytes.

Transfer from this theorem requires all of the following data:

```text
source:       a concrete binary/ternary pair of operations,
target:       the displayed root or dual module,
map:          an equivariant common affine-origin/cospan map,
preserved:    carrier, owner/address labels, and any required phase,
lost here:    physical realization and terminal amplitude,
sidecar:      a common atom or reference fixing that loss.                 (24)
```

Without `(24)`, the theorem is finite group/cohomology anatomy only.  It
does not prove a degree-four point-cap Keller map exists or fails, does not
turn a cubic resolvent into a grade-three Keller map, does not identify a
binary fraction tree or ternary triple tree with the Bass--Serre tree, and
does not canonically orient a six-vertex tournament.
