# Parabolic cusp, golden ray, Boolean planes, and the missing owner

**Synthesis date:** 2026-08-12.  This reflection routes proved objects and
open transfers; theorem files, not this narrative, are the truth sources.

## Current board

| lane | exact carrier | status | theorem boundary |
|---|---|---|---|
| parabolic branch | consecutive parameters `(t,t+1)` and the fixed Farey cusp `(1,1)` | **PROVED** | [THM-3334](../01-canon/theorems/THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md) |
| golden branch | Fibonacci parameters, three ancestry rays, affine cocycle/current no-go | **PROVED** | [THM-3339](../01-canon/theorems/THM-3339-fibonacci-three-ray-berggren-transplant-and-moving-owner-obstruction.md) |
| equal-plane collisions | Gaussian factor-choice torsors `F_2^r/<1>` | **PROVED** | THM-3334 |
| internal four-box | `W(m,n)=(n-m,m,n,n+m)` with a signed matching current | **PROVED** | THM-3339 |
| square/triangular selector | Pell-8/Markov compiler selecting square even legs on the U-spine | **PROVED** | [THM-3335](../01-canon/theorems/THM-3335-square-triangular-pell-markov-pythagorean-selector.md) |
| spine scalar intersections | square `C_t`, triangular `Q_t`, branch transplant, selector fibres | **PROVED** | [THM-3341](../01-canon/theorems/THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors.md) |
| Gaussian content curvature | primitive product, charge `H^1`, folded XOR weights, source groupoid | **PROVED** | [THM-3336](../01-canon/theorems/THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation.md) |
| prime-XOR ancestry | flat endpoint-path groupoid with source-dependent costs | **PROVED** | [THM-3345](../01-canon/theorems/THM-3345-prime-xor-ancestry-path-groupoid-and-source-dependent-berggren-cost.md) |
| LRC spectral closure | owner, phase, clock, endpoint word, global exit | **OPEN** | THM-3334/3336/3339/3345 remove no row |

## One interval, two genuinely different rays

Both constructions live on reduced fractions `m/n` in `(0,1)`, but their
dynamics have different Jordan type.

The THM-3334 spine uses

```text
g=[2 -1;1 0],          g(t,t+1)=(t+1,t+2),
(g-I)^2=0.
```

It is parabolic.  Its fractions `t/(t+1)` approach the rational fixed cusp
`1`, and consecutive vertices form a fan of Farey triangles about `(1,1)`.
Gaussian squaring turns linear parameter drift into quadratic triple growth.

The THM-3339 ray uses

```text
G=[0 1;1 1],           G(F_k,F_(k+1))=(F_(k+1),F_(k+2)).
```

It is hyperbolic up to orientation, with eigenvalues `phi` and `-1/phi`.
Its fractions approach the quadratic irrational `1/phi`.  Three Fibonacci
steps equal two Berggren parameter steps,

```text
G^3=AB,
```

except that the odd/odd residue lane must first pass through the primitive
normalizer `T`, where `T G^3 T^(-1)=CB`.  This produces the three ancestry
rays

```text
(BA)^r,             A(BA)^r,             C(BC)^r.
```

The parabolic and golden loci meet only at `(1,2)` and `(2,3)`, giving
`(3,4,5)` and `(5,12,13)`.  The next parameters `(3,4)` and `(3,5)` already
separate.  The initial coincidence is not a common infinite branch.

## The square-triangular theorem is a sparse parabolic clock

THM-3335 does not insert a fourth ancestry tree.  It selects the U-spine
indices `t` for which

```text
T_t=q^2,             x=2t+1,             s=2q,
x^2-8q^2=1.                                             (1)
```

At those and only those positive selector indices,

```text
P_t=(x,s^2,s^2+1),
W=s^2+2,
Q=2(s^2+1)+1=2W-1=x^2+2.                               (2)
```

Thus the THM-3334 plane scalar `Q` and THM-3335 skew-EW candidate scalar `W`
are losslessly affine-related on the selector.  The first rows are

```text
t=1:     P=(3,4,5),          W=6,       Q=11,
t=8:     P=(17,144,145),     W=146,     Q=291,
t=49:    P=(99,4900,4901),   W=4902,    Q=9803.
```

The third row is the `29*169=4901=70^2+1` cannonball address.  In the
Berggren tree these rows sit at U-depth `t-1`, so the selector is an
exponentially sparse clock on one parabolic branch.  Its two-state matrix has
eigenvalues `3+-2sqrt(2)`.  This is distinct from both the unipotent U-step
and the golden matrix with eigenvalues `phi,-1/phi`.

There are therefore three typed dynamics, not one overloaded “Pell tree”:

| object | state/action | limiting or selected feature |
|---|---|---|
| full U-spine | unipotent `g`, `(g-I)^2=0` | rational cusp `1`, every depth |
| square-triangular clock | `x^2-8q^2=1` | sparse depths with square even leg |
| Fibonacci/golden ray | `|n^2-mn-m^2|=1` | irrational endpoint `1/phi`, three ancestry rays |

The first three Pell-style rows in this discussion are not different
quadratic fields: `Q(sqrt(8))=Q(sqrt(2))`.  The square-`C_t` negative-Pell
selector and THM-3335 clock advance under `3+2sqrt(2)`, while the odd-parity
norm-17 `Q_t` branches advance under its square.  Their norm/parity cosets and
selected predicates differ; their ambient Pell field does not.

THM-3341 closes the remaining square/triangular table.  For `t>=1`, `C_t`
is square exactly on

```text
t=3,20,119,696,...,
```

is never triangular, and `Q_t` is never square.  The triangular `Q_t` values
are exactly two norm-17 Pell orbits beginning

```text
t=6,23,221,798,7524,... .
```

Each THM-3335 square-triangular row joins two consecutive square-`C_t` roots
as its fixed-two Markov coordinates.  Gaussian squaring then sends the middle
Berggren ray to U-depths

```text
2,19,118,695,4058,...,
```

with variable drift.  Most importantly, Boolean fibre rank remains unbounded
inside this sparse square selector.  Its first collision is
`C_696=985^2`, where the two ancestry words are `U^695` and
`UUUUUDADUDDU`.

## The four-versus-six puzzle is a typed diagram

There are two unrelated `K4` vertex sets.

1. The **internal** `K4` has the four entries of one window.  Its six edge
   products encode one Pythagorean triple only through asymmetric operations:
   one edge, twice an opposite edge, a sum, and a signed difference.
2. The **external** `K4` at `c=1105` has four distinct equal-hypotenuse
   parents.  Its six edges are coloured by the three Gaussian prime-XOR
   directions.

Both have three perfect matchings because every `K4` does.  That shared
incidence does not identify their vertices, operations, or ancestry.

On the Fibonacci window three finite shadows coexist:

- the four vertex values give the same transitive `T4` for every strict
  window and therefore forget Cassini sign, content, and ancestry ray;
- the six edge products give a transitive `T6` whose middle arc alternates;
- oriented Farey reduction gives a `C6` of all six transitive `T3` orders on
  the three matching channels.

The alternating `T6` arc is the isolated edge swap `(03 12)`.  It is odd,
whereas every edge permutation induced by `S4` is even.  So it is a ranking
sidecar, not a relabeling of the underlying `K4`.  At the quotient level each
adjacent matching reflection has four lifts.  Exactly four lift pairs fix an
owner, one for each vertex; the other twelve generate transitive `S4`.  The
six-state loop sees none of this `V4` coordinate.

## Gaussian multiplication supplies charges, not a global XOR action

THM-3336 closes the multiplication question left by the collision torsor.
Signed primitive Gaussian integers form an abelian group under multiply and
remove-content.  Its finite torsion is `C8`, generated by `1+i`; after
quotienting the four units, its exact coordinates are

```text
G = C2 direct-sum (direct sum over p=1 mod4 of Z),
tau in H^1(G;F_2),              lambda_p in H^1(G;Z).    (3)
```

Here `tau` is the odd/odd coordinate and `lambda_p` is the signed choice of a
Gaussian prime above `p`.  Odd content is the cancellation shadow

```text
v_p(k_odd(z,w))
  =(|lambda_p(z)|+|lambda_p(w)|-|lambda_p(z)+lambda_p(w)|)/2, (4)
```

an orientation-dependent coboundary, not a new `H^1` class.  This is the
explicit cohomological dictionary promised by the earlier XOR picture, but
its base is the associate multiplication group.  It is not cohomology of the
Berggren tree, the external `K4`, an LRC ancestry base, or a JC complement.

On a fixed hypotenuse fibre `X_c=F_2^r/<1>`, a nonzero XOR direction represented
by a prime subset `S` has the canonical folded weight

```text
K_c(S)={P_S,c/P_S}.                                      (5)
```

At `c=1105` these are `{5,221}`, `{13,85}`, and `{17,65}`.  They distinguish
the three perfect matchings and are invariant under all four translations,
so they choose no owner or tournament orientation.  Multiplication does not
descend through conjugation: choosing the other lift can change the output
hypotenuse.  A section turns it into source-dependent arrows, and even at
`c=65` the two allowed fixed-fibre arrows merely swap when the section swaps.
The correct object is therefore a weighted groupoid, not a global Boolean
action.

THM-3345 transports that groupoid into actual ancestry at `c=1105`.  The two
path costs in the prime matchings are

```text
p=5: {8,27},             p=13: {5,24},             p=17: {26,9}. (5a)
```

Thus even prime plus folded weight does not determine ancestry cost.  The
lawful arrow is the endpoint coboundary `P(x,y)=w_x^(-1)w_y`: it is
source-dependent but flat, and every external-`K4` loop dies in the ambient
Berggren tree.  Freezing one basepoint arrow per colour creates nonzero words
only by composing arrows with mismatched sources.  Those defects are not a
tree `H^1` class.

The first eight-parent record at `c=99905` expands the path costs to `238` and
the depth jumps to `216`, but the sharper result is infinite.  For `s>=1`, the
U-spine row at `t=25s+1` has a prime-5 partner with

```text
U^(25s) |--> DD U^(s-1) ADD,
depth jump=24s-4,                  path cost=26s+4.        (5b)
```

This is a rational transduction on that unary sublanguage.  Combining
`t=1 mod 25` with the split-prime root conditions by CRT proves that Boolean
rank and ancestry dispersion are simultaneously unbounded.  A uniform
source-reading transducer across every Gaussian fibre remains open.

This also sharpens the proposed degree-monoid analogy.  Raw Gaussian norms do
multiply.  Primitive matrix normalization does not preserve that grading:

```text
A=B=[1 -1;1 1],       Delta(A)=Delta(B)=2,
AB=[0 -2;2 0],        Delta(AB)=1.                        (6)
```

The content cocycle remembers the lost scalar.  Thus the raw family may be
degree-graded while the primitive action is cocycle-graded.  The ramified
element `1+i` is one of exactly eight signed universal grade preservers; its
odd/odd normalizer is, up to signed unit/leg gauge, precisely the third
Fibonacci-ray correction already seen in THM-3339.

## The two appearances of `-4` point in different directions

For the parabolic spine,

```text
c_t=t^2+(t+1)^2=2t^2+2t+1,        disc_t(c_t)=-4.
```

This fixed quadratic discriminant says that a prime divides some `c_t`
exactly when `-1` is a square modulo that prime.  CRT then builds arbitrarily
large Gaussian factor-choice ranks and hence unbounded equal-plane fibres.

For the known three-dimensional Keller counterexample, THM-1310 has

```text
disc_x(N)=-4 Q^2 L.
```

There the variable square class `-L`, not the visible constant `-4`, carries
the cubic `S3` monodromy.  The square factor `-4Q^2` cannot classify the
counterexample family or identify it with the Pythagorean fibres.  The golden
Pell form has discriminant `5`, adding a third, distinct arithmetic carrier.

## What the LRC lane actually inherits

THM-3334 contributes a canonical fixed-cusp Farey fan, exact Lorentz lifts,
and large Boolean side fibres.  THM-3336 adds exact content weights and charge
classes, but its determinant-gate reversals compare separately saturated
decks on different rational planes; gate failure is only loss of a sufficient
certificate.  THM-3339 contributes a branch word, a six-state channel order,
and a complete finite owner-lift obstruction.  These are useful coordinates,
but the LRC target needs a physical signed owner and phase on one typed
ancestry base.  Reduction modulo two retains the order of three channels and
destroys the integral Gram, height, endpoint history, and owner.  The exact
theta-slaved `r=5` contraction is positive evidence on a typed row, not a
bridge that restores those lost coordinates.

The branch-word question now has an exact answer.  After a lawful channel
frame and one displacement calibration, the branch matrices define a
`V4`-valued cocycle with affine image `D4`.  Its owner sequence is
`0,p,p,r,q,q`, and a moving gauge flattens it.  But the full signed current
has trivial translation stabilizer: every nonzero correction changes either
`(a,b)` to `(b/2,2a)`, flips the antisymmetric `c`, or both.  The current's
four signatures decode the missing translation rather than erase it.

The next lawful question is therefore not “is the puzzle a tournament?” or
“can the owner be flattened?” It is:

```text
Can the 24-state residue-order x signed-current-orbit bundle be transported
to one physical LRC ancestry base without forgetting owner, phase, or sign?
```

A static gauge is insufficient, and THM-3339 proves that even the derived
branch correction cannot preserve the fixed current while flattening the
owner.  A physical use must retain the current orbit as part of its state.

## Research queue opened by the synthesis

1. The square/triangular intersection table is now closed by THM-3341.  Use
   its variable-length Gaussian branch transplant to test whether ancestry
   depth dispersion in square-selector fibres is unbounded quantitatively.
2. THM-3336 closes the content cocycle and THM-3345 proves a flat
   source-dependent ancestry lift, its first two record tables, and a unary
   prime-5 transducer.  Decide whether one uniform source-reading transducer
   works across all Boolean fibres.
3. Build the exact `6*4` residue/order-current bundle, including the odd/odd
   leg-order normalizer, and test whether its `D4` holonomy has a physical
   LRC owner/phase interpretation.  The flatten-and-preserve route is closed.
4. THM-3345 proves depth and path-cost dispersion are unbounded simultaneously
   with Boolean rank.  Quantify the full cost distribution and ask whether
   other fixed split-prime toggles admit comparable regular sublanguages.
5. For LRC, retain the fixed-cusp integral ancestor and ask for one explicit
   owner/phase/current transport on the canonical typed row.  The Boolean
   rank or residue hexagon by itself is not a spectral-closure certificate.

The useful meta-pattern from incoming THM-3340 is methodological only: a
coarse quotient defect may be repaired by a delayed donor when an explicit
injection proves capacity.  In the present lane the odd/odd normalization
ray is a tempting donor, but THM-3339's stabilizer no-go says no donor can
flatten the owner while leaving the same signed current fixed.  A different
target bundle would need its own capacity theorem.  This analogy is not a
theorem dependency.
