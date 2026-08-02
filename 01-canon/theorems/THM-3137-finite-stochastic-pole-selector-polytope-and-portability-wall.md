---
id: THM-3137
title: "Finite stochastic pole-selector polytope, moment separation, and portability wall"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT / AWAITING INDEPENDENT THEOREM AUDIT.
  At support (1,3), bank I2, the probability laws on the eight physical
  one-pole virtual-prefix currents that are Hasse-positive in every degree
  5 through 9 form an exact two-facet polytope with 24 vertices, all of
  minimal support two.  A 3/5--2/5 law is feasible although its same-mean
  deterministic synthetic pole fails in degrees 8 and 9.  A three-pole law
  works in both I1 and I2 at fixed support, but a two-upset Farkas certificate
  excludes every law already at support (1,2), bank I2.  Degree 10 empties
  the original one-pole selector polytope; adding an idle state is only
  trivially feasible, and every physical two-pole law is also excluded.
  This is a finite averaged-current theorem, not an original-response
  decomposition or a stochastic pole flag.
source: root/multiscale-newton-flag/product-gamma-width3-2026-08-02
depends_on:
  - THM-3110-arbitrary-anchored-product-gamma-dominant-tail-and-low-histogram-reduction
  - THM-3115-low-degree-monomial-fibre-newton-refinement-transport
  - THM-3120-row-pole-prefix-newton-flag-positivity
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
related:
  - THM-3128-active-pole-prefix-labelled-deletion-no-positive-preimage
script: 04-computation/gmc_stochastic_pole_selector_polytope_scout.py
output: 05-knowledge/results/gmc_stochastic_pole_selector_polytope_scout.out
script_sha256: 1bca02a13fe1301fd054644a2dc79d3587dd5b01d40cbab3bb2e33730fa1d97f
output_sha256: 92e41457be736551803e6438101265e8adcb9edb7ae60786cf2058d2ca306462
hash_basis: LF-normalized bytes
---

# THM-3137 -- finite stochastic pole-selector polytope, moment separation, and portability wall

**PROVED CANDIDATE + VERIFIED-EXACT / AWAITING INDEPENDENT THEOREM
AUDIT.**

THM-3120 gives a positive scalar Newton pole flag, but the associated
one-letter virtual-prefix Young current need not lie in the nonnegative Hasse
boundary cone.  THM-3128 proves that the first active deterministic failure
cannot be repaired while its labelled-deletion image is held fixed.

There is nevertheless an exact convex repair in degrees 5 through 9: average
several failed deterministic pole selectors.  The repair is sharply finite.
It retains higher moments that the scalar row polynomial cannot see, does not
extend uniformly across supports, and is impossible at degree 10 even after
the two cheapest physical augmentations.

## 1. Averaged one-pole currents

Use the THM-3120 reduced row at support

```text
(a,b)=(1,3), bank I2.                                        (1)
```

Its physical reduced-pole multiset is

```text
(8,7,6,5,5,4,4,3,3,2,2,2,1,1,1,1).                         (2)
```

Let `Phi` be the signed bank functional and `Q` the distinguished residual
alphabet, as in THM-3115/3120.  For a physical pole value
`r in {1,...,8}`, put

```text
Phi^r(f)=sum_S epsilon_S f[S-r],          Q^r=Q-r.             (3)
```

For a partition `mu` of `N`, define the raw zero-mass Young current

```text
G_N^r(mu)=Phi^r(h_N)m_mu[Q^r]
          -Phi^r(m_mu)h_N[Q^r].                              (4)
```

For a probability vector `lambda=(lambda_1,...,lambda_8)`, define

```text
G_N(lambda)=sum_{r=1}^8 lambda_r G_N^r.                       (5)
```

This is an average of eight already-defined fixed-`Q` currents.  It is not a
claim that the original scalar row response has been decomposed into eight
positive pieces.

By THM-3127, `G_N(lambda)` is a nonnegative fine-to-coarse Hasse boundary if
and only if its mass on every coarsening upset of the partition lattice
`P_N` is nonnegative.  Every current in `(4)` has total mass zero, so these
upset inequalities are the complete finite criterion.

## 2. Exact selector polytope through degree nine

For `N=5,...,9`, the numbers of coarsening upsets, including the empty upset,
are

```text
10, 27, 47, 168, 573.                                        (6)
```

Thus there are 825 exact upset rows on the eight pole coordinates.  Exactly
815 are nonzero, with degree census

```text
8, 25, 45, 166, 571,                                         (7)
```

and all 815 have distinct primitive integer normalizations.  Of these, 812
are coordinatewise nonnegative.  The only three rows with a negative
coordinate are

```text
A=( 83547971350, 15688221032,   -569502791,  -2926894294,
   -14620631052,-44124052840, -85206825873,-116934942806),   (8)

B=( 77268972364224, 14817264695349,   -301139382544,
    -1779540733230,-11604925470240,-38159896577461,
   -74842829308296,-98601655955064),                         (9)

C=(-133125723152,-49818680888, 2688744487, 18251699828,
    112450828200,398996381528,867825875007,1262919484552).   (10)
```

The middle row is redundant.  With

```text
alpha=8315869923144/9467425097,                               (11)
```

direct exact subtraction gives

```text
B-alpha A >= 0                                               (12)
```

coordinatewise; its seventh coordinate is zero and the other seven are
strictly positive.  Hence the complete selector polytope is

```text
P = Delta_7 intersect {A.lambda>=0} intersect {C.lambda>=0}.  (13)
```

This proves the facet description: all other upset inequalities are either
coordinatewise nonnegative or follow from `A.lambda>=0` and `(12)`.

## 3. The 24 vertices and minimal support

The signs in `(8),(10)` split the pole coordinates into

```text
L={1,2},                    R={3,4,5,6,7,8}.                  (14)
```

`A` is positive on `L` and negative on `R`; `C` has the opposite signs.
Therefore no simplex vertex is feasible.  On every cross edge `{i,j}` with
`i in L`, `j in R`, there are two exact feasible boundary points:

```text
A-endpoint:
  lambda_i=-A_j/(A_i-A_j),       lambda_j=A_i/(A_i-A_j),

C-endpoint:
  lambda_i= C_j/(C_j-C_i),       lambda_j=-C_i/(C_j-C_i).    (15)
```

Exact comparison gives `lambda_i(A-endpoint)<=lambda_i(C-endpoint)` for all
12 cross edges.  Exact active-set enumeration of the simplex coordinates and
the two facet rows proves that `(15)` is the full vertex list:

```text
number of vertices = 2*2*6 = 24,
support of every vertex = 2.                                 (16)
```

The companion solves every possible rank-eight active system over
`Fraction`, retains every feasible solution, and rechecks all 815 original
upset rows on each of the 24 vertices.  Hence there are no hidden
three-support vertices.  Since no deterministic selector is feasible,
support two is minimal.

On the edge `{1,8}`, the exact interval for the pole-1 weight is

```text
3077235337/5275866162
 <= lambda_1 <=
157864935569/174505650963.                                  (17)
```

In particular

```text
lambda=(3/5)delta_1+(2/5)delta_8                              (18)
```

lies strictly inside this interval.  Its exact max-flow equals its exact
demand for every `N=5,...,9`.

## 4. Same scalar mean, different Young current

At the scalar cofactor level, a law sees only its mean:

```text
sum_r lambda_r(1-rt)=1-t E_lambda[r].                         (19)
```

It is therefore tempting to replace `(18)` by the deterministic synthetic
pole `rbar=E_lambda[r]`.  That replacement is false on Young types.  For every
power sum,

```text
E_lambda p_k[X-r] = p_k[X]-E_lambda[r^k],
p_k[X-rbar]       = p_k[X]-rbar^k.                            (20)
```

For `(18)` the first moments are

```text
E[r]=19/5,             E[r^2]=131/5,
Var(r)=294/25.                                                (21)
```

The averaged law `(18)` passes all degrees 5 through 9.  The exact current
formed from the one synthetic virtual letter `19/5`, despite having the same
scalar cofactor `(19)`, has Hasse outcomes

```text
N=5,6,7: pass,             N=8,9: fail.                       (22)
```

At `N=8`, its exact flow and demand are

```text
flow   =40036754042915950613904/1953125,
demand =40038012640374254607408/1953125,                      (23)
```

and the first unresolved type is `(5,1,1,1)`, with residual

```text
1258597458303993504/1953125.                                 (24)
```

Thus the convex repair is genuinely measure-valued.  Its higher pole moments,
already its variance, survive Young refinement although the scalar row
polynomial forgets them.

## 5. Chamber portability at fixed support

The same pole law can work in both row chambers at the fixed support `(1,3)`.
Namely,

```text
lambda_1=9/500,          lambda_3=487/500,
lambda_8=4/500.                                               (25)
```

For each of banks `I1,I2` and every `N=5,...,9`, the current `(5)` formed
with `(25)` satisfies every exact upset inequality and passes exact max flow.
Its mean pole is `751/250`.

This is an arithmetic common law on the two fixed-support current banks.  It
does not by itself construct one common physical random selector or operator
decomposition.

## 6. Cross-support Farkas wall

The fixed-support result is not portable to arbitrary supports.  At support

```text
(a,b)=(1,2), bank I2,                                        (26)
```

the physical pole values are `{1,2,3,4,5}`.  Let `lambda` be any probability
law on them.  Two necessary Hasse-upset rows, from degrees 6 and 9, are

```text
R6=(2469992,-986920,-2955376,-3435376,-2426920),             (27)

R9=(-60532076544,1282401120,16315005312,
    -34474521120,-48239160000).                              (28)
```

But their positive combination is

```text
10000 R6+R9=
(-35832156544,-8586798880,-13238754688,
 -68828281120,-72508360000),                                (29)
```

strictly negative in every pole coordinate.  If both necessary inequalities
had nonnegative pairing with `lambda`, so would `(29)`, which is impossible
for a probability vector.  Therefore no probability law on the physical
one-pole currents at `(26)` is Hasse-positive through degrees 5 through 9.

As a finite full-bank control, the rule

```text
(3/5) smallest physical pole +(2/5) largest physical pole    (30)
```

passes 1,141 of the 1,150 support/bank/degree cases in the exact THM-3120
universe and fails nine.  The near-universality of `(30)` does not evade the
strict wall `(29)`.

## 7. Degree ten empties the one-pole selector polytope

The degree-nine cutoff is sharp.  At the original fibre `(1)`, degree 10 has
42 partitions and 3,588 coarsening upsets.  Six upset rows have a negative
coordinate.  The following two nested rank-tail upsets suffice:

```text
U2={mu partition 10 : length(mu)<=8}
   =P_10\{(1^10),(2,1^8)},

U1={mu partition 10 : length(mu)<=9}
   =P_10\{(1^10)}.                                           (31)
```

Their raw response rows on pole values `1,...,8` are

```text
R=(-105427024770557952,-41110096821098400,1323165981189312,
   10490179768896384,84003185593589760,328063977223727328,
   734662852146861120,1034128690334272512),                  (32)

S=(2135515427122176,1594234834827264,-128048320760256,
   -1149608741796864,-8750331832665600,-36573564507761664,
   -89737341671630016,-137769965901950976).                  (33)
```

The positive Farkas combination

```text
R+11S=
(-81936355072214016,-23573513637998496,-85365547173504,
 -2155516390869120,-12250464565731840,-74245232361650976,
 -252447906241069056,-481340934587188224)                   (34)
```

is strictly negative in all eight coordinates.  Thus **no probability law on
the eight physical one-pole currents is Hasse-positive at degree 10**.

The particular law `(18)` first fails at degree 10.  Its exact deficit is

```text
269133385522535424/5,                                        (35)
```

at first unresolved type `(7,1,1,1)`.

## 8. The two cheapest physical augmentations also fail

First add an empty/no-deletion state.  Its `U2` and `U1` coordinates are both
zero.  Pairing `(34)` with a probability law on the empty state and the eight
one-pole states forces every one-pole weight to vanish.  Hence the lazy
selector polytope is exactly the singleton all-empty law.  There is no
nontrivial lazy repair.

Next jump directly to physical two-pole deletion.  Respecting the
multiplicities in `(2)` leaves exactly 33 unordered two-pole submultisets.
Consider

```text
U3=P_10\{(1^10),(2,1^8),(2,2,1^6)}.                          (36)
```

The companion evaluates the `U3` response on all 33 states.  Every coordinate
is strictly negative; their exact range is

```text
-2197919641631883264 <= response <= -64486707449419008.      (37)
```

Therefore no probability law supported on the physical two-pole currents is
Hasse-positive at degree 10.  The next possible convex repair must mix prefix
lengths with an additional nontrivial reference state, retain extra Schur
data, or leave the positive physical-prefix simplex.

## 9. Exact verification

The companion is independent of numerical linear programming in every
promoted assertion.  It performs the following finite exact checks:

1. enumerates every partition antichain/upset in degrees 5 through 10;
2. reconstructs every one-pole current from the THM-3110 bank definitions;
3. primitive-normalizes all 815 nonzero degree-5-through-9 upset rows;
4. proves the `B-alpha A` redundancy coordinatewise;
5. solves every candidate vertex active system over rational numbers and
   rechecks all original inequalities on all 24 vertices;
6. replays max flow for the displayed positive laws and the synthetic hostile;
7. reconstructs both cross-support Farkas rows;
8. reconstructs the degree-10 nested-tail Farkas rows, the empty coordinates,
   and all 33 physical two-pole coordinates.

Fresh normal and optimized runs must byte-match the stored 63-line transcript.
The immutable LF-normalized hashes are recorded in the frontmatter.

## 10. Scope boundary

This theorem proves an exact statement about convex averages of the finite
Young currents `(4)`.  It does **not** prove any of the following:

- that the original product-Gamma scalar response or operator is a positive
  mixture with the displayed weights;
- that a measure-valued first selector extends to a consistent stochastic
  pole flag at later prefix depths;
- that a pole law may be replaced by its mean;
- that the fixed-support chamber-common law is a canonical common-root
  selector;
- Hasse positivity above degree 9 (degree 10 is explicitly impossible in the
  physical one-pole selector space);
- arbitrary-radial NC2 or the Gaussian Moment Conjecture.

The exact first failed implication is:

```text
scalar-equivalent barycentric pole
  does not imply
type-equivalent measure-valued pole selector.                 (38)
```

The strongest positive survivor is the exact two-facet convex repair through
degree nine.  The strongest stopping result is the degree-ten Farkas wall for
one-pole, lazy one-pole, and pure two-pole probability selectors.

QED (candidate pending independent theorem audit).
