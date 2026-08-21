---
id: THM-3569
title: "Universal squarefree Danielewski two-by-three weight Darboux nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  squarefree polynomial Sigma of degree at least two, a polynomial Darboux
  pair on c^2 e=Sigma(b) cannot
  have at most two nonconstant weight pieces in one output and at most three
  in the other.  After scalar removal every possible pair must have support
  sizes (2,at least 4), (at least 4,2), or (at least 3,at least 3), hence at
  least six nonconstant homogeneous pieces.  No general Darboux pair and no
  counterexample to JC(2) is claimed.
source: codex-2026-08-20-frontier-overnight; root/elliptic_arm_counterexamples, 2026-08-20
audit: >
  Two independent hostile audits rederived the support-collision trichotomy,
  every lower and upper ladder, all squarefree-root valuation gates, and the
  sharp degree-one and repeated-root boundaries.  Neither proof uses the
  coefficients or number of roots of b(b+4).  Independent enlarged support
  and valuation censuses found no extra profile.  Ordinary, optimized, and
  stored exact transcripts agree.
depends_on: []
related:
  - THM-3404-factorized-danielewski-principal-parts-and-finite-cover-obstruction
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
script: 04-computation/jc2_danielewski_two_by_three_weight_nonentry_thm3569.py
output: 05-knowledge/results/jc2_danielewski_two_by_three_weight_nonentry_thm3569.out
script_sha256: 2f8a7de8e6b57dcd0555807e16560838a2e678b92a842531526e211685415914
output_sha256: d6afee0c5aa39239a2b4f698787becb7866a8ff6f29c88a8ffabf9d37112ae95
hash_basis: LF-normalized bytes
---

# THM-3569 -- universal squarefree Danielewski two-by-three weight Darboux nonentry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The original
quadratic specialization is the first member of the universal squarefree
statement; the proof mechanism depends only on simple arms and their degree.

THM-3561 turns a rational triple collision into an everywhere-etale
polynomial map from `A2` to a smooth Danielewski surface.  A polynomial
Darboux pair on that surface would compose with the map to give a genuine
planar Jacobian counterexample.  THM-3572 produces analogous squarefree
targets, including a cubic target with a degree-five etale noninjective map.
The present proof is independent of those construction formulas: it closes
the next fewnomial cell, two weights by three weights, uniformly for every
smooth squarefree target of degree at least two.

All rings are over `C`.  Fix a squarefree polynomial `Sigma in C[b]` of
degree `N>=2`, and put

```text
A_Sigma=C[b,c,e]/(c^2 e-Sigma(b)),
Y_Sigma=Spec A_Sigma.                                  (1)
```

The distinct roots of `Sigma` are called the arms.  Give `A_Sigma` the
symplectic Poisson bracket

```text
{b,c}=c^2,              {c,e}=-Sigma'(b),
{b,e}=-2ce.                                            (1a)
```

## 1. Graded Poisson problem and statement

Use the grading

```text
wt(b)=0,                    wt(c)=1,                    wt(e)=-2. (2)
```

The weight-`r` piece of `A_Sigma` has the unique form `c^r f(b)` on `c!=0`,
where

```text
Sigma^ceil(-r/2) divides f         when r<0.             (3)
```

For two homogeneous pieces,

```text
{c^r f,c^s g}=c^(r+s+1) W_(r,s)(f,g),
W_(r,s)(f,g)=s f'g-r f g'.                              (4)
```

Suppose `P,Q in A_Sigma` satisfy `{P,Q}=1`.  Subtract scalar weight-zero
pieces.  The claim is:

> It is impossible for one of `P,Q` to have at most two nonzero homogeneous
> pieces and the other to have at most three.

Together with the smaller-cell proof below this leaves only

```text
(#supp P,#supp Q) in
{(2,>=4),(>=4,2),(>=3,>=3)}.                            (5)
```

Thus every putative Darboux pair needs at least six nonconstant homogeneous
pieces.

We record the universal smaller-cell input, so the present theorem has no
provisional dependency.  A homogeneous scalar bracket would have weights
`r,-r-1`.  After swapping, take `r>=0` and write

```text
P=c^r f,       Q=c^(-r-1) Sigma^m g,
m=ceil((r+1)/2).                                      (5a)
```

For `r=0` the bracket is divisible by `Sigma`; for `r>=2` it is divisible
by `Sigma^(m-1)`.  For `r=1` its degree is
`deg(f)+deg(g)+deg(Sigma)-1>=1`, with nonzero leading multiplier in
characteristic zero.  Hence no homogeneous scalar bracket exists.
If one whole output had a single weight, only the unique complementary
weight of the other output could contribute to the scalar row; every other
weight would produce a disjoint nonzero row and could be deleted after its
bracket vanished.  This would leave the forbidden homogeneous scalar pair.
Thus the one-by-arbitrary-width boundary is closed as well.

The complete two-by-two cell is also empty.  Its only possible cross-matched
normal form is

```text
P=c^(-R)f+c^(T-1)F,
Q=c^(-T)g+c^(R-1)G,                                  (5b)
```

where `R,T>=1`, `Sigma^ceil(R/2)|f`, and
`Sigma^ceil(T/2)|g`.  Its two extreme rows and scalar row are

```text
Rfg'-Tf'g=0,              (R-1)F'G-(T-1)FG'=0,
(R-1)f'G+RfG'-TF'g-(T-1)Fg'=1.                       (5c)
```

At an arm the scalar row can survive only through a simple negative
coefficient, forcing `R=2` or `T=2`; weight `-1` has zero derivative
multiplier and does not survive.  If `R=2`, the first equation and simplicity
give `T=2k`, `f=A h`, `g=B h^k`, and `Sigma|h`.  The second equation gives
`F=L K^(2k-1),G=M K`, and the scalar row factors as

```text
(Kh'+2hK')
(AM-k(2k-1)LB K^(2k-2)h^(k-1)).                      (5d)
```

The first factor has degree `deg K+deg h-1>=1`, with nonzero leading
coefficient.  The case `T=2` is symmetric.  This proves the universal
homogeneous and two-by-two inputs used below.

## 2. Support collision trichotomy

It remains to exclude the exact `2 x 3` cell.  Write the weights of `P` as

```text
r_0<r_1,                     delta=r_1-r_0>0.            (6)
```

A weight-zero bracket needs `r_i+s_j=-1`.  A single contributing pair
would itself be a forbidden homogeneous Darboux pair, so at least two pairs
must contribute.  Since each weight has a unique complement, exactly two do,
and they use both weights of `P`.  Their two `Q` weights are

```text
alpha=-r_1-1,                beta=-r_0-1=alpha+delta.    (7)
```

Let `x` be the third `Q` weight.  The two cross rows of the complementary
rectangle have weights `-delta,delta`; the rows involving `x` have weights

```text
x-beta,                      x-alpha.                    (8)
```

Unless

```text
x=alpha-delta             or             x=beta+delta,  (9)
```

all four nonconstant rows are isolated and must vanish separately.  Deleting
the inert third piece then gives the universally excluded two-by-two cell.
Thus only the lower and upper arithmetic extensions in `(9)` can survive.

We will use one arm-order lemma repeatedly.  At a root of `Sigma`, let
`m=ord(f), n=ord(g)`.  If `W_(r,s)(f,g)=0`, its first possible arm coefficient
is

```text
s m-r n.                                                 (10)
```

By `(3)`, a negative-weight coefficient vanishes on every arm.  Equation
`(10)` shows that a vanishing extreme row cannot pair a negative weight with
a positive weight.  A zero weight on the other side would force its retained
coefficient to be scalar, already removed.  Hence the lower extreme has two
negative weights and the upper extreme has two nonnegative weights.

After reparametrizing by positive integers `R,T`, the only two support
patterns are therefore

```text
LOWER:
 supp(P)={-R,T-1},
 supp(Q)={-(R+2T-1),-T,R-1};                            (11)

UPPER:
 supp(P)={-R,T-1},
 supp(Q)={-T,R-1,2R+T-2}.                              (12)
```

The scalar row in both cases is the sum of the complementary pairs
`(-R,R-1)` and `(T-1,-T)`.

## 3. The lower extension fails

At an arm, if `m=ord(f)` and `n=ord(G)`, the first complementary scalar
summand has initial multiplier

```text
(R-1)m+Rn.                                            (12a)
```

It can be nonzero of order zero only when `m=1,n=0,R>=2`; divisibility
`(3)` then forces `R=2`.  At `R=1` the derivative multiplier `R-1`
vanishes, so a unit `G` does not contribute.  The other complementary
summand gives the symmetric conclusion.  Thus

```text
R=2                         or                         T=2. (13)
```

If `R=2`, the lowest extreme equation has weights `-2,-(2T+1)`.  Since
these integers are coprime, its UFD solution is

```text
f=a h^2,                    g_0=B h^(2T+1).              (14)
```

Thus `f` is not simple on an arm, so the other scalar summand must survive
and `T=2`.

Conversely, suppose `T=2` and `R>2`.  The same extreme equation gives

```text
f=a h^(R/d),                g_0=B h^((R+3)/d),
d=gcd(R,3).                                                  (15)
```

At any arm put `k_0=ord(h)` and `m=R k_0/d`.  Since `R>2`, the scalar row
must survive through the weight-`-2` middle coefficient, so that coefficient
is simple there.  Its contribution to the shared nonzero row has exact order
`m` and coefficient `R-2m`, which is nonzero for `d in {1,3}` and integral
`k_0>=1`.  The other term has order at least

```text
m+3k_0/d-1>m.                                             (16)
```

For `d=1` the strict inequality is immediate.  For `d=3`, the coefficient
`g_0` has weight `-(R+3)`, so `(3)` gives

```text
(R+3)k_0/3 >= ceil((R+3)/2),
```

and hence `k_0>=2`; this again makes `(16)` strict.  The row cannot cancel.
Hence the sole lower candidate is

```text
R=T=2,
supp(P)={-2,1},                 supp(Q)={-5,-2,1}.        (17)
```

Write the extreme-row solutions as

```text
f=a h^2,       g_0=B h^5,       F=L K,       H=M K.      (18)
```

The shared nonzero row integrates exactly to

```text
g_1=h^2(D+lambda hK),             lambda=5LB/(2a).       (19)
```

Substitution in the scalar row factors it as

```text
h(hK)' [2aM-L(2D+3lambda hK)]=1.                         (20)
```

Because `f` has negative weight, `(3)` and squarefreeness of `Sigma` give
`Sigma|h`.  The left side of `(20)` is divisible by a nonconstant polynomial,
a contradiction.

## 4. The upper extension fails

In the upper support, the lowest extreme equation has the common-power
solution

```text
f=a h^(R/d_0),          g=B h^(T/d_0),
d_0=gcd(R,T).                                           (20a)
```

The scalar arm test says either the `f`-summand is simple, which forces
`R=2,d_0=2` and therefore even `T`, or the `g`-summand is simple, which
forces `T=2,d_0=2` and therefore even `R`.  Thus the exhaustive upper
possibilities are exactly two ladders:

```text
A: R=2, T=2k,
   supp(P)={-2,2k-1},       supp(Q)={-2k,1,2k+2};        (21)

B: T=2, R=2k,
   supp(P)={-2k,1},         supp(Q)={-2,2k-1,4k}.        (22)
```

Here `k>=1`.

### 4.1 Ladder A

Put

```text
p=2k-1,       q=2k+2,       d=gcd(p,q) in {1,3}.        (23)
```

The two extreme equations have the common-power form

```text
f=a h,        g=B h^k,
F=L K^(p/d),  H=M K^(q/d),             Sigma|h.         (24)
```

Let `G` be the coefficient of the middle weight `1`.  The shared positive
row is a first-order equation for `G`.

If `d=1`, it gives

```text
G=D K+C hK^3,                                               (25)
```

where `C` is a fixed nonzero scalar.  The scalar row then has the factor

```text
h'K+2hK'.                                                 (26)
```

If `d=3`, it gives

```text
G=C hK+G_0,                    K'G_0-3KG_0'=0.            (27)
```

For `G_0=0`, the scalar row has the factor `3h'K+2hK'`.  If `G_0!=0`, UFD
comparison in `(27)` gives `K=lambda J^3,G_0=mu J`, and the scalar row has
the factor

```text
h'J+2hJ'.                                                 (28)
```

Every factor in `(26)--(28)` has positive degree: for example the leading
coefficient of `(26)` is

```text
(deg(h)+2deg(K)) lc(h)lc(K),                              (29)
```

nonzero in characteristic zero, and `deg(h)>=2` because `Sigma|h`.  None can
divide the unit `1`.  Ladder A is empty.

### 4.2 Ladder B

The extreme equations now give

```text
f=a h^k,        g=B h,        F=L K,        H=M K^(4k),
Sigma|h.                                                   (30)
```

The shared positive row integrates to

```text
G=D K^(2k-1)+C h^k K^(4k-1),                              (31)
```

with fixed nonzero `C`.  Direct substitution factors the scalar row by

```text
h'K+2hK'.                                                 (32)
```

Its positive degree again contradicts a scalar value of one.  Ladder B is
empty, completing the two-by-three exclusion.

## 5. Boundary and counterexample design

This theorem strengthens the exact THM-3561 residual system and applies
unchanged to every squarefree target in THM-3572.  In the universal scope it
says:

```text
PROVED:   homogeneous, 2 x 2, 2 x 3, and 3 x 2 polynomial weight
          cells contain no Darboux pair for every squarefree Sigma
          of degree at least two;

OPEN:     the first six-piece cells 3 x 3 and 2 x 4;
          all larger mixed cells;
          existence of any polynomial Darboux pair on Y_Sigma;
          JC(2).                                          (33)
```

The rational Darboux pair `{b,-1/c}=1` is the sharp hostile:
it passes the bracket equation by retaining a pole and therefore lies
outside `A_Sigma`.

The degree threshold is sharp.  If `Sigma=ub+v` has degree one, then

```text
A_Sigma=C[c,e],                 {c,e}=-u,
{c,-e/u}=1.                                             (33a)
```

Thus the degree-one target is `A2` and already has a one-by-one homogeneous
Darboux pair.  If `Sigma` is not squarefree, a multiple root `beta` gives a
line `c=0,b=beta` on which `c^2`, `Sigma'(b)`, and `2ce` all vanish.  Every
Poisson bracket vanishes there, so **no** polynomial Darboux pair exists at
all.  That stronger no-go is singular and Poisson-degenerate; it is not a
smooth near-counterexample and is not the mechanism proved above.  Positive
characteristic remains outside scope because it can kill the derivative and
leading-degree multipliers.

For counterexample search, coefficient count alone is now too loose.  The
cheapest exact cells are the convolution polygons `3 x 3` and `2 x 4`, and
their nonzero output rows must be retained together: deleting an apparent
extra piece silently recreates the forbidden two-by-two rectangle.

## 6. Exact companion

Reproduce with

```bash
python3 04-computation/jc2_danielewski_two_by_three_weight_nonentry_thm3569.py
python3 -O 04-computation/jc2_danielewski_two_by_three_weight_nonentry_thm3569.py
```

The companion exhausts finite support-collision controls; checks squarefree
samples in degrees two through seven; verifies the sharp degree-one and
nonsquarefree boundaries, `(19)--(20)`, and every upper factorization for
`1<=k<=8`, including all `d=3` and cube-homogeneous branches.  These are exact
identity controls.  The all-degree conclusion is the symbolic support,
valuation, UFD, and first-order-equation argument above, not extrapolation
from the sampled degrees or eight ladder rows.

**QED.**
