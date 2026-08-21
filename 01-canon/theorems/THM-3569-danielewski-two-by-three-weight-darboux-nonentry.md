---
id: THM-3569
title: "Danielewski two-by-three weight Darboux nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the smooth
  Danielewski target c^2 e=b(b+4) of THM-3561, a polynomial Darboux pair
  cannot have at most two nonconstant weight pieces in one output and at
  most three in the other.  After scalar removal every possible pair must
  have support sizes (2,at least 4), (at least 4,2), or (at least 3,at
  least 3), hence at least six nonconstant homogeneous pieces.  The result
  closes the complete two-by-three cell only; it neither supplies nor
  excludes a general Darboux pair and proves no JC(2) conclusion.
source: codex-2026-08-20-frontier-overnight
audit: >
  An independent hostile audit rederived the support-collision trichotomy,
  arm-order gates, both arithmetic extensions, every exceptional
  factorization, and the exact companion.  The ordinary and optimized
  replays agree.  The squarefree-arm, characteristic-zero, polynomial, and
  support-size boundaries are explicit.
depends_on:
  - THM-3561-rational-keller-danielewski-polynomial-completion
related:
  - THM-3404-factorized-danielewski-principal-parts-and-finite-cover-obstruction
script: 04-computation/jc2_danielewski_two_by_three_weight_nonentry_thm3569.py
output: 05-knowledge/results/jc2_danielewski_two_by_three_weight_nonentry_thm3569.out
script_sha256: 98e062a15921b9cb0e5b17afdd5bf206d11970631c3891230c1b86b0d064dfa3
output_sha256: 9f364247a3ae9114db6d8d143391facd9adb98593fb5e22ded34cefdfc430922
hash_basis: LF-normalized bytes
---

# THM-3569 -- Danielewski two-by-three weight Darboux nonentry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3561 turns a rational triple collision into an everywhere-etale
polynomial map from `A2` to a smooth Danielewski surface.  A polynomial
Darboux pair on that surface would compose with the map to give a genuine
planar Jacobian counterexample.  THM-3561 excludes homogeneous pairs and the
complete two-by-two weight cell.  The present theorem closes the next
fewnomial cell, two weights by three weights, in every degree.

All rings are over a characteristic-zero field.  Put

```text
S=b(b+4),
Y=Spec k[b,c,e]/(c^2 e-S).                              (1)
```

The two roots of `S` are called the arms.

## 1. Graded Poisson problem and statement

Use the THM-3561 grading

```text
wt(b)=0,                    wt(c)=1,                    wt(e)=-2. (2)
```

The weight-`r` piece of `k[Y]` has the unique form `c^r f(b)` on `c!=0`,
where

```text
S^ceil(-r/2) divides f             when r<0.             (3)
```

For two homogeneous pieces,

```text
{c^r f,c^s g}=c^(r+s+1) W_(r,s)(f,g),
W_(r,s)(f,g)=s f'g-r f g'.                              (4)
```

Suppose `P,Q in k[Y]` satisfy `{P,Q}=1`.  Subtract scalar weight-zero
pieces.  The claim is:

> It is impossible for one of `P,Q` to have at most two nonzero homogeneous
> pieces and the other to have at most three.

Together with THM-3561 this leaves only

```text
(#supp P,#supp Q) in
{(2,>=4),(>=4,2),(>=3,>=3)}.                            (5)
```

Thus every putative Darboux pair needs at least six nonconstant homogeneous
pieces.

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
the inert third piece then gives the two-by-two cell excluded by THM-3561.
Thus only the lower and upper arithmetic extensions in `(9)` can survive.

We will use one arm-order lemma repeatedly.  At a root of `S`, let
`m=ord(f), n=ord(g)`.  If `W_(r,s)(f,g)=0`, its first possible arm coefficient
is

```text
s m-r n.                                                 (10)
```

By `(3)`, a negative-weight coefficient vanishes on both arms.  Equation
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

At an arm, the scalar row can be nonzero only if the coefficient of weight
`-R` has a simple zero or the coefficient of weight `-T` has a simple zero.
Divisibility `(3)` and the nonzero derivative multipliers in `(4)` force

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

At an arm put `k_0=ord(h)` and `m=R k_0/d`.  Simplicity of the weight-`-2`
middle coefficient makes the first term of the shared nonzero row have
exact order `m` and nonzero coefficient `R-2m`.  The other term has order

```text
m+3k_0/d-1>m.                                             (16)
```

For `d=3`, divisibility `(3)` forces `k_0>=2`, so the strict inequality
still holds.  The row cannot cancel.  Hence the sole lower candidate is

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

Because `f` has negative weight, `(3)` and squarefreeness of `S` give
`S|h`.  The left side of `(20)` is divisible by a nonconstant polynomial,
a contradiction.

## 4. The upper extension fails

The same arm test, now combined with the lowest common-power equation,
leaves exactly two ladders:

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
F=L K^(p/d),  H=M K^(q/d),             S|h.             (24)
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

nonzero in characteristic zero, and `deg(h)>=2` because `S|h`.  None can
divide the unit `1`.  Ladder A is empty.

### 4.2 Ladder B

The extreme equations now give

```text
f=a h^k,        g=B h,        F=L K,        H=M K^(4k),
S|h.                                                       (30)
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

This theorem strengthens only the exact THM-3561 residual system.  It says:

```text
PROVED:   homogeneous, 2 x 2, 2 x 3, and 3 x 2 polynomial weight cells
          contain no Darboux pair;

OPEN:     the first six-piece cells 3 x 3 and 2 x 4;
          all larger mixed cells;
          existence of any polynomial Darboux pair on Y;
          JC(2).                                          (33)
```

The rational Darboux pair `{b,-1/c}=1` from THM-3561 is the sharp hostile:
it passes the bracket equation by retaining a pole and therefore lies
outside `k[Y]`.  Repeated-arm Danielewski equations can change `(3)` and the
arm-order lemma; positive characteristic can kill the derivative
coefficients.  Neither boundary is covered.

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

The companion exhausts finite support-collision controls, verifies `(19)--
(20)`, and checks every upper factorization for `1<=k<=8`, including all
`d=3` and cube-homogeneous branches.  These are exact identity controls.  The
all-degree conclusion is the symbolic support, valuation, UFD, and
first-order-equation argument above, not extrapolation from eight rows.

**QED.**
