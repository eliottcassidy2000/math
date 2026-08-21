---
id: THM-3579
title: "Equal-step three-by-three Danielewski Darboux nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  squarefree polynomial Sigma of degree at least two, no polynomial Darboux
  pair on c^2 e=Sigma(b) can have both nonconstant weight supports equal to
  length-three arithmetic progressions with the same step.  All five possible
  scalar-row alignments are excluded.  No general Darboux nonentry and no
  counterexample to JC(2) is claimed; unequal-step three-by-three cells and
  the two-by-four cell remain open.
source: kps-s188 / delegated arithmetic-progression support attack, 2026-08-21
audit: >
  ACCEPT.  The Laurent-chart compiler, five scalar-row split, arm-unit gate,
  extreme common-power descent, central cube first integral, reflected
  degree gates, terminal coprimality, and sharp boundary hostiles were
  independently checked.  Normal and optimized exact transcripts are
  byte-identical.
depends_on:
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
related:
  - THM-3561-rational-keller-danielewski-polynomial-completion
  - THM-3572-squarefree-danielewski-affine-modification-and-two-bracket-collapse
script: 04-computation/jc2_equal_step_three_by_three_danielewski_nonentry_thm3579.py
output: 05-knowledge/results/jc2_equal_step_three_by_three_danielewski_nonentry_thm3579.out
script_sha256: ad5ec528af81693a0f71fe703ac2ec332f50bb088bc1821786f08d136bb9cd5d
output_sha256: 71f4bdef61bb7d7f66a2fc7bb0d720abe1b8118b32f1d36235912f78bbd8025c
hash_basis: LF-normalized bytes
---

# THM-3579 -- equal-step three-by-three Danielewski Darboux nonentry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This closes the
additively minimal part of the first `3 x 3` support frontier left by
THM-3569 and THM-3572.  It does not close arbitrary three-by-three supports,
produce a Darboux pair, or settle the planar Jacobian conjecture.

All rings are over `C`.  Let `Sigma in C[b]` be squarefree of degree at least
two and put

```text
A_Sigma=C[b,c,e]/(c^2 e-Sigma(b)).                    (1)
```

Give this ring the symplectic Poisson bracket

```text
{b,c}=c^2,             {c,e}=-Sigma'(b),
{b,e}=-2ce                                                (2)
```

and the grading

```text
wt(b)=0,               wt(c)=1,               wt(e)=-2. (3)
```

On `c!=0`, a homogeneous piece of weight `u` is uniquely `c^u f(b)`.  Its
regularity in `(1)` is equivalent to

```text
Sigma^ceil(-u/2) divides f               when u<0.       (4)
```

After constants are removed, suppose that both outputs have exactly three
nonzero homogeneous pieces and that their supports are arithmetic
progressions with the same positive step `d`:

```text
P=c^r F(b,c^d),             Q=c^s G(b,c^d),             (5)

F=f_0+f_1 t+f_2 t^2,        G=g_0+g_1 t+g_2 t^2,
t=c^d.                                                       (6)
```

Every retained weight-zero coefficient is assumed nonconstant.  A constant
one can be subtracted without changing the bracket, reducing to a smaller
support cell already covered by THM-3569.

The theorem is

```text
{P,Q} != 1.                                               (7)
```

## 1. The Laurent compiler and five scalar rows

For homogeneous pieces define

```text
W_(u,v)(f,g)=v f'g-u f g'.                              (8)
```

Then

```text
{c^u f,c^v g}=c^(u+v+1) W_(u,v)(f,g).                  (9)
```

Equivalently, direct differentiation of `(5)` gives the Laurent-chart
compiler

```text
{P,Q}=c^(r+s+1) B_(r,s,d)(F,G),                        (10)

B_(r,s,d)(F,G)
 =s F_b G-r F G_b+d t(F_b G_t-F_t G_b).               (11)
```

The coefficients of `B` have convolution multiplicities

```text
1,2,3,2,1.                                             (12)
```

If `(7)` failed, there would be a unique scalar-row index
`kappa in {0,1,2,3,4}` such that

```text
r+s+1+kappa d=0,
B_(r,s,d)(F,G)=t^kappa.                                (13)
```

Thus the `kappa` row equals one and the other four rows vanish.  We exclude
all five possibilities.

## 2. Three local algebra lemmas

### 2.1 Extreme common powers

At a root of `Sigma`, a vanishing extreme Wronskian cannot pair a negative
weight with a positive weight.  Indeed, if the negative coefficient has
positive arm order and the other has nonnegative order, the first possible
Wronskian coefficient has a fixed nonzero sign.  A zero-weight boundary can
vanish only when its coefficient is constant, which has already been
removed.

Consequently, whenever both extreme rows vanish, the two low endpoint
weights are negative and the two high endpoint weights are nonnegative.
For example, if the low weights are `-R,-T`, then

```text
R f g'=T f'g.                                          (14)
```

Writing `delta=gcd(R,T)`, unique factorization gives

```text
f=A h^(R/delta),          g=B h^(T/delta)              (15)
```

for nonzero constants `A,B` and a nonconstant polynomial `h`.  Every arm is
a root of `h`.  The analogous statement at the high extreme will be used
without repetition.

### 2.2 The simple-arm scalar gate

A complementary scalar pair has weights `(-N,N-1)`, up to exchanging the
two outputs.  At an arm `beta`, let the negative coefficient have order `m`
and the nonnegative coefficient order `n`.  By `(4)`,

```text
m>=ceil(N/2),                    n>=0.                 (16)
```

The first possible scalar term has order and multiplier

```text
m+n-1,                 (N-1)m+Nn.                    (17)
```

It can be a nonzero unit exactly when

```text
N=2,                         m=1, n=0.                (18)
```

For `N=1,m=1,n=0`, the order is zero but the multiplier vanishes.  This is
the boundary that blocks the tempting weight `-1/0` channel.

### 2.3 Two degree-rigid operators

For nonzero polynomials `H,u` with `deg H+deg u>0`, put

```text
E_H(u)=H'u+2Hu',              J_H(u)=Hu'+2H'u.         (19)
```

Their leading terms are

```text
deg E_H(u)=deg J_H(u)=deg H+deg u-1,

lc E_H(u)=(deg H+2deg u)lc(H)lc(u),
lc J_H(u)=(deg u+2deg H)lc(H)lc(u).                   (20)
```

In characteristic zero these coefficients do not vanish.  In particular,
neither operator can equal one if `deg H>=2`; nor can `J_H(u)=1` when
`Sigma|u` and `deg Sigma>=2`.

## 3. The two one-channel rows

If `kappa=0` or `kappa=4`, the scalar row has one summand.  It would itself
be a homogeneous Darboux pair.  For completeness, the arm gate reduces such
a pair to weights `(-2,1)`.  If its negative coefficient is `H`, then
`Sigma|H`, and its bracket has positive degree with nonzero leading
coefficient.  It cannot equal one.  Thus

```text
kappa notin {0,4}.                                     (21)
```

This is also the homogeneous boundary recorded in THM-3569.

## 4. The lower off-centre row `kappa=1`

Write the low endpoint weights as `-R,-T`.  Equation `(13)` gives

```text
R+T=d+1,                                               (22)
```

and the two supports are

```text
supp(P)={-R,T-1,R+2T-2},
supp(Q)={-T,R-1,2R+T-2}.                              (23)
```

The scalar row consists of the two complementary pairs involving one low
endpoint and the opposite middle coefficient.  By `(18)`, one must have
`R=2` or `T=2`.  Suppose `R=2`.  The low extreme `(15)` says

```text
f=A h^(2/delta),       g=B h^(T/delta),
delta=gcd(2,T).                                        (24)
```

For the first scalar channel to survive, `f` must be simple at every arm.
Hence `delta=2`, `T=2k`, and `h` is simple at the arms.  Swapping the outputs
handles the other possibility.  We have reduced the entire row to

```text
supp(P)={-2,2k-1,4k},
supp(Q)={-2k,1,2k+2},                 k>=1,            (25)

f=A h,                    g=B h^k,      Sigma|h.       (26)
```

Let `a,q` be the two middle coefficients.  The scalar row is

```text
f'q+2fq'-2k a'g-(2k-1)ag'=1.                         (27)
```

Substituting `(26)` gives the exact factorization

```text
E_h(Aq-kB h^(k-1)a)=1.                                (28)
```

Since `deg h>=deg Sigma>=2`, equation `(28)` contradicts `(20)`.  Therefore
`kappa=1` is impossible.

## 5. The upper off-centre row `kappa=3`

Let the high endpoint weights be `U,V`.  Equations `(13)` and the extreme
sign lemma give

```text
U,V>=0,                    U+V=d-1.                   (29)
```

The scalar row uses the middle/high pairs.  Applying `(18)` shows, after
possibly swapping the outputs,

```text
(U,V)=(d-2,1).                                         (30)
```

For `d=1` this is impossible.  For `d=2`, one high weight is zero; the high
extreme equation is `F'G=0`, so its coefficient is constant and removable.
That leaves a two-by-three pair, excluded by THM-3569.

Put `n=d-2>=1`.  The remaining supports are exactly

```text
supp(P)={-n-4,-2,n},
supp(Q)={-2n-3,-n-1,1}.                               (31)
```

The high extreme Wronskian gives

```text
F=L K^n,                       G=M K.                 (32)
```

Let `a,q` be the middle coefficients of weights `-2,-n-1`.  The scalar row
becomes

```text
J_K(Ma-nL K^(n-1)q)=1.                                (33)
```

Both `a` and `q` have negative weight, so `(4)` gives

```text
Sigma divides Ma-nL K^(n-1)q.                         (34)
```

If the polynomial in `(34)` is zero, the left side of `(33)` is zero.  If it
is nonzero, its degree is at least two, and `(20)` again forbids a scalar.
Thus `kappa=3` is impossible.

## 6. The central triangle `kappa=2`

Put

```text
r=-R,                 s=-T,              R,T>=1.      (35)
```

Then

```text
R+T=2d+1                                                   (36)
```

is odd.  The three complementary scalar channels have weights

```text
(-R,R-1),
(p,-p-1),                    p=(T-R-1)/2,
(T-1,-T).                                             (37)
```

The endpoint channel `(-R,R-1)` could survive only for `R=2`.  But then
`T` is odd, so `gcd(R,T)=1`, and `(15)` makes its negative coefficient a
square in `h`, not simple at an arm.  The `T=2` endpoint is symmetric.  The
middle channel must therefore be the unique `(-2,1)` channel from `(18)`:

```text
p=1               or               p=-2,
|R-T|=3.                                               (38)
```

The genuine geometric symmetry here is only reflection: if `{P,Q}=1`, then
`{Q,-P}=1`.  Use it to take

```text
T=R+3,                    d=R+1.                       (39)
```

The supports and coefficients may now be written

```text
P=c^(-R) f+c a+c^(R+2)F,
Q=c^(-R-3)g+c^(-2)q+c^(R-1)G.                        (40)
```

The middle scalar channel is `(a,q)`.  At every arm, `(18)` makes `q`
simple and `a` a unit.

### 6.1 The lower bridge forces the `3m/3m+3` ladder

Let

```text
delta=gcd(R,R+3)=gcd(R,3),
alpha=R/delta,                 beta=(R+3)/delta.       (41)
```

The low extreme gives

```text
f=A h^alpha,                   g=B h^beta.             (42)
```

At an arm let `ord(h)=ell`.  The row immediately below the scalar row is

```text
R f q'-2f'q-(R+3)a'g-ag'=0.                           (43)
```

Because `q` is simple and `a` is a unit, the first two terms in `(43)` have
exact order `alpha ell`, with nonzero initial multiplier

```text
alpha(delta-2ell).                                    (44)
```

This cannot vanish because `delta in {1,3}` is odd.  The last two terms have
exact order `beta ell-1`.  Cancellation in `(43)` therefore requires

```text
(beta-alpha)ell=(3/delta)ell=1.                       (45)
```

Hence

```text
delta=3,                ell=1,               R=3m     (46)
```

for some `m>=1`, and `h` is simple at every arm.  In particular,

```text
f=A h^m,                         g=B h^(m+1).          (47)
```

### 6.2 The cube first integral fixes the middle coefficient

Substituting `(47)` into `(43)` and cancelling `h^(m-1)` gives

```text
mA(3h q'-2h'q)
 =(m+1)B h(3h a'+h'a).                                (48)
```

Put

```text
lambda=(m+1)B/(mA).                                   (49)
```

The identity

```text
3h(ha)'-2h'(ha)=h(3ha'+h'a)                           (50)
```

shows that `w=q-lambda ha` satisfies

```text
3h w'-2h'w=0.                                         (51)
```

In the fraction field,

```text
(w^3/h^2)'=(w^2/h^3)(3h w'-2h'w)=0.                  (52)
```

Thus `w^3/h^2` is constant.  If that constant were nonzero, a simple arm of
`h` would force

```text
3 ord(w)=2,                                            (53)
```

which is impossible.  Therefore

```text
q=lambda h a.                                         (54)
```

This is precisely where the squarefree-arm sidecar defeats the formal cube
resonance.

### 6.3 The upper bridge gives the terminal coprimality contradiction

Since `R=3m`, the high extreme has coprime weights `R+2,R-1` and gives

```text
F=L K^(R+2),                       G=M K^(R-1).        (55)
```

Substituting `(54)--(55)` in the row above the scalar row yields

```text
0=K^R [
 M(R-1)(a/K)'-lambda L(R+2)(haK^2)'
].                                                     (56)
```

Integration in `C(b)` and multiplication by `K` give

```text
a [M(R-1)-lambda L(R+2)hK^3]=C K                     (57)
```

for some constant `C`.

The bracketed polynomial is nonconstant.  It is coprime to `K`, because
modulo `K` it is the nonzero constant `M(R-1)`.  If `C=0`, `(57)` would make
that nonconstant polynomial vanish identically.  If `C!=0`, unique
factorization would make it divide a constant.  Both alternatives are
impossible.  This excludes `kappa=2` and completes the proof.

**QED.**

## 7. Sharp hostiles and symmetry ledger

The following boundaries are load-bearing.

1. **Cube resonance without a simple arm.**  The formal equation in `(51)`
   has the nonzero solution

   ```text
   h=b^3,                    w=b^2.                    (58)
   ```

   Thus the first integral alone is insufficient; the simple squarefree arm
   is what kills its nonzero constant.

2. **Degree one.**  For `h=b,u=1`, one has `E_h(u)=1`.  More geometrically,
   if `Sigma=ub+v`, then `A_Sigma=C[c,e]` and

   ```text
   {c,-e/u}=1.                                             (59)
   ```

   The hypothesis `deg Sigma>=2` is sharp.

3. **The rational pole.**  On the Laurent chart,

   ```text
   {b,-1/c}=1.                                             (60)
   ```

   This is the sharp polynomiality hostile: `-1/c` is not in `A_Sigma`.

4. **The selected symmetric near miss.**  The THM-3572 supports

   ```text
   {-3,-1,1},                  {-2,0,2}                 (61)
   ```

   are not a survivor.  At an arm with local parameter `z`, the extreme
   common powers make the central scalar row

   ```text
   2z(3AMz-BL),                                           (62)
   ```

   hence zero on the arm.

5. **Triangle symmetry is combinatorial, not geometric.**  The central row
   has three complementary channels, so it carries an abstract `S3`
   bookkeeping symmetry.  Arm divisibility and the linear order of weights
   do not preserve all those permutations.  Only the reflection
   `(P,Q)->(Q,-P)` is used in the proof.  No `S3` action on `Y_Sigma` is
   inferred.

## 8. Exact scope and remaining frontier

There is a useful additive-combinatorial corollary.  For three-element
integer sets `A,B`, the sharp sumset inequality and its equality case say

```text
|A+B|>=5,
|A+B|=5 iff A and B are arithmetic progressions with the same step. (62a)
```

Taking `A=supp(P)` and `B=supp(Q)`, this theorem removes the equality case.
Thus every still-live `3 x 3` Darboux cell has at least six potential output
weight rows.  This counts convolution addresses before cancellation; it does
not assert that all six bracket coefficients are nonzero.  In particular,
the surviving frontier has genuine additive defect rather than a hidden
reparametrization of the minimal arithmetic-progression cell.

The theorem proves

```text
PROVED:
  every equal-step 3 x 3 weight cell is empty for every squarefree
  Sigma of degree at least two, for every scalar-row alignment;

OPEN:
  unequal-step 3 x 3 cells;
  the 2 x 4 and 4 x 2 cells;
  all larger mixed cells;
  existence of any polynomial Darboux pair on Y_Sigma;
  JC(2).                                                  (63)
```

The proof is a Darboux-support nonentry theorem on a non-`A2` target.  It is
not a planar Jacobian counterexample and does not say that every three-by-
three support has equal step.

## 9. Exact verification contract

Reproduce with

```bash
python3 04-computation/jc2_equal_step_three_by_three_danielewski_nonentry_thm3579.py
python3 -O 04-computation/jc2_equal_step_three_by_three_danielewski_nonentry_thm3579.py
```

The optimization-safe exact companion checks:

- the Laurent compiler `(10)--(11)` and convolution word `1,2,3,2,1`;
- squarefree split controls in degrees two through seven;
- the complete simple-arm scalar gate `(16)--(18)`;
- the degree and leading coefficient formulas `(20)`;
- the entire row-one and row-three ladder identities for parameters through
  twelve;
- the central support and arm-order censuses, the exact bridge identities
  for `1<=m<=8`, and the cube first integral;
- terminal coprimality controls for `1<=m<=8`; and
- all four sharp hostiles above.

The finite ranges are exact identity and hostile controls.  The theorem is
the all-degree valuation, UFD, differential, and leading-degree proof in
Sections 2--6; it is not extrapolated from those ranges.
