---
id: THM-2873
title: "Two-ray factorial-response TP3 curvature"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  For every nonzero nonnegative two-ray cone
  W=x d_p+y d_(p+h), with h=1 or h=2, the consecutive three-column
  factorial-response determinant is strictly positive at every row n>=1.
  Each of its seven polarized coefficients has a common positive factorial
  factor and a residual polynomial whose n=z+1 expansion has strictly
  positive coefficients.  The cutoff is sharp: at n=0 the singleton d_4
  has negative curvature.  This gives local TP3 curvature for both sparse
  cones in the shared-high-tooth cubic-null family.  On the canonical
  THM-2846 cubic-null rectangle, a separate exact interval certificate makes
  the endpoint holonomy strictly negative throughout and excludes that one
  cell from the shared quartic line.  The general mixed secant comparison
  and shared multipole-line branch remain open.
source: codex/two-ray-factorial-response-tp3-2026-07-28
depends_on:
  - THM-2830-disjoint-positive-adjacent-cone-factorial-moment-three-detection
  - THM-2846-arbitrary-positive-cone-moment-three-transverse-boundary
  - THM-2866-positive-factorial-difference-semiring-and-cubic-pascal-response-ladder
  - THM-2872-four-slot-shared-multipole-quartic-norm-and-response-secant-reduction
related:
  - THM-2841-all-order-adjacent-difference-factorial-tensor-positivity
script: 04-computation/gmc_two_ray_response_tp3_thm2873.py
output: 05-knowledge/results/gmc_two_ray_response_tp3_thm2873.out
script_sha256: 03e613db2c0bf51f90c1fd09370aa96f691f4ecbb61492eaaa3f23b84a15ff31
output_sha256: 521c47f95d3c6b4e4aece53e0027ab842fe33f792017b40adc868b461a49813a
hash_basis: LF-normalized bytes
---

# THM-2873 -- two-ray factorial-response TP3 curvature

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

Put

```text
L(s^k)=k!,                   f_j=s^j/j!,
d_j=f_(j+1)-f_j.                                      (1)
```

For a nonzero nonnegative adjacent-difference cone `W`, define

```text
A_n(W)=L(d_nW),        B_n(W)=L(d_nW^2),
C_n(W)=L(d_nW^3),                                    (2)

Delta_n(W)=det [
  A_n       B_n       C_n
  A_(n+1)   B_(n+1)   C_(n+1)
  A_(n+2)   B_(n+2)   C_(n+2)
].                                                    (3)
```

The theorem proves that for

```text
W=x d_p+y d_(p+h),       x,y>=0,       (x,y)!=(0,0),
p>=0,                    h in {1,2},                   (4)
```

one has

```text
Delta_n(W)>0                         for every n>=1.   (5)
```

Since `A_n(W)>0`, divide row `r` in `(3)` by `A_(n+r)` and put

```text
R_n=B_n/A_n,                         S_n=C_n/A_n.      (6)
```

Then `(5)` is exactly the strict positive-curvature identity

```text
(R_(n+1)-R_n)(S_(n+2)-S_(n+1))
 -(S_(n+1)-S_n)(R_(n+2)-R_(n+1))>0.                 (7)
```

THM-2830 and THM-2866 separately make `R_n` and `S_n` strictly increasing.
Thus the consecutive secant slopes of the two-ray response curve turn
strictly counterclockwise.

The cutoff in `(5)` is sharp.  Already

```text
Delta_0(d_4)<0.                                       (8)
```

The mechanism is not raw total positivity of the positive-semiring output
measures.  Those measures already have a negative `2 x 2` minor for `W=d_1`.
The sign is created only after the exact Pascal prefix/tooth weighting.

## 1. Closed factorial tensor

For nonnegative indices `a_1,...,a_k`, set

```text
b_i=a_i+1,                    B=sum_i b_i.             (9)
```

Since

```text
d_a=s^a(s-a-1)/(a+1)!,
```

direct expansion gives

```text
L(prod_i d_(a_i))
 =1/(prod_i b_i!) *
  sum_(ell=0)^k (-1)^ell e_ell(b_1,...,b_k)(B-ell)! . (10)
```

This is also the marked-cycle/forbidden-predecessor inclusion--exclusion
formula behind THM-2841, but only the elementary identity `(10)` is needed
here.

Write

```text
(u)^[r]=u(u+1)...(u+r-1),             (u)^[0]=1.      (11)
```

Fix `m in {1,2,3}`, `0<=e<=m`, and a row shift `r in {0,1,2}`.  Let

```text
M_(m,e,r)
 =L(d_(n+r) d_p^(m-e) d_(p+h)^e).                    (12)
```

Use the block-size list

```text
b_(m,e,r)=(
 n+r+1,
 p+1 repeated m-e times,
 p+h+1 repeated e times
),                                                    (13)
```

and put

```text
T=n+r+mp+eh,

E_(m,e,r)^h
 =(n+mp+1)^[r+eh]
   sum_(ell=0)^(m+1)
    (-1)^ell e_ell(b_(m,e,r))
    (T+1)^[m+1-ell].                                  (14)
```

Factoring `(n+mp)!` from `(10)` yields the load-bearing normalization

```text
M_(m,e,r)
 =(n+mp)! E_(m,e,r)^h
  /[(n+r+1)!(p+1)!^(m-e)(p+h+1)!^e].                 (15)
```

Every `E_(m,e,r)^h` is an ordinary polynomial in `n,p`; there are no
symbolic factorials left.

## 2. The seven residual polynomials

The binomial expansion of `(4)` gives

```text
L(d_(n+r)W^m)
 =sum_(e=0)^m binom(m,e)x^(m-e)y^e M_(m,e,r).         (16)
```

For `0<=q<=6`, define the residual polynomial

```text
C_(h,q)(n,p)
 =[X^(6-q)Y^q]
  det_(r=0,1,2; m=1,2,3) [
    sum_(e=0)^m binom(m,e)X^(m-e)Y^e E_(m,e,r)^h
  ].                                                  (17)
```

In every determinant term contributing to the coefficient in `(17)`,
the three column powers are `1,2,3`, every row shift occurs once, and
exactly `q` of the six labelled atoms use the high ray.  Therefore all
factorials removed in `(15)` are common to that coefficient.  Equations
`(15)--(17)` give the exact identity

```text
[x^(6-q)y^q] Delta_n(W)

 =(n+p)!(n+2p)!(n+3p)! C_(h,q)(n,p)
  /[
    (n+1)!(n+2)!(n+3)!
    (p+1)!^(6-q)(p+h+1)!^q
   ].                                                 (18)
```

Every factor outside `C_(h,q)` is strictly positive.

Now set

```text
n=z+1,                         z>=0.                  (19)
```

Exact sparse expansion of all fourteen residuals gives:

```text
h=1:
  number of nonzero monomials for q=0,...,6
  =69,83,98,114,131,149,168;

h=2:
  number of nonzero monomials for q=0,...,6
  =69,98,131,168,209,254,303.                         (20)
```

Every coefficient in every polynomial

```text
C_(h,q)(z+1,p)                                        (21)
```

is strictly positive.  The minimum coefficients, in order
`q=0,...,6`, are

```text
12,72,180,240,180,72,12                               (22)
```

for both gaps.  Thus `(18)` makes all seven coefficients of
`Delta_n(xd_p+yd_(p+h))` positive when `n>=1`.  At least one endpoint
coefficient survives when `(x,y)!=(0,0)`, proving `(5)`.

The complete coefficient dictionaries are generated, checked, and
canonically hashed by the exact companion.  Their hashes are

```text
h=1:
c97819308c4ea3817b15fdc63197ad58b4fe9e5a141cbe0f5bb5c989e072ce60

h=2:
a234806f39364c4788826d6001c2cd82f8cbd17e0b033d9ca1c9cbdf2c9ad926. (23)
```

This coefficient expansion is the proof certificate, not a bounded
evaluation of `(5)`.

## 3. Sharp singleton boundary

For a singleton `W=d_j`, write the oriented response curvature as

```text
Gamma_n(d_j)
 =Delta_n(d_j)/
  [A_n(d_j)A_(n+1)(d_j)A_(n+2)(d_j)].                (24)
```

Direct simplification of `(10)` gives

```text
Gamma_n(d_j)
 =12 P(j,n)(2j+n)!(3j+n)!
  /[
    (n+1)(n+2)(n+3)(j!)^3
    (j+n+1)!(j+n+2)!
   ],                                                 (25)
```

where

```text
P(j,n)=
  j^5(8n-4)
 +j^4(28n^2+50n-5)
 +j^3(38n^3+150n^2+173n+53)
 +j^2(25n^4+152n^3+330n^2+302n+99)
 +j  (8n^5+65n^4+203n^3+305n^2+221n+62)
 +n^6+10n^5+40n^4+82n^3+91n^2+52n+12.              (26)
```

After `n=z+1`, `(26)` has `27` nonzero coefficients, all positive, with
minimum `1`.  At the omitted row `n=0`,

```text
P(j,0)=-4j^5-5j^4+53j^3+99j^2+62j+12.               (27)
```

Its values at `j=0,1,2,3` are

```text
12,217,748,1143,
```

but for `j=k+4`,

```text
P(k+4,0)
 =-(4k^5+85k^4+667k^3+2305k^2+3002k+140)<0.         (28)
```

In particular,

```text
Gamma_0(d_4)=-4527600.                                (29)
```

Thus `n>=1` in `(5)` is exact even on a single ray; an all-row TP3
statement is false.

## 4. Why raw output-measure TP is not the proof

The positive multiplication semiring of THM-2866 expands every response
column as a nonnegative Pascal mixture.  It is tempting to demand total
positivity of the three unweighted output-index measures before applying
the Pascal matrix.  That implication fails at the first nontrivial ray.

For `W=d_1`, the power-one output measure has unit mass at index `1`.
The pair law

```text
d_1^2=2f_2+6d_3
```

makes the power-two Pascal measure have mass `2` at each prefix index
`0,1` and mass `6` at index `3`.  The output rows `0,1` therefore have
minor

```text
0*2-1*2=-2.                                           (30)
```

So neither raw output TP2 nor raw output TP3 can prove `(5)`.  Formula
`(18)` shows what repairs the failure: the Pascal prefix weights cause
exact cancellations, and the surviving Newton coefficients are positive.
In the marked-cycle model, these survivors are the iterated
predecessor-shielding residues; a uniform injection for all six labels is
still open.

## 5. Shared-high-tooth consequence and canonical endpoint exit

The sparse cubic-null family used in the shared-multipole investigation is

```text
U=d_a+x d_(a+2),            V=d_(a+1)+y d_(a+2),
a>=1,                       x,y>0.                    (31)
```

Thus `U` is a gap-two cone and `V` a gap-one cone.  The theorem proves
strict local TP3 curvature for both response curves at every response row
`n>=1`.

This is a real input to the response-secant reduction of THM-2872, but it
does **not** compare the two mixed transport secants.  In the notation of
that reduction, the still-needed sign is

```text
beta kappa_U < alpha kappa_V,                          (32)
```

equivalently the strict sign of the quartic endpoint-holonomy determinant.
Separate convexity of the two response curves does not imply `(32)`, and
when `a=1` one transport endpoint also touches the sharp row `n=0`
boundary `(29)`.  No closure of the shared cubic--quartic multipole line is
claimed.

### 5.1 The canonical THM-2846 cell exits at the endpoint

The specific cubic-null cell of THM-2846 is

```text
U=d_1+x d_3,                    V=d_2+y d_3.           (33)
```

That theorem proves a unique common zero of the two cubic remainders in
the rational rectangle

```text
2636/10000 <= x <= 2637/10000,
23418/1000000 <= y <= 23421/1000000.                 (34)
```

Define, with no suppressed normalization,

```text
g_0=L(U^2),          g_1=L(UV),          g_2=L(V^2),

A_i=L(U^(4-i)V^i),                       0<=i<=4,      (35)
```

and the endpoint-holonomy determinant

```text
J=
 (2A_1g_0-A_0g_1)g_2^2
 -(2A_3g_2-A_4g_1)g_0^2.                              (36)
```

Thus

```text
J=g_0^2g_2^2(r_1^(L)-r_1^(R))                         (37)
```

for the two endpoint quotient determinations in THM-2872.

Exact expansion makes `J` a `30`-term polynomial in `x,y`.  Because both
coordinates in `(34)` are positive, each monomial is coordinatewise
increasing.  Bounding a positive coefficient at the lower/lower and
upper/upper corners, and reversing those choices for a negative
coefficient, gives the rigorous rational enclosure

```text
-1965120367409971587404977893331001634459
 ---------------------------------------------------
 3125000000000000000000000000000000

 <= J <=

-239202401274182466677656578205832473701
 --------------------------------------------------
 390625000000000000000000000000000

 <0.                                                   (38)
```

This is deliberately a crude separable monomial interval, not a decimal
root evaluation.  It holds uniformly throughout the whole rectangle
`(34)`, hence at THM-2846's exact algebraic cubic-null point.

At that point THM-2872 identifies

```text
beta kappa_U=(3/2)r_1^(L),
alpha kappa_V=(3/2)r_1^(R).                            (39)
```

Equations `(37)--(39)` prove the strict endpoint exit

```text
beta kappa_U < alpha kappa_V.                          (40)
```

Therefore the canonical THM-2846 cubic line is not a shared quartic line;
the midpoint defect never needs to be tested.  This strengthens the
mechanism attached to that one hostile from “fourth moment nonzero” to a
strict, rectangle-stable endpoint-holonomy mismatch.  It does not extend
`(40)` to every shifted cell in `(31)` or to an arbitrary four-slot
cone-cutting plane.

## 6. Universal TP3 lead: finite-exact only

For a general positive cone `W=sum mu_jd_j`, full polarization of `(3)`
has a symmetric six-label kernel: choose one labelled atom for column one,
two for column two, and the remaining three for column three, then sum the
`60` ordered set partitions.  Its universal positivity for `n>=1` remains
open.

The exact companion nevertheless records a stronger finite signal.  For
all

```text
0<=a_1<=...<=a_6<=6                                  (41)
```

there are `924` multisets.  The polarized kernel has degree at most
`a_1+...+a_6` in `n`.  Every one of all `17,556` Newton coefficients at
`n=1` is strictly positive; the minimum is `720`.  Equivalently, every
tested kernel has a strictly positive expansion in

```text
binom(n-1,r).                                         (42)
```

This is `FINITE-EXACT` evidence for the all-cone TP3 conjecture, not its
proof.

## 7. Exact companion and status

The exact companion:

1. derives `(10)` by integer inclusion--exclusion;
2. constructs all fourteen residuals `(17)` symbolically;
3. verifies every coefficient in `(21)` is positive;
4. checks `(18)` against direct factorial tensors in `168` exact cells;
5. derives and checks the singleton formula `(25)` in `54` exact cells;
6. certifies the complete `n=0` sign boundary `(27)--(29)`;
7. checks the raw-output hostile `(30)`;
8. constructs `(35)--(36)` directly from exact tensors and certifies the
   rational interval `(38)` term by term; and
9. performs the separately scoped `924`-multiset Newton scout.

Reproduce with

```text
python3 04-computation/gmc_two_ray_response_tp3_thm2873.py
python3 -O 04-computation/gmc_two_ray_response_tp3_thm2873.py
```

Both modes must byte-match

```text
05-knowledge/results/gmc_two_ray_response_tp3_thm2873.out.
```

The proof and exact companion are complete.  Promotion to proved canon
still requires an independent audit of the factorial normalization, all
determinant multiplicities, the positivity certificates, the interval
orientation in `(38)`, and the scope boundary in Section 5.
