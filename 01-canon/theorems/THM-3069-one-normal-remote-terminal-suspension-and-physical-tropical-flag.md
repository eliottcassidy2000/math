---
id: THM-3069
title: "One-normal remote-terminal suspension and physical tropical flag"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Adjoining one
  sufficiently remote terminal factorial slot to
  any fixed lower support suspends its first-window resultant: the new
  physical resultant is the kth power of the lower resultant times one
  explicit Beta carrier, up to an exponentially small error.  For four
  slots the lower sidecar is always nonzero by THM-2824, so every fixed low
  triple becomes first-window SFC(4) after its fourth slot is sufficiently
  remote.  Iteration gives actual physical lexicographically remote support
  towers in every slot number.  This is not a uniform moving-low theorem.
source: kind-pasteur-2026-08-01-remote-terminal-suspension
audit: >
  Two independent hostile audits rederived the determinant-one coordinate
  split, the exact Gamma layers and unique remote survivors, both sharp
  exponential gaps, the separated Artin norm and sign, one-variable
  covariance, the four-slot consequence, the recursive exponent ledger, the
  two noncommuting six-slot ordered cones, and the oriented finite-bank Beta
  holotopy.  Both replayed normal, optimized, and stored companions and
  matched the declared LF hashes; the immutable-file audit also checked every
  proved dependency and the fixed-low/simultaneous-scale boundary.
depends_on:
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
  - THM-3063-terminal-suspension-transverse-resultant-and-five-slot-tail-holotopy
related:
  - THM-2843-four-slot-projective-divisibility-and-resolvent-reduction
  - THM-2918-fixed-low-triple-high-gap-cube-root-trichotomy
  - THM-3047-formal-corner-width-product-gamma-moment-and-strict-hankel-positivity
  - THM-3054-affine-moving-lower-tropical-beta-gamma-tail-holotopy
script: 04-computation/gmc_one_normal_remote_terminal_suspension_thm3069.py
output: 05-knowledge/results/gmc_one_normal_remote_terminal_suspension_thm3069.out
script_sha256: bfdfe1621a5db35cb38a38233d74339bfe879b363a9536f9a089eaf972d683bc
output_sha256: 49590ed336c9af8dbc6b73b3c6f17c85cc97fc4fdb48127773c0d8d7d03cca64
hash_basis: LF-normalized bytes
---

# THM-3069 -- one-normal remote-terminal suspension

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3054 identifies a formal terminal face, but does not prove that the
actual physical first-window resultant has that face.  THM-2918 constructs
the three local cubic-null sheets after a remote atom is added, but explicitly
does not exclude escaping or boundary branches.  THM-3063 instead moves a
terminal pair with fixed gap and has two normal directions.

The complementary physical degeneration has only one normal direction.  In
that coordinate, one isotropic scaling contracts every physical coefficient
except the lower first-window system and one pure terminal monomial.  The
resultant then separates exactly.  This proves a global projective closure,
not only continuation of selected local sheets.

## 1. Physical family and lower sidecar

Let

```text
L(s^j)=j!,                         f_j=s^j/j!,          (1)
```

and fix

```text
k>=3, n>=0, 0<=a_0<...<a_(k-2)<C.                    (2)
```

In homogeneous coordinates `(x_0,...,x_(k-2))` on `P^(k-2)`, put

```text
Y_C=sum_(i=0)^(k-2) x_i(f_(n+a_i)-f_(n+C)),          (3)

F_(r,C)=L(Y_C^r)/L(f_n^r),                 2<=r<=k,
R_C=Res(F_(2,C),...,F_(k,C)).                         (4)
```

Thus `R_C` is the intrinsic normalized physical first-window resultant.
No formal corner or selected determinant is substituted for it.

Set

```text
w=sum_(i=0)^(k-2)x_i,
y_i=x_i,                                  0<=i<=k-3,
x_(k-2)=w-sum_(i=0)^(k-3)y_i.                         (5)
```

The change from `x` to `(y,w)` has determinant one.  If `c=a_(k-2)`, then

```text
Y_C=E(y)+w B_C,

E(y)=sum_(i=0)^(k-3)y_i(f_(n+a_i)-f_(n+c)),
B_C=f_(n+c)-f_(n+C).                                  (6)
```

Define the fixed lower forms and their resultant by

```text
H_r=L(E^r)/L(f_n^r),                       2<=r<=k,
S=Res_y(H_2,...,H_(k-1)).                             (7)
```

Here `S` is the first-window resultant on the fixed `k-1` lower slots, in
the determinant-one coordinates from `(5)`.  The general theorem assumes

```text
S!=0.                                                  (8)
```

No sign assumption on `S` is needed for nonvanishing.

## 2. Exact remote filtration and its two gaps

Put

```text
U_(k,C)=L(f_(n+C)^k)/L(f_n^k)
       =(kn+1)_(kC)/(n+1)_C^k>0,

lambda_C=U_(k,C)^(-1/k).                              (9)
```

Split a monomial contribution to `F_(r,C)` by

```text
j   = number of factors drawn from w B_C,
ell = number of those j factors equal to -f_(n+C).   (10)
```

The remaining factors have fixed degrees.  Stirling's Gamma-ratio expansion,
uniform over the finite coefficient bank, gives

```text
||lambda_C^j F_(r,j,ell,C)||
 <=poly(C)[ell^ell/k^j]^C,                            (11)
```

where `0^0=1`.  This estimate is coefficientwise in the fixed `(y,w)`
coordinates.  There are two model classes.  First,

```text
j=ell=0                                               (12a)
```

is not an error: it is exactly the fixed form `H_r`.  Second, the unique
remote survivor is

```text
r=j=ell=k.                                            (12)
```

Its contribution is exactly

```text
lambda_C^k L((-w f_(n+C))^k)/L(f_n^k)=(-1)^k w^k.   (13)
```

Every cell outside `(12a)` and `(12)` has a uniform gap.  Since `ell<=j`,

```text
ell^ell/k^j <=(ell/k)^ell.                            (14)
```

The function `x log(x/k)` is convex.  Its maximum on the integers
`1,...,k-1` is at an endpoint, and Bernoulli's inequality gives

```text
1/k<=((k-1)/k)^(k-1).                                (15)
```

Consequently the **isotropic coefficient gap** is

```text
rhohat_k=((k-1)/k)^(k-1)<1.                          (16)
```

It bounds every nonleading coefficient base in `(11)` and is sharp for this
physical whole-system filtration, attained exactly at

```text
j=ell=k-1.                                            (17)
```

The resultant has an additional character constraint.  Under `w->t w`,
every coefficient of `w`-degree `j` acquires `t^j`, while covariance makes
the resultant acquire exactly `t^(k!)`.  Thus every resultant monomial has
total `w`-weight `k!`.  Expand that monomial further by the number `ell` of
remote factors in each coefficient.  Its exponential base is a product of
terms `ell^ell`, with `ell<=j`, against the leading denominator `k^(k!)`.

The unique maximum **face** uses `(k-1)!` copies of the pure degree-`k`
terminal coefficient and the complete lower resultant `S^k`; the latter may
itself contain many coefficient monomials.  For every monomial off that face,
relax the actual Newton support and redistribute one or more parts of size
`k` into parts at most `k-1`.  Convexity of `x log x` makes the largest loss
the first split

```text
k -> (k-1)+1.
```

Therefore the sharper **scalar resultant-Newton gap** is

```text
rho_k=(k-1)^(k-1)/k^k=rhohat_k/k<rhohat_k.           (18)
```

If `L>=1` pure top cells are removed, the same convex redistribution gives
ratio at most `(rho_k L)^L<=rho_k` for `1<=L<=k-1`; increasing `L` by
`k-1` adds the strictly smaller factor
`((k-1)/k)^[k(k-1)]`.  This is the partition-sharp relaxed bound used in
THM-3054's formal face, now applied to the full expansion of the actual
physical resultant.  The actual Newton support may only improve it.

After substituting `w->lambda_C w` in every form, `(6)--(17)` give the
whole-system contraction

```text
(F~_(2,C),...,F~_(k,C))
 =(H_2,...,H_(k-1),H_k+(-1)^k w^k)
   +O(poly(C)rhohat_k^C)                              (19)
```

coefficientwise on one fixed projective space.  There is no root selection,
anisotropic line replacement, or termwise comparison of resultant monomials.

## 3. Exact separated resultant

For any scalar `a`, the separated system satisfies

```text
Res_(y,w)(H_2,...,H_(k-1),H_k+a w^k)
 =S^k a^[(k-1)!].                                    (20)
```

Here is a scheme-safe proof.  Assumption `(8)` makes

```text
A=C[y]/(H_2,...,H_(k-1))
```

an Artin complete intersection of length

```text
D=2*3*...*(k-1)=(k-1)!.                              (21)
```

Its positive-degree ideal is nilpotent, so multiplication by `H_k` is
nilpotent.  On the affine chart `w=1`, the relevant norm is therefore

```text
det_A(a I+m_(H_k))=a^D.                              (22)
```

The boundary `w=0` contributes the lower resultant `S` through the normal
degree `k`, hence `S^k`.  Equivalently, nonvanishing makes the left side of
`(20)` independent of every coefficient of `H_k`, after which the standard
separated product formula applies.  Calibration on

```text
(y_0^2,y_1^3,...,y_(k-3)^(k-1),w^k)                 (23)
```

fixes the scalar and sign in `(20)` to `+1`.  This retains all local
intersection multiplicities.

There are `k-1` variables and the form degrees are `2,...,k`.  Common
linear covariance of the homogeneous resultant thus gives exactly

```text
R~_C=lambda_C^(2*3*...*k)R_C
    =U_(k,C)^[-(k-1)!]R_C.                            (24)
```

Since `(k-1)!` is even for every `k>=3`, the scalar `a=(-1)^k` from `(13)`
contributes no sign.  Polynomial continuity of the fixed-degree resultant,
`(19)--(24)`, and `S!=0` first prove the leading coefficient and eventual
nonvanishing with the weaker gap `rhohat_k`.  The total-`w` character and
partition bound preceding `(18)` sharpen the scalar error to

```text
R_C=S^k U_(k,C)^[(k-1)!]
    (1+O(poly(C)rho_k^C)).                            (25)
```

In particular `R_C!=0` for all sufficiently large `C`, and

```text
sign(R_C)=sign(S^k)                                   (26)
```

eventually in the displayed standard orientation.

## 4. Exact base, power, and full symbol

Gauss multiplication gives

```text
U_(k,C)=k^(kC) product_(j=1)^(k-1)
          [(n+j/k)_C/(n+1)_C],                       (27)

U_(k,C)~kappa_(k,n) k^(kC) C^[-(k-1)/2],

kappa_(k,n)=Gamma(n+1)^(k-1)/
             product_(j=1)^(k-1)Gamma(n+j/k)>0.      (28)
```

Thus `(25)` becomes

```text
R_C~S^k kappa_(k,n)^[(k-1)!]
     k^(k! C) C^[-(k-1)!(k-1)/2].                    (29)
```

All coefficients in the finite resultant are Gamma-ratio expressions.  The
leading pure terminal coefficient in `(13)` is exact, so `(25)` has a full
exponentially separated Gamma-ratio symbol.  This is stronger than a bare
eventual sign assertion.

## 5. Unconditional remote fourth slot

Take `k=4` and fix any three lower offsets

```text
0<=a<b<c<C.                                           (30)
```

Now

```text
S_3=Res_y(H_2,H_3).                                   (31)
```

THM-2824 proves `S_3>0` for every arbitrary three-slot factorial support in
the standard orientation.  Therefore

```text
R_C=S_3^4 U_(4,C)^6
    (1+O(poly(C)(27/256)^C)),                         (32)

R_C~S_3^4 kappa_(4,n)^6 4^(24C) C^(-9)>0.           (33)
```

The whole transformed coefficient system contracts only at `27/64`; the
resultant character is what improves `(32)` to `27/256`.

Hence every fixed low triple, with arbitrary gaps, becomes first-window
SFC(4) after adjoining one sufficiently remote fourth slot.  Unlike
THM-2918, this excludes escaping and boundary projective branches as well as
the three local cube-root sheets.

There is an immediate THM-3063 consequence.  Its six-slot transverse sidecar
is an arbitrary four-slot first-window resultant.  That sidecar is therefore
automatically nonzero whenever its fourth fixed low slot is sufficiently
remote from its fixed lower triple.  The later moving terminal pair still has
the fixed-gap and fixed-low quantifiers of THM-3063.

## 6. Recursive physical tropical flags

Write

```text
R_j(a_0,...,a_(j-1))
```

for the normalized physical first-window resultant on `j` slots.  Formula
`(25)` says, with the first `j-1` offsets fixed and the last offset `C`
remote,

```text
R_j=R_(j-1)^j U_(j,C)^[(j-1)!]
    (1+O(poly(C)rho_j^C)).                            (34)
```

Starting from the positive two-slot quadratic resultant, choose each next
offset only after the preceding finite support has been fixed.  Induction
then constructs arbitrarily long supports for which every first-window
prefix resultant is positive.  In nested, lexicographically ordered limits,
the exponent ledger telescopes to

```text
R_k~R_2^[k!/2] product_(j=3)^k U_(j,a_(j-1))^[k!/j]. (35)
```

The powers `k!/j` are exactly the resultant coefficient multidegrees.  Thus
`(35)` is an actual physical tropical flag, not merely the formal corner of
THM-3054.

There are already two distinct nested six-slot cones.  Fix a low triple with
sidecar `S_3>0` and a gap `h>=1`; write `U_(r,T)` for `(9)` with
`k=r,C=T`.

1. First send a five-slot close pair `C,C+h` remote by THM-3063, freeze it,
   and then adjoin a one-normal terminal `R`.  In that ordered limit,

   ```text
   R_6~S_3^120 (5^h-4^h)^720
        U_(4,C)^180 U_(5,C)^144 U_(6,R)^120.          (35a)
   ```

2. First adjoin a remote fourth low slot `D` by this theorem, freeze it, and
   then send the six-slot close pair `C,C+h` remote by THM-3063.  The other
   ordered limit is

   ```text
   R_6~S_3^120 U_(4,D)^180 (6^h-5^h)^720
        U_(5,C)^144 U_(6,C)^120.                      (35b)
   ```

The powers in both formulas are forced by composition of resultant
multidegrees.  The formulas are not two descriptions of one simultaneous
cone: the scale order, the normal dimension at the first contraction, and
the surviving gap determinant differ.

The order of limits is load-bearing.  Formula `(34)` fixes every lower
support before sending the new terminal offset to infinity.  It supplies no
uniform threshold when several lower offsets move, no simultaneous ratio
cone, and no conclusion for an arbitrary support of the same width.

## 7. Beta carrier and finite-order holotopy

For every `S!=0`, orient once by

```text
epsilon=sign(S^k),
B_C=|S|^k U_(k,C)^[(k-1)!]>0.                        (36)
```

Each ratio in `(27)` is the moment sequence of

```text
Beta(n+j/k,1-j/k),                    1<=j<=k-1.      (37)
```

Consequently `(B_C)` is, up to a positive geometric scaling, a strict
full-support Hausdorff and Stieltjes product.  An affine subsequence
`C_c=sc+b` remains such a moment sequence by pushforward.

Put

```text
A=(k-1)!(k-1)/2,                  N=r(r-1)/2.         (38)
```

After the standard positive row and column factors are removed, a fixed
generalized carrier minor has first confluent term

```text
[product_(m=0)^(r-1) (A)_m/m!]
 V(u)V(v)c^(-2N)>0.                                  (39)
```

The relative error in `(25)` is exponentially smaller than `(39)`.  Hence,
for every fixed finite bank of generalized Hankel minors and

```text
T_c^(theta)=(1-theta)B_(C_c)+theta epsilon R_(C_c),
0<=theta<=1,                                          (40)
```

all minors in the bank are positive for every sufficiently large `c`,
uniformly in `theta`.  This is a literal carrier-to-physical finite-order
holotopy.  It is not one all-order Stieltjes tail.

## 8. Exact evidence and scope

The exact companion checks:

1. the covariance exponents, even sign, sharp isotropic coefficient gap, and
   sharper resultant-character partition gap for `3<=k<=12`;
2. the exact coordinate/layer partition on twelve physical forms;
3. four arbitrary three-slot sidecars and both signs of the exact separated
   identity `Res(H_2,H_3,H_4 +/- w^4)=S_3^4`;
4. a standard-orientation monomial control;
5. twelve intrinsic physical resultants from four unrelated low triples and
   three terminal gaps, all positive after division by the proved universal
   THM-2942 Macaulay flag;
6. exact Gauss/Beta identities and twelve strict carrier Hankel minors.

The largest-gap normalized physical ratios in the four rows are exactly
represented and print as

```text
0.999402940060125,
0.998769112178072,
0.994451207023886,
0.982389816402540.                                    (41)
```

Run

```text
python 04-computation/gmc_one_normal_remote_terminal_suspension_thm3069.py
python -O 04-computation/gmc_one_normal_remote_terminal_suspension_thm3069.py
```

Both modes must equal the stored transcript after LF normalization.

This theorem fixes every lower offset before moving one terminal slot.  It
does not prove arbitrary four-slot SFC, a uniform remote threshold, a
simultaneous moving-low limit, a shifted moment window, one all-order
Stieltjes tail, a wall-stripped Macaulay core, arbitrary-radial GMC(2), NC2,
LRC, or a Jacobian statement.

**QED.**
