---
id: THM-3192
title: "Reciprocal coefficient-jet transfer and z-adic Pluecker return"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Reciprocal reversal turns the inhomogeneous quadratic factorial-moment
  recurrence into an exact polynomial three-state transfer.  Its z-adic
  Smith types are (1,z,z) and, on exterior squares, (z,z,z^2); the hidden
  wedge returns in the second layer.  Every fixed top-jet quotient is
  homogeneous once its order is at most the current moment index.  On the
  offset-six PRS, the maintained arithmetic walls H,J,K are exactly two
  Toeplitz Pluecker charts and one wedge chart of these reciprocal jets, up
  to explicit p-adic units for p>=197.
source: root/multiscale-newton-flag/2026-08-02
audit: >
  The exact companion checks the reciprocal gauge, determinant, compound
  matrix, Smith witnesses, residue kernel/conormal, five first-layer flag
  words and their single rational cancellation gates, and the second-layer
  return.
  It compares 45 reversed coefficients and 35 homogeneous jet truncations
  with direct multinomial moments, then verifies the three symbolic H/J/K
  chart identities, five symbolic p-unit tests, and 54 direct integer
  A/B/R/S top coefficients at p=5,7,11.
  Exact path-cancellation and algebraic chart hostiles retain complementary
  coordinates.  Normal and optimized replay agree with the stored output.
  An independent immutable audit rederived the reciprocal gauge, Smith and
  exterior layers, arbitrary first-layer words, truncation boundary, all
  three H/J/K unit-chart identities, the U_K denominator, and both hostiles;
  it matched the LF hashes and accepted the scope.
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3176-six-step-prime-resonance-third-euclidean-newton-separation
  - THM-3186-full-exterior-continuant-path-convolution-and-cancellation-wall
related:
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
  - THM-3183-factorial-hecke-lattice-square-and-oriented-wedge-continuant
script: 04-computation/factorial_reciprocal_jet_pluecker_return_thm3192.py
output: 05-knowledge/results/factorial_reciprocal_jet_pluecker_return_thm3192.out
script_sha256: 2bc5197fb7589823b259a92ca3e2da17b5c342cf1e60752619c00740ced4bb30
output_sha256: f53735e7b6df34b71a5e427a0df0d8f1aa8295f8211f9b5454dad830a038ecea
hash_basis: LF-normalized bytes
---

# THM-3192 -- reciprocal coefficient-jet transfer and z-adic Pluecker return

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3186 gives the complete oriented path convolution for the bare scalar
exterior transfer, but it deliberately stops before the coefficient-degree
PRS.  The missing map is reciprocal reversal.  It turns top coefficients
into ordinary low-order jets, makes the boundary source disappear in every
fixed sufficiently shallow jet quotient, and exposes the offset-six factors
`H,J,K` as literal subresultant/Pluecker coordinates.

## 1. Exact reciprocal gauge

Retain THM-3124's resonant normalization

```text
q(t)=d-t+v t^2,                 M_n(v)=L(q^n),             (1)
```

and put

```text
F_n(z)=z^n M_n(1/z),            E_n(z)=(dz)^n.             (2)
```

The scalar state `Y_n=(M_n,M_(n-1),d^n)^t` of THM-3183 is changed by

```text
D_n(z)=diag(z^n,z^(n-1),z^n),
Ytilde_n=D_n(z)Y_n(1/z)=(F_n,F_(n-1),E_n)^t.              (3)
```

For every `n>=1`, direct conjugation gives

```text
Ytilde_(n+1)=J_n(z)Ytilde_n,

J_n(z)=D_(n+1)S_n(1/z)D_n^(-1)
      =[A_n  B_n                         C_n]
       [1    0                           0  ]
       [0    0                           dz ],              (4)

A_n=2(n+1)(2n+1),
B_n=n(n+1)z(z-4d),
C_n=(d-n-1)z.                                          (5)
```

Thus this is not a fitted recurrence or analogy with the scalar transfer: it
is that transfer in an exact moving reciprocal gauge.

Equivalently, the first coordinate of `(4)` is

```text
F_(n+1)=A_nF_n+B_nF_(n-1)+(d-n-1)zE_n.                   (6)
```

Since `E_n=d^n z^n`, the last term is precisely the reversed boundary source
`(d-n-1)d^n z^(n+1)`.

## 2. z-adic Smith layers and oriented return

Work over the DVR `Q(d,n)[[z]]`, with `d n(n+1)` inverted.  The determinant
and one valuation-one minor are

```text
det J_n=-dn(n+1)z^2(z-4d),
det J_n[{0,2},{0,2}]=A_n d z.                              (7)
```

Because `A_n` is a unit and `(7)` has orders one and two, respectively,

```text
Smith_z(J_n)=(1,z,z).                                      (8)
```

In the ordered wedge basis `(01,02,12)`,

```text
Lambda^2J_n=
[-B_n  -C_n   0      ]
[ 0     A_ndz B_ndz  ]
[ 0     dz    0      ].                                   (9)
```

Hence

```text
Smith_z(Lambda^2J_n)=(z,z,z^2).                            (10)
```

Every entry of `(9)` is divisible by `z`.  Its first normalized residue is

```text
(z^(-1)Lambda^2J_n)|_(z=0)=
[4dn(n+1)  n+1-d  0]
[0          A_nd  0]
[0          d     0].                                     (11)
```

This matrix has rank two, right kernel `e_12`, and left conormal

```text
(0,1,-A_n).                                                (12)
```

The whole first layer has a closed flag transport.  Put

```text
pi(x_0,x_1,x_2)=(x_0,x_1),
iota_n(u,v)=(u,v,v/A_n),                  pi iota_n=I,

G_n=[4dn(n+1)  n+1-d]
    [0          A_nd ].                                  (12a)
```

If `R_n` denotes the matrix in `(11)`, then

```text
R_n=iota_n G_n pi,
R_m...R_n=iota_m(G_m...G_n)pi                 (m>=n).    (12b)
```

Thus every first-layer word has common kernel `e_12`; its image chart is the
plane `x_2=x_1/A_m` determined only by the final index.  The internal
holonomy is upper triangular, and its diagonal ratio is

```text
[G_n]_(00)/[G_n]_(11)=2n/(2n+1),                         (12c)
```

independent of `d`.  All parameter-sensitive first-layer transport is
therefore confined to the off-diagonal exit entry `n+1-d`.

More explicitly, for a consecutive word from `n` through `m`, write

```text
G_m...G_n=[A_(m,n) B_(m,n);0 C_(m,n)].                    (12d)
```

Then

```text
B_(m,n)=d^(m-n)(U_(m,n)-d V_(m,n)),                       (12e)

U_(m,n)=sum_(j=n)^m
  prod_(i=j+1)^m 4i(i+1) (j+1)
  prod_(i=n)^(j-1) 2(i+1)(2i+1),

V_(m,n)=sum_(j=n)^m
  prod_(i=j+1)^m 4i(i+1)
  prod_(i=n)^(j-1) 2(i+1)(2i+1).                          (12f)
```

For positive integral `n<=m`, both `U_(m,n)` and `V_(m,n)` are positive
integers.  Apart from the common factor `d^(m-n)`, every fixed first-layer
word therefore has exactly one rational off-diagonal cancellation gate
`d=U_(m,n)/V_(m,n)`.  It is a positive weighted average of
`n+1,...,m+1`, hence lies in that interval.  This single-gate law belongs to
the linear reciprocal carrier; the multiple `H,J,K` walls arise only after
the nonlinear PRS chart projections in Section 4.

The missing column is not dead.  It returns exactly one layer later:

```text
(z^(-2)Lambda^2J_n e_12)|_(z=0)
   =-4d^2n(n+1)e_02.                                      (13)
```

Equations `(11)--(13)` are the reciprocal coefficient-jet version of the
hidden-to-visible return in THM-3186.  They also name the complementary chart
which survives when the first selected wedge coordinate cancels.

## 3. Fixed top jets are genuinely homogeneous

For `h>=0`, let

```text
Jet_h=Q(d)[z]/(z^(h+1)).                                  (14)
```

The coefficient identity

```text
[z^j]F_n=[v^(n-j)]M_n                    (0<=j<=n)         (15)
```

shows that `Jet_h(F_n)` is exactly the top `h+1` coefficient jet of `M_n`,
written from highest degree downward.

If `h<=n`, the boundary term in `(6)` has degree `n+1>h` and vanishes in
`Jet_h`.  Therefore, for every consecutive run beginning at
`N>=max(1,h)`,

```text
(F_(n+1),F_n)^t
 ==[A_n B_n;1 0](F_n,F_(n-1))^t             in Jet_h,
                                                        n>=N.   (16)
```

The qualification `h<=N` is load-bearing.  The initial jets still contain
all boundary history accumulated before `N`; `(16)` says only that their
subsequent propagation is homogeneous.  It does not erase the source at
deeper jet order.

## 4. Literal projection to the offset-six PRS charts

Now use THM-3176's offset-six family

```text
p>=197 prime,             d=p+6,
A=A_(p+4),                B=A_(p+5),
(A,B)->(A,R)->(R,S)->(S,T).                               (17)
```

Normalize every coefficient below by `(2p)!` and write

```text
a_j=[v^(p+j)]A/(2p)!,       b_j=[v^(p+j)]B/(2p)!,
r_j=[v^(p+j)]R/(2p)!,       s_j=[v^(p+j)]S/(2p)!.         (18)
```

Thus the reciprocal top jets used below are

```text
a=(a_4,a_3,a_2),      b=(b_5,b_4,b_3),
r=(r_3,r_2,r_1),      s=(s_2,s_1).                        (19)
```

For triples `f=(f_0,f_1,f_2)` and `g=(g_0,g_1,g_2)`, define the selected
Toeplitz chart

```text
P_2(f,g)=det[f_0 0   g_0;
             f_1 f_0 g_1;
             f_2 f_1 g_2].                               (20)
```

For pairs, put

```text
P_1(f,g)=f_0g_1-f_1g_0.                                  (21)
```

These are the actual pseudo-division coordinates.  Indeed,

```text
f_0^2g-[f_0g_0+P_1(f,g)z]f=P_2(f,g)z^2+O(z^3).           (22)
```

Thus `P_2` is the first surviving coefficient after cancelling two leading
terms, while `P_1` is the connection coefficient in the quotient.

Let

```text
H=24p+109,
J=256p^4-27648p^3-365600p^2-1528800p-2096649,
K=5120p^5-810240p^4-14447872p^3-92004672p^2
    -256323456p-265142241.                                (23)
```

Direct exact simplification of the reciprocal jets gives

```text
P_2(a,b)=U_H r_3,
P_2(r,a)=U_J s_2,
P_1(s,r)=U_K K,                                           (24)
```

where

```text
U_H=256 prod_(i=1)^4(p+i)^2
          (2p+1)^2(2p+3)^2(2p+5)^2(2p+7),

U_J=64 prod_(i=1)^4(p+i)^2
         (p+5)(2p+1)^2(2p+3),

U_K=64 prod_(i=1)^5(p+i)^2
         (p+6)(2p+1)(2p+7)/(2p-1).                       (25)
```

The inherited leading rows are

```text
r_3=-8 prod_(i=1)^5(p+i)(2p+1)(2p+3)H,
s_2= 4 prod_(i=1)^5(p+i)(2p+7)J.                         (26)
```

Every displayed factor in `(25)--(26)` other than `H,J,K` is a `p`-adic
unit for `p>=197`.  Therefore the three maintained arithmetic chart ideals
are exactly

```text
(P_2(a,b))=(H),          (P_2(r,a))=(J),
(P_1(s,r))=(K)                                      over Z_(p). (27)
```

This is the promised source-to-PRS map.  `H` and `J` are not mysterious
extra scalars attached after the recurrence, and `K` is not merely a
connection name: they are successive reciprocal-jet Pluecker coordinates.

### Proof of `(24)--(27)`

Formula `(15)` gives every entry in `(18)` directly from the exact
multinomial coefficient sum.  Substitute those entries into THM-3176's
definitions of `R,S`.  Expanding the determinants `(20),(21)` and factoring
over `Q(p)` gives `(24),(25)`.  The same substitution gives `(26)` and the
two top cancellations

```text
D_3r_3-P_1s_2=0,
D_3r_2-P_1s_1-P_0s_2=0,
P_0=4(p+6)(2p+1)K,                                       (28)
```

which independently confirms the `K` connection coordinate.  Finally the
inequalities `p>=197` make every other factor a unit in `Z_(p)`, proving
`(27)`.

## 5. Two sharp chart boundaries

First, take the actual reciprocal transfer with

```text
d=5,                 z=105/4,                 n=1,2,3.    (29)
```

Acting on `e_12`, the three-step exterior product is

```text
(60z^4(z-20)(4z-105),
 3000z^4(z-20)(z^2-20z+140),
 7500z^4(z-20))^t.                                       (30)
```

At `(29)` its `e_01` coordinate is zero and both complementary coordinates
are nonzero.  Every local `J_1,J_2,J_3` is invertible there.  This is exactly
THM-3186's two-path cancellation after reciprocal gauge: the selected path
amplitude dies, not the exterior state.

Second, as a sharp algebraic chart control only, specialize the universal
offset-six formulas to the root `p=-109/24` of `H`.  Then

```text
P_2(a,b)=0,
P_1((a_4,a_3),(b_5,b_4))
 =7207252894016784826007415625/4437222213480873984 !=0.    (31)
```

This value of `p` is not a prime and is not used in the offset-six closure.
It proves the formal boundary is a selected Toeplitz-chart wall rather than
rank death of the reciprocal pair.

## 6. Consequence and remaining obstruction

The previously missing connection in THM-3183/THM-3186 is now explicit:

```text
factorial recurrence
  -> reciprocal top-jet transfer
  -> Toeplitz/wedge Pluecker projection
  -> exact H,J,K PRS charts.                              (32)
```

The transfer preserves more than any one selected chart, and `(30),(31)`
show why an atlas is necessary.  What remains for arbitrary fixed offset is
not source typing: it is to iterate `(20)--(22)` through an unbounded PRS,
control simultaneous chart walls, and prove that some chart separates by the
empirical depth `floor(s/2)`.  This theorem does not prove that law, arbitrary
support `SFC`, `NC(2)`, or `GMC(2)`.

## 7. Exact evidence

Run

```text
python 04-computation/factorial_reciprocal_jet_pluecker_return_thm3192.py
python -O 04-computation/factorial_reciprocal_jet_pluecker_return_thm3192.py
```

and compare LF-normalized bytes with the declared output.  The companion
uses exact integer, rational, and symbolic arithmetic only.  It pins all six
dependency blobs; verifies `(4)--(13)`; checks 45 direct reciprocal
coefficients, five first-layer flag words and their convolution gates, and
35 truncations; reconstructs
`(24)--(28)`; proves five
displayed factors are `p`-adic units throughout `p>=197`; and compares 54
symbolic `A/B/R/S` top coefficients with independent direct multinomial
values at `p=5,7,11`.  It also verifies every coordinate and local determinant in
`(30)` and the complementary chart in `(31)`.  There is no floating point,
random sampling, imported executable, or assertion-sensitive test.

**QED.**
