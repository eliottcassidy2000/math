---
id: THM-2922
title: "Diameter-five Macaulay--Newton atlas"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  For every n>=0, first-window SFC(4) holds on all six
  translated four-slot supports of diameter exactly five.  The same
  fixed degree-seven Macaulay row chart used at diameter four has
  determinant degree at most 254 after exact width-five denominator
  clearing.  Six complete positive Gregory--Newton certificates,
  together with separately signed exceptional depths, make that minor
  nonzero everywhere.  The four-moment bound and the positive-depth
  five-monomial two-charge Gaussian bound eight are exact.
source: codex-gmc-holotopy-extension-2026-07-29
depends_on:
  - THM-2173-sparse-projective-factorial-moment-floor
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
related:
  - THM-2849-four-slot-first-window-macaulay-box
  - THM-2908-consecutive-four-slot-projective-resultant-closure
  - THM-2917-all-three-slot-diameter-four-factorial-detection
  - THM-2921-diameter-four-nonconsecutive-macaulay-newton-closure
script: 04-computation/gmc_diameter_five_macaulay_newton_atlas_thm2922.py
output: 05-knowledge/results/gmc_diameter_five_macaulay_newton_atlas_thm2922.out
script_sha256: 9f146073cf2f953bbbfcaaf20f1f58dbb25a105316d97a486e80db2be62b9b2e
output_sha256: 585240e0015005f946f27dcf4774ee50065811decd9e509b6774a5cad9b6164d
constructor_dependency_sha256: 42e9b5ceddd677d1f2601a9d5d668c9437281596b65999ddcb8549d4e0b9bf64
hash_basis: LF-normalized bytes
---

# THM-2922 -- diameter-five Macaulay--Newton atlas

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

## 1. Statement

Let

```text
L:C[s] -> C,                          L(s^j)=j!.         (1)
```

Fix an integer `n>=0` and one of the six offset sets

```text
(0,1,2,5), (0,1,3,5), (0,1,4,5),
(0,2,3,5), (0,2,4,5), (0,3,4,5).                     (2)
```

For every nonzero

```text
H=sum_(b in B) c_b s^(n+b),                            (3)
```

at least one of

```text
L(H),                 L(H^2),                 L(H^3), L(H^4) (4)
```

is nonzero.  These are exactly the six translated four-slot support
types of diameter five.

On every support in `(2)`, there is also a polynomial with all four
coefficients nonzero whose first three moments vanish and whose fourth
does not.  Thus the four-moment bound is exact with full support.

## 2. Mean elimination and width-five denominators

Normalize

```text
f_a=s^a/a!,                          L(f_a)=1.           (5)
```

and eliminate the mean through

```text
H=x(f_n-f_(n+5))
 +y(f_(n+a)-f_(n+5))
 +z(f_(n+b)-f_(n+5)),                                (6)
```

where `B=(0,a,b,5)`.  Divide the order-`m` moment by
`L(f_n^m)`.  As in THM-2921, its symmetric tensor entry is

```text
T_m(d_1,...,d_m;n)
 =(mn+1)_(sum_j d_j)/prod_j(n+1)_(d_j),               (7)
```

and the ordinary coefficient of a monomial with exponent triple
`alpha` is its signed top-replacement sum times the load-bearing
multinomial factor

```text
m!/(alpha_0!alpha_1!alpha_2!).                         (8)
```

Put the resulting forms in `R=Q[x,y,z]` and call them `Q,C,F`.
Exact polynomial division for every coefficient in all six families
proves that the common denominators

```text
D_2=prod_(j=1)^4(n+j),

D_3=[prod_(j=1)^4(n+j)^2](n+5),

D_4=[prod_(j=1)^4(n+j)^3](n+5)^2                     (9)
```

clear `Q,C,F`, respectively.  The coefficient degrees of

```text
Q~=D_2 Q,                 C~=D_3 C,                 F~=D_4 F
```

are at most

```text
4,                           9,                        14. (10)
```

Every factor in `(9)` is nonzero for integer `n>=0`, so the scaled and
unscaled forms have the same projective common zeros.

## 3. One fixed Macaulay chart

Use the degree-seven Macaulay map

```text
Phi_n:R_5 direct_sum R_4 direct_sum R_3 -> R_7,
Phi_n(A,B,D)=A Q~+B C~+D F~.                           (11)
```

In the THM-2849 monomial order its stored transpose has `46` rows and
`36` columns.  For every family in `(2)`, select the identical rows

```text
0,...,19;
21,...,29,35;
36,...,41.                                             (12)
```

The resulting maximal minor uses `20Q+10C+6F`.  If its determinant is
`P_B(n)`, then `(10)` gives

```text
deg P_B <=20*4+10*9+6*14=254.                         (13)
```

This is the same Pluecker chart that closes all three nonconsecutive
diameter-four families in THM-2921.  The new content is that it remains
a global chart across all six width-five support shapes.

## 4. Exact Gregory--Newton atlas

The companion evaluates the determinant exactly at the `255`
consecutive depths required by `(13)`.  The complete certificates are:

| offsets `B` | GN base | signs below base | scaled-form digest | GN-vector digest |
|---|---:|:---:|---|---|
| `(0,1,2,5)` | `3` | `---` | `f84fa9fe3a35770a0e1fc08bd3b8314117b7067035a9b11da6def17bbfb49c89` | `c9860215b2cbbe320a3fbd443725bf03c4d8e4ba56e2287dbd0be220557d9fb8` |
| `(0,1,3,5)` | `3` | `---` | `b92c83a1321adad0c53a470619d2daa152256dd73deb4d6f6d22d5c35c6021f3` | `ab09c820e9c8d93e358006fdd6f43d4d3fff4cb3116740235f2c45c6baed8924` |
| `(0,1,4,5)` | `3` | `---` | `420b9dc0ab9886781b42f3ab3ee71e710d26f15faa6447404ea719b0e575382b` | `6daeb2905271faa1acb8998300052dbcc92589e0a5517b1ad887dba58051d3e1` |
| `(0,2,3,5)` | `1` | `-` | `3933ea0b4b9c91fd49aa075ab59cc7365f68e3ab6eec48a711a7c5a08e4db5f7` | `47722854412a3e10572ee76bc95d4b2a9f53c352a946ab1161261ebeb7476806` |
| `(0,2,4,5)` | `1` | `-` | `4f00cb220d2bfadd3723bbc9f94c810c35428d51e196d9b9e59a49b870238077` | `f6f989f14e29c4ea49fb26be9b09c8c21a9da1f030712ae49c9f61fc9fcba3a2` |
| `(0,3,4,5)` | `0` | none | `b0b0e7c137979818c428cf6409ce35b4abbb6504159ed66c744952d575242d89` | `86f4f56722e75e41dfcd467eed5aaff46c9bb3fe1d907d120320a172fe50bdeb` |

Every one of the `255` Gregory--Newton coefficients in every row is
strictly positive.  Since `deg P_B<=254`,

```text
P_B(n)=sum_(k=0)^254 Delta^k P_B(r) binom(n-r,k),
                                                        n>=r, (14)
```

where `r` is the displayed base.  Hence every determinant is positive
from its base onward.  The first three determinants are separately
negative and nonzero at `n=0,1,2`; the next two are negative and
nonzero at `n=0`; the last is positive from base zero.  Therefore

```text
P_B(n)!=0                  for every B in (2), n>=0.   (15)
```

## 5. Projective and Gaussian consequences

The nonzero maximal minor gives `rank Phi_n=36`, so

```text
(Q~,C~,F~)_7=R_7.                                      (16)
```

A common projective zero would annihilate all of `R_7`, including a
seventh power of a nonzero coordinate.  This is impossible.  Restoring
the eliminated mean proves `(4)`.

THM-2173 supplies, in each four-dimensional support envelope, a nonzero
`H` whose first three moments vanish.  Every coordinate face has only
three slots, so PROVED THM-2824 excludes that witness from every face.
It therefore has full support, and `(4)` makes its fourth moment
nonzero.

For `n>=1`, write `H=s h(s)`, set `s=ZW` for a standard complex
Gaussian, and put

```text
P=alpha W+Z h(ZW),                       alpha!=0.      (17)
```

Charge balance gives

```text
E[P^(2j+1)]=0,
E[P^(2j)]=binom(2j,j)alpha^j L(H^j).                  (18)
```

Thus the corresponding genuinely two-charge five-monomial envelope has
exact uniform detection depth eight.  The full-support sharp witness
has moments one through seven zero and eighth moment

```text
70 alpha^4 L(H^4)!=0.                                  (19)
```

## 6. Conjectural continuation boundary

The exact width-four and width-five calculations expose a simple
candidate continuation.  For a general support

```text
(0,a,b,M),                       0<a<b<M,              (20)
```

the observed common denominator pattern is

```text
D_(m,M)
 =[prod_(j=1)^(M-1)(n+j)^(m-1)](n+M)^(m-2).          (21)
```

If `(21)` continues to clear every coefficient, the scaled order-`m`
row degree is at most

```text
(m-1)M-1,                                              (22)
```

and the fixed `20Q+10C+6F` minor has degree at most

```text
20(M-1)+10(2M-1)+6(3M-1)=58M-36.                     (23)
```

An independent exact continuation scout has now tested the complete
diameter-six atlas.  It finds

```text
D_(m,6)=[prod_(j=1)^5(n+j)^(m-1)](n+6)^(m-2),
row degrees=(5,11,17),                 minor degree<=312,
```

and the same row chart is nonzero on all ten diameter-six types.  The
four families `(0,1,b,6)` have positive `313`-term Newton vectors from
base `5` and five negative exceptional depths; the three families
`(0,2,b,6)` have them from base `2` and two negative exceptional depths;
the three families whose first positive offset is at least `3` are
Newton-positive from base zero.  Seventy direct four-variable modular
minor checks also pass.  This is **FINITE-EXACT continuation evidence,
not a theorem dependency or a diameter-six claim of THM-2922**.

Equations `(21)--(23)` remain **conjectural for general `M` (now first
untested at `M=7`)**.  Even a proof of denominator clearing would not
show that this one minor stays nonzero, much less that one shifted
Gregory--Newton vector is positive.  The next structural target is to
derive the observed base law and positivity uniformly in `M`, retaining
alternate minors if this Pluecker chart crosses a wall.

## 7. Scope and exact verification

The theorem is an infinite translated-depth result but only for:

```text
four slots; diameter exactly five; first four moments. (24)
```

It does not prove arbitrary SFC(4), a shifted moment window, SFC(5),
arbitrary-charge GMC(2), or the general formulas `(21)--(23)`.
Arbitrary three-slot detection is used only as the already proved
coordinate-face sidecar; it is not the novel payload.

The companion hash-pins the corrected THM-2921 constructor, rebuilds
`(9)` independently, proves all coefficient divisions inside `Z[n]`,
checks every exceptional low depth, and verifies all `1530` positive
Gregory--Newton coefficients.  A direct four-variable multinomial
constructor also reproduces the selected determinant modulo `1000003`
at seven hostile depths in each family, for `42` exact checks.  Original
and scaled depth-zero determinants agree by the exact row-scaling
invoice.

Run

```text
python 04-computation/gmc_diameter_five_macaulay_newton_atlas_thm2922.py
python -O 04-computation/gmc_diameter_five_macaulay_newton_atlas_thm2922.py
```

Normal and optimized executions byte-match the stored output with the
declared LF-normalized hashes.

**QED (candidate pending independent hostile audit).**
