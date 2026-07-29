---
id: THM-2924
title: "Diameter-six Macaulay--Newton atlas"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
  HOSTILE AUDIT.  For every n>=0, first-window SFC(4) holds on all ten
  translated four-slot supports of diameter exactly six.  The fixed
  degree-seven Macaulay chart from THM-2921/2922 has determinant degree
  at most 312 after exact width-six denominator clearing.  All 3,130
  shifted Gregory--Newton coefficients are positive, all exceptional
  low depths are strictly nonzero, and an independent direct
  four-variable constructor agrees at 70 hostile cells.  Full-support
  sharpness and the positive-depth Gaussian detection bound eight are
  exact.
source: codex-gmc-uniform-width-extension-2026-07-29
depends_on:
  - THM-2173-sparse-projective-factorial-moment-floor
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
related:
  - THM-2849-four-slot-first-window-macaulay-box
  - THM-2908-consecutive-four-slot-projective-resultant-closure
  - THM-2921-diameter-four-nonconsecutive-macaulay-newton-closure
  - THM-2922-diameter-five-macaulay-newton-atlas
script: 04-computation/gmc_diameter_six_macaulay_newton_atlas_thm2924.py
output: 05-knowledge/results/gmc_diameter_six_macaulay_newton_atlas_thm2924.out
script_sha256: b448ba5980ab8f0e538d27bef4581d66891d4bc7b31efe5c29c65d6367396b43
output_sha256: 01af047c189a78cfdf8f15e286981eb2a1ba44c03fd5e2d193d66790bf168626
constructor_dependency_sha256: 9f146073cf2f953bbbfcaaf20f1f58dbb25a105316d97a486e80db2be62b9b2e
hash_basis: LF-normalized bytes
---

# THM-2924 -- diameter-six Macaulay--Newton atlas

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT
HOSTILE AUDIT.**

## 1. Statement

Let

```text
L:C[s] -> C,                         L(s^j)=j!.          (1)
```

Fix an integer `n>=0` and any

```text
B=(0,a,b,6),                    0<a<b<6.                (2)
```

For every nonzero

```text
H=sum_(d in B) c_d s^(n+d),                               (3)
```

at least one of

```text
L(H),                  L(H^2),                  L(H^3), L(H^4) (4)
```

is nonzero.  These are precisely the ten translated four-slot support
types of diameter exactly six.

For every support in `(2)`, there is also a polynomial with all four
coefficients nonzero whose first three moments vanish and whose fourth
does not.  Thus four is the exact detection bound with full support.

## 2. Mean elimination and the width-six denominator

Put

```text
f_j=s^j/j!,                         L(f_j)=1,             (5)
```

and eliminate the mean by writing

```text
H=x(f_n-f_(n+6))
 +y(f_(n+a)-f_(n+6))
 +z(f_(n+b)-f_(n+6)).                                  (6)
```

After dividing the order-`m` moment by `L(f_n^m)`, its symmetric tensor
entry at offsets `d_1,...,d_m` is

```text
T_m(d_1,...,d_m;n)
 =(mn+1)_(sum_j d_j)/prod_j(n+1)_(d_j).                 (7)
```

The coefficient of an ordinary monomial with multiplicity vector
`alpha` is its signed top-replacement sum multiplied by the essential
ordinary-monomial factor

```text
m!/(alpha_0!alpha_1!alpha_2!).                          (8)
```

Call the resulting forms `Q,C,F` in `R=Q[x,y,z]`.  Exact polynomial
division, for every coefficient in all ten families, proves that

```text
D_2=prod_(j=1)^5(n+j),

D_3=[prod_(j=1)^5(n+j)^2](n+6),

D_4=[prod_(j=1)^5(n+j)^3](n+6)^2                       (9)
```

clear `Q,C,F`.  The scaled forms

```text
Q~=D_2 Q,                 C~=D_3 C,                 F~=D_4 F
```

have exact maximal coefficient degrees

```text
5,                          11,                         17. (10)
```

Every factor in `(9)` is positive for integer `n>=0`, so scaling
preserves projective common zeros.

## 3. The fixed Macaulay chart

Use the degree-seven map

```text
Phi_n:R_5 direct_sum R_4 direct_sum R_3 -> R_7,
Phi_n(A,B,D)=A Q~+B C~+D F~.                           (11)
```

In the fixed THM-2849 monomial order, select from its stored
`46`-by-`36` transpose the same rows used in THM-2921 and THM-2922:

```text
0,...,19;
21,...,29,35;
36,...,41.                                             (12)
```

The maximal minor uses `20Q+10C+6F`.  If its determinant is `P_B(n)`,
then `(10)` gives

```text
deg P_B <=20*5+10*11+6*17=312.                        (13)
```

Thus one Pluecker chart now survives every four-slot support of widths
four, five and six.

## 4. Exact Gregory--Newton certificate

The companion evaluates each determinant at the complete prefix needed
by `(13)`, plus one degree-exhaustion cell.  The minimal successful
Newton bases and exceptional signs group as follows:

| first positive offset `a` | families | Newton base | signs below base |
|---:|---:|---:|:---:|
| `1` | `4` | `5` | `-----` |
| `2` | `3` | `2` | `--` |
| `3,4` | `3` | `0` | none |

For every family, all `313` coefficients

```text
Delta^k P_B(r),                    0<=k<=312,           (14)
```

are strictly positive at the displayed base `r`, and

```text
Delta^313 P_B(r)=0.                                      (15)
```

Equation `(15)` independently exhausts the degree invoice.  Newton
interpolation now gives

```text
P_B(n)=sum_(k=0)^312 Delta^k P_B(r) binom(n-r,k)>0,
                                                     n>=r. (16)
```

Every displayed exceptional low value is strictly negative, hence
nonzero.  Therefore

```text
P_B(n)!=0                  for every B in (2), n>=0.   (17)
```

The exact evidence contains all `3,130` positive coefficients, proves
the symbolic denominator divisions, and compares the selected minor
against an independently constructed four-variable multinomial moment
map modulo `1000003` at seven depths in every family: `70` exact
cross-checks.

## 5. Projective and Gaussian consequences

The nonzero maximal minor makes `Phi_n` surjective onto `R_7`.  A
projective common zero of `Q~,C~,F~` would therefore annihilate every
degree-seven form, including a seventh power of one nonzero coordinate.
This proves `(4)` after restoring the eliminated mean.

THM-2173 produces a nonzero member of each four-dimensional support
whose first three moments vanish.  PROVED THM-2824 excludes such a
witness from every three-slot coordinate face, so it has full support;
`(4)` makes its fourth moment nonzero.

For `n>=1`, write `H=s h(s)`, set `s=ZW`, and put

```text
P=alpha W+Z h(ZW),                       alpha!=0.      (18)
```

Charge balance gives

```text
E[P^(2j+1)]=0,
E[P^(2j)]=binom(2j,j)alpha^j L(H^j).                  (19)
```

The associated genuinely two-charge five-monomial envelope therefore
has exact uniform Gaussian detection depth eight.  Its sharp witness
has moments one through seven zero and

```text
E[P^8]=70 alpha^4 L(H^4)!=0.                           (20)
```

## 6. Structural continuation and boundary

Widths four through six support the general denominator candidate

```text
D_(m,M)
 =[prod_(j=1)^(M-1)(n+j)^(m-1)](n+M)^(m-2),          (21)
```

the row-degree candidate

```text
deg(D_(m,M)L(H^m)) <=(m-1)M-1,                        (22)
```

and hence the fixed-chart bound

```text
deg P_B <=58M-36.                                      (23)
```

None of `(21)--(23)` is asserted here for arbitrary `M`.  In particular,
finite continuation at larger widths cannot by itself explain why one
Pluecker coordinate should remain nonzero.  The new structural target
is a valuation proof of `(21)` plus a planar-network, flagged-Schur, or
equivalent total-positivity realization of the shifted determinant.
Alternate maximal minors must be retained if this chart crosses a wall.

This theorem concerns only four slots of diameter exactly six and the
first four factorial moments.  It does not prove arbitrary SFC(4),
arbitrary-charge GMC(2), a shifted moment window, or a universal
multiplier selector.

## 7. Exact verification

The companion hash-pins the audited THM-2922 constructor, rebuilds
`(9)` independently, uses explicit `require` checks under ordinary and
optimized Python, and records LF-normalized hashes.

Run

```text
python 04-computation/gmc_diameter_six_macaulay_newton_atlas_thm2924.py
python -O 04-computation/gmc_diameter_six_macaulay_newton_atlas_thm2924.py
```

Normal and optimized executions byte-match the stored output with the
declared hashes.

**QED (candidate pending independent hostile audit).**
