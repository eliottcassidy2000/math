---
id: THM-2940
title: "Consecutive four-slot Macaulay--Newton short closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  n>=0, first-window SFC(4) holds on the
  translated consecutive support {n,n+1,n+2,n+3}.  Exact width-three
  denominator clearing and one fixed degree-seven Macaulay minor reduce
  the depth ray to 139 strictly positive Gregory--Newton coefficients.
  Arbitrary three-slot detection makes the four-moment bound sharp with
  full support; at positive depth the corresponding five-monomial
  two-charge Gaussian bound is exactly eight.
source: codex-gmc-holotopy-extension-2026-07-29
depends_on:
  - THM-2173-sparse-projective-factorial-moment-floor
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
related:
  - THM-2849-four-slot-first-window-macaulay-box
  - THM-2908-consecutive-four-slot-projective-resultant-closure
  - THM-2921-diameter-four-nonconsecutive-macaulay-newton-closure
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
script: 04-computation/gmc_consecutive_four_slot_macaulay_newton_thm2940.py
output: 05-knowledge/results/gmc_consecutive_four_slot_macaulay_newton_thm2940.out
script_sha256: 03ef73b730aea3264db9fa7f1cb9e1927ba97b6ab4b6a9b09c9a8023df194637
output_sha256: ef3431a702fd829ed8cc7a657fc511ee4ca06c52d71223a66d0c6c6703669502
constructor_dependency_sha256: 42e9b5ceddd677d1f2601a9d5d668c9437281596b65999ddcb8549d4e0b9bf64
hash_basis: LF-normalized bytes
---

# THM-2940 -- consecutive four-slot Macaulay--Newton short closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

Let

```text
L:C[s] -> C,                         L(s^j)=j!.          (1)
```

For every integer `n>=0` and every nonzero polynomial

```text
H=c_0 s^n+c_1 s^(n+1)+c_2 s^(n+2)+c_3 s^(n+3),        (2)
```

at least one of

```text
L(H),                 L(H^2),                 L(H^3), L(H^4) (3)
```

is nonzero.  Equivalently, first-window SFC(4) holds on every
translated consecutive four-slot support.

The bound is exact: for each `n` there is an `H` with all four
coefficients nonzero such that its first three moments vanish and its
fourth does not.

## 2. Mean elimination and ordinary coefficients

Put

```text
f_a=s^a/a!,                         L(f_a)=1.            (4)
```

Factorial rescaling of coordinates preserves projective nullity.  Every
vector in the mean-zero hyperplane has a unique expression

```text
H=x(f_n-f_(n+3))
 +y(f_(n+1)-f_(n+3))
 +z(f_(n+2)-f_(n+3)).                                (5)
```

It is therefore enough to prove that the homogeneous forms

```text
Q=L(H^2)/L(f_n^2),
C=L(H^3)/L(f_n^3),
F=L(H^4)/L(f_n^4)                                     (6)
```

have no common nonzero projective zero in `P^2`.

For offsets `d_1,...,d_m`, direct factorial normalization gives

```text
T_m(d_1,...,d_m;n)
 =L(prod_j f_(n+d_j))/L(f_n^m)
 =(mn+1)_(sum_j d_j)/prod_j(n+1)_(d_j).               (7)
```

If `alpha=(alpha_0,alpha_1,alpha_2)` has total degree `m`, repeat
direction `j` exactly `alpha_j` times and expand the signed `2^m`
choices replacing any selected direction by the top offset `3`.  The
ordinary coefficient of `x^alpha_0 y^alpha_1 z^alpha_2` is that signed
sum multiplied by

```text
m!/(alpha_0!alpha_1!alpha_2!).                         (8)
```

The factor `(8)` is essential; MISTAKE-323 records the false
divided-power constructor that omitted it.

## 3. Width-three denominator clearing

Exact division in `Z[n]` clears every coefficient of `(6)` with

```text
D_2=(n+1)(n+2),

D_3=(n+1)^2(n+2)^2(n+3),

D_4=(n+1)^3(n+2)^3(n+3)^2.                            (9)
```

This is also the `M=3` case of PROVED THM-2925, but the exact companion
checks every division here directly, so THM-2925 is related rather than
a dependency.  The coefficient degrees of

```text
Q~=D_2 Q,                   C~=D_3 C,              F~=D_4 F
```

are at most, and in fact attain,

```text
2,                            5,                        8. (10)
```

Every factor in `(9)` is nonzero for integer `n>=0`, so this scaling
does not change the common projective zero set.

## 4. One fixed Macaulay chart

Let `R=Q[x,y,z]` and write `R_d` for its degree-`d` part.  Use the
degree-seven Macaulay map

```text
Phi_n:R_5 direct_sum R_4 direct_sum R_3 -> R_7,
Phi_n(A,B,D)=A Q~+B C~+D F~.                          (11)
```

In canonical monomial order, its stored transpose has `46` rows and
`36` columns.  Select rows

```text
0,...,19;
21,...,29,35;
36,...,41.                                             (12)
```

The resulting square minor uses `20` quadratic, `10` cubic, and `6`
quartic rows.  Call its determinant `P(n)`.  Equations `(10),(12)`
give the rigorous degree invoice

```text
deg P<=20*2+10*5+6*8=138.                             (13)
```

## 5. Gregory--Newton positivity

The companion evaluates `P(n)` exactly at `n=0,...,138` and computes
all forward differences.  It proves

```text
Delta^k P(0)>0                    for 0<=k<=138.       (14)
```

The complete vector in `(14)` has SHA-256 digest

```text
72ba0779c61d5a0c7147e0a2cc01a092eb456dc9a7d128779aa9776c5b2468e3.
                                                               (15)
```

Since a polynomial of degree at most `138` has the exact expansion

```text
P(n)=sum_(k=0)^138 Delta^k P(0) binom(n,k),            (16)
```

equations `(14),(16)` imply

```text
P(n)>0                            for every integer n>=0. (17)
```

This is a single Pluecker chart on the entire consecutive depth ray.

## 6. Projective consequence

Equation `(17)` gives `rank Phi_n=36`, hence

```text
(Q~,C~,F~)_7=R_7.                                     (18)
```

If `[x:y:z]` were a common zero of the three forms, every element on the
left of `(18)` would vanish there.  But some coordinate is nonzero and
its seventh power belongs to `R_7`, a contradiction.  Thus `(6)` has
no common nonzero projective zero.  Together with `(5)`, this proves
`(3)`.

## 7. Sharpness and the Gaussian lift

For this fixed four-slot envelope, PROVED THM-2173 supplies a nonzero
`H` with

```text
L(H)=L(H^2)=L(H^3)=0.                                 (19)
```

It has full support: otherwise it lies on a three-slot coordinate face,
while PROVED THM-2824 detects every arbitrary three-slot polynomial
within its first three moments.  The result above forces

```text
L(H^4)!=0.                                            (20)
```

This proves exact full-support four-moment sharpness.

If `n>=1`, write `H=s h(s)`, put `s=ZW` for a standard complex Gaussian
`Z` with `W=conj(Z)`, and choose `alpha!=0`.  For

```text
P_G=alpha W+Z h(ZW),                                  (21)
```

charge balance gives

```text
E[P_G^(2j+1)]=0,
E[P_G^(2j)]=binom(2j,j) alpha^j L(H^j).               (22)
```

Thus every nonzero radial part in this two-charge envelope is detected
by one of moments `2,4,6,8`.  Applying `(22)` to `(19),(20)` yields a
five-monomial, full-radial-support polynomial whose moments one through
seven vanish and

```text
E[P_G^8]=70 alpha^4 L(H^4)!=0.                        (23)
```

The exact uniform Gaussian detection depth in this positive-depth
envelope is therefore eight.

## 8. Relation to THM-2908 and scope

Candidate THM-2908 develops a much larger moving-plane/projective
resultant atlas for the same consecutive SFC(4) headline, retaining
geometric information about its exceptional charts.  It is not used
here.  This theorem is the proved fixed-chart replacement dependency:
degree `138` and `139` Newton coefficients rather than the degree-`2804`
candidate projective eliminant.

The invariant is surjectivity of `(11)`, not the selected row address.
Here one chart remains nonzero on the whole depth ray.  At larger widths
the same chart can cross a wall; that signals a chart transition, not
projective rank loss.

The proved candidate scope is exactly

```text
slots:       four;
offsets:     (0,1,2,3);
window:      moments one through four;
depth:       every integer n>=0;
Gaussian:    n>=1, two charges, five monomials.        (24)
```

No arbitrary four-slot support, arbitrary-width SFC(4), SFC(5),
full Strong Factorial Conjecture, arbitrary-charge GMC(2), or
Jacobian-conjecture conclusion is asserted.

## 9. Exact verification

The exact companion hash-pins the independently audited
ordinary-monomial constructor of THM-2921.  It then:

1. constructs every coefficient by `(7),(8)` and proves every division
   in `(9)` inside `Z[n]`;
2. checks the exact degrees `(10)` and all `279` symbolic/numeric
   coefficient specializations at depths `0,...,8`;
3. evaluates the fixed determinant at all `139` required depths and
   verifies every coefficient in `(14)` is positive;
4. checks the original unscaled depth-zero determinant against the
   scaled determinant; and
5. independently reconstructs all three forms by direct four-variable
   multinomial expansion, reproducing the selected determinant modulo
   `1000003` at depths `0,1,2,7,31,97,197`.

The last route gives seven exact cross-constructor minor checks.  The
mixed-quadratic `2:1` hostile freezes the MISTAKE-323 boundary.  Run

```text
python 04-computation/gmc_consecutive_four_slot_macaulay_newton_thm2940.py
python -O 04-computation/gmc_consecutive_four_slot_macaulay_newton_thm2940.py
```

Normal and optimized executions byte-match the stored output with the
declared LF-normalized hashes.

An independent hostile audit rederived the width-three denominators and
row degrees, the degree-`138` invoice, the base-zero
Gregory--Newton argument, Macaulay surjectivity and projective emptiness,
the THM-2173/2824 full-support step, and the Gaussian charge formula.  It
also replayed normal, optimized, and stored output and reproduced the
declared LF-normalized hashes.  Its only requested repair was the
truth-surface qualification of still-candidate THM-2908 above; no
mathematical or exact-evidence defect remained.

**QED.**
