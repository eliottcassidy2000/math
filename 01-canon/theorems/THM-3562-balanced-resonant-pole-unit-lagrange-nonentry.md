---
id: THM-3562
title: "Balanced resonant pole-unit constancy and all-passport Lagrange nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In every
  normalized nonsplit polynomial exact-square-prefix Faber bank with
  R=3k+2>=8, no nonconstant balanced response enters, for any pole passport.
  The theorem is chart-local: unbalanced responses, nonpolynomial prefixes,
  other charts, and JC(2) remain open.
source: root-2026-08-18-planar-jacobian-counterexample-constructions
audit: >
  An independent hostile audit rederived the response factorization, the
  off-divisor zero exclusion, the resonant source divisor, the common
  pure-q pole unit, and the Lagrange partial-fraction contradiction.  It
  also checked the one-pole split boundary and the asymmetric R=8
  finite-pole control.  Ordinary and optimized replays are byte-identical to
  the stored transcript; dependency, universe, hash, documentation, and
  assertion audits pass.
depends_on:
  - THM-2796-balanced-response-stieltjes-pade-normal-form-and-one-double-zero-classification
  - THM-2827-uniform-pole-order-faber-nonresonance-atlas-and-double-pole-exclusion
  - THM-3133-common-simple-zero-faber-exclusion-and-odd-bipole-boundary
related:
  - THM-3151-resonant-odd-bipole-equality-cell-nonentry-and-degree-floor
  - MISTAKE-317
script: 04-computation/jc_balanced_resonant_pole_unit_thm3562.py
output: 05-knowledge/results/jc_balanced_resonant_pole_unit_thm3562.out
script_sha256: 2153d2ecfc65f59d24fc1be09f7af66d1637762df0f82c04dddfc66e5915160d
output_sha256: 4692b4fba5243d67b03887fd2abedab7af9aca6334c36354afe034bf4212ce24
hash_basis: LF-normalized bytes
---

# THM-3562 -- balanced resonant pole-unit constancy and all-passport nonentry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3151 excludes only the equality passport `(D,D)` by an infinity fan.
That route becomes awkward when the source coordinate `s_F` has finite
poles.  The present proof instead transports the unique local pure-`q`
coefficient to every response pole.  It forces one polynomial to take the
same nonzero value on the full pole set, where the response first integral
and the elementary Lagrange identity are incompatible.  Finite poles of
`s_F` never have to be removed.

## 1. Statement and inherited response data

Work over an algebraically closed field of characteristic zero.  Fix

```text
R=3k+2>=8,                    D=4k+3,                    k>=2.     (1)
```

Suppose that a nonconstant balanced response enters a normalized nonsplit
polynomial exact-square-prefix Faber bank

```text
Q=E_(4R-2)+sum_(j=1)^(R-2) a_j E_(4j-2).                    (2)
```

THM-3133 excludes all simple response zeros.  In THM-2796's factor normal
form this says `S=1`.  THM-2827 then forces every pole part to be divisible
by `D`.  Write

```text
p_j=D delta_j,             delta_j>=1,
W=product_j (x-beta_j),    P=product_j (x-beta_j)^delta_j,
N=sum_j p_j=2e,

V=v W^2 P^D,              M=E W,                           (3)
```

where `v!=0`, `E` is monic squarefree of degree `e`, and `E,W` are
coprime.  The balanced response first integral is

```text
2W E' - E sum_j p_j W_j=C,          W_j=W/(x-beta_j),
C!=0.                                                        (4)
```

Use the inherited source coordinates

```text
T_F=A_src^2/V,
d_F=C_0-B_src^2/(4V),
s_F=A_src B_src/(2V)-E_0,                                  (5)
```

and the exact Faber sidecars

```text
phi_Q=0,                  psi_Q in C,
K_Q=T_F H_Q,             A_src K_Q=lambda M,       lambda!=0,
H_Q in C[T_F,s_F,T_F d_F].                                 (6)
```

The claim is that these data are inconsistent with genuine nonsplitting,
for every positive pole passport `(delta_1,...,delta_h)`.

## 2. The source divisor has no hidden zeros

Let `alpha` be a zero of `A_src` away from `V`, of order `a>=1`.  Since
`V` is a unit at `alpha`, all of `d_F,s_F,T_Fd_F` are regular there and

```text
ord_alpha(T_F)=2a,
ord_alpha(K_Q)>=2a,
ord_alpha(A_src K_Q)>=3a.                                  (7)
```

The right side of the last equation in `(6)` is `lambda M`.  Because
`M=EW` is squarefree, its order at `alpha` is either zero or one, contrary
to `(7)`.  Hence every zero of `A_src` lies among the `beta_j`.

At `beta_j`, the pair of response orders is

```text
(ord_beta_j V,ord_beta_j M)=(2+D delta_j,1).                (8)
```

The exact THM-2827 resonance ray therefore gives

```text
ord_beta_j(A_src)=1+(2k+1)delta_j.                          (9)
```

There are no other zeros, so for some `a_0!=0`, and with
`tau=a_0^2/v!=0`,

```text
A_src=a_0 W P^(2k+1),              T_F=tau/P.                (10)
```

Substitute `(3),(6),(10)` into `A_src T_F H_Q=lambda EW`.
Cancellation of `W` yields the global rational identities

```text
H_Q=(lambda/(a_0 tau)) E/P^(2k),
H_Q/T_F^(2k)=cE,
c=lambda/(a_0 tau^(2k+1))!=0.                              (11)
```

This step is where the global source divisor is load-bearing.  A local
valuation equality alone would not identify the regular quotient in
`(11)`.

## 3. The pure-q face becomes one common pole unit

At each `beta_j`, THM-2827's pure-`q` analysis has a unique top face in
`H_R`.  Every retained lower row has strictly larger valuation, and exact
coefficient extraction gives

```text
H_Q/T_F^(2k) mod (x-beta_j)
  =c_R,
c_R=4 binom(R-1/2,4k+3)!=0.                                (12)
```

The residue in `(12)` is independent of the point and of `delta_j`; all
dependence on the local leading coefficient is already contained in
`T_F^(2k)`.  Comparing `(11)` and `(12)` gives one common nonzero value

```text
E(beta_j)=e_0:=c_R/c!=0                 for every j.         (13)
```

This is stronger than the valuation statement used in the old equality-cell
proof.  It remains valid when `s_F` or `d_F` has a finite pole, because the
strict pure-`q` efficiency and its coefficient are precisely the inherited
THM-2827 sidecar.

## 4. The Lagrange identity closes every pole passport

Evaluate the response first integral `(4)` at `beta_j`.  Since
`W_j(beta_j)=W'(beta_j)`, equations `(4),(13)` imply

```text
-e_0 p_j W'(beta_j)=C,
p_j W'(beta_j)=q:=-C/e_0!=0                for every j.      (14)
```

If `h>=2`, the partial-fraction expansion of the monic polynomial `W` gives

```text
1/W(x)=sum_j 1/[W'(beta_j)(x-beta_j)],
sum_j 1/W'(beta_j)=0.                                      (15)
```

The last equality is the coefficient of `x^(-1)` at infinity: `1/W` is
`O(x^(-h))`.  But `(14)` gives

```text
sum_j 1/W'(beta_j)=sum_j p_j/q=N/q!=0,                      (16)
```

contradicting `(15)` in characteristic zero.

If `h=1`, balance gives `p_1=N=2e`, so the only pole part is even.  With
`S=1`, THM-2796's exact squareclass criterion says that this response is
split, contrary to the nonsplit hypothesis.  Thus neither `h>=2` nor
`h=1` can occur.  No nonconstant balanced response enters the normalized
nonsplit chart `(2)`.  QED.

## 5. First asymmetric finite-pole hostile control

The following **FINITE-EXACT control is constructed in this session**; it
is not an inherited theorem or a chart entrant.  It demonstrates why one
must not assume `s_F` is polynomial and why the pole-unit gate adds genuine
information.

At `R=8`, take `k=2`, `D=11`, pole points `(0,1)`, and

```text
(delta_0,delta_1)=(3,1),              (p_0,p_1)=(33,11),
P=x^3(x-1),                           W=x(x-1),

V=x^35(x-1)^13,
A_src=x^16(x-1)^6,
B_src=x^18(x-1)^7,
C_0=E_0=0.                                                    (17)
```

Then the exact resonance divisor is satisfied at both poles, while

```text
T_F=1/[x^3(x-1)],
d_F=-x(x-1)/4,
s_F=1/(2x).                                                   (18)
```

Thus `s_F` has a genuine finite pole.  Let

```text
E(x)=sum_(m=0)^22 (-1)^m binom(11/2,m)x^(22-m).              (19)
```

Direct exact expansion gives

```text
2x(x-1)E'-(44x-33)E=30705345/2^41,
E(0)=930465/2^41,
E(1)=-2791395/2^41=-3E(0).                                  (20)
```

The nonzero constant in `(20)` proves that `E` is squarefree and disjoint
from `W`.  Exact expansion additionally gives that `E^2-P^11` is squarefree,
disjoint from `EP`, and has degree `21`.  Thus these data supply a balanced
response with pole passport `(33,11)`, zero partition `2^22`, and third
partition `(23,1^21)`; the high part `23=44-21` is at infinity.  The
weighted first-integral evaluations agree:

```text
33 W'(0)E(0)=11 W'(1)E(1).
```

But `E(0)!=E(1)`, so the response fails the necessary common-pole-unit
condition `(13)`.  This is a hostile control for the new obstruction, not
evidence of Keller-chart entry.

## 6. Scope, equality boundary, and exact companion

The theorem strictly strengthens THM-3151 inside the same normalized
polynomial exact-square-prefix chart: all balanced passports are excluded,
not merely `(D,D)`.  It does not cover unbalanced responses, nonpolynomial
prefixes, other nonsplit presentations, or arbitrary Keller pairs.  It does
not prove the planar Jacobian conjecture.

The exact companion uses rational arithmetic only.  It checks:

1. resonance exponents and the regular quotient through `k=2..40` and
   `delta=1..12`;
2. the direct pure-`q` coefficient and retained-row gap through `k=2..40`;
3. the full Lagrange interpolation vector on three exact node packets for
   each `h=2..16`;
4. all positive balanced `delta` passports of even total at most ten; and
5. every identity in the asymmetric `R=8` control, including the clean
   factor typing and third-fibre degree.

Run

```text
python3 04-computation/jc_balanced_resonant_pole_unit_thm3562.py
python3 -O 04-computation/jc_balanced_resonant_pole_unit_thm3562.py
```

and compare both transcripts byte-for-byte with the declared output.  The
finite universes are verification sidecars; the proof is the symbolic
argument in Sections 2--4.
