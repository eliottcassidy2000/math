---
id: THM-2969
title: "First-gap wall-stripped resultant norm-core atlas through width twenty-six"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For supports {n,n+1,n+2,n+M}, 6<=M<=26, two stable full Macaulay
  charts have an explicit common half-wall and leave the same primitive
  wall-stripped pure-resultant core.  Every core is dense,
  coefficient-positive, and PF2, hence the first four factorial moments
  cannot vanish simultaneously.  This is a finite family theorem only:
  no arbitrary radial GMC(2), arbitrary width, or fixed-prime conclusion
  is claimed.
source: codex-gmc-first-gap-norm-core-atlas-2026-07-30
audit: >
  An independent hostile auditor accepted the M=23,...,26 extension,
  rechecked the simultaneous q/c walls at M=25, and required four packaging
  repairs: distinguish paired-depth controls from determinant evaluations,
  extend the +6 baseline range through target 26, broaden the seam sentence
  through M=26, and move the no-extrapolation boundary to M>=27.  All four
  repairs are installed in this theorem.
depends_on:
  - THM-2925-general-width-terminal-pole-cancellation-and-macaulay-degree-law
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
  - THM-2943-width-seven-eight-two-chart-macaulay-resultant-closure
  - THM-2960-local-smith-jet-fitting-barcode-and-negative-depth-chamber-atlas
  - THM-2964-general-pure-factorial-moment-resonance-ladder
related:
  - THM-2957-first-gap-width-fifteen-twenty-modular-depth-ladder
  - THM-2959-first-gap-width-twenty-one-twenty-four-modular-depth-continuation
script: 04-computation/gmc_first_gap_wall_stripped_norm_core_atlas_thm2969.py
output: 05-knowledge/results/gmc_first_gap_wall_stripped_norm_core_atlas_thm2969.out
script_sha256: 5be5df0fd058f436593339f78cd6240a47975082d22929c135ba811458753bf5
output_sha256: 7867d6dfc1d08053dc6563c95c5822a087bfab763edea2fba6eed183454eef6a
hash_basis: LF-normalized bytes
---

# THM-2969 -- first-gap wall-stripped resultant norm-core atlas

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

The reusable identities `(4)--(13)` are finite exact.  Widths
`M=21,...,26` extend the two-chart seam beyond THM-2960.  The next
wider-width question begins at `M=27` and is not used here.

This note isolates the common nonlinear carrier left by the two stable
Macaulay charts for the normalized first-gap supports

```text
(0,1,2,M),                    6 <= M <= 26.             (1)
```

Here `L:C[s]->C` is the factorial functional `L(s^j)=j!`; a translated
support `(1)` means the exponent set `{n,n+1,n+2,n+M}`.

It also records the strongest all-polynomial statement available from the
calibrated charts.  The width atlas is finite.  Nothing below extrapolates
the factorization, positivity, or degree formula to `M>=27`.

## 1. Two charts and their common flag

Use the notation of PROVED THM-2942/2943.  After mean elimination and the
positive denominator clearing, let `Q_M,C_M^form,F_M` be the ternary forms
of degrees two, three, and four.  To avoid confusing the cubic form with
the core defined below, write

```text
q_M = [x0^2] Q_M,                 c_M = [x0^3] C_M^form,
R_M = Res(Q_M,C_M^form,F_M).                              (2)
```

Let `P_M,A_M in Z[n]` be the original and stable-mutated full Macaulay
determinants.  The universal Pluecker factorization is

```text
P_M = q_M^6 c_M K_M R_M,
A_M = q_M^5 c_M L_M R_M,                                (3)
```

where `K_M,L_M` are the two explicit selected-basis flags of THM-2942.
Polynomials are considered up to multiplication by a nonzero rational
constant; `pp(f)` denotes the primitive associate with positive leading
coefficient.

For every width in `(1)`, exact gcd arithmetic proves

```text
H_M := pp(gcd(q_M K_M,L_M))
     = (n+M) prod_(r=3)^floor(M/2) (n+r),                (4)

G_M := pp(gcd(P_M,A_M))
     ~ q_M^5 c_M R_M H_M.                               (5)
```

Thus the extra common flag is not an unidentified nonlinear carrier: it
is the explicit half-wall `(4)`.  The identity is a statement over
`Q[n]`, not a finite-field inference.

## 2. Smith walls, the complementary seam, and the pure core

For `1<=r<=M`, put

```text
beta_M(r) =
  6,   r=1;
  21,  r=2;
  26,  3<=r and 3r<=M;
  24,  3r>M and 2r<=M;
  20,  2r>M and 3r<=2M;
  19,  3r>2M,                                             (6)

epsilon_M(r)=1_{(M,r) in {(11,9),(12,5),(21,17)}}.       (7)
```

Here `(11,9)` and `(21,17)` are the quartic pure-coefficient resonances
proved in THM-2960/2964; `(12,5)` is the unique matrix-only sporadic in
this finite range.  Define

```text
B_M=(2n+1)^5 prod_(r=1)^M (n+r)^(beta_M(r)+epsilon_M(r)),

E_M=(2n+1)(n+M)^2
    *prod_(r=2)^(M-1) (n+r)^(3 if 2r<=M else 4).          (8)
```

The finite exact atlas proves that `H_M | B_M E_M` and, with

```text
W_M=B_M E_M/H_M,                                         (9)
```

both divisions below are exact and give associated primitive
polynomials:

```text
N_M := pp(G_M/(q_M^5 c_M B_M E_M))
     = pp(R_M/W_M).                                      (10)
```

The companion transcript calls this polynomial `C_M`; the letter `N` is
used here only to keep it visually distinct from the cubic moment form.

Consequently the common-gcd normalization and the pure-resultant
normalization leave exactly the same wall-stripped nonlinear core.  There
is no residual selected-basis flag hidden in `N_M`.
For `M=21,...,26`, this is also a finite exact two-chart verification of
the seam `E_M`; THM-2960's canonical two-chart sidecar stopped at width
twenty even though its one-chart Smith atlas continued farther.

## 3. Positivity, PF2, and the exact degree law

Write the coefficients of `N_M`, in either order, as
`a_0,...,a_d`.  For every `6<=M<=26`, exact integer arithmetic proves

```text
a_i>0                                  for every i,
a_i^2 >= a_(i-1) a_(i+1)              for 1<=i<d.       (11)
```

Thus the core is dense, strictly coefficient-positive, and PF2
(log-concave).  In particular `N_M(n)>0` for every real `n>=0`.  Since
`W_M` is a product of positive linear factors on that ray, the pure
resultant `R_M` has no nonnegative real zero in this finite bank.
By the defining property of the ternary resultant, the denominator-cleared
quadratic, cubic, and quartic forms therefore have no common projective
zero.  Restoring the eliminated mean gives the unified finite consequence

```text
L(H),L(H^2),L(H^3),L(H^4) are not all zero             (11a)
```

for every nonzero `H` supported on
`{n,n+1,n+2,n+M}`, every integer `n>=0`, and every `6<=M<=26`.

The degrees are

```text
M:       6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26
deg N: 121 144 164 184 205 226 244 268 288 308 329 351 369 392 412 431 453 475 493 516 536,
                                                                  (12)
```

equivalently

```text
deg N_M = 23M-2 floor(M/3)-2 floor(M/2)-floor(2M/3)-3-e_M,
e_M=1_{M in {11,12,21}}.                               (13)
```

This formula is not fitted to the table.  Summing the four chambers in
`(6)` gives

```text
deg B_M=19M+2 floor(M/3)+4 floor(M/2)+floor(2M/3)-20+e_M,
deg E_M=4M-floor(M/2)-4,             deg H_M=floor(M/2)-1.
```

Subtracting `deg W_M` from the inherited exact law
`deg R_M=46M-26` gives `(13)`.

### 3.1 The simultaneous width-25 pure walls

Width 25 contains both a simple quadratic wall and a simple cubic wall:

```text
ord_(-17) q_25=1,      ord_(-17) c_25=0,
ord_(-19) q_25=0,      ord_(-19) c_25=1.              (13a)
```

The exact local invoices, in the order
`(q,c,B,E,H,W,raw gcd,removed divisor,core)`, are

```text
n=-17: (1,0,19,4,0,23,28,28,0),
n=-19: (0,1,19,4,0,23,24,24,0).                       (13b)
```

Thus both resonances are absorbed by the separately removed
`q_M^5 c_M` factor together with `B_M E_M`.  Neither is a new Smith
sporadic and neither survives in `N_25`.

After adding back `e_M`, the six-width step is exactly `124`; the exact
+6 check now runs from sources `M=6,...,20` to targets
`M=12,...,26`.  The three corrections in `(13)` are not numerical
noise: two are the pure quartic ladder and one is the matrix Smith sporadic.

## 4. The general calibrated-chart divisor lemma

The reusable part is not finite-width.  Let `U,V in Z[n]` be primitive,
nonzero, and coprime over `Q[n]`.  Then

```text
gcd(U(t),V(t)) divides |Res_n(U,V)|                  (14)
```

for every integer `t`, with the convention that the positive gcd is
taken.  Indeed, the adjugate of the Sylvester matrix gives

```text
A(n)U(n)+B(n)V(n)=Res_n(U,V)                           (15)
```

for some `A,B in Z[n]`; evaluation at `t` proves `(14)`.

Apply this to the calibrated charts

```text
U_M=pp(P_M/G_M),                    V_M=pp(A_M/G_M).    (16)
```

The exact atlas checks `gcd(U_M,V_M)=1` for every width `(1)`, so `(14)`
gives a finite prime-support bound for simultaneous divisibility of the
**calibrated** chart values at every integer depth.  At `M=6`, the
nonzero resultant in `(14)` has

```text
2380 bits, 717 decimal digits,
SHA256(decimal expansion)
 =68785d6e7d64e3a73786911729c1adad14650349b5672e69818b309a9902d3a1.
                                                                  (17)
```

This is the strongest uniform arithmetic consequence presently justified.
It is a bound for common primes of `U_M(t),V_M(t)`, not a fixed-prime
divisibility theorem for either original minor separately.

Conceptually, `Res(U_M,V_M)` is the determinant/Fitting order of the
Sylvester cokernel.  Thus the natural arithmetic continuation is a Smith or
mod-`p` analysis of that finite cokernel, not a coefficient-positivity
argument for `N_M`.  For a proposed prime gate `p`, the cheapest exact test
is whether `gcd(U_M,V_M)` remains one after reduction modulo `p`.

## 5. Sharp stopping boundaries

### 5.1 No short six-width shift recurrence

Exact gcds give

```text
gcd(N_13(n),N_7(n+s))=1                 for -13<=s<=13,
gcd(N_12(n),N_6(n+s))=1                 for -12<=s<=12. (18)
```

Thus the clean degree increment `124` does not come from divisibility by
a shifted earlier core, even in these two smallest six-width comparisons.
This refutes that proposed mechanism only; it does not rule out a more
complicated recurrence.

### 5.2 Positivity is archimedean, not a prime gate

Strict positivity and PF2 cannot by themselves yield a fixed prime bank.
The squarefree polynomial

```text
(n+1)(n+2)=n^2+3n+2
```

has dense positive, strictly log-concave coefficients, but its value at
`n=p-1` is divisible by the arbitrary prime `p`.  Hence `(11)` proves
nonvanishing on the physical ray; it does not control prime divisors of
core values.  Any next-prime Macaulay gate needs an additional congruence
or calibrated-chart input.

## 6. Exact universe and replay

The companion reconstructs both full degree-`58M-36` determinants for
all twenty-one widths, checks each interpolation at the next three integer
depths, and verifies `(3)`--`(13)` over `ZZ/QQ`.  This gives

```text
21 widths,
42 exact determinant interpolants,
63 outside-grid paired-depth controls,
126 individual outside-grid determinant evaluations,
21 common-gcd/resultant/core identities,
21 dense-positivity and full PF2 scans.                (19)
```

The four added widths contribute exactly 12 paired-depth controls, i.e.
24 individual determinant evaluations.

The full `M=6,...,26` record digests are

```text
core_record_digest=29bcdeade4601f818730957a6a326fa83c1cc1cb17df6189cda6c1082bbab255,
flag_record_digest=a73fcfc2ed9a41ba86176c8e6409ac3f86d2f69df2f17f7584f74ae055e4859b.
```

It also checks the two shift hostiles, polynomial coprimality of every
calibrated pair, the explicit `M=6` resultant datum, and sampled instances
of `(14)` at depths `0,...,50`.

The canonical replay is

```text
python 04-computation/gmc_first_gap_wall_stripped_norm_core_atlas_thm2969.py --output .scratch/gmc_first_gap_wall_stripped_norm_core_atlas_thm2969.normal.out
python -O 04-computation/gmc_first_gap_wall_stripped_norm_core_atlas_thm2969.py --output .scratch/gmc_first_gap_wall_stripped_norm_core_atlas_thm2969.opt.out
```

Both modes must exit zero and byte-match in UTF-8/LF.  The canonical paths are

```text
04-computation/gmc_first_gap_wall_stripped_norm_core_atlas_thm2969.py
05-knowledge/results/gmc_first_gap_wall_stripped_norm_core_atlas_thm2969.out
```

The two replay transcripts and the stored output are byte-identical
UTF-8/LF files of 5,871 bytes.  Their frozen LF hashes are

```text
script=5be5df0fd058f436593339f78cd6240a47975082d22929c135ba811458753bf5,
output=7867d6dfc1d08053dc6563c95c5822a087bfab763edea2fba6eed183454eef6a.
```

The exact scope is `(1)`: supports
`{n,n+1,n+2,n+M}`, integer `n>=0`, and `6<=M<=26`.  The only GMC
consequence is finite-family nonvanishing of at least one of the first four
factorial moments.  The theorem gives no arbitrary-width core formula, no
short recurrence, no fixed prime bank for uncalibrated minors, and no
closure of arbitrary radial GMC(2).

The known hostile boundary is the next quartic pure wall
`(M,r)=(31,25)`, proved by THM-2964.  Any extension through width 31 must
modify and re-audit the current correction support; extrapolating the
`M<=26` formula verbatim is forbidden.

**QED (finite-exact scope).**
