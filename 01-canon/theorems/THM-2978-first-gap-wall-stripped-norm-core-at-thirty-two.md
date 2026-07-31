---
id: THM-2978
title: "First-gap wall-stripped norm core at thirty-two"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. For support
  {n,n+1,n+2,n+32}, the two THM-2969 Macaulay charts
  have the same explicit flag and primitive wall-stripped pure-resultant
  core. The degree-660 core is dense, coefficient-positive, PF2, and has no
  integer wall root. Width 32 has no new pure integer-wall resonance or exceptional
  Smith correction. Width 33 and arbitrary radial GMC(2) are outside scope.
source: codex-gmc-first-gap-norm-core-m32-2026-07-30
depends_on:
  - THM-2969-first-gap-wall-stripped-resultant-norm-core-atlas
  - THM-2973-first-gap-wall-stripped-norm-core-continuation-through-thirty-one
  - THM-2964-general-pure-factorial-moment-resonance-ladder
related:
  - THM-2942-macaulay-extraneous-flag-factor-and-pluecker-mutation
  - THM-2943-width-seven-eight-two-chart-macaulay-resultant-closure
script: 04-computation/gmc_first_gap_wall_stripped_norm_core_at_thirty_two_thm2978.py
output: 05-knowledge/results/gmc_first_gap_wall_stripped_norm_core_at_thirty_two_thm2978.out
script_sha256: 4787aac1159f8ee9bdf7ad27df9654c3ed46ab32060d38965ba691b662c132cb
output_sha256: 32e0ac71ffe557868238a9e1cf8e019864f6260675093a9cb45bfeb2e24c7f76
hash_basis: LF-normalized bytes
---

# THM-2978 -- first-gap wall-stripped norm core at thirty-two

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## Statement

Let `n>=0` be an integer and let `H` be a nonzero polynomial supported on

```text
{n,n+1,n+2,n+32}.                                         (1)
```

Use the original and stable-mutated full Macaulay charts of THM-2969,
translated to normalized support `(0,1,2,32)`. Their specialized
determinants have common flag

```text
H_32(n)=(n+32) product_(r=3)^16 (n+r).                    (2)
```

Let `q_32,c_32,f_32,B_32,E_32` be the pure quadratic, pure cubic, pure
quartic, baseline local Smith, and seam factors in the notation of THM-2969.
Exact division of the common chart gcd by

```text
q_32^5 c_32 B_32 E_32                                    (3)
```

leaves, up to a nonzero rational unit, a primitive polynomial `N_32(n)`.
Then:

1. `deg N_32=660`, its primitive coefficient tuple has SHA-256

   ```text
   3dc8b52265fd9b96d92851cc1f69df9f987c681843f61b6893ac68506edd323e;
   ```

2. every coefficient is strictly positive, there are no internal zeros, and
   the coefficient sequence is PF2;
3. none of `q_32,c_32,f_32,N_32` has a root `-r` for
   `1<=r<=32`, and no exceptional width-32 Smith correction is required;
4. the two calibrated chart polynomials are coprime in `Q[n]`;
5. the association `R_32 ~ W_32 N_32` over `Q[n]`, where
   `W_32=(B_32 E_32)/H_32`, implies that the first four factorial moments of
   `H` cannot all vanish.

Thus THM-2973's finite first-gap nonvanishing range extends from `M=31` to
`M=32`. No statement for `M=33`, arbitrary width, a fixed-prime gate, or
arbitrary radial coefficients is asserted.

## Proof

### 1. No new pure or exceptional wall

THM-2964 classifies the pure order-`k` integer walls by

```text
M=(k+1)d+1,        r=kd+1,        (k-1)d even.            (4)
```

For `M=32`, none of `k=2,3,4` satisfies `(4)`. Direct exact division of the
three denominator-cleared coefficients confirms

```text
q_32 roots: none,       c_32 roots: none,
f_32 roots: none        on {-1,...,-32}.                  (5)
```

The full local division succeeds with THM-2969's baseline Smith factor:
there is no additional sporadic correction at width 32. Separately, exact
division of the primitive core by every `n+r`, `1<=r<=32`, finds no residual
integer wall. Positivity alone would not prove this last assertion; it is a
distinct exact check.

### 2. Full-chart and wall invoices

The companion imports the immutable THM-2969 engine by LF hash, reconstructs
both degree-`1820` determinant interpolants from exact integer evaluations,
and verifies each chart at three outside-grid depths. The common flag `(2)`
has degree `15`. The exact degree invoice is

```text
deg B=693,      deg E=108,      deg H=15,      deg W=786,
deg R=1446,     deg raw-gcd=1679,
deg removed divisor=1019,       deg N=660.                (6)
```

The complementary factor has half-wall order six and integer-wall orders

```text
r=1:6, r=2:24, r=3..10:28, r=11..16:26,
r=17..21:24, r=22..31:23, r=32:20.                       (7)
```

Every order is nonnegative, so `(B_32 E_32)/H_32` is polynomial. Write
`P~Q` when two polynomials are associated over `Q[n]`. As in THM-2973, the
common-factor relations are

```text
G~q^5 c R H~q^5 c B E N,       R~W N.                   (8)
```

Consequently `(3)` removes precisely the prescribed wall divisor and leaves
the primitive pure-resultant core `N_32`.

### 3. Positivity, PF2, and calibrated charts

Primitive normalization gives the degree and digest in the statement. An
exact scan of all `661` coefficients finds every coefficient positive, and
every adjacent minor

```text
a_j^2-a_(j-1)a_(j+1)
```

is nonnegative. Hence the coefficient sequence is PF2. The independently
calibrated original and mutated chart polynomials both have degree `141`,
and their exact gcd has degree zero.

This coprimality is a formal chart-calibration statement. It is neither the
moment implication nor a fixed-prime assertion.

### 4. Resultant-to-moment implication

For real `n>=0`, coefficient positivity gives `N_32(n)>0`. Formula `(7)` and
the positive half-wall show `W_32(n)>0`; hence `(8)` gives `R_32(n)!=0`.
By the defining property of the denominator-cleared ternary resultant, the
quadratic, cubic, and quartic factorial-moment forms have no common
projective zero. Restoring the mean eliminated in their construction shows
that

```text
L(H), L(H^2), L(H^3), L(H^4)
```

cannot all vanish, exactly as in THM-2969 and THM-2973.

## Exact evidence

The canonical reproduction is

```text
python 04-computation/gmc_first_gap_wall_stripped_norm_core_at_thirty_two_thm2978.py --output .scratch/thm2978.normal.out
python -O 04-computation/gmc_first_gap_wall_stripped_norm_core_at_thirty_two_thm2978.py --output .scratch/thm2978.opt.out
```

The companion checks the two full determinant interpolants, six outside-grid
determinant evaluations, the flag and degree invoices, the complete pure and
core integer-wall censuses, every complementary-wall order, primitive
normalization, strict positivity, PF2, calibrated chart coprimality, and the
associated common-factor relations inherited from the audited THM-2969 engine.

The stored transcript is `1,230` UTF-8/LF bytes. Its frozen LF hashes are

```text
script  4787aac1159f8ee9bdf7ad27df9654c3ed46ab32060d38965ba691b662c132cb
output  32e0ac71ffe557868238a9e1cf8e019864f6260675093a9cb45bfeb2e24c7f76
```

The independent hostile audit rederived the pure-wall census from THM-2964
and by direct division, checked the complete degree and complementary-wall
invoices, separately verified the core integer-root scan, positivity/PF2,
formal calibrated-chart gcd, associated common-factor relations, and the
resultant-to-moment implication. Its optimized replay byte-matches the final
normal and stored outputs and recovers both frozen hashes.

Width `33`, transition maps in `M`, arbitrary interior support `(0,a,b,M)`,
and arbitrary radial GMC(2) remain open.

**QED.**
