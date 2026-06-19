---
id: HYP-2632
title: LRC(14) repeated-residue character kernel - the two-large tail reduces to a small chi_7/affine integer table
status: PARTIALLY-CONFIRMED
source: codex-2026-06-19-S23
depends_on:
  - HYP-2630
  - HYP-2631
  - HYP-2626
  - HYP-2624
  - HYP-2617
related:
  - HYP-2614
  - HYP-2619
  - THM-538
  - OPEN-Q-108
---

# HYP-2632 - LRC(14) Repeated-Residue Character Kernel

## Claim

The next sharp target after HYP-2630 is now a finite signed character kernel,
not another coimage census.  For the `d=9` support-six coefficient, define

```text
S_d(a_1,...,a_6)
 = sum_{r in (F_7^*)^6, a.r=0} C_d(r).
```

Then the repeated-residue packets surviving after height-2 wall deletion have
a small integer table in units

```text
U = 0.00955648353590534...
```

```text
4+2 packet (1,1,1,1,a,a):
  a=2,4      chi_7(a)=+1   S_9 = -25 U
  a=3,5,6    chi_7(a)=-1   S_9 = -18 U
  a=0                       S_9 =  -4 U
  a=1                       S_9 =   0

4+1+1 packet (1,1,1,1,a,b), with a,b in {0,2,3,4,5,6}:
  6 high entries:  +8 U
  6 low entries:   +1 U
  3 zero entries:   0
```

The observed HYP-2630 QR/NQR split is therefore not a floating-point census
artifact.  It is a literal Legendre component of a signed finite kernel:

```text
2*S_9(1,1,1,1,a,a)/U = -43 - 7*chi_7(a),  a=2,3,4,5,6.
```

For `4+1+1`, the earlier HYP-2630 multiplicative signature was not complete.
The hidden extra coordinate is affine:

```text
a+b = 2 mod 7  ->  S_9(1,1,1,1,a,b)=0.
```

Specifically the zero lane is `(0,2)`, `(3,6)`, `(4,5)`.  Short-signature
collisions force this extra coordinate:

```text
(-1,-1, 1,-1): (3,6)->0, (5,6)->8U
( 1,-1,-1,-1): (2,6)->U, (4,5)->0
```

Off this zero lane, the high/low split is again Legendre.  Define

```text
Q(a,b) = ab*(1+3(a+b))-1 mod 7.
```

Then, for unordered `a,b in {0,2,3,4,5,6}` with `a+b != 2`,

```text
S_9/U = 8  iff  chi_7(Q(a,b)) = +1,
S_9/U = 1  iff  chi_7(Q(a,b)) = -1.
```

So the two-large theorem must retain the `chi_7` phase, the affine zero lane,
and this secondary Legendre selector.

## Computation

Script:

- `04-computation/lrc14_repeated_character_kernel_codex_s23.py`
- output: `05-knowledge/results/lrc14_repeated_character_kernel_codex_s23.out`

The script verifies the finite additive-Fourier identity over all `159`
projective support-six coimage classes:

```text
S_d(a) = sum_{a.r=0} C_d(r)
       = (1/7) sum_{t in F_7} C_hat(t a).
```

Numerical verification:

```text
max identity error: 1.100e-14
max imaginary drift: 2.858e-15
additive transform cache entries: 955
```

## Exact Findings

The `4+2` table is:

```text
a  chi(a)  S_9/U  status
0       0     -4  tail
1       1      0  zero
2       1    -25  tail
3      -1    -18  tail
4       1    -25  tail
5      -1    -18  tail
6      -1    -18  tail
```

The `4+1+1` table over unordered pairs in `{0,2,3,4,5,6}` is:

```text
pair    S_9/U  chi_7(Q)
(0,2)      0      -1
(0,3)      1      -1
(0,4)      1      -1
(0,5)      1      -1
(0,6)      1      -1
(2,3)      8      +1
(2,4)      8      +1
(2,5)      8      +1
(2,6)      1      -1
(3,4)      8      +1
(3,5)      1      -1
(3,6)      0      -1
(4,5)      0      -1
(4,6)      8      +1
(5,6)      8      +1
```

Signed repeated-tail ledger:

```text
4+2 signed packet sum:      -108 U
4+2 absolute packet mass:    108 U
4+1+1 signed packet sum:      54 U
4+1+1 absolute packet mass:   54 U
combined signed sum:         -54 U
combined absolute mass:      162 U
absolute/net ratio:             3
```

This is the first place where the repeated-residue tail starts to look like a
usable signed proof object.  The absolute packet mass already overstates the
net finite kernel by a factor of `3` before any analytic reciprocal summation
is applied.

## Proof Target

The new theorem-shaped target is:

```text
height-2 wall deletion
-> finite additive-Fourier kernel
-> chi_7 + affine zero lane + Legendre selector Q(a,b)
-> two-large reciprocal hyperplane sum.
```

The analytic estimate should be a cotangent/Dedekind or summation-by-parts
bound for the two-large reciprocal hyperplane sums weighted by the signed
integer table above.  A proof that only uses the HYP-2630 four-character
signature is too coarse, because it misses the affine zero lane `a+b=2` and
the selector `Q(a,b)=ab*(1+3(a+b))-1`.

## Tournament Analysis

Candidate vertices included runners, gaps, fixed circle sections, section
boundaries, wall-crossing events, speed residues, cover arcs, additive Fourier
modes, quadratic characters, Jacobi/cross-ratio signatures, and proof
obligations.

Chosen Hamiltonian path:

```text
additive_fourier_kernel
> affine_line_a_plus_b
> quadratic_character_phase
> height2_wall_deletion
> euler_copy_capacity
> raw_coimage_enumeration
> raw_runner_vertices
```

Fingerprints:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3_cycles = 0
SCC_sizes = [1,1,1,1,1,1,1]
```

The quotient preserves the signed `d=9` support-six coimage coefficient after
height-2 wall deletion.  It destroys exact witness times, raw support
identities, and runner labels.

## Status

Partially confirmed by exact finite enumeration and additive-Fourier
verification.  LRC(14) remains open.  Next step: attach this signed integer
kernel to the actual two-large reciprocal hyperplane sums and prove a uniform
tail bound that exploits the `-108U + 54U` cancellation rather than the
`162U` absolute mass.
