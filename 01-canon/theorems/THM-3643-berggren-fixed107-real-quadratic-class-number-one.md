---
id: THM-3643
title: "Berggren fixed-107 real quadratic class number one"
status: >
  PROVED + FINITE-EXACT + RIGOROUS-ARB; PENDING INDEPENDENT HOSTILE AUDIT.
  The ordinary ideal class number of Q(sqrt(1225041)) is exactly one.
  This is a quadratic descent sidecar only: it does not prove that the
  fixed-107 Mordell curve has rank exactly two.
source: kps-s189 / fixed-107 Mordell descent sidecar, 2026-08-21
depends_on:
  - THM-3620-berggren-fixed-107-mordell-rank-two-and-local-collision
related:
  - THM-3640-berggren-positive-cube-slope-atlas-through-401
script: 04-computation/berggren_fixed107_real_quadratic_class_number_one_thm3643.py
output: 05-knowledge/results/berggren_fixed107_real_quadratic_class_number_one_thm3643.out
script_sha256: 95878dc28b5910e8330bbdd0b67db4759218712d4ab95d5d3c348d5c96e008ca
output_sha256: b607cfb55dcf18b7b9eb112ec041dcb0c6128784e91e2a12b5979e2a02c8cea9
semantic_sha256: f7fefe14d31997d3e6cd2744695fcc1da368bca544a52e7fdc1d8609d028f001
hash_basis: raw LF bytes for files; canonical JSON for semantic ledger
---

# THM-3643 -- fixed-`107` real quadratic class number one

**PROVED + FINITE-EXACT + RIGOROUS-ARB; PENDING INDEPENDENT HOSTILE
AUDIT.**  Put

```text
D=1225041=3*408347,       K=Q(sqrt(D)).                 (1)
```

Then the ordinary ideal class number of `K` is

```text
h_K=1.                                                   (2)
```

The word *ordinary* is load-bearing.  This field has no unit of norm `-1`,
so its narrow class number is twice `(2)`.

## 1. Discriminant and complete unit audit

Trial division proves that `408347` is prime.  Hence `D` is squarefree; since
`D=1 mod 4`, it is the fundamental discriminant of `K`.  The exact continued
fraction is

```text
sqrt(D)=[1106; overline(a_1,...,a_1372)],
a_1372=2212,
SHA256(a_1,...,a_1372)=58a9d02a...423220bb.            (3)
```

The period length `1372` is even.  Standard Pell theory therefore gives both
of the following facts:

1. `x^2-Dy^2=-1` has no integral solution;
2. the convergent through `a_1371` is the least positive solution of
   `x^2-Dy^2=1`.

The companion reconstructs this solution from the continued-fraction
recurrences and verifies its Pell identity exactly.  Its coordinates have
`2361` and `2351` bits, respectively.

There is one ring-of-integers trap to exclude.  Since
`O_K=Z[(1+sqrt(D))/2]`, a genuinely half-integral unit would correspond to odd
`r,s` satisfying

```text
r^2-Ds^2=+4  or  -4.                                   (4)
```

But `D=1 mod 8`, and for odd `r,s` the left side of `(4)` is `0 mod 8`, while
the right side is `4 mod 8`.  Thus `(4)` is impossible.  The Pell unit

```text
epsilon=x+y sqrt(D)>1                                  (5)
```

is consequently a fundamental generator of the full unit group modulo
signs, and the regulator is `R_K=log(epsilon)`.

## 2. Finite character formula

Let `chi=chi_D` be the primitive real quadratic character modulo `D`.  Because
`D=1 mod 4`, it is even.  Set

```text
H=sum_{1<=a<D/2} chi(a) log(2 sin(pi a/D)).             (6)
```

The finite logarithmic formula for an even primitive character and the real
quadratic analytic class-number formula are

```text
L(1,chi)=-2H/sqrt(D),
L(1,chi)= 2 h_K R_K/sqrt(D).                            (7)
```

The factor `2` in the first identity comes from pairing `a` with `D-a`.
Equating `(7)` gives the especially sharp integer identity

```text
h_K=-H/R_K.                                             (8)
```

No asymptotic truncation occurs: `(6)` has exactly `(D-1)/2` terms.  The
character is evaluated as

```text
chi_D(a)=(a/3)(a/408347),                               (9)
```

with zero on nonunits.  Its half-range counts are exactly

```text
chi=-1: 204173,       chi=0: 204174,       chi=1: 204173. (10)
```

## 3. Rigorous interval gate

At `160`-bit Arb precision the complete finite sum and the exact Pell unit
give

```text
H   =[-1636.857453696630438860368892008879746470313
       +/- 8.07e-40],
R_K =[ 1636.8574536966304388603688920088797464703130480
       +/- 4.92e-44],
-H/R_K=[1.00000000000000000000000000000000000000000
         +/- 5.23e-43].                                (11)
```

The final ball is contained in `(0.999999999999,1.000000000001)` and contains
`1`.  Since `(8)` is a positive integer, `(2)` follows.

## 4. Reproduction and boundary

Run

```text
python3 04-computation/berggren_fixed107_real_quadratic_class_number_one_thm3643.py
python3 -O 04-computation/berggren_fixed107_real_quadratic_class_number_one_thm3643.py
```

Both modes print the stored transcript byte for byte.  The code uses no
floating-point value in the algebraic decisions; all logarithms and
trigonometric values in `(6)` are Arb balls with directed rounding.

The theorem removes an arithmetic ambiguity in a possible `3`-descent on the
Mordell curve of THM-3620.  It does **not** compute either `3`-Selmer group,
does not prove `rank E(Q)=2`, and does not enlarge the bounded slope atlas of
THM-3640.
