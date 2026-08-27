---
id: THM-4237
title: "Multiplier-six binary adjacency and prime-power factorial closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. In THM-3474's exact
  prime-power reset compiler with multiplier a=6, the positive candidate
  degree intersection is empty exactly when the binary expansion of p^k has
  adjacent one-bits.  Hence every prime p=3 mod 4, p>=7, and every odd k
  gives a new coprime exact-quadratic factorial-moment window at
  exponents 6p^k-1,6p^k,6p^k+1.  A nonempty candidate intersection is not
  asserted to construct a factor.
source: root/cross-frontier-bridge/2026-08-26
depends_on:
  - THM-3474-factorial-binary-submask-polygon-and-prime-power-reset-families
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
related:
  - THM-3260-ternary-cantor-bessel-cartier-residual-exclusion
  - THM-3475-factorial-all-divisor-digit-polygon-and-pair-ledger-compiler
scripts:
  - 04-computation/factorial_multiplier_six_binary_adjacency_thm4237.py
  - 04-computation/factorial_multiplier_six_binary_adjacency_independent_audit_thm4237.py
outputs:
  - 05-knowledge/results/factorial_multiplier_six_binary_adjacency_thm4237.out
  - 05-knowledge/results/factorial_multiplier_six_binary_adjacency_independent_audit_thm4237.out
script_sha256:
  - c92715b3a956fafe6575a26f504f7103b110aec8a3e067323da3411150e6e906
  - 1493abcbc754a6f60cccedec8301397f6646e22463dcaf2196aadfcd9835f6cb
output_sha256:
  - 126a0d484c27175dc33b15db6e45f1fb4864f1f11a6446cd84001d52db93b356
  - 4d3a849e32e5a1a37cdbafcd3016e74fd87afa04ff9caa22fa2cfe11fb31f9e3
hash_basis: raw LF bytes
audit: >
  Adversarial audit PASS. For a=6, THM-3474 gives reset range T=3 at p=7
  and T=5 thereafter; carry-free splitting eliminates odd multipliers and
  identifies the even pair with adjacent binary support. Empty degree
  barcode and zero coordinate capacity give rational coprimality, while the
  explicit v=(N+1)C/A scaling transfers THM-3124's resonance root to the
  same polynomial pair. Both independent normal/optimized replays byte-match
  their frozen outputs and hashes; surviving-address hostiles construct no
  factor.
---

# THM-4237 -- multiplier-six binary adjacency and prime-power factorial closure

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Statement

Let `p>=7` be prime, `k>=1`, and put

```text
H=p^k,                         N=6H.                         (1)
```

Say that `Adj(H)` holds when the binary expansion of `H` has two adjacent
one-bits, equivalently

```text
H & (H>>1) != 0.                                             (2)
```

THM-3474 intersects its exact `p`-adic reset barcode with the complete
binary-submask barcode of the other resonant row.  For multiplier six, write
the surviving multiplier set as

```text
C(p,k)={t:1<=t<=T and tH is a binary submask of 6H},
T=min(5,(p-1)/2).                                           (3)
```

Then the exact classification is

```text
C(p,k)=empty                         if Adj(H),
C(7,k)={2}                           if not Adj(H),
C(p,k)={2,4}                         if p>=11 and not Adj(H). (4)
```

Consequently, whenever `Adj(p^k)` holds,

```text
gcd_Q(A_(6p^k-1)^(6p^k+1), A_(6p^k)^(6p^k+1))=1.           (5)
```

Here, as in THM-3474,

```text
L(x^r)=r!,                  A_n^d(v)=L((d-x+v x^2)^n).      (6)
```

In particular, if

```text
p=3 (mod 4),                       k odd,                   (7)
```

then `p^k ≡ 3 (mod 4)`, its final two binary digits are `11`, and `(5)` holds.
For every exact quadratic

```text
q(x)=A+Bx+Cx^2,                    ABC!=0,                  (8)
```

at least one of

```text
L(q^(6p^k-1)),        L(q^(6p^k)),        L(q^(6p^k+1))     (9)
```

is therefore nonzero under `(7)`.

The equivalence `(4)` classifies when this particular exact cross-place
compiler closes the multiplier-six row.  Its nonempty side is only a degree
address and does not construct a common factor.

## 2. Binary submasks are carry-free splittings

For nonnegative integers `x<=y`,

```text
x is a binary submask of y    iff    x & (y-x)=0.           (10)
```

Indeed, submask containment makes `y-x` the complementary bit set of `x`
inside `y`.  Conversely, disjoint supports make the addition
`x+(y-x)=y` carry-free, so every bit of `x` occurs in `y`.

Apply `(10)` with

```text
x=tH,                         y=6H.                         (11)
```

The submask condition in `(3)` is therefore equivalent to

```text
(tH) & ((6-t)H)=0.                                           (12)
```

This changes a static bit-containment question into an exact two-summand
carry question while retaining the candidate degree `tH`.

## 3. The multiplier-six classification

Because `H` is odd, if `t` is odd then both `tH` and `(6-t)H` are odd.  Their
binary supports both contain position zero, so `(12)` fails.  Thus no odd
multiplier survives.

The reset range is

```text
T=3                         for p=7,
T=5                         for p>=11.                     (13)
```

The only possible multiplier is consequently `t=2` when `p=7`, and the only
possible multipliers are `t=2,4` when `p>=11`.  Both test the same unordered
pair

```text
{tH,(6-t)H}={2H,4H}.                                      (14)
```

Multiplication by two and four shifts the binary support of `H` by one and
two positions.  Hence

```text
(2H)&(4H)!=0
  iff (supp_2(H)+1) intersects (supp_2(H)+2)
  iff H has adjacent one-bits.                              (15)
```

If `Adj(H)` holds, `(12)` fails for every candidate in `(13)`.  If it does
not hold, `2H` and `4H` are disjoint, so `t=2` survives, as does `t=4` when
it lies in the reset range.  This proves `(4)`.

## 4. Transfer through the exact reset compiler

For `a=6`, THM-3474 applies because

```text
a is even,                  2<=a<p,                         (16)
```

for every prime `p>=7`.  It proves that any positive common-factor degree of
the two rows in `(5)` must have the form `tH`, `1<=t<=T`, and must also be a
binary submask of `aH=6H`.  It separately proves that the coordinate-root
capacity is zero.  Thus `(4)` and `Adj(H)` leave no positive degree address,
which proves the rational coprimality `(5)`.

THM-3124's symbolic resonance reduction says that three consecutive zero
moments at `(9)` first force `B/A=-1/d`, where `d=N+1`. Put

```text
u=C/A,                 v=du,
q=(A/d)(d-x+v x^2),    L(q^n)=(A/d)^n A_n^d(v).          (16a)
```

Thus the same `v` would be a common complex root of the rational polynomial
pair in `(5)`. A common complex root of two rational polynomials has a
positive-degree minimal polynomial common to both, contradicting their
rational gcd `1`. This proves `(9)`.

Finally, `(7)` gives

```text
p^k ≡ 3^k ≡ 3 (mod 4),                                     (17)
```

because `k` is odd.  The two low binary bits are therefore both one, so
`Adj(H)` holds without inspecting any higher bit.  This proves the announced
all-height-in-odd-exponent family.

## 5. Equality and failure boundary

The binary criterion is sharp for the **compiler**.

- At `(p,k)=(23,2)`,

  ```text
  H=529=1000010001_2.                                      (18)
  ```

  There are no adjacent one-bits, and `(4)` leaves `t=2,4`, hence candidate
  degrees `1058,2116`.  This shows that even exponents cannot simply be added
  to `(7)`.  It does not exhibit a factor or a bad moment window.
- At `(p,k)=(17,1)`, `H=10001_2` gives the same two surviving multipliers.
  Thus no statement for all primes `p=1 (mod 4)` follows.
- Conversely, `(p,k)=(13,1)` has `H=1101_2`, and `(p,k)=(7,2)` has
  `H=110001_2`; both close by adjacency.  Therefore `(7)` is sufficient, not
  necessary.

The mechanism preserves the complete THM-3474 degree address and loses no
candidate inside that compiler.  What it does not preserve is factor
existence: a surviving local degree address is only necessary.

## 6. Verification and nonclaims

The primary companion checks `(4)` for every odd
`1<=H<2^16`, separately in reset ranges `T=3,5`, and checks every prime
`7<=p<500`, `1<=k<=25`.  Its declared universes contain `65,536` odd-height
cells and `2,300` prime-power cells.  It freezes the four boundary controls
above and the exact hostile candidate degrees.

The independent companion imports no primary code.  It represents integers
as explicit sets of occupied binary positions, exhausts every odd
`1<=H<2^14`, and independently rebuilds `1,273` prime-power intersections
for `p<350`, `k<=19`.  It separately checks the odd low-bit overlap and the
even shift-adjacency equivalence.

Reproduce with

```bash
python3 -B 04-computation/factorial_multiplier_six_binary_adjacency_thm4237.py
python3 -B -O 04-computation/factorial_multiplier_six_binary_adjacency_thm4237.py
python3 -B 04-computation/factorial_multiplier_six_binary_adjacency_independent_audit_thm4237.py
python3 -B -O 04-computation/factorial_multiplier_six_binary_adjacency_independent_audit_thm4237.py
```

The computation audits the elementary bit lemma and its specialization; the
unbounded claims follow symbolically from `(10)`--`(17)` and THM-3474.

No claim is made here for:

- a common factor when `(4)` is nonempty;
- arbitrary multipliers `a`, arbitrary three-slot supports, or translated
  supports;
- even exponents outside the exact adjacency condition;
- `SFC(1)`, `SFC(3)`, or the multivariable Factorial Conjecture.
