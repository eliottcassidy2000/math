---
id: THM-3580
title: "Berggren positive-cube slope-atlas completion through 101"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT.  In the primitive parity-correct
  slope universe 3<=n<=101, exactly four slopes admit the odd Pell data
  required by the positive two-cube compiler: (14,23), (26,29), (26,47),
  and (98,101).  THM-3547 already gives an infinite positive ray for each.
  Of its seven unresolved prime-screen survivors, three die by one-step
  p-adic descent and four die after a complete generalized-Pell class
  enumeration and finite fundamental-unit orbit calculation modulo 2n.
source: kps-s188, 2026-08-21
depends_on:
  - THM-3547-positive-two-cube-slope-atlas-through-101
related:
  - THM-3370-berggren-two-cube-biquadratic-norm-collision
  - THM-3375-berggren-positive-two-cube-pell-ray
  - THM-3376-positive-two-cube-slope-atlas-through-29
script: 04-computation/berggren_positive_cube_slope_completion_101_thm3580.py
output: 05-knowledge/results/berggren_positive_cube_slope_completion_101_thm3580.out
script_sha256: 09cd67a5bf0c2c0368e3e9ca4350d2cbb4c856ccc3bca4afa21254aaef26e184
output_sha256: c8dcdf02f7454a677e7bdc31a34c4909338b8d0b2084a646b44646d4d9af2f0b
hash_basis: LF-normalized bytes
---

# THM-3580 -- completion of the positive-cube slope atlas through 101

**PROVED + FINITE-EXACT + VERIFIED-EXACT.**  This closes the seven rows left
open by THM-3547.  It is a bounded slope classification, not an asymptotic
for sums of two cubes and not a classification for arbitrary denominator.

## 1. Inherited finite universe

THM-3547 exhausts the `528` reduced slopes

```text
3<=n<=101,       n odd,       m even,
n/2<m<n,         gcd(m,n)=1.                            (1)
```

Its exact screen by every prime at most `199` excludes `517` rows.  The
eleven survivors are

```text
(14,23), (26,29), (26,47), (38,47), (38,53),
(50,53), (50,71), (62,95), (74,95), (74,101),
(98,101).                                               (2)
```

The same theorem supplies explicit positive odd Pell rays for

```text
(14,23), (26,29), (26,47), (98,101).                   (3)
```

Only the other seven rows require work here.

## 2. The smaller norm equation and the exact parity address

The compiler equation is

```text
3W^2=n^2(4m^2-n^2)u^2+4(2m^2+2mn-n^2).                (4)
```

For every row in `(2)`, put

```text
K=(4m^2-n^2)/3,
C=4(2m^2+2mn-n^2)/3,
h=nu.                                                   (5)
```

Then `(4)` is the generalized Pell equation

```text
W^2-Kh^2=C.                                             (6)
```

All eleven `K` are positive odd nonsquares.  Since `n` is odd, the original
compiler requires precisely

```text
W odd,              h == n (mod 2n).                   (7)
```

Conversely `(7)` writes `h=nu` with `u` odd.  Thus the residue `(7)` is not
a heuristic filter: it is an iff translation of the parity/divisibility
sidecar lost when `(4)` is replaced by the smaller norm equation `(6)`.

## 3. Three one-step p-adic descents

Three rows fail before continued fractions are needed:

| `(m,n)` | `K` | `C` | prime `p` |
|---|---:|---:|---:|
| `(38,47)` | `1189` | `5668` | `13` |
| `(50,53)` | `2397` | `9988` | `11` |
| `(50,71)` | `1653` | `9412` | `13` |

In every row,

```text
v_p(C)=1,                 (K/p)=-1.                    (8)
```

If `(6)` had an integer solution, reduction modulo `p` would give
`W^2=Kh^2`.  Nonsquareness of `K` forces `p|W` and `p|h`; then the left side
of `(6)` is divisible by `p^2`, contradicting `(8)`.  This explains why the
old mod-`p` screen missed these rows: its residue pair `(W,h)=(0,0)` was a
false local survivor.  The missing coordinate was primitive `p`-adic order,
not a larger prime.

## 4. Complete generalized-Pell class reduction

For the remaining rows, use the standard Lagrange--Mollin--Matthews `P,Q,a`
continued-fraction reduction.  For completeness, the companion implements
the algorithm rather than treating a computer-algebra answer as an oracle:

1. for every square divisor `f^2|C`, put `c=C/f^2`;
2. enumerate every root `z^2=K (mod |c|)`;
3. run the exact `P,Q,a,A,B,G` recurrence for the full aperiodic-plus-periodic
   length of `(z+sqrt(K))/|c|`;
4. when `|Q|=1`, record the corresponding solution of norm `c`, using the
   negative-Pell correction exactly when required.

The resulting finite list contains one representative of every solution
class of `(6)`.  If

```text
epsilon=P+Q sqrt(K),       P^2-KQ^2=1,                 (9)
```

is the fundamental positive unit, all solutions arise from these
representatives by multiplication by `+/-epsilon^j` and conjugation.

Reduce this action modulo `2n`:

```text
(W,h) |-> (PW+KQh, QW+Ph) (mod 2n).                   (10)
```

Every orbit is finite and is enumerated until its exact starting pair
returns.  Signs and conjugation do not create a missed phase: because `n` is
odd, `-n=n (mod 2n)`, while the parity of `W` is sign-invariant.

The complete blocked ledger is

| slope | number of Pell classes | orbit periods mod `2n` | phases satisfying `(7)` |
|---|---:|---|---:|
| `(38,53)` | `9` | `2` in every class | `0` |
| `(62,95)` | `1` | `60` | `0` |
| `(74,95)` | `1` | `10` | `0` |
| `(74,101)` | `6` | `34` in every class | `0` |

Thus none of these four equations can have `h=nu` with odd `u,W`.
The stored transcript prints every class representative, fundamental unit,
orbit period, and target phase.  An independent library implementation of
the same LMM algorithm agrees exactly on all eleven rows; normal and
optimized runs are byte-identical.

For comparison, the admissible rows have nonempty target phases:

| slope | class count | unit-orbit period(s) | admissible class count |
|---|---:|---|---:|
| `(14,23)` | `3` | `11` | `2` |
| `(26,29)` | `6` | `10` | `2` |
| `(26,47)` | `9` | `46` | `6` |
| `(98,101)` | `6` | `34` | `4` |

This is a positive control on the parity address; THM-3547's explicit cone
certificates already show that each of these four slopes contains an
infinite positive orbit.

## 5. Classification

Combining the `517` inherited prime exclusions, the three descents of
Section 3, and the four complete orbit exclusions of Section 4 gives

```text
{primitive parity-correct slopes in (1) admitting (4)}
  = {(14,23),(26,29),(26,47),(98,101)}.                (11)
```

Every slope in `(11)` has infinitely many points with

```text
x>y>0,             x^3+y^3=(2r+1)^2+2,                (12)
```

by the inherited invariant-cone recurrences.  Every other slope in `(1)`
has no parity-correct integer point at all.

## 6. What changed and what remains

The old screen retained a scalar residue witness modulo each small prime.
The completion restores two pieces of information it had discarded:

```text
valuation depth at a prime     -> the three p-adic descents,
unit-orbit phase modulo 2n     -> the four parity-class exclusions.       (13)
```

This is the same representation-selection lesson that recurs elsewhere in
the repo: local support is not a globally owned orbit.  Here it is an exact
arithmetic statement, with no transfer to LRC ownership or Jacobian
boundary data.

The next unbounded question is whether the admissible slope set has an
arithmetic parametrization or a density.  Extending another finite prime
screen without carrying the complete Pell class and `2n` orbit from the
start would deliberately recreate the seven false survivors closed here.

## 7. Reproduction

```text
python 04-computation/berggren_positive_cube_slope_completion_101_thm3580.py
python -O 04-computation/berggren_positive_cube_slope_completion_101_thm3580.py
```

The companion uses integer arithmetic throughout.  Its optional SymPy call
is an independent equality control; the in-file LMM implementation and all
asserted conclusions do not depend on that import.
