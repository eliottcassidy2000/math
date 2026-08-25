---
id: THM-4040
title: "Sun 2-4-6-8 target CRT-class singleton census"
status: >
  FINITE-EXACT + VERIFIED-EXACT. Among all 943,194,644 positive targets at
  most 1,001,999,999,999,999 in the discovery class 459490 modulo 1062347,
  the canonical Sun 2-4-6-8 representation count vanishes exactly once, at
  896315812331399. Thus the certified hole is least and unique in this one
  fixed class through the displayed height. This is not global leastness and
  does not classify holes in other classes or above the finite bound.
source: codex exact class census + T. Adamczewski transcript cross-check, 2026-08-24
depends_on: []
related:
  - THM-4026-sun-two-four-six-eight-binomial-counterexample
  - THM-4027-sun-two-four-six-eight-universal-modular-solubility
  - THM-4036-sun-2468-energy-and-support-exponent
script: 04-computation/sun_2468_exact_residue_class_block_scan.py
output: 05-knowledge/results/sun_2468_target_crt_class_singleton_census_thm4040.out
embedded_cpp_sha256: cedc808e7665d6ae98b687c8b0c152a9ebb3680afa4785ef6066e58f27fbba92
independent_audit_script: 04-computation/sun_2468_target_crt_class_independent_audit_thm4040.py
independent_audit_output: 05-knowledge/results/sun_2468_target_crt_class_independent_audit_thm4040.out
independent_audit_report: 05-knowledge/results/sun_2468_target_crt_class_independent_audit_thm4040.md
source_gist: https://gist.github.com/tadamcz/0c578c8b2b3fb92fe8584bc0725187e3
---

# THM-4040 -- the target CRT class has one hole through 1.002 quadrillion

**FINITE-EXACT + VERIFIED-EXACT.** Let `a(n)` be the canonical OEIS-domain
representation count

\[
a(n)=\#\left\{(w,x,y,z):
n=\binom w2+\binom x4+\binom y6+\binom z8,
\ w\ge2,x\ge3,y\ge5,z\ge7\right\}.                 \tag{1}
\]

Put

```text
M=11*13*17*19*23=1,062,347,
R=459,490,
H=1,001,999,999,999,999,
N=896,315,812,331,399.                                (2)
```

Then

\[
\boxed{1\le n\le H,\ n\equiv R\pmod M,\ a(n)=0
       \quad\Longleftrightarrow\quad n=N.}           \tag{3}
\]

There are exactly

```text
943,194,644                                             (4)
```

targets in the finite progression in `(3)`. In particular, `N` is the
smallest counterexample in this progression, and it is the only one there
below `1.002*10^15`.

## 1. What this class means

The target residues are

```text
R mod (11,13,17,19,23) = (9,5,14,13,19).              (5)
```

These are the low-density local coordinates used by the discovery search;
CRT combines them into `(2)`. The class is a search prior, not a necessary
condition on a hole. THM-4027 proves that every residue modulo every modulus
is represented, including every residue modulo `M`. Thus `(3)` is a finite
height statement about one progression and cannot be promoted to a modular
classification of all counterexamples.

The public discovery transcript scanned this class in large blocks. The
companion here was independently written from the representation equation,
uses repaired exact integer decisions, and replays the full class from one.
The agreement with the historical block totals is a cross-check rather than
a logical dependency of `(3)`.

## 2. Exact block engine

Fix an inclusive target interval `[L,U]`. For the three higher roles, the
engine builds the strictly increasing exact lists

```text
C(x,4), x>=3;    C(y,6), y>=5;    C(z,8), z>=7,       (6)
```

stopping before `U`. Equality to `U` cannot participate because the
triangular role is at least one. It then visits every higher triple with

\[
S=\binom x4+\binom y6+\binom z8<U.                   \tag{7}
\]

The remaining role is triangular. With `i=w-1>=1`,

\[
\binom w2=T_i={i(i+1)\over2}.                         \tag{8}
\]

For a fixed triple `(7)`, every admissible representation in the block has

```text
max(1,L-S) <= T_i <= U-S,
T_i == R-S (mod M).                                   (9)
```

The triangular residues have period `2M`, since

\[
T_{i+2M}-T_i=M(2i+2M+1).                             \tag{10}
\]

The engine enumerates all roots `i mod 2M` once, then intersects their
arithmetic progressions with the exact index interval in `(9)`. Its initial
long-double square-root seeds are not decisive: exact unsigned 128-bit
comparisons increment or decrement each seed until it is the true triangular
ceiling or floor. Every resulting mark is checked to lie in the selected
target class and addressed by `(n-base)/M`.

This process is exhaustive and duplicate-free: every representation fixes
one higher triple and one positive triangular index, while every enumerated
pair reconstructs one representation. Binomial products and triangular
comparisons use unsigned 128-bit intermediates; throughout the certified
height all stored coefficients fit in 64 bits. The target counters saturate
at `65,535`, but saturation can only change a large positive multiplicity.
Therefore the zero predicate, every count below the saturation threshold,
and in particular every block minimum displayed below are exact.

## 3. Contiguous finite universe and output

The nine half-open input blocks in the script form the disjoint partition

```text
[1, 1,002,000,000,000,000).                           (11)
```

Their exact summaries are:

| block | inclusive interval | targets | feasible higher triples | representation marks | minimum `a(n)` | first minimizer | zeros |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | `[1,1,999,999,999,999]` | `1,882,624` | `99,211,900` | `51,852,255` | 3 | `636,132,780,743` | 0 |
| 2 | `[2,000,000,000,000,11,999,999,999,999]` | `9,413,120` | `263,480,161` | `283,941,711` | 2 | `9,480,152,433,497` | 0 |
| 3 | `[12,000,000,000,000,101,999,999,999,999]` | `84,718,082` | `845,239,561` | `2,789,161,302` | 2 | `99,934,198,737,404` | 0 |
| 4 | `[102,000,000,000,000,201,999,999,999,999]` | `94,131,202` | `1,225,867,610` | `3,272,771,677` | 3 | `122,533,388,547,887` | 0 |
| 5 | `[202,000,000,000,000,301,999,999,999,999]` | `94,131,202` | `1,525,588,515` | `3,321,577,932` | 4 | `202,687,664,107,388` | 0 |
| 6 | `[302,000,000,000,000,401,999,999,999,999]` | `94,131,202` | `1,782,451,689` | `3,390,881,393` | 3 | `350,011,776,754,814` | 0 |
| 7 | `[402,000,000,000,000,601,999,999,999,999]` | `188,262,404` | `2,220,090,999` | `6,842,995,836` | 4 | `406,999,014,660,698` | 0 |
| 8 | `[602,000,000,000,000,801,999,999,999,999]` | `188,262,404` | `2,594,651,419` | `6,944,628,451` | 4 | `677,862,454,957,862` | 0 |
| 9 | `[802,000,000,000,000,1,001,999,999,999,999]` | `188,262,404` | `2,927,893,995` | `6,984,950,518` | 0 | `896,315,812,331,399` | 1 |

The target counts sum to `(4)`. The only printed zero lane is

```text
LOW target=896315812331399 count=0.                    (12)
```

This proves `(3)`. As an independent target-level control, the block engine
visits exactly `2,755,643,831` feasible higher triples at `N`, the same
universe certified by THM-4026's disjoint exact sieves.

## 4. Positive and hostile controls

Before the large census, the wrapper performs a separate direct Cartesian
enumeration in Python on

```text
[1,10000], target class 20 mod 33.                    (13)
```

It compares every one of the `303` class counts against the C++ engine. All
counts agree; the minimum is exactly two at `9590`, and there are no zeros.
This control simultaneously exercises the canonical zero atoms, the positive
triangular boundary, multiple triangular roots, interval addresses, and
nonzero low multiplicities.

At the opposite scale, an additional one-off scan on

```text
[92,000,000,000,000,101,999,999,999,999]
```

exactly reproduces the historical ten-trillion control

```text
triples=845,239,561, marks=315,855,003, targets=9,413,120,
minimum=2 at 99,934,198,737,404, zeros=0,              (14)
```

while the aggregate block 3 and blocks 4--9 reproduce the transcript's
corresponding block totals. THM-4026's
independent prime-bank and exact-square methods verify the unique zero in
block 9 and its represented neighbours. These controls have different
purposes: the direct small check tests every multiplicity, the historical
replay tests block aggregation, and THM-4026 tests the one negative target by
different sieves.

## 5. Global leastness and classification firewall

Equation `(3)` does **not** show that `N` is the smallest counterexample over
all positive integers. The undischarged global interval remains

```text
2,000,000,000,001 <= n <= 896,315,812,331,398,         (15)
```

containing `894,315,812,331,398` integers after the cited public prefix
search. This theorem examines only about one target in every `1,062,347`.
Nor does it classify holes above `H` or in another residue class modulo `M`.
Allowing redundant zero indices changes positive multiplicities, but not the
existence or zero set, because every zero atom has the canonical representative
used in `(1)`.

The valid conclusions are precisely:

1. **FINITE-EXACT:** `N` is the least hole in the fixed class `(5)`.
2. **FINITE-EXACT:** `N` is the only hole in that class through `H`.
3. **OPEN:** global leastness, all other holes, finiteness or infinitude of
   the hole set, and its density.

THM-4036 proves represented support has logarithmic exponent one in every
fixed arithmetic progression, while THM-4027 rules out an exceptional tail
in any progression. Those asymptotic constraints are compatible with this
finite singleton census; neither upgrades it to a global classification.

## 6. Replay

On a 64-bit machine with a C++20 compiler, run

```text
python -B 04-computation/sun_2468_exact_residue_class_block_scan.py \
  --target-class-census --jobs 4 \
  --build-dir .scratch/sun2468-target-class-census
```

The wrapper compiles the embedded exact engine once, runs the direct Python
control, replays the nine blocks, checks every summary and the singleton hole
list, and returns zero only after printing `RESULT=PASS`. **QED.**
