---
id: THM-4026
title: "Sun 2-4-6-8 binomial conjecture is false"
status: >
  REFUTED + FINITE-EXACT + INDEPENDENTLY AUDITED. In the minimal-domain
  convention of OEIS A306477, a(896315812331399)=0. A regenerated 31-prime
  quadratic-residue certificate leaves no candidate among exactly
  2,755,643,831 feasible triples; an independent two-composite sieve checks
  3,796,359 survivors by restoring integer square root and also finds zero.
source: T. Adamczewski gist + root exact certificate + independent exact-sieve audit, 2026-08-24
depends_on: []
related:
  - THM-4027-sun-two-four-six-eight-universal-modular-solubility
  - THM-4028-sun-two-four-six-eight-average-order-criticality
  - MISTAKE-363
script: 04-computation/sun_2468_counterexample_modular_certificate_thm4026.cpp
output: 05-knowledge/results/sun_2468_counterexample_modular_certificate_thm4026.out
independent_audit_script: 04-computation/sun_2468_counterexample_independent_exact_sieve_thm4026.cpp
independent_audit_output: 05-knowledge/results/sun_2468_counterexample_independent_exact_sieve_thm4026.out
disjoint_bank_audit_script: 04-computation/sun_2468_counterexample_thm4026_independent.cpp
disjoint_bank_audit_output: 05-knowledge/results/sun_2468_counterexample_thm4026_independent.out
joint_period_local_audit_script: 04-computation/sun_two_four_six_eight_counterexample_independent_audit_thm4026.py
joint_period_local_audit_output: 05-knowledge/results/sun_two_four_six_eight_counterexample_independent_audit_thm4026.out
joint_period_local_audit_report: 05-knowledge/results/sun_two_four_six_eight_counterexample_period_local_audit_thm4026.md
eisenstein_row_sieve_script: 04-computation/sun_2468_eisenstein_inert_row_sieve_thm4026.py
eisenstein_row_sieve_output: 05-knowledge/results/sun_2468_eisenstein_inert_row_sieve_thm4026.out
---

# THM-4026 -- Sun's 2-4-6-8 conjecture is false

**REFUTED + FINITE-EXACT + INDEPENDENTLY AUDITED.** Put

\[
N=896315812331399.
\]

Then there are no integers

\[
w\ge2,\qquad x\ge3,\qquad y\ge5,\qquad z\ge7
\]

such that

\[
N=\binom w2+\binom x4+\binom y6+\binom z8.       \tag{1}
\]

Consequently the universal 2-4-6-8 conjecture posed by Zhi-Wei Sun is false.
This is a proof of nonexistence for one explicit integer, not a claim that
`N` is the least counterexample.

## 1. Source, conventions, and inheritance

The target and a NumPy verifier were published in the
[source gist](https://gist.github.com/tadamcz/0c578c8b2b3fb92fe8584bc0725187e3).
The original problem and the representation-count sequence are recorded at
[OEIS A306477](https://oeis.org/A306477) and in
[Sun's MathOverflow question](https://mathoverflow.net/questions/323541/positive-integers-written-as-binomw2-binomx4-binomy6-binomz8-with).
The source verifier was reproduced before the two methods below were built.

OEIS uses nonnegative shifted variables and writes

```text
C(a+2,2)+C(b+3,4)+C(c+5,6)+C(d+7,8).                  (2)
```

Thus `(2)` is exactly `(1)` with one canonical representative for each zero
atom. Sun's phrase “integers greater than one” also permits additional
indices below the lower order, such as `C(2,4)=0`; replacing any such zero by
`C(3,4)`, `C(5,6)`, or `C(7,8)` preserves existence. The conventions are
therefore existence-equivalent, although their raw multiplicity counts are
not the same.

The closest repo mechanism is the falling-factorial/Pascal coordinate of
THM-2412. The canonical hostile is MISTAKE-363: at primes dividing `r!`, one
must evaluate the integer `C(k,r)` and may not invert `r!` modulo the prime.
Both certificates below use exact integer binomials. THM-4001's factorial
response decoder is a corrected near miss: several symmetric responses can
recover the four summand values of an already-found representation, but the
single scalar `N` neither supplies those responses nor retains the four role
labels. The least-used useful sidecar is the complete role-labelled residue
distribution developed in THM-4027.

## 2. Exact finite universe

Every summand is nonnegative and the triangular summand is at least one.
Exact integer comparison with `N` gives

| term | top-index range | count | last admitted value | first excluded value |
|---|---:|---:|---:|---:|
| `C(w,2)` | `2..42,339,481` | `42,339,480` | `896315804504940` | `896315846844421` |
| `C(x,4)` | `3..12,112` | `12,110` | `896266258399820` | `896562324551740` |
| `C(y,6)` | `5..932` | `928` | `895693597430352` | `901490966993008` |
| `C(z,8)` | `7..281` | `275` | `871896500955975` | `897353333100675` |

The `y` and `x` ceilings shrink after the higher terms are chosen. There are
exactly `248,160` admissible `(z,y)` rows and

```text
2,755,643,831                                             (3)
```

feasible `(x,y,z)` triples satisfying

```text
C(x,4)+C(y,6)+C(z,8) <= N-1.                            (4)
```

This pins the universe. Applying the full `12,110`-entry `x` array to each of
the `248,160` admissible `(z,y)` rows would visit `3,005,217,600` lanes and
include `249,573,769` invalid ones. The still larger Cartesian product of the
three individual boxes has `3,090,472,000` lanes.

For a triple in `(4)`, set

```text
m=N-C(x,4)-C(y,6)-C(z,8) >= 1,
D=8m+1.                                                   (5)
```

The exact equivalence

```text
m=C(w,2) for some w>=2  <=>  D=(2w-1)^2                 (6)
```

turns the missing coordinate into a square test.

### An Eisenstein-norm sidecar

An incoming disjoint-prime-bank audit exposed a second exact coordinate. Put

```text
A=x^2-3x+1,       B=2w-1,       S=C(y,6)+C(z,8).
```

Direct expansion gives

```text
24*C(x,4)+1=A^2,              8*C(w,2)+1=B^2.          (6a)
```

Here `A>=1` and `B>=3` are odd. If `omega^2+omega+1=0`, then a
four-term representation is therefore equivalent to

\[
6(N-S)+1={A^2+3B^2\over4}
=N_{\mathbf Z[\omega]}\!\left({A+B\over2}+B\omega\right),       \tag{6b}
\]

together with the thin-image condition

```text
4A+5=(2x-3)^2.                                         (6c)
```

For a rational prime `ell=2 mod 3`, the Eisenstein norm form is anisotropic
modulo `ell`; hence an odd valuation of `ell` in `6(N-S)+1` eliminates that
entire `(y,z)` row. This is a useful row-level sieve, not a replacement for
the certificate: unrestricted norm representability loses positivity and
the square constraint `(6c)`. The incoming verifier records `(6b)` but its
decisive algorithm still uses the discriminant masks from `(6)`.

As a finite control, the routed companion applies the odd-valuation test for
all `87` inert primes through `1,000` to the complete `248,160` `(z,y)` rows.
It leaves `61,188`, the exact fraction `5099/20680`. This quantifies a cheap
row filter; those survivors are not asserted to satisfy the norm equation or
to extend to a four-term representation.

## 3. Primary regenerated modular certificate

The primary companion regenerates the quadratic residues for the 31 primes

```text
3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,
61,67,71,73,79,83,89,97,101,103,107,109,113,131,137.    (7)
```

For every `(z,y)` row and every still-valid `x`, it retains the bit for `x`
only when `D` from `(5)` is a square residue modulo the next prime. An integer
square must survive every mask, so every deletion is one-way sound. The
intersection is empty on every row:

```text
admissible (z,y) rows        248,160
feasible (x,y,z) triples   2,755,643,831
survive all 31 primes                  0.                (8)
```

Thus `(8)` already proves nonexistence; no floating-point square root is used
on the negative result. The program computes all binomials with guarded exact
division in an unsigned 128-bit intermediate and all residues from those
integers. Its only square-root routine is a restoring integer algorithm used
for positive controls that survive `(7)`.

## 4. Independent exhaustive audit

The independent program has a different organization and rejection carrier.
It uses the two coprime odd composite moduli

```text
M1=3*5*7*11*13=15015,
M2=17*19*23*29=215441.                                  (9)
```

For odd `M`, triangular residues have period `M`, since

```text
C(k+M,2)-C(k,2)=M(2k+M-1)/2.                           (10)
```

The audit builds the full triangular residue tables modulo `(9)`, scans all
triples in `(3)`, and applies a restoring binary integer square root to every
survivor. Its exact stage counts are

```text
feasible triples             2,755,643,831
survive M1                      78,900,865
survive both M1 and M2           3,796,359
integer-square hits                      0.              (11)
```

For every one of the `3,796,359` final candidates it asserts

```text
s^2 <= D < (s+1)^2
```

and tests `s^2=D`. Normal, optimized, single-threaded, serial, and a deep
residue-audit build agree after timing/thread metadata are removed. This
method shares the discriminant identity `(6)`, which is unavoidable, but not
the prime-by-prime bitset certificate or its control flow.

Two later incoming replays give further redundancy. A dependency-free Python
bitset route scans the full `3,090,472,000`-lane product of the individual
top-index boxes, leaves `324` mask survivors through prime `89`, rejects `31`
by residual admissibility and the other `293` by `isqrt`, and leaves zero
survivors after prime `137`. A separate C++ implementation uses
two disjoint 15-prime banks; each bank followed by exact square tests finds
zero representations, and their combined 29-prime cover leaves zero modular
survivors. These are independent implementations of the same necessary-square
mechanism, not additional conceptual proofs.

## 5. Positive and hostile controls

The primary certificate also exhausts both neighbours in the same convention:

```text
a(N-1)=89,      a(N)=0,      a(N+1)=67.                (12)
```

In particular,

```text
N-1 = C(33663667,2)+C(9433,4)+C(16,6)+C(9,8),
N+1 = C(40920205,2)+C(6138,4)+C(22,6)+C(13,8).         (13)
```

The independent audit finds the same witnesses and also checks the boundary

```text
0: no representation,       1=C(2,2)+C(3,4)+C(5,6)+C(7,8).  (14)
```

This catches the otherwise easy error of allowing a zero triangular term.
Only the two adjacent controls are promoted here; no wider isolation or
minimality claim is made.

## 6. What the certificate preserves and loses

The map `(x,y,z)->D` preserves exact equality to a triangular term through
`(6)`. The modular masks preserve a necessary square predicate and the exact
bounded index address. They destroy the square root itself and, before the
finite box is imposed, all magnitude information. THM-4027 proves that the
full four-term sum is soluble modulo every modulus, so no fixed congruence
class can replace the bounded certificate. The missing sidecar is precisely
the archimedean box `(4)`.

There is also a lawful Pascal-carry coordinate. The greedy rank-eight
combinadic expansion is exactly

```text
N=C(281,8)+C(279,7)+C(234,6)+C(212,5)
  +C(188,4)+C(136,3)+C(43,2)+C(15,1).                 (15)
```

Thus `N` lies on a full descending Pascal flag even though it misses Sun's
role-labelled four-atom family supported at ranks `2,4,6,8`. The identity
`C(t,k)=C(t-1,k)+C(t-1,k-1)` gives a split relation and, in reverse, a carry,
but it does not preserve atom count or the one-atom-per-role constraint.
Whether a confluent bounded carry normal form detects this hole is an **OPEN**
structural reformulation, not part of the certificate.

The source gist's NumPy `sqrt` plus integer square-back is sound at this fixed
scale (`D<2^53`), but it is not the load-bearing method here. Both canonical
paths use integer-only decisive arithmetic and obey MISTAKE-363.

## 7. Scope and reproduction

This theorem refutes the universal conjecture. It does **not** show that `N`
is the first counterexample: the previously published check only reached
`2*10^12`, far below `N`, and no complete intervening scan is imported here.
The formal-proof work described in the source transcripts is useful
provenance but is not a dependency because that large Lean artifact is not
distributed in the gist.

From the repository root, compile and run the primary certificate with

```text
g++ -std=c++20 -O3 -Wall -Wextra \
  04-computation/sun_2468_counterexample_modular_certificate_thm4026.cpp \
  -o sun2468-thm4026
./sun2468-thm4026
```

and the independent audit with

```text
g++ -std=c++20 -O2 -fno-tree-vectorize -Wall -Wextra -Wconversion -Wshadow \
  04-computation/sun_2468_counterexample_independent_exact_sieve_thm4026.cpp \
  -o sun2468-thm4026-audit
./sun2468-thm4026-audit
```

Both return zero only after their internal controls pass. **QED.**
