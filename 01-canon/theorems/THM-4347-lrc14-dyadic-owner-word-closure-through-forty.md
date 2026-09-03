---
id: THM-4347
title: "LRC(14) dyadic owner-word closure for every eleven-body through height forty"
status: >
  PROVED RELATIVE TO THM-2066/4070/4075 + FINITE-EXACT + DUAL OWNER-WORD/
  LITERAL-LIFT PACKET AUDIT. Every eleven-element positive body H with
  max(H)<=40 closes 2H union {a,b} for every two distinct positive odd tails
  a,b. THM-4070 closes non-divisor-complete bodies, THM-4075 closes the
  divisor-complete bodies through height 30, and a new exact census closes
  every divisor-complete exact-maximum layer 31,...,40. The new census scans
  2,257,174,140 bodies, retains 73,285,588 divisor-complete bodies including
  22,765 nonprimitive ones, and finds an empty complementary-owner relation
  at some clock 15<=N<=43 for every retained body. Tail magnitudes are
  unrestricted. This is a bounded half-body theorem, not arbitrary entry,
  not a closure of other two-adic root types, and not LRC(14).
source: root + arbitrary_entry / LRC14 continuation session, 2026-09-02
depends_on:
  - THM-2066-dyadic-seam-owner-word-crt-atlas
  - THM-4070-lrc14-d2-mod14-two-bank-affine-ray-firewall
  - THM-4075-lrc14-divisor-complete-dyadic-owner-word-closure-through-30
related:
  - THM-615-folding-identity-and-AP-even-part-confinement
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4326-lrc14-rank-two-wall-graph-complete-typed-universe-closure
  - THM-4330-lrc14-affine-two-adic-root-types-and-anchored-pool-entry-sieve
  - THM-4332-lrc14-fixed-pool-single-constraint-implication-rigidity
  - THM-4338-lrc14-all-pair-rank-four-cubic-majorant-and-uniform-one-appender
script: 04-computation/lrc14_arbitrary_entry_owner_word_through40_probe_arbitrary_entry_20260902.cpp
output: 05-knowledge/results/lrc14_arbitrary_entry_owner_word_through40_probe_arbitrary_entry_20260902.out
reverse_output: 05-knowledge/results/lrc14_arbitrary_entry_owner_word_through40_reverse_probe_arbitrary_entry_20260902.out
separation_script: 04-computation/lrc14_arbitrary_entry_gate_separation_probe_arbitrary_entry_20260902.py
separation_output: 05-knowledge/results/lrc14_arbitrary_entry_gate_separation_probe_arbitrary_entry_20260902.out
script_sha256: e0298461737adc15d2d84d55f86298ad4419fe5c9e74ff0fb97427d8371ebade
output_sha256: 9b1f34c9caa52eceb9cd2db67c3cb12b5a00ca5befd7764fa2d2a6e7d746d401
reverse_output_sha256: e871318601799b0028de6c80eb6ceeb2d4c19c6ee910362a2d2a969580c4924f
separation_script_sha256: 277bb7e517cd32f90c775f50c7143faf88c531e8957fdf632b04ed2a0af4e261
separation_output_sha256: 58a22fd47489981076e95b86f2ccbbe5a72599f000e758c345e053068823cea8
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT WITH SCOPED DUAL SEMANTIC AUDIT. Strict C++20 -O2 and -O3
  ascending runs are byte-identical. A descending-clock replay changes the
  selected-clock histogram but reproduces every raw/filter/primitive count
  and zero survivors. On every newly encountered packet the program compares
  the THM-2066 nearest-integer owner-word test with a separate literal test of
  every unordered odd residue pair with repetition on both physical lifts;
  all 4,168 ascending and 18,128 descending packet audits agree. The small
  exact-Fraction companion passes normal, optimized, and hash-seeded replay.
---

# THM-4347 -- dyadic owner-word closure through forty

**PROVED RELATIVE TO THM-2066/4070/4075 + FINITE-EXACT + SCOPED DUAL
SEMANTIC AUDIT. LRC(14) REMAINS OPEN.**

## 1. Statement and inherited partition

Let `H` be any eleven-element set of distinct positive integers with

```text
max(H)<=40,                                             (1)
```

and let `a,b` be any two distinct positive odd integers.  Then

```text
S(H;a,b)=2H union {a,b}                                 (2)
```

has a time `x in R/Z` such that

```text
||vx||>=1/14                   for every v in S(H;a,b). (3)
```

There is no upper bound on `a` or `b`, and no primitivity hypothesis on `H`.
The height in `(1)` is the maximum of the normalized eleven-element
**half-body**, not the maximum speed of the full row.

The proof uses the exact exhaustive partition

```text
some d in {2,...,14} divides no h in H;                 (4)

or

for every d in {2,...,14}, some h in H is divisible by d. (5)
```

THM-4070 proves `(3)` in branch `(4)` by explicit small-denominator two-lift
banks.  In branch `(5)`, call `H` divisor-complete.  THM-4075 proves `(3)`
for every divisor-complete body with `max(H)<=30`.  Sections 2--3 below close
the ten new exact-maximum layers `31,...,40`.  This proves the statement.

The closest mechanism is THM-2066's labelled owner word.  The canonical
combined-gate hostile is

```text
H_*={1,2,3,4,8,9,10,11,12,13,14},                     (6)
```

which is divisor-complete but fails both THM-4330 scalar entry tests at tail
ratio `(1,9)`.  The corrected near misses are change of the distinguished
reference, deletion without restoration of the deleted constraint, and
MISTAKE-238's attempt to transport an empty full-row safe set through a map
defined only on the nonempty quotient-core safe set.  The load-bearing
sidecars here are the canonical core-residue labels and the odd tails' two
physical lift masks.

## 2. Faithful owner-word certificate

For an integer clock `N`, write

```text
|z|_N=min_(k in Z)|z-kN|
```

and form the labelled core-safe packet

```text
A_N(H)={0<=r<N:14|hr|_N>=N for every h in H}.           (7)
```

At `r in A_N(H)`, the two full-row candidate times are

```text
x_(r,j)=(r+jN)/(2N),                  j=0,1.            (8)
```

For every odd residue `z mod 2N`, THM-2066 retains whether `z` is dangerous
on a lift in `(8)` and, whenever it is dangerous over every packet point,
the binary nearest-integer parity word identifying the killed lift.  Let
`R_N(H)` be the relation of unordered odd residue pairs, with repetition,
whose two literal masks cover both lifts over every `r in A_N(H)`.  This is
equivalently THM-2066's complementary eligible-owner relation.

THM-2066 proves the exact implication

```text
A_N(H) nonempty and R_N(H)=empty
  ==> S(H;a,b) satisfies (3) for every distinct positive odd a,b. (9)
```

Repeated residue classes are retained in `R_N(H)` because distinct physical
tails can be congruent modulo `2N`.  No tail-height truncation occurs: the
literal danger masks at the rational times `(8)` depend exactly on the tails
modulo `2N`.

The new computation proves that every divisor-complete body in the ten new
layers has a certificate `(9)` at some

```text
15<=N<=43.                                              (10)
```

## 3. Complete new finite universe

For an exact maximum `M`, the program enumerates every body uniquely as

```text
H=H_0 union {M},       H_0 in binom({1,...,M-1},10),   (11)
```

then applies exactly the thirteen divisibility predicates `(5)`.  It does
not filter on gcd.  The exact results are:

| `M` | searched | divisor-complete | primitive | nonprimitive | uncertified |
|---:|---:|---:|---:|---:|---:|
| 31 | 30,045,015 | 522,173 | 522,173 | 0 | 0 |
| 32 | 44,352,165 | 893,833 | 893,307 | 526 | 0 |
| 33 | 64,512,240 | 2,260,900 | 2,260,900 | 0 | 0 |
| 34 | 92,561,040 | 1,366,695 | 1,365,734 | 961 | 0 |
| 35 | 131,128,140 | 1,623,006 | 1,623,006 | 0 | 0 |
| 36 | 183,579,396 | 9,005,059 | 9,000,516 | 4,543 | 0 |
| 37 | 254,186,856 | 3,959,314 | 3,959,314 | 0 | 0 |
| 38 | 348,330,136 | 4,644,109 | 4,638,715 | 5,394 | 0 |
| 39 | 472,733,756 | 17,211,708 | 17,211,708 | 0 | 0 |
| 40 | 635,745,396 | 31,798,791 | 31,787,450 | 11,341 | 0 |
| **total** | **2,257,174,140** | **73,285,588** | **73,262,823** | **22,765** | **0** |

The three exact-max-40 bodies first certified at the last clock `43` in the
ascending bank are

```text
{11,14,19,23,25,26,29,34,36,37,40},
{12,18,19,22,25,26,28,31,34,37,40},
{19,22,25,26,28,29,31,34,36,37,40}.                  (12)
```

Thus `43` is a genuine endpoint for this ordered bank.  It is not claimed to
be a globally minimal clock for any body, nor is `(10)` claimed beyond
height `40`.

## 4. Dual semantic audit

The primary predicate constructs THM-2066's eligible tails and binary owner
words, then tests whether a word and its packet complement are both realized.
For every newly encountered packet, a separate predicate avoids nearest
integers entirely.  It enumerates all unordered odd residues with repetition
modulo `2N`, evaluates

```text
14 |z(r+jN)|_(2N) < 2N,             j=0,1,             (13)
```

and tests literal coverage of both lifts.  The two predicates agree on every
one of the `4,168` packet/clock states encountered by the ascending census.

Strict C++20 `-O2` and `-O3` builds give byte-identical ascending transcripts.
A descending-clock run reaches `18,128` packet states, has different
selected-clock histograms as expected, and reproduces every searched,
divisor-complete, primitive, nonprimitive, and zero-survivor total.  Layer
accounting independently checks

```text
searched=binom(M-1,10),
divisor_complete=primitive+nonprimitive,
sum(selected-clock histogram)+uncertified=divisor_complete. (14)
```

## 5. Separation from scalar entry gates

The exact-Fraction companion verifies that `(6)` has only six fixed-pool
hits under every lawful positive-rational body refactorization, with the
maximum attained at scale `10`, while

```text
mu(G_(H_*))=20411/360360 < 4/63=mu(C_(1,9)).           (15)
```

Thus the projective fixed-pool and pair-adaptive Haar gates both fail for
the row `2H_* union {1,9}`.  Nevertheless

```text
A_23(H_*)={4,9,10,13,14,19},                           (16)
```

and its literal owner relation is empty for **every** odd residue pair.  For
the displayed tails, `x=2/23` has clearance `2/23`.  This proves that the
owner-word operation is a strict faithful bypass, not merely a repackaging
of either scalar gate.

Two additional controls pin the discarded coordinates:

1. Deleting speed `11` from `{1,...,11}` exposes the safe time `1/11`, but
   the restored speed has clearance zero there.  Deletion reverses the useful
   safe-set containment unless the deleted constraint is checked at the same
   phase.
2. In `V={0,2,...,22} union {1,9}`, reference `0` has tails `(1,9)` and fails
   the pair-mass gate, while reference `22` has tails `(13,21)` and passes it.
   A convenient reference does not transfer loneliness to the requested
   runner.

## 6. Connection contract and scope

```text
source:       collision-free degree-two anchored row 2H union {a,b}
target:       a finite labelled packet and two literal lift masks
map:          y=r/N followed by x=(r+jN)/(2N)
preserved:    distinguished anchor, every body constraint, both tail residues,
              dyadic sheet, and strict-danger/closed-safety convention
destroyed:    tail magnitude above its exact periodic residue at the chosen N
sidecar:      canonical r label, N, and owner bit at every packet point
decisive test: R_N(H)=empty for one N in the declared clock bank.
```

This theorem has the following strict boundaries.

- It begins after THM-4330 has produced a collision-free degree-two anchor.
  It does not change reference and does not close degree-twelve minority
  anchors or two-adic types `k>=3`.
- The new exact work covers only divisor-complete exact maxima `31,...,40`.
  Non-divisor-complete rows and the range through `30` are inherited from
  THM-4070 and THM-4075 respectively.
- It is a direct owner-labelled bypass.  It does not map an arbitrary body
  into THM-4338's fixed thirty-label pool.
- It imposes no tail-height bound, but it does impose the half-body bound
  `max(H)<=40`.  No arbitrary-entry or all-height consequence follows.
- LRC(14) remains open.

## 7. Reproduction

From the repository root:

```bash
clang++ -std=c++20 -O2 -Wall -Wextra -Werror -pedantic \
  04-computation/lrc14_arbitrary_entry_owner_word_through40_probe_arbitrary_entry_20260902.cpp \
  -o /tmp/lrc14_arbitrary_entry_through40

/tmp/lrc14_arbitrary_entry_through40 | diff -u \
  05-knowledge/results/lrc14_arbitrary_entry_owner_word_through40_probe_arbitrary_entry_20260902.out -

/tmp/lrc14_arbitrary_entry_through40 --reverse-clocks | diff -u \
  05-knowledge/results/lrc14_arbitrary_entry_owner_word_through40_reverse_probe_arbitrary_entry_20260902.out -

python3 -B \
  04-computation/lrc14_arbitrary_entry_gate_separation_probe_arbitrary_entry_20260902.py \
  | diff -u \
  05-knowledge/results/lrc14_arbitrary_entry_gate_separation_probe_arbitrary_entry_20260902.out -
```

The descending run is a reordered replay, not a byte comparison with the
ascending transcript.  **QED.**
