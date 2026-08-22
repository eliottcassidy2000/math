---
id: THM-3366
title: All-sector pointwise complement-clock completion
status: >
  PROVED + FINITE-EXACT, with CITED input; the specialized k=1 census, unified
  all-k pool-14 census, and refined k=2 and k=3 compositions are INDEPENDENTLY
  HOSTILE-AUDITED.
  On every THM-2928 body-address quotient, a pointwise cover of the
  unsupported open cells by r integer danger combs compiles any hypothetical
  k-aligned, (7-k)-drift cover into a global cover by at most 7+r nonzero
  clocks when k>=1, and at most 8+r clocks when k=0.  Cited LRC through
  twelve nonzero speeds therefore closes every row with r<=5 for k>=1 and
  r<=4 for k=0.  An exact strict-endpoint set-cover census over clocks 1..14
  closes 17,561 k=0 rows, 19,273 k=1 rows, 19,198 k=2 rows, and 19,053 k=3
  rows in the raw support ledgers.  A separately replayed 1..28 extension
  closes 17,575, 19,301, 19,226, and 19,081 rows respectively.  These are
  support-row terminals, not an assertion that every surviving row is
  physically realizable.  Exact keywise composition with the current k=2
  post-d6 ledger deletes 298 live rows and 71,575,318 occurrences, leaving
  4,056 rows and 200,069,517,203 occurrences.  Exact keywise composition
  with the current k=3 one-spike ledger deletes 7 live rows and 7,648
  occurrences, leaving 1,897 rows and 2,548,893,834 occurrences.  Other
  refined-ledger decrements require the same explicit intersection rather
  than subtraction of overlapping screens.
source: codex-kps-s174-2026-08-14
depends_on:
  - THM-2928-critical-seven-comb-grid-tensorization-and-drift-tariff
  - LRC(<=13)
related:
  - THM-3363-d14-complement-clock-small-lrc-terminal
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
script: 04-computation/lrc14_allk_universal_complement_clock_scan_kps_s174.py
output: 05-knowledge/results/lrc14_allk_universal_complement_clock_scan_kps_s174.out
script_sha256: 372f1b0d2bd8c1cc453080e9ab55880480352c5277a9a216739a7262d198efba
output_sha256: 1104900cd805a2af05e3c7252b7d7bacd417b36b0f015a90d1976b6913fc91c5
allk_audit_script: 04-computation/lrc14_allk_complement_clock_independent_audit_20260814.py
allk_audit_output: 05-knowledge/results/lrc14_allk_complement_clock_independent_audit_20260814.out
allk_audit_script_sha256: a43610abcb1d591f4b606677e3e6a9a32897b38fa2fd591a281ee8c4dc07d823
allk_audit_output_sha256: 31c096e528a148cd717f4d9c7485c93cdec6175b33b83bae9d2be627249d743c
k1_companion_script: 04-computation/lrc14_k1_universal_complement_clock_scan_kps_s173.py
k1_companion_output: 05-knowledge/results/lrc14_k1_universal_complement_clock_scan_kps_s173.out
k1_companion_script_sha256: 65f1e5983b6d9ee4faf0d74f276c20821f85031f33a100d61c03942da3a8971b
k1_companion_output_sha256: 80d37ede35c05d7b95779be31a4bef487066869af3c40c1f5f923707df3fa5cf
composition_script: 04-computation/lrc14_k2_refined_complement_clock_composition_kps_s175.py
composition_output: 05-knowledge/results/lrc14_k2_refined_complement_clock_composition_kps_s175.out
composition_script_sha256: 414b3777cc44f2e059d5bc4258ab555fc055de8c0f8b3ebfd99ce0c90c7da14c
composition_output_sha256: 5fe0c922cf89754a9d557f3d535ab65e71995fef17fd91f123cf5673dcfa6dfd
k3_composition_script: 04-computation/lrc14_k3_refined_complement_clock_composition_kps_s176.py
k3_composition_output: 05-knowledge/results/lrc14_k3_refined_complement_clock_composition_kps_s176.out
k3_composition_script_sha256: 27e4bff52705189bf8ff73db42d76d4e2fc94c44330d295f166f8d4217cb1804
k3_composition_output_sha256: 4cb8f95113123007af9fb5a1f58b3b5373dd4637615b50be7168c6bf578b696b
k3_audit_script: 04-computation/lrc14_k3_refined_complement_clock_independent_audit_20260814.py
k3_audit_output: 05-knowledge/results/lrc14_k3_refined_complement_clock_independent_audit_20260814.out
k3_audit_script_sha256: 1204bff675273efb73241cbc511e84dd2cbc11bae4ce5dcc5caa258928f1064f
k3_audit_output_sha256: 0bf0d3b3c6b530ed3b37b3cf3d6c9dd38edebf331e2f6c7fcbcca1741a5e15c0
hash_basis: LF-normalized bytes
---

# THM-3366 -- all-sector pointwise complement-clock completion

This theorem turns missing body-address cells into actual lower-dimensional
runner clocks.  Unlike the support-mass tariff, the compiler is pointwise and
has no numerator-height parameter.

## 1. Inherited quotient and the missing-cell object

Write

```text
T=R/Z,                    D_w={x in T: ||wx||<1/14}.             (1)
```

Fix a literal six-body set `F subset {1,...,14}` and put

```text
L=14 lcm(F).
```

THM-2928 partitions the exact body-safe carrier into `1/L` cells.  For every
divisor `D|L`, let `S_D subset Z/DZ` be its projected safe-address support.
If `A` is the set of `k` aligned integer multipliers, write

```text
R_A=T minus union_(a in A) D_a,
Y_D(A)=union_(r in S_D) (r+R_A)/D.                            (2)
```

For `k=0`, use `A=empty` and `R_empty=T`.  If the other `p=7-k` tails have
reduced denominators with lcm `D`, THM-2928 proves that every physical cover
forces a pointwise quotient cover

```text
Y_D(A) subset union_(i=1)^p D_(b_i)                          (3)
```

by positive integer quotient clocks `b_i`.

Define the unsupported open-cell set

```text
U_D(S_D)=union_(r notin S_D) (r/D,(r+1)/D).                  (4)
```

The endpoints are deliberately absent from `(4)`.

## 2. Exact complement inclusion

For `k>=1`,

```text
T minus Y_D(A)
 subset U_D(S_D) union union_(a in A) D_(Da).                 (5)
```

Indeed, let `x` lie in a supported open cell and write `x=(r+u)/D`.
If `x` is outside `(2)`, then `u` is outside `R_A`, so for some `a in A`,

```text
||Da x||=||a(r+u)||=||au||<1/14.                             (6)
```

Thus `x in D_(Da)`.  Every grid point `x=r/D` lies strictly in every
`D_(Da)` because `Da x=ar` is integral.  Points in unsupported open cells
are exactly the first term on the right of `(5)`.

When `k=0`, there is no aligned clock to own the grid points, but the single
integer clock `D` does:

```text
T minus Y_D(empty) subset U_D(S_D) union D_D.                (7)
```

This is the only difference between the zero-aligned and positive-aligned
sectors.

## 3. Complement-clock compiler

Suppose positive integer clocks `c_1,...,c_r` cover `(4)` pointwise:

```text
U_D(S_D) subset union_(j=1)^r D_(c_j).                       (8)
```

If `(3)` also held, equations `(5)` and `(8)` would give a global circle
cover by

```text
{b_1,...,b_p} union {Da:a in A} union {c_1,...,c_r}.          (9)
```

There are at most

```text
p+k+r=7+r                                                     (10)
```

distinct nonzero clocks.  Repetitions only lower the count.  For `k=0`, use
`(7)` and include `D`, giving at most

```text
p+1+r=8+r.                                                    (11)
```

The cited Sungkawichai--Trakulthongchai computation proves LRC through
twelve nonzero speeds.  Any set of at most twelve positive clocks therefore
has a time at which every circular distance is at least `1/13`, which is
strictly larger than `1/14`.  Such clocks cannot cover `T` by the open danger
combs `(1)`.  Hence:

```text
k>=1 and r<=5  ==> the body/divisor support row is impossible;
k=0  and r<=4  ==> the body/divisor support row is impossible. (12)
```

This closes every aligned shape and every quotient numerator over the row.
No bound on those numerators appears in the proof.

## 4. Finite strict-endpoint set cover

For the primary finite certificate take

```text
C_14={1,2,...,14}.                                           (13)
```

For a clock `c`, the boundary of `D_c` in `[0,1]` consists of the rational
points

```text
(14j-1)/(14c),  (14j+1)/(14c)                               (14)
```

that lie in the unit interval.  Let `P_14` be the sorted union of all these
points together with `0,1`.  It has `206` points and `205` open atoms.
Every Boolean danger word for the fourteen clocks is constant on each open
atom.

For each body/divisor row, the exact target mask has two kinds of bits.

1. Every point of `P_14` lying strictly in `(4)`, except a `D`-grid point,
   is a bit.  This retains all strict-equality failures of the open combs.
2. Every open arrangement atom having nonempty intersection with `(4)` is a
   bit.  Since every candidate danger word is constant there, one midpoint
   determines the whole atom.

A recursive bitset set-cover solver chooses an uncovered bit, branches over
every candidate clock that covers it, and memoizes `(remaining,depth)`.  It
tries depths in increasing order through five.  Thus it finds a cover of
size at most five if and only if one exists in `(13)`; for `k=0`, only answers
of size at most four are accepted.  Grid points omitted from the target are
owned by `(5)` or `(7)`, so no endpoint is discarded.

The exact support universe and cutoffs imported from THM-2928 are

```text
k                 0        1        2        3        4       5      6
p                 7        6        5        4        3       2      1
cutoff        139/154  106/117  887/990  125/143   26/31  375/478 39/61.
                                                                    (15)
```

The primary pool-14 census is:

| `k` | input rows | closed rows | remaining rows | closed denominator-shape occurrences | remaining occurrences |
|---:|---:|---:|---:|---:|---:|
| 0 | 27,210 | 17,561 | 9,649 | 113,853,177,071,279 | 1,390,988,884,871,570 |
| 1 | 27,240 | 19,273 | 7,967 | 5,116,993,586,195 | 33,837,732,004,565 |
| 2 | 27,163 | 19,198 | 7,965 | 182,987,150,531 | 768,558,739,704 |
| 3 | 26,970 | 19,053 | 7,917 | 5,887,257,171 | 15,470,456,930 |
| 4 | 13,778 | 5,964 | 7,814 | 47,788,532 | 250,467,350 |
| 5 | 10,976 | 3,334 | 7,642 | 528,176 | 2,538,098 |
| 6 | 6,237 | 570 | 5,667 | 570 | 5,667 |

The `k=4,5,6` sectors were already closed globally by THM-2928; their rows
are controls, not new LRC progress.  The `k=0,1,2,3` lines are the new raw
support-row terminals.  The `k=1` count includes the one `D=14` row and its
`26` denominator shapes already closed by THM-3363.  Therefore the pool-14
increment beyond THM-3363 is `19,272` rows and
`5,116,993,586,169` occurrences, not the unadjusted total.

The `k=1` line is also produced by the separately implemented specialized
companion S173.  Its exact event sweep, strict-open atom cover and census were
independently reconstructed; the hardened stored run pins its dependency,
row, histogram, occurrence and semantic hashes.  This is an independent
check of that sector, not a second theorem to be added to the ledger.

For `p>=1`, the occurrence weight of a body/divisor row is the exact number
of unordered denominator multisets with lcm `D`:

```text
a_p(D)=sum_(e|D) mu(D/e) binom(tau(e)+p-2,p).                (16)
```

This is divisor-lattice Mobius inversion: `tau(e)-1` divisors of `e` are
greater than one, and their size-`p` multisets number the displayed binomial.
The script checks `(16)` by literal enumeration for `D<=40`, `1<=p<=7`, and
also reproduces the entire THM-2928 input row and occurrence ladders.

## 5. Pool-28 extension

The same exact algorithm was replayed with candidate clocks

```text
C_28={1,2,...,28}.                                          (17)
```

Here the arrangement has `780` points and `779` atoms.  Ordinary and
optimized runs agree and give:

| `k` | closed rows | remaining rows | closed occurrences | remaining occurrences |
|---:|---:|---:|---:|---:|
| 0 | 17,575 | 9,635 | 113,911,192,737,582 | 1,390,930,869,205,267 |
| 1 | 19,301 | 7,939 | 5,164,491,031,504 | 33,790,234,559,256 |
| 2 | 19,226 | 7,937 | 184,622,491,967 | 766,923,398,268 |
| 3 | 19,081 | 7,889 | 5,931,229,159 | 15,426,484,942 |
| 4 | 5,991 | 7,787 | 48,610,947 | 249,644,935 |
| 5 | 3,360 | 7,616 | 536,545 | 2,529,729 |
| 6 | 570 | 5,667 | 570 | 5,667 |

This shows strong diminishing returns after clock `14`: the larger pool adds
only `14,28,28,28` rows in the four open support sectors.  It is a stronger
finite certificate, but it does not prove that clocks above `28` are useless.

## 6. Exact composition with the refined `k=2` ledger

The non-subtraction guard is numerically significant.  The companion

```text
04-computation/lrc14_k2_refined_complement_clock_composition_kps_s175.py
```

reconstructs the post-septimal, post-located-`d=6` ledger of THM-2928 and
intersects it row by row with the pool-14 terminal set.  It independently
reproduces the inherited residual

```text
(shapes,occurrences,rows,bodies,divisors)
=(26,899,164,786, 200,141,092,521, 4,354, 2,966, 147).       (18)
```

Only `298` of those live refined rows are complement-clock terminals.  They
carry exactly `71,575,318` refined occurrences.  Deleting their actual keys,
not the raw support weights, gives

```text
(26,899,164,065, 200,069,517,203, 4,056, 2,966, 144).       (19)
```

The occurrence split by transverse count `c=1,...,5` is

```text
(3,524,756,
 46,320,209,782,
 112,778,172,674,
 39,482,116,896,
 1,485,493,095).                                            (20)
```

The C++ floor/exception engine is compiled independently under `-O2` and
`-O3`; both give the frozen engine transcript hash.  Ordinary and optimized
Python composition runs agree.  Equation `(19)` is the honest current `k=2`
ledger decrement supplied by this theorem.

## 7. Exact composition with the refined `k=3` ledger

The companion

```text
04-computation/lrc14_k3_refined_complement_clock_composition_kps_s176.py
```

reconstructs the final one-spike ledger of THM-2928 and independently
reproduces

```text
(shapes,occurrences,rows,bodies,divisors)
=(398,241,574, 2,548,901,482, 1,904, 1,823, 107).           (21)
```

The exact row-key intersection is nearly orthogonal to the raw support
census: only `7` live refined rows are pool-14 complement-clock terminals,
carrying `7,648` occurrences.  Their bodies are precisely

```text
{2,6,8,10,14} union {u},       u in {1,3,5,7,9,11,13},     (22)
```

and every terminal has `D=L/2`.  For `u=1,3,5,7` the divisor is `5,880`;
for `u=9,11,13` it is respectively `17,640`, `64,680`, `76,440`.  Thus the
seven computational hits form one odd-parameter divisor-chain family rather
than seven unrelated accidents.

Deleting those exact keys gives

```text
(398,241,574, 2,548,893,834, 1,897, 1,823, 107).           (23)
```

All `7,648` deleted occurrences lie in transverse count `c=3`, changing that
slice from `1,065,317,472` to `1,065,309,824`; the other three slices and the
shape union are unchanged.  Ordinary and optimized Python runs agree byte
for byte.  The sparse intersection is also a strategic negative result:
after the one-spike screens, enlarging an unstructured integer-clock pool is
unlikely to be the main route through the remaining `k=3` ledger.  The
reusable positive signal is instead the rigid core and `D=L/2` law in `(22)`.

An independent hostile audit reconstructs both inputs without the primary
composition engine: exhaustive subset enumeration recovers all `19,053`
pool-14 terminals, while a literal current-lcm dynamic program and primal
Boolean-coupling vertices recover the `1,904` one-spike rows.  Their exact
`(body,D)` intersection is again the seven rows and `7,648` occurrences in
`(22)`.  Body-only, divisor-only, and `(L,D)` keys delete respectively `1,904`,
`1,880`, and `33` rows, while raw or double subtraction gives wrong occurrence
counts.  These hostiles make the keywise, non-subtractive composition
load-bearing rather than cosmetic.

## 8. Scope and non-subtraction guard

The theorem closes the listed **support rows**.  It does not claim:

* that any remaining support row is realized by physical tail clocks;
* that a bounded clock pool detects every possible complement cover;
* that the raw occurrence deletions can be subtracted from later septimal,
  spike, phase, or translated-residue ledgers; or
* LRC(14).

Later screens overlap this one.  Their sharper residual must be obtained by
intersecting exact row or occurrence keys.  Sections 6 and 7 perform that
operation for the current `k=2` and `k=3` ledgers.

## 9. Independent unified-census audit

The independent companion rebuilds all `11,180` distinct pool-14 targets by
exhaustively enumerating clock subsets, rather than calling the primary
recursive solver.  It separately reconstructs all `251,536` body/divisor rows,
the strict endpoint and open-atom target bits, and the occurrence weights by a
divisor--Mobius formula checked on `280` brute-force cases.  Every input count,
terminal count, occurrence count, and least-cover-size histogram in all seven
sectors agrees exactly.

Three hostile witnesses show that the conventions are active: midpoint-only
coverage misses `16` strict endpoint bits at `(F,D)=({1,2,3,5,9,10},315)`;
omitting grid ownership changes `20` target bits at
`({1,2,3,4,5,13},3640)`; and that same row needs five added clocks, while the
`k=0` row `({1,2,3,4,5,6},420)` needs four.  This audit certifies the unified
pool-14 census, not the separate refined-`k=3` key intersection.

## 10. Reproduction

```bash
python 04-computation/lrc14_allk_universal_complement_clock_scan_kps_s174.py
python -O 04-computation/lrc14_allk_universal_complement_clock_scan_kps_s174.py
python 04-computation/lrc14_allk_complement_clock_independent_audit_20260814.py
python -O 04-computation/lrc14_allk_complement_clock_independent_audit_20260814.py
python 04-computation/lrc14_k1_universal_complement_clock_scan_kps_s173.py
python -O 04-computation/lrc14_k1_universal_complement_clock_scan_kps_s173.py
python 04-computation/lrc14_allk_universal_complement_clock_scan_kps_s174.py --pool-max 28
python -O 04-computation/lrc14_allk_universal_complement_clock_scan_kps_s174.py --pool-max 28
python 04-computation/lrc14_k2_refined_complement_clock_composition_kps_s175.py
python -O 04-computation/lrc14_k2_refined_complement_clock_composition_kps_s175.py
python 04-computation/lrc14_k3_refined_complement_clock_composition_kps_s176.py
python -O 04-computation/lrc14_k3_refined_complement_clock_composition_kps_s176.py
python 04-computation/lrc14_k3_refined_complement_clock_independent_audit_20260814.py
python -O 04-computation/lrc14_k3_refined_complement_clock_independent_audit_20260814.py
```

All decisions use integer arithmetic or `Fraction`.  Runtime checks remain
active under optimized Python.

**End of proof.**
