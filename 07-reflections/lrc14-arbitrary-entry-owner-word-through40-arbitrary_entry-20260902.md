# LRC(14) arbitrary entry: the labelled owner-word bypass through height 40

**Status:** `FINITE-EXACT` research result relative to the proved
THM-2066/4070/4075 dyadic-seam interfaces, with a symbolic finite-chart
no-go and exact separation controls.  This is not a new arbitrary-entry map
into THM-4338, not an all-height theorem, and not LRC(14).

## 1. Inheritance pass and concept board

The object is an already anchored collision-free degree-two row

```text
S(H;a,b)=2H union {a,b},
|H|=11, H subset Z_(>0), a,b distinct positive odd.       (1)
```

The closest proved mechanisms are:

- THM-4070, which closes `(1)` if some `d in {2,...,14}` divides no
  member of `H`;
- THM-2066, whose canonical core-safe residues and binary tail-owner words
  form a lossless finite certificate at a chosen clock;
- THM-4075, which applies that certificate to all divisor-complete cores
  through maximum `30`;
- THM-4330, whose projective fixed-pool and pair-adaptive Haar gates are
  strong but deliberately scalar;
- THM-615, which closes the AP body `H={1,...,11}` for every odd tail pair.

The canonical hostile for the current gates is THM-4330's anchored AP row
with tails `(1,9)`.  A sharper separation body is

```text
H_*={1,2,3,4,8,9,10,11,12,13,14}.                     (2)
```

It is divisor-complete, so it survives THM-4070, but it fails both displayed
THM-4330 entry tests at the worst tail ratio.  The corrected near misses are
the change-of-observer error in THM-4330, the reversed deletion direction,
and MISTAKE-238's attempt to transport the empty full-row safe set through a
homeomorphism whose domain starts only at the nonempty quotient core.  The
least-used relevant sidecar is the **canonical residue label** inside the
core packet: without it one cannot say which dyadic lift an odd tail owns.

The live concept board was:

1. projective fixed-pool refactorization;
2. pair-adaptive Haar mass;
3. AP/folding argmax suppliers;
4. deletion and restoration of the deleted constraint;
5. multiple collision-free references;
6. labelled owner words and literal lift masks.

The sixth object is the one that survives every type audit.

## 2. New finite-exact statement

For every exact maximum

```text
M in {31,32,...,40},                                    (3)
```

the new census exhausts every eleven-element set

```text
H subset {1,...,M}, max(H)=M                            (4)
```

which is divisor-complete in the exact THM-4070 sense

```text
for every d=2,...,14, some h in H satisfies d|h.         (5)
```

For every such `H`, some clock

```text
15 <= N <= 43                                           (6)
```

has empty THM-2066 complementary-owner relation `R_N(H)`.  Therefore, for
**every** pair of distinct positive odd tails `a,b`, the row `(1)` has a
closed `1/14`-safe time.

There is no primitivity filter.  The census includes all `22,765`
nonprimitive bodies in the new layers as well as all `73,262,823` primitive
ones.

Combining this with THM-4070 for bodies failing `(5)` and THM-4075 for
maximum at most `30` gives the finite consequence

> Every eleven-element positive body `H` with `max(H)<=40` closes
> `2H union {a,b}` for every two distinct positive odd tails.

The **newly enumerated part** is only the divisor-complete exact-max layers
`31,...,40`; the statement for all bodies through `40` uses the proved
symbolic non-divisor certificate.  Nothing here bounds an arbitrary body's
maximum by `40`.

## 3. Exact universe and filter

For exact maximum `M`, the raw universe is

```text
binom({1,...,M-1},10) union {M},                        (7)
```

so it has `binom(M-1,10)` labelled-by-value bodies.  The exact counts are:

| `M` | raw universe | divisor-complete | primitive | nonprimitive | survivor |
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

The three bodies first certified at the endpoint clock `43` in the ascending
bank are

```text
{11,14,19,23,25,26,29,34,36,37,40},
{12,18,19,22,25,26,28,31,34,37,40},
{19,22,25,26,28,29,31,34,36,37,40}.                   (8)
```

Thus the endpoint is genuine for this ordered bank, although no claim of
global minimality is made.  All three `(8)` bodies have Haar mass strictly
above `4/63`; they are hard for the owner-clock ordering but easy for the
pair-adaptive mass gate.  Clock hardness and mass hardness are orthogonal.

## 4. Why this is a faithful entry map

At a clock `N`, retain the canonical packet

```text
A_N(H)={0<=r<N: 14|hr|_N>=N for every h in H}.          (9)
```

For an odd residue `z mod 2N`, retain the actual strict-danger mask on both
lifts

```text
x_(r,j)=(r+jN)/(2N),             j in {0,1}.            (10)
```

Equivalently, when `z` is dangerous above every `r in A_N(H)`, retain the
nearest-integer parity word telling which lift it kills.  Two tails obstruct
the clock exactly when their words are complementary on every labelled
packet point.  Empty `R_N(H)` says no odd residue pair, with repetition
allowed modulo `2N`, can kill both lifts everywhere.  Allowing repeated
residue classes is necessary because distinct physical odd tails may be
congruent modulo `2N`.

The connection contract is:

```text
source:       anchored physical row 2H union {a,b}
target:       labelled packet A_N(H) and two literal lift masks
map:          y=r/N, then x=(r+jN)/(2N)
preserved:    distinguished anchor, all body constraints, both physical
              tail residues mod 2N, sheet labels, and strict/closed boundary
destroyed:    tail magnitudes beyond their exact periodic residue at this N
restoration:  none needed at the chosen clock; evaluation is periodic mod 2N
sidecar:      canonical r labels, clock N, and owner parity on every r
decisive test: no unordered odd residue pair with repetition covers both
               literal lifts over all r in A_N(H).
```

This bypasses the fixed pool rather than entering it.  It is faithful because
it evaluates the physical row at actual rational times and preserves the
tail-owner coordinate that Haar mass discards.

## 5. Exact separation of the tempting operations

### 5.1 Projective pool and scalar Haar mass

For `(2)`, the best lawful positive-rational refactorization has only six
pool hits, attained at scale `10`; THM-4330 needs nine.  Exact wall
integration gives

```text
mu(G_(H_*))=20411/360360 < 4/63=mu(C_(1,9)).            (11)
```

Thus both THM-4330 gates fail for tails `(1,9)`.  At clock `23`, however,

```text
A_23(H_*)={4,9,10,13,14,19},                           (12)
```

and the complete literal odd-residue-pair audit has no covering pair.  This
closes every odd tail pair.  For the displayed `(1,9)` control, `x=2/23`
has clearance exactly `2/23`.

The AP body `{1,...,11}` gives a second separation:

```text
projective pool hits = 6 at scale 10,
mu(G_H)=10931/194040 < 4/63,
x=5/24 has full-row clearance 1/12 for tails (1,9).     (13)
```

Here the correct alternative supplier is THM-615's AP folding theorem.

### 5.2 Deletion has the wrong order unless the phase is restored

If `J subset H`, then

```text
G_H subset G_J,                                        (14)
```

so a lower bound or witness for the deleted body does not lower-bound or
witness the full body.  The smallest exact control is

```text
J={1,...,10}, H={1,...,11}, x=1/11:                    (15)
min_(j in J)||jx||=1/11, but ||11x||=0.
```

Deletion becomes useful only with the deleted constraint evaluated at the
same phase, precisely the sort of address retained by `(9)--(10)`.

### 5.3 A convenient reference is not the requested observer

Take

```text
V={0,2,4,...,22} union {1,9}.                          (16)
```

The collision-free anchors `0` and `22` have the same half-body
`{1,...,11}`.  At anchor `0` the primitive tails are `(1,9)` and the
pair-adaptive threshold is `4/63`, so `(13)` fails the mass gate.  At anchor
`22` the tails are `(13,21)` and the threshold is `2/49`, so the same body
passes.  This is an exact same-configuration/opposite-verdict hostile:
re-referencing changes the runner whose loneliness is proved.

### 5.4 A finite literal chart cannot imply arbitrary singleton labels

There is also a simple symbolic stopping lemma.  Let a fixed body `B` have a
safe-set component containing an interval of length `L>0`.  Every component
of the singleton safe set `G_{ {h} }` has length `6/(7h)`.  Therefore

```text
G_B subset G_{ {h} }  implies  h <= 6/(7L).            (17)
```

Consequently any finite family of fixed, unscaled, positive-width charts can
pointwise imply only finitely many singleton labels.  This is the general
mechanism behind THM-4332's exact rigidity.  Its concrete pool hostile is

```text
199/21280 in G_P intersect D_1.                        (18)
```

Statement `(17)` does not forbid projectively rescaling a chart, using a
mass comparison, or retaining owner words.  It forbids only the tempting
finite unscaled literal-containment entry map.

## 6. Audit and reproduction

The main program performs the complete `(7)` enumeration, with no random
sampling and no gcd filter.  For every newly encountered packet it computes
the THM-2066 owner-word predicate and, independently inside the same binary,
enumerates every unordered odd residue pair modulo `2N` with repetition and
tests both literal lifts `(10)`.  The predicates agree on every packet.

The ascending run audits `4,168` distinct packet/clock states.  A descending
clock replay audits `18,128` states, changes the selected-clock histograms as
it should, and returns the same raw/filter/primitive/nonprimitive totals and
zero survivors.  Strict `-O2` and `-O3` ascending builds produce byte-identical
transcripts.

Artifacts:

```text
04-computation/lrc14_arbitrary_entry_owner_word_through40_probe_arbitrary_entry_20260902.cpp
05-knowledge/results/lrc14_arbitrary_entry_owner_word_through40_probe_arbitrary_entry_20260902.out
05-knowledge/results/lrc14_arbitrary_entry_owner_word_through40_reverse_probe_arbitrary_entry_20260902.out
04-computation/lrc14_arbitrary_entry_gate_separation_probe_arbitrary_entry_20260902.py
05-knowledge/results/lrc14_arbitrary_entry_gate_separation_probe_arbitrary_entry_20260902.out
```

Raw-LF SHA-256 values are:

```text
e0298461737adc15d2d84d55f86298ad4419fe5c9e74ff0fb97427d8371ebade  owner-word C++
9b1f34c9caa52eceb9cd2db67c3cb12b5a00ca5befd7764fa2d2a6e7d746d401  ascending output
e871318601799b0028de6c80eb6ceeb2d4c19c6ee910362a2d2a969580c4924f  descending output
277bb7e517cd32f90c775f50c7143faf88c531e8957fdf632b04ed2a0af4e261  separation Python
58a22fd47489981076e95b86f2ccbbe5a72599f000e758c345e053068823cea8  separation output
```

Reproduce from the repository root:

```bash
clang++ -std=c++20 -O2 -Wall -Wextra -Werror -pedantic \
  04-computation/lrc14_arbitrary_entry_owner_word_through40_probe_arbitrary_entry_20260902.cpp \
  -o /tmp/lrc14_arbitrary_entry_through40
/tmp/lrc14_arbitrary_entry_through40
/tmp/lrc14_arbitrary_entry_through40 --reverse-clocks

python3 -B \
  04-computation/lrc14_arbitrary_entry_gate_separation_probe_arbitrary_entry_20260902.py
```

## 7. Next sharp problem

Blindly extending one maximum at a time is now less informative than
intersecting the gates.  The next search object should be a body `H` with all
of the following coordinates retained:

```text
divisor-complete through 14;
failure of every lawful projective pool refactorization;
mu(G_H)<mu(C_(p,q)) for its actual primitive tail ratio;
nonempty complementary-owner relation at every tested clock;
the original anchor and the two actual tail residues.
```

The fact that all three clock-43 bodies in `(8)` are mass-easy is positive
signal for a **hybrid exclusion theorem**: perhaps owner-word hardness forces
enough safe mass, while mass-hard bodies admit an early owner clock.  The
cheapest decisive next test is to enumerate this intersection, not either
failure set separately, beginning at exact maximum `41` and recording the
first body that is simultaneously mass-hard and owner-hard.

Arbitrary all-height entry, a map into THM-4338's fixed pool, and LRC(14)
remain open.
