# Rank-two completion of the fixed-pool arbitrary two-outsider theorem

**PROVED RELATIVE TO THM-4231 + FINITE-EXACT + CLEAN-ROOM INDEPENDENTLY
AUDITED.  NO PHYSICAL ENTRY AND NO LRC(14).**

Let

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.
```

For every two distinct positive outsiders `q,r notin P` and every
`B in C(P,9)`, this packet proves

```text
mu(G_(B union {q,r})) >= 4/63.                         (1)
```

Here `G_S={t in R/Z:min_(s in S)||st||>=1/14}`.  The new finite computation
closes THM-4231's exact patched remainder of 181,194 pairs; THM-4231 already
closes every other outsider pair.  Thus (1) is an arbitrary-pair theorem for
this fixed pool, not just an endpoint-prefix statement.

## 1. The two-defect graph lemma

For a fixed pair `(q,r)`, partition the circle at every wall of
`P union {q,r}` on the exact integer grid `D`.  On a pair-safe open cell let
`F subset P` be its pool-failure set.  Retain only cells with `|F|<=2`, and
aggregate their exact widths as

```text
w0              for F=empty,
w_i             for F={i},
w_ij=w_ji>=0    for F={i,j}.
```

The wall set is finite.  The computation sums the widths of the open cells
between consecutive distinct walls and does not assign mass to the wall
points themselves.  Those finitely many points are Haar-null, so this
finite-null-wall convention leaves every measure statement unchanged (and
is conservative even though equality is safe at the threshold).

For a nine-label body `B`, put

```text
L2(B)=w0 + sum_(i notin B) w_i
          + sum_({i,j} subset P\B) w_ij.               (2)
```

Every cell counted in (2) is literally safe for `B union {q,r}`, hence
`L2(B)/D <= mu(G_(B union {q,r}))`.  Write

```text
W=w0+sum_i w_i+sum_(i<j)w_ij,
a_i=w_i+sum_(j!=i)w_ij.
```

Then exact inclusion-exclusion on the weighted graph gives

```text
L2(B)=W-sum_(i in B)a_i+sum_({i,j} subset B)w_ij.      (3)
```

If `a_(1)>=...>=a_(30)` are the sorted weighted degrees, nonnegativity of
the internal-edge correction implies the uniform all-body bound

```text
L2(B) >= W-sum_(k=1)^9 a_(k).                          (4)
```

Equivalently, choosing `B` covers every singleton loop and ordinary edge
incident to `B`; the naive degree sum double-counts internal covered edges,
which (3) adds back.  The observable `w_ij` is symmetric.  This is an
undirected weighted graph, not a tournament.

## 2. Exact THM-4231 proof surface

The exporter executes THM-4231's maintained patched-residual postprocessor
and freezes

```text
E_rem count       181,194
E_rem FNV         3874fecac4ecbd8a
E_rem SHA256      9dfbf0a8948bf23016ae40f40f9118020c9429ac60421681a9286fc4d34041a1
maximum endpoint  769.
```

The frozen ledger is in the postprocessor's lexicographic `(q,r)` order.
Its SHA256 pins the bytes, its FNV pins the same ordered pair stream, and the
verifier separately checks uniqueness, outsider typing, order, and equality
of the degree-ledger pair set with this source set.

The later THM-4287 universe is an exact subset:

```text
|U|=22,647, FNV=df5374d4aca67677,
SHA256=14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317,
|E_rem\U|=158,547.
```

This arithmetic is audited only as a provenance control.  The proof of (1)
uses the complete THM-4231 remainder directly and therefore does not depend
on the later carrier descent.

## 3. Complete exact result

On all 181,194 remainder pairs, the wide signed-128-bit degree scan finds

```text
strict degree-bound rows  181,087
coarse nonpositive rows       107
highest coarse exception      (50,554), uniquely at endpoint 554
full degree-ledger FNV         f8da82c5a6d732ed
full degree-ledger SHA256      5b008404482b6e23006c0c1cb97407c22ddf34eb7cf4e09fd5f604206d8a0356.
```

Thus (4) proves (1) simultaneously for all `C(30,9)=14,307,150` bodies on
181,087 pairs, without enumerating bodies.  For the 107 remaining pairs, a
flat complete all-body enumeration and a structurally separate exact
branch-and-bound optimizer agree on every minimizing body and margin.  All
107 minima are strict.  The branch upper bound is

```text
current score + sum of the largest needed current marginals;
```

it is rigorous because it omits only nonnegative pair penalties among future
choices.  The entire branch replay visits 55,104 nodes and prunes 49,891.

The 107 exceptions are not failures of (1).  They are exactly the rows where
the degree-only relaxation loses too much internal-edge overlap.  The first
such relaxation failure is `(50,554)`, while its exact rank-two minimum is
strictly positive.

## 4. Sharp normalized finite control

Raw integer ticks from different grids must not be ordered.  Intrinsic
comparisons use `ticks/D`, where

```text
ticks=63 L2(B)-4D,
rank-two certified surplus=(L2(B)/D)-4/63=ticks/(63D),
actual Haar surplus=mu(G)-4/63 >= ticks/(63D).          (5)
```

Among the 107 exactly optimized rows, the smallest normalized margin is

```text
(q,r)=(50,70), B=031c7400
ticks=245,428,469,244, D=91,205,797,082,400,
ticks/D=973922497/361927766200.
```

Exactly three degree-positive rows have a coarse ratio no larger than this:
`(50,212)`, `(50,274)`, and `(100,110)`.  Complete O2/O3 all-body
optimization gives their exact minima respectively as

```text
4,207,251,055,549,752 / 4,833,907,245,367,200 at 053c6400,
11,006,069,470,557,714 / 12,495,194,200,288,800 at 11187401,
63,178,284,254,904 / 91,205,797,082,400 at 04f06408.
```

All three ratios are larger.  Therefore `(50,70),031c7400` is the genuine
global normalized minimizer of the rank-two lower bound on the frozen
181,194-pair universe.  Its exact rank-two certified surplus over `4/63` is

```text
973922497/22801449270600 > 0.                           (6)
```

The actual Haar surplus is at least (6).  Equality is neither needed nor
claimed: rank-at-least-three cells were deliberately discarded by the graph
projection.

## 5. Overflow and independent controls

Exactly one row needs more than a signed 64-bit wall grid:

```text
(713,719), D=9,351,275,651,380,222,560.
```

The earlier signed-64-bit scout is therefore not a valid proof path for the
181,194-row theorem.  The frozen implementation uses signed 128-bit grids,
widths, degrees, scores, and ticks (the C++ source's historical `i64` alias is
explicitly bound to `__int128_t`).  On the 22,647-row THM-4287 control it
reproduces the prior CSV byte-for-byte while using the corrected two-limb
FNV `611b1ba5c25594dd`; the stale pre-wide FNV `fe6111d297a72a3d` is
explicitly rejected.

O2 and O3 full screens agree byte-for-byte.  O2/O3 flat enumeration and
O2/O3 branch-and-bound agree on all 107 exceptions.  A separate agent's
clean-room wide wall graph matches every field on all 181,194 pairs and its
independent minimization matches every one of the 107 exact minima.  That
implementation is frozen under `independent_cleanroom/`: its event-state
merge graph builder is distinct from the primary sorted-wall midpoint path,
its inner 29-file manifest has SHA256
`1fdd6025bae13f9bca08eb4fb97a82f149d845522e6dd6d9ac7d70629af79302`,
its broad screen has SHA256
`57826d842ee968251696be4757cf94b0ff57f55f3cc88f8020e2a27ccbd19258`,
and its primary/clean-room crosscheck freezes `181194/181194` fieldwise graph
matches and `107/107` exact minimum/body matches.

## 6. Consequence and firewall

THM-4231 proves (1) off `E_rem`; Sections 1--5 prove it on `E_rem`.
Therefore (1) holds for every distinct positive outsider pair.  Applying
THM-4150 to `H=B union {q,r}` also proves every universal odd-tail family

```text
2H union {a,b},   a,b distinct positive odd integers.  (7)
```

This is a carrier-free completion of the fixed-pool arbitrary two-outsider
Haar theorem and its associated structured thirteen-speed odd-tail families.
It does **not** map an arbitrary normalized LRC(14) row into this pool, prove
that every physical row has such a representation, or prove LRC(14).

The graph projection preserves exactly the rank-zero, singleton, and pair
failure widths needed by (2).  It destroys cell addresses, cyclic order,
all failure masks of rank at least three, and carrier/owner information.
Those discarded coordinates cannot be used later without a sidecar.
