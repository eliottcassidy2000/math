---
id: THM-3502
title: "Rule 30 four-fifths staircase entropy and sixteen-fifths compiler"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  The physical
  zero-right-tail Rule 30 staircase language has an
  exact height-23 overlap-graph certificate proving eta<4/5.  In the charged
  word-RAM model of THM-3480 this improves the nonrectangular macroblock query
  tariff from 13/4 to 16/5.
source: root/rule30-sharp-unlocks/height-23-staircase-certificate/2026-08-16
audit: >
  PASS (2026-08-16), independent hostile proof and replay audit.  The auditor
  rebuilt factor closure and the injective overlap-path map, checked the
  direction and exact prefactor of the positive-vector inequality, exhausted
  the floor-sensitive compiler budgets, and retained the spurious-path,
  finite-only, and fixed-state-permutation hostiles.  Ordinary, optimized, and
  stored transcripts are byte-identical.
depends_on:
  - THM-3480-rule30-staircase-transducer-entropy-and-nonrectangular-macroblock-compiler
  - THM-3491-rule30-seven-four-staircase-compiler
related:
  - THM-3456-left-permutive-trace-bijection-and-rule30-seed-boundary
  - THM-3463-rule30-mealy-section-suffix-parity-current-and-complexity-boundary
  - THM-3471-rule30-motzkin-strip-circuit-and-innovation-carry-spectrum
script: 04-computation/rule30_staircase_entropy_four_fifths_thm3502.py
output: 05-knowledge/results/rule30_staircase_entropy_four_fifths_thm3502.out
script_sha256: 5df8496fa715632e37b0498b69a6b2b1c0032d064de6e238f9c0ee32cbb255d7
output_sha256: 7f49456be4d79891586e377524f77d42b531a9129c12e2aa4dd41476ae568f29
hash_basis: raw bytes
---

# THM-3502 -- Rule 30 four-fifths staircase entropy and sixteen-fifths compiler

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-3480 replaces the two raw halos of a Rule 30 rectangle by one physically
reachable diagonal staircase state.  THM-3491 bounds the entropy of that state
by a height-21 factor graph and obtains the charged lookup tariff `13/4`.  The
present theorem pushes the exact factor language two levels farther.  A
positive integer certificate on the height-23 overlap graph gives

```text
eta <= log_2(87/50) < 4/5.                            (1)
```

At the charged table scale, the two elementary inequalities

```text
87^5 < 16*50^5,             3^5 < 2^8               (2)
```

align the state address and the direct ternary marker array.  In the exact
word-RAM model inherited from THM-3480, the proposed main-query count is

```text
Q(n)=(16/5+o(1)) n^2/log_2(n)^2.                     (3)
```

This is an upper-bound compiler for arbitrary finite rows with a physical zero
right exterior, and hence for the row sequence used to obtain one distinguished
center bit.  It is not a lower bound, a universal-row entropy statement, a
single-seed randomness statement, or a solution of any Rule 30 prize.

## 1. Inherited model and typed statement

Use Rule 30 in the form

```text
(F x)_j=x_(j-1) xor (x_j or x_(j+1)).                (4)
```

THM-3480 gives the exact three-state right-to-left Mealy quotient

```text
A=00,             B=01,             C=1*.            (5)
```

For a height-`h` serial cascade, let `Phi_0,Phi_1` be the two transitions
induced by one fresh top-row bit.  Physical zero-right-tail calibration starts
every stage at `A`; put

```text
S_h=Orb_(Phi_0,Phi_1)(A^h),       N_h=|S_h|.         (6)
```

Height concatenation injects `S_(h+k)` into `S_h x S_k`, so `log_2 N_h` is
subadditive and

```text
eta=lim_(h->infinity) log_2(N_h)/h                   (7)
```

exists.  The theorem is:

1. an exact positive integer certificate proves, for every `h>=22`,

   ```text
   N_h < 2^9(87/50)^h;                               (8)
   ```

2. consequently (1) holds; and
3. in THM-3480's charged uniform word-RAM, the staircase chunk table with the
   all-`n` parameters in Section 5 has (3) and retains

   ```text
   T(n)=O(n^2/log^2 n),       S(n)=O(n/log n).        (9)
   ```

The universe in (6) ranges over every finite binary driving word scanned from
the physical zero right exterior.  It is smaller than the universal ternary
cube `{A,B,C}^h`, but larger and differently typed than the one fixed
single-seed orbit.

## 2. Factor closure and the injective overlap-path map

Write `R_r=S_r` as a language of ternary height words.  If `q in S_h`, every
contiguous block of `r` stages in `q` lies in `R_r`.  Indeed, that subcascade
starts at `A^r` and is driven by the binary carrier stream emitted into its
first stage by the preceding part of the cascade.  This is physical factor
closure; it does not say that arbitrary locally compatible factors glue to a
global physical state.

Define `G_23` by

```text
vertices: R_22,
edge u -> v: some r in R_23 has prefix u and suffix v. (10)
```

The adjacency is zero-one: a length-23 word is determined by its length-22
prefix and last symbol, equivalently by its overlapping prefix and suffix.

For every `h>=22`, map a physical word `q in S_h` to its ordered list of
consecutive length-22 factors.  Factor closure makes every adjacent pair an
edge of `G_23`, because their union is a physical length-23 factor.  The
overlapping factors reconstruct every symbol of `q`, so the map is injective.
Thus, if `A_23` denotes the adjacency matrix and `n=h-22`,

```text
N_h <= 1^T A_23^n 1.                                (11)
```

Only injection is used.  Section 6 gives an explicit locally compatible path
which is not globally physical.

## 3. Exact height-23 orbit and integer certificate

Use THM-3491's two-mask encoding.  For `q=(q_0,...,q_(h-1))`, put

```text
E_i=[q_i!=A],       C_i=[q_i=C],
code_h(q)=sum_i E_i 2^i + 2^h sum_i C_i 2^i.         (12)
```

Exact breadth-first orbit construction gives

```text
h       14     15     16      17       18       19
N_h  18619  32885  57741  100901   175680   304714

h       20       21        22         23
N_h  526563   906525   1555372    2660178.           (13)
```

Hence `G_23` has exactly

```text
1,555,372 vertices and 2,660,178 edges.              (14)
```

Order its vertices by increasing `code_22`.  Starting with the all-one vector,
synchronously iterate

```text
v_i^(r+1)=max(v_i^r, ceil(50 sum_(i->j)v_j^r/87)).   (15)
```

The first fixed-point test succeeds on iteration `3502`, after `3501` strict
updates.  The resulting positive integer vector satisfies

```text
min(v)=6,       max(v)=1313,       sum(v)=389727182,
-86 <= 50(A_23 v)_i-87v_i <= 0.                     (16)
```

In particular,

```text
boxed: 50 A_23 v <= 87v.                             (17)
```

The canonical payload is the ascending sequence of pairs

```text
(code_22(q) as unsigned 64-bit little endian,
 v_q        as unsigned 64-bit little endian).       (18)
```

Its SHA-256 digest is

```text
a4ce932f067e0ff21e90ccd3b4216be4dbbba9440e657b4cebc7ff834b302e33. (19)
```

This digest commits to the vector; every row of (17) is independently checked
by the companion and the digest is not used as a substitute for that check.

## 4. Entropy consequence and exact constants

Because the all-one vector is coordinatewise at most `v/6`, positivity and
(17) turn (11) into

```text
N_h
 <= 1^T A_23^(h-22) 1
 <= (389727182/6)(87/50)^(h-22).                     (20)
```

The prefactor in (8) is an exact cleared-integer comparison:

```text
389727182*50^22
 = 9,291,820,096,969,604,492,187,500,000,000,000,000,000,000,000

 < 14,349,764,342,370,835,133,956,480,823,902,911,608,090,938,368
 = 6*2^9*87^22.                                      (21)
```

Equations (20)--(21) prove (8).  Taking logarithms and dividing by `h` gives
the first inequality in (1).  The clean dyadic coarsening is

```text
87^5=4,984,209,207
    <5,000,000,000=16*50^5,                          (22)
```

so `(87/50)^5<2^4` and `log_2(87/50)<4/5`.

## 5. Charged sixteen-fifths compiler

Retain exactly THM-3480's uniform random-access word-RAM:

```text
w=ceil(log_2(n+2)),                                  (23)
```

with unit-cost Boolean word operations, shifts, addition, and addressed word
loads/stores; a constant program; no advice or oracle; and every table and
membership structure built at runtime.

For sufficiently large `w`, set

```text
d=ceil(log_2 w),       L=w-8d-8,
m=floor((L-9)/8),      h=5m,
b=floor((L-9)/2).                                    (24)
```

Small inputs use direct packed simulation.  For the large regime, `m>=5`, so
(8) and (22) give

```text
N_h < 2^9(87/50)^(5m) < 2^(4m+9).                   (25)
```

Write `L-9=8q+r`, with `0<=r<8`.  Then `m=q` and
`b=4q+floor(r/2)`, whence

```text
4m+9+b=8q+9+floor(r/2) <= 8q+9+r=L.                 (26)
```

Thus one dense staircase-state ID and one `b`-bit fresh chunk fit in a word,
and the complete table has

```text
N_h 2^b < 2^L                                       (27)
```

entries.  The direct ternary marker array is charged as in THM-3480.  The
second exact inequality in (2) gives

```text
3^h=(3^5)^m=243^m<256^m=2^(8m)<=2^(L-9).            (28)
```

Since `2^d>=w` and `2^w=O(n)`, equations (24), (27), and (28) imply

```text
2^L=O(n/w^8),       3^h=O(n/w^8).                   (29)
```

THM-3480's deterministic construction costs

```text
O(3^h+hN_h+N_h2^b)                                  (30)
```

word operations and the corresponding marker/table storage.  The middle term
is also lower order because `b` eventually exceeds `log_2 h`; hence (29)
charges every preprocessing term below the main scan and below the active
packed-row storage.

At each macro time, scan the explicitly zero-padded active row from right to
left, starting the height cascade at `A^h`.  The exact diagonal table consumes
`b` fresh top bits, returns the successor staircase state, and emits the
corresponding `b` bits `h` rows lower.  THM-3480's induction proves that the
successive calls emit exactly the physical Rule 30 row.

From (24),

```text
h=(5/8+o(1))w,       b=(1/2+o(1))w,
hb=(5/16+o(1))w^2.                                  (31)
```

Summing the active widths over the `n/h` macrosteps gives

```text
Q(n)
 =n^2/(hb)+O(n/b+n/h)
 =(16/5+o(1))n^2/w^2.                               (32)
```

The fewer-than-`h` final rows use ordinary packed Rule 30 updates.  Equations
(29)--(32) prove the proposed bounds (3) and (9).  Relative to THM-3491's
`13/4`, the lookup constant improves by the factor

```text
(13/4)/(16/5)=65/64.                                 (33)
```

This is a table-query comparison within the declared model, not a leading
constant or lower bound for arbitrary algorithms.

## 6. Exact hostiles and finite lower signal

The overlap graph is strictly larger than the physical language already at
one extra symbol.  Exact height-24 construction gives

```text
N_24=4,535,965.                                      (34)
```

In the top-to-bottom height order, the length-24 word

```text
W=ABCBCBCBCBCAAAAAAAAAAAAA                           (35)
```

has two-mask code

```text
code_24(W)=22,884,124,670.                            (36)
```

Its prefix and suffix are respectively

```text
ABCBCBCBCBCAAAAAAAAAAAA,
BCBCBCBCBCAAAAAAAAAAAAA.                             (37)
```

Both words in (37) lie in `R_23`, but `W` does not lie in `S_24`.  Thus the two
overlapping physical edges define a path in `G_23` which is not the image of a
physical length-24 state.  This explicitly blocks path-surjectivity, equality
in (11), and any entropy lower bound inferred from the overlap graph.

The inherited positive signal does extend one level.  Put

```text
U=CABC,        V=BCBC.                               (38)
```

Exact orbit membership proves

```text
{U,V}^m subset S_(4m)       for every 1<=m<=6.       (39)
```

The companion checks all

```text
2+4+8+16+32+64=126                                  (40)
```

words.  If (39) held for every `m`, it would prove `eta>=1/4`; the present
computation supplies no induction, encoder, synchronizing section, or
all-height gluing theorem.  Therefore (39) is **FINITE-EXACT ONLY** and implies
no lower entropy bound.

The other inherited hostile also remains.  For each fixed staircase state,
the fresh `b`-bit input-to-output chunk map is a permutation.  The compiler
compresses only the physical diagonal boundary states, never the fresh data
bits and never the universal input rank.

## 7. Preservation, loss, and complexity boundary

| source -> target | preserves | destroys / admits | required sidecar |
|---|---|---|---|
| arbitrary finite row -> `S_h` | every physical zero-right-tail staircase | bi-infinite right phases and their boundary state | zero exterior and scan direction |
| `S_h` -> `G_23` path | every consecutive length-23 factor, injectively at the whole-word level | admits paths such as (35) | exact physical orbit, not local factors alone |
| `A_23` -> positive vector | a certified exponential upper rate | exact physical count and every lower bound | ordered payload and rowwise inequality |
| staircase/chunk -> table entry | exact bottom chunk and successor state | intermediate spacetime rows | `h,b`, padding, and dense state ID |
| all finite rows -> fixed seed | a valid upper-bound algorithm | typicality, hardness, balance, nonperiodicity | the distinguished initial row |
| finite block code -> six levels | 126 exact physical witnesses | all later concatenations | an all-height encoder or induction |

The zero-tail restriction is physical for finite packed rows, but it is not a
probability model.  The entropy in (7) counts boundary states reachable from
all finite driving rows; it is not Shannon entropy and says nothing about the
distribution of states along the Rule 30 single-seed orbit.  Likewise, the
word-RAM upper bound is model-specific.  Charging bit operations differently,
allowing advice, changing random-access conventions, or asking for a lower
bound requires a new theorem.

## 8. Exact replay and pending audit

The intended replays are

```bash
python3 04-computation/rule30_staircase_entropy_four_fifths_thm3502.py
python3 -O 04-computation/rule30_staircase_entropy_four_fifths_thm3502.py
```

The companion contains no assertion-dependent gates.  It checks:

1. the exact physical `S_23` orbit and the counts in (13);
2. an independent serial-stage Mealy transition against the packed two-mask
   transition on all `5,320,356` directed reachable state/input pairs;
3. prefix/suffix construction, zero-one adjacency, the deterministic recurrence
   (15), every row of (17), certificate statistics, and payload digest;
4. the full cleared-integer prefactor (21), entropy inequality (22), marker
   inequality (28), and floor-sensitive word-budget samples;
5. the complete `S_24` orbit, the physical prefix and suffix but nonphysical
   whole word in (35)--(37); and
6. all 126 finite block-code witnesses in (39).

Ordinary and optimized transcripts match the stored output byte for byte.  The
independent hostile audit rederived the factor-path injection, checked every
inequality direction and the charged `h,b` budget, and retained the
spurious-path and fixed-state-permutation controls.

No Rule 30 prize, literature novelty, center
density, center nonperiodicity, fixed-seed lower bound, or LRC statement is
claimed.
