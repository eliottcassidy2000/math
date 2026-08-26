---
id: THM-4201
title: "Exactly-one newcomer projected-repair depth stratification"
status: >
  PROVED RELATIVE TO THM-4150/4178/4179/4188 + VERIFIED-EXACT PRIMARY
  POOL-WALL AND INDEPENDENT JOINT-WALL AUDITS; LRC(14) OPEN. On the 69
  exactly-one-original resonance-anchor rows, native depth five has no
  seven-cover in 65 rows. The four exceptional A3 rows at q=25/105/210/256
  have exactly 2/1/3/4 seven-covers; native depth six has no seven-cover and
  exact projected transversal numbers 10/9/10/10. Off those resonances the
  q=50 depth-seven deck transports. This independently recovers THM-4190's
  1,071,961-body slice; THM-4191 already subsumes its body-safety consequence.
source: open-frontiers-incoming-20260826b / incoming-repair-depth extension
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4178-q50-divisor-complete-anchor-triple-exchange
  - THM-4179-q50-seventh-deletion-primitive-anchor-completion
  - THM-4188-all-newcomer-zero-original-anchor-hierarchy-and-resonance-filtration
related:
  - THM-4174-six-deletion-completion-of-divisor-complete-newcomer-haar-transfer
  - THM-4175-haar-failure-atom-deletion-tomography-and-anchor-exchange
  - THM-4190-every-newcomer-exactly-one-original-direct-body-completion
  - THM-4191-complete-full-pool-newcomer-haar-transfer
script: 04-computation/lrc14_exactly_one_newcomer_projected_repair_depth_stratification_thm4201.cpp
output: 05-knowledge/results/lrc14_exactly_one_newcomer_projected_repair_depth_stratification_thm4201.out
independent_audit_script: 04-computation/lrc14_exactly_one_newcomer_projected_repair_depth_stratification_independent_audit_thm4201.cpp
independent_audit_output: 05-knowledge/results/lrc14_exactly_one_newcomer_projected_repair_depth_stratification_independent_audit_thm4201.out
script_sha256: d431de4833211f563659931b324535ea75075f93da54305ba17b7647acc07cff
output_sha256: 788a6e922c4e31b388537f2c3b7d7f360ba87cfa8a01a8714bb5cfc48946f62d
primary_dependency_sha256: 25a6978484c7ab122fdc6c8e1593cfa2ad3468f7184a156045ea7e6cb2efc45d
independent_audit_script_sha256: 4dc9c2281eb58dac93b44cfef870b638dd17c907028946ced4ff95d62aa3bcdf
independent_audit_output_sha256: 3d7a5fedf37b4a6e1e80be87fa6677f8d7d27c088272cbf8a952562a0272b902
independent_dependency_sha256: 58817d5f5e1a8cc07384f3ea82a1feb221af37ab0907afde890ab4fbdd949137
hash_basis: raw LF bytes
primary_audit: >
  PASS. The THM-4188 pool-wall/safe-prefix path forms inclusion-minimal
  projected decks and recursively rejects every cover through budget seven
  on all load-bearing rows. O2, O0, and Apple-Clang-17 UBSan streams match
  after LF normalization; semantic ledger 1b82498b55f45bd3.
independent_audit: >
  PASS. The separate THM-4188 joint-wall path retains every distinct projected
  edge and literally scans all 480,700 seven-sets per row. It independently
  recovers the ten depth-five blockers and depth-six transversal numbers.
  O2, O0, and Apple-Clang-17 UBSan streams match after LF normalization;
  semantic ledger
  dd68c95b98fb9926.
---

# THM-4201 -- exactly-one newcomer projected-repair depth stratification

**PROVED RELATIVE TO THM-4150/4178/4179/4188 + VERIFIED-EXACT IN TWO
INDEPENDENT WALL IMPLEMENTATIONS; LRC(14) REMAINS OPEN.** The new content is
the exact projected repair-depth stratification. Its body-safety corollary
independently recovers THM-4190's slice and is subsumed by THM-4191's stronger
full-pool theorem. It does not prove physical entry or arbitrary-body entry.

## 1. Statement and exact target universe

Retain

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290},
A_0={120,126,143},
A_1={40,143,252},
A_2={80,143,252},
A_3={143,240,252}.
```

For `i in {1,2,3}` define the exact optional ground

```text
V_i=P\(A_i union {120,126}),             |V_i|=25.       (1)
```

For every `K in binom(V_i,7)` put `B_(i,K)=A_i union K`.  Then
`B_(i,K)` is a ten-set, contains exactly the one original label `143`, is
primitive, and contains a multiple of every integer `2,...,14`.

> **Theorem.** For every positive integer `q notin P`, every
> `i in {1,2,3}`, and every `K in binom(V_i,7)`,
>
> ```text
> mu(G_(B_(i,K) union {q})) >= 4/63.                    (2)
> ```
>
> Consequently, for every positive integer `c` and every two distinct
> positive odd integers `a,b`, some `x in R/Z` satisfies
>
> ```text
> min_(v in 2c(B_(i,K) union {q}) union {a,b}) ||vx||
>     >= 1/14.                                          (3)
> ```

The union over the three anchor presentations consists exactly of the
ten-subsets of `P` which contain `143` and `252`, contain neither `120` nor
`126`, and contain at least one of `{40,80,240}`.  Inclusion--exclusion gives

```text
3 binom(25,7)-3 binom(24,6)+binom(23,5)=1,071,961.      (4)
```

Thus the certificate promotes all `1,071,961` THM-4179 exactly-one-original
cores from `q=50` to every newcomer. This is an independent structural proof
of a subfamily already covered by THM-4191's stronger full-pool theorem

```text
binom(30,10)=30,045,015.                                (5)
```

The `1,071,961` bodies are therefore not an additional body count.

### Type-critical role of `V_i`

Removing `{120,126}` from the optional ground in `(1)` is both a counting
restriction and a load-bearing proof restriction.  The theorem claims only
the exactly-one-original slice, so its body choice `K` never uses `120` or
`126`. Repair deletions are *not* forbidden from using those two labels.
Instead, their
incidence with a possible `K` is forgotten by projection.

In particular, this theorem does **not** claim an all-newcomer result for
every `K in binom(P\A_i,7)`.  That stronger 27-vertex statement was proved at
`q=50` by THM-4179 but was not audited here away from `q=50`.

## 2. Inheritance pass and exact projected repair object

The closest proved mechanisms are THM-4179's `q=50` depth-seven completion
and THM-4188's exact all-newcomer edge-inclusion filtration.  The canonical
hostiles are THM-4179's depth-six seven-blockers and the four native
depth-five rows below.  The corrected near miss is that a shallow cover is a
body obstruction; THM-4179 maximal-deletion duality and the exact positive
body margins below refute that implication.  The least-used sidecar is the
literal 25-vertex body-choice ground, retained through repair projection.

Use THM-4188's global repair layer

```text
E_d(q)={R in binom(P,d):mu(G_((P union {q})\R))>=4/63}.
```

For an exactly-one anchor define the projected restricted hypergraph

```text
H_d(i,q)={R intersect V_i:R in E_d(q), R intersect A_i=empty}. (6)
```

For any `K subset V_i`,

```text
K hits every edge of H_d(i,q)
iff
K hits every R in E_d(q) disjoint from A_i.             (7)
```

Indeed, `K subset V_i`, so `K intersect R=K intersect (R intersect V_i)`.
Projection may discard `120,126` and any other non-`V_i` label, but it loses
no incidence observable to the declared body choices.  An empty projected
edge, if one occurred, would itself be a uniform disjoint repair rather than
a cover constraint.

Let

```text
c_7(H)=#{K in binom(V_i,7):K is a transversal of H}.     (8)
```

This exact-seven statistic is the direct proof obligation: every target body
has exactly seven optional labels.  Since `|V_i|=25`, `c_7(H)=0` also rules
out every transversal of size at most seven, because any smaller transversal
extends to a seven-set.

## 3. Proof outside the 23 resonances

THM-4188 proves for every positive `q notin P`

```text
E_7(50) subset E_7(q) iff q notin Q_7,                  (9)

Q_7={6,22,24,25,48,70,72,96,100,105,110,128,130,140,
     186,192,206,210,220,256,260,294,366}.              (10)
```

THM-4178 gives `tau(E_6^(A_1)(50))=8`.  For any seven-set
`K subset V_1`, its disjoint depth-six repair `R` can be extended by one label
outside `A_1 union K union R` to a depth-seven repair, because `6+7<27`; deletion
monotonicity preserves the threshold.  THM-4179 directly proves
`tau(E_7^(A_2)(50))=tau(E_7^(A_3)(50))=8`.  Hence every target `K` misses a
member of `E_7(50)` disjoint from its anchor.  If `q notin Q_7`, `(9)` keeps
that same labelled repair at `q`.  This proves `(2)` outside `Q_7` by safe-set
monotonicity.

## 4. Exact native closure on the 23 resonances

> **Certificate theorem.** Across the `23*3=69` resonance-anchor pairs,
> `c_7(H_5(i,q))=0` in exactly 65 rows. The remaining rows are
> `(A_3,q)` for `q=25,105,210,256`, with respectively `2,1,3,4` exact
> seven-covers. In all four rows `c_7(H_6(3,q))=0`, and the exact projected
> transversal numbers are `10,9,10,10`.

The independent joint-wall audit gives the following complete table.  The
middle triple is the number of distinct projected edges
`|H_5(1,q)|/|H_5(2,q)|/|H_5(3,q)|` after literal deduplication; the final
triple is the corresponding exact `c_7` count. The audits reproduce
THM-4188's `|E_5(q)|`; THM-4188 supplies the zero-equality ledger.

| `q` | `|E_5(q)|` | projected unique edges `A_1/A_2/A_3` | `c_7`, `A_1/A_2/A_3` |
|---:|---:|:---|:---|
| 6 | 50,160 | 21,728 / 18,647 / 16,208 | 0 / 0 / 0 |
| 22 | 54,396 | 27,213 / 23,730 / 22,003 | 0 / 0 / 0 |
| 24 | 51,188 | 23,419 / 19,593 / 18,080 | 0 / 0 / 0 |
| 25 | 34,684 | 15,819 / 14,056 / 11,786 | 0 / 0 / 2 |
| 48 | 77,327 | 35,912 / 31,703 / 33,781 | 0 / 0 / 0 |
| 70 | 73,996 | 35,292 / 31,337 / 29,866 | 0 / 0 / 0 |
| 72 | 59,802 | 28,200 / 26,184 / 22,324 | 0 / 0 / 0 |
| 96 | 21,457 | 8,734 / 7,098 / 8,083 | 0 / 0 / 0 |
| 100 | 26,410 | 11,721 / 9,079 / 8,249 | 0 / 0 / 0 |
| 105 | 37,139 | 16,606 / 14,735 / 13,228 | 0 / 0 / 1 |
| 110 | 44,484 | 17,529 / 14,127 / 16,574 | 0 / 0 / 0 |
| 128 | 31,275 | 11,325 / 8,970 / 8,612 | 0 / 0 / 0 |
| 130 | 105,407 | 50,140 / 47,062 / 48,227 | 0 / 0 / 0 |
| 140 | 59,351 | 25,665 / 21,389 / 20,669 | 0 / 0 / 0 |
| 186 | 61,671 | 28,351 / 24,431 / 23,275 | 0 / 0 / 0 |
| 192 | 43,980 | 18,108 / 15,208 / 16,532 | 0 / 0 / 0 |
| 206 | 65,683 | 31,832 / 27,009 / 28,588 | 0 / 0 / 0 |
| 210 | 26,533 | 9,604 / 8,515 / 6,844 | 0 / 0 / 3 |
| 220 | 77,964 | 35,132 / 30,779 / 30,776 | 0 / 0 / 0 |
| 256 | 25,765 | 10,176 / 8,015 / 7,164 | 0 / 0 / 4 |
| 260 | 61,262 | 26,322 / 23,085 / 22,240 | 0 / 0 / 0 |
| 294 | 68,909 | 29,666 / 25,249 / 25,007 | 0 / 0 / 0 |
| 366 | 57,331 | 25,835 / 21,550 / 19,839 | 0 / 0 / 0 |

Thus native depth five closes 65 of the 69 `(q,A_i)` rows.  Only four rows,
all on `A_3`, require escalation.  Their complete depth-six census is:

| `q` | `|E_6(q)|` | `|H_6(3,q)|` | `c_7(H_6(3,q))` | exact `tau(H_6(3,q))` |
|---:|---:|---:|---:|---:|
| 25 | 320,668 | 109,010 | 0 | 10 |
| 105 | 325,808 | 113,647 | 0 | 9 |
| 210 | 261,505 | 79,772 | 0 | 10 |
| 256 | 273,944 | 85,739 | 0 | 10 |

The four depth-five transversal numbers are exactly seven.  Their complete
seven-cover lists and direct body controls are:

| `q` | every depth-five seven-cover `K` | exact direct `63m-4D_q` |
|---:|:---|---:|
| 25 | `{85,88,95,132,145,168,193}` | 653,037,585,685,254 |
| 25 | `{85,88,95,145,168,193,264}` | 653,523,957,351,594 |
| 105 | `{85,88,95,132,145,168,193}` | 137,058,830,954,262 |
| 210 | `{16,85,88,95,145,168,193}` | 133,790,709,734,640 |
| 210 | `{85,88,95,132,145,168,193}` | 143,619,922,477,458 |
| 210 | `{85,88,95,145,168,193,264}` | 146,192,452,975,926 |
| 256 | `{16,85,88,95,132,168,193}` | 2,208,945,420,479,952 |
| 256 | `{16,85,88,95,168,193,264}` | 2,197,395,006,204,888 |
| 256 | `{16,88,95,132,168,170,193}` | 2,139,814,748,215,242 |
| 256 | `{16,88,95,168,170,193,264}` | 2,133,799,951,864,014 |

Here `D_q=lcm(D,14q)` is the joint denominator and `m` is the exact direct
joint-wall mass numerator of `G_(A_3 union K union {q})`.  Every displayed
delta is positive.  Hence these are failures of the depth-five certificate,
not unsafe bodies.  Depth six sharply raises the projected transversal number
to `10,9,10,10`.

For each resonance row and target `K`, select depth five unless the row is one
of the four exceptions, and depth six there.  The chosen projected deck has
`c_7=0`, so some repair `R` is disjoint from both `A_i` and `K`.  Therefore

```text
B_(i,K) union {q} subset (P union {q})\R,
mu(G_(B_(i,K) union {q}))
 >= mu(G_((P union {q})\R)) >=4/63.                    (11)
```

This proves `(2)` on `Q_7` and completes the all-newcomer claim.  THM-4150
and Haar invariance under multiplication by the positive content `c` prove
`(3)`.

## 5. Exact audit contract

The load-bearing exact universe is:

```text
23 newcomer resonances x 3 anchors x binom(25,7)
  =33,168,300 depth-five body candidates,

4 exceptional rows x binom(25,7)
  =1,922,800 depth-six body candidates.                (12)
```

The primary path includes THM-4188's pool-wall/safe-prefix geometry,
forms inclusion-minimal projected decks, and runs a memoized recursive cover
search through budget seven (and through the exact depth-six transversal
number on the four hostile rows).  Its semantic ledger is

```text
1b82498b55f45bd3.                                      (13)
```

The independent path includes THM-4188's explicit joint-wall refinement,
deduplicates projected edges without deleting supersets, literally enumerates
all `480,700` seven-sets in every row, and uses its separate recursive solver
only for the non-load-bearing exact depth-six `tau` sharpening.  Its semantic
ledger is

```text
dd68c95b98fb9926.                                      (14)
```

Both paths reproduce the THM-4188 native global edge counts.  All native
`E_5/E_6` threshold-equality counts on the 23 labels are zero, as are the
inherited `q=50` `E_5/E_6/E_7` equality counts.  Projection changes only the
incidence representation and introduces no new mass comparison.

Execution status on Windows MinGW GCC 15.2.0:

```text
primary O3: PASS, ledger 1b82498b55f45bd3
primary O2: PASS, ledger 1b82498b55f45bd3
primary O0: PASS, ledger 1b82498b55f45bd3
independent O2: PASS, ledger dd68c95b98fb9926
independent O0: PASS, ledger dd68c95b98fb9926
primary Apple Clang 17 UBSan: PASS, ledger 1b82498b55f45bd3
independent Apple Clang 17 UBSan: PASS, ledger dd68c95b98fb9926
```

The original Windows MinGW host could not link `-lubsan`; that environment
limitation is preserved as provenance, but the later macOS Apple-Clang-17
replays close the sanitizer audit for both implementations.

Replay from repository root:

```bash
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_exactly_one_newcomer_projected_repair_depth_stratification_thm4201.cpp \
  -o /tmp/lrc4201-primary-o2
/tmp/lrc4201-primary-o2

g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_exactly_one_newcomer_projected_repair_depth_stratification_independent_audit_thm4201.cpp \
  -o /tmp/lrc4201-independent-o2
/tmp/lrc4201-independent-o2

g++ -std=c++20 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_exactly_one_newcomer_projected_repair_depth_stratification_independent_audit_thm4201.cpp \
  -o /tmp/lrc4201-independent-o0
/tmp/lrc4201-independent-o0

clang++ -std=c++20 -O1 -g -fsanitize=undefined \
  -fno-sanitize-recover=undefined -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_exactly_one_newcomer_projected_repair_depth_stratification_thm4201.cpp \
  -o /tmp/lrc4201-primary-ubsan
/tmp/lrc4201-primary-ubsan

clang++ -std=c++20 -O1 -g -fsanitize=undefined \
  -fno-sanitize-recover=undefined -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_exactly_one_newcomer_projected_repair_depth_stratification_independent_audit_thm4201.cpp \
  -o /tmp/lrc4201-independent-ubsan
/tmp/lrc4201-independent-ubsan
```

The maintained wrappers locally suppress only the artificial `return-type`
diagnostic created by macro-renaming the complete inherited THM-4188 programs'
special `main` functions. All other warnings remain errors; the inherited
entry points are never called.

## 6. Connection contract and firewalls

```text
source:       THM-4179 q=50 exactly-one completion plus THM-4188 q-deformation
target:       exact shallow-certificate depth on the exactly-one slice
map:          q50 E7 inheritance off Q7; projected native E5/E6 decks on Q7;
              missed repair -> safe-set monotonicity -> THM-4150
preserved:    q, anchor labels, exact 4/63 threshold, 25-vertex body universe,
              divisor pins, primitivity, positive content, distinct odd tails
destroyed:    repair incidence on labels unavailable to K (especially 120,126),
              component order, canonical phase, selected repair after existence
sidecar:      anchor-labelled projection R -> R intersect V_i and full incidence
positive:     c7=0 on 65 native E5 rows and all four escalated E6 rows
hostile:      exactly 2/1/3/4 E5 seven-covers at q=25/105/210/256 on A3;
              all ten blocker bodies nevertheless have positive direct margin
decisive:     35,091,100 exact body candidates plus two wall geometries
```

Scope firewalls:

1. `q notin P` is essential for distinct newcomer labels and the core count.
2. The all-newcomer claim is exactly the `K subset V_i` slice.  It does not
   extend to body choices containing `120` or `126`.
3. The result is a sufficient Haar certificate for doubled bodies and two
   distinct positive odd tails.  It is not a necessary criterion.
4. The body-safety corollary is not a new disjoint count: THM-4190 proves the
   same slice directly and THM-4191 subsumes it for every ten-body in `P`.
5. Arbitrary bodies, physical entry, mixed/even tails, and LRC(14) remain open.
