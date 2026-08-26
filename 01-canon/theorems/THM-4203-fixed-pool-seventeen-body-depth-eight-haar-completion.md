---
id: THM-4203
title: "Fixed-pool seventeen-body depth-eight Haar completion"
status: >
  PROVED RELATIVE TO THM-4150/4188 + VERIFIED-EXACT PRIMARY POOL-WALL AND
  INDEPENDENT JOINT-WALL AUDITS; THIRD REFEREE MATCH; LRC(14) OPEN. Every
  nonempty subset of the fixed thirty-label pool of cardinality at most
  seventeen has complete 1/14-safe-set Haar mass at least 4/63. The q=50
  depth-eight repair hypergraph has exact transversal number fourteen. Its
  57,410 terminal
  covers through budget seventeen expand to exactly 66,468 transversals of
  sizes fourteen through seventeen, and every one has positive direct Haar
  margin. THM-4191 then gives every eleven-subset of every one-newcomer pool.
  Physical entry, two or more uncontrolled labels, and LRC(14) remain open.
source: lrc-full-pool-entry / open-frontiers-incoming-20260826b
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4188-all-newcomer-zero-original-anchor-hierarchy-and-resonance-filtration
  - THM-4191-complete-full-pool-newcomer-haar-transfer
related:
  - THM-4201-exactly-one-newcomer-projected-repair-depth-stratification
  - THM-4207-two-newcomer-sharp-depth-transition-base-surplus-composition-and-variable-pool-chart-number
  - THM-4211-fixed-fifty-cofinal-two-newcomer-haar-tail-and-eighteen-label-chart
  - THM-4214-two-newcomer-pascal-complete-eleven-body-haar-charts
script: 04-computation/lrc14_fixed_pool_seventeen_body_depth8_haar_completion_thm4203.cpp
output: 05-knowledge/results/lrc14_fixed_pool_seventeen_body_depth8_haar_completion_thm4203.out
independent_audit_script: 04-computation/lrc14_fixed_pool_seventeen_body_depth8_haar_completion_independent_audit_thm4203.cpp
independent_audit_output: 05-knowledge/results/lrc14_fixed_pool_seventeen_body_depth8_haar_completion_independent_audit_thm4203.out
script_sha256: 6bb1b294bc2e9c894eaa83764feeb3f262b36f74c5f524874b7cafbd0fb9d668
output_sha256: c5dca81099ce7f6c7a54d5733ff1cb77299c19afb86cd91b8fad9e86f8c27601
independent_audit_script_sha256: 0e69c7556f3689dbfde995f6055356c6ba45081ea66fb8962352b3dab148d486
independent_audit_output_sha256: 358d79793f344350283586b9a36dc8bde59343975f85e09d0b415b16f5d6a14a
hash_basis: raw LF bytes
primary_audit: >
  PASS. The inherited THM-4188 fixed-pool prefix-cell integrator builds q=50
  E7/E8 exactly, exhausts the budget-seventeen hitting-set tree, completes
  its terminal covers upward, checks every direct margin, and separately
  scans all C(30,11) eleven-bodies. Its semantic ledger is
  68ff828dae998c28.
independent_audit: >
  ACCEPT. A separate explicit q=50 joint-wall refinement rebuilds E7/E8,
  repeats the complete terminal and upward-closure censuses, checks direct
  margins in the independent fixed geometry, and agrees on every count,
  witness, minimum, and fingerprint. Its semantic ledger is
  d9e1bd77881a64a7.
---

# THM-4203 -- fixed-pool seventeen-body depth-eight Haar completion

**PROVED RELATIVE TO THM-4150/4188 + VERIFIED-EXACT IN TWO MAINTAINED WALL
IMPLEMENTATIONS AND A THIRD REFEREE PATH; LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance pass

Retain the fixed pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.
```

For a nonempty finite set `S` of positive integers, put

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14}.
```

> **Theorem.** For every nonempty `B subset P` with `|B|<=17`,
>
> ```text
> mu(G_B)>=4/63.                                         (1)
> ```
>
> Consequently, for every nonempty such `B`, every positive integer `c`, and
> every two distinct positive odd integers `a,b`, some `x in R/Z` satisfies
>
> ```text
> min_(v in 2cB union {a,b})||vx||>=1/14.               (2)
> ```

At `|B|=11`, `(2)` is the thirteen-nonzero-relative-speed form relevant to
LRC(14). The layers `|B|=12,...,17` give special configurations with fourteen
through nineteen nonzero relative speeds at the same `1/14` threshold; they
are not general lonely-runner statements at those cardinalities.

The closest proved mechanism is THM-4191, which transported every ten-pool
body after adjoining one arbitrary newcomer. THM-4201's exact depth escalation
is the incoming signal: a shallow repair-deck cover is a reason to test the
next deletion layer, not evidence that the covered body is unsafe. The
canonical hostile is the depth-seven eleven-cover in Section 4. The corrected
near miss is `certificate failure => unsafe body`, refuted here by direct
positive margins. The least-used inherited object is the unrestricted
depth-eight global repair deck.

## 2. Repair hypergraph and sharp transversal depth

For `d in {7,8}`, define the q=50 repair layer

```text
E_d(50)={R in binom(P,d):
         mu(G_((P union {50})\R))>=4/63}.                (3)
```

The newcomer `50` is retained; every repair deletion lies in `P`. The primary
pool-wall and independent joint-wall implementations agree exactly on

```text
|E_7(50)|=  821,737,     threshold equalities=0,
|E_8(50)|=4,088,855,     threshold equalities=0.         (4)
```

> **Exact depth-eight transversal theorem.**
>
> ```text
> tau(E_8(50))=14.                                      (5)
> ```

The budget-seventeen search branches on every vertex of the first uncovered
edge. Therefore every cover extension is represented on some branch; a
visited-mask table removes duplicate states without removing a distinct
chosen set. It visits

```text
nodes=4,991,151,             dead leaves=1,428,692,      (6)
```

and its smallest terminal cover has size fourteen. For the matching upper
bound, the explicit set

```text
C={8,10,16,80,85,88,95,143,145,168,176,193,240,252}     (7)
```

meets every E8 edge. Thus the exhaustive lower bound and `(7)` prove `(5)`.

The primary path also separately enumerates all

```text
binom(30,11)=54,627,300                                  (8)
```

eleven-subsets by Gosper masks and finds a disjoint E8 edge in every case,
using `4,506,326,556` literal edge checks. This scan is redundant for `(5)`
but is a useful independent control on the LRC(14)-relevant layer.

## 3. Complete bounded transversal census

The hitting-set recursion stops immediately when its chosen set becomes a
cover. Its output is therefore a census of terminal search covers within the
budget--not a census of all cover supersets. A terminal cover need not be
inclusion-minimal: a branch can choose a redundant vertex before it first
becomes a cover. This distinction is load-bearing because direct safe-set
mass is not upward monotone. Through size seventeen the terminal census is

```text
size 14:     35,
size 15:    817,
size 16:  8,296,
size 17: 48,262,
total:   57,410.                                         (9)
```

For every transversal `B` of cardinality at most seventeen, follow a search
branch by choosing at each node a vertex of `B` in the first uncovered edge.
Such a vertex always exists because `B` is a transversal. The branch reaches
a terminal cover `C subset B` within the budget. Conversely, every superset
of a terminal cover is a transversal. Hence a complete bounded census is
obtained by taking every one-, two-, and three-label extension needed to reach
the budget, sorting, and deduplicating. The result is

```text
size 14:     35,
size 15:    876,
size 16:  9,320,
size 17: 56,237,
total:   66,468.                                         (10)
```

Both maintained implementations compute the direct fixed-pool Haar margin of
every body in `(10)`. If `D=lcm_(p in P)(14p)` and
`m_B=D mu(G_B)`, then the exact minimum is

```text
min_B (63m_B-4D)=34,853,530,580,568>0,                   (11)
```

attained at

```text
{8,10,16,42,60,80,85,88,95,143,145,168,176,193,
 240,252,264}.                                           (12)
```

The sorted `(mask,direct margin)` FNV-1a/64 fingerprints are

```text
terminal covers:            48a6fc172ff16dda,
exact size-seventeen covers: ce55f5f78c168012,
all covers through seventeen: ed6302225a9454e8.          (13)
```

These shared fingerprints audit the full labelled universes, not merely their
counts and extrema.

## 4. Proof for every body through seventeen

Fix a nonempty `B subset P` with `|B|<=17`.

If `B` is not a transversal of `E_8(50)`, choose `R in E_8(50)` disjoint from
`B`. Then

```text
B subset (P union {50})\R,
G_((P union {50})\R) subset G_B.                         (14)
```

Definition `(3)` and safe-set monotonicity give `(1)`. If `B` is a
transversal, it is one of the complete `66,468` bodies in `(10)`, and its
direct margin is positive by `(11)`. This exhausts the dichotomy and proves
`(1)`. THM-4150 applies to the complete safe set and proves `(2)`; Haar
invariance supplies every positive common content `c`.

The theorem covers `879,612,196` nonempty labelled subsets of `P` through
size seventeen. In particular, all `54,627,300` eleven-subsets of `P` form
the zero-outsider input needed by later outsider-count decompositions.

## 5. Sharp certificate hostiles and depth-nine anatomy

Depth seven cannot replace depth eight. The eleven-set

```text
W={8,80,85,88,95,143,145,168,193,240,252}               (15)
```

meets all `821,737` E7 edges. Nevertheless its exact direct margin is

```text
63m_W-4D=139,755,964,331,592>0.                          (16)
```

It is therefore safely above `4/63`; `(15)` is a certificate failure, not an
unsafe body. At depth eight the disjoint repair

```text
R={16,84,120,176,190,264,286,290}                        (17)
```

has strictly positive repair margin and breaks the cover. The exhaustive
eleven-body control identifies `(15)` as its closest body and `(17)` as its
first missed edge.

The exact depth-eight cover `(7)` also has positive direct body margin

```text
63m_C-4D=77,740,560,077,796>0.                           (18)
```

Thus `tau(E8)=14` is a sharp boundary of this repair certificate, not a
safety boundary. A finite-exact addendum makes the distinction concrete:
among the `binom(16,9)=11,440` depth-nine deletions disjoint from `C`, exactly
`108` are lawful E9 repairs and there are zero threshold equalities. This
repairs the displayed E8 cover at the next depth; it does not prove any global
E9 transversal statement and is not load-bearing for the theorem.

## 6. One-newcomer corollary and physical-entry obstruction

Let `q notin P`. Every eleven-subset of `P union {q}` has exactly one of the
following two forms:

```text
H subset P, |H|=11:             binom(30,11)=54,627,300,
H=B union {q}, B subset P,
|B|=10:                         binom(30,10)=30,045,015.
```

The first type follows from this theorem and the second from THM-4191. Hence,
jointly, the two theorems close every one of

```text
binom(31,11)=54,627,300+30,045,015=84,672,315            (19)
```

eleven-bodies in every one-newcomer pool `P union {q}`.

For a physical parity row with exactly eleven even speeds and two odd speeds,
write the even block canonically as `2cH`, where `2c` is its common even
content and `H` is primitive. There is no discretionary rescaling after this
normalization: the body tested by the Haar theorem is literally `H`.
Consequently, `(1)` and `(19)` imply the necessary obstruction

```text
counterexample in this parity branch  =>  |H\P|>=2.      (20)
```

In THM-3818's exact scale-two `11+2` branch, the unresolved odd-pair lane has
`c=1` and `H=u`, so `(20)` is literal. Scale-one rows and other component
shapes do not enter this carrier.

## 7. Audit contract, scope firewalls, and generated tasks

The primary path inherits THM-4188's exact 7,133-cell pool partition and
integrates the q=50 comb by safe prefixes. The independent path constructs
the explicit 7,213-cell q=50 joint-wall refinement and rebuilds the low-arity
atom support. They independently reproduce E7/E8, the terminal cover search,
the complete upward closure, all direct margins, and the three fingerprints
in `(13)`. Both paths include positive and negative synthetic hitting-set
controls. A third support-at-most-eight joint-wall referee independently
matched all counts, extrema, and fingerprints and hard-replayed its frozen
output. Its semantic ledger is `3799dfd0142d4014`; the scratch source/output
SHA-256 pair is
`d94f1c1c7a59f29ee5a69dae56b6131b5b214a81477b5325cf929ce893948841` /
`6d68ad09d66e431679e27e87776d5a179833712b6cdbdf5d8d009dd95b172a17`.
This referee is corroborative and not a maintained dependency.

Scope firewalls:

1. `(1)` covers subsets of this one fixed pool only.
2. `(19)` covers exactly eleven-body subsets of a one-newcomer pool; it does
   not cover larger newcomer bodies or two uncontrolled labels.
3. Repair hypergraphs are sufficient certificates, not characterizations of
   body safety; `(16)` and `(18)` are explicit hostiles.
4. Physical entry into the `11+2` parity branch, scale-one rows, arbitrary
   parity, and LRC(14) remain open.
5. The E9 addendum repairs one E8 cover and proves no global E9 statement.

The next finite boundary task is the exact size-eighteen upward completion:
enumerate terminal E8 covers to budget eighteen, generate every exact
size-eighteen transversal from a terminal search subcover, and test whether the
minimum direct margin remains positive. This is a generated experiment, not a
claim of the theorem.

The cheapest precise outside-pool extension is also finite. For

```text
E_8(q)={R in binom(P,8):mu(G_((P union {q})\R))>=4/63},
Q_8={q notin P:E_8(50) is not a subset of E_8(q)},       (21)
```

extend THM-4188's incidence-inclusion/cofinal-discrepancy engine to determine
`Q8`. Universal inclusion closes every one-newcomer body `B union {q}` with
`B subset P`, `|B|<=13`, because `(5)` supplies a disjoint E8 edge. Any finite
exception set can then be tested by its native E8 transversal number. This
does not yet use the direct completion of E8 transversals, because adjoining
`q` can erode their fixed-pool direct margins.

## 8. Replay

Primary:

```bash
g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed_pool_seventeen_body_depth8_haar_completion_thm4203.cpp \
  -o /tmp/lrc4203-primary
/tmp/lrc4203-primary | diff -u \
  05-knowledge/results/lrc14_fixed_pool_seventeen_body_depth8_haar_completion_thm4203.out -
```

Independent:

```bash
g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_fixed_pool_seventeen_body_depth8_haar_completion_independent_audit_thm4203.cpp \
  -o /tmp/lrc4203-independent
/tmp/lrc4203-independent | diff -u \
  05-knowledge/results/lrc14_fixed_pool_seventeen_body_depth8_haar_completion_independent_audit_thm4203.out -
```
