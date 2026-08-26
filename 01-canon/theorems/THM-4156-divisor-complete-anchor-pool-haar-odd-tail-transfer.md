---
id: THM-4156
title: "Divisor-complete anchor-pool Haar odd-tail transfer"
status: >
  PROVED RELATIVE TO THM-4150 + VERIFIED-EXACT + INDEPENDENT EXACT INTERVAL
  AUDIT; LRC(14) OPEN. A 30-label common safe-set pool has Haar measure
  strictly above 4/63. Requiring the three anchors 120,126,143 makes every
  selected eleven-body primitive and divisor-complete through 14. THM-4150
  therefore closes 2,220,075 bodies with every distinct positive odd tail
  pair, including 1,875,709 outside both named THM-4148/4151 gates.
source: codex-lrc14-planar-jc-breakthrough-20260825
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
related:
  - THM-366-lrc-small-denominator-divisibility-sieve
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-4075-lrc14-divisor-complete-dyadic-owner-word-closure-through-30
  - THM-4148-first-window-width-universal-odd-tail-lrc14-transfer
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4151-scale-sensitive-first-window-odd-tail-lrc14-transfer
  - THM-4154-mod-six-fixed-clock-and-haar-pool-inheritance-correction
script: 04-computation/lrc14_divisor_complete_anchor_pool_haar_transfer_thm4156.py
output: 05-knowledge/results/lrc14_divisor_complete_anchor_pool_haar_transfer_thm4156.out
independent_audit_script: 04-computation/lrc14_divisor_complete_anchor_pool_haar_transfer_thm4156_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_divisor_complete_anchor_pool_haar_transfer_thm4156_independent_audit.out
script_sha256: db264924c1234323ecb299fbe7aacc4bc4d3b46290ec754bf4ba23ddffa8103d
output_sha256: c30bfc3f0fdd3f2dfd37efb490f8e28fc36c30583a442361d2acce57f790aae1
semantic_sha256: ed193657974754e64ae4faed03bc63146f453c16290553ef5d88f198bc77bd0e
independent_audit_script_sha256: e94cec1a54f4ff36f32a66a13ad1c5b9263fe24c1b6a6fa164bc0b39eb8ff2ad
independent_audit_output_sha256: c0b2bdda38d7c08d4a3b060499c99a2ff96c4b9a4cab5c368a47cfeca80ecf1e
independent_semantic_sha256: 97954c5f353ca0977c87d1649abb646a43af4cdb636bcc211861a12129b1c415
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact Fraction arithmetic reconstructs every safe-set wall and cell,
  the divisor-owner table, reflection symmetry, a direct positive control,
  the failed one-twelfth hostile, and all 2,220,075 declared bodies. Normal,
  optimized, and hash-seeded outputs byte-match.
independent_audit: >
  ACCEPT. A no-import implementation intersects each speed's closed safe
  tooth union successively, using the global wall set only as a checksum. It
  reproduces all 7,134 walls, 150 components, the exact measure and maximum
  width, both controls, and an independently grouped family census. Normal,
  optimized, and hash-seeded outputs byte-match.
---

# THM-4156 -- a divisor-complete common-Haar pool

**PROVED RELATIVE TO THM-4150 + VERIFIED-EXACT + INDEPENDENT EXACT INTERVAL
AUDIT; LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance

Put

```text
A={120,126,143},

P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                       (1)
```

For every eight-element set `K subset P\A`, define the eleven-body

```text
H=A union K.                                           (2)
```

> **Theorem.** For every body `(2)` and every two distinct
> positive odd integers `a,b`, there exists `x in R/Z` such that
>
> ```text
> min_(v in 2H union {a,b})||vx||>=1/14.               (3)
> ```

The closest mechanism is THM-4150's sharp full-safe-set Haar transfer. The
canonical hostile is MISTAKE-511: a large pool may look new only because its
labels all evade one small divisor. Here that hostile is forced to fail:
every body contains multiples of `6`, and the old clock `x=1/12` has zero
clearance at doubled speed `240`. The corrected near miss is precisely
THM-4154's three nonzero-mod-six pools. The least-used sidecar is no longer a
residue avoidance condition but a **mandatory anchor set**, retained through
the hereditary choice rather than forgotten in an all-subsets count.

## 2. The three anchors enter the divisor-complete seam

The anchor divisibilities can be frozen with one labelled owner per modulus:

| modulus | `2,3,4,5,6,8,10,12` | `7,9,14` | `11,13` |
|:---:|:---:|:---:|:---:|
| owner | `120` | `126` | `143` |

Thus every `H` in `(2)` contains a multiple of every integer `2,...,14`.
Moreover

```text
gcd(120,126,143)=1,                                    (4)
```

so every body is primitive. The anchor set is inclusion-minimal for this
particular cover: deleting `120`, `126`, or `143` respectively loses at least
one of the displayed modulus groups.

This is exactly the necessary divisor-complete quotient seam of THM-2061,
not THM-4154's already settled complement. It also lies beyond THM-4075's
height-30 finite closure because every body contains `143`.

## 3. Exact common safe-set measure

For a finite positive set `U`, write

```text
G_U={y in R/Z:min_(u in U)||uy||>=1/14}.                (5)
```

The complete rational wall arrangement for `G_P` has

```text
7,134 distinct walls,             150 components.      (6)
```

Exact cell summation gives

```text
mu(G_P)=298133356159/4560289854120
       =4/63+8591143199/4560289854120
       >4/63.                                          (7)
```

The largest component has length `37/25520`, and the component list is
invariant under `y -> 1-y`. These statistics are sidecars, not inputs to the
transfer: only the strict inequality in `(7)` is used below.

The primary certificate forms every wall

```text
(k+1/14)/p,                 (k+1-1/14)/p               (8)
```

for `p in P`, classifies each intervening open cell by an exact midpoint,
checks both closed endpoints of every retained cell, merges adjacent cells,
and sums their rational lengths. Hence `(6)--(7)` are finite-exact statements
over the explicitly printed universe, not sampled estimates.

As a positive control, the first retained component contains

```text
y=199/21280,             x=(y+1)/2=21479/42560.        (9)
```

At `(9)` the whole doubled 30-label pool has clearance `199/2660`, while
tails `1,3` have minimum clearance `20683/42560`. Conversely, the fixed
clock from THM-4154 fails literally:

```text
||2*120/12||=0.                                       (10)
```

## 4. Haar transfer and family count

For every body `(2)`, inclusion reverses the safe-set constraints:

```text
G_P subset G_H,                 mu(G_H)>=mu(G_P)>4/63. (11)
```

THM-4150 applied to `(11)` proves `(3)` for every distinct positive odd pair.
The even body speeds are pairwise distinct, the odd tails are pairwise
distinct and cannot meet an even speed, so `(3)` is a genuine thirteen-speed
row. The structurally distinct interval audit independently reproduces
`(6)--(7)`. This completes the independent verification. **QED.**

There are `27` optional labels, so the exact family size is

```text
binom(27,8)=2,220,075.                                  (12)
```

Every body has `max(H)>=143`. Therefore THM-4148's gate fails identically:
if `m=min(H)` and `M=max(H)`, its left-minus-right expression is

```text
27(13m-M)-4mM=m(351-4M)-27M<0.                         (13)
```

The complete min/max census for THM-4151's affine gate

```text
16M<=156m+13                                           (14)
```

finds exactly `344,366` admitted bodies. Consequently

```text
2,220,075-344,366=1,875,709                            (15)
```

fail both named gates. Equation `(15)` is a comparison with those two
declared criteria only; it is not a claim that no other
canonical certificate can close an individual member.

## 5. What changed and what remains open

The source-to-target contract is

```text
source:       a common closed safe set for P plus mandatory anchors A
target:       2(A union K) plus every distinct positive odd tail pair
map:          G_P -> G_(A union K), then THM-4150's two physical lifts
preserved:    closed safe phases, exact Haar mass, odd-tail sheet structure
destroyed:    individual component addresses and body-specific extra phases
sidecar:      the labelled anchor divisor cover through 14
positive:     exact strict 4/63 surplus in (7)
hostile:      x=1/12 is killed by the mandatory doubled anchor 120
decisive test: independent interval reconstruction of the same measure. (16)
```

The pool was found by a search, but no maximality, optimality, or all-pool
heredity is claimed. The theorem requires all three anchors and exactly eight
other labels only to state the counted eleven-body family; THM-4150 itself
works at other cardinalities. This result does not give physical entry into
the dyadic odd-tail seam, handle mixed/even tails, close arbitrary
divisor-complete bodies, or prove LRC(14).

Primary replay:

```bash
python3 -B 04-computation/lrc14_divisor_complete_anchor_pool_haar_transfer_thm4156.py
python3 -B -O 04-computation/lrc14_divisor_complete_anchor_pool_haar_transfer_thm4156.py
PYTHONHASHSEED=271828 python3 -B 04-computation/lrc14_divisor_complete_anchor_pool_haar_transfer_thm4156.py
```

Independent replay:

```bash
python3 -B 04-computation/lrc14_divisor_complete_anchor_pool_haar_transfer_thm4156_independent_audit.py
python3 -B -O 04-computation/lrc14_divisor_complete_anchor_pool_haar_transfer_thm4156_independent_audit.py
PYTHONHASHSEED=314159 python3 -B 04-computation/lrc14_divisor_complete_anchor_pool_haar_transfer_thm4156_independent_audit.py
```
