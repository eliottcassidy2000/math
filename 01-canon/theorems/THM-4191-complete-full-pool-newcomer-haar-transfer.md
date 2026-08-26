---
id: THM-4191
title: "Complete full-pool newcomer Haar transfer"
status: >
  PROVED RELATIVE TO THM-4150/4188 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. For every positive newcomer q outside the fixed thirty-label pool,
  all C(30,10)=30,045,015 ten-subsets of the full pool have complete
  1/14-safe-set Haar mass at least 4/63 after q is adjoined. The q=50
  depth-seven repair deck has no ten-cover, and the native depth-six deck has
  no ten-cover at each of THM-4188's twenty-three exceptional newcomer labels.
  A literal q=50 depth-six ten-cover proves that the seventh deletion is sharp
  for this certificate. THM-4150 transfers every body to every distinct
  positive odd-tail pair after doubling. This subsumes the body-safety slices
  of THM-4174/4175/4179/4188 but not their structural filtrations, arbitrary
  body entry, or LRC(14).
source: lrc14-incoming-breakthrough-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4188-all-newcomer-zero-original-anchor-hierarchy-and-resonance-filtration
related:
  - THM-4174-six-deletion-completion-of-divisor-complete-newcomer-haar-transfer
  - THM-4175-haar-failure-atom-deletion-tomography-and-anchor-exchange
  - THM-4179-q50-seventh-deletion-primitive-anchor-completion
script: 04-computation/lrc14_complete_full_pool_newcomer_haar_transfer_thm4191.cpp
output: 05-knowledge/results/lrc14_complete_full_pool_newcomer_haar_transfer_thm4191.out
independent_audit_script: 04-computation/lrc14_complete_full_pool_newcomer_haar_transfer_independent_audit_thm4191.cpp
independent_audit_output: 05-knowledge/results/lrc14_complete_full_pool_newcomer_haar_transfer_independent_audit_thm4191.out
script_sha256: ac5a7b13396c1ca3246e3187efbcd8c0acf4334807dd4ae3081f37648a371edb
output_sha256: 859fc1ccceaab8fd58dbf4eb4da8a81b1619be5d1b602a8587b86bb896bea10b
independent_audit_script_sha256: 2da9a1d9e00b9e24b0ca95a8ea1bb219b6315135bc4a0ffe71c947202cec6688
independent_audit_output_sha256: 009f47594589a6c41d017496954511207e368a3d03425b705da76eeec1c1e7e2
hash_basis: raw LF bytes
primary_audit: >
  PASS. The THM-4188 fixed-pool prefix-cell integrator constructs q=50 E7
  and every native resonance E6 exactly. A hash-ordered edge scan and Gosper
  enumeration exhaust all 30,045,015 labelled ten-subsets in each of 24 rows,
  with zero threshold equalities and no cover. It also verifies the sharp q=50
  E6 hostile and locks the full semantic ledger.
independent_audit: >
  ACCEPT. The separate THM-4188 joint-wall geometry explicitly refines every
  newcomer wall with the fixed pool, then recursively enumerates all
  ten-subsets. It agrees on every edge count, equality count, body count,
  edge-check total, closest body, missed edge, hostile, and no-cover verdict.
  A fresh referee replay byte-matches both frozen outputs.
---

# THM-4191 -- complete full-pool newcomer Haar transfer

**PROVED RELATIVE TO THM-4150/4188 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED; LRC(14) REMAINS OPEN.**

## 1. Theorem and inheritance pass

Retain the fixed pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                        (1)
```

For a finite positive set `S`, put

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14}.                (2)
```

> **Theorem.** For every positive integer `q notin P` and every
> `B in binom(P,10)`,
>
> ```text
> mu(G_(B union {q}))>=4/63.                            (3)
> ```
>
> Consequently, for every positive integer `c` and every two distinct
> positive odd integers `a,b`, some `x in R/Z` satisfies
>
> ```text
> min_(v in 2c(B union {q}) union {a,b})||vx||>=1/14.   (4)
> ```

The closest proved mechanism is THM-4188: it transported only the
`1,491,665` primitive divisor-complete zero-original bodies, using a
`32/297/24` hierarchy of small good anchors. THM-4174/4175/4179 separately
handled selected original-anchor-count slices. The present theorem returns to
the discarded object: the complete labelled global repair deck itself. Its
transversal number is already larger than ten before any anchor or arithmetic
filter is imposed.

The canonical hostile is an actual q=50 depth-six ten-cover, not a body known
to be unsafe. The corrected near miss is to interpret a certificate failure
as danger. The least-used sidecar is the full unprojected deletion edge. The
selected method card is “correct the object before sharpening the technique”:
test all labelled bodies against the lawful deck before restricting to the
candidate class that motivated the deck.

## 2. Global repair hypergraph

For `d in {6,7}`, retain THM-4188's exact repair layer

```text
E_d(q)={R in binom(P,d):
        mu(G_((P union {q})\R))>=4/63}.                 (5)
```

The newcomer is never deleted. Regard `E_d(q)` as a labelled hypergraph on
the thirty vertices of `P`. A ten-body `B` is a transversal precisely when
it meets every lawful repair edge. If it is not a transversal, choose
`R in E_d(q)` with `R intersect B=empty`. Then

```text
B union {q} subset (P union {q})\R,
G_((P union {q})\R) subset G_(B union {q}).             (6)
```

Safe-set monotonicity and `(5)` prove `(3)`. No arithmetic property of `B`
is used. The exact connection contract is

```text
source:       lawful global repair edges R subset P
target:       every ten-subset B subset P
observable:   labelled incidence B intersects R
preserved:    q, repair mass, full edge, threshold
destroyed:    wall addresses and the chosen safe phase after integration
sidecar:      the full failure-atom mass used to certify R
decisive test: tau(E_d(q))>10.                           (7)
```

## 3. Base deck and all nonresonant newcomers

At `q=50`, exact pool-wall integration gives

```text
|E_7(50)|=821,737,                 threshold equalities=0. (8)
```

Both audits enumerate the complete universe

```text
|binom(P,10)|=binom(30,10)=30,045,015                  (9)
```

and find a disjoint E7 edge for every body. Any transversal with fewer than
ten vertices extends inside the thirty-vertex ground set `P` to a ten-vertex
transversal. Thus absence of a ten-vertex transversal also rules out every
transversal of size at most ten. Therefore

```text
tau(E_7(50))>10.                                      (10)
```

The primary scan performs `896,394,953` ordered edge checks. Its closest
body before a miss is

```text
B={80,85,88,95,143,145,168,193,240,252},
R={8,16,42,63,132,170,290}.                            (11)
```

THM-4188 proves the exact inclusion law

```text
E_7(50) subset E_7(q)       iff       q notin Q_7,     (12)

Q_7={6,22,24,25,48,70,72,96,100,105,110,128,130,140,
     186,192,206,210,220,256,260,294,366}.             (13)
```

For `q notin Q_7`, the q=50 edge missed by `B` remains a lawful q-repair.
Equations `(6)` and `(10)--(12)` prove `(3)` for every nonresonant newcomer,
including the complete cofinal tail.

## 4. All 23 resonances

For every `q in Q_7`, both audits independently construct the native
depth-six deck and exhaust all `30,045,015` ten-bodies. The exact edge counts
are

| `q` | `|E_6(q)|` | `q` | `|E_6(q)|` | `q` | `|E_6(q)|` |
|---:|---:|---:|---:|---:|---:|
| 6 | 389365 | 22 | 421289 | 24 | 390537 |
| 25 | 320668 | 48 | 475101 | 70 | 469924 |
| 72 | 418692 | 96 | 247722 | 100 | 284863 |
| 105 | 325808 | 110 | 365850 | 128 | 301247 |
| 130 | 539067 | 140 | 416495 | 186 | 432488 |
| 192 | 364769 | 206 | 451190 | 210 | 261505 |
| 220 | 479384 | 256 | 273944 | 260 | 422570 |
| 294 | 454315 | 366 | 422399 |  |  |

Every row has zero threshold equalities and

```text
tau(E_6(q))>10                         for q in Q_7.    (14)
```

Here the same extension-to-ten argument used after `(9)` converts the exact
ten-body exhaustion into the stated transversal inequality. Thus `(6)` proves
`(3)` at every exceptional newcomer. Together with
Section 3, this exhausts every positive `q notin P`.

## 5. Sharp certificate boundary

The seventh deletion at the base is necessary for this hypergraph proof. At
`q=50`, the exact ten-body

```text
C={80,85,88,95,145,168,193,240,252,286}               (15)
```

meets every one of the `85,324` edges of `E_6(50)`. Hence

```text
tau(E_6(50))<=10.                                      (16)
```

Both implementations verify `(15)` against the complete labelled deck. The
first failed implication is

```text
no depth-six certificate  !=>  unsafe body.            (17)
```

Indeed Section 3 proves that the same body misses a lawful depth-seven edge.
The strongest survivor is exact: depth seven closes the base, while native
depth six closes all twenty-three newcomer labels where base-edge transport
fails and a native repair is needed here.

## 6. Transfer and complete family count

Equation `(3)` satisfies THM-4150's only body hypothesis. That theorem accepts
every nonempty finite positive body; primitivity and divisor-completeness are
not required. It gives `(4)`, and Haar invariance under multiplication by `c`
supplies every positive common content.

The original-label decomposition is

```text
zero of {120,126,143}:       C(27,10)   =  8,436,285,
exactly one:               3*C(27,9)   = 14,060,475,
exactly two:               3*C(27,8)   =  6,660,225,
all three:                   C(27,7)   =    888,030,
                                          ----------
                                           30,045,015. (18)
```

Thus every newcomer supplies all `30,045,015` scale-one ten-body cores, not
only the earlier selected slices. For `N>=290`, the outsider labels through
`N` supply exactly

```text
(N-30)*30,045,015.                                     (19)
```

The earlier theorems remain useful: they identify primitive/divisor-complete
subfamilies, minimal anchors, exact depth staircases, and deformation
resonances. This theorem supersedes only their restricted **body-safety
counts**, not that structure.

## 7. Two-path exact audit

The primary path reuses THM-4188's fixed-pool prefix-cell integration, orders
edges by a fixed SplitMix64 key, and enumerates all ten-subsets with Gosper
masks. It locks the complete semantic ledger at

```text
FNV1A64_LE=53664fb38b90d1c0.                           (20)
```

The independent path reuses THM-4188's separate explicit joint-wall
geometry, uses the same fixed edge order only to make first-miss diagnostics
comparable, and recursively generates all ten-subsets. Its independently
encoded ledger is

```text
FNV1A64_LE=e9ae70b15ef96dfb.                           (21)
```

The paths agree on every load-bearing edge count, equality count, body count,
no-cover verdict, edge-check total, closest body, missed edge, and the q=50
hostile. A separate referee replay of both final sources exits cleanly and
byte-matches the frozen outputs.

## 8. Scope firewall

This is a complete body theorem only inside one fixed thirty-label pool and
with exactly one outside-pool newcomer in the eleven-speed body. It does not
put a hypothetical LRC counterexample into that form, handle two or more
uncontrolled body labels, close other parity/entry branches, or prove
LRC(14). It also does not assert that depth seven is globally minimal among
all conceivable certificates; `(16)` proves sharpness only for this nested
repair-deck mechanism.

## 9. Exact artifacts and replay

Primary:

```text
04-computation/lrc14_complete_full_pool_newcomer_haar_transfer_thm4191.cpp
sha256 ac5a7b13396c1ca3246e3187efbcd8c0acf4334807dd4ae3081f37648a371edb

05-knowledge/results/lrc14_complete_full_pool_newcomer_haar_transfer_thm4191.out
sha256 859fc1ccceaab8fd58dbf4eb4da8a81b1619be5d1b602a8587b86bb896bea10b
```

Independent:

```text
04-computation/lrc14_complete_full_pool_newcomer_haar_transfer_independent_audit_thm4191.cpp
sha256 2da9a1d9e00b9e24b0ca95a8ea1bb219b6315135bc4a0ffe71c947202cec6688

05-knowledge/results/lrc14_complete_full_pool_newcomer_haar_transfer_independent_audit_thm4191.out
sha256 009f47594589a6c41d017496954511207e368a3d03425b705da76eeec1c1e7e2
```

Replay:

```bash
g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_complete_full_pool_newcomer_haar_transfer_thm4191.cpp \
  -o /tmp/lrc4191-full-primary
/tmp/lrc4191-full-primary | diff -u \
  05-knowledge/results/lrc14_complete_full_pool_newcomer_haar_transfer_thm4191.out -

g++ -std=c++20 -O2 -DNDEBUG -pthread \
  -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_complete_full_pool_newcomer_haar_transfer_independent_audit_thm4191.cpp \
  -o /tmp/lrc4191-full-independent
/tmp/lrc4191-full-independent | diff -u \
  05-knowledge/results/lrc14_complete_full_pool_newcomer_haar_transfer_independent_audit_thm4191.out -
```
