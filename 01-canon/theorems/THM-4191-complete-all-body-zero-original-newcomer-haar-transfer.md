---
id: THM-4191
title: "Complete all-body zero-original newcomer Haar transfer"
status: >
  PROVED RELATIVE TO THM-4150/4188 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. For every positive newcomer q outside the fixed thirty-label pool,
  every one of the C(27,10)=8,436,285 zero-original ten-bodies has complete
  1/14-safe-set Haar mass at least 4/63. The q=50 depth-seven projected repair
  deck has no ten-cover, and all twenty-three exceptional newcomer labels are
  closed by their native depth-six decks. A literal q=50 depth-six ten-cover
  proves that the extra deletion is genuinely needed by this certificate.
  THM-4150 transfers every certified body to every distinct positive odd-tail
  pair. This strengthens THM-4188's 1,491,665 primitive divisor-complete
  subfamily but does not prove arbitrary-body entry or LRC(14).
source: lrc14-incoming-breakthrough-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4188-all-newcomer-zero-original-anchor-hierarchy-and-resonance-filtration
related:
  - THM-4174-six-deletion-completion-of-divisor-complete-newcomer-haar-transfer
  - THM-4175-haar-failure-atom-deletion-tomography-and-anchor-exchange
  - THM-4179-q50-seventh-deletion-primitive-anchor-completion
script: 04-computation/lrc14_all_body_zero_original_newcomer_haar_transfer_thm4191.cpp
output: 05-knowledge/results/lrc14_all_body_zero_original_newcomer_haar_transfer_thm4191.out
independent_audit_script: 04-computation/lrc14_all_body_zero_original_newcomer_haar_transfer_independent_audit_thm4191.cpp
independent_audit_output: 05-knowledge/results/lrc14_all_body_zero_original_newcomer_haar_transfer_independent_audit_thm4191.out
script_sha256: 8bc1213c253cd450bd797af3e8e926dd69adc145af7eb85459b2b3acf0bbb7f2
output_sha256: 9f2ae05d3c362190533514f9fc919d8b324e7e3d64230b79ab40b086e8b5493c
independent_audit_script_sha256: ce956bbcbe3bdd940b294ad3c1cde7cf953677c77e92625c8f6494894886f0e2
independent_audit_output_sha256: 5b09a62c4746839d0ca7e98f0e4297be929cb64b4b4d0e7587748dc47320293e
hash_basis: raw LF bytes
primary_audit: >
  PASS. The THM-4188 fixed-pool prefix-cell integrator constructs q=50 E7
  and every native resonance E6 exactly. Gosper enumeration exhausts every
  labelled ten-subset of the twenty-seven-label ground set and records the
  first disjoint repair. It certifies 24 complete 8,436,285-body rows, zero
  threshold equalities, the sharp q=50 E6 hostile, and a locked semantic
  FNV-1a ledger.
independent_audit: >
  ACCEPT. The separate THM-4188 joint-wall geometry explicitly refines every
  newcomer wall with the fixed pool, then recursively enumerates ten-subsets.
  It agrees on every edge count, equality count, body count, cover verdict,
  edge-check total, closest body, and missed-edge diagnostic. Its independent
  semantic ledger is locked separately.
---

# THM-4191 -- complete all-body zero-original newcomer Haar transfer

**PROVED RELATIVE TO THM-4150/4188 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED; LRC(14) REMAINS OPEN.**

## 1. Theorem and inheritance pass

Retain THM-4188's pool and original labels

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290},
A_0={120,126,143},                         U=P\A_0.       (1)
```

For a finite positive set `S`, put

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14}.                (2)
```

> **Theorem.** For every positive integer `q notin P` and every
> `B in binom(U,10)`,
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

THM-4188 proved `(3)` only for the `1,491,665` primitive divisor-complete
members of `binom(U,10)`. Its closest mechanism was an inclusion-minimal good
anchor inside each such body. The discarded objects are the other
`6,944,620` labelled bodies: although they were irrelevant to that normalized
arithmetic slice, they remain legitimate vertices of the repair-transversal
problem and legitimate inputs to THM-4150.

The corrected near miss is to project away `A_0` **after** constructing a
lawful global repair. Constructing a repair only on `U` would forget whether
the deleted pool labels really certify Haar mass. The canonical hostile is
the exact q=50 depth-six cover in Section 4. The least-used sidecar is the
complete labelled deletion deck rather than its minimal-good-anchor quotient.

## 2. Projected repair duality

For `d in {6,7}`, retain the global repair layer

```text
E_d(q)={R in binom(P,d):
        mu(G_((P union {q})\R))>=4/63}.                 (5)
```

The newcomer is never deleted. Project a repair to the zero-original ground
set by

```text
pi(R)=R intersect U,                 H_d(q)=pi(E_d(q)). (6)
```

For every `B subset U` and `R subset P`, one has the exact labelled
equivalence

```text
B intersects pi(R)=empty  iff  B intersects R=empty.   (7)
```

Thus a ten-body blocks every lawful depth-`d` repair exactly when it is a
ten-vertex transversal of the projected deck `H_d(q)`. If it is not a
transversal, take `R in E_d(q)` disjoint from `B`. Then

```text
B union {q} subset (P union {q})\R,
G_((P union {q})\R) subset G_(B union {q}),             (8)
```

and safe-set monotonicity proves `(3)`. Projection destroys the identities of
deleted original labels only after `(5)` has preserved the actual repair
mass; the full edge remains the required sidecar.

## 3. Base deck and nonresonant newcomers

At `q=50`, exact pool-wall integration gives

```text
|E_7(50)|=821,737,                threshold equalities=0. (9)
```

Both audits enumerate the complete universe

```text
|binom(U,10)|=binom(27,10)=8,436,285                   (10)
```

and for every body find an edge of `E_7(50)` disjoint from it. Equivalently,

```text
tau(H_7(50))>10.                                       (11)
```

The primary path performs `359,519,024` ordered edge checks. Its closest row
before finding a missed edge is

```text
B={8,80,85,88,95,168,193,240,252,286},
R={10,15,120,126,143,145,290}.                          (12)
```

THM-4188 proves the exact all-newcomer inclusion law

```text
E_7(50) subset E_7(q)       iff       q notin Q_7,      (13)

Q_7={6,22,24,25,48,70,72,96,100,105,110,128,130,140,
     186,192,206,210,220,256,260,294,366}.              (14)
```

For `q notin Q_7`, the q=50 repair missed by `B` therefore remains a lawful
q-repair. Equations `(7)--(8)` prove `(3)` for every nonresonant newcomer,
including the complete cofinal tail; no new asymptotic estimate is needed.

## 4. The 23 native resonance rows and the sharp hostile

For every `q in Q_7`, both audits independently build the native depth-six
deck and exhaust all `8,436,285` bodies. The exact edge counts are

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

Every row has zero threshold equalities and no ten-cover:

```text
tau(H_6(q))>10                         for q in Q_7.     (15)
```

This closes `(3)` at every resonance by `(7)--(8)`. The largest ordered scan
for one body is `3,248` edges at `q=210`; it still ends with a disjoint lawful
repair.

The extra base deletion is not cosmetic. At `q=50`, the body

```text
{80,85,88,95,145,168,193,240,252,286}                 (16)
```

meets every one of the `85,324` edges of `E_6(50)`. Hence

```text
tau(H_6(50))<=10.                                      (17)
```

Both implementations verify `(16)` directly. The first failed implication
of a depth-six proof is a genuine labelled transversal, not low edge count or
an unsafe body. The strongest survivor is precisely the depth-seven theorem.

## 5. Transfer and exact family counts

Equation `(3)` meets THM-4150's only body hypothesis. That theorem applies to
every nonempty finite positive body; primitivity and divisor-completeness are
not required. It yields `(4)`, and Haar invariance under multiplication by
`c` supplies every positive common content.

Every newcomer has exactly

```text
binom(27,10)=8,436,285                                 (18)
```

distinct zero-original scale-one cores. THM-4188's `1,491,665` good bodies
remain the distinguished primitive divisor-complete subcount, but no longer
the safe-family count. For `N>=290`, the newcomer labels through `N` supply

```text
(N-30)*8,436,285                                       (19)
```

such cores. Combining the pairwise-disjoint all-three-original,
exactly-two-original, and zero-original slices of THM-4174/4175 and this
theorem gives

```text
888,030+6,660,225+8,436,285=15,984,540                 (20)
```

certified cores per newcomer. At `q=50`, THM-4179's proved exactly-one slice
adds `1,071,961`, for a four-slice total of

```text
17,056,501.                                            (21)
```

Equation `(21)` is only a q=50 synthesis here; a uniform exactly-one theorem
requires its own proof.

## 6. Independent audit contract

The primary implementation reuses THM-4188's fixed-pool prefix cells and
integrates each newcomer comb over them. It enumerates ten-subsets with
Gosper masks and locks the complete positive ledger at

```text
FNV1A64_LE=14e15a5a00ad0764.                           (22)
```

The independent implementation reuses THM-4188's separate joint-wall
geometry: each newcomer wall is explicitly merged with the fixed-pool wall
arrangement. It recursively generates the ten-subsets and locks its distinct
encoding at

```text
FNV1A64_LE=18bb68e547da2a9e.                           (23)
```

They agree on every load-bearing edge count, equality count, body count,
no-cover verdict, edge-check total, closest body, and missed edge. They also
agree on the positive q=50 E7 control and the negative q=50 E6 hostile.

## 7. Scope firewall

This theorem proves a complete labelled slice inside one fixed thirty-label
pool. It does not show that every arbitrary eleven-body has mass `4/63`,
provide entry of a hypothetical LRC counterexample into this pool, close the
exactly-one-original slice uniformly, handle other parity branches, or prove
LRC(14). The connection contract is

```text
source:       lawful global repair edges in P
target:       every ten-subset of U
map:          R -> R intersect U
preserved:    q, full repair mass, labelled disjointness, threshold
destroyed:    identities of deleted original labels after certification
sidecar:      the full unprojected repair edge
decisive test: exact projected transversal number through size ten. (24)
```

## 8. Exact artifacts and replay

Primary:

```text
04-computation/lrc14_all_body_zero_original_newcomer_haar_transfer_thm4191.cpp
sha256 8bc1213c253cd450bd797af3e8e926dd69adc145af7eb85459b2b3acf0bbb7f2

05-knowledge/results/lrc14_all_body_zero_original_newcomer_haar_transfer_thm4191.out
sha256 9f2ae05d3c362190533514f9fc919d8b324e7e3d64230b79ab40b086e8b5493c
```

Independent:

```text
04-computation/lrc14_all_body_zero_original_newcomer_haar_transfer_independent_audit_thm4191.cpp
sha256 ce956bbcbe3bdd940b294ad3c1cde7cf953677c77e92625c8f6494894886f0e2

05-knowledge/results/lrc14_all_body_zero_original_newcomer_haar_transfer_independent_audit_thm4191.out
sha256 5b09a62c4746839d0ca7e98f0e4297be929cb64b4b4d0e7587748dc47320293e
```

Replay:

```bash
g++ -std=c++20 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_all_body_zero_original_newcomer_haar_transfer_thm4191.cpp \
  -o /tmp/lrc4191-primary
/tmp/lrc4191-primary | diff -u \
  05-knowledge/results/lrc14_all_body_zero_original_newcomer_haar_transfer_thm4191.out -

g++ -std=c++20 -O2 -DNDEBUG -pthread \
  -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_all_body_zero_original_newcomer_haar_transfer_independent_audit_thm4191.cpp \
  -o /tmp/lrc4191-independent
/tmp/lrc4191-independent | diff -u \
  05-knowledge/results/lrc14_all_body_zero_original_newcomer_haar_transfer_independent_audit_thm4191.out -
```
