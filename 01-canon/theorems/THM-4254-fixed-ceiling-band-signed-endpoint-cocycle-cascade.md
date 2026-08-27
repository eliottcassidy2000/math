---
id: THM-4254
title: "Fixed ceiling-band signed endpoint-cocycle cascade"
status: >
  PROVED RELATIVE TO THM-4228/4231/4233/4238/4242/4252
  + FINITE-EXACT + INDEPENDENTLY AUDITED. For every current post-THM-4252
  residual pair whose second endpoint lies from 755 through 768, and every
  nine-body in the fixed thirty-label pool, adjoining the pair leaves
  Haar-safe mass at least 4/63. There are exactly 59 such pairs. An unlabelled
  4,675-mask carrier supports 56,419 pair-labelled prefix occurrences and
  replaces 59 full superset transforms at proof-check time. Removing the band
  leaves 181,064 residual edges, maximum endpoint 754.
  THM-4252 is proved. This theorem does not prove an arbitrary pair, a ray,
  or LRC(14).
source: root/cross-frontier-bridge/2026-08-26
depends_on:
  - THM-4228-common-gcd-two-outsider-periodic-observable-haar-ray
  - THM-4231-arbitrary-pair-cofinal-depth-six-haar-repair-and-exact-outsider-lift
  - THM-4233-pair-specific-primitive-observable-oscillation-haar-charts
  - THM-4238-forty-small-label-uniform-r590-haar-tail-closure
  - THM-4242-fixed-fifty-direct-r590-tail-and-twenty-three-label-chart
  - THM-4252-exact-signed-endpoint-cocycle-residual-edge-closure
related:
  - THM-4245-primitive-observable-component-floor-and-cofinal-gate-redundancy
primary_script: 04-computation/cascade_pair_exhaustive_primary.cpp
primary_script_sha256: c4acbbc18eb0d9b3bb3105efe090ae4bb08ccc15ac8c230b4af14dec0db00627
primary_replay_directory: 05-knowledge/results/lrc14_endpoint_cascade_thm4254/replay_band
primary_replay_manifest: 05-knowledge/results/lrc14_endpoint_cascade_thm4254/SHA256SUMS
primary_replay_manifest_sha256: 531d6b48cd50b999d7dfb8c209ba4bd7f5b04f85dd46ce7d3934a360e3005a1d
primary_replay_manifest_basis: >
  59 LF-terminated lines generated inside the replay directory by sorted
  SHA-256 over the relative names ./q_r.out; the manifest is 4,720 bytes.
primary_normalized_aggregate_sha256: 84aefdd53bd018b89d2d5c039ff711d74798807fc246cea79d6f18b68023a145
primary_dependency_script: 04-computation/dependency_lrc14_all_newcomer_zero_original_anchor_hierarchy_thm4188.cpp
primary_dependency_script_sha256: 25a6978484c7ab122fdc6c8e1593cfa2ad3468f7184a156045ea7e6cb2efc45d
compressed_checker: 04-computation/cascade_prefix_union_exact_verifier.cpp
compressed_output: 05-knowledge/results/lrc14_endpoint_cascade_prefix_union_thm4254.out
compressed_checker_sha256: 43ecbcd022e88ca269609a44222042d374d4390c36fa0c1eb347dc31bdad3382
compressed_output_sha256: a018cfd6447147202c539bbad475cfa52ccbba1aca4a0ab81d4fa048ec31f1fe
independent_audit_script: 04-computation/lrc14_endpoint_cascade_direct_wall_body_audit.cpp
independent_audit_output: 05-knowledge/results/lrc14_endpoint_cascade_direct_wall_body_audit_thm4254.out
independent_audit_script_sha256: bbfc55a3181882cf6456900951658f95634680fb30f68ad763ffa830812c66e7
independent_audit_output_sha256: 6fb7ecf203e4966c1bda9c10fc8a393dd73a284af8a2b5d7cab74b0ad4e517c7
residual_postprocess_script: 04-computation/cascade_residual_postprocess_thm4254.py
residual_postprocess_output: 05-knowledge/results/lrc14_endpoint_cascade_residual_postprocess_thm4254.out
residual_postprocess_script_sha256: f2c29193538de9f81a4f7c5c4649107b7dfc3058f4787578643e8df4a94d1b0b
residual_postprocess_output_sha256: 5428e430b001a25fc932bf32b96bfe976b806dda195803bd2782dfa72fb35785
audit: >
  PASS / ACCEPT. THM-4252 is proved. The discovery compiler thresholds all
  5,852,925 repairs and scans all 14,307,150 bodies for each pair. The compact
  checker verifies only 56,419 pair-labelled prefix occurrences directly from
  signed endpoint atoms and repeats all 844,121,850 body cases without a zeta
  transform. A clean-room companion instead builds literal joint walls,
  integrates each occurrence directly, and recursively enumerates every body.
  Within each checker family, O0/O3/ASan+UBSan outputs are byte-identical.
---

# THM-4254 -- fixed ceiling-band signed endpoint-cocycle cascade

**PROVED RELATIVE TO
THM-4228/4231/4233/4238/4242/4252 + FINITE-EXACT + INDEPENDENTLY AUDITED;
LRC(14) REMAINS OPEN.**

## 1. Statement and exact universe

For a finite positive set `A`, write

```text
G_A={x in R/Z:min_(a in A)||ax||>=1/14},   alpha=4/63,
```

and retain the fixed labelled pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,264,286,290}.   (1)
```

Let `C` be the following 59-pair set:

| second endpoint | first endpoints |
|---:|---|
| 768 | 616, 721, 744, 750, 765, 766 |
| 767 | 616, 721 |
| 766 | 616, 721, 744, 750, 765 |
| 765 | 616, 704, 721, 744, 750 |
| 764 | 616, 704, 721 |
| 763 | 616, 704, 721, 732, 744, 750 |
| 762 | 616, 704, 721, 732, 744, 750 |
| 761 | 616, 704, 721, 726, 732 |
| 760 | 616 |
| 759 | 616, 704, 721 |
| 758 | 616, 704, 721, 726 |
| 757 | 616, 698, 704, 721, 726, 732 |
| 756 | 616 |
| 755 | 616, 698, 704, 721, 726, 732 |

Then for every `(q,r) in C` and every `B in binom(P,9)`, one has

```text
mu(G_(B union {q,r})) >= 4/63.                           (2)
```

Within the fixed band, no edge was selected or omitted by outcome.
Independently recomputing the current THM-4231/4238/4242 residual, deleting
THM-4252's three proved edges, and applying the semantic filter
`755<=second endpoint<=768` gives exactly this set, with

```text
COUNT 59
FNV-1a-64 b3d54b78babbcaec
SHA256(lines q,r\n) 6b54d8fa3b408325fc309bec3ed769f5e56ce370fa34fa7ad1bb6d7ed4cafc36.  (3)
```

## 2. Exact endpoint cocycle and denominator cancellation

Fix `(q,r)=g(u,v)` with `gcd(u,v)=1`, and put

```text
N=14uv,   A=G_u intersect G_v,   beta=mu(A)=S/N,
D=lcm(14s:s in P)=18,241,159,416,480.                  (4)
```

On the primitive wall cells

```text
0=e_0<e_1<...<e_s=N,
```

let `A_j` be the indicator of `A`, `d_j=e_(j+1)-e_j`, and define

```text
C_0=0,   C_(j+1)=C_j+d_j(NA_j-S).                      (5)
```

For `H(x)=integral_0^x(1_A-beta)`, direct affine integration gives

```text
H(e_j/N)=C_j/N^2.                                      (6)
```

For a fixed-pool tick `t`, let `p_t` be the representative of `gt mod D` in
`[0,D)`, choose `j` with `D e_j<=N p_t<=D e_(j+1)`, and set

```text
Z_g(t)=D C_j+(N p_t-D e_j)(N A_j-S).                   (7)
```

The adjacent formulas agree at a wall, and

```text
H(gt/D)=Z_g(t)/(D N^2).                                (8)
```

For `R in binom(P,8)`, write `U_R=G_(P\R)`. Its
**positive-length circular components** have nonwrapping lifts
`[a_i/D,b_i/D]`; isolated safe points are excluded and do not affect Haar
mass. Put

```text
m_R=sum_i(b_i-a_i),
K_R=sum_i(Z_g(b_i)-Z_g(a_i)).                          (9)
```

The chain rule, including its factor `1/g`, yields

```text
mu(U_R intersect G_q intersect G_r)
  =(g N S m_R+K_R)/(g D N^2).                         (10)
```

Therefore `R` is active exactly when

```text
63(g N S m_R+K_R) >= 4gDN^2.                          (11)
```

The implementation stores the correctly cancelled numerator

```text
M_R=g S m_R+K_R/N
```

on denominator `gDN`, and tests `63M_R>=4gDN`. Here `K_R/N` is integral
because the program's endpoint value is `Z_g(t)/N`, the numerator of `H(gt/D)`
on denominator `DN`. Thus `(11)` and the implemented predicate differ by one
cancelled common factor `N`, not by a changed threshold.

## 3. Repair hypergraph and consequence direction

For a fixed pair define

```text
E(q,r)={R in binom(P,8):R satisfies (11)}.              (12)
```

The discovery calculation evaluates every one of

```text
binom(30,8)=5,852,925                                  (13)
```

repairs. If a body `B in binom(P,9)` is disjoint from `R in E(q,r)`, then

```text
B subset P\R,
G_((P\R) union {q,r}) subset G_(B union {q,r}).         (14)
```

The inclusion points from the more constrained repair set to the target body
set. Hence the active repair's mass lower-bounds the target mass, proving
`(2)`. Equivalently, every labelled nine-set fails to be a transversal of
`E(q,r)`.

## 4. Exact exhaustive discovery and compressed certificate

For each of the 59 pairs, the primary program:

1. reconstructs the 7,133 fixed-pool cells;
2. builds the primitive endpoint cocycle and checks every cell against a
   separately evaluated direct pair-safe prefix integral;
3. applies the labelled ordinary-colex superset transform to all 5,852,925
   repairs (`152,170,690` exact additions per pair);
4. thresholds by `(11)` and freezes the active-deck ledger; and
5. scans all `14,307,150=binom(30,9)` bodies.

All 59 scans have zero failures and zero activation equalities. Across the
band there are `258,307,684` active-repair occurrences. The 59 verbatim
typed transcripts live in
`05-knowledge/results/lrc14_endpoint_cascade_thm4254/replay_band`; the
adjacent `SHA256SUMS` manifest commits every file. Their deterministic
pair-order concatenation after stripping only the two timing fields has
SHA-256 `84aefdd53bd018b89d2d5c039ff711d74798807fc246cea79d6f18b68023a145`. Active repairs are
ordered by

```text
(SplitMix64(mask xor 0x4245422842334245),mask).         (15)
```

For each pair, the first `K(q,r)` active repairs already cover every body.
The prefix sizes range from 438 to 3,205 and sum to 56,419. In every row the
frozen worst body intersects the first `K-1` repairs and is disjoint from the
last repair. This proves prefix minimality **only within order `(15)`**; no
globally minimum deck is claimed.

Taking distinct masks in first-occurrence order over the sorted 59 pair rows
produces only 4,675 masks:

```text
UNIQUE PREFIX UNION COUNT 4675
FNV-1a-64 over little-endian u64 masks ce4e76ec11df057c.                (16)
```

This is a **pair-labelled union**, not a common-active deck. A mask may be
active for one row and inactive for another. The compact checker reads each
row's own labelled prefix, recomputes all 56,419 masses by direct summation of
the endpoint atoms, and repeats all

```text
59*binom(30,9)=844,121,850                             (17)
```

body cases. It performs 153,346,842 atom-membership tests and no zeta
transform. The pair-labelled prefix-incidence ledger, which retains `q,r,K`
and duplicates across rows, is `acc8347addf27ac3`.

Thus the 4,675-mask union is a compact carrier for 59 explicitly typed
certificates; it is not falsely asserted that every carrier mask works for
every pair.

## 5. Independent literal-wall audit

The clean-room checker does not evaluate `(5)--(11)`. For each pair it builds

```text
L=lcm(14s:s in P union {q,r})
```

and all literal safe/unsafe walls for the 32 speeds. On each resulting open
cell it evaluates pool failures and `q,r` safety directly at an exact
midpoint. For every emitted repair it sums literal cell widths and tests
`63*mass>=4L`. It then recursively enumerates every body and reproduces every
order-minimality witness.

Across all 59 rows it checks 590,771 joint cells, 56,419 direct repair
integrals, 844,121,850 body cases, and 25,170,340,269 body/repair incidences.
There are zero failures, and every final direct mass reduces to the same
fraction as the primary cocycle mass.

The largest literal grid is

```text
L(698,757)=4,819,186,629,718,100,640 < 2^63-1.
```

Although each wall fits signed 64-bit, `left+right` can overflow. The audited
source promotes before addition; `(698,757)` is the explicit arithmetic
hostile. O0, O3, and ASan+UBSan transcripts are byte-identical.

## 6. Proof-graph consequence

The post-THM-4242 residual has 181,126 edges. THM-4252 removes three, leaving

```text
COUNT 181123
FNV 6ec03ed4c4dc841b
SHA256 9a9b6fbe14db00e9d7f8f08ecddaa1e3d263fd063c6b3c003e18c210b3334ef8.
```

Removing the 59 proved rows of `C` leaves

```text
COUNT 181064
FNV 8f550dcc2e552962
SHA256 0167652b41139bfd00c52236338fdd50e3be604641fe03e71eb66c68ee497d35,
MAX ENDPOINT 754,
TOP LAYER (616,754),(698,754),(704,754),(721,754),
CUTOFF 755.                                                   (18)
```

The dependency graph is

```text
THM-4228 + THM-4233  -> exact endpoint transfer
THM-4231 -> THM-4238 -> THM-4242 -> THM-4252 -> inherited residual
                                                \
                                                 -> THM-4254.
```

THM-4254 uses THM-4252 as a proved dependency; detaching it would be invalid.

## 7. Hostile boundary and scope

As an exploratory hostile, the frozen 4,675-mask union was tested on all 297
post-THM-4252 residual edges with second endpoint 733 through 754. It closes
all of them. At endpoint 732, the first compression failure is `(542,732)`:
its pair-active union subdeck has 3,227 repairs and leaves exactly two bodies
uncovered. The first is mask `0x151a5400`; among union repairs disjoint from
that body, the best activation numerator is still

```text
-86,086,766,768,380,560
```

on denominator `50,659,493,860,723,587,840`.

This is a boundary of the **frozen compact union**, not of the endpoint method:
the full 5,852,925-repair computation for `(542,732)` does close every body,
with 2,614,902 active repairs and a 4,518-repair prefix. These exploratory
successes are not added to statement `(2)` or residual `(18)` without their
own frozen universe and independent promotion audit.

Nothing here proves a dilation ray, arbitrary-pair entry, a common-active
deck, global prefix optimality, literal safety outside `C`, or LRC(14).
