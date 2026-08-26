---
id: THM-4166
title: "Two-deletion repair-graph Haar odd-tail transfer"
status: >
  PROVED RELATIVE TO THM-4150/4156 + VERIFIED-EXACT GLOBAL CENSUS +
  INDEPENDENT SEQUENTIAL-FRACTION/GLOBAL-LATTICE/C++ AUDITS; LRC(14) OPEN.
  The exact two-deletion repair graph qualifies 1,032 and only 1,032
  newcomer labels for the fixed THM-4156 anchored pool. They close
  916,446,960 primitive divisor-complete eleven-bodies, or 918,667,035
  together with the old THM-4156 family. The last qualifier is q=8,265;
  an exact analytic supergraph excludes every q>=49,494 and the complete
  intervening census proves that no qualifier re-enters after 8,265.
source: codex-lrc14-jc-sharp-fronts-20260825b
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
related:
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-4148-first-window-width-universal-odd-tail-lrc14-transfer
  - THM-4151-scale-sensitive-first-window-odd-tail-lrc14-transfer
  - THM-4154-mod-six-fixed-clock-and-haar-pool-inheritance-correction
  - THM-4160-anchored-haar-deletion-cover-and-content-tower
script: 04-computation/lrc14_two_deletion_repair_graph_thm4166.py
output: 05-knowledge/results/lrc14_two_deletion_repair_graph_thm4166.out
independent_audit_source: 04-computation/lrc14_two_deletion_repair_graph_thm4166_independent_audit.cpp
independent_audit_output: 05-knowledge/results/lrc14_two_deletion_repair_graph_thm4166_independent_audit_cpp.out
global_census_script: 04-computation/lrc14_two_deletion_repair_graph_thm4166_global_census.py
global_census_output: 05-knowledge/results/lrc14_two_deletion_repair_graph_thm4166_global_census.out
global_census_independent_source: 04-computation/lrc14_two_deletion_repair_graph_thm4166_global_census_independent.cpp
global_census_independent_output: 05-knowledge/results/lrc14_two_deletion_repair_graph_thm4166_global_census_independent_cpp.out
script_sha256: e5f91aaaff3876568dbcace7623bdf5c03c8643cc450c5d2945d74ff449d1e55
output_sha256: 276be2254d9e9311fb00a60d7e8397af5d179943aab134fe5530ace319860c86
independent_audit_source_sha256: bc3c3547ebdc6e9191dd5f390a01d759d2223b7faf5eeffa2d065b29621e61b1
independent_audit_output_sha256: 9dda2b76fe4701c03359da9ae4dc1405da1939cf79354d1f99cc977e2653dbe5
global_census_script_sha256: 4cc8b067458a365c3d7f13fe9af1b0168157f66d573eeeaa99b211f9f1a0fa5a
global_census_output_sha256: b9a06d4cb1c66bf9f41ef07ef5aa224c216ead6882a8694161844737ddf36708
global_census_independent_source_sha256: 17acfe56f3d121858877ff8eeb701b2a70bde75c52a1d2e9079ddd4b931f75f6
global_census_independent_output_sha256: 3ff8ab357211c5f51cc9f73001bc5e77c375298baee9e032903b15ee32ed1cb1
q_le_200_semantic_ledger_sha256: 13404c6c2986bd5a14bad57519eb30792adc8e8ea33146b7c68969e33b255394
global_semantic_word_xor_mul64: 995aa971af1069e4
hash_basis: raw LF bytes
primary_audit: >
  PASS. Sequential exact Fraction intersections independently reconstruct all
  351 leave-two banks and every Gamma_q through q=200. They recover the 53
  qualifiers, every exact alpha/tau value, both strict threshold margins,
  and positive and hostile controls. Normal, optimized, and hash-seeded
  outputs byte-match.
independent_audit: >
  ACCEPT. A C++20 global arrangement of all pool and q<=200 walls uses
  arbitrary-precision rational atom lengths. For every one of the 175 labels
  it exhausts all C(27,20)=888,030 possible twenty-vertex independent sets.
  Its complete adjacency ledger byte-matches the Fraction ledger. O0, O2,
  and UBSan outputs byte-match under warnings-as-errors builds.
global_census_audit: >
  PASS. A Fraction-built 7,134-wall arrangement embeds exactly in a checked
  45-bit common lattice. Vectorized signed-int64 arithmetic evaluates all
  49,463 labels through the proved cutoff, derives the cutoff itself, and
  freezes the complete qualifier list, alpha/tau histogram, and semantic
  fingerprint. Normal, optimized, and hash-seeded outputs byte-match.
global_census_independent_audit: >
  ACCEPT. Independent C++20 signed-128 arithmetic integrates the newcomer
  comb over all 7,133 fixed-pool cells, scans every one of the 49,463 labels,
  reproduces both exact margin extrema, the full qualifier list and the same
  semantic fingerprint. A separate depth-seven vertex-cover recursion agrees
  with maximum-independent-set duality on all 49,463 graphs. O0, O2, and
  UBSan outputs byte-match under warnings-as-errors builds.
---

# THM-4166 -- two-deletion repair graphs

**PROVED RELATIVE TO THM-4150/4156 + VERIFIED-EXACT GLOBAL CENSUS + FOUR
CROSS-CHECKING AUDITS; LRC(14) REMAINS OPEN.**

## 1. Statement

Retain the THM-4156 anchors and pool

```text
A={120,126,143},

P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290},

O=P\A.                                                   (1)
```

For a finite positive set `S`, write

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14}.                 (2)
```

For a positive integer `q notin P`, define a graph `Gamma_q` on the 27
vertices `O` by

```text
{r,s} in E(Gamma_q)
  iff mu(G_((P union {q})\{r,s}))>=4/63.                 (3)
```

The comparison in `(3)` is deliberately nonstrict. Let

```text
Q_2={q>=1:q notin P and tau(Gamma_q)>7},                 (4)
```

where `tau` is minimum vertex-cover number.

> **Theorem.** The set `Q_2` has exactly `1,032` elements. Its complete
> increasing list is frozen in both global census outputs. Its largest
> member is `8,265`; in particular no `q>=8,266` belongs to `Q_2`.
>
> For every `q in Q_2`, every `K in binom(O,7)`, every positive integer
> `c`, and every two distinct positive odd integers `a,b`, there is an
> `x in R/Z` such that, for
>
> ```text
> H=A union {q} union K,
> ```
>
> one has
>
> ```text
> min_(v in 2cH union {a,b})||vx||>=1/14.                (5)
> ```

There are exactly

```text
1,032 binom(27,7)=916,446,960                            (6)
```

new scale-one bodies. They are pairwise distinct because each has the unique
label `q` outside `P`. They are disjoint from the old THM-4156 family, so

```text
916,446,960+binom(27,8)=918,667,035                      (7)
```

bodies are closed by the union of the two results. Every scale-one body is
primitive because it contains `A` and `gcd(A)=1`; the anchors also provide a
multiple of every divisor `2,...,14`, exactly as in THM-4156. Equation `(6)`
counts scale-one bodies, while `(5)` retains every positive common content.

## 2. The deletion-hypergraph lemma

The graph criterion is the `d=2` case of an elementary operation that also
identifies the next frontier. Fix `q notin P` and an integer `d>=1`. Define a
`d`-uniform hypergraph on `O` by

```text
E_q^(d)={R in binom(O,d):
          mu(G_((P union {q})\R))>=4/63}.                (8)
```

Call `K in binom(O,7)` **d-deletion certified** if some `R in E_q^(d)` is
disjoint from `K`. For such an `R`,

```text
A union {q} union K subset (P union {q})\R,

G_((P union {q})\R) subset G_(A union {q} union K).      (9)
```

Thus its body-safe mass is at least `4/63`.

A set `K` fails this certificate exactly when it meets every hyperedge, that
is, when it is a transversal. Consequently

```text
every K in binom(O,7) is d-deletion certified
  iff tau(E_q^(d))>7.                                   (10)
```

Indeed, if a transversal has size at most seven, extend it inside the
27-element set `O` to a seven-set. Conversely, a seven-set meeting every
hyperedge is itself a transversal. This proves both directions, including
the empty-hypergraph boundary.

For `d=1`, the hyperedges are the singleton repair labels of THM-4160 and
their transversal number is their cardinality. For `d=2`, `(8)` is precisely
the graph `(3)`. Since `|O|=27`, graph duality gives

```text
tau(Gamma_q)=27-alpha(Gamma_q),                          (11)
```

so `(4)` is equivalently `alpha(Gamma_q)<20`.

## 3. Why two deletions strictly extend one deletion

Let `D_1(q)` be THM-4160's one-deletion repair set. If `r in D_1(q)`, then
for every `s in O\{r}` monotonicity gives

```text
G_((P union {q})\{r}) subset
G_((P union {q})\{r,s}).                                (12)
```

Hence `Gamma_q` contains the full star centered at `r`. If
`1<=|D_1(q)|<=26`, the union of these full stars has vertex-cover number
exactly `|D_1(q)|`: in the star-union graph all noncentres form an independent
set, while any independent set containing a centre has size one. Therefore

```text
tau(Gamma_q)>=|D_1(q)|.                                 (13)
```

All ten THM-4160 qualifiers have at least eight one-deletion repairs and are
therefore automatically in `Q_2`; the global census checks all ten directly.
Thus the two-deletion family subsumes the one-deletion family. The gain beyond
THM-4160 is

```text
(1,032-10)binom(27,7)=907,566,660 bodies.                (14)
```

The corrected one-deletion hostile `q=235` exposes the strict gain:

```text
|E(Gamma_235)|=238,   alpha(Gamma_235)=10,
tau(Gamma_235)=17.                                      (15)
```

Distributed pair repairs succeed where no eighth star centre exists.

## 4. Exact analytic cutoff

For `{r,s} in binom(O,2)`, put

```text
U_(r,s)=G_(P\{r,s}),
M_(r,s)=mu(U_(r,s)),
c_(r,s)=number of positive-length components of U_(r,s). (16)
```

Let `g(t)=1_[1/14,13/14]({t})`. Its mean is `6/7`, and the periodic primitive
of `g-6/7` has oscillation `6/49`. Exact integration over a union of `c`
intervals therefore gives

```text
|mu(U_(r,s) intersect G_q)-(6/7)M_(r,s)|
  <=6c_(r,s)/(49q).                                     (17)
```

The left side of `(3)` is exactly the intersection in `(17)`. Call
`{r,s}` a stable-limit edge when

```text
(6/7)M_(r,s)>4/63,       equivalently M_(r,s)>2/27.      (18)
```

Exact reconstruction finds 39 such edges and no equality in `(18)`. Their
graph `L` has

```text
alpha(L)=21,                 tau(L)=6.                   (19)
```

For each of the other 312 pairs define

```text
B_(r,s)={6c_(r,s)/49}/[4/63-(6/7)M_(r,s)].              (20)
```

If `q>B_(r,s)`, `(17)` makes that pair strictly fail `(3)`. Thus the actual
repair graph is a subgraph of the nested analytic supergraph

```text
Gamma_q subset S_q
 :=L union {{r,s}: floor(B_(r,s))>=q}.                   (21)
```

At the decisive transition, the edge `{170,240}` has

```text
c_(170,240)=148,
M_(170,240)=3480322961/47256889680,
4/63-(6/7)M_(170,240)=60562157/165399113880,
B_(170,240)=2997437002560/60562157,
floor(B_(170,240))=49,493.                              (22)
```

The exact supergraphs on the two sides are

```text
q=49,493: |E(S_q)|=53, alpha=19, tau=8,
q=49,494: |E(S_q)|=52, alpha=20, tau=7.                  (23)
```

At `q=49,494`, one seven-cover is

```text
{85,88,95,168,193,252,290}.                             (24)
```

The supergraphs decrease with `q`. Hence `(21)` and `(24)` give

```text
tau(Gamma_q)<=7                    for every q>=49,494.  (25)
```

This is a global theorem, not a search truncation. The analytic supergraph
still has `tau=8` at `49,493`; the exact graph there already has 38 edges and
`tau=6`.

## 5. Complete exact census below the cutoff

The remaining universe is

```text
1<=q<=49,493,                 q notin P,                 (26)
```

containing 49,463 labels. The fixed-pool walls have the integer sidecar

```text
7,134 walls,       7,133 open cells,
L=lcm{14p:p in P}=18,241,159,416,480.                   (27)
```

Every wall lies on `L^(-1)Z`. Exact midpoint classification of the cells by
failed optional labels gives

```text
no optional failure:150,  one failure:328,  two failures:518,
ignored (anchor failure or >=3 optional failures):6,137. (28)
```

Only the first three classes can contribute after two optional deletions.
For a wall `z=p/d` in lowest terms, where `d|L`, write `qp=kd+r` with
`0<=r<d`. The exact numerator `J_q(z)` of the safe-comb integral through `z`,
on denominator `14qL`, is

```text
J_q(z)=12kL+
  0,                         if 14r<=d,
  (14r-d)L/d,                if d<14r<13d,
  12L,                       if 13d<=14r.               (29)
```

Subtracting `(29)` at consecutive walls gives a nonnegative exact cell
contribution. For fixed `q`, accumulate these as `b_q` for zero failures,
`u_(q,r)` for one failure `r`, and `v_(q,r,s)` for the pair `{r,s}`. Then

```text
mu(G_((P union {q})\{r,s}))
 ={b_q+u_(q,r)+u_(q,s)+v_(q,r,s)}/{14qL},               (30)

{r,s} in E(Gamma_q)
 iff 9[b_q+u_(q,r)+u_(q,s)+v_(q,r,s)]>=8qL.             (31)
```

Thus every graph edge in `(26)` is one checked integer comparison. The
Fraction/lattice script checks all signed-int64 bounds before vectorization;
the independent C++ implementation uses signed 128-bit arithmetic.

The exact vertex-cover histogram is:

|`tau`|0|1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|count|45|127|124|596|793|6241|38003|2502|377|435|81|45|32|19|11|13|4|5|3|5|2|

Exactly 1,032 rows have `tau>7`. The full qualifier list and all
`(q,|E|,alpha,tau)` rows are frozen in both global outputs under the following
wordwise 64-bit XOR/multiply recurrence:

```text
h_0=0xcbf29ce484222325,
h_(j+1)=((h_j XOR u_j)*0x100000001b3) mod 2^64,
word-hash=995aa971af1069e4.                             (32)
```

For each `q` in increasing order, the unsigned integer words are
`q,|E|,alpha,adjacency[0],...,adjacency[26]`. There is deliberately no byte
serialization; this is a reproducibility sidecar, **not** standard bytewise
FNV-1a-64 (MISTAKE-513).

The maximum cover number is 20, attained exactly at `q=380,386`. The final
positive and immediate hostile are

```text
q=8,265: |E|=52, alpha=19, tau=8,
q=8,266: |E|=45, alpha=21, tau=6.                       (33)
```

The census has no re-entry through `49,493`; `(25)` excludes every larger
integer. Among all 17,361,513 finite edge comparisons in `(26)`, there is no
threshold equality. The smallest strict surplus and deficit are

```text
+109205/929387072269656 at (q,r,s)=(42,798,145,168),
-3737/5186443642380      at (q,r,s)=(878,95,193).        (34)
```

The signs are relative to `4/63`. These extrema are sidecars; the nonstrict
convention in `(3)` remains part of the theorem.

## 6. Independent `q<=200` validity slice

The sequential Fraction implementation reconstructs all 351 leave-two banks
and intersects each separately with every `G_q` for

```text
1<=q<=200,                   q notin P.                  (35)
```

It finds 53 qualifiers. Its smallest positive edge margin is

```text
195043643/362999072387952 at (q,r,s)=(199,80,264),       (36)
```

and its closest failure is

```text
2193407/22921456898340 at (q,r,s)=(191,42,290).          (37)
```

The independent C++ program instead forms the global arrangement of every
pool and `q<=200` wall and sums arbitrary-precision atom lengths. For each of
the 175 graphs it exhausts all

```text
binom(27,20)=888,030
```

candidate twenty-vertex independent sets, for 155,405,250 tests in total.
By `(11)`, a graph qualifies exactly when none survives. The complete
adjacency/alpha ledgers byte-match, with SHA-256

```text
13404c6c2986bd5a14bad57519eb30792adc8e8ea33146b7c68969e33b255394. (38)
```

This slice is a hostile independent check on the prefix-integral census.

## 7. Haar transfer, content, and QED

Fix `q in Q_2` and `K in binom(O,7)`. Equations `(10)--(11)` provide an edge
`{r,s}` disjoint from `K`. By `(3)` and `(9)`,

```text
mu(G_(A union {q} union K))>=4/63.                       (39)
```

THM-4150 applied to `(39)` supplies `(5)` at content one. For every positive
integer `c`, the surjective circle endomorphism `m_c(y)=cy` satisfies

```text
G_(cH)=m_c^(-1)(G_H),             mu(G_(cH))=mu(G_H).    (40)
```

Applying THM-4150 to `cH` proves `(5)` at every content. The body speeds are
even and the tails odd, so no body--tail collision occurs; the tails are
distinct by hypothesis. Equations `(6)--(7)` give the counts. **QED.**

## 8. Reproduction and scope

Primary replay:

```bash
python3 -B 04-computation/lrc14_two_deletion_repair_graph_thm4166.py
python3 -B -O 04-computation/lrc14_two_deletion_repair_graph_thm4166.py
PYTHONHASHSEED=271828 python3 -B \
  04-computation/lrc14_two_deletion_repair_graph_thm4166.py
```

Global Fraction/lattice replay:

```bash
python3 -B 04-computation/lrc14_two_deletion_repair_graph_thm4166_global_census.py
python3 -B -O 04-computation/lrc14_two_deletion_repair_graph_thm4166_global_census.py
PYTHONHASHSEED=314159 python3 -B \
  04-computation/lrc14_two_deletion_repair_graph_thm4166_global_census.py
```

For each C++ source, replace `SOURCE` and `BINARY` below by the corresponding
checked-in path and a temporary output path:

```bash
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror SOURCE -o BINARY
BINARY
g++ -std=c++20 -O0 -Wall -Wextra -Wpedantic -Werror SOURCE -o BINARY
BINARY
g++ -std=c++20 -O1 -g -fsanitize=undefined \
  -fno-sanitize-recover=undefined -Wall -Wextra -Wpedantic -Werror \
  SOURCE -o BINARY
BINARY
```

The source-to-target contract is

```text
source:       exact two-deletion Haar repair hyperedges on the THM-4156 pool
target:       every anchored seven-choice body and every distinct odd pair
map:          missed graph edge -> safe-set inclusion -> THM-4150 transfer
preserved:    anchors, exact 4/63 threshold, content, divisor completeness
destroyed:    component addresses and body-specific safe phases
sidecar:      vertex-cover duality, safe-comb discrepancy, common wall lattice
positive:     q=8265 has tau=8; q=235 repairs the one-deletion near miss
hostile:      q=8266 has tau=6; the q=49494 supergraph has a seven-cover
decisive test: global exact alpha/tau census plus analytic tail cutoff. (41)
```

This theorem classifies one precise two-deletion mechanism for one fixed
pool. It does not prove that these are all bodies accessible to THM-4150,
establish optimality of the pool, treat arbitrary tail parity, or prove
LRC(14). The `d=3` hypergraph in `(8)` is the immediate next operation:
graph vertex-cover duality is then replaced by a genuine three-uniform
transversal problem.
