---
id: THM-4160
title: "Anchored Haar deletion cover and content tower"
status: >
  PROVED RELATIVE TO THM-4150/4156 + VERIFIED-EXACT + INDEPENDENT GLOBAL-WALL
  FRACTION AUDIT + INDEPENDENT CHECKED-INTEGER C++ AUDIT; LRC(14) OPEN. An
  exact one-deletion cover adds ten and only ten singleton labels to the
  THM-4156 anchored family, closing 8,880,300 new primitive eleven-bodies and
  11,100,375 bodies in total with every distinct positive odd tail pair.
  This one-deletion mechanism admits no forced newcomer set of size 2--8.
  Exactly 9,703,274 bodies remain beyond both named THM-4148/4151 gates at
  every positive common body content. FINITE-EXACT at content one, exactly
  9,036,418 cores lie outside every THM-4158 wrapped-carrier alphabet,
  including 7,267,924 of the new singleton-newcomer cores; this last
  separation is not dilation invariant.
source: codex-lrc14-jc-sharp-fronts-20260825b
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
related:
  - THM-366-lrc-small-denominator-divisibility-sieve
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-4148-first-window-width-universal-odd-tail-lrc14-transfer
  - THM-4151-scale-sensitive-first-window-odd-tail-lrc14-transfer
  - THM-4154-mod-six-fixed-clock-and-haar-pool-inheritance-correction
  - THM-4158-three-band-wrapped-carrier-odd-tail-lrc14-transfer
script: 04-computation/lrc14_anchored_deletion_cover_thm4160.py
output: 05-knowledge/results/lrc14_anchored_deletion_cover_thm4160.out
independent_audit_script: 04-computation/lrc14_anchored_deletion_cover_thm4160_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_anchored_deletion_cover_thm4160_independent_audit.out
integer_audit_source: 04-computation/lrc14_anchored_deletion_cover_thm4160_independent_audit.cpp
integer_audit_output: 05-knowledge/results/lrc14_anchored_deletion_cover_thm4160_independent_audit_cpp.out
content_tower_script: 04-computation/lrc14_anchored_deletion_cover_thm4160_content_tower_audit.py
content_tower_output: 05-knowledge/results/lrc14_anchored_deletion_cover_thm4160_content_tower_audit.out
anchor_audit_script: 04-computation/lrc14_anchored_deletion_cover_thm4160_anchor_minimality_audit.py
anchor_audit_output: 05-knowledge/results/lrc14_anchored_deletion_cover_thm4160_anchor_minimality_audit.out
carrier_overlap_script: 04-computation/lrc14_anchored_deletion_wrapped_carrier_overlap_thm4160.py
carrier_overlap_output: 05-knowledge/results/lrc14_anchored_deletion_wrapped_carrier_overlap_thm4160.out
carrier_overlap_independent_audit_script: 04-computation/lrc14_anchored_deletion_wrapped_carrier_overlap_thm4160_independent_audit.py
carrier_overlap_independent_audit_output: 05-knowledge/results/lrc14_anchored_deletion_wrapped_carrier_overlap_thm4160_independent_audit.out
script_sha256: 6cafff20874f4854bcf6f95303a607e13844a53bc066bfaca6ae2d9b150feea8
output_sha256: cc80874f70b021abc25fb6fc90ba2e9fe8b4e3b4aca82d00b09cfe84011c387f
semantic_sha256: 7306501ecd70b0c4796963ea61c67cd9b2ec55e5189dcbcee83dc0c8ee0f4e3d
independent_audit_script_sha256: 7829055606a1dcb2d12f3d5a0315d1be8bb1cd1e281473ebe53930268f4a2280
independent_audit_output_sha256: c55ba8d63608bbd56fee419d7afd3ccd3ebd7dc3a1dc841464b0305c7e36c7cf
independent_singleton_semantic_sha256: 91fac18ab0ec1a7e0a9b679c4859963c4905252e3e9d9158760ee7d3f05ad7a8
independent_pair_semantic_sha256: 87e7775750729c679a0c87cb8771c9b602a4f7e8b29683d4d408919116d8f2cc
integer_audit_source_sha256: 419439d1c2aa7fcf6dd75ced8d4d0880e992b74d83eedb4368bd55f8c2cce266
integer_audit_output_sha256: b1b3aeeccbcb26546bc3732560ad30d8d0f9ca4da2f1170cd1df4f1ec0a31738
content_tower_script_sha256: 35ed92277fb229e349054b4af43b292ae249f4497dd87c53f3a77546cfcaaedf
content_tower_output_sha256: d84d5aa5db0309396485ce81f50f86c47e25cd234d01924433565434bd35f098
content_tower_semantic_sha256: b3c557dc73505fbe46a3c3586612af64e23af1d2956b719ec1f6d219da139a5e
anchor_audit_script_sha256: 59c214e97ed94e8c2f49da4052ada5949125669c3abcc31a5a0fa2fad6e93f62
anchor_audit_output_sha256: af29bca2ed5ec371c97f049e6f95cdbacabb5a093b0b1050629fafebb39a294e
anchor_audit_semantic_sha256: 78c304f62c71a4c08881081a4b94270dd4b53df495bd29fbf9f72f2eb1b995d3
carrier_overlap_script_sha256: 954179cb6a87a254bf6435a6d32f4e7e4a51f05a714cd380622a88375a29d884
carrier_overlap_output_sha256: 9479922dc562f7d873cf3de78622957aeb3bffb252a4d89b93f6d4771b644ff3
carrier_overlap_semantic_sha256: 3e1a60dbe4e29904d8eb342e952c634b5c1565cce8a2dd98be6802504ad65e7d
carrier_overlap_independent_audit_script_sha256: dc94fe602dfed2f3e64b160a117acf8fc26e8ffd63ddfb863ef4187ef606ecf6
carrier_overlap_independent_audit_output_sha256: 86a9746211c33d97e6188f58159ecdc8ee34152a612f7091c78ccd31e21fdbeb
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact Fraction interval intersections reconstruct all 27 leave-one
  banks, the complete rank cutoff, singleton/pair/triple candidate universes,
  every qualifying repair set and strict margin, both hostile controls, the
  pigeonhole thresholds, and the exact family count.
independent_audit: >
  ACCEPT. A no-import global 7,134-wall midpoint classification and centered-
  primitive integer evaluator scan the complete q<=27,816 universe, recover
  573 positive singleton masks and 45 repairwise pair edges, and reject every
  repaired triple. Normal, optimized, and hash-seeded outputs byte-match.
integer_audit: >
  ACCEPT. A structurally independent C++20 interval lattice uses checked
  int64 endpoint products and arbitrary-precision measure comparisons. O0,
  O2, and UBSan builds byte-match; all warnings are errors under
  -Wall -Wextra -Wpedantic.
sidecar_audits: >
  PASS. Literal and grouped-binomial content censuses agree exactly, including
  maximum histograms and bounded-height constants. A separate complete anchor
  census proves the fixed-pool height-143 minimality statement.
carrier_overlap_audits: >
  PASS. A literal 11,100,375-body carrier-mask census and an independent
  inclusion-exclusion/subset-DP audit reproduce the content-one overlap,
  deduplicate all multiple-carrier witnesses, and exhibit a dilation hostile.
  Normal, optimized, and hash-seeded outputs byte-match.
---

# THM-4160 -- anchored Haar deletion covers

**PROVED RELATIVE TO THM-4150/4156 + VERIFIED-EXACT + TWO STRUCTURALLY
INDEPENDENT AUDITS; LRC(14) REMAINS OPEN.**

## 1. Statement

Retain the THM-4156 anchors and pool

```text
A={120,126,143},

P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290},

O=P\A.                                                   (1)
```

Put

```text
N={5,66,182,298,336,340,380,386,528,572}.                (2)
```

> **Theorem.** Let either
>
> ```text
> H=A union K,                         K in binom(O,8),   (3)
> ```
>
> or, for some `q in N`,
>
> ```text
> H=A union {q} union K,               K in binom(O,7).  (4)
> ```
>
> For every positive integer `c` and every two distinct positive odd
> integers `a,b`, there is an `x in R/Z` such that
>
> ```text
> min_(v in 2cH union {a,b}) ||vx|| >= 1/14.             (5)
> ```

There are exactly

```text
binom(27,8)                 = 2,220,075 old bodies,
10 binom(27,7)              = 8,880,300 new bodies,
total                        =11,100,375.                (6)
```

The ten classes in `(4)` are disjoint: `N` is disjoint from `P`, and each
body has one labelled newcomer. Multiplication by the same positive `c` is
injective on body labels. Thus no newcomer can duplicate an anchor or old
body label after common dilation. The physical body labels `2cH` are even,
whereas the tails are odd, so no body--tail collision occurs either.

Every scale-one body in `(3)--(4)` is primitive because it contains `A` and
`gcd(A)=1`; every body contains a multiple of each divisor `2,...,14`. The
common dilation in `(5)` records its content and does not assert that the
dilated body remains primitive.

## 2. The deletion-cover lemma

For a finite positive set `S`, write

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14}.                 (7)
```

Let `A` be a subset of a finite pool `P`, put `O=P\A`, and fix `k`. Let `Q`
be a nonempty finite set disjoint from `P`, with `t=|Q|<=k`, and define

```text
D(Q)={r in O:mu(G_((P union Q)\{r}))>=4/63}.             (8)
```

If `|D(Q)|>k-t`, then every `K subset O` with `|K|=k-t` omits some
`r in D(Q)`: otherwise `D(Q) subset K`. For that omitted repair,

```text
A union Q union K subset (P union Q)\{r},
G_((P union Q)\{r}) subset G_(A union Q union K).        (9)
```

Thus the body in `(9)` has Haar mass at least `4/63`, and THM-4150 supplies
every distinct positive odd tail pair after doubling. The comparison in `(8)`
is nonstrict: THM-4150 includes equality. This is a hereditary quantifier over
**every** `K`, not an average or a choice of `K` after `r`.

For `(1)`, `k=8`. The independent literal census checks all 8,880,300 pairs
`(q,K)` in `(4)` and finds an omitted repair every time. Its minimum numbers
of omitted repairs, in the order `(2)`, are

```text
4,8,2,4,2,1,9,11,5,4.                                  (10)
```

## 3. A finite cutoff from safe-comb discrepancy

Let `g(u)=1_[1/14,13/14]({u})`. Its mean is `6/7`. The periodic primitive
`F'=g-6/7`, normalized by `F(0)=0`, falls from `0` to `-3/49` on
`[0,1/14]`, rises to `3/49` on `[1/14,13/14]`, and returns to zero. Hence

```text
osc(F)=6/49.                                             (11)
```

For every interval `[alpha,beta]`, exact integration gives

```text
mu([alpha,beta] intersect G_q)-(6/7)(beta-alpha)
   ={F(q beta)-F(q alpha)}/q.                            (12)
```

If `U` is a union of `c` intervals of total mass `M`, summing `(12)` yields

```text
|mu(U intersect G_q)-(6/7)M| <= 6c/(49q).               (13)
```

For `r in O`, put `U_r=G_(P\{r})`, let `M_r=mu(U_r)` and let `c_r` be its
component count. Exact reconstruction finds

```text
delta_r=4/63-(6/7)M_r > 0                               (14)
```

for all 27 repairs. By `(13)`, repair at `r` forces

```text
q <= B_r=(6c_r/49)/delta_r.                              (15)
```

The complete descending list is:

|rank|`r`|exact `B_r`|floor|
|---:|---:|---:|---:|
|1|252|`273617391247200/9836515057`|27816|
|2|85|`76221987561720/13439926789`|5671|
|3|88|`88227648606240/15641623853`|5640|
|4|95|`125082235998720/23089795681`|5417|
|5|168|`2931614906220/593415569`|4940|
|6|145|`89344454284800/18162161711`|4919|
|7|193|`623790943776/133195981`|4683|
|8|240|`115701068298816/26358573847`|4389|
|9|264|`10974763495080/2613894587`|4198|
|10|176|`17198807449824/4350517789`|3953|
|11|290|`3566798135901/913013672`|3906|
|12|286|`300979130371920/81651228889`|3686|
|13|8|`19299798132615/5236528393`|3685|
|14|170|`84877231570560/23235046577`|3652|
|15|16|`56144867294880/16218075763`|3461|
|16|42|`38110993780860/11279974859`|3378|
|17|80|`84877231570560/25158741857`|3373|
|18|132|`297070310496960/92205621463`|3221|
|19|190|`601958260743840/191486558657`|3143|
|20|63|`3762239129649/1241068982`|3031|
|21|10|`120391652148768/45713236603`|2633|
|22|15|`37133788812120/14242896053`|2607|
|23|40|`14658074531100/5821304243`|2518|
|24|84|`293161490622000/116882514613`|2508|
|25--27|20,30,60|`293161490622000/118997676883`|2463|

The last three labels give the only tie. If `|Q|=t`, every `q in Q` must by
monotonicity have at least `9-t` singleton repairs before the deletion-cover
inequality can hold. The `(9-t)`-th largest bound therefore gives

```text
t:       1     2     3     4     5     6     7      8
q<=   4389  4683  4919  4940  5417  5640  5671  27816. (16)
```

Eight repairs imply directly that the integer `q<=4389`. The hostile
`q=4390` can have at most seven by `(15)`; exact interval intersection finds
zero. Since `B_252<27817`, every `q>=27817` has no repair at all.

## 4. Complete singleton classification

The exact universe for one newcomer is all 4,359 integers

```text
1<=q<=4389,                    q notin P.                (17)
```

Exactly the ten labels in `(2)` have at least eight repairs. Their full repair
sets and smallest strict Haar surpluses are:

|`q`|`|D({q})|`|`D({q})`|minimum surplus over `4/63`|
|---:|---:|:---|:---|
|5|11|`85,88,95,145,168,176,193,240,252,264,290`|`1228523/3265513680`|
|66|15|`42,80,85,88,95,145,168,170,176,193,240,252,264,286,290`|`261933509/4560289854120`|
|182|9|`85,95,145,168,176,193,240,252,290`|`504443/38976836360`|
|298|11|`85,88,95,145,168,176,193,240,252,264,290`|`1762326833/113247198043980`|
|336|9|`85,88,95,145,168,176,193,252,264`|`58030039/1520096618040`|
|340|8|`85,88,95,145,168,240,252,264`|`2854953/7795367272`|
|380|16|`8,16,85,88,95,132,145,168,170,176,193,240,252,264,286,290`|`273319799/1302939958320`|
|386|18|`8,16,42,85,88,95,132,145,168,170,176,190,193,240,252,264,286,290`|`35335481/760048309020`|
|528|12|`8,85,88,95,145,168,170,193,240,252,286,290`|`52627243/138190601640`|
|572|11|`85,88,95,145,168,176,193,240,252,264,290`|`426566897/829143609840`|

There is no equality repair. The cardinality boundary is `q=340`, with
exactly eight repairs. The strongest label is `q=386`, with 18. Its weakest
repair is `r=190`, where

```text
mu(G_((P union {386})\{190}))
 =144877112923/2280144927060
 =4/63+35335481/760048309020.                            (18)
```

The clean just-failing hostile is

```text
q=235,      D({235})={85,88,145,168,193,240,252},       (19)
```

with seven repairs. An exact `>=` audit also finds no `q` addable to the
entire pool `P`: such a `q` would repair all 27 deletions, hence lie in `(17)`,
but none does. There is no hidden equality case.

Applying Section 2 with `t=1` proves `(4)--(5)` at content one; `(3)` is
THM-4156. Finally

```text
G_(cH)=m_c^(-1)(G_H),                 m_c(y)=cy,          (20)
```

and the surjective circle endomorphism `m_c` preserves Haar measure. Therefore
THM-4150 applies at every positive common body content. **QED.**

## 5. Exact obstruction to stacking newcomers in this mechanism

For pairs, every member must have seven singleton repairs. The complete
candidate set is

```text
{5,66,182,235,298,336,340,380,386,528,572}.             (21)
```

Among its 55 pairs, 27 have zero common repairs and 28 have one; every
nonzero repair is `r=252`. No pair reaches the required seven. The complete
repairwise scan over all positive labels is stronger:

```text
max{|Q|:Q nonempty, r in D(Q)}
                    = 0  for r in {10,15,20,30,40,60,63,84},
                      2  for r=252,
                      1  for every other r in O.         (22)
```

Exactly 45 unordered pairs survive any repair, all solely at `r=252`. A
strict positive control is

```text
Q={5,66}, r=252,
mu(G_((P union Q)\{252}))=347277745/5333672344>4/63.     (23)
```

No triple survives even one repair. The independent pair graph has 46
triangles, and direct exact intersection rejects all 46, without equality.
Equivalently the checked-integer search tests all 19,768 monotone triple
extensions of the 45 surviving pair branches and rejects them all. Hence

```text
max_(|Q|=1)|D(Q)|=18,
max_(|Q|=2)|D(Q)|=1,
D(Q)=empty for every |Q|>=3.                             (24)
```

The thresholds for `|Q|=1,...,8` are `8,7,...,1`, so `(24)` proves that the
ten singleton classes are the **complete classification of this one-deletion-
cover mechanism**. This does not obstruct two simultaneous old deletions,
body-specific geometry, or another Haar lower bound.

## 6. Content-stable separation from the named gates

For a body `H`, put

```text
D_H=16 max(H)-156 min(H).                               (25)
```

At common content `c`, THM-4151's affine first-window gate is `cD_H<=13`.
A literal census and a different grouped min/max binomial reconstruction agree:

|family|`D_H<=0`|`1<=D_H<=13`|`D_H>=14`|
|:---|---:|---:|---:|
|old THM-4156|344,366|0|1,875,709|
|ten new classes|1,052,735|0|7,827,565|
|total|1,397,101|0|9,703,274|

For the new family, the largest nonpositive defect is `-20` and the smallest
positive defect is `88`; for the old family they are `-20` and `192`.
Therefore exactly `9,703,274` bodies remain beyond THM-4151 at **every**
positive content. All 11,100,375 bodies also fail the stated THM-4148 gate,
so the same 9,703,274 fail both named gates at every content. This compares
only those two declared sufficient criteria.

Every body here omits the label `7`, so this family is disjoint from
THM-4158's **named `m=7` specialization**, whose anchor set requires `7`.
The full wrapped-carrier overlap can also be classified at content one.

**FINITE-EXACT content-one overlap.** Put

```text
C={H:there exists m>=1 with H subset P_m},                (25a)
```

where `P_m` is THM-4158's exact alphabet, and call `H` canonical when
`H subset P_min(H)`. Exact literal carrier masks and an independent
inclusion-exclusion/subset-DP route agree:

|content-one family|canonical|in `C`|outside `C`|
|:---|---:|---:|---:|
|old THM-4156|422,222|451,581|1,768,494|
|ten new classes|1,298,364|1,612,376|7,267,924|
|total|1,720,586|2,063,957|9,036,418|

Every core contains `120`, while `min(P_m)=m`; hence `m<=120` is an
exhaustive carrier search. Intersecting one 120-bit membership mask per label
counts each body once even when several `m` witness it. There are `1,664,667`
such multiple-witness bodies, so summing carrier incidences would overcount.

All `1,397,101` nonpositive-defect bodies are canonical and lie in `C`.
Among the `9,703,274` positive-defect bodies, `666,856` lie in `C` and the
remaining `9,036,418` lie outside every content-one wrapped carrier. Thus the
last number is a genuine content-one separation from all THM-4158 alphabets,
not merely from its named `m=7` specialization.

This separation is deliberately **not** promoted through the content tower.
For the hostile core

```text
H={20,30,40,42,60,63,120,126,143,168,264},               (25b)
```

`H` lies in no `P_m`, whereas `2H subset P_35`. The all-content statement in
the first half of this section therefore remains only the separation from the
named THM-4148/4151 gates.

The output freezes all maximum histograms. If height means `max(cH)<=X`,
primitivity of the cores makes content representations unique, and

```text
#{cH:max(cH)<=X}=C X+O(11,100,375),
C=sum_H 1/max(H).                                       (26)
```

The exact constants for all/stable/beyond bodies are respectively

```text
2300406906171178901/64712684596560,
25289558064375623/5392723716380,
23493320110572605/761325701136.                         (27)
```

## 7. Anchor-lowering boundary

The complete finite universe

```text
1<=a<b<c<=143, gcd(a,b,c)=1,
each d=2,...,14 divides at least one of a,b,c            (28)
```

contains exactly

```text
(70,72,143), (72,140,143), (120,126,143).               (29)
```

Thus maximum anchor `143` is forced. Keeping `O` fixed, the corresponding
30-label pool masses are

```text
3242643641/80005085160,
90567155849/2280144927060,
298133356159/4560289854120,                             (30)
```

and only the last reaches `4/63`. Hence `A` is the unique primitive divisor-
complete triple of height at most 143 compatible with this fixed optional
pool's common-Haar threshold. This is bounded fixed-pool minimality, not
global optimality.

## 8. Audit boundary and replay

The primary Python path builds `U_r` by exact safe-tooth intersections. The
independent Python path classifies one common 7,134-wall arrangement by exact
midpoints, recovers all 27 leave-one banks simultaneously, scans all
`q<=27816`, and builds exact repairwise pair graphs.

The C++20 path reconstructs intervals again and uses arbitrary-precision
integer measure comparisons. Every rational endpoint is at most
`14*27816=389424` before reduction; every int64 multiplication is checked.
The theoretical endpoint cross-product bound `151651051776` is far below
`2^63-1`; the observed maximum is printed in the certificate.

Primary replay:

```bash
python3 -B 04-computation/lrc14_anchored_deletion_cover_thm4160.py
python3 -B -O 04-computation/lrc14_anchored_deletion_cover_thm4160.py
PYTHONHASHSEED=271828 python3 -B 04-computation/lrc14_anchored_deletion_cover_thm4160.py
```

Independent Fraction replay:

```bash
python3 -B 04-computation/lrc14_anchored_deletion_cover_thm4160_independent_audit.py
python3 -B -O 04-computation/lrc14_anchored_deletion_cover_thm4160_independent_audit.py
PYTHONHASHSEED=314159 python3 -B 04-computation/lrc14_anchored_deletion_cover_thm4160_independent_audit.py
```

Checked-integer replay:

```bash
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_anchored_deletion_cover_thm4160_independent_audit.cpp \
  -o /tmp/lrc4160-o2
/tmp/lrc4160-o2
g++ -std=c++20 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_anchored_deletion_cover_thm4160_independent_audit.cpp \
  -o /tmp/lrc4160-o0
/tmp/lrc4160-o0
g++ -std=c++20 -O1 -g -fsanitize=undefined \
  -fno-sanitize-recover=undefined -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_anchored_deletion_cover_thm4160_independent_audit.cpp \
  -o /tmp/lrc4160-ubsan
/tmp/lrc4160-ubsan
```

Sidecar replay:

```bash
python3 -B 04-computation/lrc14_anchored_deletion_cover_thm4160_content_tower_audit.py
python3 -B 04-computation/lrc14_anchored_deletion_cover_thm4160_anchor_minimality_audit.py
```

Content-one wrapped-carrier overlap replay:

```bash
python3 -B 04-computation/lrc14_anchored_deletion_wrapped_carrier_overlap_thm4160.py
python3 -B 04-computation/lrc14_anchored_deletion_wrapped_carrier_overlap_thm4160_independent_audit.py
```

The source-to-target contract is

```text
source:       one-deletion Haar repair sets on the THM-4156 pool
target:       anchored eleven-bodies with every distinct positive odd pair
map:          omitted repair -> safe-set inclusion -> THM-4150 transfer
preserved:    exact Haar threshold, anchors, divisor completeness, content
destroyed:    body-specific extra safe phases and component addresses
sidecar:      safe-comb discrepancy, repair graph, carrier masks, content census
positive:     q=386 has 18 repairs; {5,66} survives at r=252
hostile:      q=235 has only seven singleton repairs; no triple has any
decisive test: complete finite cutoff plus two independent exact geometries. (31)
```

This theorem strengthens a structured dyadic odd-tail family only. It does not
enter the dyadic seam for arbitrary divisor-complete bodies, treat mixed/even
tails, prove maximality beyond this one-deletion mechanism, or prove LRC(14).
