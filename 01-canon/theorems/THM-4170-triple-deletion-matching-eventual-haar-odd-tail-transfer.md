---
id: THM-4170
title: "Triple-deletion matching eventual Haar odd-tail transfer"
status: >
  PROVED RELATIVE TO THM-4150/4156 + VERIFIED-EXACT FRACTION/INTEGER FINITE
  CENSUS + INDEPENDENT CHECKED-INTEGER C++ AUDITS; LRC(14) OPEN. Every
  positive newcomer q outside P except exactly 61 labels passes the
  triple-deletion transversal certificate; the last exception is q=924.
  Hence every q>=925 closes all 888,030 anchored seven-subset bodies with
  every distinct positive odd tail pair at every common body content.
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
  - THM-4160-anchored-haar-deletion-cover-and-content-tower
  - THM-4166-two-deletion-repair-graph-haar-odd-tail-transfer
  - THM-4172-multideletion-support-tomography-and-same-parity-johnson-holonomy
script: 04-computation/lrc14_triple_deletion_matching_eventual_haar_thm4170.py
output: 05-knowledge/results/lrc14_triple_deletion_matching_eventual_haar_thm4170.out
independent_audit_source: 04-computation/lrc14_triple_deletion_matching_eventual_haar_thm4170_independent_audit.cpp
independent_audit_output: 05-knowledge/results/lrc14_triple_deletion_matching_eventual_haar_thm4170_independent_audit_cpp.out
finite_census_script: 04-computation/lrc14_triple_deletion_matching_eventual_haar_thm4170_finite_census.py
finite_census_output: 05-knowledge/results/lrc14_triple_deletion_matching_eventual_haar_thm4170_finite_census.out
finite_census_independent_source: 04-computation/lrc14_triple_deletion_matching_eventual_haar_thm4170_finite_census_independent.cpp
finite_census_independent_output: 05-knowledge/results/lrc14_triple_deletion_matching_eventual_haar_thm4170_finite_census_independent_cpp.out
script_sha256: d6c7413e94a4e92a096e96eb7586287f58a033a58f7a0fb27fdd43e34cb49c25
output_sha256: 4a48470d4af3b4bc2e35874ddacaa42ebcf296670680521252f4127e411f1ab1
semantic_sha256: d6ae7826c902ee3dd166a5b0961d29b3f8fd3829ca1c1d256a6677a24da46024
independent_audit_source_sha256: 0fed5f5e2b8968309b86f4196901ff4396b99a576853fabfa11803e1132994ee
independent_audit_output_sha256: 91f19492b174d78747f5fcf3d65fe0af60251651831c17abd507b5b9b53faf01
finite_census_script_sha256: 932847febcaa4aea0f564fe458c82b52fa02e5bdc693b55f14ab8a1d061e5709
finite_census_output_sha256: 86d38ff3cdedb7414852ceaa736e35216ac7b2f17ed33c19cdb535b746bff80d
finite_census_independent_source_sha256: 85395ddef97adcfc6aed3b0c9485b517da06e0bc2417561156f3a8e8791a15ea
finite_census_independent_output_sha256: 046b2952e253bdbbf43c37ac03921b51e9c28af6dbb017eda1fc5acc3b4eb63b
finite_qualifier_q_word_xor_mul64: 02784121a66537ac
hash_basis: raw LF bytes
primary_audit: >
  PASS. Exact Fraction arithmetic constructs the 7,134 global walls,
  classifies every cell by its optional failure mask, reconstructs all 2,925
  triple-deletion banks, and independently integrates the q=9,699/9,700
  controls. Normal, optimized, and hash-seeded outputs byte-match.
independent_audit: >
  ACCEPT. A C++20 common-denominator cell ledger uses signed 128-bit measure
  comparisons, exhausts every strict-limit and boundary hyperedge, and
  independently rejects every seven-vertex cover at the two controls. O0,
  O2, and UBSan outputs byte-match under warnings-as-errors builds.
finite_census_audits: >
  PASS/ACCEPT. A vectorized exact common-lattice Python census and a separate
  C++20 grouped-cell implementation classify all 9,669 newcomer labels below
  9,700. They agree on all 61 certificate failures, every minimum-transversal
  histogram entry, zero equalities, q=924/925 controls, and the specified
  wordwise XOR/multiply qualifier ledger. Python normal/-O/seed and C++
  O0/O2/UBSan outputs
  respectively byte-match.
portability_audit: >
  PASS. Both Python paths force LF stdout, both C++ paths force binary stdout
  on Windows, and exact-path attributes keep both C++ sources LF. Fresh
  MinGW C++20 builds and Python replays preserve all four frozen output hashes.
---

# THM-4170 -- triple-deletion matching eventual Haar transfer

**PROVED RELATIVE TO THM-4150/4156 + TWO CROSS-CHECKING ANALYTIC AUDITS +
TWO CROSS-CHECKING FINITE CENSUSES; LRC(14) REMAINS OPEN.**

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
B={3,6,22,24,25,46,48,50,55,64,70,72,75,83,93,96,100,
   103,105,110,122,127,128,140,147,153,158,166,172,173,
   183,186,192,206,210,220,256,260,270,282,294,306,313,
   320,325,332,346,366,372,384,416,440,462,512,516,520,
   550,567,744,768,924}.                               (2)
```

For a positive integer `q notin P union B` and `K in binom(O,7)`, put

```text
H_(q,K)=A union {q} union K.                             (3)
```

> **Theorem.** For every positive integer `c`, every two distinct positive
> odd integers `a,b`, every positive `q notin P union B`, and every
> `K in binom(O,7)`,
> there is an `x in R/Z` such that
>
> ```text
> min_(v in 2cH_(q,K) union {a,b}) ||vx|| >= 1/14.       (4)
> ```

In particular, `(4)` holds for every integer `q>=925`. The set `B` is the
exact failure set of the triple-deletion transversal **certificate**, not a
set of bodies proved unsafe.

For each `q`, this gives exactly

```text
binom(27,7)=888,030                                      (5)
```

primitive divisor-complete core bodies. Thus `(4)` is an infinite labelled
family. The cofinal assertion combines a complete finite census with a
proved discrepancy tail; it is not a finite extrapolation.

## 2. The deletion-hypergraph lemma

For a finite positive set `S`, write

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14}.                 (6)
```

Let `A subset P`, put `O=P\A`, let `Q` be finite and disjoint from `P`, and
fix integers `d,s>=0`. Define the `d`-uniform repair hypergraph on `O` by

```text
E_d(Q)={R in binom(O,d):
        mu(G_((P union Q)\R))>=4/63}.                    (7)
```

If `tau(E_d(Q))>s`, then every `K in binom(O,s)` misses an entire repair
edge: `K` is not a transversal, so some `R in E_d(Q)` has `R intersect
K=empty`. Consequently

```text
A union Q union K subset (P union Q)\R,
G_((P union Q)\R) subset G_(A union Q union K).          (8)
```

The target safe set has Haar mass at least `4/63`, and THM-4150 closes every
distinct positive odd tail pair after doubling. A matching of `s+1` repair
edges is a transparent sufficient certificate: an `s`-set can meet at most
`s` pairwise-disjoint edges. This is sufficient, not an equivalence.

Repair certification is also monotone upward in deletion arity. If every
`s`-body `K` misses some `d`-repair `R` and `d<d'<=|O|-s`, extend `R` inside
`O\K` to a `d'`-set `R'`. Then `(P union Q)\R'` is a subset of
`(P union Q)\R`, so safe-set inclusion reverses and `R'` remains a repair
disjoint from `K`. This is a body-by-body implication; it need not produce
one fixed matching at arity `d'`.

## 3. Strict-limit triple banks

Let `g(u)=1_[1/14,13/14]({u})`. Its mean is `6/7`. The normalized periodic
primitive of `g-6/7` has oscillation `6/49`; hence, if `U` is a union of
`c_U` intervals with mass `M_U`, then

```text
|mu(U intersect G_{q})-(6/7)M_U| <= 6c_U/(49q).         (9)
```

For `R in binom(O,3)`, set

```text
U_R=G_(P\R),   M_R=mu(U_R),   c_R=#components(U_R),
s_R=(6/7)M_R-4/63.                                    (10)
```

Because `G_((P union {q})\R)=U_R intersect G_{q}`, every strict-limit triple
with `s_R>0` is a repair whenever

```text
q>=B_R=(6c_R/49)/s_R.                                  (11)
```

The complete exact universe has `2,925` triple banks, `1,335` strict-limit
triples, and zero limit equalities `s_R=0`. The theorem uses the following
eight pairwise-disjoint triples.

|`R`|`c_R`|`M_R`|`s_R`|`B_R`|`ceil(B_R)`|
|:---|---:|:---|:---|:---|---:|
|`(8,84,252)`|148|`53623011809/701583054480`|`4961689987/2455540690680`|`44500410884160/4961689987`|8969|
|`(15,88,95)`|178|`140987531863/1824115941648`|`17603497445/6384405795768`|`139153987548576/17603497445`|7905|
|`(40,85,170)`|166|`4641131947/59611632080`|`6087298409/1877766410520`|`38168476426080/6087298409`|6271|
|`(63,176,193)`|164|`23261449/302928780`|`7400521/3180752190`|`63874697040/7400521`|8632|
|`(10,145,290)`|162|`12058235549/157251374280`|`1229956807/550379809980`|`10917738271440/1229956807`|8877|
|`(42,132,264)`|156|`354200884867/4560289854120`|`49204909241/15961014489420`|`304887950246880/49204909241`|6197|
|`(16,168,286)`|168|`4772688913/62044759920`|`1591026937/651469979160`|`13401668142720/1591026937`|8424|
|`(80,190,240)`|164|`13415915189/175395763620`|`1270909207/613885172670`|`12327816528720/1270909207`|9700|

They use 24 distinct optional labels, leaving only `20,30,60`. Their largest
activation ceiling is `9700`, so all eight are repair edges for every
`q>=9700`.

The number `9700` is sharp only for these eight displayed **discrepancy
inequalities**: the last `B_R` lies strictly between `9699` and `9700`. It is
not claimed to be the first safe newcomer, the first exact matching, or a
global minimax over other matchings. Indeed, the direct exact `q=9699`
control is already positive.

## 4. Complete finite census below the discrepancy tail

For every positive `q<9700` outside `P`, both audits construct the complete
hypergraph `E_3({q})` and decide whether it has a cover of size at most seven.
The explicit universe has `9,669` labels. The exact result is

```text
tau(E_3({q}))>7:  9,608 labels,
tau(E_3({q}))<=7:    61 labels, exactly the set B in (2),
threshold equalities: 0.                                (12)
```

The minimum-transversal histogram on the 61 failures is

|`tau`|0|1|2|3|4|5|6|7|
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|count|4|6|5|8|5|9|10|14|

Of the 9,608 qualifiers, 8,837 have the immediate lexicographic matching-of-
eight certificate; the remaining 771 pass the exhaustive no-seven-cover
branch. The ordered qualifying `q` ledger has the reproducibility sidecar

```text
h_0=0xcbf29ce484222325,
h_(j+1)=((h_j XOR q_j)*0x100000001b3) mod 2^64,
h_final=02784121a66537ac.                                (13)
```

Here the qualifying `q_j` are taken in increasing order and each is one
unsigned integer word. This is a wordwise XOR/multiply recurrence, not
standard bytewise FNV-1a-64 (MISTAKE-513).

The last certificate failure and the next label form a sharp finite control:

|`q`|hyperedges|greedy matching|exact result|
|---:|---:|---:|:---|
|924|258|4|`tau=6`, cover `{16,85,88,145,168,252}`|
|925|1705|9|`tau>7`|

At `q=925`, the first eight greedy edges are

```text
(8,10,85), (15,16,88), (20,30,252), (40,80,95),
(42,63,145), (60,132,264), (84,168,170), (176,190,193).
                                                               (14)
```

They are pairwise disjoint and every exact Haar margin is strictly positive;
the frozen outputs give all eight fractions. Separate direct Fraction
integration reproduces both the `q=924` cover ledger and the `q=925` matching
measures. Thus `924` is the last failure of this **triple-deletion
certificate**, not a body proved unsafe.

## 5. Pigeonhole closure, content, and uniqueness

Fix a positive `q notin P union B` and `K in binom(O,7)`. If `q<9700`, the
finite census gives `tau(E_3({q}))>7`; if `q>=9700`, the displayed eventual
matching gives the same conclusion. In either case the deletion-hypergraph
lemma supplies a repair `R` disjoint from `K`, and hence

```text
mu(G_(H_(q,K)))>=4/63.                                  (15)
```

THM-4150 proves `(4)` for `c=1`. For every positive integer `c`, the
surjective circle endomorphism `m_c(y)=cy` gives

```text
G_(cH_(q,K))=m_c^(-1)(G_(H_(q,K))),
mu(G_(cH_(q,K)))=mu(G_(H_(q,K))).                       (16)
```

THM-4150 applied to `cH_(q,K)` proves `(4)` at every common content.

Every core is primitive because `gcd(120,126,143)=1`, and the anchors contain
a multiple of every divisor `2,...,14` as recorded in THM-4156. Because
`q notin P`, every core has exactly one newcomer; equality of two cores first
forces equal newcomers and then equal optional subsets. Across contents,
`gcd(2cH_(q,K))=2c`, so equal physical body sets first force equal contents
and then equal cores. Even body labels cannot collide with odd tails. This
proves `(4)--(5)` and the uniqueness claims. **QED.**

The arity monotonicity gives a clean corollary: the same bodies are certified
by every `d`-deletion repair hypergraph with `4<=d<=20`. After finding a
disjoint triple, extend it inside the 20-label complement `O\K`. This is a
second certificate for the same `888,030` bodies per `q`, not an additional
family count.

There is also an eventual separation from THM-4158 at **content one**. Every
core `H_(q,K)` contains `120`. If it were contained in a THM-4158 alphabet
`P_m`, then `m=1` is impossible because `max(P_1)=44`, while for `m>=2` the
first label of `P_m` is `m`, so `m<=120`. The third-band endpoint is increasing
in `m`, and therefore

```text
max(P_m)<=floor(41(12*120+1)/16)=3692.                 (16a)
```

Consequently, for every `q>=3693` and every `K in binom(O,7)`, the core
`H_(q,K)` lies in no THM-4158 alphabet. Since every such `q` qualifies by the
theorem, this gives `888,030` content-one cores per newcomer outside the full
union of THM-4158 carrier alphabets. This is a mechanism-level alphabet
comparison only. It does not claim disjointness after a common dilation and
does not produce additional LRC bodies beyond `(5)`. THM-4174 later proves
that the same `q>=3693` separation holds at every common content; that stronger
statement is not used in this theorem.

## 6. Exact controls and audit architecture

The common arrangement has `7,134` rational walls and `7,133` open cells.
Classifying cells by optional failures gives

```text
failure sizes 0,1,2,3,>3-or-anchor: 150,328,518,678,5459. (17)
```

The Fraction audit reconstructs all banks from failure buckets and separately
integrates the safe primitive at the boundary:

|`q`|repair hyperedges|equalities|chosen eight active?|
|---:|---:|---:|:---:|
|9699|1351|0|yes|
|9700|1356|0|yes|

The C++ audit embeds the cells in common denominator `18241159416480`, uses
signed integer prefix differences, and tests all `2,925` triples. It finds no
seven-vertex cover at either boundary. Reduced pool endpoints are at most
`14*290=4060`; the unreduced midpoint product bound `32967200` fits signed
64-bit arithmetic, and Haar comparisons use signed 128-bit integers.

The activation ledger has 665 strict-limit triples with discrepancy ceiling
at most `9699` and 666 at most `9700`; the previous ceiling event is `9687`.
These are side data only. No global matching-minimax claim is made.

The finite Python audit uses a different exact representation. It embeds the
Fraction-built walls in the same lattice, vectorizes the `9,669 x 2,925`
measure comparisons in bounded signed-int64 blocks, and then runs a separate
recursive transversal search. Direct Fraction primitives rebuild `q=924`
and `q=925`. The independent C++ audit instead streams each newcomer comb
through the 7,133 cells and accumulates base, 27 singleton, 351 pair, and
2,925 triple failure buckets. Both use the exact identity

```text
mass(R,q)=bucket(empty)+sum_(empty!=F subset R) bucket(F,q), (18)
```

then compare by signed integer cross-multiplication. Both freeze the complete
61-label failure ledger and a valid minimum cover for every failure.

Primary replay:

```bash
python3 -B 04-computation/lrc14_triple_deletion_matching_eventual_haar_thm4170.py
python3 -B -O 04-computation/lrc14_triple_deletion_matching_eventual_haar_thm4170.py
PYTHONHASHSEED=4170 python3 -B 04-computation/lrc14_triple_deletion_matching_eventual_haar_thm4170.py
```

Independent replay:

```bash
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_triple_deletion_matching_eventual_haar_thm4170_independent_audit.cpp \
  -o /tmp/lrc4170-o2
/tmp/lrc4170-o2
g++ -std=c++20 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_triple_deletion_matching_eventual_haar_thm4170_independent_audit.cpp \
  -o /tmp/lrc4170-o0
/tmp/lrc4170-o0
g++ -std=c++20 -O1 -g -fsanitize=undefined \
  -fno-sanitize-recover=undefined -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_triple_deletion_matching_eventual_haar_thm4170_independent_audit.cpp \
  -o /tmp/lrc4170-ubsan
/tmp/lrc4170-ubsan
```

Finite-census replay:

```bash
python3 -B 04-computation/lrc14_triple_deletion_matching_eventual_haar_thm4170_finite_census.py
python3 -B -O 04-computation/lrc14_triple_deletion_matching_eventual_haar_thm4170_finite_census.py
PYTHONHASHSEED=925 python3 -B 04-computation/lrc14_triple_deletion_matching_eventual_haar_thm4170_finite_census.py
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_triple_deletion_matching_eventual_haar_thm4170_finite_census_independent.cpp \
  -o /tmp/lrc4170-finite-o2
/tmp/lrc4170-finite-o2
g++ -std=c++20 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_triple_deletion_matching_eventual_haar_thm4170_finite_census_independent.cpp \
  -o /tmp/lrc4170-finite-o0
/tmp/lrc4170-finite-o0
g++ -std=c++20 -O1 -g -fsanitize=undefined \
  -fno-sanitize-recover=undefined -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_triple_deletion_matching_eventual_haar_thm4170_finite_census_independent.cpp \
  -o /tmp/lrc4170-finite-ubsan
/tmp/lrc4170-finite-ubsan
```

## 7. Scope boundary

THM-4166's two-deletion graph has last qualifier `q=8265`. Triple deletion
instead has only 61 exceptional newcomer labels and no exception after
`q=924`; the eventual matching supplies the proof beyond the finite census.
The contract is

```text
source:       strict-limit triple-deletion banks of the THM-4156 pool
target:       anchored eleven-bodies with every distinct positive odd pair
map:          tau>7 or eight disjoint repairs -> missed triple -> transfer
preserved:    exact Haar threshold, anchors, divisor completeness, content
destroyed:    component addresses and body-specific extra safe phases
sidecar:      safe-comb discrepancy and exact boundary hypergraphs
positive:     q=925 has nine greedy repair triples; all q>=925 qualify
hostile:      q=924 has the six-cover {16,85,88,145,168,252}
decisive test: complete q<9700 census plus the eight-edge discrepancy tail. (19)
```

The 61 labels in `B` are exactly the failures of this triple-deletion
transversal mechanism, not unsafe bodies and not exclusions from every known
mechanism. Labels in `P` are outside the newcomer definition. This theorem
does not assert maximality of the matching, treat arbitrary divisor-complete
bodies, treat mixed/even tail parities, or prove LRC(14).

THM-4172 supplies a cross-thread deletion-tomography warning. Symmetric
deletion layers can invert an outside-support-size histogram, but that scalar
histogram does not retain which repair triples overlap. Here the preserved
target is `tau(E_3({q}))>7`, so the missing coordinate is the labelled
triple-incidence deck (or equivalent overlap data), not another arity count.
The cheapest decisive hostile is a pair of hypergraphs with the same
support-size/deletion-layer histogram and different transversal number.
THM-4174 subsequently performs the exact `d=4,5,6` residual census and closes
the fixed-pool newcomer filtration, while confirming that tomography cannot
replace cover recursion without this incidence sidecar.
