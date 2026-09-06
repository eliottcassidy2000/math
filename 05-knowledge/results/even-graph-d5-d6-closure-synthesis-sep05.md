# Cumulative signed-cycle gaps at D=5 and D=6

**Status: PROVED CANDIDATE, pending final independent proof review and canon
promotion.** Both finite bases and the local comparison lemmas have passed
two independent exhaustive exact implementations. This note supplies the
all-order induction, including equality. It makes no external priority claim.

## Inheritance

The closest proved mechanism is
[THM-4083, cumulative D3/D4 spectral gap](../../01-canon/theorems/THM-4083-even-graph-cumulative-d3-d4-spectral-gap.md):
balanced-deletion rigidity plus weighted deletion. The canonical hostile is
the antibalanced signing, for which every even cycle is positive. Consequently
a lower bound on each individual long layer cannot work. MISTAKE-496 is the
corrected near miss: the trivial balanced character must be excluded from a
spectral-gap minimum. The least-used sidecar is the full local vector of
negative-cycle counts, already available in THM-4078's exact Fourier atlas.

[THM-4200, four-edge D5 firewall](../../01-canon/theorems/THM-4200-even-graph-four-edge-d5-frustration-firewall.md)
and [THM-4084, matching profile](../../01-canon/theorems/THM-4084-even-graph-matching-character-profile-and-d5-firewall.md)
restrict possible D5 obstructions to frustration at least five. The present
proof eliminates every remaining frustration at once; these restrictions are
useful provenance, not assumptions in the new theorem.

## Statement and normalization

For a signing of the complete graph K_n let c_k be the number of unoriented
simple k-cycles whose edge-sign product is -1. Switching changes signs along
a cut and preserves every c_k. The balanced class is characterized by c_3=0.
Set

```text
S_D = sum_(k=3)^(D+1) c_k,
A_(n,D) = sum_(k=3)^(D+1) (n-2)!/(n-k)!.
```

For D=5 and every n>=6, and for D=6 and every n>=7,

```text
min_(nonbalanced switching classes) S_D = A_(n,D).
```

In both cases equality holds exactly on the switching classes represented
by one negative edge. There are binom(n,2) such labelled classes, forming
one relabelling orbit. Explicitly,

```text
A_(n,5) = (n-2)(n^3-11n^2+41n-50),
A_(n,6) = (n-2)(n^4-17n^3+108n^2-301n+310).
```

Under [THM-4078, triangle-quotient spectrum](../../01-canon/theorems/THM-4078-even-graph-triangle-quotient-spectrum-and-boolean-noncommutation.md),
the unnormalized cumulative cycle Laplacian therefore has gap 2A_(n,D).
Its full labelled multiplicity is binom(n,2), and its relabelling-invariant
quotient multiplicity is one. If Q=sum_(k=3)^(D+1) n!/(2k(n-k)!) is the
degree, the normalized algebraic gap is 2A/Q and the lazy gap is A/Q.
These statements use the multiplicity-weighted cycle operator, not a
Booleanized quotient adjacency matrix.

## Local comparisons retain the missing layer

The following two finite statements include the balanced character:

```text
on K_6: 3c_6 >= 2c_4;                         (L6)
on K_7: 2c_7 >= c_4+c_5.                      (L7)
```

They are FINITE-EXACT lemmas, proved by exhausting all 2^10=1,024 and
2^15=32,768 labelled switching classes. Gauge every edge at vertex zero
positive. Every other edge sign is independent, and each switching class
has exactly one such representative. No frustration, degree, symmetry,
connectivity, or positivity filter is imposed.

For transparency the complete K6 profile table is below. The last column is
the number of labelled root-gauged signings, not an isomorphism count.

| c3 | c4 | c5 | c6 | multiplicity |
|---:|---:|---:|---:|---:|
|0|0|0|0|1|
|4|12|24|24|15|
|6|18|36|36|60|
|8|20|32|24|45|
|8|24|40|32|180|
|8|24|48|32|15|
|10|18|36|36|20|
|10|22|36|28|180|
|10|26|36|28|180|
|10|30|36|20|12|
|12|20|40|24|45|
|12|24|24|32|15|
|12|24|32|32|180|
|14|18|36|36|60|
|16|12|48|24|15|
|20|0|72|0|1|

Thus L6 is sharp even away from zero. One equality witness is a negative
5-cycle on five vertices and a positive apex, with profile (10,30,36,20).
Antibalanced and balanced signings both have c4=c6=0, as required.
For L7 the minimum of 2c7-c4-c5 over all nonbalanced signings is 108;
balanced gives zero. Only the weak inequality is used below.

Sum L6 over every induced six-set of a signed K_n. Each negative four-cycle
appears in binom(n-4,2) such sets and each six-cycle appears once. Hence

```text
3c6 >= (n-4)(n-5)c4,                         n>=6.    (1)
```

In particular c6>=c4 for n>=7. Similarly, summing L7 over seven-sets gives

```text
2c7 >= binom(n-4,3)c4 + binom(n-5,2)c5,      n>=7.    (2)
```

For n>=8 the two binomial coefficients are at least 4 and 3, so (2) gives
c7>=2c4+c5. The local comparisons are applied with their exact occurrence
multiplicities; no independence assumption is involved.

## Balanced deletions and the only remaining case

Write b for the number of vertices whose deletion leaves a balanced signing.
The elementary trichotomy in THM-4083 Sections 3 and 4 says that for any
nonbalanced signing of K_n with n>=4:

- If b>=2, its switching class is a single edge.
- If b=1, it has a negative star representative with part sizes r,s>=2,
  r+s=n-1, and S_D=rs A_(n,D)/(n-2)>A_(n,D).
- If b=0, every vertex deletion is nonbalanced.

For completeness, if two deletions u,v are balanced, every negative triangle
contains u,v. The product of the four triangle signs on any four vertices
is positive, forcing every triangle uvx to be negative. Triangle signs
determine the switching class, proving the first case. If only deletion v
is balanced, gauge its complement positive; split that complement according
to the signs of edges at v. A negative cycle uses one neighbour from each
part, yielding rs(n-3)!/(n-k)! negative k-cycles and the second case.

In the remaining case b=0, induction applies to each deletion. The exact
identity is

```text
sum_v S_D(H-v) = sum_(k=3)^(D+1)(n-k)c_k.               (3)
```

## D5 induction

The exhaustive bases and equality multiplicities are

```text
n                  6       7       8
minimum S5        64     205     516
labelled ties     15      21      28.
```

Every tie is a single-edge class. Suppose n>=9 and b=0. If S5<=A_(n,5),
then induction, (3), (1), and c3<=binom(n,3) give

```text
n A_(n-1,5)
 <= (n-5)S5 + 2c3+c4-c6
 <= (n-5)A_(n,5) + 2binom(n,3).
```

This is impossible, since the difference between the leftmost and rightmost
quantities factors as

```text
(n-5)(n-4)(n-3)(3n-25)/3 > 0,                n>=9.     (4)
```

Consequently b=0 is strictly above the candidate. The other two cases prove
both the all-n minimum and its exact equality classification.

## D6 induction

The independently recomputed bases are

```text
n                  7       8
minimum S6       325    1236
labelled ties     21      28.
```

Again every tie is a single-edge class. Suppose n>=9 and b=0. If
S6<=A_(n,6), induction and (2) imply

```text
n A_(n-1,6)
 <= (n-6)S6 + 3c3+2c4+c5-c7
 <= (n-6)A_(n,6) + 3binom(n,3).
```

The difference is exactly

```text
(n-4)(n-3)(2n^3-42n^2+285n-620)/2.                     (5)
```

Putting n=9+t changes the cubic factor to 2t^3+12t^2+15t+1, positive for
all t>=0. Thus b=0 is again strict, completing the D6 theorem and equality.

## Independent exact audits and boundaries

```bash
python3 -B 04-computation/even_graph_d5_six_vertex_synthesis_sep05.py
python3 -B -O 04-computation/even_graph_d5_six_vertex_synthesis_sep05.py
g++ -O3 -std=c++17 -mpopcnt 04-computation/even_graph_d5_six_vertex_synthesis_sep05_independent.cpp -o /tmp/math-d5-sep05-audit
/tmp/math-d5-sep05-audit
```

The primary path generates simple cycles by subset/cyclic-order enumeration
and computes every character by an integer Walsh transform. The independent
C++ path generates cycles by full permutations, deduplicates them, and tests
every cycle parity for every character directly. It has no transform and
imports no Python profile or theorem table. Both paths exhaust all 2,130,944
classes at n=6,7,8, check both base equality sets, and prove L6/L7. All gates
use explicit failure branches; Python optimization cannot remove them.

The source object is a switching class, the intermediate observer is its
negative-cycle vector, and the target is the cumulative spectral character.
Induced restriction preserves edge signs and cycle parity but loses ambient
incidence; binomial occurrence factors restore exactly the data needed by
deletion. The scalar S_D alone cannot perform this step. This is a counting
connection to the repository's component-incidence theme, not a map from
tournaments or LRC configurations into signed complete graphs.

The cumulative conjecture for D>=7 remains OPEN. Individual even layers
still have the antibalanced zero hostile. No all-D theorem, tournament
H>=disc inequality, LRC(14) theorem, or Boolean metagraph gap follows.
