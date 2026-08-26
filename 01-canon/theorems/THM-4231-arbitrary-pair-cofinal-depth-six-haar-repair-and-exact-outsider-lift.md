---
id: THM-4231
title: "Arbitrary-pair cofinal depth-six Haar repair and exact outsider lift"
status: >
  PROVED RELATIVE TO THM-4150/4156/4170/4191 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. For the displayed thirty-label pool, every nine-body
  becomes 1/14-Haar-safe after adjoining any two distinct integers q,r at
  least 17548. Together with the inherited zero- and one-outsider layers,
  this fills all C(32,11)=129,024,480 eleven-bodies in every such two-outsider
  chart. The adjacent 17547/17548 transition is exact only for the complete
  depth-six component-discrepancy activation filter; it is not a minimal
  literal entry threshold. Combined with THM-4227/4228, every remaining pair
  has gcd at most 3466 and lies on one of finitely many small-outsider rays or
  in one finite box. Arbitrary finite pair entry and LRC(14) remain OPEN.
source: codex-frontier-synthesis-20260826
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
  - THM-4191-complete-full-pool-newcomer-haar-transfer
  - THM-4227-two-outsider-scale-separated-depth-eight-haar-wedge
  - THM-4228-common-gcd-two-outsider-periodic-observable-haar-ray
related:
  - THM-4207-two-newcomer-sharp-depth-transition-base-surplus-composition-and-variable-pool-chart-number
  - THM-4211-fixed-fifty-cofinal-two-newcomer-haar-tail-and-eighteen-label-chart
  - THM-4214-two-newcomer-pascal-complete-eleven-body-haar-charts
  - THM-4223-cyclic-cut-cover-boolean-mobius-hierarchy-and-two-bad-owner-obstruction
primary_script: 04-computation/lrc14_arbitrary_pair_cofinal_depth6_primary_thm4231.py
primary_output: 05-knowledge/results/lrc14_arbitrary_pair_cofinal_depth6_primary_thm4231.out
independent_audit_script: 04-computation/lrc14_arbitrary_pair_cofinal_depth6_independent_audit_thm4231.py
independent_audit_output: 05-knowledge/results/lrc14_arbitrary_pair_cofinal_depth6_independent_audit_thm4231.out
primary_script_sha256: 794f0df69956e46c5c73ad6489498b1bd404b9d3643722b129f4b63b092c890a
primary_output_sha256: ddaa2fdd3c822126a5c51c48b526450f34cd22537a40421bba42524ed5c51834
independent_audit_script_sha256: a52740967b84dcef6be68e9ad362cdc51ea9552234c6f5cd8e725a87f0d947f9
independent_audit_output_sha256: 0e7dffe54514f3d46d8076b7dda0c69418ec774c6d47ba95aa2b50bb73f76cec
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. The primary midpoint-cell and lower-zeta implementation
  exhausts all C(30,6) deletions and proves the new cutoff by a path-complete
  cover DFS. A separate endpoint-toggle implementation reverse-scatters mass
  and signed component boundaries, then literally scans all C(30,9) bodies
  against independently ordered decks at both adjacent cutoffs. Normal and
  optimized Python replays byte-match both frozen outputs.
---

# THM-4231 -- arbitrary-pair cofinal depth-six Haar repair and exact outsider lift

**PROVED RELATIVE TO THM-4150/4156/4170/4191 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED; LRC(14) REMAINS OPEN.**

## 1. Statement and inheritance

Retain the THM-4156 pool

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                        (1)
```

For a finite positive label set `S`, put

```text
G_S={y in R/Z:min_(s in S)||sy||>=1/14},
alpha=4/63.                                              (2)
```

> **Cofinal arbitrary-pair theorem.** For every two distinct integers
> `q,r>=17548` and every `B in binom(P,9)`,
>
> ```text
> mu(G_(B union {q,r}))>=alpha.                          (3)
> ```

The quantifier is uniform over the entire northeast quadrant, apart from the
diagonal `q=r`; neither newcomer is fixed. This strengthens THM-4211's
fixed-`50` cofinal ray in a different direction. It does not include pairs
with one small outsider, arbitrary finite pair entry, or all thirteen-speed
instances.

The closest proved mechanism is THM-4170's component-discrepancy estimate.
The canonical hostile is THM-4207's failure of marginal one-newcomer deck
intersection. The corrected near miss is to regard a cover of a sufficient
repair deck as an unsafe body. The least-used sidecar is the exact component
count of each fixed-pool safe set, retained together with its mass.

## 2. Symmetric two-newcomer discrepancy

For `R in binom(P,6)`, define

```text
U_R=G_(P\R),        M_R=mu(U_R),        c_R=#components(U_R). (4)
```

THM-4170 equation `(9)` applies to every union of `c_R` circle intervals:

```text
mu(U_R intersect G_s)>=(6/7)M_R-6c_R/(49s).             (5)
```

Bonferroni inside `U_R` therefore gives

```text
mu(U_R intersect G_q intersect G_r)
 >=(5/7)M_R-(6c_R/49)(1/q+1/r).                         (6)
```

Set

```text
s_R=45M_R-4.                                            (7)
```

When `s_R>0`, define the integer activation

```text
kappa_2(R)=ceil(108c_R/(7s_R)).                         (8)
```

If `q,r>=kappa_2(R)`, then

```text
(6c_R/49)(1/q+1/r)<=12c_R/(49kappa_2(R))<=s_R/63.
```

Since `(5/7)M_R-alpha=s_R/63`, equation `(6)` proves

```text
mu(G_((P\R) union {q,r}))>=alpha.                       (9)
```

Thus every strict-limit deletion becomes a lawful joint repair on an explicit
uniform quadrant. No independence or equidistribution of `G_q` and `G_r` is
assumed; `(6)` uses only the two one-frequency discrepancy bounds and
Bonferroni.

## 3. Activation deck and transversal consequence

For an integer cutoff `Q`, let

```text
A_6^(2)(Q)={R in binom(P,6):s_R>0 and kappa_2(R)<=Q}.    (10)
```

If `q,r>=Q`, every edge of `(10)` is a lawful repair by `(9)`. If the
transversal number satisfies

```text
tau(A_6^(2)(Q))>9,                                      (11)
```

then every `B in binom(P,9)` misses some edge `R`. For that edge,

```text
B union {q,r} subset (P\R) union {q,r},
G_((P\R) union {q,r}) subset G_(B union {q,r}).          (12)
```

Safe-set monotonicity and `(9)` prove `(3)`. This is a sufficient
hypergraph certificate, not an equivalence between covers and unsafe bodies.

## 4. Complete primary census and adjacent transition

The exact pool arrangement has common denominator

```text
D=18,241,159,416,480,
7,134 walls, 7,133 open cells.                           (13)
```

The primary program exhausts all

```text
binom(30,6)=593,775                                     (14)
```

deletions. Exactly `140,082` have `s_R>0`, `453,693` have nonpositive
surplus, and there are no equalities. The minimum and maximum activations are

```text
3,077 at {88,143,168,193,252,286},
4,636,948,909 at {80,95,120,143,170,193}.               (15)
```

At the adjacent cutoffs,

```text
|A_6^(2)(17547)|=54,563,
|A_6^(2)(17548)|=54,566.                                (16)
```

The nine-set

```text
W={85,88,143,168,193,240,252,264,290}                  (17)
```

covers the complete `17547` deck. Three edges activate at `17548`; one edge
disjoint from `W` is

```text
R#={8,16,42,95,132,145},
M_R# ticks=1,694,858,026,164,
c_R#=206,
(45M_R#-4) ticks=3,303,973,511,460.                     (18)
```

A path-complete recursion on the complete `54,566`-edge deck chooses the
first uncovered edge and branches on each of its six vertices. It visits
`663,464` states, all terminally dead, and proves that no cover of size at
most nine exists. Therefore

```text
tau(A_6^(2)(17547))<=9,
tau(A_6^(2)(17548))>9.                                  (19)
```

Equations `(11)--(12)` now prove the theorem.

The exact labelled row fingerprints are

```text
all 140,082 strict rows:       a8b79ad77ad91a62,
Q=17547 filtered rows:         476fef92619d2c0b,
Q=17548 filtered rows:         d20636ace1522a29.         (20)
```

Each row contributes five little-endian unsigned 64-bit words
`(mask,kappa,mass_ticks,components,surplus_ticks)` to bytewise FNV-1a-64.

## 5. Independent exact audit

The referee implementation does not import the primary geometry, zeta
lookup, edge order, or cover search. It instead:

1. constructs all `7,133` cells by endpoint enter/leave toggles rather than
   midpoint classification;
2. forms mass and signed cyclic component-boundary coefficients and
   reverse-incidence scatters them to all `593,775` deletions rather than
   summing submasks separately for each deletion;
3. orders the filtered edges by an independent SplitMix key; and
4. scans every one of the `14,307,150` labelled nine-bodies at both cutoffs
   rather than solving the dual cover problem recursively.

At `17547`, the scan performs `233,058,925` incidence checks and finds
exactly one cover, namely `W`. At `17548`, it performs `233,056,301` checks
and finds zero covers; the closest body `W` misses `R#`. The referee also
checks every exact ceiling inequality

```text
7s_R(kappa_2(R)-1)<108c_R<=7s_R kappa_2(R)              (21)
```

after clearing `D`, and recovers both filtered decks directly with zero
boundary equalities. It agrees on every count and all three fingerprints in
`(20)`.

As hostile controls, an independent literal joint-wall sweep at the pair
`(17547,17548)` finds that both the newly activating repair and the old
covering body are already strictly safe:

```text
mu(G_((P\R#) union {17547,17548}))
 =1344542319402481/19682934366703608,
63mu-4=94832200881617/312427529630216>0;                (22)

mu(G_(W union {17547,17548}))
 =2384057538965009671/12317373137900310480,
63mu-4=1602002101638005831/195513859331750960>0.        (23)
```

Thus `17548` is the first cutoff at which this complete activation deck has
no nine-cover. It is not the first literal safe pair, and `(17)` is not an
unsafe-body witness.

## 6. Exact outsider lift and odd-tail consequence

Fix distinct `q,r>=17548`. Every eleven-subset of the thirty-two-label set

```text
P union {q,r}                                            (24)
```

is safe:

- with no outsider, THM-4156 proves the stronger full-pool bound and safe-set
  monotonicity handles every subset;
- with one labelled outsider, THM-4191 handles every ten-subset of `P`; and
- with both outsiders, `(3)` handles every nine-subset of `P`.

The labelled block count is exactly

```text
binom(30,11)+2binom(30,10)+binom(30,9)
 =54,627,300+2(30,045,015)+14,307,150
 =129,024,480=binom(32,11).                              (25)
```

For each two-outsider body in `(3)`, every positive integer `c` and every two
distinct positive odd integers `a,b`, THM-4150 gives some `x in R/Z` with

```text
min_(v in 2c(B union {q,r}) union {a,b})||vx||>=1/14.   (26)
```

The same transfer applies to the zero- and one-outsider bodies. These are
genuine infinite LRC(14) families of thirteen relative speeds, not a proof
for an arbitrary instance.

## 7. Boolean sidecar and transfer from the tournament work

For each cyclic pool cell `i`, retain its failed-label mask `F_i` and length
`ell_i`. For a deletion `R`, the exact formulas behind the computation are

```text
M_R=sum_i ell_i 1_(F_i subset R),
c_R=sum_i [1_(F_i subset R)
           -1_(F_(i-1) union F_i subset R)].            (27)
```

The first line is a lower Boolean zeta transform. The second retains the
cyclic-adjacency sidecar needed to turn safe-cell counts into component
counts. This is the lawful lesson from THM-4223's owner-refined Boolean
hierarchy: retain the tensor and the operation-specific adjacency, not only a
scalar total. There is no intrinsic pairwise orientation here, so no
tournament is imposed.

The connection contract is

```text
source:       cyclic pool cells (F_i,ell_i) and adjacent mask unions
target:       the depth-six joint activation deck
map:          lower-zeta mass and transition count -> (M_R,c_R,kappa_2(R))
preserved:    labels, exact mass, components, and the uniform q/r quantifier
destroyed:    component addresses and literal q/r phase alignment
sidecar:      adjacent failure-mask unions and literal boundary scans
decisive test: tau(A_6^(2)(Q))>9 on the complete deletion universe. (28)
```

## 8. Boundaries and next frontiers

This theorem proves neither arbitrary finite pair entry nor a minimal literal
cofinal threshold. It does, however, combine sharply with THM-4227. For a
pair of distinct positive outsiders, write

```text
m=min(q,r),        M=max(q,r).                           (29)
```

THM-4231 handles `m>=17548`. If `3391<=m<=17547`, THM-4227 handles the pair
as soon as

```text
M>=ceil(321902813232(m+130)/10633545731).               (30)
```

The right side is increasing and equals `535125` at `m=17547`. Consequently,
every pair not covered by THM-4227 or THM-4231 lies in

```text
m<=3390,
or
3391<=m<=17547 and M<=535124.                           (31)
```

THM-4228 separately handles `gcd(q,r)>=3467`. Thus every pair left by all
three theorems also satisfies

```text
gcd(q,r)<=3466.                                          (32)
```

The unproved pair plane is therefore a gcd-thinned finite box together with
finitely many fixed-small-outsider rays. This is a reduction, not a closure:
none of the three theorems supplies the cofinal tail on each ray with
`m<=3390`.

The main remaining LRC routes are therefore:

1. prove one cofinal tail for every outsider below `3391`, retaining its
   literal label, then exhaust the finite remainder in `(31)`;
2. combine phase-sensitive or common-gcd discrepancy with the symmetric
   Bonferroni deck without losing the literal pair labels;
3. count the intrinsic multiplicity of lawful repairs per residual body,
   rather than an order-dependent greedy separator count; and
4. leave the fixed pool through an entry or replacement theorem.

The fixed-object firewall remains essential: a compatible certificate for
moving pairs at every depth would not imply one fixed pair. Here `(6)--(12)`
keep the same literal `q,r,B` throughout and are uniform after one finite
cutoff.

## 9. Reproduction

From the repository root, each command must produce the matching frozen
output byte for byte:

```bash
python3 -B 04-computation/lrc14_arbitrary_pair_cofinal_depth6_primary_thm4231.py
python3 -O -B 04-computation/lrc14_arbitrary_pair_cofinal_depth6_primary_thm4231.py

python3 -B 04-computation/lrc14_arbitrary_pair_cofinal_depth6_independent_audit_thm4231.py
python3 -O -B 04-computation/lrc14_arbitrary_pair_cofinal_depth6_independent_audit_thm4231.py
```

The first path uses exact midpoint cells, submask zeta sums, and a
path-complete transversal search. The second uses endpoint toggles,
reverse-incidence scatter, and exhaustive body enumeration. Both run with
integer/rational arithmetic only.

**QED.**
