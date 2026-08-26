---
id: THM-4175
title: "Haar failure-atom deletion tomography and anchor exchange"
status: >
  PROVED RELATIVE TO THM-4150/4156/4170 + VERIFIED-EXACT COMPLETE C++ ALL-q
  AND INDEPENDENT PYTHON q=50 WALL-LATTICE/TOMOGRAPHY AUDITS; LRC(14) OPEN.
  Boolean deletion masses
  recover every labelled failure atom by Mobius inversion and have
  nonnegative mixed differences. Applied to the THM-4156 pool, this closes
  every one-anchor-exchange body for every newcomer q outside the pool by
  at most six total pool deletions, a bound sharp for this mechanism at
  q=50, and supplies 6,660,225 bodies per newcomer.
source: codex-frontier-synthesis-creative-20260826av
depends_on:
  - THM-4150-safe-set-haar-measure-universal-odd-tail-lrc14-transfer
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
  - THM-4170-triple-deletion-matching-eventual-haar-odd-tail-transfer
related:
  - THM-4114-ocf-mobius-positivity-tropical-defect-layer-and-opposite-ear-cut-curvature
  - THM-4160-anchored-haar-deletion-cover-and-content-tower
  - THM-4166-two-deletion-repair-graph-haar-odd-tail-transfer
  - THM-4172-tournament-multideletion-tomography-and-constant-finite-banks
  - THM-4174-six-deletion-completion-of-divisor-complete-newcomer-haar-transfer
script: 04-computation/lrc14_haar_failure_atom_anchor_exchange_thm4175.py
output: 05-knowledge/results/lrc14_haar_failure_atom_anchor_exchange_thm4175.out
independent_audit_script: 04-computation/lrc14_one_anchor_exchange_all_newcomers_thm4175_independent.cpp
independent_audit_output: 05-knowledge/results/lrc14_one_anchor_exchange_all_newcomers_thm4175_independent.out
script_sha256: d0dd60d900655147d399e2ba4aa4ba18bb37472ed43e60ca941e229d8669976b
output_sha256: 56b9bb8f0608320fedd34357ceec7b00cf1719b97d6e90c13c5304757907b63f
independent_audit_script_sha256: b1f40b89b9a4335b5ad261346c09e419893a1f8211b7c14c65db58073e74a91b
independent_audit_output_sha256: b5e048cd4b1d262e943ebbf56eb53cebd39aa7eee96328001e9940b21c5b081e
hash_basis: raw LF bytes
primary_audit: >
  PASS. A standard-library Fraction/submask implementation reconstructs the
  7,134-wall pool at q=50, retains all 2,950 full failure-mask buckets,
  checks 1,024 labelled deletion masses, their complete Boolean inversion,
  all 59,049 disjoint mixed differences, the symmetric transform and its
  inverse, both mixed and forced-anchor hostile ladders, and the literal
  body census. Normal, optimized, and hash-seeded streams byte-match.
independent_audit: >
  ACCEPT. A separate signed-128-bit C++ implementation reconstructs the
  integer wall lattice, proves three strict-tail nine-matchings, scans every
  finite residual and transversal through total deletion depth six, and
  literally enumerates the uniform divisor-complete bodies. Its 113,249,682
  exact comparisons have zero threshold equalities. O2, O0, and UBSan
  streams byte-match. The paths share the exact wall identity and threshold
  but use different languages, arithmetic representations, and reducers.
---

# THM-4175 -- failure atoms, deletion tomography, and anchor exchange

**PROVED RELATIVE TO THM-4150/4156/4170 + VERIFIED-EXACT COMPLETE C++ ALL-q
AND INDEPENDENT PYTHON q=50 AUDIT; LRC(14) REMAINS OPEN.**

## 1. A finite-measure deletion theorem

Let `(Omega,mu)` be a finite measure space, let `S` be a finite labelled set,
and let `G_s` be a measurable good event for every `s in S`. Put

```text
B_s=Omega\G_s,
C_F=(intersection_(s in F) B_s)
    intersection (intersection_(s in S\F) G_s),
w(F)=mu(C_F)                                               (1)
```

for `F subset S`. The `C_F` form the exact labelled failure-atom partition.
For a deletion set `R subset S`, define

```text
g(R)=mu(intersection_(s in S\R) G_s).                     (2)
```

> **Failure-atom deletion tomography.** For every `R,F subset S`,
>
> ```text
> g(R)=sum_(F subset R) w(F),                              (3)
> w(F)=sum_(R subset F)(-1)^(|F|-|R|) g(R).               (4)
> ```
>
> If `A intersection R=empty` and
>
> ```text
> Delta_A g(R)=sum_(H subset A)(-1)^(|A|-|H|)g(R union H),
> ```
>
> then
>
> ```text
> Delta_A g(R)=sum_(H subset R)w(A union H)>=0.            (5)
> ```
>
> Consequently the complete labelled deletion table determines every atom,
> every threshold repair hypergraph, and every mixed deletion gain.

Indeed, a point is admitted by `(2)` exactly when every label it fails lies
in `R`, which proves `(3)`. Boolean Mobius inversion gives `(4)`. Expanding
the mixed difference, the coefficient of `w(F)` is nonzero exactly when
`F\R=A`; this proves `(5)`.

There is also an exact symmetric quotient. Write `n=|S|` and

```text
W_s=sum_(|F|=s)w(F),              A_r=sum_(|R|=r)g(R).   (6)
```

Then

```text
A_r=sum_(s<=r) binom(n-s,r-s)W_s,                        (7)
sum_r A_r z^r=sum_s W_s z^s(1+z)^(n-s),                 (8)
W_s=sum_(r<=s)(-1)^(s-r)binom(n-r,s-r)A_r.              (9)
```

This quotient is intentionally lossy. For example, take `S={1,2,3}` and
`w(empty)=2/5`. Giving the three singleton atoms weights `(3/5,0,0)` or
`(1/5,1/5,1/5)` produces the same `W_s` and hence the same `A_r`. At threshold
`4/5`, however, the first one-deletion repair hypergraph has the single edge
`{1}`, while the second is empty. Thus symmetric layers do not retain labelled
overlap or transversal data.

The inheritance pass is now explicit. The closest proved algebraic mechanism
is THM-4114's Mobius positivity, while THM-4172 is the survival-side deletion
dual: there deletion must avoid a support, whereas here deletion must contain
a failure support. The canonical hostile is the symmetric example above. The
corrected near miss is therefore “the deletion-size histogram determines a
repair cover.” It does not. The least-used sidecar is the full labelled atom
support, retained below until the transversal decision has been made.

## 2. The one-anchor-exchange theorem

Retain the THM-4156 anchor set, pool, and optional ground set

```text
A={120,126,143},

P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290},

O=P\A.                                                   (10)
```

For a finite positive set `U`, put

```text
G_U={y in R/Z:min_(u in U)||uy||>=1/14}.                (11)
```

Fix an omitted anchor `a in A`, a positive newcomer `q notin P`, and a total
deletion depth `t`. A legal deletion is

```text
D in binom({a} union O,t).                              (12)
```

Project it to `rho_a(D)=D intersection O` and define the nonuniform repair
hypergraph on `O`

```text
F_(a,t)(q)={rho_a(D): D satisfies (12) and
  mu(G_((P union {q})\D))>=4/63}.                       (13)
```

An edge has size `t-1` if `D` deletes `a`, and size `t` if it retains `a`.
Retaining `a` is lawful: `a` is absent from the target body, so it is only a
harmless extra constraint in the safe superset.

For `K in binom(O,8)`, set

```text
H_(a,q,K)=(A\{a}) union {q} union K.                    (14)
```

> **One-anchor-exchange theorem.** For every choice in `(14)`, some total
> depth `t<=6` satisfies `tau(F_(a,t)(q))>8`. Consequently, for every positive
> integer `c` and every two distinct positive odd integers `r,s`, there is an
> `x in R/Z` such that
>
> ```text
> min_(v in 2cH_(a,q,K) union {r,s})||vx||>=1/14.       (15)
> ```
>
> The depth-six bound is sharp for the projected-repair mechanism: at `q=50`,
> all three depth-five hypergraphs have transversal number at most eight.

The exact criterion behind the statement is

```text
tau(F_(a,t)(q))>8
 iff every K in binom(O,8) misses some projected repair. (16)
```

The forward direction is the definition of nontransversality. Conversely, a
transversal of size at most eight can be extended to an eight-transversal.
If `rho_a(D)` misses `K`, then no member of `D` belongs to `(14)`, because the
only deleted label outside `O` can be the already omitted `a`. Hence

```text
H_(a,q,K) subset (P union {q})\D,
G_((P union {q})\D) subset G_(H_(a,q,K)).              (17)
```

Equations `(13),(17)` give Haar mass at least `4/63`; THM-4150 gives `(15)`
at content one. The surjective circle map `y -> cy` preserves Haar measure
of the safe-set pullback and gives every content. A repair disjoint from `K`
at depth below twenty can be extended inside `({a} union O)\K`; thus the
predicate in `(16)` is monotone upward in `t` through the depths used here.

## 3. Exact failure filtration

For `t=3,4,5,6`, let

```text
B_t^a={q positive:q notin P and tau(F_(a,t)(q))<=8}.   (18)
```

The complete depth-three failure sets are

```text
B_3^120={3,6,22,24,25,31,46,48,50,55,59,64,70,72,75,83,
 93,96,100,103,104,105,110,116,122,127,128,130,140,147,
 153,158,166,172,173,183,186,192,206,208,210,220,244,256,
 260,270,271,282,294,306,313,320,325,332,346,361,366,372,
 378,381,383,384,416,420,437,440,462,512,516,519,520,540,
 550,567,626,650,704,722,744,768,924,1134},

B_3^126={3,6,22,24,25,31,46,48,50,55,59,64,70,72,75,83,
 93,96,100,103,104,105,110,116,122,127,128,130,140,147,
 148,153,158,166,172,173,183,186,192,206,208,210,220,244,
 256,258,260,270,271,282,294,306,313,320,325,332,346,361,
 366,372,378,381,383,384,416,420,437,440,462,512,516,519,
 520,540,550,567,626,650,704,722,744,768,924,1134},

B_3^143={3,6,22,24,25,46,48,50,55,64,70,72,75,83,93,96,
 100,103,105,110,122,128,140,147,153,158,166,172,173,183,
 186,192,206,208,210,220,256,260,270,282,294,306,313,320,
 325,332,346,366,372,384,416,420,440,462,512,516,519,520,
 550,567,704,744,768,924}.                                  (19)
```

The remaining layers are

```text
B_4^120={6,24,25,50,96,100,105,128,140,183,192,210,256},
B_4^126={6,24,25,50,96,100,105,128,140,183,192,210,256,366},
B_4^143={6,25,50,96,100,105,128,192,210,256},

B_5^120=B_5^126=B_5^143={50},
B_6^120=B_6^126=B_6^143=empty.                         (20)
```

Thus the failure counts through depths `3,4,5,6` are respectively

```text
omit 120: 82,13,1,0;
omit 126: 84,14,1,0;
omit 143: 64,10,1,0.                                  (21)
```

At the common hostile `q=50`, depth four has four minimal repairs and
`tau=1`, with cover `{193}`. The exact depth-five and depth-six ledgers are

| omitted `a` | `t=5` raw/minimal | `tau` and one cover | `t=6` raw/minimal | `tau` |
|---:|---:|:---|---:|:---|
|120|`2,467/535`|`5`, `{95,168,193,240,290}`|`59,320/13,158`|`>8`|
|126|`2,126/1,417`|`5`, `{95,168,193,240,290}`|`51,810/11,051`|`>8`|
|143|`2,476/654`|`7`, `{80,85,88,145,168,193,286}`|`60,509/14,321`|`>8`|

Edge count is again not the invariant. Thousands of depth-five repairs still
admit small covers; depth six changes their labelled overlap enough to cross
the required transversal threshold.

## 4. Why the all-`q` audit is finite

At depth three, each omitted anchor has nine repairs with pairwise disjoint
projections. In the following display, tag `A` means delete the omitted
anchor together with the listed pair; tag `P` means retain it and delete the
listed triple.

```text
a=120:
 A:(42,88),(80,252),(85,170),(95,190),(145,290),(168,286),
   (176,193),(240,264); P:(8,16,132).                  cutoff 9,376

a=126:
 A:(42,88),(60,252),(145,290),(176,193);
 P:(8,16,286),(30,85,170),(63,132,264),(80,95,190),
   (84,168,240).                                      cutoff 14,511

a=143:
 A:(15,252),(16,286),(85,170),(95,190),(132,264),(145,290);
 P:(8,168,176),(42,63,88),(80,193,240).               cutoff 7,883
                                                               (22)
```

For a fixed deletion `D`, write `G_(P\D)` as `c_D` circle intervals and let
its measure be `m_D`. The exact safe-comb discrepancy used in THM-4170 gives

```text
mu(G_(P\D) intersection G_{q})
 >=(6/7)m_D-6c_D/(49q).                                (23)
```

Every row in `(22)` has `54m_D-4>0`. The exact integer ceiling

```text
ceil(54c_D/(7(54m_D-4)))                               (24)
```

is at most the displayed cutoff. Therefore all nine disjoint projections are
repair edges at and above that cutoff, forcing `tau>=9`. The frozen output
prints every exact `m_D`, component count, and row ceiling.

Below the three cutoffs there are respectively `9,345`, `14,480`, and `7,852`
newcomers after removing the 30 pool labels. The signed-128-bit audit builds
the `7,134` rational walls and `2,950` failure-mask groups, integrates the
`q`-safe comb exactly on every cell, and tests

```text
9N(q,D)>=8qL,             L=18,241,159,416,480.        (25)
```

It exhausts all depth-three rows and then only the exact residuals at depths
four through six. The total is

```text
113,249,682 exact threshold comparisons,
0 threshold equalities.                                (26)
```

Minimal projected edges are deduplicated and inclusion-reduced before an
exact depth-eight transversal recursion. A greedy disjoint-edge packing is
used only as a valid lower-bound prune; every returned cover is checked.
Equations `(22)--(26)` prove `(19)--(21)` and hence the theorem.

## 5. Family size and the divisor-complete seam

For fixed `q notin P`, the missing anchor identifies `a`, the unique
outside-pool label identifies `q`, and the remaining eight optional labels
identify `K`. Hence `(14)` gives exactly

```text
3 binom(27,8)=6,660,225                               (27)
```

distinct scale-one bodies per newcomer. They are disjoint from THM-4156 and
THM-4174, whose bodies contain all three anchors. For `N>=290`, the exact
newcomer count with `1<=q<=N` is `(N-30)*6,660,225`.

There is a large q-independent primitive divisor-complete subfamily. Count a
body when the retained anchors together with `K` already contain a multiple
of every modulus `2,...,14`, before using `q`. Literal enumeration and
inclusion-exclusion give

```text
omit 120:
 binom(27,8)-binom(18,8)-binom(17,8)-binom(20,8)
 +binom(11,8)+binom(14,8)+binom(12,8)-binom(8,8)
 =2,029,699;

omit 126:
 binom(27,8)-binom(23,8)-binom(25,8)+binom(22,8)
 =967,956;

omit 143:
 binom(26,7)=657,800 divisor-complete,
 binom(26,7)-binom(20,7)=580,280 uniformly primitive. (28)
```

The first two rows are primitive because the retained anchor pairs have gcd
one. In the third row, divisor completeness forces `286 in K`. Exactly
`binom(20,7)=77,520` choices have only even remaining optionals and gcd two;
all others contain an odd optional and have gcd one. Thus every even newcomer
has at least, and every odd newcomer at least,

```text
3,577,935 and 3,655,455                                (29)
```

primitive divisor-complete exchange bodies respectively. These are the exact
parity-class minima: take `q=2p` or `q=p` with a prime `p>290`, so `q` supplies
no missing modulus in `2,...,14`; odd `q` breaks the gcd-two rows and even `q`
does not. Equation `(28)` counts the named uniform subfamily, not all extra
q-specific divisor completions. Its primitive cores have distinct content
towers because `gcd(2cH)=2c`.

## 6. Concrete deletion tomography at the hostile

The independent Fraction audit specializes `q=50` but retains all thirty pool
labels in the failure masks. On the same `7,133` cells it records `2,950`
full mask buckets. For the ten-label slice

```text
(8,15,42,63,85,88,120,126,143,252),                   (30)
```

it checks all `1,024` direct deletion integrals against zeta sums, recovers
every exact atom by `(4)`, checks all `3^10=59,049` mixed differences in `(5)`,
and verifies both symmetric transforms `(7),(9)`. Forty selected-set atoms
have positive mass. The same cell integration independently reproduces the
entire mixed `q=50` table above and, as a hostile sidecar, the weaker convention
that always deletes the omitted anchor. The latter needs total depth seven
only for `a=126`; allowing both legal branches in `(12)` is what sharpens the
uniform theorem to six.

This is the operational payoff of `(3)`: once exact failure atoms are retained,
every deletion mass is a submask sum. But optimization must still use the full
labelled edge family; the symmetric quotient cannot see the crossing in
transversal dispersion.

## 7. Scope, loss ledger, and replay

The source-to-target contract is

```text
source:       labelled Haar failure atoms and total-depth pool deletions
target:       every one-anchor-exchange body with every odd tail pair
map:          failure zeta sum -> projected repair -> missed edge -> THM-4150
preserved:    exact 4/63 mass, body labels, omitted-anchor sidecar, content
destroyed:    chosen repair, component address, body-specific safe phases
sidecar:      full projected edge support until the transversal decision
positive:     every q=50 depth-six family has no eight-cover
hostile:      q=50 depth five has tau 5,5,7 despite thousands of repairs
decisive test: strict-tail nine-matchings plus 113,249,682 exact comparisons.
                                                               (31)
```

The restriction `q notin P` is type-critical for the novelty count. If `q=a`,
one recovers an old all-anchor body; if `q` is a retained anchor or belongs to
`K`, the set has only ten labels; and if `q in O\K`, it is an already safe
eleven-subset of `P` with non-newcomer presentations. The theorem remains a
fixed-pool, one-anchor-exchange, dyadic-body/odd-tail result. It gives no
physical entry theorem, no arbitrary-body classification, no mixed/even-tail
closure, no necessity of the Haar threshold, and no proof of LRC(14).

Replay the Fraction audit with

```bash
python3 -B 04-computation/lrc14_haar_failure_atom_anchor_exchange_thm4175.py
python3 -B -O 04-computation/lrc14_haar_failure_atom_anchor_exchange_thm4175.py
PYTHONHASHSEED=4175 python3 -B \
  04-computation/lrc14_haar_failure_atom_anchor_exchange_thm4175.py
```

Replay the complete all-newcomer C++ audit with

```bash
g++ -std=c++20 -O2 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_one_anchor_exchange_all_newcomers_thm4175_independent.cpp \
  -o /tmp/lrc14-thm4175
/tmp/lrc14-thm4175

g++ -std=c++20 -O0 -Wall -Wextra -Wpedantic -Werror \
  04-computation/lrc14_one_anchor_exchange_all_newcomers_thm4175_independent.cpp \
  -o /tmp/lrc14-thm4175-o0
/tmp/lrc14-thm4175-o0

g++ -std=c++20 -O1 -g -fsanitize=undefined \
  -fno-sanitize-recover=undefined \
  04-computation/lrc14_one_anchor_exchange_all_newcomers_thm4175_independent.cpp \
  -o /tmp/lrc14-thm4175-ubsan
/tmp/lrc14-thm4175-ubsan
```

All declared replay streams byte-match their frozen outputs. **QED.**
