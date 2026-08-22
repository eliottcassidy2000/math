# LRC(14) `z1=216`: intrinsic ruler cost prefixes and the status-cut taxonomy

**Status: FINITE-EXACT DIRECT PROGRESS, PROJECTED-ATLAS SCOPE ONLY.**  Two
bounded, hostile-audited probes close `23` rows in the maintained projected
`k=3`, `z1=216` necessary wall atlas.  They change the exact projected ledger

```text
373184 -> 373179 -> 373161,
z1=216 wall rows: 380 -> 375 -> 357,
projected cap: 216 (unchanged).                              (1)
```

This is not a physical-cover classification.  It does not address arbitrary
`k<=1`, the rung, endpoint origin, phase/current, or LRC(14).

Companions:

- [`lrc14_j7_k3_z216_first_nontrivial_ruler_cost_prefix_closure_scout_20260803.py`](../04-computation/lrc14_j7_k3_z216_first_nontrivial_ruler_cost_prefix_closure_scout_20260803.py)
  and its [frozen output](../05-knowledge/results/lrc14_j7_k3_z216_first_nontrivial_ruler_cost_prefix_closure_scout_20260803.out);
- [`lrc14_j7_k3_z216_second_nontrivial_ruler_cost_prefix_closure_scout_20260803.py`](../04-computation/lrc14_j7_k3_z216_second_nontrivial_ruler_cost_prefix_closure_scout_20260803.py)
  and its [frozen output](../05-knowledge/results/lrc14_j7_k3_z216_second_nontrivial_ruler_cost_prefix_closure_scout_20260803.out).

## Inheritance pass

**Closest proved engine.**  [THM-3139 -- projected-k3 z225 terminal and z224
screen double-layer descent](../01-canon/theorems/THM-3139-projected-k3-z225-terminal-and-z224-screen-double-layer-descent.md)
provides the exact ray/status upper-relaxation screen.  [THM-3281 --
projected-k3 z216 three natural wall-family screen descent](../01-canon/theorems/THM-3281-projected-k3-z216-three-natural-wall-family-screen-descent.md)
shows that complete `(gcd(216,L),L)` ruler fibres are lawful units of work.

**Canonical hostile.**  [THM-3264 -- projected-k3 z216 low-cost gcd8
seventeen-row terminal descent](../01-canon/theorems/THM-3264-projected-k3-z216-low-cost-gcd8-seventeen-row-terminal-descent.md)
and the independently audited [costly gcd-eight
closure](../05-knowledge/results/lrc14-z216-costly-gcd8-closure-opus-20260803.md)
show that the screen need not be empty: row `64` leaves eight masks, while
costly rows `238,370` leave `115,19`.  Thus a zero residual here cannot be an
unnoticed “every screen dies” convention.

**Corrected near miss.**  MISTAKE-331 and MISTAKE-333 forbid semantic hashes
that bind a solver-selected Farkas basis or its normalization-dependent
contradiction magnitude.  Both new companions verify every returned dual over
exact rationals but persist only the deterministic infeasible instance,
canonical `row[:19]`, branch counts, and solver-free cut data.

**Least-used relevant sidecar.**  The canonical atlas prints a component
count `r` which is not part of the four-field row returned by the screen.  The
pair

```text
(g,L,r),       g=gcd(216,L),       cost=L*r                 (2)
```

is useful as a work-order sidecar.  It is not a safety statistic: no row is
excluded because its cost is small, and no large-cost row is declared hard or
dangerous.

## Concept board

| live concept | representation / operation | preserved predicate | loss / cheapest test |
|---|---|---|---|
| projected wall row | labelled body `E`, ruler `L`, high gate | membership in the maintained necessary atlas | loses physical entry, phase, owners and current; run the complete upper screen |
| intrinsic ruler family | fibre of `(g,L)` | every row remains individually labelled and screenable | family equality supplies no common physical action; compare all rows, not one representative |
| work invoice | `sum_E L(E)r(E)` | deterministic scheduling order | preserves no safety implication; hostile is a costly row that still closes |
| status state | four binary status coordinates, marginals, capacity function and load tail | exact feasibility of the common status table | scalar verdict hides the cut; reconstruct the monotone event |
| physical restoration | endpoint/owner/phase/current sidecars | actual cover consequence | absent here; do not reverse the quotient |

The useful new comparison is between the last two rows: most solver statuses
are small monotone Boolean tail cuts.  The family lens selects which rows to
run; the cut lens explains why they die.  Neither restores physical LRC data.

## Exact universe and family queue

The pinned atlas has

```text
z1=216: 480 rows = 447 wall + 33 order,
gcd strata: 8^19, 18^1, 24^135, 36^15, 72^310.             (3)
```

After the audited gcd-eight, gcd-eighteen, order-row, and THM-3281 closures,
the live wall universe is exactly `380` rows in `39` complete `(g,L)`
families.  Its intrinsic invoice census is

| gcd `g` | families | rows | total `sum Lr` |
|---:|---:|---:|---:|
| 24 | 12 | 91 | 319,446,960 |
| 36 | 2 | 11 | 156,008,160 |
| 72 | 25 | 278 | 1,609,301,232 |
| **total** | **39** | **380** | **2,084,756,352** |

The first companion freezes all `39` ranks.  The beginning is

| rank | `(g,L)` | rows | `sum Lr` |
|---:|---|---:|---:|
| 1 | `(72,11088)` | 1 | 243,936 |
| 2 | `(24,18480)` | 1 | 443,520 |
| 3 | `(24,5880)` | 3 | 517,440 |
| 4 | `(72,27720)` | 1 | 609,840 |
| 5 | `(72,32760)` | 1 | 655,200 |
| 6 | `(72,3528)` | 16 | 1,361,808 |

The bounded policy is structural rather than a raw row-cost slice: take the
cost-ordered prefix of **whole ruler families** through the next family having
more than one row.  Singleton fibres before that stopping family are included
because they precede it in the same intrinsic order.

## First prefix: two sentinels and a triangle family

Ranks `1--3` contain five rows and cost `1,204,896`.  The three `L=5880`
bodies have a literal triangle description:

```text
common core = {4,6,10,14},
outer triple = {1,2,12},
bodies = core union each two-subset of the outer triple.     (4)
```

Their exact screen, together with the two lower-cost singleton fibres, is

| family | rows | states | crude | exact status | residual |
|---|---:|---:|---:|---:|---:|
| `gcd72/L11088` | 1 | 8 | 8 | 0 | 0 |
| `gcd24/L18480` | 1 | 22 | 21 | 1 | 0 |
| `gcd24/L5880` | 3 | 12 | 9 | 3 | 0 |
| **total** | **5** | **42** | **38** | **4** | **0** |

### The four exact statuses are elementary union cuts

Let `x_P` be the mass assigned to binary status pattern
`P subseteq {0,1,2,3}`.  At a load threshold `t`, suppose the
capacity-good patterns are exactly those meeting a coordinate set `S`.
Then their mass is at most the union bound

```text
sum_(j in S) marginal_j.                                  (5)
```

The load histogram simultaneously requires a tail mass `T`.  Each selected
status has `T` strictly larger than `(5)`:

| atlas row | threshold `t` | coordinates `S` | required `T` | marginal bound | gap |
|---:|---:|---|---:|---:|---:|
| 50 | 3 | `{1,3}` | 896 | 528 | 368 |
| 18 | 2 | `{0,1,3}` | 486 | 252 | 234 |
| 18, alternate | 3 | `{0,1,3}` | 272 | 252 | 20 |
| 140 | 2 | `{0,1,2}` | 608 | 420 | 188 |
| 299 | 4 | `{0,1,2,3}` | 80 | 60 | 20 |

Thus the first prefix does not merely replay four optimizer verdicts.  Its
only non-crude states have explicit Hall/union-bound contradictions.  The
triangle is a coherent family at selection level; its four divisor tuples and
coordinate subsets differ, so no equivariant common certificate is claimed.

## Second prefix: the complete `L=3528` lcm fibre

After those five rows, the cost queue starts with two singleton fibres and
then the sixteen-row family

```text
gcd(216,L)=72,       L=3528=14*252,       lcm(E)=252.       (6)
```

The second companion takes exactly this next prefix.  Its results are

| family | rows | states | crude | exact status | residual |
|---|---:|---:|---:|---:|---:|
| `gcd72/L27720` | 1 | 29 | 28 | 1 | 0 |
| `gcd72/L32760` | 1 | 95 | 79 | 16 | 0 |
| `gcd72/L3528` | 16 | 626 | 372 | 254 | 0 |
| **total** | **18** | **750** | **479** | **271** | **0** |

All `271` exact Farkas certificates are independently rebuilt and checked.
The solver-free status packet gives a sharper anatomy:

```text
248  literal coordinate-union cuts,
  1  coordinate-union cut after deleting patterns forced off by a zero marginal,
 11  two-fan cuts,
 11  remaining weighted exact-Farkas states.                (7)
```

For a two-fan, the capacity-good event is

```text
B union (A intersect (C union D)).                          (8)
```

Its mass is at most

```text
m_B + min(m_A, m_C+m_D).                                   (9)
```

In all eleven cases the histogram tail exceeds `(9)`.  The zero-marginal
case first discards every pattern containing the zero coordinate, after which
the ordinary union bound applies.  Consequently `260/271` status states have
short monotone-event proofs.  The remaining eleven are exactly verified by
the inherited multirow Farkas checker, but “weighted core” means only that the
two implemented elementary templates do not close them; it is not a proof of
minimal proof complexity.

The weighted-core row distribution is

```text
134:1, 121:1, 17:1, 27:2, 56:1, 138:3, 298:2.             (10)
```

This is the most reusable side result: the common-status LP is largely a
four-bit monotone-tail cut algebra.  A finite catalogue of monotone Boolean
capacity events may replace most generic dual searches and isolate the truly
coupled thresholds.

## Direction, controls, and loss ledger

The logical direction is only

```text
exact necessary upper relaxation empty
    => corresponding projected wall row empty.              (11)
```

No converse from screen feasibility is used.  The source-to-target contract
is

```text
source: labelled projected wall row E
map:    E -> (g,L,r) family -> exact ray states -> status tail event
kept:   every inherited projected filter and the one-way exclusion predicate
lost:   physical speed entry, endpoint origin, owner, phase, current, chronology
needed sidecar for LRC: a lawful physical lift retaining those lost coordinates
decisive local test: complete exact screen with a known residual control.
```

The hostile suite includes:

- prior `gcd24/L2352` row `14`, which exercises `35` exact statuses and closes;
- prior gcd-eight row `64`, which retains eight residual masks of quotient
  type `8`, proving the engine does not force empty output;
- the costly independent audit's rejected negative Farkas multiplier;
- exact checking of all `4+271` selected status certificates; and
- canonical digests excluding all solver-selected dual vectors and
  normalization magnitudes under MISTAKE-331/333.

The code contains no `assert`-dependent truth gates and no floating-point
mathematics.  Normal and optimized runs match their stored transcripts with
empty standard error.

## Honest stopping boundary and next tests

After `(1)`, exactly

```text
357 z1=216 wall rows in 33 complete ruler families          (12)
```

remain.  The next intrinsic prefix is the singleton
`gcd72/L72072` followed by the three-row `gcd24/L25872` family, with invoices
`1,729,728` and `2,121,504`.  Ranked next moves are:

1. run that four-row prefix as the next cheap direct ledger test;
2. classify the eleven weighted states in `(10)` by two-threshold laminar or
   submodular cuts before expanding the family queue;
3. compile the complete four-bit monotone capacity-event catalogue, with
   trigger, marginal bound, and counterexample for each event;
4. treat this as a projected-atlas compiler only until an endpoint/owner/phase
   lift proves that the quotient operation is physical.

The method card used was “audit and close sections under their next native
operation”: a proposed intrinsic family was closed under the full exact screen
before trusting its scalar cost.  The niche that emerged is the monotone
status-cut catalogue.  The wildcard was whether the Farkas rows were merely
coordinate unions; the hostile census repaired that guess to the four-part
taxonomy `(7)`.

## Reproduction and hashes

```text
python3 04-computation/lrc14_j7_k3_z216_first_nontrivial_ruler_cost_prefix_closure_scout_20260803.py --processes 3
python3 -O 04-computation/lrc14_j7_k3_z216_first_nontrivial_ruler_cost_prefix_closure_scout_20260803.py --processes 2

python3 04-computation/lrc14_j7_k3_z216_second_nontrivial_ruler_cost_prefix_closure_scout_20260803.py --processes 3
python3 -O 04-computation/lrc14_j7_k3_z216_second_nontrivial_ruler_cost_prefix_closure_scout_20260803.py --processes 2
```

```text
first script   e3718e71d2f46ea691e269aaf7e05c292b930c4d6a9b6d6d89a835bc38df75fe
first output   6dba76e65300371399376aede7e79b9a6fa865b6ebb58475a14ef8e262328be1
second script  bd2236e4e3fe6c47984e4ecf10cd0b87f2d4b937a8d5c0ebb6da70e81884dd5e
second output  a5d68cb56dfeb3f44e49133451c3b2445ca387370a0a21998c0a6be5531ed8c9
```
