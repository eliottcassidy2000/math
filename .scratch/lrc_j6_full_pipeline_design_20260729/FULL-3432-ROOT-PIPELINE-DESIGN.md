# Full 3432-root j6 proof-pipeline design

**Status:** scratch architectural and finite-exact audit, 2026-07-29.  No
uniform theorem is claimed here.  The all-root hitting atlas is **PROVED**
as `THM-2896-seven-body-adaptive-six-cover-hitting-gate-atlas.md`; the
singleton-complement parity descent and four-root closure is **PROVED** as
`THM-2895-singleton-complement-parity-descent-and-four-root-j6-closure.md`;
and the rank-selective cap activation used below is **PROVED** as
`THM-2897-partition-cap-tropical-convolution-and-alternating-pair-ladder.md`.
The child `q3+B2` composition is **PROVED** as
`THM-2900-flag-conditioned-rank-selective-partition-closure.md`.

## 1. Exact input size

The all-root top-forty discovery run gives:

```text
roots                                      3,432
least gate sizes K                         2 through 25 (no K=24)
sum K                                     41,415
mean K                                    12.0673076923
sum K^2                                  535,645
other-gate coverage records              494,230 = sum K(K-1)
outside-gate top-five records             207,075 = 5 sum K
total reusable apex candidate records     701,305
fixed K=12 roots                            1,976
roots needing K>12                         1,456
maximum least K                               25
unique K=25 root             (1,8,10,11,12,13,14)
maximum root tail-first                      2,201
maximum retained root speed                    574
```

The complete K distribution is

```text
(2:1, 3:3, 4:7, 5:33, 6:60, 7:140, 8:207, 9:285,
 10:375, 11:436, 12:429, 13:385, 14:333, 15:245,
 16:188, 17:112, 18:76, 19:55, 20:30, 21:20,
 22:9, 23:2, 25:1).
```

At base horizon 1600, constructing one complete singleton ledger for every
apex means

```text
41,415 * (1600-14) = 65,684,190
```

exact carrier/speed evaluations, plus only those tail extensions required
by the strict discrepancy seals.  The gate labels all lie below 575, so
their 494,230 cross-coverages should be extracted from the same scan rather
than recomputed.

This profile construction, and later child-residual top-three scans, are
the likely dominant work.  The Boolean activation search is secondary.

## 2. One reusable profile per root/apex

For a root `E`, let `A_E` be its sealed least top-K hitting set.  For every
`a in A_E`, construct the literal carrier

```text
C_a = G_E minus D_a = G_(E union {a}),       h_a=|C_a|.
```

Store:

1. its exact interval list, mass and component count;
2. `c_a(b)` for every `b in A_E minus {a}`;
3. the globally sealed five largest `c_a(w)` over `w notin A_E`;
4. the discrepancy cutoff and exact scalar controls.

The sorted union of (2) and (3), at most `K+4 <= 29` entries, evaluates

```text
S_a(P) = sum of the five largest allowed entries
```

for **every** marked prefix `P subset A_E`.  Only labels actually in `P`
may be removed.  In particular, a low current rank is never permission to
remove a future label.  Five outside-gate entries suffice because no sixth
outside-gate entry can enter the overall top five.

Process one root end-to-end and emit a hashed ledger.  Do not retain all
65 million exact fractions in memory.  Cache interval carriers and
coverage arrays by their literal sorted family, extending the array lazily
to the largest requested cutoff.

## 3. Scalar closure as an exact closed-state graph

Let

```text
cl(P) = least fixed point above P of
        P |-> P union {a : S_a(P)<h_a}.
```

This is an extensive, monotone, idempotent closure operator.  Search only
its closed masks:

```text
start = cl(empty)
edge  X --a--> cl(X union {a}),     X=cl(X), a notin X.
```

Every edge costs one nonscalar (paid) branch.  Breadth-first search gives
the exact minimum scalar-bootstrap seed number.

Proof of exactness: for any ordered seed list `(a_1,...,a_m)`, repeated
application of the displayed edge ends at `cl({a_1,...,a_m})`.  Conversely
the labels of every path are a seed list with that endpoint.  Hence shortest
path length to `A_E` equals minimum seed size.  Closed masks deduplicate all
seed subsets having the same scalar consequence.

On the committed four-root battery, early-stopping BFS reproduced minima
`(4,3,5,3)` with

```text
9,960 closed states
34,929 edges
10,292 closure calls
```

instead of the committed raw census's `34,661` seed subsets and `574,940`
activation-predicate checks.  The individual state/edge counts were

```text
K=19: 1,956 /  5,681
K=20:   213 /    402
K=21: 7,694 / 28,672
K=13:    97 /    174.
```

The unique K=25 root has exact scalar minimum `6`; its BFS used

```text
68,407 closed states, 323,776 edges, 68,407 closure calls.
```

This is tractable and is a better hostile scaling datum than extrapolating
the easy roots.

## 4. Proof order: the non-negotiable invariant

A scalar seed is a workload object, not a branch proof.  The theorem-bearing
state is an **actual processed prefix** `P` together with its explicit
order.

For a paid transition at closed state `P`:

1. prove the branch for `a` on the literal suffix excluding exactly
   `P union {a}`;
2. append `a`;
3. append a recorded fair scalar activation order from
   `P union {a}` to `cl(P union {a})`.

Simultaneous scalar rounds are soundly linearized: every member is certified
at the pre-round prefix, and earlier insertions only enlarge that prefix.

Thus the closed-state predecessor path is already the correct interleaved
proof order.  It is stronger and safer than putting an unordered minimum
seed into the prefix all at once.  No paid branch may use a later paid apex
as an exclusion.

The cheap BFS may first ignore whether its paid edges possess stronger cap
or parity certificates.  It supplies a scheduling seed, not the final
cost label.  Saturate every actual state by the monotone THM-2897 activation
below, then pay for parity only at the resulting joint fixed point.  A
failed candidate path is not a mathematical obstruction: search alternative
paths, then larger paths.

## 5. Rank-selective pair and matching activation before parity

For an apex `a` at actual prefix `P`, let

```text
q5_a(P) = fifth largest singleton coverage on C_a,
B2_a(P) = exact two-label union cap on C_a,
M22_a(P) = exact maximum weight of two vertex-disjoint pair edges,
```

over labels outside `P union {a}`.  THM-2897 proves the monotone activation

```text
q5_a(P) + 2 B2_a(P) < h_a.                                (RP)
```

Indeed, in any five-set choose its minimum-singleton member, bounded by the
**fifth** order statistic, and partition the other four labels into two
pairs.

The stronger matching repair

```text
q5_a(P) + M22_a(P) < h_a                                  (M)
```

retains the requirement that the two pair blocks use four distinct labels.
For a target `T=h_a-q5`, THM-2897 reduces the decision `M22<T` to the
finite `L`-heavy threshold graph, where `L=T-B2`.  This graph is intrinsic,
undirected, and weighted; ties are retained.  It is not a tournament.

Do not confuse the three adjacent quantities:

```text
scalar top-five: q1+q2+q3+q4+q5 < h_a;
matching repair:                   q5+M22 < h_a;
rank-pair fallback:                q5+2B2 < h_a;
H4 finiteness:                           q1 < 3h_a/7.
```

The singleton in `(RP)` is `q5`, not `q1`; using `q1` needlessly discards
distinctness.  Both `q5_a(P)` and `B2_a(P)` are nonincreasing under prefix
deletion, so `(RP)` belongs in the closure operator:

```text
repeat scalar top-five activations to a fixed point;
try the finite matching repair on eligible hostile states;
compute/cache exact B2 on remaining apex states;
append every passing rank-pair apex and repeat both layers;
buy an H4 parity edge only at the joint fixed point.
```

For an exact minimum-parity search, the closed masks must be closed under
scalar, matching, and rank-pair predicates.  The scalar-only BFS minima
remain valid scheduling upper bounds, not minimum parity counts.

## 6. Paid edge: singleton-complement parity certificate

At an apex carrier `C` with five labels remaining and allowed universe
excluding its actual prefix:

1. Compute the exact singleton maximum `B1=q1`.
2. Require `B1 < 3|C|/7`.
3. Put

   ```text
   tau=(|C|-B1)/4,       H4={w:c_C(w)>=tau}.
   ```

   The strict inequality makes `tau>|C|/7`, so discrepancy seals `H4`
   finitely.  Every hypothetical five-cover contains at least two H4
   labels.
4. For every unordered `L in C(H4,2)`, construct the literal nonempty
   residual `R=C minus D_L`.  It needs three more labels.
5. Apply, in this order:

   ```text
   inherited parent top-three bound;
   inherited parent B2(C)+B1(C) bound (one B2 amortized per paid edge);
   fresh globally sealed child top-three bound;
   fresh child q3(R)+B2(R) bound (THM-2900);
   fresh B2(R)+B1(R), only if child top-three fails.
   ```

6. A survivor is descended with `(k,s,ell)=(3,2,2)`: require its own
   `q1<5|R|/7`, enumerate the finite heavy H2 edges, subtract each
   literally, and use the longest-component singleton horizon plus exact
   tooth containment.

The parent inheritance is sound because, for
`R=C minus (D_x union D_y)`,

```text
c_R(w)<=c_C(w),        U_R(Q)<=U_C(Q).
```

Exact residual mass is still required.  The parent caps may range over the
larger pre-L allowed universe; that only weakens the bound.

### Exact inheritance audit

On the original 25 branches / 784 H4 pairs:

```text
parent top-three closes                       86
parent B2+B1 adds                             13
parent union                                 99  (12.6%)
fresh child profiles                        685
hard recursive rows                           5, all closed
parent paid pair unions                      206 across 23 parent caps.
```

On the exact 15-minimum-seed order / 1,464 H4 pairs:

```text
parent top-three closes                       28
parent B2+B1 adds                              9
parent union                                 37  (2.5%)
fresh child profiles                      1,427
hard recursive rows                           2, both closed
parent paid pair unions                      191 across 15 parent caps.
```

So inheritance is valid and cheap, but it does not remove the dominant
early-prefix child-profile workload.

### Stage child B2 only after top three

On the 15 paid branches, fresh child top-three closes `1461/1464` pairs.
Only **three** rows need a child B2 computation; its union closes one more,
and the remaining two close recursively.  A production implementation
must not compute B2 on all 1,464 rows as the exploratory helper does.

Likewise, replace the fixed child horizon 2500 by an adaptive top-three
seal.  After a provisional exact `q3>|R|/7`, scan only through

```text
ceil((99/70) r_R / (7(q3-|R|/7))) - 1
```

and verify the final strict tail inequality.  Start with a modest head,
extend, re-rank, and seal.  This directly attacks the dominant workload.

## 7. Hostile K=25 structural pilot

The exact minimum scalar-bootstrap seed path for

```text
E=(1,8,10,11,12,13,14)
```

uses apices

```text
23, 27, 19, 46, 18, 17
```

The clean proof schedule uses parity only at

```text
23, 27, 19, 46, 18,
```

after which scalar closure adds `(168,182)`.  At the resulting actual
prefix, apex `17` inherits the smaller fixed-rank certificate with excluded
set `(23,27,19,46,18,17)`.  Its THM-2897 `L`-heavy threshold graph has
endpoint cutoff `1844` and exactly three edges:

```text
{25,37}: 33906683/440740300
{25,45}: 3985831/53603550
{37,50}: 65302143/881480600.
```

The two edge-pairs reaching `T=h-q5` share vertex `25` or `37`.  The unique
disjoint pair in this threshold graph is `{25,45},{37,50}` and lies below
`T` by `69186919/39666627000`.  This is a threshold-graph gap, not a claim
that it is the exact global value of `h-(q5+M22)`.  The cutoff theorem
nevertheless proves `q5+M22<h`, so matching certifies apex `17`; scalar
closure then exhausts the gate.

The five paid apices satisfy `B1<3h/7`.  Their sealed H4 data on the
actual prefixes are

```text
a=23: Htail= 920, |H4|=13, pairs=78
a=27: Htail=1192, |H4|=13, pairs=78
a=19: Htail= 898, |H4|=10, pairs=45
a=46: Htail= 943, |H4|=11, pairs=55
a=18: Htail= 791, |H4|= 7, pairs=21
```

The full probe closes `277/277`: child top-three closes `276`, `B2+B1`
closes the last, and no recursive row survives.  No rank-pair activation is
used in this proof.  The repository-ready scratch verifier is

```text
python3 .scratch_thm2898_20260729/
  lrc14_j6_k25_five_parity_matching_closure_thm2898.py
```

The next shard after that is all `62` roots with `K>=20`, containing only
`1,289` apex profiles (3.11% of the full profile universe) while covering
every largest-gate geometry.

## 8. Concrete all-root run

Shard by `(stratum, K range)` and process each root as follows:

```text
ROOT(E):
  construct and globally seal least hitting gate A
  construct the K reusable apex profiles
  build the scalar top-five closure
  enrich it by finite q5+M22<h and then q5+2B2<h activations
  enumerate closed masks in BFS order
  for each candidate paid edge (P,a):
      try B1<3h/7 and estimate |H4| choose low-work edge first
      run full literal H4 -> child top3 -> rare B2 -> H2 -> singleton proof
      if PASS:
          record paid certificate
          append deterministic scalar cascade to cl(P union {a})
  stop at A; emit complete path/certificate ledger
```

For a proof-only run, stop at the first certified path.  For an exact
minimum-certified-path claim, exhaust every smaller BFS layer.

Every shard ledger must record universe, actual prefixes, all literal
families, strict margins, equality cases, tail cutoffs, H4 pairs, child
survivors, recursive heavy edges, singleton horizons, controls, and a
canonical digest.  Ordinary and optimized runs must agree byte-for-byte.

## 9. Precise obstructions to expect

The pipeline can fail for structurally different reasons; they must not be
merged:

1. `q1>=3h/7`: the first H4 finite-core certificate is unavailable.  A
   later prefix or different paid apex may repair it because q1 decreases
   under exclusions.
2. `q1<3h/7` but `tau` is close to `h/7`: H4 is finite but too large.
3. An H4 pair leaves an empty residual: this is a real small cover witness,
   not a successful recursive row.
4. Child top-three, q3+B2, and B2+B1 all fail: send only that row deeper.
5. Recursive `q1>=5m/7`: the H2 parity step is unavailable; use exact
   triple maximization or choose another paid edge/order.
6. A terminal singleton contains the literal residual: the branch remains
   open and the witness must be preserved.
7. A minimum scalar path fails parity: only that proof path failed.  Search
   another closed-state path before diagnosing the root.

The central underlying object is therefore not a tournament.  It is a
directed graph of scalar/matching/rank-pair-closed proof states, whose
expensive edges carry undirected heavy hypergraphs and literal residual
flags.  The direction comes from proof-prefix growth; the matching,
H4, and H2 relations themselves remain symmetric.

## 10. Theorem-ready proof outline for the K=25 root

1. THM-2896 gives the exact 25-label gate
   `(23,27,19,46,18,17,25,34,38,100,63,156,29,125,130,44,37,50,92,168,72,32,110,54,182)`
   and proves that every hypothetical six-cover meets it.
2. Order the gate by the locked combined schedule.  Assign a cover to its
   unique earliest gate apex.  The other five labels lie outside the actual
   marked prefix.
3. At every scalar step, the exact sum of the five largest allowed
   singleton coverages of the literal apex residual is strictly below its
   mass, so subadditivity excludes that branch.
4. The first five hostile apices `(23,27,19,46,18)` have `q1<3h/7`.
   THM-2895 forces two labels from each finite H4 core.  All `277` literal
   pair residuals fail to admit a three-cover: `276` by exact child
   top-three and the last by exact `B2+B1`.  No recursive row remains.
5. Scalar closure adds `(168,182)`.  At apex `17`, THM-2897's finite
   matching cutoff leaves the three-edge `L`-heavy threshold graph above.
   Its two hostile edge-pairs share a vertex, so no two disjoint edges
   reach `T=h-q5`; hence `q5+M22<h`.  Adding `17` triggers the locked scalar
   cascade that activates every remaining gate apex.
6. Every gate apex branch is excluded, contradicting the THM-2896 hitting
   conclusion.  Hence this one seven-body root has no six-label external
   cover.  No statement about the other 3431 roots or LRC(14) follows
   without their corresponding ledgers.
