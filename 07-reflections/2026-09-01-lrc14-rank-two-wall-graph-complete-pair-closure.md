# LRC(14): the rank-two wall graph closes the fixed-pool pair chart

**Status:** PROVED RELATIVE TO THM-4231 + FINITE-EXACT + CLEAN-ROOM
INDEPENDENT AUDIT. **LRC(14) remains open.**

This session began as another endpoint descent with the inherited 3,925-mask
carrier. The decisive move was to stop asking whether that carrier could be
repaired at the next endpoint and ask what part of the literal wall geometry
already sufficed for *every* nine-body. Retaining only pair-safe cells whose
pool-failure set has cardinality at most two turns the body problem into one
weighted graph on the thirty pool labels. That graph closes not merely the
next endpoint but all `181,194` pairs in THM-4231's exact finite remainder.

The resulting theorem is strong but sharply typed: it completes the
arbitrary two-outsider chart over one fixed pool and, through THM-4150, gives
its universal doubled-body/two-odd-tail families. It does not map a general
thirteen-speed row into that chart.

## 1. Inheritance pass and concept board

The closest proved mechanism was THM-4231: it already closed every outsider
pair outside one exact finite remainder. The canonical hostile was the
`q=50` fibre, which repeatedly survived carrier descent and contains the
eventual normalized rank-two minimum `(50,70)`. The corrected near miss was
the proposed THM-4325 endpoint staircase: its computations were positive,
but endpoint-by-endpoint carrier bookkeeping was no longer the right proof
object. The least-used sidecar was the *cardinality of the pool-failure mask
on each literal wall cell*.

The live board was:

| Lane | Object | Decisive question |
|---|---|---|
| Anchor | THM-4231's `181,194` outsider pairs | Can one certificate cover every labelled nine-body on every residual pair? |
| Niche | rank-zero/singleton/pair wall cells | Does their mass alone exceed `4/63`? |
| Wildcard | normalized margin ordinals and tournament analogies | Which scale-free order is intrinsic, and is the relation actually directed? |

The niche overtook the anchor because it supplied a complete proof object.
The wildcard contributed two guardrails. First, THM-4286's normalized-margin
ordinal and MISTAKE-532 forced comparison by `ticks/D`, not by raw tick
numerators. Second, the pair observable is symmetric, so the correct object
is an undirected weighted graph. Orienting it as a tournament would destroy
the exact overlap correction.

## 2. The graph compiler

Fix outsiders `q,r` and partition the circle at all walls of the pool and the
two outsiders. On a pair-safe open cell `C`, let `F(C)` be the set of pool
labels that fail there. Aggregate cell widths for `F=empty`, `F={i}`, and
`F={i,j}` as `w0`, `wi`, and `wij`. For a nine-body `B`, the retained mass is

```text
L2(B)=w0+sum_(i notin B)wi+sum_(i<j, i,j notin B)wij.
```

All discarded cells have nonnegative width, hence `L2(B)/D` is a lower bound
for the full Haar-safe mass. If

```text
T=w0+sum_i wi+sum_(i<j)wij,
di=wi+sum_(j!=i)wij,
```

then

```text
L2(B)=T-sum_(i in B)di+sum_(i<j, i,j in B)wij.
```

This is weighted maximum nine-vertex coverage in disguise. The last term is
the exact repair for edges counted twice by the degree sum. Since it is
nonnegative,

```text
L2(B) >= T - (sum of the nine largest weighted degrees).
```

That one inequality certifies all `C(30,9)=14,307,150` bodies on a pair. A
nonpositive degree relaxation says only that the overlap correction matters;
it is not evidence of an unsafe body.

The more general compiler is already visible. If widths `w_F` are retained
for `|F|<=d`, then

```text
Ld(B)=sum_(F intersect B=empty, |F|<=d) w_F
```

is a weighted rank-`d` hypergraph complement mass. Selecting `B` maximizes
the weight of hyperedges it meets. The sum of incident degrees always
overcounts covered hyperedges and therefore gives a valid coarse lower
bound; exact optimization restores their intersection multiplicities. Rank
two is special because the complete correction is just the internal-edge
sum.

## 3. Exact closure and the hostile boundary

The maintained THM-4231 postprocessor exports the exact remainder with

```text
count                 181,194
ordered-pair FNV      3874fecac4ecbd8a
SHA-256               9dfbf0a8948bf23016ae40f40f9118020c9429ac60421681a9286fc4d34041a1
maximum endpoint      769.
```

The wide primary degree screen splits it as

```text
strict coarse-positive pairs    181,087
coarse-nonpositive pairs            107
highest coarse exception        (50,554).
```

The degree bound proves all bodies on the first `181,087` pairs. Complete
enumeration of `107*C(30,9)=1,530,865,050` cases proves every remaining
rank-two mass strict. A separately written branch-and-bound search agrees on
every minimum and least minimizing mask. A clean-room event-state sweep
independently reconstructs all `181,194` graphs and all `107` exact minima.

The scale-free minimum requires one more hostile test. Among the
coarse-positive rows, exactly three degree bounds could lie below the best
exact exceptional ratio: `(50,212)`, `(50,274)`, and `(100,110)`. Exact
optimization rules all three out. Therefore the global retained-mass minimum
on the finite remainder is

```text
(q,r)=(50,70),
B={80,85,88,95,143,145,168,193,240},
body mask=031c7400,
D=91,205,797,082,400,
L2=5,794,739,949,188,
63L2-4D=245,428,469,244,

L2/D
 =1448684987297/22801449270600
 =4/63+973922497/22801449270600.
```

This is an exact minimum of the **retained rank-two lower bound**. The full
safe mass is at least this value and can be larger because cells of failure
rank at least three were deliberately discarded. Calling the displayed
fraction the exact full Haar surplus would repeat a quotient-loss error.

Exactly one residual grid exceeds signed 64-bit range:

```text
(713,719), D=9,351,275,651,380,222,560.
```

The proof path therefore uses signed 128-bit arithmetic throughout. The old
signed-64 scout is a discovery artifact, not evidence for the theorem.

## 4. What was preserved and what was destroyed

```text
source:       pair-safe wall cells with labelled pool-failure masks
target:       every nine-body over the fixed pool
map:          keep only masks of rank 0, 1, or 2 and aggregate their widths
preserved:    exact retained Haar mass and every singleton/pair overlap
destroyed:    cyclic address, owners, cell order, and all rank >=3 masks
sidecar:      exact grid D, pair identity, and normalized margin ticks/D
hostile:      (50,70), body 031c7400
decisive test: 63L2(B)-4D>0 for the exact worst body on every residual pair.
```

The lost coordinates are harmless for this measure lower bound but cannot be
silently reused for arrival, ownership, or entry. That is why the graph
closes the fixed-pool chart while leaving the main conjecture open.

## 5. Consequence and remaining obstruction

THM-4231 already proves the fixed-pool two-outsider theorem off its finite
remainder, and the graph closes the remainder. Thus for every distinct
positive outsiders `q,r`, every nine-body `B` from the pool satisfies

```text
mu(G_(B union {q,r})) >= 4/63.
```

Haar invariance under multiplication also preserves this inequality for
`c(B union {q,r})`. THM-4150 then proves every row

```text
2c(B union {q,r}) union {a,b},
```

where `c` is positive and `a,b` are distinct positive odd integers. These are
genuine structured thirteen-speed families.

The main obstruction is now cleaner: **entry, not pair closure**. A general
normalized LRC(14) row need not have eleven even speeds whose halves consist
of nine labels in this fixed pool plus two outsiders. Nothing in the wall
graph creates such a parity normalization or pool embedding.

## 6. Generated next tasks

1. **Entry atlas.** Search for a small family of scale-equivalent pools whose
   `9+2` charts cover a provable normalization class of arbitrary rows. The
   output must be a typed map, not a numerical resemblance.
2. **Pool optimization.** Treat
   `inf_(q,r) min_(|B|=9) L2(B)/D` as the robust pool objective. Local changes
   should first be tested against the exact hostile `(50,70),031c7400` and
   then against the full finite ledger.
3. **Rank-`d` hypergraph compiler.** Implement `Ld` with exact maximum
   coverage and a dual certificate. Rank three may trade a larger packet for
   enough robustness to use smaller or more flexible pools.
4. **Analytic pair tail from graph statistics.** Seek bounds on `T` and the
   top weighted degrees that hold directly in `q,r`, replacing the finite
   pair ledger by an arithmetic inequality.
5. **Formal kernel.** Isolate the finite-null-wall lemma, graph identity,
   degree direction, and branch upper bound for Lean. The numerical ledgers
   would remain external exact certificates, but the logical compiler would
   become independently checkable.

The important shift is procedural: when a carrier is failing one layer at a
time, inspect the literal failure-rank distribution before enlarging the
carrier. A low-rank mass skeleton can be weaker pointwise yet stronger as a
uniform proof object.
