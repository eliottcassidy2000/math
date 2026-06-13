# Lonely Runner Endpoint Protection

**Session:** codex-2026-05-30-S359
**Status:** exploratory progress; exact rational probe added

## Starting Point

The previous Lonely Runner pass reframed the reduced conjecture as a finite
open-cover problem.  For `k` nonzero speeds `V` and threshold `1/(k+1)`, each
speed forbids

```text
F_v = {t in R/Z : ||v t|| < 1/(k+1)}.
```

The conjecture says the open union `F(V)=union_v F_v` is never all of the
circle.

There are three cases:

```text
positive gap       -> open interval of lonely witnesses
boundary-only      -> tight examples; witnesses are endpoints
full open cover    -> counterexample
```

The new object is the endpoint-protection graph.

## Endpoint Protection

Every forbidden interval endpoint has exact rational form

```text
t = (m +/- 1/(k+1))/v.
```

Call an endpoint protected if it lies strictly inside another forbidden
interval.  Then:

```text
full open cover  => every endpoint protected.
```

Conversely, if total forbidden length is `1` and some endpoint is unprotected,
then that endpoint is a boundary lonely witness.  This makes a counterexample
a finite incidence object:

```text
endpoints -> intervals that strictly protect them.
```

The script `04-computation/lonely_runner_endpoint_protection_s359.py` computes
this incidence graph exactly over `Fraction`.

Stored output:

```text
05-knowledge/results/lonely_runner_endpoint_protection_s359.out
```

## Protection Core Peeling

The session added an explicit endpoint-core peel.  An interval is allowed to
protect an endpoint only if both endpoints of that protecting interval are
still present.  Then repeatedly delete every endpoint with no surviving
protector.

```text
active endpoints E
active intervals I with both endpoints in E
delete e in E if no interval in I strictly contains e
iterate to a fixed point
```

A full open cover would have every forbidden endpoint strictly protected, so
it would survive as a nonempty core.  Thus a theorem that all primitive speed
sets peel to the empty core would prove LRC in this reduced form.

The bounded scan found no nonempty cores, including inside the necessary
counterexample subfamily where at least one speed is divisible by `k+1`:

```text
k=3, max_speed=20: primitive 997,  mod-filter 563,  nonempty cores 0
k=4, max_speed=16: primitive 1745, mod-filter 1066, nonempty cores 0
k=5, max_speed=13: primitive 1281, mod-filter 819,  nonempty cores 0
k=6, max_speed=11: primitive 462,  mod-filter 252,  nonempty cores 0
```

The hardest empty-core examples are not shallow artifacts.  The longest peel
in this scan used 17 rounds, for example

```text
(4,5,6,7,10,11): removed profile
(9,16,4,2,4,8,2,2,4,4,4,4,2,6,4,4,2)
```

This makes protection peeling a stronger object than simply finding an
initial exposed endpoint: some examples have high protected ratios and long
dependency chains, but the self-supporting core still collapses.

## Main Observations

### 1. Known tight examples are not close to all-protected

For all standard tight examples checked, there are explicit unprotected
endpoints.  Examples:

```text
(1,2,3,4):          4 unprotected forbidden endpoints
(1,3,4,7):          4 unprotected forbidden endpoints
(1,3,4,5,9):        2 unprotected forbidden endpoints
(1,4,5,6,7,11,13):  4 unprotected forbidden endpoints
```

So the known tight stratum is not an almost-counterexample in the protection
graph.  It is a rigid boundary skeleton.

### 2. The boundary skeleton is the unit group modulo `k+1`

In every known tight example checked, the unprotected endpoints are exactly

```text
{a/(k+1) : 1 <= a <= k and gcd(a,k+1)=1}.
```

The lonely witnesses are these nonzero unit residues `a/(k+1)`.  S362 later
clarified the earlier shorthand: `0` belongs to the Dirichlet-pigeonhole orbit,
but it is a forbidden center, not a forbidden endpoint.

This holds for:

```text
initial k=3,4,5,7
sporadic {1,3,4,7}
sporadic {1,3,4,5,9}
sporadic {1,4,5,6,7,11,13}
sporadic {1,2,3,4,5,7,12}
```

This is much sharper than "the first witness is `1/(k+1)`."  Tight examples
appear to collapse to the quotient `Z/(k+1)Z`, and the surviving boundary
set is the unit group.

### 3. The first counterexample filter is quotient-divisibility

If no speed is divisible by `k+1`, then `t=1/(k+1)` is immediately lonely:

```text
v not == 0 mod (k+1)
=> ||v/(k+1)|| >= 1/(k+1).
```

So every counterexample must contain at least one speed divisible by `k+1`.
This is elementary, but it is now part of the protection-graph story: a
counterexample must first protect the unit-boundary skeleton.

### 4. Near-tight positive-gap examples still have many unprotected endpoints

The near-tight scans rank positive-gap sets by `max_gap/threshold`.
Even the tightest small positive-gap examples retain visible unprotected
endpoint sets.  For example:

```text
k=6, speeds=(1,4,5,6,7,11)
gap/threshold = 0.023810
protected_ratio = 0.923077
unprotected endpoints = 5
```

This is the closest thing in the scan to "almost all protected," and it still
has a small explicit survivor set.

## Emerging Proof Shape

The endpoint-protection formulation suggests a possible proof strategy:

```text
1. Any counterexample must contain a speed divisible by n=k+1.
2. Divisible speeds try to protect the unit-boundary skeleton a/n.
3. Protection of unit endpoints creates new dependency endpoints elsewhere.
4. A denominator/interval potential should decrease along protection chains.
5. Peeling then forces the endpoint core to vanish.
```

This is analogous to endpoint-transfer private pivots in the tournament repo:

```text
unprotected endpoint = private witness
all-protected core   = collision block
proof                = leaf-peeling / no-core theorem
```

## New Conjectural Handles

### Unit-boundary skeleton

Tight instances may be exactly those where the protection graph collapses to
the unit-boundary skeleton.  This is not literally a classification yet, but
it is the strongest small-data pattern.

### Protection-peeling

A full open cover would require a nonempty all-protected endpoint core.  The
new finite theorem to pursue is:

```text
primitive speed set => endpoint-protection core is empty.
```

The scan supports this through the bounded boxes above.  The hard examples
suggest the proof cannot be only "find one exposed endpoint"; it needs a
monotone charge, denominator height, interval endpoint order, or transfer
argument that explains why every protection cycle eventually leaks.

### Divisibility descent

Speeds divisible by `k+1` are necessary for a counterexample.  If many are
divisible, divide them by `k+1` and use lower-dimensional LRC on that subset;
if few are divisible, the non-divisible speeds may be unable to protect all
unit endpoints.  This is the same split as the earlier failed induction, but
now stated in endpoint-protection language.

## Next Experiments

- Optimize endpoint protection so larger mod-filter scans are feasible.
- For every endpoint, record protecting speed chains and build a directed
  protection hypergraph with removal depths.
- Search for a denominator-height potential that strictly decreases after
  compressing protection chains.
- Run larger mod-filter scans restricted to sets containing a multiple of
  `k+1`, since other sets expose `1/(k+1)` immediately.
- Test whether boundary-only examples beyond the known ones always have the
  unit-boundary skeleton.
