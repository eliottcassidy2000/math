---
id: HYP-1859
status: OPEN
source: codex-2026-05-31-S391
related:
  - THM-357
  - THM-365
  - THM-366
  - THM-367
  - HYP-1841
  - HYP-1842
  - HYP-1844
  - HYP-1857
  - HYP-1858
  - HYP-1866
---

# HYP-1859: n=16 needs a global dyadic debt-flow inequality

## Statement

At Lonely Runner denominator `n=16`, a full proof should use a global
dyadic endpoint-debt flow rather than a purely local private-endpoint lemma.

The target inequality is:

```text
primitive 15-speed all-protected endpoint system
+ at least one 16-gate
=> positive dyadic endpoint divergence,
```

where the divergence is measured over owner layers and protector `v2` layers.
Equivalently, after THM-366 forces a `16`-gate and THM-367 gives the exact
local endpoint-count law, the remaining branch should fail because the
fifteen available speeds cannot simultaneously:

1. close the old unit skeleton;
2. pay the nine-residue local bill at the first pure dyadic owner;
3. keep every descendant dyadic endpoint protected; and
4. avoid either a positive gap or a peelable endpoint leaf.

## Evidence

THM-367 proves the pure dyadic local count formula.  For owner `u=2^k`, a
protector residue `p=2^j q mod 16u` protects exactly:

```text
super-gate p=0 mod 16u:                   2u endpoints
j >= k, non-super:                         0 endpoints
k-j >= 3:                                  2^(k-2) endpoints
k-j = 2 and q in {+/-1,+/-3} mod 16:       2^(k-1) endpoints
k-j = 1 and q in {+/-1} mod 16:            2^k endpoints.
```

This turns the vague S389/S390 debt language into a precise local capacity
law.

For maximal-owner lower protection:

```text
u=2,4,8  cannot be closed by lower protectors;
u=16     has exact lower-cover number 9;
u=32,64  also have exact lower-cover number 9 in the bounded solver.
```

The important surprise is that only `u=16` has the clean private-endpoint
certificate.  For `u>=32`, lower protectors can close the local endpoint layer
self-similarly, and the obvious private-leaf proof disappears.  The obstruction
must therefore be global: a local closure at one dyadic height exports a debt
pattern to another height, or consumes protector residues needed elsewhere.

Odd-payload scans for `u=16w` with `w=1,3,5,7,9,15` show that the maximum
capacity of each `v2` layer is exactly `w` times the pure `u=16` capacity:

```text
v2(p)=0: max 4w
v2(p)=1: max 4w
v2(p)=2: max 8w
v2(p)=3: max 16w.
```

Counts within a layer vary once the odd payload is present, so the likely
formal statement is not a literal tensor product.  It should be a capacity
upper bound plus a residue-class compatibility constraint.

## Proof Route

The conjectural proof should build a weighted endpoint-flow ledger:

```text
owner layer        demand = number of uncovered strict endpoints
protector layer    capacity = THM-367 count, or odd-payload extension
flow edge          p protects endpoint of u with strict slack
divergence         demand - incoming capacity after quotient constraints
```

The weights should penalize protection that uses a maximal-speed super-gate,
because a maximal owner cannot be protected by a larger exact multiple.  They
should also remember half-turn parity: odd protectors spend their budget
one-sidedly, while deeper even protectors fold through a dyadic quotient.

The hoped-for theorem is a finite certificate:

```text
Every n=16 all-protected endpoint core has positive dyadic divergence.
```

Combined with THM-365, this would rule out labelled endpoint cycles in the
gated branch.  Combined with THM-366, it would prove the `n=16` Lonely Runner
case.

HYP-1866 reframes this target as a product inequality:

```text
ArchGap * Debt_2 >= c(16) > 0.
```

The raw dyadic ladders already show two exact plateaus, `34/33` and `35/33`,
with a `35/34` phase tax when the exposed layer changes.  This suggests that
the global debt-flow divergence should be normalized as a 2-adic debt size,
not merely as an endpoint count.

## Sources

- THM-367.
- `04-computation/lrc_n16_dyadic_endpoint_formula_s391.py`
- `05-knowledge/results/lrc_n16_dyadic_endpoint_formula_s391.out`
- HYP-1857.
- HYP-1858.
- HYP-1866.
