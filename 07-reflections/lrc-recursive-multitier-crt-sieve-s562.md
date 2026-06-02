---
source: codex-2026-06-02-S562
status: reflection + exact packet audit
tags: [LRC, CRT, sieve, recursive-sieve, p-adic, endpoint-debt, n14, n17, n18, Tournament-Analysis]
---

# Recursive multi-tier CRT sieves

The useful shift is to stop picturing a CRT sieve as a stack of independent
prime filters.  A flat filter asks:

```text
does some q have no speed divisible by q?
```

That is THM-369's coarse sieve, and it is powerful.  But once all those gates
are covered, the problem has not disappeared.  The obstruction has moved.  A
multi-tier CRT sieve should therefore recurse on the residual object created
by a covered gate.

S561 already found the first algorithmic version for `k+1=2q`: the mod-2 tier
is forced, the mod-q tier is cheap, and the hard set is not searched but
generated as an apex-stuck or ratio-cover residual.  S562 asks what that
residual becomes if we keep pushing it down local denominator trees.

## The node

The right vertex is not a runner.  For this pass the vertex is the whole packet

```text
{1} union {scale*q : 1 <= q < n, q != skip}.
```

The recursive state is:

```text
(n, skip, v_p(scale) over the active primes).
```

For `n=14`, the packet keeps the odd payload `7` and moves in the dyadic
coordinate:

```text
scale 7 -> 14 -> 28
v(scale): (2:0,7:1) -> (2:1,7:1) -> (2:2,7:1).
```

For `n=18`, the packet keeps the square payload `3^2`:

```text
scale 9 -> 18 -> 36
v(scale): (2:0,3:2) -> (2:1,3:2) -> (2:2,3:2).
```

For `n=17`, the packet keeps the prime payload `17`, but the same dyadic
endpoint translation is visible:

```text
scale 17 -> 34 -> 68
v(scale): (2:0,17:1) -> (2:1,17:1) -> (2:2,17:1).
```

That last line is the pleasant surprise.  The dyadic tier can be introduced by
the residual endpoint denominator even when `2` is not a base factor of `n`.
So multi-tier means product-tree recursion on residuals, not just CRT factors
of the original denominator.

## The exact invariant

The S562 audit computes exact S356 interval reports.

```text
n=14: gap*boundary = 5/11
n=17: gap*boundary = 225/136
n=18: gap*boundary = 1
```

Each dyadic lift has the same local form:

```text
gap_factor = 1/2
boundary_factor = 2
product_preserved = True
```

That does not prove LRC.  It says something more targeted: once a residual
packet exists, further local lifts do not erase the obstruction.  They trade
visible Archimedean gap for exposed endpoint/frontier debt.

This is exactly the debt-export language of HYP-1844 and HYP-1866, now phrased
as a recursive CRT sieve.  HYP-1868 says raw boundary count should eventually
be replaced by a normalized product-building frontier mass.  S562 only tests
the raw invariant because that is the exact quantity already available from
S356.

## Why doubled primes still matter

The doubled-prime cases `2q` are not merely awkward composites.  They are the
first place where the cheap prime-field method loses the mod-2 freedom and the
apex becomes a zero-divisor.  But after S561, that failure is not a wall.  It
is a tier split:

```text
forced parity -> mod-q oracle -> generated apex/ratio residual.
```

The doubled prime contributes the first dyadic hinge.  The important doubled
prime is not just a number-theoretic term in Lemoine/Goldbach language; it is
the first even local coordinate where residual debt can be exported instead of
deleted.  That matches the earlier "doubled primes as parity hinge" reading,
but now the hinge has a computational shape.

## Tournament Analysis

Vertices:

```text
whole residual packets / proof obligations
```

Observable:

```text
(gap/th, gap*boundary, valuation depth, boundary count)
```

Switch:

```text
smaller visible gap wins; if tied, prefer lower conserved debt and deeper
translation.
```

Fingerprint:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3_cycles=0
sccs=[1,1,1,1,1,1,1,1,1]
hamiltonian_paths=1
```

This is not a failure of Tournament Analysis.  It says the bare valuation
recursion is a transitive ledger.  To get meaningful cycles we probably need
the labels that the quotient discards: endpoint owners, pressure leaves,
wall-crossing order, or signed cross-prime couplings.

## Assumption challenge

Considered vertex sets:

```text
runners, gaps, fixed sections, section boundaries, wall-crossing events,
residues, cover arcs, Fourier/Gabor modes, matroid circuits, CRT channels,
p-adic product-tree nodes, skipped gate labels, endpoint owners, and proof
obligations.
```

Chosen quotient:

```text
vertices = residual gate-packet proof obligations.
```

Preserved:

```text
sieve completeness, skipped label, valuation translation, exact gap ratio,
boundary debt, and raw frontier product.
```

Destroyed:

```text
runner identity, exact interval adjacency, owner signs, pressure leaves, and
cross-prime cancellation data.
```

Challenged assumption:

```text
CRT sieve tiers are fixed once the base denominator n is fixed.
```

S562 suggests tiers are created recursively by residual denominators.  The base
CRT factorization is only the first floor of the building.

## Handoff

The next proof object should be a weighted recursive sieve state:

```text
ResidualPacket(
    denominator,
    skipped_label,
    valuation_vector,
    owner_multiset,
    product_frontier_mass,
    pressure_leaf_status,
)
```

Then a descent lemma can try to prove:

```text
local witness exists
or residual has positive product-frontier mass
or residual exports to children with conserved positive mass.
```

That would turn "multi-tier CRT sieve" from a fast filter into a proof route.

**Artifacts:** `04-computation/lrc_multitier_crt_sieve_s562.py`,
`05-knowledge/results/lrc_multitier_crt_sieve_s562.out`, HYP-2073.
