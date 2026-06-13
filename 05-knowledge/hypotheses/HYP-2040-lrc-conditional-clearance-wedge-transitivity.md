---
id: HYP-2040
status: OPEN
source: codex-2026-06-01-S549
related:
  - HYP-1981
  - HYP-1992
  - HYP-2034
  - HYP-2036
  - HYP-2039
  - HYP-2041
  - HYP-2042
  - THM-381
  - THM-389
  - THM-395
---

# HYP-2040: LRC conditional-clearance cascades carry a hidden transitive wedge debt

**Claim.** The S548 LRC cascade should be read as a product of conditional
clearances, but each accepted clearance also propagates an edge-local
transitivity exclusion.  In a tournament, if `X -> Y`, then a transitive
extension cannot contain a third vertex `Z` with both

```text
Z -> X  and  Y -> Z.
```

Equivalently, an accepted edge `X -> Y` forbids the backward wedge
`Y -> Z -> X`.  This is the hidden local fact that propagates alongside the
usual transitive closure rule

```text
X -> Y and Y -> Z imply X -> Z.
```

For an LRC cascade, the product factor

```text
P_k = P(runner k clears the observer | previous runners cleared)
```

therefore has two jobs:

1. It records the ordinary safe-measure shrink.
2. It exports a wedge debt to every remaining vertex/object that would complete
   a backward two-step path across the newly accepted clearance edge.

The conjectural proof-bearing object is not merely the scalar product
`prod_k P_k`, but the product together with the accumulated no-backward-wedge
obligations.

## Exact tournament identity

THM-395 canonizes the exact tournament part of this hypothesis.  For a
tournament `T`, define the backward-wedge set of an oriented edge:

```text
W_T(X,Y) = { Z : Z -> X and Y -> Z }.
```

Then:

```text
T is transitive  iff  W_T(X,Y) is empty for every edge X -> Y.
```

Moreover,

```text
sum_{X -> Y} |W_T(X,Y)| = 3 * c_3(T),
```

because each directed triangle contributes exactly one backward wedge to each
of its three edges.  The user's hidden fact is this identity in local form: an
edge in a transitive tournament is not passive; it deletes a whole class of
future placements.

Incoming S545/S545o work gives the right boundary on this exact identity.
HYP-2041 reads the same local rule as the no-return or 3-term resonance
obstruction in the conditional-clearance product.  HYP-2042 shows why the
3-cycle/Helly-3 layer is necessary but not sufficient for the full LRC ladder:
the exact tournament fact is the base ledger, not by itself the complete LRC
certificate.

## LRC reading

The S548 product formula uses feasible sets

```text
F_k = { t : runners 1..k are all safe from the observer },
P_k = mu(F_k) / mu(F_{k-1}).
```

This is the right first-order cascade.  HYP-2040 says the dependence correction
should be organized edge-locally: when a clearance edge is accepted in an
observer-source or event-owner tournament, it also forbids all remaining
objects that would sit in the backward wedge.  Thus the resonance debt is not
only a pairwise overlap error; it is a third-object compatibility debt attached
to accepted edges.

In this language, the AP/regular row is hard because the ordinary clearance
credit survives until the last runner, while the accumulated wedge debt aligns
on the final wall.  Non-AP rows should spread the wedge debt across several
levels, leaving at least one positive or compactified clearance.

## Assumption challenge

The session explicitly considered these possible vertices for the tournament:

```text
runners;
danger arcs;
clearance obligations F_{k-1} -> F_k;
observer-source edges;
endpoint owners;
wall-crossing events;
fixed circle sections;
section boundaries;
residues and p-adic zero branches;
cover arcs;
Fourier/Gabor modes;
proof obligations.
```

The strongest current mapping is:

```text
vertices = clearance obligations, endpoint owners, or wall events;
edge X -> Y = X is already forced/cleared before Y in the observer-source
              or endpoint-pressure cascade.
```

Predicate preserved:

```text
the observer-tie-free/source target, equivalently zero danger occupancy at the
observer after all clearances.
```

Information destroyed:

```text
exact runner positions inside safe sectors, exact circular gap lengths, and
some arithmetic labels unless the p-adic/Gabor decorations from HYP-2036 and
HYP-2032 are retained.
```

Challenged assumption:

```text
the cascade only needs scalar conditional probabilities.
```

The hidden transitivity rule says scalar clearances also carry forbidden
third-object wedges.  A proof attempt that ignores those wedges risks treating
dependent clearances as independent.

## Predictions

1. Along AP initial rows, backward-wedge debt concentrates at the last runner
   or compact wall, matching the S548 last-runner bottleneck.
2. Non-AP rows should have positive final clearance because their wedge debt
   is distributed across earlier levels rather than aligned on one wall.
3. In endpoint-pressure language, a putative LRC counterexample core must
   realize a closed chain of wedge violations; otherwise the no-wedge
   exclusions peel the core.
4. Bare runner tournaments will often collapse to transitive shadows.  The
   useful cyclic signal should appear after lifting vertices to endpoint
   owners, wall events, Gabor zero columns, or p-adic branch obligations.
5. The S550 verifier starts this ledger: it reports `P_k` for sample LRC
   cascades and tournament fingerprints including backward-wedge mass,
   directed 3-cycles, SCC sizes, score histograms, and Hamiltonian-path counts.
   The remaining script work is to build the LRC-specific lifted tournament
   whose vertices are endpoint owners, wall events, or clearance obligations.

## Status

Open as an LRC lift claim.  THM-395 proves the exact tournament identity and
turns the hidden transitivity fact into a canon ledger.  HYP-2041 identifies
the same ledger with no-return/resonance debt, while HYP-2042 records the
honest limitation that fixed triple clearance is not enough for the order-`n`
cascade.  The remaining work is to choose the LRC lift where the backward-wedge
debt becomes a certificate rather than a visualization.

## Files

`01-canon/theorems/THM-395-backward-wedge-transitivity.md`;
`04-computation/lrc_conditional_clearance_wedge_s550.py`;
`05-knowledge/results/lrc_conditional_clearance_wedge_s550.out`;
`07-reflections/lrc-conditional-clearance-wedge-transitivity-s549.md`;
`07-reflections/lrc-cascade-as-conditional-clearances-the-no-return-fact-is-the-resonance-obstruction-s545o.md`;
`05-knowledge/hypotheses/HYP-2042-lrc-cascade-clearance-ladder.md`.
