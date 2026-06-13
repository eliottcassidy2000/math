---
id: THM-391
name: lrc-zero-branch-star-core-peeling
status: PROVED
date: 2026-06-01
session: codex-2026-06-01-S548
depends_on:
  - THM-359
  - THM-369
  - HYP-2036
---

# THM-391: A single LRC zero-branch star has empty endpoint core

> **Renumbered from THM-390 by monad-reviewer-2026-06-01 (QC).**
> codex-S547 and codex-S548 independently both claimed THM-390 — a
> namespace collision of two distinct, both-PROVED theorems. The first
> claimant (S547's `lrc-padic-zero-branch-cover-core`, committed in
> fa44a9d) keeps THM-390; this S548 star-peeling result is renumbered to
> THM-391. The mathematics is unaffected. See MISTAKE-052.

## Statement

Let `n >= 2` and `2 <= q <= n`.  Let `C` be a finite subset of the nonzero
q-grid points

```text
C subset {1/q, 2/q, ..., (q-1)/q} subset R/Z,
```

and let `S` be a finite multiset of positive speeds, each divisible by `q`.
For every center `c = u/q` in `C` and every speed `s` in `S`, form the open
danger interval

```text
I(c,s) = (c - 1/(n s), c + 1/(n s)).
```

Use strict endpoint protection: an endpoint is protected only when it lies in
the interior of another active interval.

Then no nonempty subfamily of these intervals is an endpoint-protection core.
Equivalently, the THM-359 endpoint peeling process deletes every interval.

More sharply, if `m_s` is the multiplicity of speed `s` in `S`, then the peel
layers are explicit.  Speeds peel in increasing speed order:

```text
layer(s) = |C| * m_s.
```

In the LRC p-adic application, take `q=p^d` and `C` to be the unit points
`u/q`, `(u,q)=1`.  A covered prime-power zero branch kills the THM-369 unit
witnesses, but it cannot itself hold a local endpoint core.

## Proof

Write

```text
r_s = 1/(n s).
```

Since every `s` is divisible by `q`, we have `s >= q`, hence

```text
r_s <= 1/(n q) <= 1/(2q).
```

Distinct q-grid centers are separated by at least `1/q`.  Therefore intervals
with different centers are disjoint as open intervals, except for possible
endpoint touching in the equality case.  Endpoint touching is not strict
protection, so an endpoint of an interval centered at `c` cannot be protected
by an interval centered at a different q-grid point.

Now suppose a nonempty endpoint-protection core existed.  Choose an interval in
the core with maximal radius, equivalently with minimal speed among the
intervals in the core.  At its own center, every other interval has radius at
most this radius.  An equal-radius interval has the same endpoints and does
not strictly contain them; a smaller-radius interval misses those endpoints.
Intervals at different centers cannot strictly contain those endpoints by the
separation argument above.

Thus the chosen interval has both endpoints unprotected, contradicting the
definition of a protection core.  Hence no nonempty core exists.

For the layer formula, the same argument applied to the whole active family
says exactly the intervals with maximal active radius are removed first.  These
are precisely the intervals belonging to the smallest active speed `s`, one for
each center and each multiplicity of `s`.  Removing them leaves the same
q-grid star with larger speeds.  Induction gives the stated peel layers.

## Consequences

1. The S546b computation is now a theorem: the local prime-power zero-branch
   cores are empty for structural reasons, not only in the audited rows.
2. HYP-2036's product-tree branch ledger is a sieve-survival ledger, not a
   local counterexample-core ledger.
3. A genuine LRC endpoint core, if one exists, must mix at least one additional
   label: descendant endpoints, off-grid event centers, endpoint owners,
   cross-prime product-tree coupling, or compactified wall data.
4. HYP-2037/HYP-2038 entropy quantities measure global spread and criticality;
   THM-391 says they cannot be converted into a local zero-branch proof without
   retaining the exported descendant labels.

## Verification

The exact verifier checks the theorem and the peel-layer formula over bounded
q-grid stars, including non-prime-power q-grid examples and the recent n=14/n=18
LRC-style stars:

```text
04-computation/lrc_zero_branch_star_theorem_s548.py
05-knowledge/results/lrc_zero_branch_star_theorem_s548.out
```

The bounded verifier reports:

```text
bounded exact cases checked: 3255
max intervals in bounded check: 21
all bounded cores empty and all peel layers match the formula
```

The proof above is independent of the verifier.

**QC verification (monad-reviewer-2026-06-01):** Re-derived from definitions. The
key bound `r_s = 1/(ns) ≤ 1/(nq) ≤ 1/(2q)` (using `s ≥ q`, `n ≥ 2`) with center
separation `≥ 1/q` forces intervals at distinct q-grid centers to be at most
endpoint-touching; strict protection then rules out cross-center protection. In any
nonempty subfamily the minimal-speed (maximal-radius) interval has both endpoints
unprotected — same-center smaller radii miss them, equal radii do not strictly
contain them — so no nonempty core survives. Checked the boundary case `n=2, s=q`
(adjacent stars touch at a shared endpoint, not strict). Proof CONFIRMED.
