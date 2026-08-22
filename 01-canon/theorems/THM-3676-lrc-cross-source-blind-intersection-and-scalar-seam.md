---
id: THM-3676
title: "LRC cross-source blind intersection and scalar seam"
status: >
  PROVED; PENDING INDEPENDENT HOSTILE AUDIT.  For the three possible
  selected blocker sources in one owner-pivot row, two distinct sources
  already cut the common all-packet blind space to zero when the six unit
  residues have nonzero sum, and to the global scalar line when that sum is
  zero.  Consequently eight explicit packet charts detect any one common
  phase-cone residue measure off the stated hostile.  Actual owner currents
  change at handoff, and no common-current transport, scalar-line exclusion,
  covering-row consequence, or LRC(14) conclusion is proved.
source: kps-s195 / cross-source owner-packet seam, 2026-08-21
depends_on:
  - THM-3671-lrc-all-packet-blind-subspace-and-phase-cone-gate
related:
  - THM-2305-canonical-blocker-word-handoff-hypergraph
  - THM-2337-expiration-word-residue-invisibility-and-first-bockstein-sidecar
  - THM-3671-lrc-all-packet-blind-subspace-and-phase-cone-gate
---

# THM-3676 -- two owner roles leave only the scalar seam

**PROVED; PENDING INDEPENDENT HOSTILE AUDIT.**  This is a finite linear
statement in one fixed mod-thirteen row.  Its conditional measure corollary
requires literally the same grouped residue measure in both source charts.
It does not identify the distinct semantic currents produced by two owner
strata.

## 1. Three source roles in one ambient relation kernel

Work over `F=F_13`.  Let

```text
U={0,1,2,3,4,5},                 B={1',2',3'}       (1)
```

be the six guard/unit labels and the three blocker labels in one fixed
THM-2309 row.  Write the speed-residue vector as `w in F^9`, with

```text
w_u!=0 for u in U,              w_j=0 for j in B,  (2)

K=w^perp,                       sigma_U=sum_(u in U)w_u.  (3)
```

For a selected source `j in B`, let `{a,b}=B minus {j}`.  For every ordered
pair of distinct grafts `k,l in U`, the THM-3666 target-dual plane has
zero-target fibre

```text
L^j_(k,l)={r in K:r_a=r_k and r_b=r_l}.            (4)
```

Put

```text
W_j=intersection_(k!=l)L^j_(k,l).                  (5)
```

THM-3671 computes this space for one source.  Repeating its equality-graph
argument in the common ambient coordinates gives the especially convenient
form

```text
W_j=span(e_j),                    if sigma_U!=0,

W_j=span(1,e_j),                 if sigma_U=0,      (6)
```

where `1` is the vector with all nine coordinates equal to one.  Indeed all
coordinates outside `j` must equal one scalar `d`, the source coordinate is
free, and the only remaining kernel equation is `d sigma_U=0`.

## 2. Exact cross-source intersection

Let `j,j' in B` be distinct.  If `sigma_U!=0`, (6) gives two distinct
coordinate lines, so

```text
W_j intersection W_(j')={0}.                       (7)
```

If `sigma_U=0`, write an element of the intersection in the two ways

```text
alpha 1+beta e_j=gamma 1+delta e_(j').             (8)
```

Comparison at any unit coordinate gives `alpha=gamma`.  The remaining
identity `beta e_j=delta e_(j')` forces `beta=delta=0`.  Hence

```text
W_j intersection W_(j')=span(1).                   (9)
```

Intersecting the third source changes neither answer.  Thus the pure-source
hostile of THM-3671 is genuinely source-relative: one owner change kills it.
The scalar hostile in the zero-sum branch is genuinely source-independent.

The four explicit charts

```text
(k,l)=(0,1),(2,3),(4,5),(1,0)                     (10)
```

already cut out `W_j` for each source.  Therefore eight charts, four for
each of two distinct sources, cut out exactly (7) or (9).

## 3. Conditional common-measure detector

Let `mu` be one nonzero finite complex measure on `K`.  Assume every nonzero
coefficient lies in one strict open half-plane through the origin.  For each
of the eight charts push `mu` through its two target differences as in
THM-3671.  Sums inside a nonempty fibre cannot cancel under the half-plane
hypothesis.  Therefore, if every pushforward were supported only at target
zero, every point of `supp(mu)` would lie in the intersection (7) or (9).

It follows that some one of the eight charts has a nonzero target aggregate
under either of the following sharp alternatives:

```text
sigma_U!=0 and 0 notin supp(mu);                    (11)

sigma_U=0 and supp(mu) intersection span(1)=empty.  (12)
```

In particular, all-coordinate-unit support supplies the premise in (11),
but not in (12): every nonzero scalar vector has all nine coordinates
nonzero.  The scalar seam is therefore not an artefact of forgetting unit
support.

## 4. Exact remaining bridge

THM-2305 gives positive-measure blocker-word handoffs, but the current on an
incoming subset of `E_j` and the current used for an outgoing subset of
`E_j` need not agree and can live on disjoint physical subsets.  Equations
(7)--(12) cannot be applied to those two currents without a transport map.
Likewise THM-2337's first word jet is decomposition-coloured and is not a
function of the uncoloured relation address; it does not presently exclude
the scalar line in (9).

The exact next interface is consequently:

```text
owner handoff
  -> retain one nonzero grouped current in two source-role charts
  -> use (7), or exclude the scalar residue line in (9)
  -> obtain a nonzero target aggregate.                         (13)
```

This theorem proves only the middle finite-algebraic gate.  It does not
prove the first arrow, a phase cone for the actual complex coefficients, or
the scalar-line exclusion.  No scalar row is excluded and LRC(14) remains
open.  **QED.**
