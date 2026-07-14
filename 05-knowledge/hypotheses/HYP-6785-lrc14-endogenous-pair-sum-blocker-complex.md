---
id: HYP-6785
title: The endogenous pair-sum blocker complex is exact at fixed S; a scale-normal blocker-pressure theorem is open
status: EXACT EQUIVALENCE PROVED via THM-668 / structural pressure hypothesis open
source: codex-2026-07-14-S2
related: [THM-366, THM-523, THM-668, THM-718, THM-722, THM-753, THM-758, THM-760, HYP-2443, HYP-2972, HYP-6280, HYP-6780, HYP-6815]
script: 04-computation/lrc14_endogenous_blocker_complex_codex_S2.py
result: 05-knowledge/results/lrc14_endogenous_blocker_complex_codex_S2.out
---

# HYP-6785 — Endogenous pair-sum blocker complex

## Exact construction

Let S={v_0,...,v_12} and fix a target delta=A/B.  For each unordered pair of
indices i<=j, put q=v_i+v_j.  For every multiplier p with 1<=p<q, including
nonunits, retain:

- the event t=p/q;
- the protected endpoint set E={i,j};
- the blocker set

  B(i,j,p)={h not in E: ||v_h p/q|| < delta},

provided the protected endpoints themselves satisfy
||v_i p/q||,||v_j p/q|| >= delta.

The exact object is a two-sorted incidence complex: one sort is proof
obligations (i,j,p), the other is runners, and incidence means “this runner
blocks this obligation.”  Its hypergraph shadow has runners as vertices and
obligation-labelled blocker sets as edges.  For proof search, the obligation
nodes—not the runner vertices—are the natural primary objects.

## Exact equivalence

THM-668 proves that the maximum of min_v ||vt|| occurs at a pair-sum event
p/(v_i+v_j).  Therefore

> M(S)>=delta if and only if the endogenous blocker complex contains an empty
> blocker edge.

This is an exact finite predicate for each fixed S.  It includes all pair-sum
rulers and all multipliers, so it does not conflict with THM-566's refutation
of a fixed external denominator window.

The construction is the exact completion of the older dynamic twist-ladder
blocker picture in HYP-2443/HYP-2972:

~~~text
old object:
  selected external twist ladder -> blockers

completed object:
  every endogenous THM-668 pair-sum event
    + protected endpoint owners
    + simultaneous blocker hyperedge.
~~~

## Structural hypothesis

For primitive covering thirteen-speed families at delta=1/14, one of the
following should hold:

1. an empty blocker edge exists;
2. a safe peel or a witness-sheet dodge reduces the family;
3. a bounded good-period certificate exists after scale normalization; or
4. repeated pair-sum rulers force a classified AP/GW/deep-well/near-dilate
   packet.

Equivalently, a scale-normal primitive covering blocker complex should not
cover every protected obligation unless its additive realization belongs to a
finite tight packet class.

The missing result is a local-to-global blocker-pressure theorem, not the
finite equivalence above.

## New scale lane closed

THM-760 proves that cP union {w}, with c>=2 and gcd(c,w)=1, has

M(cP union {w}) >= min(M(P), 1/2-1/(2c)).

For a twelve-speed core this is at least 1/13.  Thus the unbounded
codimension-one scaled-core lane exposed by HYP-6780 is already discharged.
The open scale residual has at least two exceptional residue classes or no
twelve-speed common-factor core.

## Computational atlas

The exact companion script:

- verifies THM-760 on 80 deterministic random twelve-core instances;
- checks the c=26,52,104 near-dilate ray and a safe-killer variant exactly;
- computes M from every pair-sum event;
- reports blocker-edge size histograms at 1/14, 14/183, and 1/13;
- reports pair-sum multiplicity; and
- runs the required Tournament Analysis as a deliberately lossy diagnostic.

Representative exact results:

| family | exact M | maximum pair-sum multiplicity | qualitative complex |
|---|---:|---:|---|
| AP {1,...,13} | 1/14 | 6 | tight, repeated rulers |
| GW | 1/14 | 5 | tight, repeated rulers |
| deep well | 14/183 | 6 | two empty edges at 14/183 |
| second rung | 7/89 | 5 | strict slack |
| near dilate c=26 | 1/13 | 6 | many empty edges; sheet proof |
| spread blocker example | 406/1669 | 2 | very loose, abundant empty edges |

This supports, but does not prove, an additive dichotomy:

> High pair-sum multiplicity forces a structured ruler/sheet packet; low
> multiplicity forces blocker sparsity, decorrelation, or a small good period.

## Tournament Analysis

**Vertices:** runners, for the diagnostic quotient only.

**Pairwise observable:** exclusive blocker responsibility across obligations
where neither competing runner is protected.

**Switch/gauge:** orient toward the runner with greater exclusive
responsibility.

**Tie Hamiltonian path:** increasing speed.

**Fingerprints reported:** score histogram, directed 3-cycles, SCC sizes,
Hamiltonian-path count, and edge flips between 14/183 and 1/13.

The tournament destroys simultaneous coverage intersections, exact margins,
and ruler labels.  It is not equivalent to the LRC predicate.  Alternate
vertex sets explicitly considered were runners, distinct rulers, residues,
wall events, witness sheets, matroid-like circuits, and proof obligations.
Proof obligations are the exact vertices; sheets are the correct vertices for
THM-760.

## Immediate proof pulls

1. Double-count private blockers against protected endpoint pairs.
2. Form the protected set-cover LP and study its dual obstruction weights.
3. Prove an inverse theorem from high pair-sum multiplicity to AP/symmetric
   cores.
4. Prove low multiplicity yields empty edges or a bounded good period.
5. Attach scale/gcd residue fibers and derive a well-founded quotient descent.
6. Relate safe-peel irreducibility to minimum degree or circuits in this
   complex.

The q<=25 arm is now being audited for uniformity by
05-knowledge/hypotheses/HYP-6820-q25-and-n12-uniformity-audit.md
(codex-2026-07-14-S3; renumbered from HYP-6810 by opus-S299 after the
collision with opus-S298's earlier assembly-paper reservation).  Its output
should attach good-period certificates to normalized blocker states.

The full ranked backlog and historical guardrails are in
00-navigation/LRC14-FRONTIER-AND-AVENUES-2026-07-14.md.

## Incoming affine and endpoint-sidecar connection

HYP-6815 represents V(c)=cP+R exactly as the slope-c geodesic fiber of the
two-torus threshold function min_i ||p_i u+r_i t||.  Its R=0 cylinder is
HYP-6780 scaling; nonzero R retains the exceptional offsets discarded by a
pure scale quotient.  The blocker incidence complex should be attached to
the strip strata of this suspension.

The exact endpoint-sidecar audit in
04-computation/lrc14_endpoint_tournament_sidecar_audit_structure_S1.py
supports the same controlled-forgetting rule: neither its runner tournament
nor its unweighted endpoint tournament preserves the theorem-facing
predicates.  Covering needs the divisor mask; the THM-755 decision needs the
projective cap ratio; Bernoulli discrepancy needs signed endpoint phases; and
exact M needs a peak witness/value sidecar.  These become explicit fields of
the proposed scale-normal blocker automaton.
