---
id: HYP-2748
title: LRC14 Delsarte/Tanner carriers are signed dessins with Ferrers quotients, not weakly regular LDPC-like graphs
status: OPEN guardrail; exact finite audit
source: codex-2026-06-21-S73
depends_on:
  - HYP-2740
  - HYP-2746
  - HYP-2745
  - HYP-2742
  - HYP-2741
  - HYP-2739
  - HYP-2726
related:
  - HYP-2738
  - HYP-2727
  - HYP-2723
  - HYP-2602
  - OPEN-Q-108
---

# HYP-2748: Signed Tanner/Dessin Weak-Regularity Audit

## Claim

The user's Tanner graph / unit-distance / weakly regular graph / Doyle-Holt /
Belyi prompt should be read as a carrier audit for the THM-534 Delsarte duals,
not as permission to import a sparse LDPC or unit-distance theorem.

For each binding Delsarte dual, define a bipartite carrier

```text
check vertices = nonzero moment rows r,
atom vertices  = missed-depth atoms t=0..6,
edge r--t      iff r <= t,
edge sign      = sign(y_r) in the binomial-basis certificate.
```

This carrier preserves the exact per-atom dominance predicate
`g(t) >= scale*[t=0]`.  It destroys the generated LRC speed word, the relation
code, and the L7 residue odometer.  The useful invariant is the signed
half-arc structure: support automorphisms never mix positive and negative
edge orbits.  The weakly regular and dessin quotients are address layers only.

## Exact Audit

Script:

```text
04-computation/lrc14_tanner_dessin_weak_regular_codex_s73.py
```

Stored output:

```text
05-knowledge/results/lrc14_tanner_dessin_weak_regular_codex_s73.out
```

Finite fingerprints:

```text
K8 : edges=25, density=5/7,  check-degree hist={7,6,5,4,3},
     atom-degree hist={1,2,3,4,5x3}, automorphisms=6,
     edge-orbits=15, mixed-sign-orbits=0,
     dessin black degrees=[3,4,5,6,7],
     white degrees=[1,2,3,4,5,5,5],
     face lengths=[4,8,13].

K9 : edges=22, density=11/14, check-degree hist={7,6,5,4},
     atom-degree hist={1,2,3,4x4}, automorphisms=24,
     edge-orbits=10, mixed-sign-orbits=0,
     dessin black degrees=[4,5,6,7],
     white degrees=[1,2,3,4,4,4,4],
     face lengths=[4,6,12].

K11: edges=18, density=6/7, check-degree hist={7,6,5},
     atom-degree hist={1,2,3x5}, automorphisms=120,
     edge-orbits=6, mixed-sign-orbits=0,
     dessin black degrees=[5,6,7],
     white degrees=[1,2,3,3,3,3,3],
     face lengths=[7,11].
```

For all three carriers:

```text
dominance check: true
degree quotient: equitable
weak-regularity verdict: not biregular; common-neighbor counts form a chain
mixed-sign support orbits: 0
```

So the carrier has a clean Ferrers/equitable quotient but not a weakly regular
or LDPC-like local structure.  The Doyle-Holt analogy survives only in the
half-arc sense: the support has symmetries, but the sign/orientation classes
are not interchangeable.

## External-Source Reading

- Samorodnitsky's Delsarte LP optimum paper reinforces a guardrail: an LP
  optimum is a global relaxation certificate and need not construct the
  extremal object.
- Li's `Unique Optima of the Delsarte Linear Program` is directly useful:
  the extend/puncture parity phenomenon for quasicodes matches HYP-2740's
  K8/K9/K11 parity target.
- The antipodal-contact Delsarte method (`math/0306405`) points to spectral
  positivity over orthogonal polynomial bases, not to a sparse Tanner graph.
- The MathOverflow Tanner/Delsarte question asks whether Delsarte extremizers
  guide Tanner or LDPC constructions.  In this LRC carrier they guide a signed
  finite certificate, but the graph is dense and Ferrers-chain-like.
- The IACR ePrint PDF in the prompt was blocked by Cloudflare in this shell,
  so this hypothesis does not rely on its contents.
- Incoming HYP-2746 tests Tanner graphs on the relation code `Lambda(E)` and
  finds the same warning from a different carrier: girth, expansion, and
  absolute enumerators do not give the LRC cap; signed weight distribution
  does.  HYP-2748 is the parallel audit for the THM-534 Delsarte dual carrier.
- Incoming HYP-2745 completes the general-prime residue discrepancy program;
  HYP-2748 therefore keeps Belyi/dessin language as an address layer rather
  than a replacement for the signed Delsarte proof stack.

## Proof Order

The remaining LRC bottleneck is not an L7 discrepancy constant: HYP-2739 proves
the exact residue numerator, and HYP-2741 repairs the finite-`f1` convergence
rate.  The honest remaining gaps are the r>=3 pairwise lift, base-size
domination, and the upstream consec-max / generated-depth extremality problem.

HYP-2748 therefore points to this proof order:

```text
generated depth word
  -> signed Delsarte/quasicode parity
  -> HYP-2738 aggregate/consec-max ledger
  -> cap scalarization.
```

Weakly regular graphs, unit-distance graphs, and Belyi passports are useful
only if they retain the signed Delsarte orientation.  The unsigned quotient
forgets exactly the information that makes `g(t) >= scale*[t=0]` true.

## Tournament Analysis

Vertices are proof lenses:

```text
generated_depth_occupancy
> signed_delsarte_dual
> delsarte_quasicode_parity
> halfarc_sign_orbit
> weak_regular_equitable_quotient
> dessin_belyi_passport
> sparse_tanner_expander
> unit_distance_spectral_transfer.
```

Pairwise observable:

```text
(preserves_bound, exact_audit, LRC_remaining_gap, sign_structure, locality)
```

Switch/gauge: generated words and signs are retained before scalarization.
The computed tournament is transitive after strict tie-break, with no directed
3-cycles.

## Assumption Challenge

The useful vertices are not runners, arcs, unit-distance points, or Tanner
checks alone.  Moment rows and missed-depth atoms preserve the Delsarte
dominance predicate, but they destroy generated-word information.  The
quotient is acceptable as an audit layer only if it is handed back to the
generated-depth proof stack before any final cap claim.
