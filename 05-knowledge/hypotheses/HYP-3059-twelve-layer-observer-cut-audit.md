---
id: HYP-3059
title: Observer-extension cut payload and the twelve-layer audit
status: SYNTHESIS / exact finite audit and proof-carrier abstraction; not a proof
source: codex-2026-06-26-S223
tangent: T1141
script: 04-computation/observer_extension_cut_payload_codex_s223.py
result: 05-knowledge/results/observer_extension_cut_payload_codex_s223.out
related:
  - HYP-3056
  - HYP-3055
  - HYP-3054
  - HYP-3053
  - HYP-3052
  - HYP-3051
  - HYP-3050
  - HYP-3049
  - HYP-3048
  - HYP-3047
  - HYP-3046
  - HYP-3043
  - HYP-3039
  - HYP-3031
  - HYP-3013
  - HYP-3008
  - HYP-2120
  - HYP-2121
  - THM-381
  - THM-385
  - OPEN-Q-108
---

# HYP-3059: Observer-extension cut payload and the twelve-layer audit

## Claim

The user's `12` observation is real, but the arithmetic must be corrected:

```text
48 + 12 = 60, not 56.
```

The exact first-failure arithmetic is

```text
P(5) = 48
U(6) = 56
U(6) - P(5) = 8
U(6) = P(5) + U(5) - U(4) = 48 + 12 - 4.
```

This last identity is not promoted as a theorem by itself.  It is a payload
alarm: the `12` slice is a real extension boundary, while the `4` is the
previous parent/source layer that must be quotient-corrected when a
source/sink slice is blended with the arbitrary rooted-perspective cache.

As fractions of `U(6)=56`,

```text
P(5)/U(6)      = 6/7
defect/U(6)    = 1/7
U(5)/U(6)      = 3/14
U(4)/U(6)      = 1/14
1              = 6/7 + 3/14 - 1/14.
```

## Exact audit

Building on HYP-3056's observer-cut orbit ledger, HYP-3055's duodecimal bridge,
and the HYP-3054 observer-extension cut calculus, S223 computes the relevant small tournament surfaces without brute-forcing all
labelled six-tournaments; it reaches `U(6)` through the established
`U(5)` plus incident-word extension carrier.  The incoming S218
`observer_extension_payload_codex_s218.py` ledger is the broad payload-field
base layer; this S223 audit is the exact first-defect stress test for that
field list.

The exact `12` layers near the first defect are:

```text
P(4) rooted perspectives                         = 12
U(5) unrooted tournament classes                 = 12
5 -> 6 source-word extension slice               = 12
5 -> 6 sink-word extension slice                 = 12
SC(6) self-converse classes                      = 12
K_{4,5} rectangle cycle-space redundancy (S217)  = 12
```

The finite fiber checks add several useful details:

- Source and sink child sets in `5 -> 6` each have size `12`, but overlap in
  `4` classes.
- `SC(6)` has size `12`; it intersects the source slice in `2` classes and
  the sink slice in `2` classes.
- The deletion decks of the `12` self-converse six-classes touch all `12`
  five-vertex parent classes.
- Deletion decks for all six-classes have `52` unique decks across `56`
  classes, with four size-2 collisions.
- Ordered-pair sector decks still behave as in HYP-3049: size and internal
  sector decks separate `55/56`; cross-sector orientation separates `56/56`.

The source/sink overlap is the most concrete finite witness for why
`56 = 48 + 12 - 4` is the right splice to investigate: `4` is not arbitrary
mass, but a visible overlap between the two extremal observer slices at the
same surface where rooted perspectives first lose `8` unrooted classes.

## Observer-extension / cut payload

The proposed reusable carrier is:

```text
base quotient Q
observer/cut word sigma
stabilizer Aut(Q)
sidecar C(sigma)
sink map Phi(Q,sigma)
legality exit
```

The controlled-forgetting recurrence is:

```text
retained_{r+1} = (retained_r, C_r) / Aut(retained_r)
debt_r         = kernel(C_r -> LRC predicate)
```

A quotient is allowed to forget the payload only when the debt is
reconstructed, annihilated by a dual certificate, descended to a smaller
family, recognized as AP/Goddyn-Wong equality, or named as residual debt.

This packages several previous repo motifs in one language:

- A000568/rooted perspective: rooted node type forgets incident word and
  cross-sector chirality; ordered-pair sector decks repair it.
- LRC threshold packets: safe/open status forgets endpoint owner, period deck,
  and route label; labelled packet sheaves repair it.
- Residual capacitor cuts: coarse boundary type forgets which generator/cut
  blocks a route; finite min-cut payload and owner strips repair it.
- Haar drop/add squares: row/column margins forget mixed `zeta`; rectangle
  product residues repair it.
- S217 fixed-path diagonal flow: line-count and half-tiling shadows forget
  rectangle/hourglass cycle residues; `K_{k,k+1}` cut-space sidecars repair it.
- Matrix sidecar observability: scalar invariants forget hidden route/status
  coordinates; observability columns repair them.
- Moser/fibbinary automata: language membership forgets scale, bit-position
  parity, and carry boundary; automaton state plus valuation/topology sidecars
  repair it.
- Perfect-number/divisor lanes: abundancy or Euler-product scalars forget
  prime-power lanes; divisor-lattice or Smith-normal sidecars repair them.

## Tournament Analysis

Vertices are proof carriers, not runners or tournament nodes.  S223 orients
the carrier tournament by retained cut payload, exact finite checkability, LRC
predicate relevance, and automaton/handoff value, with proof cost as the tie
gauge.

The fingerprint is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}
directed_3_cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1]
hamiltonian_paths=1
```

One Hamiltonian path is:

```text
proof_obligation_automaton
> sidecar_observability_matrix
> endpoint_owner_packet
> rectangle_hourglass_residue
> deletion_fiber_payload
> ordered_pair_sector_deck
> incident_word_orbit
> rooted_node_cache
> raw_scalar_or_count
```

## Proof target

The next useful theorem is an observer-extension cut-payload theorem:

1. Start with a quotient that preserves an LRC predicate or tournament growth
   predicate.
2. Identify the observer/cut word and its stabilizer action.
3. Compute the sidecar needed to separate every mixed route/status fiber.
4. Prove the sidecar is retained, reconstructed, dual-annihilated, descended,
   AP/Goddyn-Wong boundary equality, or named residual debt.

This theorem would subsume the first A000568/rooted-perspective defect,
HYP-3049 cross-sector orientation, HYP-3052 deletion-parent fibers, HYP-3053
rectangle/hourglass residues, and HYP-3039 controlled forgetting in one
proof-carrier protocol.

## Assumption challenge

Considered vertices included runners, rooted nodes, directed edges, ordered
pairs, incident words, deletion fibers, source/sink slices, self-converse
branch loci, rectangle/hourglass cycles, endpoint owners, automaton states,
divisor lanes, matrix columns, and proof obligations.

The selected quotient uses proof carriers with explicit observer/cut payload.
It preserves LRC safe-box route/status information and tournament extension
recoverability.  It destroys labelled runner identity, raw word order, full
half-tiling presentation, and scalar-only counts unless those are retained as
sidecars.
