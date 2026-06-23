---
id: HYP-2908
status: PROOF-TARGET / transfer guardrail; graph atom exact, LRC state-lift open
source: codex-2026-06-22-S118
tags: [lrc14, tournaments, digraphs, forbidden-H, binary-arcs, state-lift, conflict-graphs, open-q-108]
depends_on:
  - HYP-2905
  - HYP-2906
  - HYP-2907
  - HYP-2881
  - HYP-2882
  - HYP-2883
  - HYP-2677
  - HYP-2691
  - HYP-2648
  - HYP-2632
  - THM-002
  - THM-200
  - THM-201
  - THM-343
  - THM-572
related:
  - KPS-S31y
  - THM-344
  - OPEN-Q-108
results:
  - 04-computation/lrc_h7_state_lift_codex_s118.py
  - 05-knowledge/results/lrc_h7_state_lift_codex_s118.out
  - 04-computation/lean/TournamentH7/TournamentH7/LRCTournamentStateLift.lean
---

# HYP-2908: LRC14's forbidden-H7 route needs a state lift, not just a digraph

The owner's prompt is almost a proof, but only after one missing functor is
made rigorous.

The two usable tournament facts are:

```text
arc state = one of two orientations on every unordered pair
H(T) = I(Omega(T),2)
7 is forbidden as an H-value
```

The LRC14 number is exactly the tempting scalar:

```text
14 = 2 * 7.
```

But the scalar identity is not enough.  The proof target is the following
state-lift theorem.

```text
Every primitive top-balanced LRC14 counterexample produces a connected
binary packet graph G_LRC with I(G_LRC,2)=7, and G_LRC is tournament
conflict-realizable, or at least satisfies the THM-201/THM-343 closure rule
that a connected K3 conflict atom forces extra odd-cycle mass.
```

If this lift is proved, LRC14 is impossible to disprove by a counterexample:
`I(G_LRC,2)=7` forces `G_LRC=K_3`, while the tournament conflict category
forbids the `K_3` atom.

Codex S128/THM-572 now formalizes this closure endpoint in Lean.  The theorem
`LonelyRunner.not_bad_of_tournament_state_lift` says that any bad-atom
predicate is empty once each bad atom constructs a `TournamentStateLift` whose
activity-two packet value is `7` and agrees with tournament `H`.  Thus the
tournament route's remaining work is exactly the first arrow:

```text
LRC14 bad atom  ->  TournamentStateLift.
```

## Exact atom

For a connected simple graph `G` on `n` vertices,

```text
I(G,2) >= 1 + 2n
```

from the empty independent set and all singleton independent sets.  Hence
`n>=4` already gives `I(G,2)>=9`.  For `n=3`, the only connected graph with no
independent pair is `K_3`; all other connected graphs have an independent pair
and contribute at least `1+6+4=11`.  Therefore

```text
connected I(G,2)=7  <=>  G=K_3.
```

This is the rigorous version of `I(K_r,2)=1+2r` at `r=3`.  It also explains
why the obstruction is extremely low-level: from value `11` onward, non-clique
connected preimages already exist.

For tournaments, THM-002 gives `H(T)=I(Omega(T),2)`.  THM-343 is the modern
all-`n` proof that `H(T)!=7`; THM-200 records the OCF/K3 route.  If a
conflict component contributed value `7`, the component would be `K_3`, and
THM-201 forbids `K_3` as a connected component of `Omega(T)`.  Equivalently,
the `K_3` atom cannot remain isolated: three pairwise vertex-conflicting
primitive odd cycles force extra odd-cycle mass, visible in HYP-2881 as the
directed `C_5` support layer.

This is a component statement, not a subgraph statement.  THM-344 is the
guardrail: at `n=8`, `H=63` occurs with `Omega=K31`, so `Omega` contains many
`K_3` subgraphs.  The forbidden object is a whole `I=7` atom/component.

## Digraph guardrail

The word "digraph" is dangerous here.  If arcs are allowed to be present or
absent independently, then `H=7` is easy.  The S118 audit finds a 4-vertex
directed graph with exactly seven Hamiltonian paths:

```text
(0,1), (0,2), (0,3), (1,2), (1,3), (2,0), (2,1), (3,0).
```

So the proof cannot use arbitrary binary digraphs.  It must use complete
binary orientations, or the equivalent odd-cycle conflict graph of a
tournament.  This matches HYP-2883: even graphs realize `K_3` and therefore
realize `I(G,2)=7`; the tournament scar is orientation-born and disappears if
the quotient forgets the wrong state.

## LRC state lift

After HYP-2906, a counterexample reaching the finite atom cannot have one
locally huge speed:

```text
v_max <= 13 * v_second.
```

After HYP-2905/HYP-2904, the remaining non-descending object is the bounded or
top-balanced covering atom.  HYP-2908 builds on HYP-2907 and proposes that the right vertices of the
finite atom are not runners.  Candidate vertices are:

```text
sector-pair obligations
wall-crossing events
cover-arc packets
exact-period phi packets
measured sector-state words
sector-transfer DP packets
support-six relation packets
packet-sign tournament atlas states
binary side states around the observer
```

The edge relation should record incompatibility of two packets inside a single
would-be covering certificate.  The activity `2` comes from the two local side
states, just as OCF counts each independent odd-cycle packet with two signs.

The missing theorem has four parts.

1. **Packet extraction:** a primitive LRC14 counterexample has a connected
   apex-7 packet core after top-large peeling and prime/dilation reductions.
2. **Activity two:** each primitive packet has exactly two compatible side
   states in the quotient that preserves the LRC witness predicate.
3. **Exact value:** the non-lonely covering atom has weighted independence
   value exactly `7`, not just a loose lower or upper bound.
4. **Tournament-realizability / closure:** the packet graph lies in the
   tournament conflict category, or at least obeys the THM-201/THM-343 rule
   that a connected `K_3` conflict atom forces extra odd-cycle mass.

Then the proof is immediate:

```text
counterexample
  -> connected G_LRC with I(G_LRC,2)=7
  -> G_LRC=K_3
  -> THM-201/THM-343 closure forbids the atom
  -> contradiction.
```

In the Lean endpoint this same closure is packaged as:

```text
counterexample
  -> TournamentStateLift with packetValue = H(T) = 7
  -> contradiction by H_ne_seven.
```

## What this does and does not prove

This does not prove LRC14 yet.  It isolates the exact missing bridge in the
S31y slogan

```text
apex-7 over-cover  <=>  forbidden K3 conflict graph.
```

The S118 contribution is the guardrail around that slogan.  If the bridge lands
in arbitrary digraphs or even graphs, `7` is realizable and there is no
contradiction.  If it lands in tournament conflict packets, `7` is impossible.

After THM-572, the contradiction side is no longer a prose handoff.  The live
mathematical target is a construction theorem for the `TournamentStateLift`
object from the remaining `|M14|<=6` finite-core / shell-height atom.

So the route is better read as "prove that an LRC14 disproof is impossible"
rather than as a construction of a disproof.  A disproof would have to realize
the one connected graph atom that the tournament category forbids.

## Tournament Analysis / assumption challenge

Candidate vertices considered:

```text
runners
pairwise runner arcs
complete tournament orientations
arbitrary directed arcs
odd cycles
Omega conflict vertices
sector-pair obligations
cover arcs
wall crossings
exact-period phi packets
measured sector-state words
sector-transfer DP packets
support-six relation packets
packet-sign tournament atlas states
bounded-core proof obligations
```

Chosen quotient: primitive binary packets before scalarization.  This quotient
preserves conflict incidence and the activity-2 state count needed for
`I(G,2)`.  It destroys raw speed identity and raw denominator identity, which
is acceptable only if the LRC witness predicate is preserved by the packet
extraction.

Challenged assumption: "binary digraph" is enough.  It is not.  Present/absent
digraphs realize `H=7`; complete orientations and tournament odd-cycle
conflict graphs are the necessary carrier.

## Artifacts

- `04-computation/lrc_h7_state_lift_codex_s118.py`
- `05-knowledge/results/lrc_h7_state_lift_codex_s118.out`
- `07-reflections/lrc14-forbidden-h7-state-lift-codex-s118.md`
