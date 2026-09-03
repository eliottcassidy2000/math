# Attachment address before quotient: a JC/LRC synthesis

**Status: RESEARCH SYNTHESIS / TASK GENERATOR, NOT A THEOREM.** The planar
Jacobian inputs below are now proved in THM-4352/4354/4355/4356. The bounded
LRC test generated here has also been run and is recorded as a FINITE-EXACT
sidecar; it proves no new LRC(14) row safe.

## The signal

Three nearby planar-Jacobian boundary pictures separate data that a simple
incidence graph, tournament, or support mask would conflate:

1. THM-4354 has two rational middle components and eleven separately
   addressed transverse nodes. Replacing the eleven parallel edges by one
   simple adjacency changes graph genus from ten to zero.
2. THM-4352's even cusp has two tail ends entering one connected complement.
   The two attachments add one graph cycle.
3. The THM-4355/4356 `A15` packet has two tail ends entering two different
   complementary components. The attachments join the components and add no
   graph cycle. THM-4355's `A3` returns to the THM-4352 behavior.

The reusable object is therefore not “two ends” but

```text
(normalized packet, normalized complement components,
 attachment edges, endpoint-to-component address map).
```

For finite graphs `P` and `C`, add `k` attachment edges without identifying
vertices. Euler's identity gives the exact bookkeeping formula

```text
b1(P union C plus attachments)
 = b1(P)+b1(C)+k-c(P)-c(C)+c(final).                  (1)
```

Thus a connected packet with two ends contributes `+1` when the complement
is connected, but `+0` when its two ends connect two complementary
components and the final graph is connected. Formula `(1)` explains both
outcomes without analogy or parity guessing.

## What the tournament quotient loses

A tournament has at most one directed edge between two vertices. It cannot
faithfully encode the eleven `mu_11`-addressed nodes of THM-4354, and orienting
one representative does not recover the ten erased cycles. A lawful
tournament view would need each *node occurrence* as a vertex or an explicit
parallel-edge/root-address sidecar. Either choice is a different object from
the two-component tournament.

This is a concrete instance of the repository guardrail “do not manufacture
tournaments from ties.” The intrinsic binary relations available here are
branch-to-component landing and first-resolution precedence. Component
adjacency alone is symmetric and should remain a labelled multigraph.

## Exact map to the LRC representation problem

The present LRC(14) proof map already calls for an owner-labelled
relation-overlap hypergraph. The closest precise dictionary is:

| Jacobian degeneration | LRC boundary word |
|---|---|
| normalized complement component | safe component / owner sector |
| marked tail end | located endpoint or arrival occurrence |
| attachment address | owner, side, clock, and occurrence address |
| several nodes with one component pair | parallel returns with the same coarse labels |
| graph connection versus new cycle | bridge-to-child versus genuine return |

The map preserves labelled incidence and connectivity. It destroys metric
tooth length, signed phase, chronology, current magnitude, and the loneliness
predicate. Consequently `(1)` supplies a representation audit for LRC, not a
covering argument.

The sharp warning is the same in both subjects: counting two appearances of
an endpoint says nothing until their target components are known. In LRC
language, two return-looking occurrences may join previously separate safe
components rather than close a usable owner cycle. Conversely, quotienting
parallel located returns to one simple edge can erase independent cycle
constraints.

## Reversible natural addresses

THM-4354's 64 exact supports use a reversible coordinate

```text
n=1+8c+4u+2Phi+eta,                  1 <= n <= 64,
```

where `c` retains the coupled coefficient class. The odd-square label is an
equivalent display coordinate because

```text
((2n-1)^2-1)/8 = n(n-1)/2 = T(n-1).                    (2)
```

Equation `(2)` is useful because the large arithmetic-looking label retains
the small ordinal exactly. It is not a selection theorem: the inverse record
still needs the coefficient class and presence bits. THM-4355's proposed
address concatenates its `Theta!=0` and `Theta=0` strata into the same
`1..64` interval and is being audited with this inverse requirement.

For LRC the analogous design goal is not a single clever number, but a
Gödel-style reversible address whose fields are

```text
(row, owner, side, occurrence, clock, phase sheet, arrival component).
```

Only after proving an inverse on the exact finite lane should one replace the
tuple by its natural-number rank. Odd squares can be used as display labels;
the ordinal `n`, not the square, is the computational coordinate.

## Cheapest decisive LRC probe

The current frontier identifies owner/arrival as missing on the 165-row and
anchored `2+12` lanes. A bounded next experiment is:

1. choose one already finite exact lane, not the arbitrary-height problem;
2. construct the bipartite multigraph between normalized safe components and
   located endpoint occurrences;
3. retain owner, side, clock and phase sheet on every attachment;
4. partition configurations by the old support/owner quotient;
5. within each fibre, test whether the attachment-component map changes the
   next blocker, first exit, or signed current;
6. compare the cycle rank from `(1)` with the simple-graph/tournament rank.

A positive signal is a quotient fibre in which the old state is identical
but the exact attachment map changes a physical next operation. That would
identify the missing sidecar and justify enlarging the state. A negative
signal—attachment maps vary but every physical consumer agrees—would justify
compression on that lane. Merely finding more cycles is not progress unless
one of those consumers changes.

### Exact outcome

The test produced the pair `P=1287,9009` in the `h=420` safe-control bank.
Their completed owner/status quotient is identical, but 209 partial component
itineraries differ. On component 9, the next path changes from
`C3->P->D5->C2` to `C3->D4`, changing the local exit and signed-current
consumer. Every changed trace is nevertheless still missing, and the global
first missing component is unchanged. Thus the probe proves a precise
nonfactorization, not an LRC advance: partial attachment addresses are needed
by continuation/current consumers, while completed-cover status alone does
not need them on this pair. Full data and two exact implementations are in
[the finite-exact sidecar](lrc14-partial-attachment-address-nonfactorization-codex-20260902.md).

## Procedurally generated tasks

- **Anchor:** THM-4355 and the literal `U=K=0` reconstruction are complete;
  pull the closed endpoint atlas back through the source-normal row-eight
  response before spending work on another abstract wall.
- **Niche:** package `(1)` as a small attachment-rank lemma only if a second
  non-Jacobian consumer genuinely uses it.
- **Wildcard:** after the hostile pair above, instantiate the `h=420,u=3`
  collar deletion and test the remaining 828 components with an actual
  counterexample-conditioned ten-tail family.
- **Hostile:** force two marked endpoints to different components while
  keeping their unordered owner labels fixed. Any proposed `+1` cycle rule
  must reject this configuration.
- **Compression test:** verify an exact inverse before replacing a structured
  tuple by `n` or `(2n-1)^2`; otherwise the rank is only decoration.

The useful cross-problem principle is narrow: preserve the landing address
until after the consumer has been evaluated. No Jacobian-to-LRC implication
is claimed.
