---
id: HYP-2907
status: GUARDRAIL / sharp proof obligation
source: codex-2026-06-22-S118
tags: [lrc14, digraphs, tournaments, h7, forbidden-values, realizability, open-q-108]
depends_on:
  - HYP-2905
  - HYP-2906
  - HYP-2881
  - HYP-2878
  - THM-200
  - THM-343
  - KPS-S31y
  - KPS-S31z
related:
  - HYP-2885
  - HYP-2908
  - HYP-2903
  - HYP-2904
  - HYP-+2880
results:
  - 04-computation/lrc_digraph_h7_guardrail_codex_s118.py
  - 05-knowledge/results/lrc_digraph_h7_guardrail_codex_s118.out
---

# HYP-2907: binary digraphs do not inherit the tournament `H=7` obstruction

## Claim

The prompt's two facts,

```text
arcs have two states
H=7 is forbidden for tournaments
```

are not enough by themselves to disprove LRC14 or to prove that a disproof is
impossible.  They become powerful only after one additional realizability
statement:

```text
an LRC14 apex over-cover induces a tournament OCF conflict graph,
not a generic binary digraph or incomplete oriented graph.
```

Equivalently, the missing theorem in the KPS-S31y route is:

```text
apex-7 over-cover  ->  Omega = K_3 in the tournament conflict-graph image.
```

Once that implication is exact, THM-200/THM-343 block the counterexample,
because `Omega=K_3` is exactly the forbidden `H=7` object.

## Exact audit

The script separates three alphabets that the slogan "two arc states" can
conflate.

### Tournament orientation bits

For tournaments, each unordered pair has exactly one of two orientations.
Exact enumeration gives:

```text
tournaments n=3: H=7 count=0, all H odd=True
tournaments n=4: H=7 count=0, all H odd=True
tournaments n=5: H=7 count=0, all H odd=True
tournaments n=6: H=7 count=0, all H odd=True
```

This is the finite audit shadow of THM-200/THM-343.

### Generic ordered-arc digraph bits

For ordinary directed graphs, each ordered arc is present or absent.  These
are also binary arc states, but they are not tournament orientation states.
They realize `H=7` already on four vertices:

```text
arcs=[(0,1),(0,2),(0,3),(1,2),(1,3),(2,0),(2,1),(3,0)]
H=7
unordered pair-state profile={'one_forward': 2, 'both': 3, 'none': 1}
```

So generic binary digraph structure cannot carry the forbidden-value proof.
The obstruction needs tournament completeness and antisymmetry, not just
bit-valued arcs.

### Incomplete oriented graphs

Even if bidirected pairs are forbidden, completeness still matters.  Oriented
graphs on five vertices with missing pairs realize `H=7`:

```text
oriented graphs on n=5 with H=7: 1440
first arcs=[(0,3),(1,2),(1,4),(2,3),(2,4),(3,1),(3,4),(4,0)]
profile={'none': 2, 'one_forward': 6, 'one_backward': 2}
```

Thus the LRC proof cannot merely construct an oriented conflict shadow.  It
has to construct the tournament OCF conflict shadow, or a labelled packet
equivalent strong enough to import THM-200.

### Tie-free winding tournament

The audit also corrects the common grid-sampling pitfall.  If winding
tournaments are sampled on wall times where a phase is `0` or `1/2`, the object
is not a tournament and `H` can become even or zero.  On tie-free cells for the
AP7 winding tournament,

```text
cells=24
H=7 cells=0
all cell H odd=True
cell H spectrum=[1,33,105,123,137,151,175]
H at x=1/7 is 175
```

So the winding audit supports the forbidden-value story only after the
tie-free/tie-resolved tournament condition is enforced.

## Consequence for LRC14

KPS-S31y is the right conceptual target, but HYP-2907 prevents an overclaim.
The finite atom theorem should be stated as:

```text
If a primitive bounded prime-covering 13-set over-covers at threshold 1/14,
then its apex-7 sector-pair conflict data realizes Omega=K_3 in the
tournament OCF image.
```

Then THM-200 says that object cannot exist.  This would prove the bounded atom
cannot be a counterexample; HYP-2906/HYP-2904 handle the large-speed branches.

Without that exact correspondence, the facts "binary" and "7 forbidden" do
not exclude a disproof, because binary non-tournament digraphs and incomplete
oriented graphs can have seven Hamiltonian paths.  The proof obligation is
therefore not a numerical identity `14=2*7`; it is a realizability theorem
placing the LRC over-cover inside the tournament conflict-graph functor.

## S48 order-2 packet integration

Incoming mac-mini S48 gives a sharper geometric candidate for the same
realizability theorem.  At the AP extremal `{0,...,13}` and `t=1/14`, the
winding comparison has exactly seven tied diameter pairs:

```text
(0,7),(1,8),(2,9),(3,10),(4,11),(5,12),(6,13).
```

This is the `14=2*7` packet as order-2 antipodal symmetry with seven orbits.
It is useful because it identifies the boundary packet that a counterexample
would have to beat.

It is also exactly where the guardrail matters.  A tied diameter pair is not a
tournament arc yet.  It is a wall-time degenerate comparison that must be
resolved into one of two orientations before THM-200 can apply.  Therefore
S48 should be read as a candidate source of the labelled OCF packet:

```text
resolved apex diameter packet -> tournament Omega=K_3.
```

The unproved statement "every resolution stays at `M>=1/14`" is the Node-2
atom bound in another language.  HYP-2907 says how to make it theorem-shaped:
prove that any attempted `M<1/14` resolution realizes the forbidden
`Omega=K_3`, rather than only an order-2 symmetric binary shadow.

## S31z logical-status integration

Incoming KPS-S31z clarifies the phrase "impossible to disprove."  LRC14 is a
`Pi^0_1` assertion over integer speed sets with decidable finite checks.  For
any sound arithmetic theory, a real counterexample would be a finite
certificate and hence refutable.  Therefore:

```text
impossible to disprove LRC14  <=>  LRC14 is true.
```

So HYP-2907 should not be read as an independence or indeterminacy route.  If
the apex-over-cover-to-`Omega=K_3` implication is proved, it is a proof of
LRC14's truth.  The guardrail is methodological: the implication must land in
the tournament OCF image, because generic binary digraphs can realize `H=7`.

## Tournament Analysis

Candidate vertex sets considered:

```text
tournament orientation bits
generic ordered-arc digraph bits
incomplete oriented graphs
abstract K_3 conflict graph
tie-free winding cells
wall-time winding digraphs
LRC apex sector-pair conflicts
labelled OCF packets
```

Pairwise observable: preserves the predicate "THM-200 forbids the realized
`H=7` object."

Ranking:

```text
labelled OCF packets
> tournament orientation bits
> tie-free winding cells
> LRC apex sector-pair conflicts
> abstract K_3 conflict graph
> incomplete oriented graphs
> generic ordered-arc digraph bits
> wall-time winding digraphs.
```

Challenged assumption: "two states per arc" means the same thing for LRC
sector shadows, ordinary digraphs, and tournaments.  The audit says no.  The
binary alphabet must be the tournament orientation alphabet, or the `H=7`
forbidden theorem does not apply.
