# HYP-3100: Tournament-Contradiction Grammar

Status: SYNTHESIS / programmatic proof-frontier router; not a proof.

## Claim

The project's older proof-by-contradiction move

```text
constructed subproblem -> tournament shadow -> H in {7,21} -> impossible
```

should be treated as one member of a larger tournament-certificate grammar.
The grammar is useful only after the encoding map is explicit: a tournament
certificate can be pulled back to LRC14, H=21 closure, or a sieve sidecar only
when the relevant predicate is constant on the fibers, reconstructible from
retained sidecars, dual-annihilated, or routed to named residual debt.

## Programmatic evidence

`04-computation/tournament_contradiction_grammar_codex_s260.py` imports the
S31ah certificate engine, integrates the concurrent HYP-3099 applications and
the `TournamentH7.LRCBleedingEdgeFrontier` wrapper, and scores `24`
tournament proof moves against `12` current repo targets.  The stored output is
`05-knowledge/results/tournament_contradiction_grammar_codex_s260.out`.

The script verifies the seed certificates through the KPS engine:

- `H=7` and `H=21` are rejected by the permanent-gap table.
- even `H` values are rejected by parity.
- exact enumeration through `n<=5` has no even Hamiltonian-path counts and
  realizes `H in {1,3,5,9,11,13,15}`, with `7` already absent.
- sample Landau gates reject non-tournament score shadows.

The rebase over incoming S31ah/S65 adds three constraints that sharpen the
grammar:

- single-component H gaps are exactly `{7,21}` through the verified ladder,
  packaged as `K_3` and `K_10` nonrealizability in Omega;
- mac-mini HYP-3099 shows cap optimality is a non-transitive improvement
  tournament with bounded finite-check status, so tournament cycles can be
  diagnostic rather than fatal;
- the tempting `apex-7 = H=7` bridge is a coincidence: the genuine tournament
  lever is order-2 antipodal symmetry versus odd tournament automorphism
  groups, while baby-Hodge holes are c5/spectral and not forbidden-H events.

The generated target ranking is:

- coarse mod-14 winding: first repair by trienerment / fine-scale lift, not
  raw H.
- fine mod-7 winding: compare Paley/regular rigidity, cycle census, score,
  and skew spectrum.
- H=21 closure: combine H forbidden-value, SCC product descent, Omega-shape
  mining, and alpha-vector gates, with the `K_10` Omega case now the clean
  single-component target.
- level-7 state lift: H=7/H=21 can be terminal only after the output is a
  complete tournament certificate, not a loose directed graph.
- HYP-2963 packet normalizer: use automaton-state, matrix-observability,
  median-center, and bridge-protection tournaments to decide whether the H
  test is legal.
- sexy primes: use residue-sieve tournaments as admissibility and bookkeeping
  sidecars only; they do not break sieve parity or prove fixed-gap infinitude.

## Tournament Analysis

Vertex sets considered: runners, arcs, residues, Omega components, score
sequences, proof obligations, observer charts, sidecar columns, quotient maps,
finite sieve channels, automaton states, rooted perspectives, and branch-graph
edges.

Chosen vertex set for this pass: certificate functors and target proof
obligations.

Pairwise observable: number of target ledgers on which one technique scores
higher.

Switch/gauge: orient `A -> B` when `A` wins more target ledgers; break ties by
fewer destroyed coordinates.

Fingerprint of the selected eight-technique frontier:

```text
score histogram: {0: 1, 1: 1, 2: 1, 3: 1, 4: 1, 6: 3}
directed 3-cycles: 1
SCCs: [[7], [3], [1], [5], [2], [0, 4, 6]]
Hamiltonian path count: 3
path: Automaton-state tournament
  > H forbidden-value certificate
  > No-naked-bridge orientation
  > SCC product descent
  > Residue-sieve tournament
  > H-max rigidity
  > Trienerment lift for ties
  > Matrix observability tournament
```

The cycle among automaton-state routing, H forbidden-value testing, and
no-naked-bridge protection is the important signal: an H contradiction may need
packet-normalizer and bridge-protection checks before it is theorem currency.

## LRC14 implication

The live proof route should not ask first whether a raw residual has `H=7` or
`H=21`.  It should ask:

1. What finite-address packet or observer chart is being encoded?
2. What LRC predicate survives the tournament quotient?
3. Which sidecar is destroyed?
4. Does the destroyed coordinate have a legal discharge?
5. Only then, does the terminal complete tournament emit a forbidden H value,
   forbidden Omega shape, non-Landau score, invalid cycle census, or protected
   branch violation?

This integrates the S258/S259 observer-gluing frontier with S31ah's certificate
engine, HYP-3099's three concrete applications, the S64 symbolic-overlap
update, and the new `TournamentH7.LRCBleedingEdgeFrontier` formal wrapper.
THM-577 strengthens the cap side, and HYP-3099 explains why the remaining
cap-optimality problem can be a bounded but non-transitive finite check.  The
grammar still says the cap chart must glue to normalized arcs, level-7/CRT
descent, branch/K33 handoff, and formal witness discharge.

## Sexy-prime relation

The transfer to the sexy prime conjecture is methodological, not a proof.  A
residue-sieve tournament can track midpoint residue channels for the fixed gap
`6`, the local collapse modulo `3`, Hardy-Littlewood weight sidecars,
prime-power exceptions, parity debt, distribution debt, and terminal exits.
That is the same controlled-forgetting discipline used in LRC14, but it does
not supply a parity-breaking lower-bound sieve.

## Next pulls

- Extend the Omega-realizability miner to the remaining candidate
  `I(Omega,2)=21` components, starting from the `K_10` package and the
  known P4 exclusion.
- Feed THM-572/K33/F7 state-lift outputs through completeness, SCC product,
  and H-spectrum checks.
- Add tournament-certificate columns to the HYP-2963 packet normalizer and to
  `TournamentH7.LRCBleedingEdgeFrontier`; measure route/fiber separation before
  applying any H contradiction.
- Build a fine mod-7 winding scout that records score, cycle census, skew
  spectrum, Paley distance, and tie-lift status.
- For sexy primes, keep the tournament sidecar as a local-residue audit and
  state explicitly where parity and distribution remain external obligations.

## Links

Related: THM-200, THM-202, THM-343, THM-573, THM-577, HYP-2908, HYP-2909,
HYP-2963, HYP-3037, HYP-3047, HYP-3086, HYP-3092, HYP-3095, HYP-3096,
HYP-3097, HYP-3098, HYP-3099, T1178, LTI-239, LTT-137,
`TournamentH7.LRCBleedingEdgeFrontier`.
