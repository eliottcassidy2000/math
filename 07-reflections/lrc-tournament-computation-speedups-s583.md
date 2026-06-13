---
source: codex-2026-06-03-S583
status: roadmap after repo archaeology
tags: [LRC, Tournament-Analysis, speedups, A000568, automata, SCC, Cprime, HYP-2113]
---

# Tournament computations to leverage for LRC

The useful tournament speedups split into two types.

The first type is **certificate-safe**: it preserves an LRC predicate exactly or
preserves a proof obligation with enough labels to certify a branch.

The second type is **search-safe**: it cheaply ranks or compresses cells, but
cannot certify alone.  These are still valuable because LRC now has tiny hard
residuals; a good prefilter can decide which residual deserves exact CRT or
interval arithmetic.

## 1. Source-marked reachability is exact

THM-381 is the best tournament bridge in the repo:

```
observer lonely at t  <=>  observer is a source in the marked tournament.
```

This means source-marked target classes are real proof targets, not metaphors.
The target set size is `A000568(n-1)` after deleting the marked source.  For
algorithm design, this suggests:

- precompute source target classes;
- drive the arithmetic clock as a walk in marked class space;
- prove every primitive clock hits a source target or compactified boundary
  target.

Raw runner subtournament classes are unsafe: S511/S535 show they can mix good
and bad states.  The mark and the `1/n` observer threshold are the payload.

## 2. Use quotient stacks, not one quotient

S535 is the right warning and the right recipe.  Raw phase classes are too
coarse, but the stack

```
circular body -> source target -> gap/apex threshold fiber
-> kinetic flow fiber -> endpoint blocker fiber
```

is pure in bounded audits.  S539 strengthens the lesson: vertices do not need
to be runners.  Fixed sections, half-sections, section boundaries, voids, cover
arcs, residues, and proof obligations can all be better vertices.

Practical use: make colored quotient classes into certificate targets.  If a
cell lands in a good-only colored class, record the certificate and skip exact
interval search.  If it lands in a mixed or unknown class, keep its labels and
send it to the next gate.

## 3. `A^2` fingerprints can remove the canonicalization bottleneck

Sorted rows of `A^2` are verified complete for tournament isomorphism through
`n<=9`, with zero random collisions at `n=10`.  Even where not proved, they are
excellent cache keys with a canonical fallback.

This matters for LRC because quotient walks classify many nearby cells.  The
fast path should be:

```
emit colored/marked tournament -> compute A^2 fingerprint -> hash lookup
-> canonicalize only on collision or outside certified range.
```

For the n=14 fixed-boundary program, the 13-runner classes are larger than the
verified range, so the responsible use is deletion/fractal fingerprints or a
round/fixed subclass audit, not blind proof reliance.

## 4. SCC/good-cut machinery is endpoint protection in disguise

THM-354 says good-cut count is `n - #SCC(T)`.  THM-395 says backward-wedge debt
is exactly `3*c3`.  These are cheap, structural protection coordinates.

The LRC analogy is direct:

- tournament bad cuts are unprotected endpoints;
- backward arcs are protector intervals;
- SCCs are protection cores;
- private leaves are peelable certificate exits.

Use SCC and backward-wedge features on endpoint-owner graphs, pressure graphs,
certificate-germ graphs, and middle-state graphs.  A nontrivial SCC in one of
those proof graphs is the shape of a real residual.  A DAG should peel.

## 5. Keep middle states visible

S582/HYP-2109 adds the missing automaton layer.  A binary tournament edge should
often be the terminal projection of a pair automaton:

```
event word over L/M/R -> terminal state -> edge or tie path.
```

For LRC, `M` is the endpoint-cover wall, live residue, midpoint corridor, or
certificate gluing seam.  This is where HYP-2112, HYP-2108, and HYP-2107 act:

- all components terminal `M` around a cover circuit should violate the
  HYP-2112 gap functional or HYP-2108 midpoint positivity;
- live bounded residue `M` states should be killed by factored CRT/dominance.

S577 adds a useful downstream exit for this same middle layer.  A 3-term
summand relation `v_c=v_a+v_b` turns a fold into a shield on runner `c`'s
clock; if the additive circuit is absent, the sampled circuit-free rows sit
comfortably above `delta=1/(k+1)`.  So summand-graph circuit labels should be
kept until the certificate DAG can consume them.

Algorithmically: do not exact-check every pair.  Let stable `L/R` cells exit
early; spend arithmetic only on terminal or recurrent `M` cells.

## 6. Cheap TDA features are residual triage

The existing `tournament_tda.py` features are ideal prefilters:

- score sequence;
- `c3` and transitivity ratio;
- SCC defect;
- deletion-residue rank;
- Omega alpha packets;
- low Walsh/Krawtchouk coordinates.

These do not prove LRC by themselves.  Their job is to sort residuals:
transitive proof ledgers, SCC pressure cores, near-kills, exact residue kills,
and cyclic gluing failures should land in different bins before exact checking.

## 7. Near-transitive perturbation is the tie-wall mode

The repo's near-transitive speedup idea says that if a tournament is a base
order plus `k` upsets, local methods behave like `O(k*2^k)` rather than
Held-Karp `O(n^2*2^n)`.

In LRC, the analogue is not arbitrary ranking; it is perturbation around AP,
`V*`, and mod-7 tie walls.  The question becomes:

```
how many pair/certificate cells change when we perturb the wall?
```

If `k` stays small, the quotient update can be incremental.  If `k` explodes,
that itself is proof signal: the wall is unstable and should expose positive
measure or a cheap pair.

## Wild but concrete next build

Build a single stack runner:

```
lrc_marked_source_speedup_stack.py
```

For each speed row it should emit:

1. observer-source marked class;
2. threshold-colored section and boundary classes;
3. `A^2` fingerprint cache key;
4. SCC/good-cut and backward-wedge features on proof graphs;
5. `L/M/R` middle automaton summaries;
6. Cprime/certificate gate route with `Phi/P(S)` terminal calls and additive
   fold/shield exits;
7. residual middle circuits with endpoint-owner/residue labels.

The output should be a certificate ledger, not a yes/no verdict.  Every row
should end as one of:

- source target;
- known positive-measure gate;
- private pivot / cheap pair;
- bounded CRT contradiction;
- additive fold/shield exit;
- HYP-2112 gap-functional contradiction;
- HYP-2108 midpoint-positivity contradiction;
- small labelled residual for a human theorem.

That is the way to make tournament computation help the proof rather than just
decorate the search.
