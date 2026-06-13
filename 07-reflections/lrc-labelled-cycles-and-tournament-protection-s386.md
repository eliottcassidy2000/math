# LRC Labelled Cycles and Tournament Protection S386

**Session:** codex-2026-05-31-S386  
**Script:** `04-computation/lrc_tournament_labelled_cycle_bridge_s386.py`  
**Stored output:** `05-knowledge/results/lrc_tournament_labelled_cycle_bridge_s386.out`  
**Hypothesis:** HYP-1853

The bridge is sharper than a loose analogy:

```text
LRC endpoint protection is the circular arithmetic version of tournament
good-cut protection.
```

In a tournament, choose a Hamiltonian path

```text
v_0 -> v_1 -> ... -> v_{n-1}.
```

The cuts are the gaps between prefix and suffix.  A backward arc

```text
v_j -> v_i, i < j
```

protects every cut `k` with `i < k <= j`.  So every backward arc is literally
an interval of protected cuts on the path.

THM-354 then says:

```text
goodCuts = n - number_of_SCCs.
```

Read in the LRC dialect:

```text
good cut          -> protected endpoint
bad cut           -> unprotected endpoint
backward arc      -> protector interval
strong component  -> protection core
condensation      -> endpoint-debt / quotient-layer leak order
```

That is the same grammar as THM-365, except the LRC version lives on a circle
and its protectors carry arithmetic speed labels.

## The Important Difference

Tournament backward arcs are already valid labels: if the arc exists, the cut
interval it protects is automatically legal.

In LRC, a putative protection arrow needs the integer label

```text
|p*(n*m+eps)-a*n*u| < u.
```

That label is the whole difficulty.  S384 found abstract circular-arc
protection cycles immediately, but those are mirages unless the speed labels
realize them.

So the analogy says:

```text
Do not prove "no cycles."
Prove "no labelled cycles with the right arithmetic slack and measure balance."
```

## Tournament Parallels

Several repo threads line up with this:

1. **Endpoint transfer:** THM-356 says support matching is not enough; one
   needs private pivots, nonzero minors, or incidence rank.  LRC endpoint
   cycles likewise cannot be judged from overlap support alone.

2. **Endpoint collision hypergraphs:** HYP-1792 says collision triples are not
   parent-metagraph triangles.  LRC near-disproofs may look cyclic in the
   projected endpoint graph while having leaves in the labelled incidence
   hypergraph.

3. **Root signs:** tournament cycles are zero sums of type-A roots.  LRC
   labelled endpoint cycles are near-zero relations among characters
   `t -> v*t`, with strict slack as the boundary defect.

4. **OCF/Omega:** Omega stores which odd-cycle packets can coexist.  The LRC
   analogue should store which labelled protection arrows can coexist.  Its
   independence polynomial would not count Hamiltonian paths; it would count
   compatible endpoint-repair packets.

5. **Single-core anomalies:** H=63 single-core tournaments concentrate odd
   cycles through one vertex.  A true LRC counterexample would be a
   concentrated all-protected endpoint core.  Both ask for a finite language of
   core signatures with missing target values.

## New Picture

The labelled endpoint cycle is a boundary-state object:

```text
projection shadow:       endpoint overlap graph
incidence layer:         owner/protector labels and inequalities
residue:                 unprotected leaves, slack, quotient debt
core obstruction:        nonempty labelled cycle after peeling
```

This is exactly the repo's residue-phase-incidence split.  The LRC work should
now use the tournament toolkit:

```text
good-cut protection
endpoint-transfer private pivots
collision hypergraph peeling
root-sign packet modules
Omega-style compatibility graphs
```

## Questions

1. Can THM-354 be restated as a formal endpoint-protection theorem and used as
   a model for LRC?
2. What is the LRC condensation order replacing SCC condensation?
3. Can the first impossible labelled endpoint cycle at `n=14` be killed by a
   private endpoint, exactly like a private child column?
4. What is the LRC Omega graph whose vertices are labelled protection arrows?
5. Is the seven-ladder's tiny gap the endpoint analogue of a support matching
   with zero rank: convincing in projection, dead in incidence?

The slogan I would keep:

```text
LRC is tournament good-cut theory after the path has been wrapped into a
circle and every protector arc has been forced to pass an arithmetic label
test.
```
