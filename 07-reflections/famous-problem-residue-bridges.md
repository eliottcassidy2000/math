# Famous Problems Through Residue Calculus

**Session:** codex-2026-05-30 S356
**Mode:** exploratory synthesis with one exact probe
**Web anchors sampled:** Lonely Runner survey and recent finite-runner work;
Gilmer/Yu on union-closed sets; Oliveira-Thatte on graph reconstruction;
Nathanson on Caccetta-Haggkvist; Elsholtz-Tao and related Erdos-Straus work.

This session asked a deliberately loose question:

```text
Can the repo's projection / quotient / residue language point at a famous
or less famous open problem in a way that creates usable mathematical work?
```

The answer is yes, but not as a claimed proof.  The productive thing is a
problem grammar:

```text
choose a natural quotient
find its occupied fibers
look for gaps, boundary residues, or surviving deletion residues
feed those residues back as the next quotient
```

That grammar appears in several web-searched problems, with Lonely Runner as
the cleanest immediate target.

## Feedback Loop

### Iteration 1: Start From The Repo

The recent residue calculus sessions gave four stable motifs.

1. **Good-cut interval gas:** a set of intervals covers a path; bucket `1` is
   impossible; THM-355 turns missing buckets into zero transport rows/columns.
2. **Deletion residue:** H=63 is an exact odd-cycle kill, while THM-025 is a
   near-kill whose tiny residue has enough rank to be dangerous.
3. **Flat versus localized residues:** Paley-like objects spread mass evenly;
   core-like objects concentrate mass around one projection.
4. **Base-42 Erdos-Straus:** finite residue classes are covered by identity
   families, and the uncovered-looking classes are the real problem.

The common noun is "residue", but the useful verb is "transport": how does
mass move or fail to move through a quotient?

### Iteration 2: Search For Problems With The Same Skeleton

The web pass separated the problems into quotient types.

| Problem | Natural quotient | Residue target |
|---|---|---|
| Lonely Runner | time circle -> speed torus / forbidden arcs | uncovered time residue |
| Union-Closed Sets | family -> coordinate frequencies / entropy shadow | high-frequency coordinate |
| Graph Reconstruction | graph -> vertex-deleted deck | deck residue that distinguishes graphs |
| Erdos-Straus | prime/residue class -> identity family | uncovered congruence class |
| Caccetta-Haggkvist | digraph -> short-cycle/outdegree quotient | forced directed cycle residue |

The strongest match is Lonely Runner, because the standard formulation already
is an interval-cover problem on a circle.  No metaphor tax.

### Iteration 3: Build A Probe

`04-computation/lonely_runner_residue_probe_s356.py` pulls each forbidden arc

```text
dist(v t, Z) < 1/(k+1)
```

back to `R/Z`, using exact rational arithmetic.  For integer speeds every
boundary lies in

```text
Q(V) = (k+1) * lcm(V).
```

The probe records:

```text
component_count,
forbidden_length,
positive_gap_count,
max_gap,
boundary_witness_count,
boundary_modulus,
max_gap_boundary_residues.
```

This makes the LRC witness a quotient-gap certificate.  Initial segments
`{1,...,k}` are tight: no positive gap, but boundary witnesses such as
`1/(k+1)`.  Mixed or random sets have positive gaps, sometimes tiny.  That
splits the search space into boundary-only and positive-gap regimes.

## Bridge 1: Lonely Runner

The reduced LRC says that for any `k` nonzero distinct integer speeds, there
is a time `t` such that every runner is at distance at least `1/(k+1)` from the
origin.

Residue calculus translation:

```text
forbidden fibers = pullbacks of open arcs around integers
lonely time      = complement residue
tight example    = boundary residue, not positive gap
generic example  = positive quotient gap
```

This resembles the good-cut gas almost exactly.  Good-cut buckets are interval
unions on a path; LRC forbidden sets are interval unions on a circle.  THM-355
says empty fibers are silent rows/columns in finite transport.  The LRC
analogue would say:

```text
if all positive gaps vanish, the surviving witness must be a boundary residue
visible in Q(V).
```

That is HYP-1794.

S357 sharpened this point: for finite open interval unions, "positive gap or
boundary residue" is automatic.  The genuine LRC obstruction is an open cover
with no surviving endpoint.  This becomes HYP-1802, the endpoint-protection
hypergraph in the finite boundary quotient `Q(V)`.

What to try next:

1. Exhaust primitive speed sets for small `max(V)`.
2. Rank by `max_gap/threshold`.
3. Study the boundary-only stratum separately from the positive-gap stratum.
4. Search for transport moves on speed sets that preserve or shrink the
   boundary quotient while keeping witness residues visible.

This feels like the best bet because the computation is exact and the problem
already wants finite quotient language.

## Bridge 2: Union-Closed Sets

Union-closed sets ask for an element appearing in at least half the sets of a
nonempty union-closed family.  Gilmer's breakthrough replaced the target `1/2`
with the first constant lower bound through entropy growth under union; Yu and
others sharpened the information-theoretic optimization.

Residue calculus translation:

```text
upstairs object  = distribution on sets
projection       = coordinate-frequency vector
operation        = union transport A,B -> A union B
residue          = entropy increase not seen by low coordinate frequencies
desired gap      = a coordinate fiber with mass >= 1/2
```

Repo connection: HYP-1781 says the residue after a projection can be more
informative than the projected invariant.  Union-closed families may be a
clean non-tournament lab for that: if every coordinate frequency is below
`1/2`, union transport should leak entropy in a direction that cannot remain
hidden forever.

Speculative attack:

```text
Define deletion residues F_i^0={A in F: i notin A} and F_i^1={A in F: i in A}.
If every |F_i^1| < |F|/2, then every deletion shadow is larger than its
occupied fiber.  Try to show union closure forces a low-rank residue collision
analogous to THM-025: small survivor, high rank.
```

The repo's `rank_res` idea could become a family-rank statistic:

```text
rank_union_res(i) = max chain/antichain complexity surviving in F_i^0.
```

If every coordinate has low frequency, at least one zero-fiber should become
too large and too union-rich to avoid generating the missing high-frequency
fiber.

## Bridge 3: Graph Reconstruction

The reconstruction conjecture says a finite simple graph on at least three
vertices should be determined by its vertex-deleted deck.  Oliveira and Thatte
frame it algebraically: Kocay-style covering-number matrices give rank lower
bounds on the number of decks, and full rank would certify reconstruction.

Residue calculus translation:

```text
upstairs object  = graph G
projection       = all vertex deletions G-v
shadow           = deck
residue          = invariants surviving deletion comparisons
rank certificate = covering-number matrix rank
```

This is strikingly close to HYP-1785.  The repo already treats vertex deletion
as a projection and measures the max-loss residue of odd-cycle conflict
graphs.  Reconstruction asks for the universal version: can deletion residues
over all vertices distinguish the upstairs object?

Speculative attack:

```text
For each graph G, build a residue vector from all deleted cards:
  cycle-count residues,
  cut-rank residues,
  deck-level support multiplicities,
  small induced-subgraph covering ranks.
Then search for pairs with same deck-like low features but different high
residue, exactly the way Paley/Interval split after sharing shadows.
```

The Oliveira-Thatte rank formulation gives a serious target: residue vectors
are only useful if they generate rows of a covering-number matrix whose rank
approaches the number of isomorphism classes.

## Bridge 4: Erdos-Straus

The repo already has a base-42 thread for Erdos-Straus:

```text
4/n = 1/x + 1/y + 1/z.
```

Elsholtz-Tao emphasize the prime case and count solution averages; the local
repo notes classify many identities as residue covers.  The residue calculus
reading is:

```text
identity family = transport rule
covered class   = occupied fiber
hard class      = quotient gap
multi-r search  = boundary residue repair
```

The new feedback from LRC is to distinguish:

```text
positive gap     = no identity family covers a residue neighborhood
boundary witness = a limiting identity catches a class exactly
```

For Erdos-Straus, "boundary" is not geometric length but algebraic equality:
a parameter `r` may land exactly on a divisor condition.  The next useful
computation is a boundary-residue ledger:

```text
For p == 1 mod 12, record the first r such that p(p+r)/4 has a prime
factor 2 mod 3, then study first-r residues modulo 42, 84, 168, ...
```

That would turn the old multi-r covering observation into a quotient transport
table.

## Bridge 5: Caccetta-Haggkvist

Caccetta-Haggkvist says a digraph with minimum outdegree at least `n/k` should
contain a directed cycle of length at most `k`.  Nathanson's survey highlights
the additive-number-theory proof in Cayley and vertex-transitive cases.

Residue calculus translation:

```text
projection       = out-neighborhood growth / Cayley sumset quotient
forbidden state  = no short directed cycle
residue          = expansion mass that cannot return too late
certificate      = short-cycle boundary
```

This connects to the tournament side through cycle-packing and support
residue.  A no-short-cycle digraph is a transport system where mass keeps
escaping without closing.  Minimum outdegree says each row has enough outgoing
mass.  The conjecture asserts that the escape cannot avoid a short return.

Speculative repo move:

```text
Measure short-cycle residues under vertex deletion:
R_v^{<=k}(D) = short directed cycles avoiding v.
```

Then look for a Caccetta-Haggkvist analogue of HYP-1785:

```text
minimal counterexamples should have no exact kill and no small high-rank
near-kill; otherwise deletion or quotient transport exposes a short cycle.
```

This is weaker than a proof, but it gives a feature layer for searching small
extremal digraphs.

## What These Threads Really Represent

The repo has been repeatedly finding one structure in different clothes:

```text
An invariant becomes useful only after we ask what it forgot.
```

Score forgets cycle structure.  A deck forgets the deleted vertex labels.
Union frequencies forget the interaction under union.  LRC speed sets forget
the forbidden boundary quotient.  Erdos-Straus congruence classes forget which
identity family can repair them.  Caccetta-Haggkvist outdegree forgets return
time.

So the residue calculus thesis should be pushed one notch wider:

```text
A hard combinatorial problem often asks whether every object has a witness
outside a forbidden projection.  The next useful invariant is the residue
left on the boundary of that projection.
```

## Ranked Bets

1. **Lonely Runner quotient-gap certificate.** Best immediate target.  Exact
   rational computation, direct interval geometry, and a close match to
   good-cut gas.
2. **Graph Reconstruction deletion-rank dictionary.** Best structural match to
   HYP-1785 and existing deck/rank literature.
3. **Union-Closed entropy residue.** Best chance to connect repo language to a
   recent breakthrough method, but the entropy layer needs careful new math.
4. **Erdos-Straus quotient-cover ledger.** Strong local continuity with the
   base-42 notes; likely computationally rich.
5. **Caccetta-Haggkvist short-cycle transport.** Conceptually attractive, but
   the repo has fewer direct tools for general digraph expansion.

## Next Concrete Work

1. Extend `lonely_runner_residue_probe_s356.py` to enumerate primitive speed
   sets up to a bound and store the tightest positive gaps.
2. Add a `boundary_only` mode for tight LRC families.
3. Build a small "projection residue table" for union-closed toy families,
   using coordinate deletion fibers and entropy under union.
4. For graph reconstruction, create a tiny deck-residue matrix for graphs on
   `n <= 7` and compare its rank to deck classes.
5. For Erdos-Straus, convert the base-42 multi-r search into a first-r residue
   ledger.

The playful version: the same shadow keeps showing up on different walls.  The
serious version: quotient boundaries are becoming a reusable search object.
