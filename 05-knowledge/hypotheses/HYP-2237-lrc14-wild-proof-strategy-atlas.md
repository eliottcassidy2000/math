---
id: HYP-2237
status: OPEN strategy theorem with S661 repo-search atlas
source: codex-2026-06-05-S661
related:
  - HYP-2236
  - HYP-2235
  - HYP-2231
  - HYP-2230
  - HYP-2222
  - HYP-2171
  - HYP-2167
  - HYP-2166
  - HYP-2165
  - HYP-2164
  - THM-401
  - THM-407
---

# HYP-2237: LRC14 Wild Proof Strategy Atlas

## Claim

The next serious LRC `n=14` proof attempt should not reopen a raw speed search.
The current quotient tower already says where the proof seam is:

```text
least-positive Res_27 quotient: classified by HYP-2164
owner/certificate bridge:       bounded-clean in HYP-2165
carry fiber:                    exact seam in HYP-2167/HYP-2230
odd wall/gcd carrier:           HYP-2231/HYP-2222
```

The most promising theorem shape is a **no-leak owner-derivative theorem**:

```text
fixed odd wall
+ fixed C=27 gcd shell
+ visible carry/owner/deletion derivatives
  => AP, Vstar, nonprimitive 2AP, or strict looseness.
```

The phrase "deletion derivative" is intentional.  S658/S660's tournament-deck
lane says scalar decks collide until one records loss profiles such as
`H-H_card`, `c3-c3_card`, deleted score, and boundary/influence data.  The LRC
analogue is to perturb or delete a runner/carry coordinate and record the exact
loss in the floor certificate:

```text
Delta M,
active odd wall pairs,
C=27 gcd-shell mass,
n-clock residue word,
pair-sum apex congruences,
owner/certificate route.
```

If those local derivatives reconstruct the carry cocycle or force strictness,
the global lift/CRT seam becomes a finite reconstruction theorem instead of a
large lift search.

Incoming S660's completed deck atlas makes this sharper.  Through tournament
`n=6`, the ordinary full deck and global scalar repairs still collide, and even
the unpaired deleted-score multiset still collides.  The clean repair is the
paired derivative:

```text
(card isomorphism type, deleted vertex outdegree).
```

That pairing resolves every checked collision.  The LRC analogue should
therefore keep the pair:

```text
Res_27/certificate card  +  the carry/owner derivative attached to that card.
```

## S661 Search Atlas

`04-computation/lrc14_wild_strategy_atlas_s661.py` scans the research repo
excluding vendored/build directories, ranks fifteen proof strategies, and builds
a proof-strategy tournament.

Search scope:

```text
scanned_files = 6629
dirs = 00-navigation,01-canon,04-computation,05-knowledge,07-reflections
```

The majority tournament is transitive:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1,10:1,11:1,12:1,13:1,14:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
hamiltonian_paths=1
```

The majority leaders are:

```text
carry_deletion_derivative      outscore 14
owner_concurrency_jackknife    outscore 13
apex_sheaf_gluing              outscore 12
two_block_helly_extractor      outscore 11
three_state_middle_automaton   outscore 10
pincer_grip_ledger             outscore 9
```

The weighted strategy score gives the same top two:

```text
1. carry_deletion_derivative
2. owner_concurrency_jackknife
3. apex_sheaf_gluing
4. orbit_boundary_rigidity
5. two_block_helly_extractor
6. three_state_middle_automaton
7. pincer_grip_ledger
```

This is useful because the wilder routes do not point away from the existing
proof tower.  They point back to a concrete derivative/reconstruction theorem
over the carry-owner seam.

## Strongest Routes

### 1. Carry Deletion Derivative

For each HYP-2164 floor or near-floor atom, delete one coordinate or toggle one
carry and record exact changes in `M`, odd wall pairs, gcd-shell mass, n-clock
residue word, and owner route.

Desired lemma:

```text
All non-scalar one/two-local derivatives are strict, and the only global
derivative patterns with zero tax are scalar AP/Vstar carry cocycles or
known owner/Cprime exits.
```

This is the most direct extension of HYP-2167 and HYP-2230.

The S660 reinforcement is that the derivative must be paired to the projected
card.  An unpaired multiset of lost labels can still leak.

### 2. Owner-Concurrency Jackknife

S659's finite-field Kakeya/Falconer model says direction coverage and distance
support saturate too early.  Pinned fibers and concurrency remain decisive.

LRC transfer:

```text
hold odd-wall and C=27 support fixed;
vary owner/carry packets;
measure which hidden packets decide floor vs strict.
```

This is the finite-field "pinned distance" idea translated into LRC owner data.

### 3. Apex Sheaf Gluing

S579 already says cheap-pair certificates are local sections over unit, nonunit,
and apex charts.  The S661 upgrade is to include carry derivatives as
restriction maps.

Desired theorem:

```text
global cheap section exists,
or failed gluing creates endpoint-owner positive measure.
```

This keeps the abstraction honest: the "sheaf" is a finite table of cheap
pairs, shields, anchors, D/U/N labels, carry moves, and restriction residuals.

### 4. Orbit Boundary Rigidity

At `n=14`, the unit witness sheet is rigid, doubling exits it into a gcd
boundary, and the apex chart is self-antipodal.  S590's orbit-functor language
suggests a compact target:

```text
no boundary monodromy class survives endpoint-owner peeling,
denominator shields, pincer closure, and CRT splitting.
```

This should be tested by tracking active wall-pair labels around the
`G=<2,-1>` shell fold.

### 5. Two-Block Helly Extractor

For rows that route into the multiple/Cprime branch, do not solve the full CRT.
Extract singleton or pair determinant Helly contradictions as in S599.

This route is not the whole carry theorem, but it is likely the right discharge
for residual rows after cheap-pair and local carry derivatives fail.

### 6. Three-State Middle Automaton

The tiny L/M/R automaton says an owner side cannot flip without passing through
a visible middle state.  A counterexample would need too many proof cells to
remain terminal-middle.

Target:

```text
compile owner/carry obligations into L/M/R states;
build the middle graph;
prove every closed-middle circuit peels by owner, CRT, or endpoint labels.
```

## Wacky Reserve Routes

Four routes remain worth keeping as controlled weirdness:

- **matroid fold circuits:** distinguish observer-coupled 3-term circuits from
  observer-blind 4-term energy;
- **rooted source perspectives:** quotient source/sink classes only after
  threshold and carry labels are retained;
- **Fourier major/minor shell split:** test whether `C=27` character packets
  separate gcd-3/gcd-9 carry behavior;
- **tropical min-plus polytope:** turn active pair-sum walls into a slack cone
  and search for a degeneracy-rank proof of strictness.

These are not first moves, but they are useful escape hatches if derivative
and owner-jackknife methods stall.

## Assumption Challenge

Candidate tournament vertices considered:

```text
runners, speeds, gaps, fixed circle sections, wall-crossing events, residues,
cover arcs, endpoint owners, carry cocycles, deletion derivatives, sheaf charts,
automaton states, proof obligations, and proof strategies.
```

S661 uses proof strategies as vertices because the session's predicate is:

```text
which proof route best preserves the data needed for LRC n=14 no-leak?
```

For future computations, the vertex set should become more concrete:

```text
carry derivatives,
owner/certificate germs,
L/M/R middle cells,
or determinant component languages.
```

Raw runners are the least useful default here; they forget the side channel the
current theorem actually needs.

## Handoff

The next script should be something like:

```text
04-computation/lrc14_carry_derivative_reconstruction_s662.py
```

Minimum contents:

1. start with AP, `Vstar`, and `2AP` over the HYP-2164 floor atoms;
2. generate one-coordinate and two-coordinate carry derivatives, including
   minimal apex carries from HYP-2230;
3. compute exact `M`, active odd wall pairs, gcd-shell mass, n-clock residues,
   pair-sum apex congruences, and owner route;
4. build a derivative deck for each row;
5. test whether identical scalar/wall decks but different derivative decks ever
   disagree on floor-vs-strict.

Success would be a small reconstruction lemma:

```text
the derivative deck is sufficient on the sampled carry neighborhood,
and every non-scalar derivative has positive tax.
```

That would be a real stride toward the lift/CRT conservativity theorem.

**See:** `04-computation/lrc14_wild_strategy_atlas_s661.py`;
`05-knowledge/results/lrc14_wild_strategy_atlas_s661.out`;
`07-reflections/lrc14-wild-proof-strategy-atlas-s661.md`; HYP-2236, HYP-2235,
HYP-2231, HYP-2230, HYP-2222, HYP-2171, HYP-2167, HYP-2166, HYP-2165, HYP-2164.
