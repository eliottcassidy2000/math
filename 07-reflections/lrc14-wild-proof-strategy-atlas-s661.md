# LRC14 Wild Proof Strategy Atlas S661

The user asked for a long, weird search for the missing LRC `n=14` insight.  I
went through the current core (`HYP-2164` through `HYP-2171`), the recent carry
and odd-wall work (`HYP-2230`, `HYP-2231`, `HYP-2235`, `HYP-2236`), and the
older stranger threads: certificate sheaves, pincer/rigidity, L/M/R automata,
two-block Helly, union-closed pressure, unit-distance impairment, CH-style
forcing jackknives, coimage/Yoneda, and anti-Poisson depth polynomials.

The thing that kept surviving every tangent is:

```text
the proof is no longer about finding the right scalar;
it is about proving that the retained side-channel is sufficient.
```

HYP-2164 already killed the least-positive `Res_27` quotient.  HYP-2167 says
why that quotient is not conservative: a speed is `v=r+27k`, and the `n=14`
clock sees `r-k`.  HYP-2230 then identifies the same carry `k` as the parity
and apex obstruction.  HYP-2231 names the active wall pairs.  HYP-2235 says the
finite-field analogue is pinned/concurrency data, not raw distance support.

So the missing theorem should be named:

```text
no-leak owner derivative theorem
```

or, in a more operational form:

```text
fixed odd wall
+ fixed C=27 gcd shell
+ carry/owner/deletion derivatives
  => AP, Vstar, 2AP, or strict looseness.
```

The S661 atlas script made this more than a mood.  It scanned `6629` research
files after excluding vendored build directories, ranked fifteen proof
strategies, and built a proof-route tournament.  The tournament was transitive:
no directed 3-cycles, singleton SCCs, one Hamiltonian path.  The majority
leaders were:

```text
carry_deletion_derivative
owner_concurrency_jackknife
apex_sheaf_gluing
two_block_helly_extractor
three_state_middle_automaton
pincer_grip_ledger
```

That is a good sign.  The weird methods did not create an incoherent swarm.
They organized themselves into a work queue.

## The Core Insight

The side-channel should be studied by **derivatives**, not only by sections.

S658/S660's deck work is the clue.  Tournament deck scalars collide, but
sorted loss profiles such as

```text
H-H_card, c3-c3_card, deleted score
```

separate the checked collisions.  The deleted vertex boundary is the missing
owner label.

LRC `n=14` has the same shape.  The least-positive section forgets carry.  The
way to recover carry may be to ask what changes under deletion or small carry
perturbation:

```text
Delta M,
active wall pair,
gcd-shell mass,
n-clock word,
pair-sum apex congruence,
owner route.
```

If a derivative deck reconstructs the carry cocycle on the floor neighborhood,
then the lift theorem becomes a reconstruction/no-leak lemma.  This feels like
the cleanest missing insight from the session.

Incoming S660 made this crisper after the S661 atlas was written.  Its
completed tournament deck derivative computation shows that through `n=6`,
global scalar repairs do not fix full-deck collisions, and even the unpaired
deleted-score multiset still leaks.  The clean repair is paired:

```text
(card isomorphism type, deleted vertex outdegree)
```

That is exactly the warning LRC needs.  We should not merely collect a multiset
of carry losses.  We need to attach each carry/owner derivative to the
`Res_27` card or certificate chart it modifies:

```text
projected proof card + boundary/carry derivative.
```

The pairing, not the raw lost-label inventory, is the candidate no-leak object.

## How The Tangents Fit

Kakeya/Falconer says scalar support saturates too early.  In `F_5^2`, every
one-line-per-direction row has full distance support, while pinned fibers and
concurrency still vary.  LRC translation: wall/pair-sum support can be fixed
while owner/carry packets still decide the proof.

Certificate sheaves say cheap pairs are local sections.  The next table should
not only ask whether a section exists; it should ask how a section restricts
after a carry move.

Pincer calculus says pair-sum pinches are only contact until the observer,
threshold, endpoint owner, and escape labels are attached.  A pincer proof must
say where failed pinches export their labels.

Rigidity/orbits say symmetry is not enough.  `n=14` has a rigid unit sheet, a
gcd boundary under doubling, and an apex chart.  The proof should show boundary
monodromy is owned.

The L/M/R automaton says middle states are where proof data hides.  A raw
tournament edge is too binary; the middle graph should be attacked directly.

Two-block Helly says some CRT residuals are not global at all.  They collapse
to singleton or pair determinant walls.  This should be the discharge after
the derivative/sheaf routes leave a multiple branch.

Union-closed pressure is a wilder but plausible repair: element frequencies
collide, but set pressure separates the first small collision.  LRC proof
obligations may need an analogous pressure on obligation packets.

The CH/forcing tangent gives the protocol name: build scalar twins, perturb
the generic side channel, audit absoluteness, and promote the labels that make
the predicate invariant.

## What To Do Next

The next computation should be narrow and ruthless:

```text
lrc14_carry_derivative_reconstruction_s663.py
```

It should start with AP, `Vstar`, and `2AP`, then generate one- and two-carry
derivatives, including minimal apex carries.  For each row, compute exact `M`,
wall-pair label, gcd-shell mass, n-clock word, pair-sum apex congruences, and
owner route.  Then ask:

```text
Can two rows have the same scalar/wall deck but different derivative deck?
Can two rows have the same derivative deck but different floor/strict status?
Do all non-scalar derivatives pay positive tax?
```

If the answer is clean in the local neighborhood, the theorem target becomes:

```text
local derivative sufficiency + global carry coherence
```

rather than "search more lifts."

The wildest reserve move I still like is tropical/min-plus: the active pair
walls form a slack cone, and AP/Vstar should be the degenerate zero-slack cells.
But that is a reserve move.  The main road is now much less mystical:

```text
differentiate the proof certificate and see whether the carry can still hide.
```

Artifacts:

- `04-computation/lrc14_wild_strategy_atlas_s661.py`
- `05-knowledge/results/lrc14_wild_strategy_atlas_s661.out`
- `05-knowledge/hypotheses/HYP-2237-lrc14-wild-proof-strategy-atlas.md`
