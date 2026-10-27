# Independent adversarial audit — Session 75 black-core geodesics

**Verdict: PASS for the stated mathematical results and exact certificates.**
No mathematical correction is required. The regular theorem is reconstruction
of the user's inherited statement, not a novelty claim. The even extension
really requires the descending black-path input, not mere black activity.
The all-colour distance theorem and the black core theorem have different
targets; they must remain separate.

Audit scope: the frozen producer packet in `sources/producer/`, its 15 pinned
repository sources, a fresh independent exact implementation, and the
structural proof arguments below. Repository pin:
`96e35642110104cc547516a63f37eebcf930f707`.
No shared checkout, theorem, hypothesis, or user-facing document was edited.

## 1. Inheritance and adversarial posture

Closest proved mechanism: labelled endpoint/cut current under reversal of a
whole directed Hamilton path. Canonical hostiles: pure-blue C3 and R5,
`C3 -> C3`, `C3 -> 1 -> C3`, and a source–middle–sink endpoint. Corrected
near miss: averaging a cycle's score changes or merely finding a black path
cannot select a descending black path. Least-used sidecar: the successor map
on one retained Hamilton cycle of the deleted core, which produces two
*different named* starts.

Concept board: actual anti-automorphisms; endpoint current; directed cuts;
deleted-core strength; two named starts; two different target sets.
Anchor: regular exact black distance. Niche: even balanced extension.
Wildcard hostile: one-back-arc join of two pure-blue R5 blocks.

The operative palette is tested on every arc. Neither Hamilton-path counts,
closure of the path into a cycle, nor endpoint-score compatibility is used
as a sufficient condition for blue. Endpoint-score incompatibility is used
only as the valid sufficient obstruction to an anti-automorphism.

## 2. Source and correction-lineage audit

All 15 SHA-256 entries in the producer `SOURCE_LEDGER.json` match `git show`
at the stated pin. Their bytes are bundled at the same relative paths under
`sources/`. Eight frozen producer files are also bound in
`SOURCE_BINDINGS.json`; the portable checker verifies all 23 files first.

The load-bearing descending result is precisely:

- `05-knowledge/results/black_descent_20261020c/cycles.md`, Sections 2–6:
  every strong black-active tournament has a **black closed Hamilton path**
  with `score(last) >= score(first)`.
- `05-knowledge/results/blue_quotient_20261020b/structure.md`, Sections 2–3:
  black closed-path existence, blue-cut cosets, and the elementary even
  all-blue-cycle repair.
- `05-knowledge/results/pure_blue_atom_classification_20261020/UNIFORM_PREFIX_REROUTING_CLASSIFICATION.md`:
  repaired odd-circulant obstruction, using explicit failed arc equations.

I read the arguments rather than accepting their PROVED labels. I also
searched the pinned mistakes ledger for the exact descent/pure-blue/endpoint
claims. Its later cycle-changing repair explicitly supersedes the old
fixed-cycle obstruction. The historical affine-only test in THM-4730 is not
reused. THM-4730's current frontmatter distinguishes that INCOMPLETE historical
argument from its independent corrected completion; this distinction is
retained here. The filtered pure-blue exclusions in HYP-9348 remain intact.

The prefix proof can be used here without its older Lu-crossing reduction:
`structure.md` starts from any Hamilton cycle, handles a black cut directly,
and if all cuts are blue obtains the odd circulant by composing reversals;
its even case has its own elementary replacement. Thus this audit's actual
dependency chain does not silently require the historical crossing proof.
Strong-tournament Hamiltonicity is also accompanied by the elementary
longest-cycle insertion argument in the inherited source.

Startup: `agents/start_session.py` worked read-only. The processor lookup was
blocked because `.machine-id` is absent; none was invented. The producer
ledger's blanket statement that both commands require that file is at least
over-broad for this audit run; this is procedural metadata, not a proof gap.

A cached remote advance to `3613d304395` was observed during review. It adds
only honest RESERVED stubs THM-4811..4813 / HYP-9421, including an overlapping
arbitrary-strong terminal-distance topic. No reserved result was used and no
ID was claimed. This audit remains explicitly pinned, not silently rebased.

## 3. Structural argument: regular odd order

For a regular tournament on `n=2m+1`, summing outdegrees within a shore gives
exactly `c(A)=|A|(n-|A|)/2`. Consequently every cut begins at least `m`, and
every cut with both shores non-singleton begins at least `2m-1`.

The signed path-crossing sum telescopes to

`c_after(A)-c_before(A) = 1_A(last)-1_A(first)`.

Hence every state before move `m` is strong, independently of the chosen
paths. A source/sink endpoint needs at least `m` moves because one named
vertex's score changes by at most one per move. This is a genuine lower
bound to the requested core target, not only to a particular construction.

After the first black move ending at `v`, deleting `v` leaves a strong
core: before the move its minimum semidegree is at least `m-1`, every cut
has at least `m-1` arcs, and restricting the move to this deletion is again
a Hamilton-path reversal, changing each cut by at most one. Thus at least
`m-2>=1` remains. This exact inequality explains the order-seven threshold.

The remaining `m-1>=2` inneighbors of `v` have distinct successors on one
fixed core cycle. Rotating that *same* cycle and appending `v` supplies two
actual candidate starts, so one differs from the first start. At time one
its endpoint-score sum is `2m+1`, excluding blue. Thereafter two already
used distinct starts imply no vertex has been used as start more than
`j-1` times after `j>=2` moves; every later path ending at `v` likewise has
endpoint-score sum at least `2m+1`. Blackness is therefore justified anew
at every step; it is not inherited from an endpoint-closing slogan.

At move `m`, `v` is a source, every other vertex has score at least one,
and every nonsingleton cut is positive. If the deleted core were nonstrong,
its terminal SCC would yield a zero full-tournament cut: if its size were
at least two this is a forbidden nonsingleton cut, and if its size were one
it is a forbidden score-zero sink. **This is the necessary strong-core
argument.** Mere arrival at a singleton cut would not suffice.

The regular conclusion is therefore PROVED: exact black distance `m` to
`{1 -> S, S -> 1 : S strong on n-1 vertices}`, for `2m+1>=7`.
C3 and R5 are genuine pure-blue failure boundaries, not omitted small tests.

## 4. Structural argument: even balanced order

For `n=2m`, all scores `m-1` or `m`, both minimum semidegrees are `m-1`.
Writing `s` for the smaller shore,

`c(A) >= s(m-1)-s(s-1)/2`.

The concave quadratic is at least `m-1`; on `2<=s<=m` it is at least
`2m-3`, for `m>=4`. Thus preterminal strength and the terminal nonsingleton
cut margin are valid for `m-1` moves.

The descending black-path selector is essential. If its last score is
`m`, keep it; otherwise its two endpoint scores must both be `m-1`, so
conversing the tournament and reversing the word makes the last score
`m`. An arbitrary black path from score `m` to `m-1` would not justify
this step. The proof explicitly avoids that invalid substitution.

After deleting this chosen `v`, initial semidegrees are at least `m-2`;
smaller shores have size at most `m-1`, giving every cut at least `m-2`.
After the first move this is at least `m-3>=1`. There are `m-2>=2`
inneighbors of `v`, so the same-cycle successor construction again gives
the different second start. All subsequent endpoint sums are at least
`2m>n-1`. At the terminal, every other score is at least one; the regular
terminal-SCC proof applies verbatim. The lower bound is `delta=m-1`.

The proof is structural for `n>=8`, with independently complete exact heads
at `n=4,6`. Therefore the full even conclusion is supported, not an
extrapolation from sampled orders.

### Audit of the inherited even descending lemma

The source's Section 5 has the correct cases. Blue cuts on a Hamilton cycle
form either the empty set or a coset of the subgroup induced by even
rotations. At even order a nonzero size is `2d` with `d` an odd divisor of
`n/2`. More than two but fewer than all blue cuts give an odd-prime orbit
rotation; the block lemma applies with at least four orbits. Its additive
boundary condition follows from the explicitly proved Cauchy–Davenport
inequality. Equal orbit endpoint scores force blackness at even order.
The all-blue case uses the separately checked two-orbit construction.

In the remaining hypothetical obstruction there are exactly two opposite
blue cuts, and all others have strictly increasing scores. Starting at a
minimum-score vertex of score `delta`, blue reversal gives the terminal
maximum `Delta=n-1-delta`; each half-block increases strictly, so
`delta<=floor(m/2)`. The last eligible inneighbor of the start identifies a
terminal block of size

`2 <= t <= delta+1 < m`.

Thus that block lies wholly within the second increasing-score half.
Its initial vertex `x` has at least `delta` prefix inneighbors, whereas at
most `delta-1` of their distinct successors can beat the old terminal `u`:
`u` already has its predecessor as an inneighbor inside the terminal block.
This pigeonhole count provides both new boundary arcs for insertion.
The new last vertex `y` has `delta<=s(y)<Delta`. Consequently its path is
descending, while `s(y)+delta<n-1` excludes blue. The `n=4` separate
alternating-score argument is also valid. No arbitrary edge-to-cycle
extension theorem is assumed.

## 5. All-colour distance and the nonregular hostile

The all-colour distance from a strong tournament to the nonstrong locus is
exactly `lambda=min_A c(A)`. Lower bound: the eventual zero cut decreases
by at most one per move. Upper bound: retain a minimizing shore, remain
strong before `lambda` steps by the same lower bound, then cut a Hamilton
cycle at a complement-to-shore arc. The resulting word starts in the
shore and ends outside, decreasing that chosen cut by one. This proves
arrival in exactly `lambda` steps. No colour assertion is inserted.

The two-R5 one-back-arc source has `lambda=1`, `delta=2`. Independent full
first-step enumeration finds 657 Hamilton words, 624 black. None of the
657 reaches the strong-core target. A black nonstrong exit is

`(6,7,8,9,5,0,2,3,4,1)`.

Its actual failed anti equation is the image of `0 -> 1`, which would
require `6 -> 5`. The bundled producer two-black route is independently
replayed to a singleton plus strong core. Thus black nonstrong distance
is one and black core distance is two. This separates the target sets;
it is not evidence against either separately formulated conjecture.

The repeated-start seven-vertex hostile has colours exactly black, blue,
blue. Its final source is vertex 0, final sink vertex 6, and its middle
`{1,2,3,4,5}` is regular and strong (R5). It therefore lies outside the
core target despite having singleton extreme cuts. The two-start guard
protects both colour and the actual consequence object.

## 6. Independent computation, universe, and receipts

`audit.py` uses only the Python standard library and imports neither
producer program. Every labelled tournament is visited by a Gray walk
through all unordered-pair orientations; one arc and its two degree
counters change at a time. Balanced objects are filtered only by the
stated score range. Operations then use sets of ordered arcs rather than
the producer's bitwise move formulas. Anti-automorphisms are tested on
all actual directed arcs; strongness uses forward/backward reachability.

Fresh independently found routes give:

| universe | retained labelled objects | route length | actual R/K moves checked |
|---|---:|---:|---:|
| all balanced n4 | 24 | 1 | 24 |
| all balanced n6 | 2640 | 2 | 5280 |
| all regular n7 | 2640 | 3 | 7920 |
| total | 5304 | — | 13224 |

Every route verifies Hamilton directedness, actual nonblue status,
blackness after the move, exact endpoint score update, **every proper cut
current**, literal nonpath-complement K conversion, retained named labels,
strong intermediate states, distinct first two construction-gauge starts,
and final deleted-core strength. Structural-range constructors additionally
check the deleted core after the first move. The endpoint bound certifies
minimality. The full stream is SHA-256 committed per order in the receipt;
a literal witness is printed for each order.

An independent exhaustive descending-path test covers all strong labelled
sources through order six:

| n | strong sources | pure-blue exceptions | descending black closed witnesses |
|---|---:|---:|---:|
| 3 | 2 | 2 | 0 |
| 4 | 24 | 0 | 24 |
| 5 | 544 | 24 | 520 |
| 6 | 22320 | 0 | 22320 |

Every no-witness source is then checked against **all** of its Hamilton
words and found pure-blue. This is independent finite support for, not a
substitute for, the all-order proof review.

The all-colour minimizing-cut construction is independently run on every
strong labelled tournament in the same universes. At n6 the lambda counts
are 19680 with lambda one, 2640 with lambda two; the n5 counts are 520/24.
This does not claim a full all-order black-distance theorem.

All 210 saved producer moves, including larger balanced examples through
n13 and the n10 separation, are separately replayed and compared with
their literal saved K carriers. Pure-blue controls independently enumerate
all permutations: C3 has 3 Hamilton words, R5 has 15, and the two ordered
C3 hostiles each have 9, all blue. These samples are controls, not a large-
order census.

Reproduce from this directory:

```sh
python3 -B audit.py > normal.json
python3 -B -O audit.py > optimized.json
cmp normal.json optimized.json
```

The checker uses explicit exceptions rather than Python `assert`, and the
normal/optimized receipts are byte-identical. Each receipt binds the
checker and source manifest. `CHECKSUMS.sha256` binds the final packet.

## 7. Remaining scope and stopping point

This audit establishes the producer's regular and even-balanced results,
the all-colour nonstrong distance theorem, and both hostiles. It does not
prove arbitrary strong-source black distance equals minimum semidegree,
arbitrary strong-source black nonstrong distance equals lambda, all-order
black connectivity, arbitrary module lifting, or distances to forest or
transitive targets. Those frontiers remain separate and open in this packet.

Both R and literal labelled K routes are covered by the alternating-
converse identity; converse preserves the two targets, lambda, and delta.
The named-start invariant is recorded in the construction R gauge (or its
explicit reversed/converse gauge), not falsely asserted for the alternating
K word starts. Quotient statements concern invariant target sets and can
be lifted by carrying relabelings; no labels are erased inside a certificate.

Coherent milestone reached: proof audit, source authentication, independent
complete finite heads, normal/optimized replay, and explicit strong-core
hostile checks. No shared edits or publication were needed for this assigned
read-only audit. The parent retains responsibility for subsequent canon
promotion and collision-free publication.
