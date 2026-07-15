# LRC(14) continued-fraction frontier — what the old threads actually preserve

This is a historical and forward-looking map after THM-778.  Continued
fractions have repeatedly found real structure in the LRC workspace, but the
phrase “use continued fractions” has referred to several mathematically
different operations.  Keeping those operations separate is the main lesson.

## The accumulated record

| thread | exact contribution | what it did not preserve / prove |
|---|---|---|
| THM-536, Sturmian seven-sector walk | Rewrites the consecutive cluster as partial sums of a mechanical word modulo seven; proves the `k<=6` vanishing and subset domination. | The joint seven-sector cover is aggregate. Per-block extremality, spread monotonicity, and translation invariance all fail. |
| HYP-3732/3738/4078, Ostrowski ladder | Identifies the exact restricted-family spine `[0;n-1,k]=k/((n-1)k+1)` and its three-gap geometry. | The court case `CASE-convergent-not-covering-min` refutes the early universal covering-min reading. A beautiful rung is not automatically the global extremum. |
| HYP-3718 and clean-base reflections | Explains semiconvergents, distance-one landers, and the rotation-orbit tournament at a chosen modulus. | Several early global-optimality sentences are historical overclaims. The base and orbit still need owner, metric, and realizability data. |
| THM-622, Farey-cell void | Turns a spectral gap into an exact numerator/depth restriction inside a Farey cell. | It is a reduction to a finite covering obstruction, not the obstruction itself. |
| THM-637, Farey roof | Proves the largest-gap function for `AP_k` is linear on each Farey cell and gives exact tail/mean ledgers. | The formal proof needs only Farey neighbours, not a continued-fraction slogan; extrapolation from AP cells to arbitrary families remains separate. |
| THM-565, apex-ruler sampling | Gives the exact `Vc +/- arcCount` lattice-sampling error and connects a positive slow good set to a finite witness. | The proof is interval counting. Calling it “three-gap” does not supply the missing positive-measure floor. |
| THM-736, far-peel Farey arcs | Computes the deep-well safe set as twelve exact Farey arcs with measure `H_12/91`. | It is a closed-form instance, subsumed by broader peel results; it does not classify arbitrary cores. |
| HYP-3791, convergent resonance lattice | Detects the persistent `13Z` resonance skeleton around `[0;13,14]` and the empirical redundancy sign. | The multi-far uniform bound remains heuristic; a convergent names where correlation lives, not its globally useful sign. |
| HYP-2772, Stern–Brocot tree sum | Proves an important negative result: the unsigned per-ratio discrepancy sum diverges level by level. | Absolute summation over the Farey atlas is unavailable; a signed/covolume mechanism is required. |
| THM-745, Euclidean residual tower | Isolates an exact sawtooth/Euclidean residual identity. The tempting CF fluctuation then vanishes from exposed segments by mirror pairing and the no-wrap lemma. | A visible Euclidean tower may cancel before reaching the proof functional. One must identify the surviving signed observable first. |
| HYP-3114, irrational approximation sidecar | Gives the correct interval-margin transfer: a convergent is useful only after a positive witness interval and robust radius are known. | Named constants or unusually good approximants cannot create the interval they are supposed to hit. |
| THM-780, phase-pigeonhole measure bridge | A strictly deeper witness gives an explicit uniform positive-measure safe set by an anchored joint-phase return; no height or continued fraction is needed. | Positive measure is metric substrate, not a certificate that a prescribed convergent denominator lands in the required component.  Component geometry and targeting remain separate. |
| THM-722/HYP-6280, leader ledger | Interprets the deep-well Ostrowski rung as a stopping time of a metric handoff staircase. | The average conservation law is nearly a factor two too weak; the cutting/lander pattern, not its total mass, is the missing information. |
| THM-773, prime token polynomial | Gives the exact finite sheet fibre and proves endpoint order/inverse steps/carry are missing from the metagraph node. | It deliberately left the endpoint schedule as an undecomposed field. |
| THM-778, centered endpoint skew product | Supplies that missing field: pairwise midpoint clocks are centered mechanical words; centered Beatty ranks reconstruct all simultaneous walls; the Euclidean parity cocycle recursively addresses the word; the word drives the `F_7` token fibre. | This is an exact language and reduction, not a proof that every persistent cover tears. |
| THM-779/783/784/788, r=8 token-walk chain | Characterizes blocking exactly, proves zero-sum visitor laws, refutes every absolute raw-wall cap by fast refinement, and contracts runs to decorated active-period normal form `E_0,V_1,...,V_A,E_A`. Joined to THM-778, this is a metric schedule-versus-fibre problem. | Fixed-index de-phasing is false (MISTAKE-148). The open theorem is a varying-index bound on active visitor packets plus core-component incidence, not a wall-count bound. |

## The conceptual correction

The useful continued-fraction object is not a scalar approximation quality. It
is one of four typed carriers:

```text
value address       [0;13,k] or a Farey parent cell,
metric roof         exact interval endpoints and widths inside that cell,
symbolic schedule   a mechanical/Christoffel cutting word,
recursive transport Euclidean substitutions acting on a labelled fibre.
```

Confusing these carriers caused most of the old overclaims.  A value address
does not imply global optimality.  A roof does not identify owners.  A word
does not determine the metric interval unless its scale and intercept survive.
A recursive address does not determine the LRC future unless the fibre action
is carried with it.

THM-778 now gives a precise four-layer object on the prime-seven lens:

```text
metric wall time x and centered rank R
    -> Euclidean/continued-fraction owner word
        -> labelled token translation in F_7^r
            -> polynomial cover observation and metagraph shadow.
```

The ordinary next-event tournament is always transitive.  Its isomorphism
class is therefore constant even while the eight-owner example traverses 948
labelled Hamiltonian paths.  Continued fractions belong to the labelled path
stalk, not to the iso node.

## Precisely what must be retained

For exact prime-sheet continuation, a continued-fraction packet must retain:

```text
owner labels and local event indices,
common scale gcd as well as the reduced ratio,
the centred half-lattice phase bit,
odd/odd simultaneous-wall blocks,
metric wall positions or an equivalent centered rank packet,
inverse steps w_a^(-1) mod 7,
current owner-to-sheet assignment,
endpoint deletion convention,
global carry,
the core-safe base interval.
```

Discarding the gcd loses repetition count.  Discarding the phase bit loses the
midpoint coset under Euclidean shears.  Discarding tie blocks changes strict
wall coverage.  Discarding inverse steps makes the same owner word act by a
different fibre transport.  Discarding the metric base repeats the winding
atlas mistake: a correct symbolic winding may miss every relevant component.

## New frontier avenues

### CF-01 — Euclidean block transfer on the redundancy polynomial

For `r>=8`, write `D_x=(X^7-X)Q_x` in a covered chamber.  A continued-fraction
partial quotient creates a long run of one owner inside a pair word.  Compute
the exact prefix action of that run on `Q_x`, not merely its net translation
modulo seven.  The desired lemma says that a sufficiently long Euclidean block
either exposes a noncovered prefix, hits a simultaneous wall, or returns with
a strictly simpler redundancy factor.

### CF-02 — The ten-wall absent-owner word

For the HYP-6835 survivor, THM-778 reduces 1,205 walls to ten covered ones with
palindromic owner word

```text
162,108,108,206,197,197,206,108,108,162.
```

The masks are now attached:

```text
25773,32153,31115,14635,615,
30093,31115,615,14233,6035.
```

They are not palindromic even though the owner word is.  The all-5,040 audit
explains this precisely: reflection does not descend to the 25-mask quotient
without the owner-labelled lift.  Next compute the Euclidean block between
successive lifted return masks.  Every pair has at least three intervening wall
blocks, so THM-779 identifies all ten as isolated one-wall blocking runs, not a
single persistent run.  The adjacent redundancy root already has the exact
period-five word `((1->0),(4->6),(2->4),(0->2),(6->5))^2`; ask whether the
between-wall block lengths
`57,301,3,24,329,24,3,301,57` reduce the base to five types and their mirror.
Ask whether that five-block first-return half classifies the mod-seven
degeneracy packets that could sustain a genuinely consecutive survivor under
scale/Farey mutation.

### CF-03 — Mod-seven S-adic automaton

Christoffel words are built by substitutions directed by partial quotients.
Let those substitutions act on the finite token fibre, with owner `a` acting
by `-w_a^{-1}`.  Minimize the resulting automaton by continuation equivalence.
The state must expose exactly which continued-fraction digits matter modulo
seven and which require a prefix-cover sidecar despite having the same net
translation.  THM-779 supplies the acceptance condition: `SURJ` on pieces and
the collision-pair hop constraint at walls.

### CF-04 — Centered Farey concatenation law

Standard Christoffel words concatenate under Farey mediants.  Derive the exact
centered version, including the parity cocycle and tie-block correction.  This
would make a Farey mutation of speed ratios into an explicit edit of the LRC
wall word, rather than an analogy about nearby fractions.

### CF-05 — Attach substitutions to metagraph edges

Use the merged metagraph node only as a base address.  Label its fibre edges by
the Euclidean substitutions that carry one owner-labelled sheet tiling to the
next.  Compare this substitution order with the blue/black and spine/rib/sea
address orders.  The goal is an objective index of a *transported* tiling, not
another static class rank.

### CF-06 — Minimal pairword observation graph

All pairwise centered words reconstruct the global wall order but are highly
redundant.  Find the smallest graph of owner pairs whose word packets still
reconstruct every wall/tie block for a named family.  This is a comparison-
graph/Menger problem with endpoint events as vertices; failure pairs would
identify genuinely higher-order scheduling information.

### CF-07 — Signed Farey sums after the divergence warning

HYP-2772 rules out unsigned Stern–Brocot summation.  Revisit the tree with the
actual signs supplied by token translations, mirror rank, or the THM-745
telescoping observable.  Prove convergence only after restricting to the
relation lattice realized by one speed family.

### CF-08 — Leader-ledger landers in centered ranks

Express every lander that cuts a sum-ruler climb as a centered Beatty-rank
constraint.  The hoped-for statement is that cutting all climbs below
`14/183` forces the AP endpoint lattice and hence violates covering.  This is
the most direct bridge from the old Ostrowski stopping-time insight to the
present endpoint word.

### CF-09 — Pairwise CF versus multidimensional CF

Do not import Jacobi–Perron or another multidimensional continued fraction
until pairwise ranks are shown insufficient.  THM-778 proves they are already
lossless for chronology.  A higher-dimensional algorithm is justified only if
it compresses the core-safe/polynomial continuation predicate, not merely the
speed vector.

## The closest concrete next theorem

The most focused next step is CF-02/03:

> Use the period-five redundancy-root word and the five exact interval-block
> types `(57,301,3,24,329)`, compile their Euclidean substitutions, and test
> whether `(Euclidean remainder state, phase bit, absent owner, hop target,
> mask, redundancy root, carry)` is continuation-complete on the exact
> first-return movie.

This would turn the current historical synthesis into a finite automaton whose
fields each have a mathematical reason to survive, and aim the corrected
THM-788 exit target at active zero-sum visitor packets rather than an empirical
or raw wall count.
