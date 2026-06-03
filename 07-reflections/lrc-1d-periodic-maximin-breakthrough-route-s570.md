---
source: codex-2026-06-03-S570
status: breakthrough-route synthesis
tags:
  - lonely-runner
  - maximin
  - one-dimensional
  - pinch
  - spectral-gap
  - endpoint-core
  - prime-gears
  - transversal
  - proof-strategy
---

# LRC as a 1D Periodic Maximin: Where This Repo Has Leverage

The important specialization is not merely "LRC is one-dimensional."  It is a
very rigid one-dimensional periodic maximin:

```text
M(S) = max_t min_i ||v_i t||
```

The objective is the lower envelope of finitely many periodic tents.  That means
its maxima are not mysterious high-dimensional events.  They are active-set
breakpoints.  In this repo, that fact has already been converted into several
nearly compatible proof tools:

```text
active-set exactness      -> pinch lemma / pair-sum moments
off-stratum clearance     -> minimax margin / spectral gap
tight-stratum structure   -> n-clock / CRT prime gears
boundary obstruction      -> endpoint protection core
recursive simplification  -> even-fold, transversals, multisieve
```

The possible breakthrough is to prove that all non-tight strata have a uniform
margin, then classify the remaining tight stratum by a tiny arithmetic clock.

## The Crucial Shape

In a generic high-dimensional optimization problem, the active set can be large
and unstable.  Here, the lower envelope is on a circle.  A non-apex maximum is
pinned by two tents, one rising and one falling.  S557/HYP-2059 proves the
consequent exact form:

```text
t* = m / (v_a + v_b)
M(S) = r / s, where s = (v_a + v_b) / gcd(v_a, v_b)
```

So every possible best time is a pair-pinch time.  This is the first real
resource: the continuum has collapsed to an arithmetic active-set catalogue.

For LRC at threshold `1/n`, a counterexample would need:

```text
M(S) < 1/n
```

The pinch lemma immediately says its optimal binding reduced pair-sum must have
`s > n`.  A tight family with `M(S)=1/n` must have `s` divisible by `n`, and the
floor case has `s=n`.

That is a scenario-specific gift.  It turns "find one good time" into:

```text
find a clearing pair-pinch, or prove the obstruction forces all optimal pairs
past the threshold denominator
```

## The Breakthrough Route

The route I would bet on is:

```text
Prove an off-stratum margin theorem.
Then classify the equality stratum.
Then show the endpoint-core obstruction cannot survive inside that stratum.
```

More explicitly:

```text
Theorem A, margin:
  If S is not in the tight arithmetic/resonance stratum, then
  M(S) >= 2/(2n - 1) > 1/n.

Theorem B, tight classification:
  If M(S) = 1/n, then every active optimal pair has reduced pair-sum n
  and the n-clock / CRT gear pattern is AP-like or one of the known finite
  sporadic tight forms.

Theorem C, no sub-tight core:
  If M(S) < 1/n, the endpoint protection cover has a nonempty labelled core.
  But such a core forces either a missed antipodal pair, a private endpoint,
  an even-fold descent, or a prime-gear clear setting.  Contradiction.
```

The repo already has partial ingredients for all three.

## Repo Resources That Matter

**1. S557 / HYP-2059: active-set exactness.**  This is the core theorem for this
specific maximin scenario.  The maximum is pinned by a straddling pair and the
value is `r/s` with `s` a reduced pair-sum.  This is much stronger than a
generic finite-checking claim.

**2. S562 / HYP-2075: pair-sum moduli are complete.**  The multi-sieve result
says pair-sum pinch clocks find exact witnesses in the tested pipeline.  The
"apex obstruction" disappears when the right active-set coordinates are used.
Fresh upstream stress tests by monad-compute extend this evidence: roughly
39,900 configs across `n=14..17` had zero PINCH misses.  The refinement matters:
pair-sum moduli are a complete witness family in the tested sense, but not
always the minimal-denominator witness family; small division moduli can clear
earlier.

**3. S552/S573: minimax margin data, corrected.**  Exact residue and small-box
enumerations showed `M1=1/n` and suggested a clean second value
`M2=2/(2n-1)`.  S573 corrects the global version: lifted integer rows can sit
strictly in `(1/n,2/(2n-1))`, for example `M=5/33` at `n=7`.  The margin is
still structural as a unit-shell edge, but the global theorem target is now a
clock-blocker ledger proving `M(S)>=1/n`, not a universal second-value gap.

**4. S553: antipodal transversals mod `2n-1`.**  This gives the most promising
proof of the margin.  Missing an antipodal pair mod `2n-1` immediately gives
the witness `M(S) >= 2/(2n-1)`.  The residual collapses to perfect
transversals, a `2^(n-1)` family whose AP member is the all-lower choice.

**5. S564 / HYP-2081: worry/ignore split and prime gears.**  Positive-measure
instances are ignorable.  The rare measure-zero residual is resonance-maximal,
and for it the n-clock `j/n` factors into CRT prime gears.  At `n=14`, the time
clock is the coupled parity gear and 7-gear.

**6. S556: the forced-multiple tension.**  A counterexample must contain a
multiple of `n`, since otherwise `t=1/n` clears.  But the searches suggest
multiples of `n` make configurations loose, not tight.  Turning this tension
into a quantitative lemma for sets with `n | v` is one of the cleanest
near-term targets.

**7. S357/S361/S430/HYP-1900: endpoint-core duality.**  If no primal witness is
found, the failure must be an all-protected endpoint cover.  This supplies the
dual certificate that a maximin proof usually lacks.

## A More Precise Breakthrough Candidate

Try to prove this statement for `n=14` first:

```text
If S is a primitive 13-speed set and M(S) < 2/27,
then the residues of S modulo 27 form a perfect antipodal transversal,
and the flip-set is one of the known tight flip-sets.
```

Then prove:

```text
No perfect antipodal transversal with a forbidden flip-set can realize
a nonempty labelled endpoint protection core below 1/14.
```

This uses the repo's exact tools in order:

```text
pinch lemma        -> any bad maximum has a pair denominator s
mod 27 witnesses   -> missing antipodal pair gives 2/27 clearance
transversal check  -> only flip-set residual remains
n-clock gears      -> tight residual must clear on j/14 unless a full gear aligns
endpoint core      -> any alleged sub-tight residual exports a finite dual object
pressure/peeling   -> eliminate the finite dual object
```

The reason this is better than broad finite checking is that every stage
shrinks by a theorem-shaped quotient:

```text
continuum time     -> pair-pinch moments
all speed sets     -> mod 27 antipodal transversals
all transversals   -> flip-set residuals
all residuals      -> endpoint protection cores
all cores          -> owner/pressure cycles or peels
```

## What To Compute Next

The companion script now started this direction:

```text
04-computation/lrc_witness_or_core_s570.py
05-knowledge/results/lrc_witness_or_core_s570.out
```

It returns a primal route (`n_clock`, `pair_sum`, `antipode`, boundary, or core)
plus endpoint-core data.  In the bounded run, pair-sum candidates recovered the
exact tropical maximum in every tested primitive box, n-clock handled the
resonance-maximal tight rows, and all endpoint cores peeled empty.

The next script should not merely search for witnesses.  It should classify
which theorem gate catches each instance:

```text
04-computation/lrc_maximin_breakthrough_audit_s570.py
```

For each primitive `n=14` set in a bounded range, record:

```text
M(S), t*, active straddling pair, reduced pair-sum s
whether M(S) >= 2/27
residue type mod 27: missed pair, zero residue, or perfect transversal
flip-set if transversal
whether S contains a multiple of 14
n-clock clear settings j/14
endpoint-cover status: positive gap, wall witness, empty peel, nonempty core
```

The output should answer:

```text
Are all below-2/27 sets forced into the known tight flip-sets?
Do all sets with a multiple of 14 have a uniform margin above 1/14?
When pair-pinch fails early, does endpoint peel always empty?
```

## Tournament Analysis / Assumption Challenge

The vertices should not be runners.  For this maximin scenario the right
vertices are active constraints:

```text
pinch pairs
antipodal residue pairs mod 2n-1
flip-set coordinates
n-clock gear settings
endpoint obligations
```

Suggested tournament:

- **Vertices:** candidate active pinch pairs `(a,b)` with reduced sum `s`.
- **Pairwise observable:** which pair gives the larger certified minimum over
  all runners at its best residue.
- **Switch/gauge:** orient toward the pair with larger max-min clearance; tie by
  smaller `s`, then by endpoint-core peel advantage.
- **Tie Hamiltonian path:** increasing reduced pair-sum, then lexicographic pair
  label.
- **Fingerprints:** score histogram, SCCs of mutually blocked pairs, directed
  cycles among small pairs, edge flips after adding endpoint-owner labels.

Predicate preserved: whether a small active pair can certify `M(S) >= 1/n`, and
whether failure exports a labelled core.

Information destroyed: continuous chamber order and exact wall timing, unless
endpoint-event labels are attached.

Challenged assumption: the runner set is not the proof object.  The proof object
is the active-set/endpoint-cover complex produced by a one-dimensional lower
envelope.

## Short Version

The repo's possible breakthrough is:

```text
Use the 1D maximin geometry to reduce time to pair pinches.
Use mod (2n-1) antipodal witnesses to prove a positive off-stratum margin.
Use CRT n-gears to classify equality.
Use endpoint cores to rule out sub-tight covers.
```

That fits the specific type of problem exactly.  It attacks the pathology that
a 1D periodic maximin can have: not lack of samples, but boundary equality with
too many resonances.
