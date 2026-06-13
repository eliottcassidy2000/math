---
id: HYP-2101
status: SUPPORTED by bounded S579 apex-lift sheaf audit; unbounded gluing proof open
source: codex-2026-06-03-S579
related:
  - HYP-2102
  - HYP-2100
  - HYP-2099
  - HYP-2098
  - HYP-2097
  - HYP-2096
  - HYP-2095
  - HYP-2063
  - HYP-2045
  - HYP-2088
  - HYP-2083
  - THM-397
---

# HYP-2101: n=14 should be proved by gluing cheap-pair certificates across the mod-7 apex lift

## Claim

The next useful n=14 proof object is a certificate sheaf over the tight
mod-7 tie-wall, not another class table.

The base is the 14-gon boundary seen at the tight time:

```text
six nonzero mod-7 antipodal lanes:  +/-1, ..., +/-6
one self-antipodal apex lane:        0
```

A local section is an HYP-2095 certificate germ:

```text
small pair (a,b), reduced sum <= 14,
pinch or n-clock time,
paired-shield status,
endpoint-anchor status,
D/U/N obligation owners,
unit-spine and slack-owner labels.
```

The conjectural theorem is:

> Every measure-zero n=14 realization on the tie-wall has a global cheap-pair
> section after reducing and lifting through the mod-7 apex seam. If the
> section cannot be glued, the failed gluing forces a positive-measure interval.

Equivalently:

```text
measure-zero -> unblocked small pair,
block all local sections -> positive measure.
```

This is HYP-2095 localized to the HYP-2098/HYP-2099 n=14 boundary.

## Why This Extends The Current Picture

HYP-2100 says the unit-spine exchange lemma should be sieved by cheap pairs
first: in its bounded one-unit-lift audit, all `13169` full D/U/N covers had
an unblocked small-pair witness and the no-cheap residual was empty.
HYP-2102/Opus S571's lift-lemma checkpoint adds the other side of the apex seam:
bounded rows containing a multiple of `n` were all loose/positive-measure for
`n=4,6,8,10,12,14`, while a crude measure lower bound failed and left a
thin-arc-cover residual.
HYP-2098 says the 64 self-converse classes are a boundary overcount and the
tight speed family appears to collapse to `{AP,V*}`. HYP-2099 says the 64
classes are still valuable, but only as a surface on which owner labels must be
reattached. Opus S570 says every one of the `8191` gcd-1 transversals mod `27`
has a cheap unblocked pair, with AP the unique floor-tight transversal. S578
says the canonical unit-spine composite fibre has the same behavior: AP and
`V*` are the only floor rows through slack `42`, both using `(1,13)` at
`1/14`, while the block-all controls are positive measure.

The sheaf reading treats these as one statement:

```text
transversal flip lattice       = global sections already glue;
unit-spine exchange lifts      = cheap-pair sieve fires before exchange;
multiple-of-n apex rows        = positive-measure side of failed gluing;
canonical V* composite lift    = first nontransversal apex lift still glues;
block-all rows                 = failed sections become positive measure;
64 fixed classes               = boundary base, not certificate content.
```

So the proof target is not "prove 64 rows one by one." It is:

```text
prove local cheap sections survive every allowed apex/nonunit lift,
or prove the obstruction to survival is already positive measure.
```

## Bounded S579 Audit

`04-computation/lrc_n14_apex_lift_sheaf_s579.py` makes the sheaf finite over
the HYP-2100 unit-spine lift site.

```text
objects:       unit-shell lift subsets over a fixed four-slack row
restrictions:  lower one lifted unit-shell representative
sections:      ledger_failure, cheap_d14, cheap_side, positive_measure, residual
```

Here `cheap_d14` is the denominator-`14=2*7` apex chart.  It is central, but not
complete.

For the all-slack one-lift scan:

```text
rows=45045
full covers=13169
route_hist={ledger_failure:31876, cheap_d14:7943, cheap_side:5226}
restriction_pair_hist={cheap:13119, positive_measure:50}
restriction_residual=0
```

For named two-lift/local stress:

```text
AP:       full=240, cheap_d14=187, cheap_side=53, restriction_residual=0
V*:       full=336, cheap_d14=216, cheap_side=120, restriction_residual=0
open-gap: full=240, cheap_d14=168, cheap_side=72, restriction_residual=0
controls: zero full covers
```

The two-lift restriction fibres show why this must be a sheaf rather than a
single apex witness:

```text
AP:       {cheap,cheap}:150, {cheap,ledger}:67, {ledger,ledger}:1
V*:       {cheap,cheap}:250, {cheap,ledger}:57, {ledger,ledger}:1
open-gap: {cheap,cheap}:150, {cheap,ledger}:67, {ledger,ledger}:1
```

The `{ledger,ledger}` rows are union-only sections: each one-lift restriction
is not a full D/U/N cover, while the two-lift union is full and has a cheap
witness.  Therefore ledger-failure patches must be valid local sections.

The audit also corrects the apex language.  Witnesses involving the apex shell
`7` are rare; the apex is the denominator-14 chart, not usually the witness
runner itself.

## The Clocks That Matter

This framing keeps only five clocks.

```text
mod-7 tie-wall clock:
  the six antipodal lanes plus apex seam; decides boundary gluing.

mod-27 unit-visibility clock:
  the summand graph acted on by units; detects missed unit shells.

n-clock floor seam:
  the AP/V* wall witness at j/14; hosts the straddle cheap pair.

pair-sum pinch clock:
  the HYP-2095 witness extractor.

endpoint-anchor clock:
  the THM-397 dual ledger when cheap pairs fail.
```

Reset length, unlabelled round class identity, bare half-turn tournaments, and
raw Hamiltonian-path counts are diagnostics here. They matter only after they
carry one of the labels above.

## Candidate Theorems

### 1. Apex-Lift Preservation

Let a local cheap-pair section exist on the mod-7 tie-wall away from the apex.
Lifting a representative through a nonunit mod-27 shell either:

```text
preserves a cheap pair with the same reduced sum,
or creates a new pair-sum/n-clock witness,
or introduces an endpoint anchor that opens positive measure.
```

This is the missing bridge between the transversal census and `V*`-type
composite cousins.

### 2. Section-Obstruction Positivity

If every local small-pair section is shielded or anchored, then the shield and
anchor labels form a cover-circuit. For n=14, every such cover-circuit should
have a leaf, a private D/U/N pivot, or an endpoint owner that peels. Therefore
it cannot be measure-zero.

This imports HYP-2095 and S574's obligation-hypergraph logic.
Opus S571 warns that the proof should not rely on a crude global measure bound:
the residual is a thin-arc cover, so the positivity argument probably needs
endpoint-owner and private-pivot labels.

### 3. Tie-Wall Limit Functor

The two perturbation directions at a tight time may yield different round
classes, as `V*` already shows. The invariant should be the tie-wall limit
plus certificate labels, not either perturbation class alone.

This would repair the boundary-subtle containment from HYP-2098.

### 4. Strict-Off-Section Lemma

Any self-converse fixed boundary class that does not admit the AP/V* cheap
section should be strictly lonely. In the sheaf language, "not AP/V*" means
the local section either glues with positive slack or fails in a positive-
measure way.

## Summand Graph Connection

The base of the sheaf is the odd summand graph at `C=27`, reduced through the
mod-7 tie-wall:

```text
addition       -> antipodal shells {a,27-a};
multiplication -> unit inverses that expose missed shells;
nonunits       -> hidden lift layers;
apex           -> the self-antipodal residue-0 seam at the tight 14-gon.
```

The older summand-graph lesson was "addition makes shells, multiplication
makes them visible." The new extension is "local visibility certificates must
glue across the apex." Nonunit shells are not enemies by themselves; they are
places where a visible section may need to be lifted or transferred.

## Tournament Analysis

The tournament vertices should be certificate germs and gluing failures, not
runners and not naked classes.

```text
Vertices:
  local cheap-pair germs, endpoint-anchor germs, unit-spine owner germs,
  and apex/nonunit lift failures.

Pair observable:
  same pair family?,
  same reduced sum?,
  apex-owner conflict?,
  endpoint-anchor conflict?,
  private D/U/N pivot preserved?,
  positive-measure witness produced?

Switch:
  obstruction beats section;
  endpoint-anchor beats paired shield;
  preserved cheap section beats class identity;
  lower lift depth beats higher lift depth.

Tie Hamiltonian path:
  transversal sections, V* apex lift, other gcd-3 lifts, gcd-9 lift,
  endpoint-cover failures.
```

Expected fingerprints are mostly transitive ledgers. Nontrivial cycles should
mean competing certificate transports around the apex seam; those cycles are
the only "cohomology" this proof route needs.

## Assumption Challenge

This hypothesis rejects three assumptions:

```text
1. fixed round classes are the proof vertices;
2. transversals and V* are separate cases;
3. blocking all cheap pairs leaves a hard residual.
```

The quotient preserves the predicate HYP-2095 needs: existence or failure of an
unblocked small pair. It destroys exact round-class identity and most of the
Hamiltonian-path complexity, intentionally. The staircase computations are a
warning here: near-regular objects can have huge path counts, but the proof may
live in a smaller local incidence law.

## Next Concrete Test

S579 gives the first local lift-site audit.  A useful S580 computation would
attach the same section data to the HYP-2099 fixed-boundary classes:

```text
rows:
  AP, V*, transversal flip representatives, and minimal gcd-3/gcd-9 lifts

columns:
  mod-7 lane owners,
  mod-27 summand shell,
  best cheap pair,
  shield/anchor blockers,
  D/U/N private pivots,
  positive-measure or measure-zero status.
```

The desired output is not another large census. It is a list of restriction
maps: which cheap-pair certificates survive when a speed is moved through the
apex/nonunit lift, and exactly where a failed move creates positive measure.

## Honest Status

This is not a proof of n=14.  S579 now has bounded computation on the local
unit-shell lift site, but the unbounded fixed-class sheaf remains open.
The value is organizational plus computational: it names a bridge lemma that
absorbs HYP-2100's cheap-pair-first exchange sieve, HYP-2102's multiple-of-n
positive-measure signal, HYP-2098's tie-wall correction, HYP-2099's owner
scaffold, Opus S570's transversal certificates, HYP-2096's unit-spine normal
form, HYP-2095's paired-or-anchored split, and the summand-graph
addition/multiplication dictionary.

## Files

- `04-computation/lrc_n14_apex_lift_sheaf_s579.py`
- `05-knowledge/results/lrc_n14_apex_lift_sheaf_s579.out`
- `07-reflections/lrc-apex-lift-certificate-sheaf-s579.md`
