# HYP-2099: n=14 fixed-round certificates need a two-layer scaffold: 64 class boundary plus unit-spine owner labels

**Session:** codex-2026-06-03-S578  
**Status:** SUPPORTED by S578 scaffold/normal-form audit; bridge lemma open.

## Claim

The useful n=14 proof object is not a naked 64-row self-converse round-class table.  It is a two-layer certificate scaffold:

1. a fixed-boundary round scaffold preserving the 64 self-converse dihedral classes;
2. a labelled speed-level owner scaffold restoring HYP-2095 unblocked-pair data and HYP-2096 unit-spine/four-slack data.

The class layer is necessary because HYP-2094/HYP-2097 isolate the fixed boundary.  The owner layer is necessary because HYP-2095 certificates, THM-397 endpoint owners, and HYP-2096 slack debt are not determined by the unlabelled class.  Opus S570 adds the complementary transversal view: all `8191` gcd-1 antipodal transversals mod `27` have unblocked small pairs, and AP is the unique floor-tight transversal.  S578 adds the canonical unit-spine view, where V* reappears as the composite non-transversal floor cousin.

## Evidence From S578

`04-computation/lrc_n14_fixed_round_certificate_s578.py` builds the scaffold without expensive full A000568 canon by using the round generator's labelled `d`-vectors and dihedral anti-witnesses.

**Fixed-boundary scaffold.**

- labelled round `d`-vectors: `4096`
- labelled vectors with a dihedral anti-witness: `820`
- round classes: `316`
- fixed classes: `64`
- non-fixed controls: `252` classes, or `126` converse pairs
- score-span histogram over the 64 fixed classes: `{0:1, 2:20, 4:21, 6:12, 8:7, 10:2, 12:1}`
- SCC-count histogram: `{1:63, 13:1}`
- dihedral anti-count histogram: `{1:63, 13:1}`

So the fixed boundary is almost entirely strong and almost entirely has a unique dihedral anti-witness.  The regular class is the exceptional high-symmetry case.

**Speed-level owner scaffold.**

S578 then returns to the HYP-2096 canonical unit spine

`(1,2,4,5,7,8,10,11,13)`

and scans four slack runners among multiples of `3` through `42`, now adding the HYP-2095 unblocked-small-pair and positive-measure tests:

- slack rows: `1001`
- full D/U/N quotient covers: `531`
- rows below `1/14`: `0`
- floor rows: `2`
- open-gap rows: `1`
- full-cover rows with an unblocked small pair: `529/531`
- block-all rows positive-measure: `2/2`
- measure-zero full-cover rows with unblocked pair: `2/2`
- route histogram: `{clock-ledger failure:470, cheap unblocked pair:529, block-all but positive-measure:2}`

The two floor rows are exactly the AP slack `(3,6,9,12)` and `V*` slack `(3,6,9,24)`.  Both have the cheap unblocked pair `(1,13)` at `1/14`.  The two block-all full-cover controls are positive measure:

- slack `(3,9,12,42)`, `M=2/23`
- slack `(3,12,27,42)`, `M=3/32`

This matches HYP-2098's claim that the operative tight worry is the speed family `{AP,V*}`, not all 64 self-converse classes.

It also refines Opus S570.  The transversal flip-lattice handles all gcd-1 antipodal transversals and makes AP the unique tight point.  The unit-spine slack scaffold handles the first composite-27 non-transversal cousin: V* is floor-tight, still has the cheap pair `(1,13)`, and differs from AP only in the four-slack owner layer.

## Interpretation

HYP-2098 corrected the naive target: the 64 self-converse classes are a boundary overcount, and tightness seems to collapse to about four perturbation classes from `{AP,V*}`.  HYP-2099 keeps the 64 scaffold but demotes it to a routing surface.

The proof route becomes:

1. Work on the tie-wall / fixed-boundary limit, not a single generic round class.
2. For each candidate realization, test the HYP-2095 cheap route: find an unblocked small reduced-sum pair.
3. If every small pair is blocked, prove the row is positive measure by the paired-or-anchored shield/anchor ledger.
4. Use HYP-2096 to normalize forced unit-shell owners to a nine-runner spine, leaving only four composite slack runners to carry nonunit and endpoint debt.
5. Show all measure-zero rows in the normalized fibre reduce to the AP/V* cheap pair, or else cannot remain measure-zero.

For transversals, Opus S570 verifies this with `8191/8191` unblocked-pair certificates.  For the canonical composite slack fibre, S578 verifies the same dichotomy through slack `42`.

## Bridge Lemma

The remaining bridge is:

> Every n=14 tight-boundary realization either lifts to the canonical unit-spine owner scaffold without changing its HYP-2095 cheap-pair status, or failure of that lift already exposes an unblocked small pair / positive-measure interval.

Equivalently, the 64 fixed classes must be connected to labelled speed owners in a way that preserves the dichotomy:

`measure-zero => unblocked small pair`,  
`block-all => positive measure`.

S579/HYP-2101 gives the local language for this bridge: over each labelled
speed-owner fibre, use an apex-lift certificate sheaf whose objects are
unit-shell lift subsets and whose restriction maps lower one unit representative
at a time.  The denominator-14 apex chart is central but not complete; side
charts and ledger-failure local sections have to be part of the gluing data.
The bounded S579 audit found zero restriction residuals in the tested lift
site.

## Tournament Analysis

For the speed-level part of S578, vertices are canonical full-cover slack rows over the forced unit spine.  The pair observable is:

`(unblocked-pair flag, exact M, private D/U/N count, max speed)`.

The switch orients toward harder rows: no cheap pair first, then smaller exact `M`, fewer pivots, and smaller row.  The sampled 24-row fingerprint is transitive with no directed 3-cycles, so this is still a proof ledger.  Cycles should be sought only after arbitrary unit representatives, endpoint-owner fibres, or tie-wall perturbation directions are reattached.

## Assumption Challenge

This session explicitly rejects two tempting assumptions:

- The 64 self-converse classes are not themselves the proof certificates.
- A class-only table cannot decide HYP-2095's unblocked-pair route.

The quotient preserves fixed-boundary structure, score span, SCC shape, and dihedral anti-witness data.  It destroys small-pair blockers, endpoint owners, speed normal forms, and tie-wall perturbation direction.  The speed-level unit-spine scaffold restores those labels, but still destroys arbitrary large unit-shell representatives.  That is the exact next bridge.

## Honest Status

This is not a proof of n=14.  It is a certificate architecture plus a bounded normal-form verification.  The unbounded closure from HYP-2098 remains central: we still need to know whether tight speed families beyond `{AP,V*}` exist at large speeds, and whether arbitrary unit-shell representatives can be normalized without losing the cheap certificate.
