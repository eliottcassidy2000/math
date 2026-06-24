# LRC14 Apex-Aperture Comb Trichotomy - Codex S154

The main correction from this pass is that HYP-2961's `|Q|>=7` apex-multiple
family should not be treated as live in the ordinary proof frontier.  THM-571
already closes it modulo the accepted LRC<=13 input.  The branch left after
THM-568/THM-571 is `1<=|M14|<=6`.

After fetching concurrent work, this is HYP-2967 rather than HYP-2963,
HYP-2965, or HYP-2966: incoming HYP-2963 is the labelled-packet audit,
HYP-2964 is the Moon-core proof skeleton, HYP-2965 is the boundary-gap packet
bridge, and HYP-2966 is the NORK pinch-template atlas.  The aperture-comb
trichotomy is a local reducer inside that Moon core.

I wrote an exact local certifier.  At a denominator-14 unit apex `a/14`, each
non-multiple core speed is closed-safe.  If one side of the apex is not blocked
by an inward-moving boundary contact, there is an explicit rational core
aperture `(0,U)`.  If the danger combs of the `14`-multiples do not cover that
aperture, the midpoint of an uncovered gap is a strict LRC14 witness.

Full audit over the S151 AP-neighborhood banks:

```text
live_low_m14_rows = 18909
aperture-comb-certified = 12548
all-apertures-first-order-blocked = 4661
all-apertures-comb-saturated = 1700
```

The local lemma proves a large chunk of the remaining low-multiple rows with
explicit rational witnesses.

The guardrail is just as important: the named covering rows that S151 discharges
by exact Haar mass are first-order blocked at every denominator-14 apex.  So a
pure local-apex proof is false as a complete strategy.  Those rows need off-apex
Haar/source-kernel movement or state-lift labels.

The new proof shape is:

```text
qdiv/THM571 discharged
or aperture-comb certified
or full unit-support AP/GW-skeleton core
or comb-saturated tiny aperture.
```

Full unit-support means the core hits all six unit residues mod 14, so it is
AP/GW-skeleton-like.  Comb saturation gives rational inequalities between the
smallest multiple-of-14 comb and the core boundary speed; that feels finite or
scale-separated.

Assumption challenge: I considered runners, 14-multiple speeds, residual core
speeds, apex units, one-sided apertures, danger-comb teeth, packet families, and
proof gates.  The six-unit tournament is always the same transitive isomorphism
class, which means the scalar aperture observable is a reducer, not a complete
invariant.
