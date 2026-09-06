# A recovered eleven-body realizes the joint overlap certificate

**Status: FINITE-EXACT physical geometry; PROVED family consequence relative
to the independently audited five-profile envelope; INDEPENDENTLY AUDITED.** The core is recovered
prior work, not a newly discovered speed set. The conclusion compares the
new joint sufficient certificate with its two scalar relaxations. It does
not claim the family was unsafe, previously unknown, or missed by every
other LRC tool. LRC(14) remains **OPEN**.

The [root audit](creative_20260906_root_audit.md) independently reconstructs
the entire safe set by a different method: 150 rational wall points and
149 intervening cells. All components, thresholds and the full-row phase
agree exactly.

## Inheritance and recovery

The new supplier is the independently audited
[creative_20260906_inert_pareto.md](creative_20260906_inert_pareto.md),
with [independent review](creative_20260906_inert_pareto_audit.md).
It preserves the joint pair of cross-comb mass and maximum literal component
width rather than optimizing those two coordinates separately.

Exact-set and exact-fraction searches recovered the eleven-body

    H={1,3,4,5,6,7,8,9,10,11,13}

in the earlier finite output
[lrc14_census_dichotomy_synthesis_kps.out](lrc14_census_dichotomy_synthesis_kps.out),
line 11, with approximate safe mass `0.04238` and clearance maximum `2/17`.
Other historical outputs also contain related `11/728` widths. Those outputs
supply provenance, not the proof below: the whole closed safe set is rebuilt
exactly in the standalone source. The contribution here is to reuse this
known physical body as a realization of the new joint certificate.

The corrected near miss is deleting isolated safe points when replacing an
ordinary measure calculation by a connected-component statement. The six
isolated points below are retained. The live concepts are the closed body
safe set, its mass and width, the coupled tail profile, and the physical
half-lift. Their map is exactly the doubled-body consumer of the envelope;
there is no asserted map from an arbitrary selected split to an actual
THM-3818 decoder partition.

## 1. Complete closed geometry

For `G_H={y mod1: ||hy||>=1/14 for every h in H}`, the complete connected
components, including zero-length ones, are:

| Closed component | Length |
|---|---:|
| `{1/14}` | `0` |
| `[15/182,13/154]` | `2/1001` |
| `{3/14}` | `0` |
| `{5/14}` | `0` |
| `[29/70,41/98]` | `1/245` |
| `[85/182,27/56]` | `11/728` |
| `[29/56,97/182]` | `11/728` |
| `[57/98,41/70]` | `1/245` |
| `{9/14}` | `0` |
| `{11/14}` | `0` |
| `[141/154,167/182]` | `2/1001` |
| `{13/14}` | `0` |

Thus

    M=mu(G_H)=5939/140140,
    L=max component length=11/728.                       (1)

These formulas follow by listing every open danger tooth
`((14k-1)/(14h),(14k+1)/(14h))` for `h in H`, clipping to `[0,1]`, and taking
the complement of their union. Two open teeth are merged only for a
positive overlap; teeth meeting at an excluded endpoint leave that isolated
safe point. Zero and one lie inside danger and create no circular joining.
The source gives this exact finite union certificate and independently
checks the clearance at every listed endpoint and midpoint and at every
intercomponent gap midpoint.

## 2. A physical gain for retaining both profile coordinates

Both scalar sufficient gates from the envelope fail at odd tail scale one:

    M < 20/469 < 4/91,             L < 1/49.              (2)

Nevertheless `M` is at least each of the first four frontier masses

    2/49, 138/3325, 12/287, 78/1855,

and `L>=2/469`, the width of the final `(1,67)` profile. Hence all five
mass-or-width disjunctions hold. This is a physical realization of the joint
certificate passing while each of its independently maximized scalar gates
fails; the older sufficient mass gate `4/91` fails as well. A failed
sufficient test is only uncertified, never unsafe.

**Family consequence.** For every coprime `0<p<q` with `gcd(pq,6)=1` whose
sum has only prime factors `2 mod3`, each with exponent at most two, and
for every positive odd integer `g`, the thirteen-speed row

    2H union {gp,gq}                                      (3)

has a common phase of clearance at least `1/14`.

**Proof.** At `g=1`, the five disjunctions hold by (1)--(2) and the comparisons
above. For `g>=1`, the mass is unchanged and `gL>=L`, so the disjunctions
continue to hold. Apply the proved doubled-body consumer of the envelope.
All thirteen labels are distinct: the body has eleven distinct even labels
and the tails are distinct odd labels. The row is primitive because the
body contains two and its gcd is two, while both tails are odd. No height
bound is needed for this consumer. QED.

For the extremal mass pair `(p,q)=(1,67)` and `g=1`, the literal row is

    {1,2,6,8,10,12,14,16,18,20,22,26,67}.

At `x=9/34` its clearance is exactly `2/17>1/14`. This directly verifies
one full thirteen-speed consequence. It is not the phase claimed uniformly
for every tail pair in (3).

The displayed example is **not** an actual `11+2` decoder partition:
the physical body label ten and tail label one give the crossing atlas
ratio `1:10`, whose sum eleven is an allowed inert prime. The body/tail
split in this note belongs to the physical doubled-body consumer. It must
not be relabeled as actual decoder equality entry.

## 3. Exact reproduction

[Standalone source](../../04-computation/creative_20260906_physical_control.py)
and [frozen output](creative_20260906_physical_control.out) retain all twelve
components, six isolated points, exact mass and width, every threshold
comparison, and the full-row safe phase. The finite universe is this single
explicit body and completion, with phase reflection and the isolated-point
boundary as controls. The infinite family follows from the proved envelope,
not from a finite tail scan.

    python3 -B 04-computation/creative_20260906_physical_control.py
    python3 -B -O 04-computation/creative_20260906_physical_control.py

Normal and optimized outputs are byte-identical and pass **47 always-active
gates**. Source SHA256:
`05f0d3136d271123c204b70d0dd5ddad433dfe83e501334ef4916e79f3285533`.
Output SHA256:
`3d2cbe79ea2154e4628d731e47584377c48ab37899030f6b85cab7f0fa0569b6`.
