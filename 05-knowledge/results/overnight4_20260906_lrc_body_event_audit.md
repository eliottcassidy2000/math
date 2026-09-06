# Independent audit of the actual body-to-event interface

**Status: INDEPENDENT AUDIT PASS.** The elementary event-interface proof
and all declared finite controls in
[the producer report](overnight4_20260906_lrc_body_event.md) are accepted.
No source correction was needed. The actual-entry obligation remains open;
the event criterion is sufficient and is explicitly not necessary.

## Analytic scope

At `y=(14k+3)/(14w)`, a ternary-unit tail w has no strictly bad physical
lift: its minimum of the three clearances is exactly `1/14`. Each other
ternary-unit tail can spoil at most one lift, since the three phases are
separated by `1/3` and the strict danger interval has length `1/7`.
Choosing a missing owner therefore gives a weak-safe physical phase
`x=(y+j)/3`. The divided body is preserved exactly because
`(3c)x=cy+cj`. This proves the claimed actual-witness implication, including
its endpoint semantics; an arbitrary lift need not be safe.

For `g=gcd(c,w)`, `q=w/g`, `h=c/g`, the event values modulo `14q` are
`h(14k+3)`. Since `gcd(h,q)=1`, they run over the complete coset `3h mod14`,
with multiplicity g. The bad representatives satisfy `-q<n<q`; counting
that progression gives the displayed two-floor formula. The count is at
most `g*ceil(q/7)`. For `q>=2` this is less than w. For `q=1`, the sole
value is bad exactly when `14|h`, giving the exact full-cover boundary
`14w|c`. The strict `q-1` endpoint is necessary.

Consequently, the four divisible frequencies in the ten-body divisor gate
are harmless when none is divisible by `14w`. Each of the other at most
six frequencies blocks at most `ceil(w/7)` points under its coprimality
assumption. The strict union bound is valid. The seven-residue argument
proves its all-integer cutoff 37; filtering out multiples of three gives
the uniform ternary-unit cutoff 31. The failures at 36 and 29 concern the
count statistic, not sharpness of the underlying event noncoverage claim.

The phase-loss example at w=5 is correctly typed: the gcd and h modulo
fourteen determine each marginal count but not the union of its affine
event words. Keeping c modulo `14w` repairs that loss. The positive
divisor gate is an inherited event-interface specialization, not a general
entry theorem or a newly claimed LRC family.

## Independent exact paths

The [audit source](../../04-computation/overnight4_20260906_lrc_body_event_audit.py)
does not import the producer. It constructs the sorted rational event
points and finds those in each open danger interval
`((14m-1)/(14c),(14m+1)/(14c))` by binary search. Both boundary comparisons
are strict. This replaces the producer's modular-residue membership test.

Every full event word agrees on all **4,920** pairs: `1<=w<=61`, w a ternary
unit, and `1<=c<=120`. The exact count formula, gcd bound, and single-owner
iff pass on every row. The whole-word semantic digest, not merely a count
digest, is
`bbeeb45e66237a65b66f0837833787a2bb2275abd5ae5f1eb1c7bc8a306d3671`.

For the w=43 positive control, the interval path independently recovers
all sixteen survivors:

```text
1,2,7,8,10,11,20,22,25,26,28,29,37,40,41,42.
```

Its marginal blocker sum is 37 and actual union size is 27. The phase
`y=17/602` is body-safe. Direct physical danger intervals give safe lift
labels `{1,2}`; the stated label-one witness `x=619/1806` has exact minimum
clearance `1/7`, which is stronger than required. Every small-clock divisor
from two through fourteen occurs among its thirteen physical speeds.

The all-event-cosets hostile also passes independently: the frequencies
`14w` block both endpoint signs for tails 1,5,11, while `x=1/13` is safe
for the full physical row with minimum clearance `1/13`. This directly
preserves the nonnecessity boundary of event incidence.

## Full packets and the norm-five boundary

For `C={1,...,10}`, the full strict simultaneous return set of radius
`delta=3/154` is the centered interval of radius `3/1540`. Indeed c=1
first gives `|u|<delta`; `10delta<1/2` prevents any wrap for c=10 and
forces `|u|<delta/10`. The converse pays all smaller speeds. Both chosen
anchors have exact deep margin `1/11`, so equal packet size alone cannot
distinguish them.

The audit reconstructs the closed body component around `2/11` directly
from rational safe intervals as `[5/28,13/70]`. The full translated packet
`(277/1540,283/1540)` lies strictly inside it. On physical lifts j=0,1,2,
the entire packet is respectively contained in danger teeth
`(frequency,tooth integer)=(1,0),(5,2),(11,8)`. This is a direct physical
interval certificate, independent of the producer's effective nearest-
integer owner computation. Translating the same packet to `1/11` makes
tail five inactive on all three lifts. The anchor, rather than packet
measure or body lattice alone, is necessary information.

For the norm-five control `(10,11,16)`, the unit-residue defect restriction
and `|delta|<15/14` force zero defect. The complete multiplier list is
`k=+-1,+-2` because the tight roof cutoff is `63/28`. I independently use
the effective intervals centered at nearest integers `(2k,2k,3k)`, compute
their literal pairwise intersections, and cap by the shortest interval.
They give

```text
E=(17/176,9/140,3/55), physical mass=331/6160,
min E - physical mass=1/1232.
```

The minimizing omitted coordinate is the second at k=1 and the third at
k=2. All projection columns, physical mass, and four-carrier count agree
with the [native-audited H63 table](overnight4_20260906_lrc_parityfree_probe_head63.tsv).
Thus the fixed norm-four projection identity cannot be transferred to
norm five; a scalar central-section integral does not repair that switch.

## Reproduction and limits

[Frozen audit output](overnight4_20260906_lrc_body_event_audit.out) retains
**14,805** explicit gates. Normal and optimized outputs are byte-identical.

```text
python -B 04-computation/overnight4_20260906_lrc_body_event_audit.py
python -B -O 04-computation/overnight4_20260906_lrc_body_event_audit.py

SHA-256, LF bytes:
source 195e7c0df5226a58785ebe1a545bb11e9d4b7fc6910c637e1c20f36aee751632
output 78c8f65e436fa88a352516c5f8c28173fd917debc43b1a2e496bc9d7ba71416f
native-audited head c3d33fdd136245aafe512b04963a6eb6f1b5db6f1a572a3e8535ef59d01a09fa
```

The interval replay is a complete verification of the declared bounded
marginal universe and the individual controls. It is not a census of
all bodies. No body Haar floor, universal event incidence, choice-free
deep-anchor argument, arbitrary entry, or LRC(14) closure follows.
