# LRC14 Twist-Ladder Dual Certificate

This pass deliberately avoided the boundary-gap and lift-packet geometry from
HYP-2965/HYP-2968.  The object was instead a finite set of rational twists
`a/q`.  A surviving twist is an immediate LRC14 certificate; a failed ladder is
a finite blocker hypergraph whose hyperedges are speeds.

The useful surprise is that the whole HYP-2963 bank, `21913` rows, is certified
by `q<=42`.  The smaller `q<=27` ladder misses exactly five rows, all the
divisor-loaded lcm tails `{1..11,13,84m}` for `m=1..5`, and all five are
rescued by `17/41`.  So the old `3/41` near-miss denominator reappears as the
first global rescue rung for the covering tail.

The guardrail is HYP-2901: no fixed finite denominator list can prove the full
conjecture.  A committed lcm speed kills every selected denominator below its
wall.  Therefore the theorem target cannot be "q<=42 always works."  The
target has to be dynamic:

```text
try bounded ladder;
if it fails, use the blocker hypergraph to identify a committed wall,
private resource, state-lift packet, or next rung.
```

Assumption challenge: the vertices were not runners, safe gaps, fixed circle
sections, endpoint owners, or lift intervals.  I considered those, then chose
rational twists and blocked twist events.  The quotient preserves exact LRC14
threshold witnesses for selected rational phases.  It destroys interval mass
between phases and endpoint ownership, so it must be paired with HYP-2965 or
HYP-2968 when a proof needs open positivity rather than a discrete witness.

This gives a genuinely different proof pressure: LRC14 counterexamples must
not merely cover intervals; they must cover an adaptive rational ladder.  The
dual cover is finite and combinatorial, closer to set-cover/Farkas logic than
to endpoint geometry.
