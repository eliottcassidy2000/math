# The one reframing that bit: an existence argument for an existence problem

**Session:** mac-mini-2026-06-17-S4. **Result:** THM-526 (arc-width lemma), HYP-2578 (reframings ledger).

I spent a session throwing eight creative lenses at the last LRC(14) inequality — the
statement, after three reductions, that every covering 13-set has gap `M(S) ≥ 1/14`. Seven of
the eight bounced off. One bit. The pattern in *which* one bit is the lesson.

## What didn't work, and the shape of why

The Farey / continued-fraction lens organized the witnesses beautifully — the optimal time
`τ*` is a best-approximant, the binding pair `(a,b)` is the reflection `b ≡ −a (mod D)`, the
covering constraint excludes exactly the convergent-denominator sawtooth minima — and yielded
*nothing*. Three-distance gave nothing. The LP/averaging dual gave a clean fence but no bound.
The binding-crossing tournament's independence number turned out to be an *anti*-certificate
(the hard sets have the *largest* α). The Diophantine clearing-count argument proved the count
of good crossings is even, which is charming and useless until you can force it nonzero.

The common failure has a single shape. `M(S) ≥ 1/14` is an **existence** statement — *find a
time when every runner is far*. Almost every elegant tool I reached for is a **structural** or
**upper-bound** tool: three-distance bounds gaps from above; the reflection `b ≡ −a` holds for
*any* pair summing to the denominator and so carries zero content; the tournament/Farey lenses
describe *where* a witness would live without producing one. You cannot prove "there exists a
good τ" by describing the geometry of good τ's. The reflection `‖aτ*‖ = ‖bτ*‖` is the same
tautology as `D ≤ 14k ⟺ M ≥ 1/14` wearing a nicer hat. Structure is not a witness.

## What the winning lens did differently

The arc-width lemma is the only one of the eight that *constructs the witness*. The core `A`,
being strictly loose, has a level-`1/14` safe set that is a union of arcs; let `W` be the
widest. The parked runner `14m` has danger teeth of width `1/(98m)` spaced seven times that
apart. If `W` exceeds a single tooth, the widest arc simply *cannot* be swallowed — it pokes
into a gap, and that poke is the lonely time. No optimization, no knowing `M` in advance: a
point you can name, `τ₀ = (2k+1)/(2·14m)`, where the parked runner sits at distance ½ and the
core is safe because we're inside its arc. The conclusion `M(S) ≥ 1/14` is *derived*, and that
is exactly the non-tautological move every other angle failed to make.

It bit because it answered an existence question with an existence argument — pigeonhole, not
description. And it discharges the entire large-`m` direction with an *explicit* threshold
(sharper than the workflow's parallel derivation by a factor of seven, because catching a
*safe* point is easier than catching the *safest* point). The reward is a genuine theorem: an
infinite sub-family of covering 13-sets — the canonical hard cores `{1..13}\{j} ∪ {14m}` among
them — for which LRC(14) is now proved unconditionally.

## The honest boundary, and where it points

The lemma stops exactly where the conjecture gets hard. The safe-arc width `W(A)` is not
uniform — push a core's speeds up and its widest arc shrinks toward zero, so the threshold runs
away. What survives is precisely **large-speed cores in small-`m` resonance**, the regime where
the parked runner *is* the binding partner. That is the same stone every lens has now bared from
its own side: the measure side (THM-522) calls it bounding the lcm, the binding-pair side
(THM-524) calls it the clearing crossing, and this side calls it bounding `W(A)` below. The
distilled target is sharp now — **clearing-crossing existence**: for the binding pair
`D = flank + w`, show the crossing index `j = (flank·num mod D)` is forced `≥ D/14`. It is, once
more, an existence statement, and the session's lesson is which kind of argument to bring to it:
not a description of where the good crossing would be, but a construction that produces one.
