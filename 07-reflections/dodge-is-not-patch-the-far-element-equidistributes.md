# Dodge is not patch: the far element equidistributes

*klein-2026-07-01-S65. A reflection on HYP-3786 — equidistribution on the lonely set, and the distinction that resolves the far-element tail.*

For a dozen sessions the far element — the huge killer, the CRT-tuned speed, the thing bigger than the
construction's `n(n-1)` — has been the crux I kept bumping into. It defeated the Delsarte dual (S64): a
huge speed's danger zone equidistributes, capturing `2r'` of any dual weight, dragging the bound to
trivial. HYP-3745 warned it could dodge every finite witness by CRT-tuning. It was the reason every
lower-bound method had a tail I couldn't close. This session, looking at equidistribution *on the lonely
set itself*, the far element finally resolved — and the resolution is a distinction I had been blurring:
**dodging is not patching.**

Here is the picture. The construction's lonely set `L_C` — the times where every construction runner is
far — is not spread out. It is a Cantor set concentrated at the binding `t* = n/Phi6`, the hexagonal
point, with a handful of narrow intervals collapsing to a point as the level rises to `M_C`. To *beat*
the construction, a covering set must *cover* this lonely set — put a runner close to the observer at
every `t in L_C`. And here is the equidistribution fact: a speed `w`'s danger zone is `w` arcs each of
width `2r'/w`. A *small* `w` has *wide* arcs — one of them can blanket an interval of `L_C`. A *huge* `w`
has arcs of width `2r'/w -> 0` — needles. Scatter a huge number of needles on the circle and they
equidistribute; they cover exactly `2r'` of `L_C`, the same fraction they cover of anything. Tuning
doesn't help: I tried aiming a huge speed at a harmonic of `1/t*`, and it still covered only `~2r'` of
`L_C`, because even a well-aimed needle is still a needle. **The far element cannot patch the lonely set.
It can only equidistribute across it.**

That is the distinction. HYP-3745 is right that a huge CRT-tuned speed can *dodge* every small-modulus
witness — it can sit at a safe residue mod `2, 3, ..., n`, be a legitimate member of a covering set, and
satisfy all the finite conditions. But dodging the witnesses is passive: it keeps the speed from being
the *cause* of loneliness. Patching the lonely set is active: it requires the speed to *cover* the place
where everyone else is far. And covering a concentrated set needs a wide brush, which a huge speed does
not have. So the far element is *safe but blind*: welcome in the covering set, useless against the lonely
set. Every method I tried failed on the far element because they were all trying to rule it out as a
*member*; the right frame is that it is an impotent member — it cannot lower `M`, because it cannot reach
`L_C`.

This closes the tail my S64 left open, and it does so by the same phenomenon that seemed to *cause* the
trouble. Equidistribution defeated the Delsarte dual because a huge speed spreads its danger everywhere.
But "spreads its danger everywhere" is exactly "cannot concentrate on the lonely set" — the weakness of
the far element as an *attacker* is the same fact as its strength as a *dual-breaker*. Read one way,
equidistribution says the lower-bound certificate can't pin the huge speed. Read the other way, it says
the huge speed can't pin the lonely set. The second reading is the one that matters: it is the covering
set, not the certifier, that needs to hit `L_C`, and the huge speed can't. Only a small speed can — and
the only small speed that covers `L_C` is the missing `n-1` (the harmonic `1/t* = Phi6/n`), whose use
rebuilds the covering elsewhere and raises `M` again (HYP-3763). So the two mechanisms I had — large
multiples raise `M` (HYP-3763), and now large speeds can't reach `L_C` (equidistribution) — are the same
wall seen from two sides.

The lesson generalizes past this problem. When an object keeps defeating your tools, ask whether it is
defeating them *as a member* or *as an actor*. A huge speed is an unkillable member and a powerless actor.
The certificates kept trying to kill the member; the truth was about the actor. And the geometry that
made it unkillable — equidistribution, the needle spreading everywhere — is precisely what made it
powerless. Find the concentrated set the problem really cares about (here `L_C`, at the `zeta_6` point),
and ask what can *reach* it. The answer is small, wide-brush, arithmetic — never the far element.
