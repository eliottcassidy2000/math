# Consolidated structure analysis of the collapse family (S621)

Last session (S617) I wrote, with some pride, that the collapse family "turned out, on an
exhaustive small search, to be *exactly* the additive chains." This session I gave that sentence
the search it deserved, and it broke — in the most instructive way.

The first thing that changed was the framing. "Collapse" (`p₀ = 0`) is not a combinatorial
curiosity; it is the statement that the loneliness gap is pinned at exactly `1/m` — the width-`2δ`
arcs cover the circle with not a point to spare. That is the *exactly-tight* / *barely-lonely*
configuration of the actual Lonely Runner literature (Kravitz; the Perarnau–Serra survey). So the
sub-problem the user pointed at is a classical object wearing the covering-depth costume, and the
right questions are the classical ones: how many are there, how big can they be, where do they
achieve the gap.

All three have clean answers, and none of them is "additive chain."

**How many, how big.** The collapse family is *finite* per `n`, and tiny: `1, 2, 2, 1, 3, 1` for
`n = 3..8`. Every member obeys `v_max ≤ 2n-1`, with the sporadics hitting the bound exactly
(`(1,3,4,7), (1,3,4,5,9), (1,4,5,6,7,11,13)`). The count is range-stable — push the search to
`8n` and nothing new appears — and it is not in OEIS. That `2n-1` is, up to the usual
normalization shift, the `2n-2` tight-instance barrier the literature already knew as the wall
every proof of LRC runs into. It is satisfying to meet that wall from the covering side, as a
plain bound on how large a speed can be before its thin comb of arcs leaves a gap nobody can fill.

**Where.** This is the part that pays. Every tight set — the AP and all the sporadics — achieves
the gap on *exactly the same* set of times: `{k/m : gcd(k,m)=1}`, the `(ℤ/m)*` orbit, `φ(m)/2` of
them. The sporadics are witness-equivalent to the AP. That is not a new object; it is THM-403's
witness orbit, the very `(ℤ/m)*` of the perspective key, now shown to govern not just the AP but
the whole extremal family. A tight instance is precisely an integer lift that keeps the gap pinned
across the entire primitive-residue orbit. Rigidity = orbit-type, read on the clock.

**Why additive chains were a mirage.** The subadditivity mechanism is real: `‖(a+b)t‖ ≤ ‖at‖ +
‖bt‖`, so a relation `a+b=c` chains three arcs and a tight boundary pins the third. But that pins
*one* arc. Covering the *whole* circle needs every arc pinned simultaneously across the whole
witness orbit — incomparably stronger. The numbers make the gap visceral: at `n=7` there are 3
collapse sets, 364 strict additive chains, and 5858 sets whose largest speed is a sum of two
others. Additive richness is a *symptom* of collapse (the clean necessary condition is just
`max(S) ∈ S+S`, which every tight set satisfies), never the cause. I had mistaken the smoke for
the fire. The honest statement is one-directional: collapse ⟹ resonance-rich, and not back.

The meta-lesson is about exhaustiveness. "Exhaustive search" with the wrong bound (`max ≤ 12`)
returned a family small enough to pattern-match to additive chains, and the pattern was a
coincidence of the truncation. The moment the true bound `v_max ≤ 2n-1` is in hand, the family is
*more* constrained, not less, and the constraint is arithmetic-geometric (the witness orbit), not
additive-combinatorial. A characterization that is necessary-but-not-sufficient feels like a near
miss; here the sufficiency failed by three orders of magnitude. Better to know.

What remains is to prove what the computation sees: the `2n-1` bound, the `(ℤ/m)*` witness
characterization, and the even-`m` uniqueness dichotomy (`n=3,6,8` give only the AP — the
2-adic-seam fingerprint again). Those are HYP-2196. The covering-depth master object stands; the
additive-chain story is corrected; the perspective key turns out to be holding the door.
