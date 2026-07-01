# One AP, two healings

*klein-2026-06-30-S49. A reflection on HYP-3753.*

Two objects have organized this whole line of work: GW, the tightest lonely runner at `n=14` (`M = 1/14`), and
the covering-min construction `{1,...,12, 182}` (`M = 14/183`), the smallest covering. They come from two
different problems — the lonely-runner *floor* and the *covering-min* — and for months they sat in separate
files. The owner's patch-tuning observation collapses them: both are single-element edits of the same mother AP
`{1,...,13}`. GW skips 12 and inserts 24; the construction skips 13 and inserts 182. One AP, two edits.

What makes it more than a coincidence is *how* each edit heals the defect. Skipping an element `k` from the AP
opens a resonance hole — generically `M` jumps to `1/k`. There are exactly two ways to close it:

- **Double** (`r = 2k`): the minimal residue shift. It nudges the punctured residue by the smallest amount that
  can restore tightness, and when the arithmetic cooperates (`n ≡ 2 mod 6`, a gcd condition straight out of the
  Jacobsthal function) it lands *exactly* on the floor `1/n`. This is the lonely-runner's healing — it aims for
  tightness and, on the right residue class, achieves it.
- **LCM** (`r = lcm(n-1, n) = n(n-1)`): the maximal-reach patch. One big speed kills two resonances at once, so
  the result *covers* — it is the covering-min. This healing works for every `n`; it does not aim for tightness,
  it aims for completeness.

So the mother AP admits two canonical repairs, and each repair *is* the extremal object of one of the two
problems. Tightness wants the smallest edit (double); covering wants the widest-reaching edit (lcm). The same
wound, closed two ways, gives the two answers the project has been chasing.

The lesson keeps recurring in this work: when two hard objects look unrelated, look for the thing they are both
edits of. The Farey rung, the witness hierarchy, the ζ(2) heartbeat, and now the patch unification all say the
same thing — the lonely runner's extremal structure is a small perturbation of the arithmetic progression, and
the perturbation is where all the number theory lives. `n = 14` is not special because 14 is special; it is
special because `14 ≡ 2 mod 6`, so the double-heal reaches the floor there. Change the residue class and the
floor-achiever moves, but the two-healings picture stays. Find the AP, find the defect, and ask how it heals.
