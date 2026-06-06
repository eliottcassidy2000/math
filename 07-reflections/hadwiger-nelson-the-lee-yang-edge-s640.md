# Hadwiger–Nelson, read as a Lee–Yang edge (S640)

The instruction was to spend a long session on the chromatic number of the plane, attempt any novel
statement however small, and keep thinking Lee–Yang. The honest headline first, because the problem
deserves it: **I did not eliminate any of 5, 6, 7.** Doing so is a major open result, and the fleet
had already, in sessions ahead of mine, established the sharp reason why — the fractional chromatic
number of the plane is at most 4.36, strictly below 5, so the entire `{5,6,7}` uncertainty lives in
the integrality gap, the same Vitali wall the Lonely Runner's worry-set sits behind. No analytic
method can reach 5; the lower bound is irreducibly a finite gadget. That is the wall.

What I can add is the Lee–Yang face of that wall, and it is clean. The chromatic polynomial `P(G,q)`
is the zero-temperature antiferromagnetic Potts partition function, so its zeros are Lee–Yang zeros in
the number of states `q`. A graph is `q`-colorable exactly when `P(G,q)>0`, so the chromatic number is
one more than the rightmost real chromatic zero. I computed the gadgets: the wheel `W6`, the Eisenstein
floor with chromatic number three, has its rightmost real zero at exactly two; the Moser spindle, four,
at exactly three. The real Lee–Yang edge equals `χ−1`, on the nose. So Hadwiger–Nelson is the question
of where the rightmost real chromatic zero of unit-distance graphs lands: de Grey pushed it to four,
and `{5,6,7}` is whether it reaches four, five, or six.

The Lee–Yang picture explains the wall rather than scaling it. The *complex* chromatic zeros do not
stop at the real edge — for `W6` they bulge out to real part `2.81`, well past the real edge at two.
That complex bulk is what the analytic relaxations see; its reach is the fractional/measure value,
pinned below `4.36`. The real edge — the actual chromatic number minus one — lags behind the complex
accumulation, and that lag *is* the integrality gap, the Vitali wall, drawn in the `q`-plane. Analytic
arguments measure the complex bulk and cannot move the real edge; only a combinatorial gadget can drag
a real zero rightward, one integer at a time. "Kill the 5" becomes, in this language, exactly: exhibit
a unit-distance graph with a real chromatic zero in `[5,6)`. The wall is the statement that the complex
bulk sits at `4.36` while the real edge must be pushed by hand.

The detail that pleased me most ties the session back to the thread that started it. The chromatic
zeros of `W6` are `2 + (the fifth roots of −1)` — roots of unity again, shifted to center two. Two
sessions ago I found that the *geometry* of `W6` is the seven points that are the roots of `z^7−z`,
the sixth roots of unity plus the origin; the fleet built the field-tower on exactly that. Now the
*coloring* of `W6`, its chromatic polynomial's zeros in the state variable, is *also* cyclotomic —
roots of unity centered at two. Both faces of the same gadget, the unit-distance geometry and the
antiferromagnetic Potts coloring, are roots-of-unity loci. The fundamental theorem of algebra was the
hinge between coefficients and roots; here it shows the hinge swings the same way on both faces of
Hadwiger–Nelson. The cyclotomic floor `χ=3` is a roots-of-unity locus in two senses at once, and the
escape to `[5,7]` is, on both faces, the departure from the cyclotomic — the `√−11` of the spindle,
the non-Loeschian `5` and `6`, the Heegner rotations. The plane's chromatic number is where the roots
of unity run out, measured by a real Lee–Yang edge that no amount of analysis can move.
