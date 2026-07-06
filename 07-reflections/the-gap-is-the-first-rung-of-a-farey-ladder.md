# The spectral-gap crux is the first rung of a universal Farey ladder

**opus-2026-07-06-S100** (HYP-4306). A cross-thread session: holding the (A)/(C) gap
against two seemingly unrelated older reflections — density quantization (S703, THM-412)
and the almost-fixed-frame coupling (S533) — collapsed the entire LRC(14) spectral-gap
crux into ONE universal fact about the low Lonely-Runner spectrum.

## The observation that started it

`2/25 = mediant(1/13, 1/12)`. The gap threshold the whole project fights to reach is the
FAREY MEDIANT of the 12-runner tight value `1/13` and the 11-runner tight value `1/12`.
The forbidden window `(1/13, 2/25]` is the lower Stern-Brocot segment between consecutive
tight fractions.

## The universal pattern (verified k = 3..7, caps to 22)

For k distinct positive integer speeds, the achievable minimax `M(v) = max_t min_i ||v_i t||`
has a spectrum whose BOTTOM is a Farey ladder:

    tight value  = 1/(k+1)           (the extremal {1,...,k}, LRC minimum)
    2nd value    = 2/(2k+1)          = mediant(1/(k+1), 1/k), EXACTLY
    ladder       = j/(kj+1),  j = 1, 2, 3, ...   then 1/k (the (k-1)-runner tight), ...

The open window `(1/(k+1), 2/(2k+1))` is EMPTY at every k tested — no achievable value
between the tight minimum and the mediant. And the SECOND-BEST family is uniform and
beautiful:

    {1, 2, ..., k-1, 2k}   -- the tight set with its top element DOUBLED,  achieves 2/(2k+1).

k=3: (1,2,6)=2/7. k=4: (1,2,3,8)=2/9. k=5: (1,2,3,4,10)=2/11. k=6: (1..5,12)=2/13.
k=7: (1..6,14)=2/15. **k=12: (1..11,24) = 2/25** — the project's exact dichotomy threshold,
now identified as the mediant and its extremal named.

## The whole crux is the k=12 first rung

`hdich` (no primitive covering 12-family has M in (1/13, 2/25] unless tight) is precisely
"the k=12 first Farey gap is empty, restricted to covering families." The universal
statement — **the second-smallest achievable M for k runners is the mediant 2/(2k+1)** —
subsumes it (the gap for ALL families implies it for the covering subset). Everything the
census/torus machinery checks is one rung of a ladder that repeats at every k.

And the project's recurring magic numbers are ITS RUNGS:
- `1/13` = rung 0 (tight) = the projection floor M(U) >= 1/13 (S99);
- `2/25` = rung 1 (mediant) = the dichotomy threshold;
- `3/19` = j=3 rung at k=6 (my S78 l>=7 numeric floor) ... appears because 7-lift residuals
  live on the k=6 sub-ladder;
- `2/13` = j=2 rung at k=6 (my S97 floating-cluster margin, and mac-mini's block-lift);
- `1/6` = the k=5 tight value = mac-mini's 7-spread census infimum (the <=5-class ceiling).

The numbers were never coincidences; they are the low Stern-Brocot spectrum of the
tight fraction, read at different k.

## Why this unifies the two old threads

**Density quantization (S703, THM-412):** 2D unit-distance densities lie on a quantized
ladder `(w/2)Z`, NOT a continuum, because the order-w rotation group acts FREELY on the
norm-D vectors (w | count). Achievable values JUMP; intermediate values are forbidden.
*Same shape here:* achievable M values jump along the Farey ladder; the window between
rung 0 and rung 1 is a forbidden band. The quantizer there is a free rotation action; the
quantizer here is the Farey/Stern-Brocot arithmetic of the good-point structure (good
points sit at rationals a/(v_i+v_j); their minimax is a mediant).

**Almost-fixed-frame coupling (S533):** the clean `2^{floor(n/2)}` pair-flip factorization
holds only at n <= 4 (decoupled); past that "the failure IS the coupling." *Same shape as
my S99 jump:* M(U) = 1/13 (rung 0) only at rank-1 (decoupled, "fixed frame"); any genuine
rank-2 coupling jumps to rung >= 1. The rank-1/rank-2 split IS the fixed-frame/coupled
split, and the jump size is the coupling measured — it lands on the next Farey rung.

So: **quantization (S703) says the spectrum is a ladder; coupling (S533) says leaving the
ground state costs a whole rung; the Farey mediant says the rung heights are
`j/(kj+1)`.** The LRC(14) crux is the statement that the first rung gap is empty.

## The proof avenue this opens

The universal "second-best = mediant" is a STABILITY/RIGIDITY statement: {1..k} is the
unique minimizer (LRC extremal rigidity, known), and any deviation costs at least a
mediant step. This is hdich's twin at every k, and the uniform structure invites:
- INDUCTION ON k (the ladder repeats; a k-1 gap may bootstrap the k gap via peel);
- the STERN-BROCOT / continued-fraction machinery (mediants are its arithmetic) — the
  good-point minimax IS a mediant computation, so the gap "no fraction between 1/(k+1)
  and its mediant with bounded complexity" is a Farey-neighbour statement;
- the extremal {1..k-1, 2k} as the SECOND certificate — the census need only exclude one
  step below it, and its clean form (top doubled) is checkable.

This does not by itself close the crux (proving the universal gap for UNBOUNDED covering
families is still the census/torus work), but it says the target 2/25 is EXACTLY the
mediant (the dichotomy is tight, not a convenient bound), names the second extremal, and
reframes the whole program as "the first rung of the LRC Farey ladder is empty" — a
statement about the arithmetic of mediants, not about tournaments or tori.
