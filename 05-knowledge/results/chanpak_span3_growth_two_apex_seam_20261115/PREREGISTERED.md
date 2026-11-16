# Pre-registered predictions for the runs still in progress (committed 2026-11-15 before any of their outputs was read)

Session klein-2026-11-15-ribbon. Four computations were launched and are running at the time of this commit:
`hratio mats3.txt 9` (the exact one-vertex growth ratio maxima of literal span-3 words of lengths 8 and 9),
`two_apex 11` (every two-apex transitive tape with spine 11), `two_apex_small 12 3` and `two_apex_small 13 3` (two-apex
tapes whose second apex beats at most three spine vertices, spines 12 and 13), and `image_orbit.py 6000` (the orbit of
image subspaces of span-3 prefix products, capped at 6,000 subspaces). Predictions:

| # | prediction | expectation from |
|---|---|---|
| Q1 | the maximum of H(Wl)/H(W) over words of length 8 and of length 9 is 39, attained by W = 0^(n-3) 7 6 0 with l = 6 (HYP-9802's P105) | lengths 5-7 give 39 with exactly these maximisers; H(W) = 1 there |
| Q2 | the maximum of the two-fragment part Z_l(W)/H(W) at lengths 8, 9 is 38 (so H(Wl) = 1 + 38 at the maximiser: one appended path, no prepended path, 38 bridged two-covers) | lengths 6, 7 |
| Q3 | two-apex tapes with spine 11 (all 2^23 tapes) attain every odd N < 2^15 except 7 and 21, and the first miss above 21 lies in [45,000, 70,000] | spine 10: every odd N < 2^14 except 7, 21, first miss 27,141 (max H = 69,485); the small-overlap spine-11 run already read gives every odd N < 2^15 except 7, 21 with first miss 50,023 |
| Q4 | with the second apex beating at most 3 spine vertices, spine 12 attains every odd N < 2^16 except 7, 21, and spine 13 every odd N < 2^17 except 7, 21 | spine 11 reached 2^15 (first miss 50,023 = 1.5 * 2^15); the maximum count grows about 2.8x per spine step and the covered range about 2x |
| Q5 | the values needing overlap exactly 3 remain a small positive fraction (between 0.1 % and 2 %) below 2^16 at spine 12, and overlap at most 3 keeps sufficing below 2^16 | spine 11: 105 of 16,382 values below 2^15 need overlap 3 (0.6 %); overlap 2 covered 98 % at spine 10 |
| Q6 | orientation A -> B suffices for all but a handful (< 10) of the odd values below 2^16 at spine 12 | spine 11: one value (32,507) below 2^15 needs B -> A |
| Q7 | the image orbit is NOT finite within the cap: the number of distinct image subspaces of span-3 prefix products exceeds 6,000 without closing (so an exact rank-tagged automaton on image subspaces is infinite, and the rank potential must be coarsened) | the nilpotent/expanding structure of the 119-state transfer; the orbit of a positive semigroup on subspaces is rarely finite |
