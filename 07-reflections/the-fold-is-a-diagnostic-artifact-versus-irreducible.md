# The fold is a diagnostic: coordinate-artifact versus metric-irreducible

*klein-2026-07-01-S79. A reflection on HYP-3812 — carrying the tournament fold-thinking across to the
lonely runner, and learning as much from where it fails as from where it works.*

The owner asked me to apply the tournament quarter-fold to LRC-14 and to do the complement-fold. The
instinct was to hope the same magic would happen: last session an obstruction that looked real in one
coordinate system — the self-complementary classes resisting low-dimensional covers — simply dissolved when
I moved to the fold-adapted (half-address) coordinates. It was a coordinate artifact, not an obstruction. So
the natural hope was that the lonely-runner covering-min, which has been the open crux for months, might
dissolve the same way under the right fold.

The folds transfer cleanly at the level of structure. The lonely runner has its complement: the antipode
`t -> 1-t`, under which the loneliness profile, the danger pattern, and the lonely measure are all
invariant, so the whole problem halves onto `[0, 1/2]` — the two binding times `t*` and `1-t*` fold to one,
the Verblunsky coefficients are real, the parity lemma is the fixed-point count. And there is even a quarter:
the complement `-1 mod Phi6` factors, by the Chinese remainder theorem over the Eisenstein primes `3` and
`61` of `Phi6 = 183`, into two partial complements, and those generate a Klein four-group in the units of
`Z/Phi6` — the exact algebraic twin of the tournament's `<sigma, flip>`. So the fold is there, both halves
of it.

But the obstruction does not dissolve, and the way it refuses is the real content. I measured how lonely the
construction can be at each modulus separately. At `2, 3, 6, 7, 14` it is not lonely at all — the small
speeds cover those; at `61` it is only shallowly lonely; and the full covering-min value appears **only** at
the composite modulus `Phi6 = 3 * 61`. The construction covers every prime factor and binds only at their
product. So the covering-min is metric-irreducible: it lives at the deep composite modulus, and the CRT
quarter-fold, which sees only the prime factors, cannot reach it. Where the tournament's hardness was a
mirage that a change of basis erased, the runner's hardness is genuine and sits precisely at the place no
factorization can touch. This is the same fact as the old theorem that the singular series is not an Euler
product, but seen from the covering side and made vivid: the binding is at `Phi6`, not at `3` or `61`, so
there is no shortcut through the primes.

That turns a hoped-for dissolution into a diagnostic, which is more valuable than the dissolution would have
been. The fold-thinking now sorts obstructions into two kinds. A coordinate-artifact obstruction dissolves
under the fold — it was never in the object, only in the axes — and such problems are easy once you find the
natural coordinates. A metric-irreducible obstruction survives every fold, because it lives at a scale the
symmetry cannot factor, and such problems are genuinely hard and must be attacked at that scale. The
tournament SC-cover was the first kind; the lonely-runner covering-min is the second. Knowing which kind you
face is knowing whether to hunt for a change of basis or to settle in at the irreducible scale and build the
certificate there. For the runner, the answer is the latter: the proof must live at the full `Phi6`, which
is exactly where the Chebyshev two-point dual sits — the alternation on the slowest runner and the killer,
at `t*` and its complement. The fold told me not to look for a shortcut through the primes, and that
negative is a map: it says the whole difficulty is concentrated at one composite modulus and its
two-point boundary, and nowhere else.

The methodological keeper: when a technique works on one problem, port it to the next not to win but to
measure. If the obstruction dissolves, you have found the coordinates; if it refuses, you have found the
scale. Either way the fold has done its job — it has told you what kind of hard the problem is. The runner
is metric-hard at `Phi6`, and that is the most precise thing I can say about why it has resisted, and the
most useful thing to know before spending the next effort.
