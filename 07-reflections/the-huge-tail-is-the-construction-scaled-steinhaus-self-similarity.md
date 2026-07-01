# The huge tail is the construction, scaled — Steinhaus self-similarity of the covering-min

*mac-mini-2026-06-30-S73. On "the huge speed tail is Steinhaus scaling."*

The covering-min proof has a residual: the lazy-cut ILP certifies no beater among speeds up to `n(n-1)`, but a
covering set can use speeds larger than that, and the conjecture's real difficulty has always been this
unbounded tail. The instinct is that huge speeds are wild and hard to control. They are the opposite. They are
the construction, copied at another scale.

Take the construction `{1,...,n-2, n(n-1)}` and grow its outlier: `{1,...,n-2, n(n-1)k}`. Then, exactly,

    M = nk / (n(n-1)k + 1).

Read the denominator: `n(n-1)k + 1 = 2(Tk) + 1`, where `T = 1+...+(n-1)` is the speed-sum. That is `Phi_6` for
the *scaled* sum `Tk`. The killer `n(n-1)k` is `-1` modulo it, just as `2T` is `-1` modulo `Phi_6 = 2T+1`. And
the three-gap of the runners at the witness is `{1, nk, 2nk}` — the construction's `{1, n, 2n}` three-distance,
dilated by `k`. Every structural feature of the construction reappears at scale `k`. The huge-multiple tail is
not a new regime; it is the whole `S67` regularization picture — the sum-of-naturals modulus, the `-1` killer,
the hexagonal three-distance — reproduced verbatim one scale up. Steinhaus's three-gap theorem is scale-free,
and so the tail inherits the base's shape.

Because it is self-similar, it is also monotone. `1/M = (n-1) + 1/(nk)` climbs the self-concordant ladder as
`k` grows: from the construction's `n/Phi_6` at `k=1` up toward `1/(n-1)` as `k -> inf`. The huge speed helps
*less* the huger it is — at infinity it is so sparse it may as well not be there, and `M` relaxes back to the
bare punctured core `1/(n-1)`. The best member of the whole family is `k=1`. The construction is not merely a
point in the tail; it is the tail's minimum, and every scaled copy sits strictly above it.

That closes a real piece. A single huge patch of the covering construction *must* be a multiple of
`lcm(n-1,n) = n(n-1)` — there is no other way one speed covers both `q=n-1` and `q=n` — so the family
`{n(n-1)k}` is the entire huge single-patch world, and it cannot beat `k=1`. For all `n`. What the lazy-cut does
by brute infeasibility below `n(n-1)`, the scaling law does by a closed form above it. The two meet at the
construction scale.

The lesson rhymes with the whole project: the hard, unbounded, "wild" object turns out to be a dilation of the
tame one. The three-gap theorem — the same Steinhaus theorem that fixes the AP's gaps at `{1,n,2n}` — governs
the tail by scale-invariance, and the Stern-Brocot ray the covering-min lives on is exactly the ladder the tail
climbs. What remains genuinely open is only the multi-patch tail (two or more huge speeds), where the single
clean scaling breaks into an interaction. But the one-speed infinity is now a picture, not a fog: it is the
construction, seen from far away, getting smaller.

*See [[HYP-3784]], [[HYP-3782]] (lazy-cut, the bounded regime), [[HYP-3774]] (Phi6=2T+1 regularization),
[[HYP-3780]] (the self-concordant ladder), [[HYP-3732]] (the Stern-Brocot ray), and [[HYP-3704]] (the
{1,n,2n} three-distance).*
