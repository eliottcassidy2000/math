# The clean answer is the loose one

*klein-2026-06-30-S53. A reflection on HYP-3764, the covering-min open edge, and a pattern the project keeps re-learning.*

There is a shape this research keeps falling into, and it is worth naming because it has cost several
sessions and two logged mistakes. You find a beautiful covering set — the construction `{1,…,n-2,
n(n-1)}`, whose gap is `n/Φ₆(n)`, the sixth cyclotomic value, the hexagonal lattice, the Eisenstein
integers, a clean closed form. It sits right where a good answer should sit, a hair below `1/(n-1)`.
You verify it is the tightest thing your search can find. You call it the covering-min. And it is not.

This session pinned *why* with an embarrassingly simple observation. Every covering set's gap is a Farey
rung above the floor: `M = r/(r(n-1)+1)`, and this is **increasing** in the rung `r`. The floor `1/n`
is rung 1; the loose ceiling `1/(n-1)` is rung infinity; and the construction is rung `n` — almost all
the way up, near the loose end. The construction is not the tightest covering. It is nearly the
*loosest*. Tightness lives at *small* rungs, and small rungs are spread sets: `{1,3,4,5,7,11,18,32}`
beats the construction at `n=9`, `{1,2,3,5,6,7,8,9,30}` at `n=10`. The clean, dense, cyclotomic object
is the easy upper bound; the extremum is a scattered thing with a killer at `3n` and a binding pair at
some prime nobody would guess.

Why does the clean object masquerade as the answer? Because it is the only one a *local* search finds.
The construction is the fixed point of "add the densest core, patch the two loose resonances with their
lcm." It is what you converge to from almost any starting set, because it is the local basin every
greedy move drains into. The true covering-min sits in a different basin — a specific spread
configuration — that you reach only if you already know where to look. This session made that concrete
and, I think, decisive: I ran three searches (random, hill-climb, structured drop-≤2), and *every one*
of them missed a beater that is known to exist. Random missed all of them. Hill-climb found `n=7,8,10`
and missed `n=9,11`. The targeted drop-≤2 found `n=9,10` and missed `n=11`. Three methods, three
different blind spots, each certifying "no beater here" at some `n` where a beater provably lives.

That is the real content of the "open edge." The repository believes the construction becomes the
covering-min at `n≥12` — a transition, driven by the radius-1 band getting over-constrained. Maybe. But
the *evidence* for that belief is exactly the kind of search that just failed at `n=11`. You cannot use
a method that misses the `n=11` beater to certify there is no `n=14` beater. The transition might be
real; it might be the horizon where our searches stop being able to see. Right now those two are
observationally identical, and honesty demands we say so. My own bet — and it is a bet — is that it is
the horizon: the rung ladder `2,2,4,4,3` has no reason to leap to `n` at `n=12`, and "the answer got
clean exactly when the search got hard" is the signature of a mirage, not a theorem.

The saving grace, and the reason this does not threaten the conjecture it looks like it should, is that
the *whole question is orthogonal to the Lonely Runner Conjecture.* LRC(n) asks only whether the
covering-min stays above `1/n` — whether the smallest rung is at least 2. And it is, unconditionally,
because covering sets have `M > 1/n` strictly (rung 1 is the floor, forbidden). So `covering-min ∈
[2/(2n-1), n/Φ₆]`, both a whisker above `1/n`, and LRC holds with margin `1/(n(2n-1))` no matter which
rung is extremal. The thing we cannot compute (the exact rung `a(n)`) is not the thing we need (rung
`≥ 2`). The razor-thin edge is real, and it is genuinely open, and it does not matter for the
conjecture — three facts that are easy to conflate and worth keeping apart.

The general lesson, the one behind MISTAKE-087 and -088 and now this: when an optimization over the
integers hands you a clean closed form, suspect that you have found the *boundary of your search*, not
the *extremum of the problem*. The clean answer is usually the loose one. The real extremum is where
the number theory is — irregular, spread out, with no closed form, sitting one `1/n²` step past the
elegant object you were admiring. `width(G_n)` abandoned `C(n-2,⌊·⌋)` at `n=7`; the covering-min
abandons `n/Φ₆` somewhere too, and probably always did. Find the clean object, then distrust it, then
go looking in the mess just beyond it.
