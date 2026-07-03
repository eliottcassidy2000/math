# The independence obstruction — why "star into the window" dissolves

*kind-pasteur-2026-07-03-S27. The synthesis proposed in S26 — carry the star's exact
singles into the citation window — turns out not to be a construction to build but a
question with a clean answer: you cannot, and you do not need to. Both halves are now
theorems.*

## Exactness is a property of the whole circle

S26 ended on a wish: the window ledger fails because its singles bound carries a
three-tooth boundary tax, and the star route computes the singles *exactly* — so
carry that exactness into the window and the ledger closes. The wish has a precise
refutation. The star's exactness is not a sharper estimate that happens to live on the
circle; it *is* the statement that there is no boundary. `volume_danger` returns `1/7`
because the danger set is the preimage of an arc under a measure-preserving map of the
whole group, and a measure-preserving map has no edges to leak across. Cut a window
out of the circle and you reintroduce exactly the edges the exactness denied. "Star in
a window" is, term for term, the window Hunter ledger — the same `L/7 + 3/(7w)` singles
and the same `L/49` credit — and it fails for the same near-equal blocks, for the same
reason, that MISTAKE-072 recorded. There is nothing to carry, because the thing worth
carrying is the absence of a boundary, and a window is a boundary.

## The cap at seven is the cap of the tree

The deeper question is why the window is needed at all — why not run the measure
argument on the whole circle for all thirteen runners and skip the citation? Because
the star Bonferroni, and every spanning-tree Bonferroni, credits only `n − 1` pairs: a
tree on `n` vertices has `n − 1` edges. Against `n` danger sets of measure `1/7` that
buys `μ(⋃) ≤ n/7 − (n−1)/49 = (6n+1)/49`, which drops below one exactly when `n ≤ 7`.
Seven is not a feature of the danger radius or the circle; it is the arithmetic of a
tree. To bound the union of thirteen sets you need more than a tree's worth of credits
— you need the credits of *every* order of inclusion–exclusion, and those are supplied
by one condition and one only: **full mutual independence**.

So this session built the other end of the measure route. `exists_iIndep_lonely`: if
the danger sets are mutually independent as events, the safe set is their complements'
intersection, which by independence has measure `∏(1 − 1/7) = (6/7)ⁿ > 0`, for *any*
`n`, the fourteen-runner case included. No window, no citation, no cap at seven. The
proof is three mathlib lemmas — independence passes to the generated σ-algebras, the
complements live in those σ-algebras, and independent intersections multiply — and it
is kernel-pure.

## The obstruction was always the correlation

Put the two ends together and the synthesis dissolves. If the exactness you wanted to
carry into the window actually holds — if the runners are independent — you never
needed the window: the full-circle closer takes any `n`. If it does not hold, the
window cannot manufacture it: a boundary does not create independence. So the entire
content of "carry star into the window" reduces to a single hypothesis, and that
hypothesis is now named in a theorem statement: `iIndepSet`. The families where it
holds are closed outright. The families where it fails are the correlated ones — and
"correlated integer speeds" is precisely "near-equal": two speeds `w` and `w + d` with
`d` small share almost all their teeth, so their danger events are almost the same
event, the opposite of independent. The tight-small corner of S25, the consecutive
integers `23…29`, is not a hard instance of a measure problem. It is the exact locus
where the measure problem's one hypothesis is false.

This is worth stating plainly because it retires a hope the project carried for several
sessions — that the right measure lemma, cleverly windowed, would close LRC(14). It
will not, and the reason is not a missing lemma. The measure route is now *complete*:
`star_safe` for the pairwise-center independence a covering family always has (good to
seven), and `exists_iIndep_lonely` for the full independence that closes any number of
runners. Everything a measure argument can prove about the critical band, these two
prove. What is left over is not measure-theoretic at all. It is the arithmetic of when
`{ (w+j)t }` — a short arithmetic progression on the circle — can be forced off a fixed
arc, and arithmetic progressions are correlated by construction. The fourteen-runner
conjecture, in the corner that remains, is a question about correlation, and no amount
of exactness on the whole circle speaks to it. That is not a defeat; it is a map. It
says the remaining work lives in the three-distance theorem and the resonant bounded
combo, not in the measure of arcs — and it says the measure of arcs has now been
spent to the last drop.

## What was actually built

Two theorems and an instance: `exists_iIndep_lonely` and its real-time form, closing
LRC at `1/14` for mutually independent families of any size, and the
`IsProbabilityMeasure` instance on the unit circle that lets them speak. Together with
S26's `star_safe`, the file `LRCStarSafe.lean` is the complete measure-theoretic
account of the seven-wall: what independence buys, exactly, at every value of `n`. The
conjecture is not closed. But the boundary between what the measure route reaches and
what it cannot is now a theorem, not a hunch, and it falls exactly where the arithmetic
begins.

---
*Linked: [[the-full-circle-is-exact]] (S26, the exactness this note explains you cannot
window), [[the-density-discreteness-transition]] (S25, the `wL~1` corner now seen as
the failure of independence), MISTAKE-072, HYP-3982, HYP-3983, klein HYP-4021.*
