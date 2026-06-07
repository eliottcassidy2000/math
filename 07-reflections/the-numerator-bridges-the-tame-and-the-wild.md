# The numerator polynomial bridges the tame and the wild

*monad-explorer-2026-06-07 (deep-research / analytic lane, 10th session). Builds on THM-438
ADDENDUM-7 (column denominator `(1−x)^{2m−1}` PROVED) → ADDENDUM-8 (the binomial reframing).*

## Two collapses, perpendicular

The cycle-rank triangle `t(k,m)` of `(★★)` is wild in the raw. Its diagonal `t(m,m)=A088368(m)`
grows like `e·m!`; its unsigned row sums `1,4,23,160,1262,10944` are uncatalogued; the s-expansion
coefficients `R_s(m,e)` and the top residues `1,3,20,181` are in no OEIS sequence. Everywhere you
take an honest *count*, you get a factorial-scale object with no name.

And yet the table collapses to something tame along **both** of its axes, and along each axis the
collapse is an *alternating sum that annihilates the factorial*:

- Read along **rows** (alternate over cycle rank `m`): `Σ_m (−1)^m t(k,m) = (−1)^k C_k`. This is
  `(★★)` itself — the Catalan law. The factorial-scale unsigned row sum collapses to a Catalan
  number. ADDENDUM-5 showed this collapse is *genus-blind*: it lives ACROSS cycle rank, not inside
  the genus of any one graph.

- Read along **columns** (alternate over the number of lines `e`): `Σ_e (−1)^e R_s(m,e) = 0` for
  `m ≥ 2`. This is `deg P_m = m−2`. The factorial diagonal `R_s(m,m) = A088368(m) ~ e·m!` is one of
  the terms — and it *cancels*. The collapse here lives ACROSS line count, WITHIN a fixed cycle
  rank — exactly the "one level down" involution ADDENDUM-7 predicted would be where the real
  pairing first becomes visible.

So `(★★)`'s Catalan and the column's degree drop are the same phenomenon read in perpendicular
directions: **a factorial unsigned count, alternately summed along an axis, collapses to a tame
value.** The project's master slogan — *the Catalan is a cancellation, not a count* — was always a
statement about one axis. It holds on the other one too.

## The binomial reframing makes the axis visible

The thing that turned this from slogan into structure was noticing that `[x^k](x/(1−x))^e =
C(k−1,e−1)`, so the column GF is, coefficient by coefficient,

```
   t(k,m) = Σ_{e=m}^{2m−1} R_s(m,e) · C(k−1, e−1).
```

This is the column's Ehrhart/`h`-vector form. Two things fall out for free. First, the small-`k`
zeros (`t(k,m)=0` for `1≤k≤m−1`) need no argument at all — the binomials simply have no support
there. Second, the *only* surviving small value is at `k=0`, and `t(0,m) = −Σ_e(−1)^e R_s(m,e) =
−Q_m(−1)`. So the entire content of `deg P_m=m−2` is the single statement **`t(0,m)=0`**, equivalently
**`T_m(x)→0` as `x→∞`**. The handoff stops being an inequality about polynomial degrees and becomes
a vanishing at one point.

## Reading the numerator from infinity

Here is the image the reframing gives. The denominator `(1−x)^{2m−1}` is the Euler-characteristic
ceiling (ADDENDUM-7). What is left after the ceiling is subtracted — the numerator `P_m(x)`, degree
`m−2` — turns out to be a *bridge between a tame end and a wild end*, and you can watch it being built
by expanding the column at `x=∞` instead of at `x=0`.

At `x=0` (the bottom coefficient of `P_m`) sits `A088368(m) ~ e·m!` — the all-pairings overcount, the
wild factorial. At `x=∞` (the top coefficient, the lead) sits `2^m − 1` — Mersenne, tame, rational
generating function `y/((1−y)(1−2y))`, and morally the number of nonempty subsets of the `m`
independent cycles. The numerator polynomial `P_m` *interpolates from the wild factorial bottom to the
tame Mersenne top*, and its degree `m−2` is exactly the room a single column has to make that
crossing. Rationality of the column is the statement that the crossing is polynomial of the right
length — not too long (the Eulerian ceiling caps the denominator at `2m−1`), not too short (the
numerator must reach from `e·m!` down to `2^m−1`).

Stacking the columns, the two handoffs become the first two terms of one asymptotic series. With
`s = x/(1−x) = −1 − 1/x − 1/x² − …`, the section `x=∞` is `s=−1`, and

```
   [x^0]  U(x,y) = V(−1,y)        = −y                      (only cycle rank 1 survives)
   [x^{−1}] U(x,y) = −V_s(−1,y)   = −y/((1−y)(1−2y))         (the Mersenne leads)
```

`U(x,y)` is resurgent at `x=0` (factorial growth, no algebraic equation — ADDENDUM-6). But read from
the *other* end, at `x=∞`, its leading terms are rational in `y`. The wildness is all at the origin;
infinity is tame. The numerator is what you see when you stand at infinity and look back — and from
there the factorial has already cancelled itself away, leaving Mersenne numbers and a clean pole.

## Why this is the right kind of explanation

It is, once more, the project's recurring move in a new costume: an invariant that looks arithmetic
(`2^m−1`, `2m−1`, `A088368`) is really the shadow of a single combinatorial-topological discipline —
here, *one open walk* on a rank-`m` skeleton. The single-walk (Eulerian-trail) condition is what caps
the denominator at `2m−1` (ADDENDUM-7), what makes `S_k` a free *cumulant* rather than a moment
(ADDENDUM-4), and — conjecturally — what forces the numerator's tame-to-wild bridge to have exactly
degree `m−2`. The pole order, the Catalan law, and the degree drop are three readings of the same
"there is one walk, not a disjoint union."

The honest status is unchanged: the two collapses are *verified*, not *proved*. Each is generated by
an involution we have not yet written down — the `m`-shifting one for the rows (Catalan fixed points),
the `e`-shifting one for the columns (no fixed points above rank 1). The reframing's contribution is
to say precisely what those involutions must do and where to look: not in the genus of a graph, not in
an algebraic equation for `U`, but in the term-by-term sign pattern of a single counting axis, read
from infinity.
