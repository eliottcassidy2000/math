# Pollock's conjecture is the bounded-arity currency — and its proof method is the template for two of our open cores

**Source:** kind-pasteur-2026-06-13-S3. Dispatch: think deeply about a proof
method for Pollock's conjecture and apply it to our work.

## Where Pollock already lives in the repo

Pollock's conjectures (1850) — every integer is a sum of at most 5 tetrahedral,
or at most 7 octahedral, numbers — are not a new topic for us; they are the
**bounded-arity currency** of the additive-coverage frame the repo built in S501
(`goldbach-polygonal-zeckendorf-additive-bases`). That frame names three ways an
additive basis pays for full coverage:

1. **smoothing** (Goldbach: dense, noisy primes; Hardy–Littlewood singular series;
   binary open, ternary proved because the third summand buys an averaging
   dimension);
2. **bounded arity** (Fermat/Cauchy polygonal numbers: sparse atoms `~√N`, but `k`
   summands erase the residue gaps);
3. **normal form** (Zeckendorf: `O(log N)` Fibonacci atoms, coverage by recurrence,
   uniqueness by confluent carries = the repo's `I(P_m,1)` independent-set lens).

Pollock is squarely regime 2: tetrahedral/octahedral atoms are sparse, and the
whole content is whether a *bounded number of summands* suffices to hit every
integer. So a "proof method for Pollock" is a method for the bounded-arity
currency — and the repo already owns a working instance of exactly that method.

## THM-462 is a Pollock-method theorem in miniature

THM-462 proves the cubic channel `c3` (directed-3-cycle count) is gap-free —
*every* value in `[0,M(n)]` is realized. Read as additive coverage, that is
"`c3` is a complete basis." Its engine is the **literal prototype Pollock/Waring
theorem**: the score deviation is `f(s) = m(n) + (Σ e_i²)/2`, and **Lagrange's
four-square theorem** makes `Σ e_i²` reach every value in the top window, with
induction below. The arity is 4 (four squares); the coverage is total. So we have
been doing Pollock-method proofs without the name: a sum-of-squares
representation + the four-square completeness theorem + induction.

## The new datum: completeness degrades with cycle length (THM-498)

If `c3` is Pollock-complete, what about the higher channels? This session
computed (exhaustively, two independent counters, 0 disagreements):

- `c3` (degree 3): gap-free at all `n` — complete.
- `c5` (degree 5): gap-free at `n=5`, but **first gap at `n=6`, value 10**
  (`[0,12]∖{10}`, with 9, 11, 12 all realized). `c5` is *not* score-determined,
  so this escapes the four-square argument — and exactly there, a gap appears.
- `H = I(Ω,2)` (degree `n`): forbidden values 7, then 21.

So the OCF cycle channels form a **Pollock-completeness hierarchy**: the channel
is complete precisely while it is a sum-of-squares (degree 3), and additive
completeness fails sooner as the cycle length climbs. The forbidden values
(`c5=10`, `H∈{7,21}`) are the **Pollock exceptions** — the finite residue
obstructions a circle-method proof would isolate as "singular series vanishes on
a finite set." This reframes the repo's long-standing "why is H=7,21 forbidden?"
as one instance of a general Pollock-completeness curve of the OCF.

## The payoff: the LRC deficit *is* a circle-method singular series

The deepest transfer is back to LRC(14). The covering deficit `D(q,S)` (THM-497)
has precisely the **Hardy–Littlewood shape** that a Pollock/Waring circle-method
proof produces:

```
D(q,S)  =  MAIN TERM  q·(6/7)^13   +   CHARACTER-SUM corrections   +  finite exceptional/resource set
           (the "expected"            (the "singular series" /         (the configs that survive —
            number of witnesses)        minor-arc fluctuation)           the over-correlated regime)
```

This is the same architecture as "every large `n` is representable (main term
dominates) + a finite list of small exceptions (checked by computer)." The
Pollock proof method — *bound the minor arcs / the singular-series fluctuation so
the main term dominates for all `n` past an explicit bound, then finitely check
below it* — is exactly the missing template for the THM-497 D open core (the
deficit lower bound in the over-correlated regime) and the resource bound t-0124.
Concretely: treat the over-correlation as the singular-series factor and bound it
by **incomplete multiplicative-character (Weil) estimates**, the way Waring's
problem controls minor arcs — respecting the tool-domain boundary (THM-497 D2:
LRC is multiplicative-character, not additive/code). The honest obstruction, as
last session found, is that the fluctuation is *large* (over-correlation grows
faster than `√q` for structured configs) — but that is the same difficulty
Waring's problem confronts at small numbers of variables, and the circle method's
answer (more variables = more averaging; or a sharper exceptional-set argument)
is the menu to try.

## The grounded proof method (what the 2025 proof actually does)

The strongest recent result — Basak, Dong, Saettone & Zaharescu, *Representations
as Sums of Icosahedral and Dodecahedral Numbers: Proof of Pollock's Conjectures*,
IMRN 2025 (DOI 10.1093/imrn/rnaf180) — and its direct precedent (Brady, *Sums of
seven octahedral numbers*, JLMS 2016, arXiv:1509.04316) use a **two-engine
hybrid**, the modern template for Waring-type problems on integer-valued
polynomials:

1. **Linnik's quadratic-form reduction.** Pair the values: `f(q+x) + f(q−x)` for a
   cubic `f` collapses to `2f(q) + (6aq)·x²` plus lower order; six paired terms
   become `6aq·(x² + y² + z²)` — a **ternary sum of three squares** — and a 7th
   term tunes a congruence mod `6aq`. The problem becomes counting lattice points
   on an algebraic curve with congruence constraints, controlled by **Bombieri's
   bound on exponential sums along curves** + Cobeli–Zaharescu / Hasse–Weil point
   counting.
2. **The Hardy–Littlewood circle method with an EXPLICIT power-saving error**,
   giving an asymptotic representation count positive for all `N` above a
   computable threshold (≈ `9.6·10³⁵` icosahedral).
3. **A finite descent + computer check** below the threshold.

Two things this nails down for us. (i) **THM-462's `Σe_i²`/four-square engine is
literally Linnik's ternary-form trick in miniature** — the `f(q+x)+f(q−x)` pairing
that produces a sum of squares is exactly the score-deviation `f(s)=m+(Σe_i²)/2`.
So we were already running Engine 1. (ii) The refinement lesson: Pollock's
original bounds **13 and 21 are WRONG** — the true bounds are 15 and 22, off by a
finite exceptional set (`{47,83,94,95,119}` and `{79}`) that only a computation
reveals. This is the same lesson the repo keeps relearning (the LRC ceiling
`f(13)=41` was wrong; `H`'s forbidden `{7,21}`; `c5`'s forbidden `10`): the clean
conjectured bound is generically off by an exceptional set you must compute, never
guess. *Verified:* our additive-basis DP, with the atom generator swapped to 3D
figurate numbers, **independently reproduces the entire landscape** — Waring
`{23,239}`, tetrahedral's 241 exceptions (largest 343867), and the 2025 icosa/dodeca
corrections exactly (HYP-2490). Our machinery is the finite-check engine of the
Pollock template.

## The sharpened LRC transfer (the actual toolkit)

With the method grounded, the LRC transfer (HYP-2489) sharpens: the over-correlated
regime of the LRC deficit — last session's open core — is exactly a
**count of lattice points on a curve / a multiplicative-character sum with
congruence constraints**, which is what **Bombieri-along-curves + Hasse–Weil**
control in Engine 1. So the precise instruments to try on the LRC deficit are not
just "some Weil bound" but the Linnik/Cobeli–Zaharescu point-counting machinery the
2025 Pollock proof uses, plus the explicit-error circle method for the main term
and a finite resource check below the threshold — the full Pollock template, ported.

## The one-line synthesis

Pollock's conjecture is the bounded-arity face of additive coverage; its proof
method is "main term by averaging + finite exceptional set." The repo already
runs this method (THM-462, four squares → complete `c3`); the new finding is the
completeness *hierarchy* (THM-498: `c5` breaks at 10); and the prize is that the
LRC(14) deficit wears the same circle-method clothes, so Pollock's method is the
right template for the LRC open core — not the additive/code machinery the repo
usually reaches for.

Cross-links: THM-462 (the four-square engine), THM-498 (the cycle-spectrum
hierarchy), THM-497 (the LRC covering deficit / its open core),
[[goldbach-polygonal-zeckendorf-additive-bases-s501]] (the three currencies),
HYP-2487/2488/2489.
