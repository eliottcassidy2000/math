---
source: monad-researcher-2026-06-01-S1
status: REFLECTION — methodological correction + the n=4 picture
tags: [LRC, covering, fourier, mod-4-character, measure-gap, tight-set, methodology]
---

# The safe measure is trivially nonnegative; the tight set is the whole problem

Working the S526 covering+character methodology down to n=4, one identity reframes
everything:

```text
|SAFE| = ∫_0^1 ∏_i 1[||s_i t|| >= 1/n] dt .
```

The integrand is a product of `[0,1]`-valued indicators, so **`|SAFE| >= 0` for
every speed system and every n, with no work at all.** The roots-of-unity
expansion `|SAFE| = (1-2/n)^{n-1} + (resonance corrections)` is exact and pretty,
but the inequality "corrections `>= -(1-2/n)^{n-1}`" that it seems to invite is
*just* `|SAFE| >= 0` — already free, and **not** LRC.

The reason LRC is hard is now stark: the regular polygon / arithmetic progression
has `|SAFE|` **exactly 0**. A measure bound that is tight at 0 cannot certify a
nonempty *set*. So the entire conjecture lives in the measure-zero stratum:

```text
LRC(n)  <=>  every measure-zero (tight) speed system still has a CLOSED point
             t* with ||s_i t*|| >= 1/n for all i.
```

This is the precise form of the "set-vs-measure gap" oracle named in HYP-2039.
It is not a residual sitting next to the main argument — it *is* the argument.
Everything with positive measure is free (positive measure ⇒ nonempty); all the
mathematics is in the boundary of the cases where the open danger arcs manage to
cover the circle up to measure zero.

## What n=4 shows cleanly

At n=4 the harmonic substrate is the **odd Dirichlet character mod 4**:
`g_k = -chi4(k)/(πk)`, because `sin(πk/2) = chi4(k)` — even harmonics die. The
pairwise term collapses to a closed form (THM-393),

```text
|S_a ∩ S_b| = 1/4 + chi4(a')chi4(b')/(4 a' b'),
```

the exact mod-4 sibling of the n=3 Legendre formula. But the 3-term resonance
`R3` is a genuine triple character sum on a rank-2 lattice and does **not**
factor — the same wall as n=14, in miniature.

Yet the *tight stratum* is astonishingly simple at n=4, and that is what matters:

- `{1,2,3}` is the **unique** measure-zero primitive triple (checked over 135,739
  triples, speeds ≤ 100). At n=5,6 the AP is not unique; at n=4 it is.
- Off the AP there is a hard **measure gap**: `|SAFE| >= 1/28`, equality at
  `{1,6,7}`, with the slow family `(1,4k+2,4k+3)` filling `1/28,1/22,1/20,...`
  from above.
- The AP itself is rescued by one rational point, `t=1/4`.

So LRC(n=4) is not "bound a triple character sum"; it is "prove a `1/28` measure
gap, then write down `t=1/4`." The gap is a finite-flavored, self-contained
target. The lesson generalizes: at every n the productive question is not *how
negative can the resonance correction be* (answer: down to `-(1-2/n)^{n-1}`,
trivially, at the AP), but *which speed systems are tight, and does each carry a
rational boundary witness.* The character machinery is for **locating** the
tight stratum (and `{1,2,3}` is forced there by an even speed + the all-odd
`t=1/4` lemma), not for a measure inequality that is either trivial or false-as-
stated.

The mathematics keeps pointing at the same place: the regular polygon, where the
measure vanishes and only the boundary survives. That is where to push.
