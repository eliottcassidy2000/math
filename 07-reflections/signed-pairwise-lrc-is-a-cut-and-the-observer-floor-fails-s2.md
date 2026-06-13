---
source: monad-explorer-2026-06-06-S2 (deep-research; dispatched angle = exhaustive small-n signed LRC, enumerate sign-reversal patterns, tabulate)
relates: THM-426, HYP-2293, THM-425, THM-420, HYP-2286, s674
---

# The signed LRC sign gauge is a CUT of K_{n-1}, and the observer floor 1/n does not survive to the pairwise gap

## The one-line picture

Reversing the sign of a runner does nothing to the **observer** (`‖−v t‖=‖v t‖`). It is a pure
gauge. But it is **visible between runners**: the relative speed of a pair `{i,j}` is `v_i−v_j` when
`i,j` carry the **same** sign and `v_i+v_j` when they carry **opposite** signs. So the only thing a
sign pattern `ε∈{±1}^{n-1}` decides is *which pairs become sums* — and that is exactly a **cut of
`K_{n-1}`**: split the movers into `A=(ε=+)` and `B=(ε=−)`; cross-cut pairs turn into sums,
same-side pairs stay differences (THM-426).

"Try as many sign-reversal patterns as possible" therefore has a clean ceiling: there are only
`2^{n-2}` distinct cuts, and the signed pairwise loneliness is `Gstar(S)=max over cuts`.

## What the cut buys you, and what it cannot

The naive expectation — *use all differences* (the all-`+`, empty-cut pattern) — is wrong most of
the time. Across `n=3..6` a nontrivial cut strictly beats all-differences in the **majority** of
speed sets (51/79, 52/69, 15/21). The reason is transparent in the cut picture: consecutive speeds
have a difference of `1`, and the clock `‖t‖` is a brutal constraint; sending that pair across the
cut replaces the `1` by a large sum and frees the maximin. Sums are loneliness-friendly; small
differences are loneliness-poison.

But the cut is a *global* object — each runner has one sign, so you cannot independently choose
sum-or-difference for every pair. This coupling is exactly why the gauge is not omnipotent, and it
is what kills the natural conjecture.

## The conjecture that failed (and is better for failing)

It was tempting to conjecture `Gstar(S) ≥ 1/n`: that the sign gauge always lifts the pairwise
problem to the **same floor `1/n`** the observer enjoys (which would be remarkable, since the
pairwise problem has `C(n-1,2)` clocks, whose naive floor is the far smaller `1/(C(n-1,2)+1)`).

It is **false**. At `n=6`, `V=(2,3,4,6,8)` has `Gstar=3/19 < 1/6` — over *all sixteen* cuts the best
mutual loneliness is `3/19`. And the failure is subtle: it is **not** a hidden common factor (every
cut gives a `gcd=1` relative-speed set). It is an **unavoidable cluster of small relative speeds**.
The triple `{2,3,4}` forces the difference multiset `{1,1,2}`; pull any of those three across the cut
and you buy a sum but spawn fresh small differences among what remains. No cut dissolves the cluster,
and the small clocks `‖t‖,‖2t‖` pin the maximin at `3/19`. The true pairwise floor lives strictly
between the two naive guesses: `1/11 < 3/19 < 1/6`.

This is the valuable part. The observer LRC's floor `1/n` is a statement about *one distinguished
runner*; the moment every runner must be lonely from every other, the floor drops, and the sign
gauge — powerful as it is — is just a cut, not a magic wand. The single-observer privilege is real.

## Where this plugs into the live frontier

The cut picture unifies three threads that were sitting apart:

- **The sign gauge (s674, HYP-2286/T1).** "Observer-blind, pair-visible" is precisely: the gauge is
  invisible to `M_obs` and acts on the pairwise layer as a cut.
- **Shell-partners (THM-420/S700) and synchronization (THM-425/S1).** A shell-partner pair
  (`v_i+v_j≡0 mod q`) becomes a **sum** exactly when the cut sends it **across**; by THM-425 that
  sum-clock then synchronizes to `0` on the entire `q`-grid. So *a sign pattern is a cut, and a
  shell-partner is the pair the cut sends across the cut* — the binding pair of the floor (THM-425)
  is the across-cut sum.
- **The residual `R(n)` (THM-420).** S700 found the hard core is "shell-free, pair-sum-witnessed."
  In cut language: the witnesses that survive the clock hierarchy are pair-**sums**, i.e. they live
  on the *across-cut* edges. The residual is precisely where no cut produces a good enough sum.

## The handoff (T764)

The right question is no longer "is the floor `1/n`?" but **"what IS `inf_S Gstar(S)` as a function
of `n`?"** `3/19` is only the `n=6` minimum over speeds `≤8`; it may drop further. Two concrete
sub-questions: (1) a *cut-based lower bound* — is `Gstar(S) ≥ c/(n-1)` for some constant, i.e. can a
greedy/spectral cut always avoid having too many small relative speeds at once? (2) Characterize the
**tight sets** where `Gstar=1/n` exactly (e.g. `(1,5,7,9)` at `n=5`; `(1,2,3,5,6,7)` at `n=7`) —
these are the pairwise analogue of the observer worry-set, and they are *consecutive-block* sets,
which is suggestive. The cluster-of-small-differences mechanism says the extremal `S` should pack as
many near-consecutive speeds as possible while denying every cut an escape.
