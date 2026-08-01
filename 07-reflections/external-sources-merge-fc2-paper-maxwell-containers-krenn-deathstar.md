---
source: death-star-2026-08-01-coinC2
status: >
  SYNTHESIS + external-source integration. Four owner-supplied sources merged against the
  repo's open lanes. The FC one is decisive for our FC lane: it CORROBORATES my THM-3018
  section 4b audit flag and CONFIRMS THM-3031's identification of transcendence (not
  nonvanishing) as the operative hypothesis. Provenance NOT verified by this repo: the FC
  item is a self-described research draft, not peer-reviewed.
tags: [factorial-conjecture, transcendence, E-functions, pakovich-muzychuk, lrc14, extremality,
       hypergraph-containers, krenn, external-merge]
related: [THM-3018, THM-3031, THM-3012, THM-1123, THM-1134, kps-S166, HYP-9078]
---

# Merging four external sources into the repo's lanes

Owner-supplied, 2026-08-01. Each entry states what it is, what it *changes* here, and what
is merely adjacent. Adjacent-but-unusable is marked as such rather than dressed up.

## 1. octonion/mathematics/fc -- a full two-variable Factorial Conjecture attack

**What it claims.** `For every f in C[x,y], if L(f^m) = 0 for all m >= 1 then f = 0`, with
`L(x^a y^b) = a! b!` -- i.e. **general, inhomogeneous** `FC(2)`. Eight steps: algebraic
specialisation to `Qbar`; radial-projective coordinates via Gamma integrals; Mittag-Leffler
generating functions; logarithmic symbols; **constellation monodromy (Pakovich-Muzychuk)**
to exclude nonconstant projective leading forms; **flat-top reduction** to
`f_D = a(x+y)^D`; an **E-function transcendence theorem**; final contradiction. Self-described
research draft, AI-assisted, not peer-reviewed. **This repo has not verified it.**

The transcendence theorem stated there is

```text
for nonconstant q in Qbar[x] and distinct algebraic A, B:
        int_A^B e^{q(x)} dx  is TRANSCENDENTAL (hence nonzero),
```

by `E`-function moment systems plus **Beukers' refinement of Siegel-Shidlovskii**.

### What it changes here

**(a) THM-3031 is confirmed, and its point is the operative one.** THM-3031 proved that a
`FC(2)` counterexample pins `int_0^1 e^Q` to the value **1**, so the *nonvanishing* form is
useless and the bridge needs `!= 1`; transcendence supplies it because `1` is rational. The
external theorem is stated as "transcendental (**hence** nonzero)" -- transcendence is the
load-bearing half there too. Our correction to the relayed "`!= 0`" phrasing was right, and
the parenthetical in the source is exactly where the relay lost the content.

**(b) It CORROBORATES the THM-3018 section 4b audit flag.** Section 4b tried to kill
nonconstant leading forms with a soft max-modulus/Laplace argument; I flagged it today as
`AUDIT-REQUIRED` ((G1) the type shortcut is circular; (G2) the interior-maximum oscillation
factor `e^{-R b^2/(4c)}` is unaccounted). The external program does the *same job* --
excluding nonconstant projective leading forms -- and needs **constellation monodromy with
inverse-branch asymptotics at infinity** to do it. If a two-line Laplace argument sufficed,
that machinery would be unnecessary. This is independent evidence that 4b's step is not
soft, and the flag should stand until (G1)/(G2) are closed.

**(c) A caution against the obvious repair.** The tempting fix is "cite Pakovich-Muzychuk's
solution of the polynomial moment problem". **That does not directly apply.** The classical
PMP asks for `int_a^b p^i q' dx = 0` for **all** `i >= 0`; our condition excludes `i = 0`,
and indeed `int_0^1 g^0 du = 1 != 0`. THM-3018 section 5 already recorded the related
obstruction (with `q(u) = u`, `deg W = 1` makes `W(0) = W(1)` impossible, so the Composition
Condition is unavailable). The external work uses a *positive-interval constellation* variant,
not the plain PMP theorem. Any repair of 4b must respect this gap.

**(d) The right division of labour, now visible.** THM-3018 reduces the **homogeneous** case
to the simplex moment problem. The external program targets the **general** case and reaches
the flat-top residual `f_D = a(x+y)^D`, which is precisely where the exponential integral
enters. Note that on the simplex `x + y = 1` a flat top restricts to a **constant**, killed
by `m = 1` -- so the flat-top residual is a genuinely *inhomogeneous* phenomenon and cannot
be seen from THM-3018's reduction at all. That explains why our lane and theirs did not meet.

**(e) A concrete generalisation target for us.** Their integral is over **general algebraic
endpoints** `int_A^B`, ours over `[0,1]`. THM-3031's value-`1` computation used
`int_0^1 g^0 du = 1`. Over `[A,B]` the forced level becomes `B - A`, an algebraic number, so
the minimal bridge generalises to "`int_A^B e^Q dx != B - A`". Worth stating properly; it is
the same rigidity phenomenon with a different normalisation, and it is what a
`[A,B]`-version of THM-3031 must say.

## 2. arXiv 2607.27197, Arathoon-Ball-Kvalheim: "The Maxwell Conjecture is False"

Five point charges whose electrostatic potential has **at least 24 non-degenerate critical
points**, against Maxwell's conjectured `(n-1)^2 = 16`. Explicit construction.

**Transfer, and it is a live one.** Two of our lanes currently rest on unproven *extremality*
guesses:

* LRC(14) `r = 6`: THM-1134 closed all 792 cores for the **consecutive** killer shape; the
  remaining gap is whether consecutive is **globally extremal** among killer shapes.
* AMM 12592: whether the `gamma = 3/5` / floor-profile constructions are extremal.

Maxwell's bound was natural, classical, believed, about **critical points of a potential**,
and **false** -- refuted by explicit construction rather than by theory. Our LRC safe-set
analysis is also a critical-point/level-set count for a sum of periodic potentials. The
methodological content: when an extremality conjecture is the last step, **hunt for the
counterexample with structured constructions before investing in a proof**. This is already
in the LRC agent's brief and is now reinforced by a same-month precedent.

## 3. arXiv 2607.25973, Or Zamir: "k-Coloring is Faster than Computing the Chromatic Number"

Randomized `(2 - eps_k)^n` for every fixed `k`; **hypergraph containers** plus an iterable
`(k+1)`-list-colouring to `k`-list-colouring reduction.

**Honest assessment.** The *result* is not applicable here. The **container method** is
plausibly relevant to the LRC covering case, which is literally a covering/independent-set
question ("can `r` far killers cover the safe set?"), and containers are the standard modern
tool for turning such counting statements into *structure* statements -- which is the shape of
the AP-uniqueness inverse theorem that both LRC routes funnel into. But this is a **tool
suggestion, not a result**; I have not checked that the LRC hypergraph satisfies any container
hypothesis (bounded codegree / supersaturation), and it would need that before the method
means anything. Recorded as a lead, deliberately not as a claim.

## 4. Mixon's blog: Krenn's graph-colouring conjecture (3000 EUR)

Weighted edge-coloured bipartite `G = (A,B,E,k,w)` from quantum optics (GHZ states).
`k_max(n) = sup |k(E)|` over monochromatic `G` with `|B| = n`, `deg(a) = 2` for all `a in A`.
**Krenn:** `k_max(n) = floor(n/2)` for even `n >= 4`. Lower bound proved via perfect-matching
decompositions (1-factorisations of `K_4`, cycles `C_n`); **upper bound open for all n >= 4**;
`k_max(2) = infinity`.

**Fit with the repo.** Genuinely new here -- no existing lane. Closest neighbours are the
tournament/design work (TournamentH7, Paley, A000568 enumeration) and klein's covering and
perfect-matching arguments. The structure -- `deg(a) = 2` forces `E` into a union of perfect
matchings, and monochromaticity is a constraint on how those matchings can be coloured -- is a
1-factorisation extremal problem, which is squarely in that toolkit. Prize-bearing, clean, and
self-contained; a reasonable target if someone wants a bounded new lane rather than another
increment on FC/AMM/LRC.

## What I would actually do with these

1. **FC lane.** Treat item 1(b) as corroboration and keep 4b flagged. State the `[A,B]`
   generalisation of THM-3031 (forced level `B - A`) -- cheap and it aligns our statement with
   the external theorem's shape. Do **not** attempt the PM-citation repair, per 1(c).
2. **LRC lane.** Prioritise the extremality *counterexample hunt* over the extremality proof,
   with Maxwell as the precedent.
3. **Containers.** Park as a lead with the hypothesis-check named, not as a plan.
4. **Krenn.** Flag to the group as an available bounded problem; do not open it unilaterally
   while FC/AMM/LRC have live frontiers.
