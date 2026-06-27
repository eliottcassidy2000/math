# LRC(14) in three engines: lift, decorrelate, lift-to-H=7

*kind-pasteur-2026-06-27-S253. Integrating the current proof DAG (THM-523…572, the labelled-packet
families) with the witness-route floor work, the whole proof collapses to three engines — and two
problems that were being attacked separately turn out to be one.*

## The proof as it stands

Every primitive 13-set `S` is sorted by its **divisibility threshold** `q(S) = min{d : S has no
multiple of d}`:

- **`q(S) ≤ 14` — non-covering.** THM-523's `τ = 1/q` witness gives `M(S) ≥ 1/q ≥ 1/14`. The tight
  boundary (`M = 1/14`) is the census `{AP, GW}` — bounded + single-swap *proven* (HYP-2920/2921);
  the full census is the three-gap / consec-maximizes rigidity (open, but characterized as a binding
  scale on the tournament-spectrum Farey tree).
- **`q(S) > 14` — covering.** Contains a multiple of every `q ∈ {2,…,14}`, so `r := |14ℤ ∩ S| ≥ 1`.
  This is the hard core, and the literature/repo result for `≤13` runners (LRC≤13) is the accepted
  input that feeds it.

That sorting is settled. The action is all in the covering branch, and the current program splits it
by `r` and by a zoo of labelled-packet families (L1–L5, F0–F7). The synthesis below says that zoo is
three engines.

## Engine 1 — q-witness + census (the non-covering world)

`τ = 1/q` and the census. Proved except the tight-locus completeness, which is its own (genuinely
open, literature-hard) rigidity statement and is *isolated* from the rest of the proof: tight sets are
non-covering, so they never collide with the covering hard core.

## Engine 2 — lift-and-decorrelate (the covering world, all `r`)

Here is the unification. Write a covering set as `S = R ∪ 14Q`, `R` the 14-free part (`13−r` speeds),
`14Q` the `r` multiples. Substitute `u = 14t`: a multiple `14m` is safe iff `‖m·u‖ ≥ 1/14`, so **the
multiples of 14 become a sub-lonely-runner instance for `Q` (`r` runners), lonely by LRC≤13**. The
14-free `R` plays the witness "small part." Loneliness of `S` is exactly

> `R`-safe `∩` `Q`-lonely `≠ ∅`  (the **decorrelation floor**).

And this single engine covers all `r`:

- **`r ≥ 7`** is THM-571's gamma descent. With `≥9` speeds divisible by 7, the 7-lift is *explicit* —
  the decorrelation has a clean comb structure and `≥3` of 7 lifts survive. Proved (modulo LRC≤13).
- **`1 ≤ r ≤ 6`** is the open "few-apex" branch (HYP-2968). I verified this session that it is the
  *same* decorrelation floor — and, crucially, that **Bonferroni fails here**: for `{1..11,13,84}`
  (`r=1`) and `{1..10,13,84,154}` (`r=2`), `meas(R-safe) + meas(Q-lonely) − 1 = −0.13, −0.18 < 0`.
  Loneliness survives only because the two events are *quasi-independent*: `R' = meas(both) /
  (meas R · meas Q) = 0.50, 0.93 > 0`. That `R' ≥ c > 0` is exactly the witness-route Node-3 floor —
  the L2-Cauchy–Schwarz spectrum bound (HYP-2861), the `3/π²` Farey floor (HYP-2856), the
  complement-even handle (HYP-2871).

So the two things the program attacks separately — **THM-571 (`r≥7`, the 7-lift)** and **HYP-2968
(`r≤6`, lift-packet positivity)** — are one engine: the `u=14t` lift plus the decorrelation floor
`R' ≥ c`. My Node-3 floor, built for the general small-part/cluster decorrelation, is precisely the
instrument the open `r≤6` branch needs, in the exact regime (Bonferroni-negative) where nothing
cheaper works. **Prove the uniform `R' ≥ c` and both the few-apex covering branch and the general
witness floor close together.** Two open problems become one.

Better still, the few-apex case is the *favorable* end of that one problem: `Q` has `≤6` runners, so
`Q`-lonely is a large chunk (`M(Q) ≥ 1/(r+1) ≥ 1/7`); the difficulty is concentrated in `R`'s
14-free, dense arc structure — which is bounded-spread-friendly, the regime where my floor is already
rigorous.

## Engine 3 — the H=7 lift (the residual hard kernel)

What the decorrelation floor cannot reach is the zero-open-mass kernel (the K33 / L4 families): atoms
with no positive Haar-open witness at all. For those, THM-572 supplies the complementary closure — a
"bad atom" that lifts to a tournament with Rédei count `H = 7` is impossible, because `H = 7` is
forbidden for *every* tournament (THM-029 / THM-343). This is the project's home-turf result doing
real work on the conjecture: where analysis (the floor) sees zero, combinatorics (the H-spectrum gap)
sees a contradiction. The open obligation here is *constructing* the state-lift (HYP-2908), not the
closure itself.

## The improved argument, in one line

> **LRC(14) = [q-witness + census] ⊕ [`u=14t` lift + decorrelation floor `R' ≥ c`, all `r`] ⊕
> [H=7 tournament-state lift for the zero-mass kernel].** The middle engine subsumes THM-571 and the
> open few-apex branch into a single floor obligation — the uniform `R' ≥ c` — and that obligation is
> the witness-route work, in its most favorable (small-cluster, Bonferroni-negative) regime.

## What this buys, honestly

It does not prove LRC(14). The uniform `R' ≥ c` (the wide-`V` floor) is still open, and so is the
state-lift construction and the census completeness. But it **collapses the residual**: the
labelled-packet zoo and the `r`-split are not independent obstacles — the soft families are one
decorrelation floor, the hard kernel is one H=7 lift, and the tight boundary is one census. Three
obligations, not a dozen. And it tells the next session where to push: **the uniform decorrelation
floor `R' ≥ c`, attacked first on the few-apex (`r≤6`, small-`Q`) covering sets**, because closing it
there closes the open covering branch *and* advances the general witness floor — the highest-leverage
single move on the board.

— Related: [[lrc14-thread]], HYP-3121 (this synthesis), THM-571 (gamma descent), THM-572 (H=7 lift),
HYP-2968 (few-apex), HYP-2861/2916/2856/2871 (the Node-3 floor), HYP-2920/2921 (census), THM-523;
`the-tournament-spectrum-is-the-object.md`, `the-lonely-runner-is-a-random-round-tournament.md`.
