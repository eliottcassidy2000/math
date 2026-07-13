# The arc-count bound is false — far-element isolation, not diameter, is the classifier

*klein-2026-07-13-S289. Owner: prove the arc-count bound `r < 3√2·v|G'_{~v}|`. I could not — because it
is false. A census of covering 13-sets produces 938 failures at `max≤18` (ratios to 8.4), and a single
clean counterexample settles it: `{1,90,91,…,101}` is a covering set of large diameter that fails the
certificate at every peel. This is an honest negative that corrects my S288 framing: the crude `disc_v`
bound does not reduce the covering case to a combinatorial arc-count lemma. It closes exactly the
isolated-far-element covering sets — which include the covering-min extremals — and the rest reduce to the
same endpoint cancellation the density route needs.*

---

## The refutation

S288 proved `disc_v ≤ r²/(3v²)` and reduced the certificate to `r < 3√2·v|G'_{~v}|`, calling the residue
tight (ratio 0.92) and hoping the condition held combinatorially. It does not. Exhaustive census of
covering 13-subsets of `{1..N}`: **938 failures for `max≤18`**, ratios up to **8.4** — not marginal. The
decisive counterexample is structural, not a boundary case:

> `{1, 90, 91, …, 101}` — covering, diameter 100 — has `r=132` arcs and **fails at every peel** (best
> ratio 3.57). `{1,2,3,50..59}` fails at 4.48.

So the bound is not a theorem, and no combinatorial argument will make it one.

## The real classifier: isolation, not diameter or max

I first guessed the failures were "small max" (finite, checkable) or "bounded diameter" (handled by the
existing finite check). Both are wrong. `{1,90..101}` has large max and large diameter yet fails. The
actual dividing line is **whether the far element is isolated**:

- **Isolated far element** (`v ≫ max(W)`): then `r ≤ 2max(W) ≪ v`, so `r < 3√2 v|G'|` holds — the
  certificate fires. The covering-min extremals live here: the deep well (`v=182`, `max(W)=12`) and the
  residue (`v=84`, `max(W)=13`) are *forced* to have an isolated far element by the divisibility structure
  (`{1..12}` misses `13,14` ⟹ a far element divisible by `182`). That is *why* they are the covering-min,
  and why the crude bound reaches them.
- **Non-isolated far element** (`v ≈ max(W)`: compact sets, or large-diameter sets with a clustered top
  block like `{1,90..101}`): the top block is high-frequency, so `r` is large (`132`!) while `v` is only
  moderate. Ratio `≫ 1`. These are an **infinite** family — not bounded-diameter, not a finite check.

## Why this is not a detour but the same wall

The crude bound discards the endpoint cancellation via `|U(ℓ)| ≤ 2r`. For an isolated far element that
loss is affordable (`v` is huge). For a clustered top block it is fatal: `r=132`, but the *true* `disc_v`
is far smaller because `|U(mv)| ≪ 2r` (the `132` endpoints interfere). That interference is **exactly** the
endpoint cancellation the density route needs for `Q_s` (THM-729). So S288's "reduces to a combinatorial
arc-count condition" is only true for the isolated-far class; for everything else it reduces to the *same*
shared cancellation. There is no combinatorial shortcut around the analysis — consistent with the
S285–S287 picture of one hard cancellation underlying both routes.

## What stands

- **THM-732 is correct and its certification of the covering-min extremals stands** — the deep well and
  residue, the proven binding families, are rigorously certified `L>0`. That was the eighty-session
  sticking point and it is genuinely closed by elementary means.
- **The crude route closes all isolated-far-element covering sets.**
- **`r ≤ 2max(W)`** is a clean verified regularity (0/966 violations) — the honest form of "arc count is
  controlled," though it does not rescue the certificate for non-isolated sets.
- **The general covering case still needs the shared endpoint cancellation** (`|U(mv)|≪2r` = density
  `Q_s`, THM-729) for the non-isolated families.

The lesson: when a crude bound suddenly closes the hardest known family (the residue), check whether it
closes the *easy* families too before declaring victory. Here it does not — the residue is easy for the
crude bound (isolated far element) precisely because the divisibility structure that makes it extremal
also isolates its far element. The genuinely hard families for the crude bound are elsewhere, and they
sit on the same cancellation as everything else.

*Files: `04-computation/lrc14_arccount_census_klein_S289.py`, `lrc14_twolarge_check_klein_S289.py`,
`lrc14_clustered_top_klein_S289.py` (+.out). HYP-6505. Corrects
[[the-disc-v-bound-is-arc-counting-not-analysis-and-the-good-sets-have-few-arcs-klein-S288]]; the residual
is the shared cancellation of [[the-covering-middle-order-x-integral-is-a-good-set-autocorrelation-discrepancy-and-it-certifies-klein-S287]] (THM-729).*
