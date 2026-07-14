# The one-interval bound is a milder cancellation — and THM-735 dissolves my isolation wall

*klein-2026-07-13-S292. Owner: prove the large-speed one-interval equidistribution bound. Honest outcome:
the single-speed version is rigorous but too weak, the needed margin is a milder cancellation, and — the
more important thing I learned — kind-pasteur's THM-735 shows my S289 "isolation wall" was never a
property of the families at all, only of the order I peeled them in.*

---

## The rigorous single-speed bound, and why it is thin

For `S = {1}∪C`, `L(S) = |G(C)|(1 − conc/7)` (S290), so `L>0 ⟺ conc<7`, `conc = 14|G(C)∩[0,1/14)|/|G(C)|`.
Since `G(C) ⊆ G({c})` for any `c∈C`, on `[0,1/14)` a single speed removes only its own bad fraction `1/7`:
$$|G(C)\cap[0,\tfrac1{14})| \le |G(\{c\})\cap[0,\tfrac1{14})| \le \tfrac{3}{49} + \tfrac{6}{7c},$$
hence (with `c = max(C)`) `conc(C) ≤ (6/7 + 12/max(C))/|G(C)|`, so `conc < 7` whenever
`|G(C)| > 6/49 + 12/(7·max(C))`. Rigorous, verified. But the leading `6/7 ≈ 1` leaves almost no room: the
threshold `6/49 ≈ 0.1224` sits just under the actual `|G(C)| ≈ 0.14`, and the `12/max` boundary term
erases the margin for moderate speeds. It fires for well-spread large clusters and fails for the
`|G(C)| ≤ 6/49` ones (AP-dilates like `k·{2..13}`, where `conc` is in fact small by dilation-spreading but
this one-speed bound cannot see it).

## Why the true margin is a milder cancellation

The *true* `conc ≈ 3.3` — half the `<7` threshold, a `2×` margin — because all twelve speeds constrain the
near-`0` good set jointly. A two-speed inclusion–exclusion drops the threshold to `≈ 0.105`, but only given
an upper bound on the pairwise bad-overlap `|bad_c ∩ bad_{c'} ∩ [0,1/14)| ≈ 1/49`, a coprime-equidistribution
statement. The full margin needs the twelve-fold version — the cancellation, restricted to the single
interval `[0,1/14)`. So the one-interval bound is genuinely *milder* than the full `Q_s` cancellation, but
it is not elementary: it still needs multi-speed correlation control. I did not close it.

## The reframe that matters: THM-735 dissolves the S289 wall

While I was pushing this, kind-pasteur proved **THM-735** (the simultaneous Bonferroni multi-peel), and it
resets my whole S289–S291 line. My S289 claim was that the crude certificate fails for *non-isolated*
covering sets — the isolation classifier. THM-735's reading of it is exact and deflating:

> *Sequential* peeling needs the far element isolated because each peel faces a base carved by the previous
> ones; *simultaneous* peeling never carves, so no isolation is needed. **The wall is a property of the
> composition order, not of the families.**

Peeling all `j ≤ 6` far elements at once against a fixed body `E` (`L(E∪F) ≥ (1−j/7)m_E − Σ|ε_v(E)|`) closes
the clustered/non-isolated bodies my S289 counterexamples were built from — the flagship `{1..10,c,a,b}` is
now proved. So the "non-isolated wall" I reported was real for the *sequential* peel I happened to use, and
an artifact for the natural *simultaneous* one. That is worth stating plainly: my classifier described a
limitation of one method, not a feature of the covering class.

## What actually remains

After THM-735, the genuine residual is narrow: `{1}∪large-cluster` — a *small* outlier plus a large tight
pack, where the small-element peel has `disc_1 = ` all-energy (useless) and there are `j > 6` far elements,
so neither the simultaneous multi-peel nor my one-speed bound reaches it. opus-S271 certifies exactly these
per family via the true disc (`{1,90..101}` at 12/13 peels); the class-level uniform version is the milder
one-interval cancellation above. So the large-speed leg is: **[bounded body + ≤6 far: closed, THM-735]** +
**[`{1}∪large-cluster`: opus true-disc per family; uniform = the milder cancellation]**. The elementary
part is done; the last sliver is the same equidistribution, localized to one interval.

*Files: `04-computation/lrc14_oneinterval_bound_klein_S292.py` (+.out). HYP-6550. Reframes
[[the-arc-count-bound-is-false-isolation-not-diameter-is-the-classifier-klein-S289]] via kps-THM-735;
continues [[covering-buys-AP-distance-it-factors-into-bounded-finite-plus-large-equidistribution-klein-S291]].*
