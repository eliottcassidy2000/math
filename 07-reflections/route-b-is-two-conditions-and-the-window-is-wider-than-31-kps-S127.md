# Route B clearing is two conditions, not one — and the window is wider than 31

*kind-pasteur-2026-07-11-S127 cont.46. Owner: "attack the residual 8.5% as a bounded fold-class covering on
the small coprime sub-family." Attacking it surfaced the exact shape of clearing — and a correction: the
fleet's recently-narrowed `[15,31]` window is a sampling artifact.*

---

## Clearing is exactly two conditions

Building on the coprime reduction (cont.45), the exact condition for a family to clear at modulus `q` on the
`{0,±1}` band (`q ∈ [15,28]`) via *any* multiplier is:

> **(a)** `q ∤ vᵢ` for every `i` — no runner is a multiple of `q` (such a runner sits at residue `0` for
> *every* `p`, so it is always in the danger arc and `q` can never clear); **and**
> **(b)** the coprime-to-`q` sub-family **misses a unit fold-class** (so some unit multiplier avoids all of
> `{±vᵢ⁻¹}`).

And this is **exact and complete**: over 4500 spread divisor-complete families, **0** clear by any mechanism
other than (a)∧(b) at some window `q` — non-unit multipliers never add clearing. So "spread DC clears on the
window" is precisely "some window `q` is divisibility-unblocked *and* its coprime sub-family misses a
fold-class." This is the composite-modulus completion of klein's THM-718 (which counted the prime case): the
whole of Route B is this two-condition covering statement.

## The residual clears by collisions

The residual 8.5% (families where the elementary count guarantee `#coprime ≤ φ(q)/2−1` fails everywhere) all
clear anyway — by the genuine fold-class **miss**: at the clearing `q` the coprime runners *collide*, a median
of **5** of them sharing fold-classes, so even with `≥ φ(q)/2` coprime runners they occupy fewer than
`φ(q)/2` classes and leave one empty. The residual is a real anti-concentration on the small coprime
sub-family — but a robust one (max-over-window missed classes is `≥ 2` for essentially all of it).

## The window is wider than 31 — a sampling artifact corrected

opus-S238's minimal covering set `{15,16,19,23,25,29,31}` and klein-S258's finish-map window `[15,29]`/`[15,31]`
are **too narrow**. They were read off samples where families rarely have small elements coinciding with
window moduli. But condition (a) is exactly a *divisibility-blocking*: a spread DC family whose elements have
many divisors in the window blocks most of it. Adversarially:

- `v = [23,26,29,30,31,40,42,44,48,50,51,54,57]` has runners `= 23,26,29,30,31`, blocking those moduli; first
  clears at **q = 33**.
- Hill-climbing to maximise the min-clearing modulus reaches **q = 37**, then **q = 39**
  (`[36,66,133,…,216]` blocks all of `[15,36]` except `31,37` and fails (b) at `31`).

So the true bounded window is **not `[15,31]`** — it reaches at least `39`, and the safe empirical bound is
`[15,43]` (my cont.34 `[8,43]`). This does **not** threaten the conclusion: the band-edge margin
`⌈q/14⌉/q > 1/14` holds at *every* non-14 `q` (`3/37, 3/39, 4/43` all exceed `1/14`), so clearing anywhere in
the wider window still gives `M > 1/14`. What changes is the **finite-check range** — the diameter-free
residue check must run over `[15,43]`, not `[15,31]`, or it will falsely report un-clearable families that in
fact clear at `q ∈ [33,39]`. The narrowing was optimistic; the honest window is the wider one.

## The sharpened picture of the hard core

The spread DC hard core now splits cleanly by *which* condition binds:

- **(a)-hard families**: block most small moduli by divisibility (elements rich in small-window divisors).
  They need the *upper* end of the window — a `q` coprime to all runners. This is a covering-by-divisors
  question, and it is why the window must be wide.
- **(b)-hard families**: divisibility-unblocked but the coprime sub-family covers all fold-classes. This is
  the genuine anti-concentration, now on the small (median-3) coprime sub-family.

The elementary guarantee kills everything except a thin `(b)`-hard residual, and the `(a)`-hard families are
handled by widening the window, not by anti-concentration at all. Two different difficulties wearing one name
— the recurring shape of this project.

## Honest scope

No proof of boundedness here — that the min-clearing modulus is bounded over *all* spread DC families is Route
B itself, and remains open. What this session fixes is the **exactness** of the clearing condition (two
conditions, captures all), the **collision** mechanism of the residual, and the **window width** (`[15,43]`,
not `[15,31]`). The next move on the residual is the `(b)`-hard anti-concentration on the coprime sub-family;
the `(a)`-hard part is now understood as divisibility-covering and just needs the correct window.

*Files: lrc14_residual_foldcover_kps_S127.py/out, lrc14_true_window_kps_S127.py/out. HYP-6105. Corrects the
`[15,31]` window (opus-S238/klein-S258); extends the coprime reduction
[[the-coprime-reduction-shrinks-route-b-to-a-handful-of-runners-kps-S127]]; composite companion to
klein-S259/THM-718.*
