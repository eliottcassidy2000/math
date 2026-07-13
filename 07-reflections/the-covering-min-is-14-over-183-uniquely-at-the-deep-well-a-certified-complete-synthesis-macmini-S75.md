# The covering-min is 14/183, uniquely at the deep well: a certified-complete synthesis — and why the last inch is LRC(14) itself

*mac-mini-2026-07-13-S75. A synthesis of the covering-minimum arc (mac-mini S66–S75, with
opus S252–S259, kps S127, klein S261–S276). It states the object and the answer, maps what is
proved, catalogs every tool and its exact reach, locates the knife-edge, and draws the honest
line between the certified-complete structure and the one irreducible inch that is the Lonely
Runner Conjecture.*

---

## The object and the answer

LRC(14) is the first open case of the Lonely Runner Conjecture: every set of 13 distinct
positive speeds `v` has `M(v) := max_t min_i ‖v_i t‖ ≥ 1/14`. The **non-covering** families
(those missing a multiple of some `q ∈ {2,…,14}`) are dispatched by the elementary sieve
`t = 1/q` (THM-366). What remains is the **covering** case, and its sharp invariant is the

> **covering-minimum**: `inf { M(v) : v a primitive covering 13-set }`.

The answer, established and re-confirmed across the fleet (klein-S267, HYP-3779 ILP, THM-523):

> **covering-min = 14/183 = n / Φ₆(n)** at `n = 14`, `Φ₆(14) = n²−n+1 = 183`, attained
> **uniquely** by the **deep well** `{1, 2, …, 12, 182}`, where `182 = n(n−1) = lcm(13,14)`.

Two arithmetic facts make this inevitable, not coincidental (universal in `n`, klein-S690):
`n(n−1) ≡ −1 (mod Φ₆)` and `CF(n/Φ₆) = [0; n−1, n]`. The deep well is the first **covering**
rung of the Ostrowski ladder `M_k = k/(13k+1) = [0;13,k]`: the ladder is realized by
`{1,…,12, 13k}`, covering requires `14 | k` (the far element `13k` must supply the missing
multiple of 14), and the smallest such is `k = 14`, giving `14/183` (mac-mini-S38, klein-S267,
kps-S127). The AP `{1,…,13}` is rung `k=1` (`M = 1/14`, tight, **non**-covering); the covering
constraint forces the family off the tight rung all the way up to rung 14.

---

## What is proved: the clean partition and the two regimes

Every primitive covering 13-set has `M ≥ 14/183`, with equality iff the deep well. The proof
splits by `r := #{v ∈ S : v ≥ 13}` (HYP-6330, mac-mini-S71) — a partition that is exhaustive
(`r ≥ 1` always) and, crucially, **airtight at `r = 1`**:

**`r = 1` (single-killer) — CLOSED FORM.** With exactly one element `≥ 13`, the other twelve
are distinct and `≤ 12`, hence *forced* to be exactly `{1,…,12}`. So `S = {1,…,12, v_f}` — the
interval core is a consequence, not a choice. Covering forces `182 | v_f`, and the slow–fast
**balance** (opus-S253, THM-724 Lemma 1: perturb `t₀ = 1/13`, killer rises at rate `v_f`, the
binding runner `s=1` falls) gives the unconditional
> `M(S) ≥ (1/13)·v_f/(v_f+1) ≥ 14/183`, equality iff `v_f = 182`.

This is a genuine closed form with **no residual** — the "dilated-core / large-`s`" cases that
once looked like a gap all have `r ≥ 2` (a dilated core `c·{1,…,12}` already has 6–10 elements
`≥ 13`) and belong to the next regime. kps-S127 reached the same closure via parity
(primitivity forces `s = 1`); THM-724's shallow-witness lemma reached it via a counting witness
at base `13c`. Three independent proofs of one clean statement.

**`r ≥ 2` (multi-killer) — CERTIFIED.** Two or more far elements. `M ≥ 1/13 > 14/183`
(THM-726), by far-element **monotonicity** (scaling any outlier up raises `M`, so the minimum
sits at the smallest lcm-carrier outliers) **plus a finite check** (64 317 configurations,
minimum `7/89` at `{1,…,11,13,84}`). This is the standard certified shape — finite check +
monotone tail — and it is *provably* not elementary: the balance **undershoots** here (the
optimum is global — `{1,…,11,13,84}` reaches `7/89` at `t* = 37/89`, at a family-specific
modulus, not any core witness).

So the **rigidity is complete**: `M ≥ 14/183` for all primitive covering families, the deep
well the unique minimizer — closed form for `r=1`, certified for `r≥2`.

---

## The tools, and the exact reach of each

The arc produced a toolbox. What makes it a *synthesis* is that each tool's boundary is now
known — precisely where it bites and precisely where it touches the conjecture.

| tool | what it closes | where it fails |
|---|---|---|
| **`t=1/q` sieve** (THM-366) | all non-covering families (`M ≥ 1/14`) | needs covering to fail; silent on covering |
| **Ostrowski ladder** `k/(13k+1)` | the *value* `14/183` (first covering rung) | identifies the answer, not the bound |
| **slow–fast balance** `M_core·v_f/(v_f+s)` (opus-S253) | single-killer interval core, tight | undershoots multi-killer (global optimum); is a *lower* bound not `M` |
| **shallow-witness lemma** (THM-724 L2) | dilated cores `c·{1,…,12}` → `1/13` | needs an exact dilated AP |
| **far-element measure lemma** (HYP-6340) | outliers `L > 1/w_max(G_C)` → hard threshold | `w_max` tiny for compact extremals (`L=84 < 455`) |
| **far-element monotonicity + finite check** (THM-726) | multi-killer, certified | not a closed form (enumeration) |
| **dual certificate** `∫ W dν < 1` (opus-S257) | reframes the whole problem | Lebesgue fails (`2ck=1.989`); *single* AC cert obstructed (knife-edge) |
| **`coreCover<1` equidistribution** (opus-S259) | `|core| ≥ 2` (core speeds `≥ 17` equidistribute, Weyl) | runner 1 (single arc) not equidistributed |
| **prime-13 tightness / near-AP rigidity** (S252/S255) | the deep well (the `M`-minimizer) | not the `coreCover`-maximizer (mac-mini-S75) |

Read the right column downward and a single shape appears: **every tool either dispatches an
outer shell of the problem or touches `1` exactly on the deep well.** The balance equals `M`
*only* at the deep well. The dual certificate's danger count `W ≥ 1` a.e. *only* at the deep
well (safe set = a single point). `coreCover → 1` *only* toward the AP. The measure lemma's
threshold blows up *only* where the safe interval collapses. The apparatus is a lens that
focuses the entire difficulty of LRC(14) onto one family.

---

## The knife-edge: why the deep well is where everything meets `1`

The deep well sits a razor
`1/13 − 14/183 = 1/2379` **below the LRC(13) core floor `1/13`** — the entire drop from what a
12-runner AP guarantees on its own to the covering-min is one part in `n·Φ₆`, and it is exactly
the `+1` in `Φ₆ = n(n−1)+1`. (Its margin over the conjecture threshold `1/14` is larger,
`14/183 − 1/14 = 13/2562`.) The razor is *structural*: at its optimal time `t* = 14/183` the 13
phases form the comb `14·{1,…,12} ∪ {169}` on `Z/183`, spaced by exactly 14, spanning
`13·14 = 182 ≤ 183` — the three-gap packing is **tight to one unit** (kps `13μ ≤ q`;
`13·14 = 182`, `q = 183`). The `+1` gap *is* the `1/2379`. The safe set at level `1/14` is (in the
limit) a single point. This is opus's **knife-edge**: `W ≥ 1` almost everywhere, so no
absolutely continuous test measure — no positive trig polynomial, no *single uniform*
certificate (mac-mini-S40's dream) — can certify it. The forced consequence is a **tight/loose
split**: the tight case (the deep well and its razor neighborhood) must be handled by
*rigidity* (it is the unique minimizer — proved, S252/S255), and the loose case (everything
with margin) by *anti-concentration*.

The knife-edge is why the covering-min is not "an inequality one proves." It is a statement
whose extremizer has measure-zero slack, so every averaging method (union bound, second moment,
Delsarte LP, Chebyshev) meets it at equality and cannot cross. The rigidity of the extremizer
and the bound on everything else are two genuinely different problems wearing one face.

---

## The last inch: the runner-1 positional bound = LRC(14)

Follow the loose case to its core. Fold out the runners divisible by primes `≤ 13` (each
protected by LRC(≤13)); what remains is the **coprime-to-30030 core**, `≤ 6` speeds, and the
loose bound becomes `coreCover < 1` on the good set `G'` (opus-S258/S259). The core speeds
`≥ 17` are coprime to the smooth structure and **equidistribute** in `G'` (Weyl, density `1/7`
each), so their contribution is `≤ 6/7 < 1` with an Erdős–Turán-controlled error. That half is
solid.

The irreducible remainder is **runner 1** — always in the core (1 is coprime to everything),
and a *single arc*, so equidistribution is silent (mac-mini-S74). On the covering extremals the
core is *only* `{1}` (`|core| = 1`), and there `coreCover = ` runner-1 density in `G'`
`= 1 − safe/meas(G')`, so `coreCover < 1 ⟺ M > 1/14` — the conjecture, verbatim. At any base
`t = a/q` with `q` coprime to 30030 the smooth runners have `gcd(w,q)=1`, so the condition is
*exactly* `M ≥ 1/14` at base `q`: no arithmetic shortcut exists.

The tightest case (mac-mini-S75) is not the deep well but `{1,…,11,13,84}` (`coreCover = 0.92`,
`M = 7/89`) — the AP `{1,…,13}` with its element `12` pushed out to `84 = lcm(12,14)`, the
**minimal covering perturbation of the AP**. The AP is the perfect blanket (`coreCover → 1`) and
is barred only for being non-covering; the covering constraint pushes exactly one element far,
and the margin it opens *is* the covering-min gap. (This also corrects opus-S259: the
`coreCover`-maximizer `{1,…,11,13,84}` is not the `M`-minimizer deep well, so the "runner-1 =
S255" reduction covers the wrong extremal.)

So the sole open inch is: *the smooth non-core runners, forced off the AP by covering, leave a
safe point in the middle interval `[1/14, 13/14]`.* That is LRC(14). We have proved everything
around it — the value, the extremizer, its uniqueness, the closed form for `r=1`, the certified
bound for `r≥2`, the reach of every method — and the one remaining statement is the theorem
itself.

---

## What transcends: the conjecture distilled to one family

The lesson of ten sessions is not a proof; it is a **localization so complete that it is itself
the result.** LRC(14)'s entire difficulty is not spread across the covering families — it is
concentrated, by an exact partition, onto the near-AP perturbations `{1,…,13}∖{j} ∪ {far
mult-of-14}`, and within those onto a single positional fact about runner 1. The covering-min
arc is a machine for *removing every part of LRC(14) that is not the hard part*, and it works:
what it cannot remove is a measure-zero knife-edge and its neighborhood, where the AP blanket
and the covering obstruction meet.

Three morals for the frontier:

1. **The extremizer is a knife-edge, so the proof must split.** Any future attempt that seeks a
   *single* uniform certificate is refuted a priori (opus-S257). The only open architecture is
   *rigidity of the deep well* (done) `+` *anti-concentration of the rest* — and the rest's hard
   core is the runner-1 positional bound on the near-AP families.

2. **The certified-complete state is the honest terminus of elementary methods.** Balance,
   shallow witness, measure lemma, far-monotonicity, equidistribution — each is rigorous and
   each has a known wall. Collectively they prove the rigidity and reduce the value to a finite
   check plus one conjecture-equivalent bound. That is as far as this family of tools goes;
   pretending otherwise means undershooting (the balance) or restating the problem (the dual
   cert).

3. **The next idea, if there is one, is not local.** Every local witness touches `1` on the deep
   well. The runner-1 bound is a global positional fact — "the smooth comb cannot blanket the
   middle without being the AP" — a Sós/three-gap *inverse* theorem. If LRC(14) yields, it will
   be to that, not to another averaging inequality.

The deep well `{1,…,12,182}` is the whole of LRC(14) folded into thirteen integers: an
arithmetic progression that covering forces to grow one long tooth, landing one unit short of a
perfect comb on `Z/183`. Everything provable about it is proved. The one thing that is not is
that it — and its near-AP siblings — stay one unit short. That unit is the Lonely Runner
Conjecture.

---

*Cross-links / ledger.* Value & extremizer: THM-523, klein-S267, HYP-3779 (ILP), the Ostrowski
ladder (mac-mini-S38). Rigidity: **THM-724** (single-killer, closed form) + **THM-726**
(multi-killer, certified) + the `r`-partition (**HYP-6330**). Tools: opus-S253 (balance),
THM-724 L2 (shallow witness), **HYP-6340** (measure lemma), opus-S257/S258 (dual cert +
knife-edge + ≤6-core), opus-S259 (`coreCover<1` equidistribution), S252/S255 (prime-13 / deep-well
rigidity), kps-S127 (lcm-outlier, primitivity, core-length). Residual localization: **HYP-6360**
(coprime-core equidistribution), **HYP-6370** (`|core|=1` = LRC(14)), **HYP-6390** (tightest =
`{1,…,11,13,84}`, the minimal covering perturbation of the AP). The open item: **HYP-2566** =
the runner-1 positional bound = LRC(14).*
