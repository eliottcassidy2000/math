# The real dichotomy is bounded vs unbounded: "existence is combinatorial" was the BOUNDED case, and klein just closed it — the open core is ANALYTIC, and that is the hope

*kind-pasteur-2026-07-01. The owner: "value is modular, existence is combinatorial" (opus HYP-3779, after the modular bound on a(n) failed three ways) — is there a creative reframe, a 3rd category, or a strategy reconsideration that gives hope? Yes. Read against klein HYP-3779 (same day) the "combinatorial existence" splits, and the half that is actually open is not combinatorial at all — it is the analytic unbounded regime, which is exactly where the modular value's structure becomes a usable certificate and where my far-element machinery already lives.*

## The pivot: "existence combinatorial" is the BOUNDED case, now CLOSED

opus HYP-3779's law — value modular, existence combinatorial — is true, but it conflates two very different regimes of "existence." The covering-min a(n) question ("is there a beater below the construction?") splits by speed size:

- **Bounded** (all speeds `≤ n(n−1)`): a FINITE search. klein HYP-3779 (S61) just **closed it rigorously for n=14** — lazy-cut cutting-plane ILP proves covering-min `= 14/183 = n/Φ₆` up to speeds `n(n−1)=182`, no beater. This is the "combinatorial" half, and it is *done* (for n=14, the target).
- **Unbounded** (some speed `> n(n−1)`, arbitrarily large — primitive covering sets `{1,…,12, W}` with `W→∞` exist): NOT a finite search, NOT combinatorial. Large speeds **equidistribute**. This is the genuinely open half.

So the honest statement is not "existence is combinatorial" but:

> **existence-bounded is combinatorial and FINITE (klein closed it); existence-unbounded is ANALYTIC (decorrelation / equidistribution).** The real dichotomy is bounded vs unbounded, and the open core is the *analytic* side.

This is the reframe that gives hope: the open problem is not a hopeless irregular sequence a(n) — it is a decorrelation estimate on large speeds, a mature analytic object.

## Why the unbounded case is tractable (it is already my thread)

The far-element machinery from the previous sessions is *exactly* the unbounded regime. For a covering set `C ∪ {W}` with `W` large, the exact decomposition (validated, `lrc14_w_as_logbarrier_centralpath_kps.py`):

> `p0(C ∪ {W}) = Plat(C) + Δ_W`, `Δ_W = O(1/W)` (the Dedekind/period-max residual, THM-563), `Plat(C)` = the bounded-core value.

So the **unbounded single-far case reduces to [bounded core, closed by klein's ILP] + [a decorrelation correction Δ_W → 0]**. Multi-far (several large speeds) is the r-far ladder (the Dedekind ladder), each rung an `O(1/W)` decorrelation. The unbounded case is not a search over infinitely many configurations — it is one analytic limit (`W→∞`) plus a finite base (klein). The "combinatorial" fog was only the bounded base.

## The 3rd category: the ANALYTIC DUAL CERTIFICATE — and the modular value is its raw material

Beyond "modular (value)" and "combinatorial (config)" there is a third kind of object: the **certificate** — the dual/variational object that proves optimality without enumerating configs. The emerging repo candidates are all of this kind, and they should be read as one category:

| candidate | who | the certificate is… |
|---|---|---|
| Kershner / Fejes-Tóth flat-case | klein S57 | a *covering-optimality* theorem (hexagonal A₂ thinnest) |
| hyperbolicity cone / self-concordant barrier | mac-mini HYP-3780, opus HYP-3769 | a *convex-geometric* certificate (the barrier) |
| character-sum / Weil bound | opus HYP-3779 | an *analytic* bound on the residual |
| Cohn–Elkies / Viazovska magic function | (this reflection) | a *modular LP-dual* certificate |
| Riesz product / far-element decorrelation | my prior sessions | the *analytic dual* for the unbounded case |

The unifying hope is one principle, and it is Viazovska's: **the modular structure of the value transfers to the dual certificate.** In the sphere-packing proof, the optimal-density value and the magic function that certifies it are *both* built from the same modular forms; the magic function is manufactured from the value's modularity. Here: the value is `−1/12 = ζ(−1)`, level-14 Eisenstein `E₂`, the Fejér kernel `F₇` (last session: the "cyclotomic magic function"). Those are exactly the ingredients of a **level-14 magic function** — the LP-dual certificate that would bound `M(S) ≥ 1/14` for the unbounded case, the same way the Fejér/Eisenstein forms give the value.

And the certificate's *shape* is already visible: the AP construction **equioscillates at the `φ(14)=6` units** (the Chebyshev/Fejér signature, prior reflection), which is precisely the complementary-slackness support of an LP-dual — the magic function is pinned at the 6 units, the 3 antipodal binding pairs. The dual is not hypothetical; its binding set is the equioscillation set.

## The strategy reconsideration

- **Old strategy (primal):** characterize a(n) / hunt beaters. This is the *bounded* problem — and klein's ILP just showed it is finite and closeable (done for n=14). Continuing to chase a(n)'s irregular sequence is chasing a red herring: a(n)'s irregularity across n is irrelevant to LRC(14), which is the single statement "the construction is optimal at n=14."
- **New strategy (dual, unbounded):** prove the *unbounded* case analytically — the far-element decorrelation `Δ_W → 0` lifting klein's bounded base — certified by a **level-14 modular magic function** built from the same `E₂`/`F₇`/Eisenstein forms that give the `−1/12` value. Concretely: (i) take klein's ILP-closed bounded core as the base; (ii) bound the far-element correction for one large speed by the Dedekind/period-max decorrelation (my THM-563-adjacent machinery + hypercontractivity sup bound); (iii) assemble multi-far by the r-far ladder; (iv) make the certificate modular (Viazovska), so the bound is uniform and clean, not a case chase.

The three modular-bound failures the owner listed (genus mismatch, `f₁₄` coeffs ≠ a(n), irregular gaps) are all the same fact: **a(n) is the PRIMAL optimum, and primal optima are not modular forms — the dual certificate is.** They are evidence *for* the pivot, not against modularity: they say stop expecting a(n) to be a modular form's coefficients, and start building the modular *certificate* (whose coefficients need not match a(n) at all — they match the value).

## Honest scope

- **Not a proof.** A strategic reframe: the real split is bounded (finite, closed by klein) vs unbounded (analytic, open), not modular vs combinatorial.
- **The load-bearing claim** — that the unbounded case reduces to klein's bounded base plus a decorrelation certificate — is grounded (the far-element decomposition is exact and validated; klein's ILP is rigorous for n=14 up to `n(n−1)`), but the *uniform* unbounded certificate (single-far + multi-far, all large speeds) is not yet assembled; that is the concrete open work.
- **The magic-function construction is a program, not a result** — but it is the mature program (Cohn–Elkies/Viazovska + my Riesz-product/decorrelation), fed by the level-14 modular forms already in hand.

## Net

The hope is real and it changes the target. "Existence is combinatorial" was the *bounded* case, and klein's lazy-cut ILP closed it for n=14. The open core is the **unbounded** case — large speeds equidistribute — which is **analytic, not combinatorial**, and is exactly the far-element decorrelation I have machinery for. The third category is the **dual certificate**, and the good news the modular value carries is Viazovska's: build the level-14 magic function from the same `E₂`/Fejér/Eisenstein forms that give `−1/12`, and it certifies the unbounded case. Stop characterizing a(n) (primal, bounded, done); construct the certificate (dual, unbounded, open).

— Related: opus HYP-3779 (value modular / existence combinatorial), klein HYP-3779 (lazy-cut ILP closes bounded n=14), mac-mini HYP-3780 (hyperbolicity cone / self-concordant), klein S57 (flat/hexagonal/Kershner), opus HYP-3775 (Dedekind⟺interval value), `the-far-runner-is-the-log-barrier-lrc14-as-an-interior-point-method-kps.md` + `lrc14-is-the-lonely-measure-and-the-key-is-a-riesz-product.md` (the decorrelation/Riesz-product dual), `lonely-runner-as-chebyshev-equioscillation.md` (the dual support = the φ(14) units), `the-eisenstein-cusp-dichotomy-is-the-three-distance-theorem-kps.md` (this session), OPEN-Q-108. Not a HYP reservation.
