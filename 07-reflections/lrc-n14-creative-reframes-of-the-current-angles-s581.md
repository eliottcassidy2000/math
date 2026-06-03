---
source: opus-2026-06-03-S581 (remote-control / "star")
status: BRAINSTORM + two rigorous reframes (verified on AP, V*). Creative reframes of the *current* n=14 angles of attack (C′ reduction, pinch/THM-396-397, paired-or-anchored split, 64 self-converse classes, resonance-energy θ-bound). Marks rigor vs conjecture per item.
tags: [LRC, n14, creative-reframes, Cprime, pinch, danger-budget, divisor-cover, dominance-dodge, bivariate-nullstellensatz, structure-vs-randomness, 6-cube, S550, S557, S558, S572, S574, S576, S580]
---

# LRC@14: creative reframes of the current angles of attack

**Prompt (user, remote-control):** spend this session trying many creative
reframes on the current angles of attack for n=14.

Grounding (the live frontier as of S580). The program has collapsed onto **one
crux in three dresses**:

- **Pinch dress (HYP-2095, S569):** every measure-zero 13-set has an *unblocked
  small-reduced-sum pair*.
- **C′ dress (THM-398, S572):** if `14 | v` for some `v∈S`, then `M(S) > 1/14`;
  the residual is "all-short" configs where one runner's thin danger-comb must
  cover the others' safe set.
- **Congruence dress (S574):** the large-owner endpoint congruence system is never
  simultaneously feasible.

Everything below is a *different lens on this crux* (or a way to enlarge the part
of the energy/structure axis that the existing tools already cover). Conventions:
`n=14`; observer + 13 speeds `S`; `‖x‖`=dist to ℤ; `M(S)=max_t min_{v∈S}‖vt‖`;
LRC(14) ⟺ `M(S)≥1/14`; "loose" = `M(S)>1/14`. `G(S)={t:‖vt‖≥1/14 ∀v∈S}` (the
level-1/14 lonely set). All numerics verified inline.

---

## Reframe 1 — The danger **budget** is exactly 1/7 (RIGOROUS lever)

**Observation (exact).** For *any* single speed `v`, the level-1/14 danger set
`D_v = {t : ‖vt‖<1/14}` is `v` arcs of width `2/(14v)`, so

> **`μ(D_v) = 2/14 = 1/7`, independent of `v`.**

One runner can poison at most **1/7 of the circle**, no matter how large. This is
the *currency* the whole C′ residual is denominated in. Immediate:

> **Lemma 1 (1/7-budget; rigorous).** If for *some* `v∈S` the 12-runner subsystem
> has `μ(G(S∖{v})) > 1/7`, then `S` is loose.
> *Proof.* `G(S)=G(S∖{v})∖D_v`; if `μ(G(S∖{v}))>μ(D_v)=1/7` it cannot be covered
> by `D_v`, so `G(S)≠∅`. ∎

**Contrapositive (necessary condition for a counterexample).**

> **A counterexample is "irredundant": for *every* runner `v`,
> `μ(G(S∖{v})) ≤ 1/7`.** Each of the 13 runners is essential — deleting any one
> still leaves the other twelve lonely on ≤1/7 of the circle.

**Why this is leverage (synthesis with S550).** The resonance-energy bound (S550)
applied to the **12-runner subsystem at threshold 1/14** gives
`μ(G(S∖{v})) ≥ (1−2/14)^{12} − E_{12} = (6/7)^{12} − E_{12}`. Now

> `(6/7)^{12} = 0.15730… > 1/7 = 0.14286…`, gap `≈ 0.01445`.

So **Lemma 1 fires (config is loose) as soon as some 12-subsystem has
`E_{12}(S∖{v}) < (6/7)^{12} − 1/7 ≈ 0.0145`.** Hence:

> **The entire C′ residual lives in the high-energy core of *every* 12-subsystem:**
> a counterexample needs `min_v E_{12}(S∖{v}) ≥ 0.0145`.

This is sharper than the plain S550 self-bound (which needed `E_{13} < (6/7)^{13}
≈ 0.135`): we spend the 13th runner's fixed 1/7 budget to lower the energy bar
we must clear from `0.135` down to `0.0145`, on a subsystem whose energy is also
smaller (one fewer runner ⇒ fewer resonances). Two opposing knobs, but the picture
is concrete and the threshold is a *number*.

**The geometric punchline (rigidity).** At a counterexample, `G(S∖{v})` has
measure ≤1/7 and is **exactly covered** by the periodic comb `D_v` (period
`1/(14w)` when `v=14w`). A periodic comb tiling a set forces that set to be
**invariant under `t↦t+1/(14w)`**. So the 12 other runners' *joint* lonely set
must be `(1/14w)`-periodic — a severe structural constraint (their danger pattern,
of mixed periods `1/u`, would have to organize into a single fine period). This is
the "all-short alignment" core (S572/S574) restated as an **exact-tiling
rigidity**, and it is plausibly *provably* impossible.

**Concrete probes.**
1. For the S574 all-short residual configs (the ~19–39 survivors at n=14),
   compute `μ(G(S∖{v}))` for each `v`; any `>1/7` is immediately discharged.
2. Test the periodicity rigidity directly: is `G(S∖{v})` ever `(1/14w)`-periodic
   for a primitive 12-set? If never (for `w≥1`), C′ is proved on the residual.

*Verified:* `μ(D_v)=1/7` exact; `(6/7)^12=0.15730`. (n=14.)

---

## Reframe 2 — Pinch-saturation **is** a divisor-cover (RIGOROUS restatement + a new bridge)

Combine THM-396 (a single universal blocker of a small pinch must be a *sum-multiple
shield* `D|c`) with THM-397 (a collective non-shield cover needs an *endpoint
anchor*, impossible when `D≤14`). The net, stated as a clean arithmetic fact:

> **Lemma 2 (divisor-cover; rigorous given THM-396/397).** If `S` contains a pair
> `a,b` with `a+b ≤ 14` and **no element of `S` divisible by `a+b`**, then the
> pinch time `t=m/(a+b)` clears all runners and `M(S) ≥ 1/14`.

Define `S` **pinch-saturated** iff for *every* pair `a,b∈S` with `a+b≤14`, some
`c∈S` has `(a+b) | c`. Then:

> **Equivalence.** A 13-set without the Lemma-2 witness ⟺ `S` is pinch-saturated.
> So: **LRC(14) ⟸ [ every pinch-saturated 13-set is loose ].**
> (The non-saturated case is *already proved* by Lemma 2.)

This converts the analytic crux into a **self-referential divisor-covering
condition on 13 integers** — finite, combinatorial, checkable.

*Verified:* AP's unblocked pairs are exactly `(1,13),(2,12),…,(6,8)` (sum 14, no
multiple of 14 present) — i.e. the pinch witness **is** the unblocked pair, matching
the S557 straddle table. V*={1..11,13,24} likewise (it lacks 12, so 5 pairs).

**The new bridge (pinch ⇒ dominance).** Pinch-saturation is *expensive*: it forces
`S` to contain a multiple of `D` for **every realizable small pair-sum** `D≤14`.

- To even *have* small pairs, `S` must contain small speeds; but then `D=14`
  forces a multiple of 14, `D=13` a multiple of 13, `D=12`→mult of 12, … . Stacking
  these forces large speeds.
- The two escape routes both *help us*:
  - **Vacuous saturation** (all pair-sums `>14`, i.e. speeds spread out / ≥8) ⇒ the
    config is "dominant/spread" ⇒ **dominance-dodge (THM-398 Lemma B)** and the
    1/7-budget (Reframe 1) make it loose.
  - **Genuine saturation** (small pairs present but all carry a divisor-multiple)
    ⇒ `S` contains a large multiple of a small `D` ⇒ a **dominant runner** ⇒
    again dominance-dodge.

> **Unification claim (conjecture, high-value):** *pinch-saturation ⇒ `S` contains
> a runner `v > 13·max(others∖v)` (a dominant runner) ⇒ loose by THM-398 Lemma B.*
> This would **merge the pinch line (HYP-2095) and the C′/measure line (THM-398)
> into one argument** — the two halves of the proof program become one.

**Concrete probe.** Enumerate pinch-saturated 13-sets in a bounded box (or prove
none are primitive). Specifically test the implication "pinch-saturated ⇒ has a
dominant runner": search for a *counterexample to the bridge* (a saturated set with
no dominant runner). If none exists up to height `B`, the bridge is strongly
supported and worth a direct proof via the divisor-stacking count.

---

## Reframe 3 — Break `14=2·7` with a **bivariate Nullstellensatz over 𝔽₂×𝔽₇**

S558 named the single most-leveraged target: the polynomial method of
Sungkawichai–Trakulthongchai proves the hardest tuple "proper" *for free* when
`k+1` is an **odd prime** (Combinatorial-Nullstellensatz needs a **field**), and it
**fails at `k+1=14`** precisely because `ℤ/14` is not a field — this *is* the
literature's wall at n=14.

**Reframe.** Don't work over the non-field `ℤ/14`. Use the CRT splitting into a
**product of two fields**:
> `ℤ/14 ≅ 𝔽₂ × 𝔽₇`,  `t ↔ (t mod 2, t mod 7)`.

Replace single-variable CN (Alon) over `𝔽_p` by a **two-variable Nullstellensatz
over `𝔽₂ × 𝔽₇`**: build indicator polynomials `P(x,y)` with `x∈𝔽₂`, `y∈𝔽₇`,
handle the 2-part and 7-part by *separate* (legitimate, field-based) CN
applications, and recombine the non-vanishing coefficient via CRT. The "tight
tuple proper" obligation that the prime case discharges analytically becomes a
**bilinear coefficient-nonvanishing** statement — one CN over `𝔽₂` (degree budget
1 per variable — almost trivial) ⊗ one CN over `𝔽₇`.

This is the concrete shape of "the **2q-analogue** of `(1,…,k)` is proper" that
S558 asked for, and it is **exactly our `2·7` wheelhouse** (the mod-2 fold ⊗ mod-7
singleton that Theorem B(2)/oracle-S552o showed obstruct on the *same* runner).
The same `𝔽₂×𝔽₇` factorization should also be the natural home for the
**certificate-compression `n`-clock ⊗ pair-clock** quotient (S580): the four clocks
are really `2`-adic ⊗ `7`-adic data.

**Status:** speculative but *the* right target. **Probe:** write the prime-`p`
tight-tuple CN polynomial from the Apr-2026 paper explicitly, then attempt the
`𝔽₂⊗𝔽₇` factorization on the smallest doubled-prime case `k+1=2·3=6` (where LRC
is *proved*, so the method can be *calibrated against a known answer*) before
pushing to `2·7=14`. Calibrating the new tool on `n=6` is the cheap, decisive test.

---

## Reframe 4 — **Margin-and-dodge** off the *proven* LRC(13) (concrete, uses S558 correction #2)

S558 flagged that we kept using LRC(7) when **LRC(13) is now proven**. Cash it in
quantitatively via the even-fold.

Let `s` be a LRC(13) witness on `fold(S)={u : 2u∈S}` (≤12 evens, so `M(fold)≥1/13`
by *proven* LRC(13)). Pull back to `t₀ = s/2`: every **even** runner satisfies
`‖2u·(s/2)‖ = ‖u s‖ ≥ 1/13 > 1/14` — with **slack `1/13−1/14 = 1/182`**. So at
`t₀` we *already win on all evens, with margin*.

> **Reframe of the residual.** The only obstruction is the **odd** runners at
> `t₀`. We may perturb `t∈(t₀−δ, t₀+δ)` and keep evens safe as long as
> `2u_max·δ ≤ 1/182`, i.e. `δ ≤ 1/(364·u_max)`. Inside this window we must **dodge
> the odd runners' danger arcs** (each odd `r` has arcs of radius `1/(14r)`,
> spacing `1/r`).

This is a *much smaller* problem: an LRC-style dodge among ≤13 odd runners, but
only on a single short window where the evens are pre-cleared. It quantifies the
"odd coupling residual" (S554) as a **covering count on one interval**: the window
has length `2δ`; odd runner `r` blocks ≤ `⌈2δ·r⌉` of its arcs there, total blocked
measure `≤ Σ_{odd r} 2/(14r)` per period — a Borel–Cantelli/`Σ1/r`-style budget. If
the odd danger budget on the window is `<2δ`, a dodge exists.

**Probe.** For each n=14 config, compute the even-margin window `δ` from the
LRC(13) witness and the odd danger measure inside it; check the budget inequality.
This turns the antipodal/odd-coupling residual into an explicit inequality
`Σ_{odd} (window-arc-count)·2/(14r) < 2δ`, parametrized by the LRC(13) witness.

---

## Reframe 5 — **Structure vs randomness** (the organizing dichotomy)

A counterexample is squeezed from both ends of the **resonance-energy axis**:

- **Low energy (random-like speeds).** Few short resonances ⇒ `E` small ⇒ the S550
  bound (sharpened by Reframe 1's 1/7 budget, threshold `≈0.0145` on a 12-subsystem)
  gives `μ(G)>0` ⇒ **loose**.
- **High energy (structured / AP-like speeds).** Many additive relations ⇒ speeds
  are Freiman-structured (close to an AP) ⇒ **even-fold + proven LRC(13)** (Reframe 4)
  and the **divisibility sieve** (THM-369) bite, because structure = detectability.

> **The counterexample must be neither random nor structured — an empty middle.**
> This is the Tao-style structure/randomness transference applied to LRC: the two
> *already-working* tools (measure bound; fold+sieve) cover the two ends, and the
> proof is finished by making their coverage **overlap on the energy axis**.

This is the right *narrative skeleton* for the paper: stop treating measure-bound,
fold, and sieve as separate cases and instead state one **energy-threshold
partition** `[0, e_lo] ∪ [e_hi, ∞)` with `e_lo ≥ e_hi`. **Probe:** quantify `e_hi`
= the largest energy still reachable by a fold/sieve-reducible config, and `e_lo` =
Reframe-1 threshold; the open work is exactly the band `(e_hi, e_lo)`, and the goal
is to prove it empty (i.e. `e_lo ≥ e_hi`).

---

## Reframe 6 — The 64 self-converse classes as a **6-cube; prove per-edge monotonicity**

The worry-set collapses (HYP-2094) to `2^6 = 64` self-converse round classes =
flip-sets `F ⊆ {6 antipodal pairs mod 27}`, i.e. the **vertices of the 6-cube
`{0,1}^6`**, with the **AP at the origin** and (S570) every single flip strictly
increasing `M` above the floor `1/14`.

> **Reframe.** Define `Φ:{0,1}^6→ℝ`, `Φ(F)=M(S_F)`. LRC-on-the-64 ⟺ `Φ≥1/14`
> everywhere. We know `Φ(0)=1/14` exactly and `Φ(e_i)>1/14` for each unit flip.
> **Target: `Φ` is "discretely monotone from the origin"** — every coordinate flip
> *up* never decreases `M`. Then `Φ(F) ≥ Φ(0)=1/14` for all `F`, done.

This replaces a continuum statement over 64 classes by a **single per-edge sign
lemma** ("flipping antipodal pair `i` from a transversal does not lower `M`"),
plausibly supplied by the **two-gap monotonicity** dynamics (S519/THM-387). The
discrete-derivative structure on `{0,1}^6` is exactly where Fourier-on-the-cube /
monotone-Boolean machinery applies.

**Caveat (the honest gap, per S570):** the 64 cube covers the *transversal*
realizations (bounded speeds); the **structural lift** to general speed sets
realizing one of the 64 classes (e.g. V*-type non-transversal cousins) is the real
open step. So: prove the cube monotonicity *and* the lift `general realization →
transversal representative without lowering M` (this is the S576 unit-spine
exchange lemma in disguise). **Probe:** compute `Φ` on all 64 vertices and verify
monotonicity edge-by-edge from the origin; identify any non-monotone edge (would be
the precise obstruction to attack).

---

## Reframe 7 — The certificate as an **LP-dual measure on pinch times** (revive BHK averaging)

`M(S)≥1/14` is feasibility of a covering LP; by the Pinch Lemma (S557) its optimum
sits at a 2-runner-binding vertex, so the *only* times that matter are the pinch
times `{m/(a+b)}`. LP duality ⇒ a certificate is a **probability measure `ν` on
those finitely many pinch times** with `∫‖v_i t‖ dν ≥ 1/14` (suitably) for all `i`.

This is precisely the **Bohman–Holzman–Kleitman averaging method** that *proved
n=6* — and S558 flagged it as **underused for `n≥8`**. The reframe: build the n=14
certificate as an *explicit averaging measure* supported on `{m/D : D=a+b≤14}`,
rather than as a single witness time. The pinch structure tells us the support a
priori; we only need the *weights*.

> **Target:** an explicit `ν` (e.g. uniform on the `≤6` floor-pair pinch times
> `{(2k+1)/14}` plus corrections) that certifies every 13-set, with the weights
> chosen by the THM-396/397 shield/anchor data. The "certificate sheaf / four-clock
> quotient" (S579/S580) is then literally the **support of the dual measure**.

**Probe.** Solve the finite LP (max-min over pinch times) for a sample of hard
configs; read off the dual optimal `ν`; check whether a *single config-independent*
`ν` (or a small finite menu indexed by the 64 classes) certifies all — that would
be a BHK-style uniform averaging proof.

---

## Synthesis — which reframe attacks which open target

| open target (current frontier) | reframe(s) that bite |
|---|---|
| all-short / large-owner congruence residual (S574) | **R1** (1/7 budget ⇒ residual is `E_{12}≥0.0145`; exact-tiling rigidity) |
| unblocked-small-pair lemma / paired-or-anchored (HYP-2095) | **R2** (divisor-cover restatement) |
| merge pinch-line and C′/measure-line | **R2** bridge (saturation ⇒ dominant runner) |
| the `2·7` polynomial-method wall (S558 #1) | **R3** (bivariate CN over 𝔽₂×𝔽₇) |
| odd-coupling residual after even-fold (S554) | **R4** (margin-and-dodge off proven LRC(13)) |
| organize measure-bound + fold + sieve into one proof | **R5** (energy-axis partition `e_lo ≥ e_hi`) |
| prove the 64 self-converse classes lonely (S570/S576) | **R6** (6-cube monotonicity + spine-exchange lift) |
| certificate compression / four-clock completeness (S580) | **R7** (LP-dual measure = certificate support; BHK) |

**Two are rigorous now** (R1's 1/7-budget lemma; R2's divisor-cover equivalence,
given THM-396/397) and immediately shrink the residual / restate the crux. **R3
and R5** are the highest-ceiling: R3 attacks the literature's exact wall, R5 is the
paper's organizing principle. **R4 and R6** are the most concrete next computations.

**Suggested new hypotheses:**
- **HYP-2106** (R1): a counterexample is irredundant — `μ(G(S∖{v}))≤1/7` for all
  `v`; equivalently `min_v E_{12}(S∖{v}) ≥ (6/7)^{12}−1/7`.
- **HYP-2107** (R2 bridge): pinch-saturation ⇒ `S` has a dominant runner
  `v>13·max(others)` ⇒ loose. (Merges HYP-2095 with THM-398 Lemma B.)
- **HYP-2108** (R6): `Φ:{0,1}^6→ℝ`, `Φ(F)=M(S_F)`, is monotone from the AP origin.

**Artifacts:** verification inline (exact: `μ(D_v)=1/7`; `M(AP)=M(V*)=1/14` at
units/14; AP unblocked pairs `=(1,13),…,(6,8)`). Builds on S550 (energy), S557
(pinch), S558 (methodology map; LRC(13) proven), S569 (paired/anchored), S572
(C′), S574 (small-owner), S576 (unit spine), S580 (certificate). New: HYP-2106,
HYP-2107, HYP-2108.
