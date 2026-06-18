---
id: THM-525
title: "The easy-dominates-hard covering reduction of LRC(14): hard covering 13-sets are an easy LRC(12)-core C plus a runner parked in section 0; centering the parked runner ('perfect middle of section 0') is a constructive certificate that localizes the covering case onto OPEN-Q-108, and the loneliness-derived tournament maps factor through LOCAL (round) tournaments"
status: MIXED — see per-part marks. PROVED (elementary): STEP-1 q-witness universal bound, STEP-2 measure identity, survivor-sufficiency lemma, meas(T_w)=6/7, the per-arc comb inequality, the diff-wind tie-free⇒local lemma at n=5. VERIFIED (exact compute): M({1..11,13}∪{84m})=7m/(84m+5) and edge-binding (m=1..8); meas(G_C) extremal 7/858; ~105k covering 13-sets zero counterexamples; the two forbidden tournament classes at n=5. CONJECTURE / OPEN: the reduction reaches nothing weaker than OPEN-Q-108 (GAP A) + the transversality estimate (GAP G2). **NO claim of M≥1/14 is proved; LRC(14) is NOT proved.**
source: kind-pasteur-2026-06-17-S2 (user reframe: hardest configs park a runner in section 0 forever; center it in the perfect middle; easy & hard come hand in hand; use the easy cases' structure to kill the hard cases — plus the creative tournament-map hunt). Built by a workflow (3 structure investigators + 3 perspective-diverse reduction skeptics + 8-theme tournament fan-out with adversarial refutation + synthesis/completeness critic), independently re-verified.
depends_on:
  - THM-524   # binding-pair reduction + regions/sections reframe (covering = off-grid, parked runner)
  - THM-523   # covering needs a multiple of every q in 2..14 (the hard core IS the covering family)
  - THM-522   # measure-side scale-invariance + quantization
related:
  - OPEN-Q-108  # the uniform fattening lemma = THE crux; this localizes the covering case onto it
  - HYP-2573    # centering certificate (REFUTED as optimizer behavior; CONJECTURE as a device)
  - HYP-2574    # easy↔hard correspondence (VERIFIED, lcm-parking subfamily)
  - HYP-2575    # section-shuffling census (VERIFIED; sharing ANTI-correlates with hardness)
  - HYP-2576    # tournament maps factor through local/rotational structure (VERIFIED n=5)
external: proven LRC(12)/LRC(13) (Sungkawichai–Trakulthongchai 2026) — the EASY cores; Goddyn–Wong tight family; Bang-Jensen local (round) tournaments.
---

# THM-525 — Easy dominates hard: the centered-parked-runner covering reduction

**One-line.** A covering 13-set (the only place a LRC(14) counterexample can hide, THM-523) is
`S = C ∪ {w}` with `w ≡ 0 (mod 14)` parked in section 0 and `C` a 12-speed **easy core** that
proven LRC(12) already makes lonely at gap `≥1/13`. The covering case of LRC(14) reduces — by a
clean, non-circular implication — to a single uniform-measure estimate (OPEN-Q-108). The "perfect
middle of section 0" turns out to be a **constructive lower-bound device**, not what the optimizer
does. Separately, every non-trivial loneliness-derived tournament map factors through **local
(round) tournaments**, so the non-local iso-classes are forbidden — by circle geometry, not by τ.

`M(S) = max_{τ∈[0,1)} min_{v∈S} ‖vτ‖`; LRC(14) ⟺ `M(S) ≥ 1/14` for all primitive `|S|=13`.

---

## A. The reduction (the Q ⟹ P localization)

Let `P` = LRC(14); `Q` = OPEN-Q-108 = "∃ c>0 with `meas(G_C) ≥ c` for every 12-subset `C`", where
`G_C = {τ : ‖vτ‖ ≥ 1/14 ∀v∈C}` is the gap-`1/14` lonely set of `C`.

- **STEP 1 — non-covering case CLOSED (PROVED, THM-523).** If `S` omits a multiple of some
  `q ∈ {2,…,14}` then `M(S) ≥ 1/q ≥ 1/14` (witness `τ` near `a/q`). Re-verified: 0 failures over
  3647 non-covering primitive 13-sets. So WLOG `S` is **covering** and contains `w ≡ 0 (mod 14)`.
- **Bookkeeping correction (VERIFIED, a real trap).** A *13-speed* covering set must **drop ≥1
  small speed** to make room for `w`: `{1..13} ∪ {14m}` has **14** speeds (it is the `N=15`
  problem), not ours. The peel rule is well-defined: peel the smallest multiple of 14; `C` is a
  12-core with `meas(G_C) ≥ meas(G_S)`.
- **STEP 2 — measure identity (PROVED).** `L(S) = meas(G_C) − meas(G_C ∩ D_w)`, where
  `D_w = {τ : ‖wτ‖ < 1/14}` is `w`'s danger comb and `L(S) = meas(G_S)`; recall `L(S)>0 ⟺ M(S)>1/14`.
- **STEP 3 — decoupling floor (PROVED, but WEAKER than it looks).**
  `L(S) ≥ (6/7)·meas(G_C) − r(C)/(7w)`, `r(C)=#mouths of G_C`; positive once `w > w₀ = r/(6·meas(G_C))`.
  **Caveat (the key sharpening):** `w₀` grows ≈ **linearly in the largest core speed** because `r`
  does (`w₀ ≈ 54.8, 139, 1014, 9762` as `v_big = 13, 500, 5000, 50000`). So the *proved* certificate
  is **vacuous for large-speed cores even with a single parked `w`**; below `w₀` the case rests on
  direct exact computation, not the bound. (E.g. the drop-6 core: floor `>0` at `w=84` (`4/21021`)
  but `<0` at `w=14` (`−237/7007`).)
- **STEP 4 — the residual (OPEN = GAP A = OPEN-Q-108).** The `k ≥ 3` arithmetically-coordinated
  growing-speed regime is uncontrolled: peeling speeds one at a time needs *every* intermediate
  core to keep `meas ≥ c` uniformly — that **is `Q` unchanged.**

> **The reduction is a genuine, non-circular implication `Q ⟹ P`** (verified by an adversarial
> "circularity" skeptic): STEP 1 is unconditional, STEP 2 is an identity, STEP 3 uses the
> *separately-proven* LRC(12) — none assumes `P` or `Q`. **But it reaches nothing strictly weaker
> than `Q`.** "One estimate away" undersells it: that one estimate is conjecture-equivalent to the
> whole open problem (finiteness of the n=13 primitive tight locus).

**Evidence (VERIFIED, not proof):** ~33k + ~64k + 8.6k ≈ **105,000 covering primitive 13-sets,
zero with `L=0`, zero tight, zero counterexamples.** Closest coordinated approach `M = 5/53 ≈ 0.0943`
(far above `1/14`). New datum: driving two coordinated parked speeds to ∞ on a fixed core does **not**
drive `L→0` — `L` **plateaus ≈ 0.0238**, coordinated 12-cores keep `meas(G_C) ≥ ~0.095 ≫ 7/858`; the
feared "`meas(G_C)→0` as `v_max→∞`" failure mode did not materialize in any probe.

---

## B. Centering the parked runner — the "perfect middle of section 0"

The user's directive, made exact. Two distinct `τ`'s must never be conflated.

- **At the TRUE optimum, the parked runner is at the EDGE, not the middle (VERIFIED, m=1..8).**
  For `S = {1..11,13} ∪ {84m}` (84 = lcm(12,14)), `M(S) = 7m/(84m+5)` exactly (min `7/89` at m=1,
  `→1/12`). At `τ*` the parked `w=84m` is **binding** (`‖wτ*‖ = M`) with `frac(wτ*) ≈ 0.92`
  (m=1: `82/89`) — it **hugs the boundary** of its section-0 safe band, **co-binding with core
  runner 5** (`84m+5` is the optimum denominator, a THM-524 binding pair `{5, w}`). So "the
  optimizer keeps `w` in the perfect middle" is **REFUTED.**
- **The constructive centering DEVICE works (the correct reading of "perfect middle").** Choose
  `τ ∈ G_C` maximizing `‖wτ‖`: this pushes `w` toward `1/2` (its perfect middle — reaches `1/2`
  exactly for several `w`, slack up to `3/7`) while a **core runner** becomes the `1/14`-binding
  floor. The full gap is then exactly `1/14` — a valid **certificate** in every tested case, never
  below `1/14`. It is core-dependent: for the AP-drop-6 core the perfect-middle center is literally
  the dihedral-clock event `τ = 9/56 = (14·2−1)/(14·12)` (where `‖84τ‖=1/2`, core runner 12 binds);
  for the hard core `{1..11,13}` no band-center `(2j+1)/168` lands in `G_C`.
- **PROVED atoms.** *Survivor sufficiency:* `G_C ∩ T_w(1/14) ≠ ∅ ⟹ M(S) ≥ 1/14` (trivial from
  defs), where `T_w(h)={τ:‖wτ‖≥h}`. *Comb measure:* `meas(D_w)=1/7`, `meas(T_w)=6/7` for every `w`.
  *Per-arc inequality:* `cover_w(A) ≤ |A|/7 + 1/(7w)` (exhaustive-exact, even conservative).
- **The two named gaps.** **GAP A** = `Q` (above). **GAP G2** (the transversality/fattening
  estimate): even granting `Q`, prove `w`'s danger comb (measure `1/7`, concentrated near grid
  points `a/w`) **cannot contain all of `G_C`** — i.e. `G_C` avoids `O(1/w)`-neighborhoods of
  `{a/w}`. The intersection was nonempty in **every** computed case; **no general argument exists.**
- **REFUTED route (do not retry — HYP-2573 mechanism / GAP B):** the naive "LRC(12) + Lipschitz on
  one safe arc" lever. The LRC(12) safe arc has half-width `(M(C)−1/14)/v_binder`, which **shrinks
  like `1/v_max`** when the binders are large speeds (drop-6: `5/3864 ≈ 0.00129 ≪ 1/91`). One must
  bound **total measure**, not one arc.

---

## C. Section-shuffling / SDR census (VERIFIED) — and the easy↔hard correspondence

- **Section-sharing is ANTI-correlated with hardness (CORRECTS the framing).** The perfect SDR /
  AP `{1..13}` is the **tight FLOOR** (`M=1/14`), *not* the easy end; sharing a nonzero section
  *raises* `M` (1 collision→`1/13`, 2→`1/12`, …; non-strict, k=6 dips to `2/21`). The single sharp
  hardness coordinate is **`z` = #runners in section 0 = #multiples of 14**: `z=0 ⟹ M ≥ 1/13`
  (LRC(13)); a counterexample needs `z ≥ 1`, killing grid-loneliness and forcing the optimum
  off-grid. Census: 373 section-sharing patterns; 101 grid-loneable; **exactly 1 is the SDR.**
- **Easy↔hard correspondence (PROVED, lcm-parking subfamily; HYP-2574).**
  `{C ∪ {84m}}` = (easy 12-SDR core with exactly **one** nonzero hole) + (a 13th runner forced into
  section 0 by `84m ≡ 0`). **`(ℤ/14)*`-compatible:** `0` is a fixed point of unit multiplication,
  so the pattern is a `(ℤ/14)*`-invariant (verified 20000 multisets, 0 failures); covering→covering,
  core→core. Hole positions fall in 3 orbits: `{1,3,5,9,11,13}, {2,4,6,8,10,12}, {7}` (antipode 7
  unit-fixed). **No Hall theorem** at any modulus (SDR-ability is the trivial "distinct & nonzero",
  since mult-by-`a` is a permutation); the genuine obstruction is **band-avoidance (covering)** —
  a covering-design problem, not a matching problem.
- **Scope (open):** proved only for the lcm-parking subfamily; a general covering set's section-0
  runner need not be `84m`. This lens is on-grid and **blind to the off-grid optima** that decide
  LRC(14).

---

## D. The tournament ledger (the creative-map hunt) — honest verdict

**The prize as literally posed was NOT won: no map yields a *loneliness*-forbidden iso class.**
Across all non-trivial maps (switch-parity, danger-interval, residue-orbit, pinning,
character-spectral, first-return, pairwise-gap-load), the adversarial control returned
**forbidden-by-loneliness = 0**: the iso-class set realized at the lonely-optimal `τ` **equals** the
set realized at arbitrary `τ`. The choice of lonely time does not constrain the achievable
tournament beyond the map's own algebra. (This generalizes the dead overtaking/snapshot map.)

**But a genuine structural result survived adversarial refutation (VERIFIED n=5; HYP-2576):**

> **Loneliness-derived pairwise maps factor through LOCAL (round) tournaments.** The
> **difference-winding** map `i→j ⟺ frac((v_i−v_j)τ) ∈ (0,1/2)` *is* the circular/phase tournament
> on the points `a_v = frac(vτ)` on `ℝ/ℤ`; **tie-free ⇒ local (round) tournament** (0/45108
> non-local). Hence the **unique maximal-H non-local n=5 class** — score `(1,2,2,2,3)`, `c3=4`,
> `H=15` — is **unreachable** at the lonely optimum *and* at free phases. The signed
> danger-arrival map (M3) and the section-rotational map (M4) forbid the **same** class (and M3 the
> full set of 8 non-rotational classes), making it the canonical "structurally-forbidden"
> tournament. Independently reproduced (`lrc14_diffwind_local_tournament_kps-S2.py`): diff-wind at
> the lonely optimum realizes 10/12, tie-free random phases realize only the 7 local classes,
> target unreachable in both. **Cause = circle geometry (the LRC lives on `ℝ/ℤ`), NOT loneliness.**

**Dead maps (do not retry):** any single-scalar key → total order → transitive (`H≡1`):
forward-clearance, binding-slope lex, approx-quality, Stern-Brocot ancestry, three-distance gap,
first-hit-of-section-0, separation-parity (degenerate constant map). **Realize-everything controls:**
meeting-parity×speed, CF-parity, Ostrowski-depth (full A000568 at n≤5).

---

## E. Status & next step

- A (reduction): `Q⟹P` **PROVED non-circular**; STEP 1 **PROVED** (non-covering closed); STEP 3
  **PROVED but vacuous for large-speed cores**; STEP 4 = **OPEN (=Q)**. 105k-set scan = **VERIFIED
  evidence**. **No M≥1/14 proved.**
- B (centering): survivor-sufficiency + `meas(T_w)=6/7` **PROVED**; `7m/(84m+5)` + edge-binding
  **VERIFIED**; the universal centering certificate **CONJECTURE** (rests on `Q`+`G2`).
- C (sections): census **VERIFIED**; easy↔hard correspondence + `(ℤ/14)*`-compatibility **PROVED**
  (lcm-parking subfamily).
- D (tournaments): local-tournament factoring + forbidden class **VERIFIED n=5 / PROVED mechanism**;
  **CONJECTURE** universal-n; loneliness adds **zero** iso-class constraint (clean negative).
- **NEXT (the one that would close it):** a **bounded-speed reduction** (Tao Thm 1.3 / MSS-style):
  prove LRC(14) on covering sets need only be checked for `v_max ≤ V₀`; then the ~105k scan becomes
  a finite *certificate* up to `V₀`. Failing that, attack **GAP G2** directly: `G_C` (a union of
  arcs with binding-crossing endpoints `k/(v_a±v_b)`, all `v ≤` bound) cannot concentrate in
  `O(1/w)`-neighborhoods of `{a/w}`. Either one — not both — moves OPEN-Q-108 from "localized +
  strongly evidenced" to "closed."

Scripts: `04-computation/lrc14_{parked_centering,easy_dominates_hard,section_shuffle_census,
verify_reduction_*,diffwind_local_tournament,tourmap_*,refute_*}_kps-S2*.py`; outputs in
`05-knowledge/results/` (same stems `.out`).
