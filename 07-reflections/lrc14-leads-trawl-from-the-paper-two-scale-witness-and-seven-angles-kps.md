# LRC(14) leads trawl: the two-scale witness and seven angles from the paper

*kind-pasteur-2026-06-27-S31ag. A real-time trawl of leads thrown up by reading arXiv:2604.23906's
proof techniques against the project's state. Most are flagged SPECULATIVE; a few are SOLID. Shared for
the team to pick up — the goal is breadth, not closure.*

## The organizing structure: the two-scale witness `t = s/14 + r/p`
The paper's witness (Prop 4.4) is **two-scale**: a coarse part `s/(k+1) = s/14` plus a fine part `r/p`,
with `{t u_i} = {s u_i'/14 + r·i/p}`. The coarse scale is the apex denominator 14; the fine scale is a
large prime `p`. The gap lemma (Lemma 4.3): consecutive discontinuities of the coarse discretizer are
`≥ 1/(k(k+1)) = 1/182` apart, so for `p > k(k+1) = 182` the fine grid refines the coarse one. **This is
the paper's version of the project's three-gap / Node-1 lemma (THM-565), with the explicit constant
`1/182`.**

> **SOLID LEAD 1 — why covering blocks the coarse scale (and is therefore the only hard case).**
> A coarse witness `m/14` is `1/14`-safe iff `14 ∤ s_i m` for all `i`. If `m` is a unit mod 14, this is
> `14 ∤ s_i`, i.e. **`m/14` works ⟺ no `s_i ≡ 0 mod 14` ⟺ non-covering.** For covering sets the
> multiple of 14 lands on `0` at *every* `m/14`, so **no coarse witness exists** and the witness must use
> the fine scale `r/p` — which is exactly the measure-floor / equidistribution nut. This is the precise
> reason "covering = hard": it is the vanishing of the paper's good-residue set `G` mod 14. *(Confirms
> the S59 redirect from the paper's own machinery; sharpens "covering is hard" to "the coarse `G` is
> empty.")*

## Seven angles

> **SOLID LEAD 2 — the reflection-Perron block IS the project's complement (T→T^op / SC) symmetry.**
> HYP-3089's reflection symmetry of the pairwise matrix `M` (`a↔6−a`) comes from the phase involution
> `x→−x`, under which `frac(e_i x)→1−frac(e_i x)` and sector `s→6−s`. That involution is the
> **tournament complement `T→T^op`** (kps-S35, lrc-witness-is-complement-symmetric). So mac-mini's
> reflection-Perron half-block is the project's **complement-even (self-complementary) subspace**, and
> the dominant Perron mode lives in the SC part. *Actionable for CRUX 1:* the merged-metagraph / SC-spine
> tools (which already handle the complement-even part of every project invariant) should apply to the
> `S2` Perron bound. The "fixed sectors" `{3,6}` (sector 3 = antipode, sector 6 = the one mapping to the
> anchor) are the SC-fixed coordinates.

> **SOLID LEAD 4 — the lonely measure is a signed relation-lattice sum, zero exactly on the tight locus.**
> `meas(L) = Σ_{n∈Λ(E)} ∏_i ĉ(n_i)`, `Λ(E)={n: Σ n_i e_i = 0}` the relation lattice, `ĉ` the Fourier
> transform of the `6/7`-safe arc — and **`ĉ(n)=0` whenever `7|n`** (Fejér zero, apex-7). The `n=0` term
> is `(6/7)^{13}`; the signed remainder cancels it to **exactly 0 on the tight AP** (which has open
> lonely measure 0). So the floor `meas(L) ≥ c` is "bounded away from the finite tight locus `{AP,GW}`"
> — and the only surviving relations have **all coordinates coprime to 7** (the support condition). The
> tightness is the exact zero-crossing of a signed cyclotomic sum.

> **SPECULATIVE LEAD 5 — Mertens/ζ(2) is the average good-residue density.** The project's provable floor
> `3/π² = 1/(2ζ(2))` (HYP-2856) should be the **`p`-averaged** density of the paper's good-residue set:
> `avg_p |G_p|/p`. The paper computes `I(k,p,1)` per prime; averaging the good density over primes is a
> Mertens product `∏_p(1−1/p·…)` giving a `ζ(2)^{-1}` constant. *Test:* compute `|G_p|/p` for the
> covering-binding family over many primes `p` and check the average against `3/π²`. If they match, the
> project's analytic floor is literally the paper's good-density on average — a bridge between the two
> programs.

> **SPECULATIVE LEAD 6 — Fano / Hamming[7,4] covering radius.** The 7 sectors = the Fano plane points
> (mac-mini's two-Fano fact for T7). The danger combs cover `1/7` each; "all 7 sectors hit" (`S7`) is a
> **covering condition on `F_2^3` via the Hamming[7,4] code** (whose 7 coordinates are the Fano lines).
> The cap `cap_k` may be a covering-radius / MacWilliams quantity for Hamming[7,4]. *Test:* compare
> `cap_8,cap_9,cap_10` to Hamming[7,4] covering-design numbers.

> **SPECULATIVE LEAD 7 — Ostrowski / continued-fraction address of the hard configs.** The witness
> denominators form a "good ladder" `{12,14,41,53,55,65,67,79,…}` (HYP-2866); these are best-approximation
> denominators. The three-gap structure (THM-565) is the Stern–Brocot/CF skeleton. *Lead:* index the hard
> covering configs by the **continued-fraction expansion of the binding ratio** (the slope where the
> deepest sink sits, denom 41 = the recurring Farey neighbor of 1/14). The hard cases should be exactly
> those whose CF has a specific bounded-partial-quotient pattern.

> **CONCRETE LEAD 8 — TESTED: the sieve and the hard core sit at OPPOSITE ends of the `|H_7|` axis.**
> THM-573/574 fail at `|H_7| = 6` only if the `13−6 = 7` speeds coprime to 7 kill **7 distinct** lifts
> (a perfect/derangement system); any overlap leaves a survivor. **Test result
> (`lrc14_residual_H7eq6_structure_kps.py`): over 3000 primitive covering 13-sets with `|H_7| = 6`, ALL
> have `M > 1/14` (`minM = 0.12 ≫ 1/14`), ZERO tight.** So `|H_7| = 6` is *benign* — overlap essentially
> always leaves a survivor. The genuinely hard/tight configs (dilated AP, GW) sit at **`|H_7| = 1`** (the
> AP `{1..13}` and GW each have exactly one multiple of 7, namely 7; `2·AP` also has `|H_7| = 1`). So
> **the level-7 sieve closes the apex-MAJORITY end (`|H_7| ≥ 7`), which is EASY; the binding core is at
> the SMALL-`|H_7|` end (`1,2,3`), the opposite end, squarely unreached.** This is honest and sharpening:
> THM-573/574 are clean but tangential to the difficulty — the measure floor on **small-`|H_7|` covering
> configs** (dilated AP/GW + perturbations) is the true residual, and there the `c=7` lift cannot fire.
> (Confirms HYP-3084's "hard tight cases have `|H|` small.")

> **CONCRETE LEAD 9 — TESTED: second moment (pairwise) is NOT enough for the cover bound; difficulty is
> order-3/4 in the LHS extremality.** Since THM-576 makes the *cap* (RHS) pairwise-exact, I tested whether
> the *cover bound* `meas(S7)=P(N=0) ≤ cap_k` closes by **Paley–Zygmund** `P(N=0) ≤ 1−S1²/E[N²]` (pairwise
> only). **Result (`lrc14_paley_zygmund_cover_bound_kps.py`): it fails at the binding config** — for
> consec_8 it gives `0.504` vs true `P(N=0)=0.327` and `cap_8=0.381`; only `2/302` configs (incl. consec)
> need more than pairwise, but they are exactly the binding ones. So the asymmetry is sharp: **the cap is
> order-2, but the cover-bound LHS extremality needs order-3 (S3, closing k=9,10 via mac-mini's
> `L_yK9=18−13S1+8S2−3S3`) and order-4 (S4, only k=8 via `L_yK8=…−9S3+6S4`).** The whole LRC(14) binding
> difficulty is "consec maximizes a degree-3 (k=9,10) / degree-4 (k=8) moment functional of `N∈{0..6}`" —
> the reflection-Perron crux. Consec is the meas(S7)-maximizer at every binding row (verified), with
> margin to cap (0.054 at k=8). This pins the order and the extremal config exactly.

## Net
Lead 1 (covering ⟺ empty coarse `G`) and Lead 2 (reflection-Perron = complement/SC) are solid and
immediately useful. Lead 9 pins the cover-bound difficulty to order-3/4 at consec (PZ/pairwise fails). Leads 5,6 are bridges worth a quick numeric test. Lead 8 sharpens the residual and
is tested next. None closes LRC(14); together they re-tool the endgame around the paper's two-scale
witness and the project's complement symmetry.

→ HYP-3087/3088/3089, THM-573/574/565, HYP-2856 (3/π²), HYP-2724 (relation lattice), HYP-2866 (good
ladder), mac-mini-S60 (reflection-Perron), arXiv:2604.23906 (Prop 4.4, Lemma 4.3), [[lrc14-thread]].
