---
source: klein-2026-07-23-S405
status: SYNTHESIS / IDEA-GENERATION per owner directive ("mine past work, let connections pop, find useful
  insights the snippet's STRUCTURE gives our open problems"). Five ranked structural transfers from the
  eq(27) fragment into the repo's live LRC(14) machinery. Each is a hypothesis + concrete next action, NOT
  canon. The fragment's exact form is settled (klein-S404); this is about what its SHAPE lends us.
tags: [lrc14, eq27, structural-transfer, second-moment, log-energy, rapidity, adelic, certified-log, ideas]
---

# What the eq(27) snippet's structure lends the repo — five transfers

**klein-2026-07-23-S405.** The fragment is a **certified log-energy floor**: a second-moment-weighted
linear form in two rapidities, `(2457/6592)·2artanh(t_B) − 2artanh(t_A) > 1/25`, proven float-free by a
truncate-plus-geometric-tail artanh sandwich. Stripping the (bespoke) numbers, its reusable STRUCTURE is
four primitives: **(P1)** certified two-sided rational bounds on a transcendental; **(P2)** a linear form in
logs/rapidities bounded *below by a clean rational floor* (not near-zero); **(P3)** a **second-moment weight**
`2457 = 3·Σ_{k=1}^{13}k²`; **(P4)** an **odd-harmonic** kernel `2artanh(t)=Σ_{k odd}t^k/k`. Where each plugs in:

## The arithmetic bridge (why this is plausibly *our* arithmetic, not a coincidence)
`2457 = 27·91 = 3·819`, where **`91 = C(14,2) = Σ_{1}^{13}k`** is the repo's LRC(14) pair-core constant
(`91^6` height bound THM-2052; `1/91` in THM-2053/THM-965), **`819 = 9·91 = Σ_{1}^{13}k²`** the extremal
second moment, and **`729 = 27² = 8·91+1`** a tracked repo anchor (LTI-385). The pair-covariance floor
`δ_(a,b) ≥ −6/637` has `637 = 91·7`; so the snippet's weight numerator (`91·27`) and the repo's covariance
denominator (`91·7`) are BOTH `91`-multiples. If the fragment is LRC(14), its weight is *the config's second
moment* and it lives in the exact `91`-arithmetic of our certificates.

## T1 (most actionable) — the certified-artanh engine as a Lean-ready primitive
opus-S2 already built `log_lower/log_upper(P,Q)` (rational two-sided `log(P/Q)` bounds). It is the "missing
continuous layer" beside `reciprocal_block_bounds` in `SupportHarmonicFigurate.lean`. **Transfer:** any repo
certificate that reduces to `log(rational) ⋛ rational` becomes kernel-pure — THM-2000 mass orderings (opus
already did `M(6,2)>M(4,3)` float-free), and any log/rapidity-valued threshold. **Action:** land the Lean
lemma `log_ratio_{lower,upper}` (P2-style), the first float-free transcendental-comparison brick.

## T2 (deepest, highest payoff) — a log-energy relaxation *complementary* to the polynomial moment-LP
Our density/covering floors run on POLYNOMIAL moment LPs: provable `D3=0.3088`, tight `d→∞` `=νConsec=0.4425`
(THM-661; opus-S185). opus-S185 proved the *naive* 2nd moment (Paley–Zygmund on empty-arc counts) is dominated
by that LP because `ν=μ=P(W>0)`. **But the snippet's second moment is used differently** — as a multiplicative
weight on a **log-energy / logit** form (P3×P4), a functional OUTSIDE the polynomial-moment hierarchy. Entropy
appears in our LRC files only as a heuristic *correlate* (LTI-143 phase entropy; `corr(p0,residue_entropy)=+0.54`),
never as a rigorous log-energy *certificate*. **Hypothesis:** a log-energy lower bound may reach past `0.3088`
at LOWER degree/cost than the polynomial `B_d` (entropy vs L² relaxations differ; log can be tighter), possibly
sidestepping opus-S185's "coupled-region difficulty." **Action:** on the tight `consecutive` cluster, compute
the Σv²-weighted log-energy analogue of `B_d` and compare to `D3=0.3088` at matched degree. Honest caveat:
`ν=μ` means it cannot beat the tight `0.4425`; the bet is *efficiency* (fewer terms to a usable floor), not
a new ceiling. This is the concrete experiment worth running next.

## T3 (revive dormant thread) — the archimedean rapidity floor for the adelic proof-shape
THM-252 puts flip-chain rapidities in the log-prime lattice with the adelic product formula `∏_v|x|_v=1`; the
`rapidity-supersingular-adelic` shape (oracle-S19) needs a **certified archimedean lower bound** to exclude the
counterexample corner `(gap=0,debt=0)` via `|gap|_∞·|debt|_p=1`. The fragment IS a worked certified archimedean
rapidity separation (`>1/25`). **Transfer:** the artanh engine (P1) supplies exactly the "missing archimedean
floor" klein-S402 named. **Action:** for a small config, pair a p-adic covering-debt bound with an artanh-
certified archimedean rapidity floor and check the product formula closes the corner.

## T4 — odd-harmonic certification for the density route's Q_s and the tent loneliness integral
THM-729 (klein-S280): the density 2nd moment `Q_s=Σ_ℓ|U_s(ℓw)|²/ℓ²` with **odd-harmonic `sin(πℓ/7)/ℓ²`
weights**; closing the density row needs `Q_s=o(r²)` (a width-weighted 2nd-moment estimate). The loneliness
measure itself is `∫`(products of tent `‖v_i x‖`), whose Fourier support is odd-harmonic — the same class as
`2artanh`'s `Σ_{k odd}t^k/k`. **Transfer:** the snippet's truncate+geometric-tail certification gives float-free
two-sided bounds on exactly these odd-harmonic tails. **Action:** apply the sandwich to certify a clean upper
bound on the `Σsin²(πℓ/7)/ℓ²` factor and, more ambitiously, a certified `Q_s` tail.

## T5 (hint, not method) — 1/25 = 1/(2n−1) as a clean explicit target beating 1/(2n)
The fragment's target `1/25=1/(2·13−1)` beats classical `1/26=1/(2·13)`, and ML-values live in the exact family
`s/(ns+1)` (Kravitz 1912.06034: `ML(1,…,n−1,ns)=s/(ns+1)`, maxima at pair-sum times `m/(v_i+v_j)` — our
pair-sum machinery). **Insight:** a *rational-floor* certificate (P2) targeting the next Ostrowski/`s/(ns+1)`
rung is a natural intermediate milestone between our `1/13`-margin bounds (THM-2053: `M(v)≥1/13−R/2N`) and the
conjecture `1/14`. Frame explicit-floor pushes as "reach the next `s/(ns+1)` rung," certified by the artanh engine.

## Ledger
Strongest bet: **T2** (log-energy as an efficient complement to the moment-LP) — genuinely untried, testable.
Most actionable now: **T1** (Lean brick). Best "revive": **T3** (adelic archimedean floor). All rest on the
`91`-arithmetic bridge tying the snippet's weight to LRC(14). → THM-729/731 (2nd-moment vein), THM-661/opus-S185
(moment-LP), THM-252 (rapidity), THM-2052/2053/965 (`91`/`637` constants), opus-S2 (log engine), klein-S404
(exact form). Next: run the T2 experiment; land the T1 lemma.
