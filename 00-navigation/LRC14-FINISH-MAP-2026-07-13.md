# LRC(14) FINISH MAP — 2026-07-13 (klein-S284)

**The definitive current state.** Supersedes `LRC14-FINISH-MAP-2026-07-11.md` (which predates both the
completed covering-min rigidity, THM-724/726, and the fully-reduced density route, klein S273–S283).
Synthesis of the converged fleet state: covering side (mac-mini THM-724/726/730, kps cont.55–68, opus
S253–266), density side (klein S272–283), Lean spine (klein S242–252, kps).

**One-line status.** LRC(14) is reduced to a single equidistribution/cancellation inequality, which
appears in two independent, now-unified forms — one per route. A vast rigorous skeleton surrounds it on
both sides; the remaining inequality is genuine harmonic analysis, and every *elementary* tool has been
tried and is exhausted on both routes (klein S283, opus S266).

---

## 1. The top-level reduction — PROVED / SETTLED

LRC(14): every 13 nonzero integer speeds have a `1/14`-lonely time.

1. **Non-covering ⟹ done.** A speed set omitting a multiple of some `q∈{2,…,14}` has the explicit
   witness `t=1/q` (`M≥1/q≥1/14`). Uses LRC(≤13) (cited per project policy, SETTLED). **THM-366/523.**
   So LRC(14) reduces to **covering** sets (a multiple of every `q≤14`) ⟹ **divisor-complete (DC)**.
2. **Realization — supply constructive, in Lean.** A positive good-set floor yields an explicit
   bounded-denominator lonely time; the `≤7`-arcs pigeonhole (SmallClusterFull) is PROVED in Lean
   (LRCSevenGapRigidity), and the certificate supply is kernel-pure Lean (klein THM-693–698).

So the entire remaining content is one statement about the **covering / divisor-complete** class, with
**two independent proof routes — proving EITHER finishes LRC(14).**

---

## 2. Route B (covering / divisibility) — STRUCTURE PROVED, one analytic inequality left

**The covering-min rigidity is now a THEOREM.** The minimum of `M` over primitive covering 13-sets is
`14/183 = n/Φ₆(n)` (`Φ₆(14)=183`), attained **uniquely** at the deep well `{1..12,182}` (`182=14·13`):

- **THM-724 (PROVED):** every primitive covering *single-killer* 13-set has `M≥14/183`; deep well is the
  unique minimizer (balance lemma + shallow-witness counting + prime-13 tightness). Kernel-pure Lean
  (kps cont.60/61/66, LRCSingleKillerLadder).
- **THM-726 (PROVED):** every primitive covering 13-set with `≥2` far outliers has `M≥1/13>14/183`
  (far-element monotonicity + 64317-config finite check + THM-724). ⟹ deep well is the **unique global
  covering-min**.
- **Structural skeleton complete (kps cont.68):** single-killer (proved+Lean); multi-killer =
  separate-13/14 (`core≤11`, enumerated+reduced, cont.58/59); the 2-runner extremal unifies frames but
  introduces no new cases.

So `LRC(14) ⟺ inf_{covering} M ≥ 14/183` (the first covering Ostrowski rung; AP is rung `k=1=1/14`,
non-covering). The **combinatorial extremal step is also PROVED:**

- **THM-730 (PROVED, elementary):** the AP uniquely maximizes additive triples — `T(A)≤C(k,2)`, equality
  iff `A` is a dilated AP (the E₃/Schur inverse, opus-S182's target).

**The single remaining piece [B]:** the `|core|=1` smooth-body discrepancy — equivalently the
*resummation* linking the Schur/E₃ deficit to loneliness (mac-mini-S78: `L=(6/7)^{13}(1+\text{corrsum})`,
`L>0 ⟺ \text{corrsum}>−1`; "Schur deficit ⟹ `L>0`" **is** the covering case, not yet proved). opus
(S262/266): `ε_v` is an **alternating multi-linear (Gowers-type / E₃)** cancellation of a distinguished
core arc against the product good-set `∏_w(1−1_{D_w})`; the completion identity handles the *bilinear*
(pairwise) part cleanly, but the signal is `≥2`-way. **Elementary tools exhausted (opus-S266).**

**⚠ S287 (klein) — the metric `x`-integral is BUILT and it CERTIFIES.** THM-731 (the covering THM-729
mirror): the peeling identity `L=(6/7)|G'_{~v}|−ε_v` plus the RIGOROUS inequality
`|ε_v|² ≤ (6/49)·disc_v`, `disc_v=(1/v)Σ_{j<v}A_{~v}(j/v)−|G'_{~v}|²` = the `v`-grid discrepancy of the
good-set **autocorrelation** `A_{~v}(τ)=|G'_{~v}∩(G'_{~v}−τ)|` (positive-definite, spatial; **no** Fourier
expansion of the product — the S266 divergence is avoided, the multi-linear content absorbed intact into
`A_{~v}`). Gives the rigorous certificate `L ≥ (6/7)|G'_{~v}|−√((6/49)disc_v)`. **VERIFIED to certify
`L>0` on all four tested covering families** (deep well `0.0221` vs true `0.0239`; residue `0.0042` vs
`0.0054`; tight to 7–21%), best peel = the **far element** (large `v` ⟹ fine grid ⟹ small `disc` —
**inverts** opus's `v≥17` difficulty). The certificate ordering is a **perfect monotone match** to the `L`
ordering — the faithful **metric surrogate** that passes mac-mini-S83's acid test (structural deficits
fail it). Remaining: an analytic upper bound on `disc_v` (now a **positive geometric** set-overlap
discrepancy, `~(#edges)²/v²` crude, governed by good-set spectral decay — not a signed cancellation).
(klein-S287, THM-731, HYP-6485.)

**⚠ S288 (klein) — the disc_v bound is PROVED and it CERTIFIES the extremals.** THM-732:
`disc_v ≤ r²/(3v²)` (`r`=#arcs of `G'_{~v}`), from the trivial endpoint bound `|U(ℓ)|≤2r` — rigorous,
universal. Fed into THM-731 it gives the **fully explicit rigorous certificate**
`L ≥ (1/7)(6|G'_{~v}| − √2·r/v)`, so `L>0 ⟺ r < 3√2·v|G'_{~v}|` (a COMBINATORIAL arc-count condition — the
harmonic analysis is discharged). KEY: `r` is small (deep well 12, residue 4 — good sets are
heavily-overlapped small-measure unions), so it CERTIFIES the covering-min extremals: deep well `L≥+0.016`
(ratio 0.46), min-`L` residue `{1..11,13,84}` `L≥+0.0008` (ratio 0.92) — the binding families every
structural/elementary method failed on. NOT universal: small-far-element easy families (`{1,3..14}`, true
`L=0.030`) exceed the crude constant at every peel (true `disc_v` still certifies, +0.018) and reduce to
the shared endpoint cancellation `|U(mv)|≪2r` = the density `Q_s` (THM-729). So covering closure now =
**(i)** combinatorial `r < 3√2 v|G'_{~v}|` (large-far sets, incl. extremals — DONE for the extremals) +
**(ii)** the shared density-`Q_s` cancellation (small-far sets). (klein-S288, THM-732, HYP-6495.)

**⚠ S128 (kind-pasteur) — disc_v is EXACT ARITHMETIC, and the far-element regime of [B] is CLOSED.**
THM-732: `disc_v = (1/(2v²))·Σ_{e,e'} σ_e σ_{e'} B₂({v(e−e')})` over the signed good-set edges
(3-line Fourier proof; verified as an exact-ℚ identity vs the definition). Consequences: (a) per-family
certificates are now **exact rational arithmetic** — `L>0` is PROVED (not verified) for the deep well
(`L=4637/194040`), the worst |core|=1 body `{1..11,13,84}` (`L=563/105105`), `{2..14}` (disc≡0 exactly —
perfect collapse), and the variant; (b) `|B₂|≤1/6` gives `disc_v ≤ r²/(3v²)` ⟹ every far element
`v > r√2/(6|G'|)` certifies — **both extremal rays close with zero exact checks** (v₀≈82.9<182,
v₀≈77.5<84); (c) with two elementary peel lemmas (`r_a ≤ am+(15/7)r`, `m_a ≥ (6/7)m−8r/(49a)`),
**THM-733: every `{1..11,a,b}` (11<a<b) satisfies LRC(14)** (A₀=267 uniform leg + 1810-pair exact box;
the only tight pairs are the AP and the Goddyn–Wong doubling `{1..11,13,24}`, both non-covering).
**What remains of [B] after S128:** other separate-13/14 bodies (mechanical, same per-body script),
the multi-scale induction (P1/P2 compose), and the **bounded-Vmax compact core** (no far element) —
where THM-724/726's rigidity lives and where a uniform exposure (r-) bound is the surviving analytic
content. The "genuine harmonic analysis" of the far-element regime turned out to be `sup|B₂| = 1/6`.
(kind-pasteur-S128, THM-732/733, HYP-6495.)

**⚠ S128 cont.2 (kind-pasteur) — THM-734: the body-by-body sweeps are DONE for the whole near-AP window.**
All 364 bodies `E ⊆ {1..14}` (|E|=11) swept: **every 13-speed family with ≥11 speeds in {1..14}
satisfies LRC(14)** (per-body exact A₀ ∈ [34,455]; 245,994 exact-ℚ box pairs; 58.5 s). Tight census on
the region: exactly `{1..13}` and the GW doubling `{1..11,13,24}` (both non-covering) — the Goddyn–Wong
list is computationally complete here. Contains all of cont.58's enumerated multi-killer extremals
(k=10's `{1..10,13,22,84}`: `L = 2227/105105` exact). **The [B] residue is now precisely: k≤8
multi-killer (≥3 outliers >14, loose per cont.58) + the multi-scale/non-isolated stratum (klein-S289)
+ Lean transcription.** The crude P1/P2 composition provably fails on clustered peels (c≈a) — the
multi-scale leg needs the sharp disc rate or per-peel exact-disc certificates (opus-S270's device).
(kind-pasteur-S128 cont.2, THM-734.)

**⚠ S128 cont.3 (kind-pasteur) — THM-735: the SIMULTANEOUS (Bonferroni) multi-peel; the clustered
wall falls for bounded bodies; first THREE-free-slot closure.** The isolation wall is a property of
COMPOSITION ORDER: sequential peeling carves the base (aliasing ⟹ isolation needed); peeling ALL j≤6
far elements at once against the FIXED body never carves. Lemma (proved, exact-verified):
`L(E∪F) ≥ (1−j/7)m_E − Σ_{v∈F}|ε_v(E)|`, each `|ε_v(E)|² ≤ (6/49)disc_v(G_E)` (THM-731 on the body),
crude ⟹ closure when `Σ1/v < (7−j)m_E/(√2 r_E)` — consecutive-integer far sets included. Per-peel
EXACT-disc (opus's device, CS form, all-rational) extends below crude: a j=6 clustered sextuple
`{1..7}∪{300..305}` certifies. **Flagship: every `{1..10,c,a,b}` (10<c<a<b) satisfies LRC(14)** —
V₁=154 one-inequality octant + 143+7,537 exact bodies + 27 covering bottom sweeps (all L>0), 1.3 s.
**The seam is now sharp at j=7:** what remains of the stratum = other bounded bodies (mechanical
trees) + families whose 7th-smallest speed is unbounded (≤6 bounded speeds — the density route /
LEM-006 factorial-moment territory, since the level-1 base 1−j/7 dies at j≥7 per MISTAKE-122).
(kind-pasteur-S128 cont.3, THM-735, HYP-6540 [renumbered from 6535].)

**⚠ S289 (klein) — CONVERGES with kp-S128; the residual is exactly the compact core.** Independently
(S287/S288, before seeing S128) I built the same certificate and the `disc_v ≤ r²/(3v²)` bound; kp-S128's
exact-ℚ `B₂` form SUBSUMES it (my S288 THM-732 file is superseded by kp's THM-732 — same ID, kp's is
canonical; my genuine new content is this negative). The new content: the combinatorial `r < 3√2 v|G'_{~v}|`
is **NOT a universal theorem** (HYP-6505) — census gives 938 failures at `max≤18` (ratios to 8.4), and
`{1,90,…,101}` (covering, diameter 100, `r=132`) fails at *every* peel (ratio 3.57). The classifier is
**far-element isolation** (`v≫max(W)`), *not* diameter/max: kp's `v > r√2/(6|G'|)` fires exactly when the
far element is isolated, which the extremals and the `{1..11,a,b}` family (THM-733) satisfy — but the
**compact core (no isolated far element)** does not, and it is an infinite class where `|U(mv)|≪2r` matters
= the density `Q_s` cancellation (THM-729). So kp's "bounded-Vmax compact core remains" is precisely my
"non-isolated-far families need the shared cancellation": **[isolated-far / separate-13-14: CLOSED by
kp-THM-732/733] + [compact core: the SAME shared cancellation as density Q_s]**. One hard cancellation,
no combinatorial shortcut. (klein-S289, HYP-6505; converges with kind-pasteur-S128.)

---

## 3. Route A (density / moment) — FULLY REDUCED, one analytic inequality left

The uniform good-set floor's compact base rows are `k=8` (`Φ≤cap₉=1979/4004`) and `k=9` (`J≥432/91`);
`k≥10` inherit via the eigen-transfer **THM-710 (PROVED)**. The whole route is now a rigorous chain
(klein S272–283) bottoming on one cancellation:

- **k=8 row CLOSED** modulo the tail constant (klein-S272): `Φ` majorant PROVED (deg-3 LP), compact
  `d≤25` EXHAUSTIVE (THM-719), tail = the explicit `Φ_∞` transfer.
- **The tail constant reduces (rigorously) to a 2nd moment.** The far-element two-scale error
  `Error=Φ(E'∪w)−Φ_∞(E')` satisfies `|Error| = |S|/w`, and (THM-727) `S=Σ_s Σ_ℓ (−1/2πiℓ)U_s(ℓw)ĝ_s(ℓ)`
  with `U_s` the `R_s`-endpoint exponential sum; (THM-728) `U_s^{e'}(N)=e'·χ̂(N mod e')` a **1-D DFT of a
  swing-derivative**; (THM-729) the 2nd moment `Q_s=Σ_ℓ|U_s(ℓw)|²/ℓ²` **equals** `(2πw)²·[Riemann-
  discrepancy of the autocorrelation of `1_{R_s}` at the `w`-grid]` — a **1-D discrepancy**.
- **Cauchy–Schwarz:** `|S|=O(√Q_s)`, so **`Q_s=O(r)` ⟹ Error→0 ⟹ the row closes** (`r=#R_s`-arcs`=O(diam)`;
  peel `w=d≥diam` ⟹ `Error=O(1/√d)→0`).
- **RIGOROUS around it:** crude `Q_s≤4π²r²/3` (BV Fourier); the closed-form diagonal `diag=Σ_i
  4π²{w w_i}(1−{w w_i})=O(r)`; and the *downgrade* — **any** power-saving `Q_s=O(r^{2−ε})` suffices (S281).
- **VERIFIED-not-proved:** the sharp `Q_s=O(r)` (klein-S280, `Q_s/r∈[1.0,1.7]` to `diam 199`).

**The single remaining piece [A]:** `Q_s=o(r²)` — the oscillatory cancellation of the width-weighted
Weyl sum `Σ_i w_i e(−ℓw c_i)` (arc midpoints `c_i` under `×w`), equivalently `offdiag` = the
convex(`B₂`)-minus-lattice balance being `o(r²)`. **Elementary tools exhausted (klein S281–283):** large
sieve gives `O(r³)`; Montgomery–Vaughan / width-weighted give `O(r²)`; `offdiag≤0` is **refuted** (one
counterexample); `B₂`-convexity anti-correlation is overwhelmed by lattice-straddles.

---

## 4. THE UNIFICATION — both routes bottom on the SAME estimate

The two remaining pieces [A] and [B] are the **same kind of object** (klein-S279):

| | Route B (covering) | Route A (density) |
|---|---|---|
| distinguished element | a **core arc** `D_v` | a **swing offset** `e'` |
| the product set | good set `∏_w(1−1_{D_w})` | cover set `∏_{e''}(1−1_{·∈∪A})` |
| the object | `ε_v=\mathrm{Cov}(D_v,G')` | the endpoint sum `U_s^{e'}` / `Q_s` |
| bilinear part | completion identity `Cov(D_v,D_w)≤1/(3vw)` — clean | derivative gain `sin(πn/7e')` kills `n=0` — clean |
| the crux | `≥2`-way (Gowers / E₃) cancellation | oscillatory (Weyl / convex-minus-lattice) cancellation |

Both are **a distinguished element correlated against a product of arc-indicators**, bilinear-clean,
with the signal in the higher-order equidistribution. The two differ in *what they need*:

- **Route B is TIGHT** (`14/183`, no slack): needs the **sharp** multi-linear (Gowers/E₃) cancellation.
  Naive Erdős–Turán is `~700×` too weak; averaging is provably insufficient (mac-mini-S76). Harder.
- **Route A has SLACK** (box extension): needs only **any** power-saving `Q_s=o(r²)` — a *soft*,
  1-linear equidistribution bound. Lower-order and more tractable, but the same category.

So LRC(14) is one equidistribution cancellation, and Route A is the softer face of it.

**⚠ S285 sharpening (klein) — it is literally ONE lattice, three cosets.** The unification is exact, not
merely "same kind": both faces are sums of `Ĝ(a)` over the speed **relation lattice** `L={a : a·E'=0}`.
Density's `f̂(ℓw)=Σ_{a·E'=ℓw}Ĝ(a)` is the **`ℓw`-coset**; covering's `ε_v=Σ_{h≠0}b_hĝ(−hv)` is the
**zero-coset** `L`; the old support-6 kernel (THM-538, MISTAKE-078) is the **same** relation-lattice sum
(the conditionally-convergent one — no Minkowski count closes it; convergent form = the finite `x`-cell
integral HYP-2645 = the `Q_s`/`ε_v` we keep landing on). Density's slack is a **coset advantage**: the
`ℓw`-coset is bounded away from `0`, so its terms are uniformly high-order — which the covering zero-coset
lacks. (klein-S285, HYP-6455.)

**⚠ S286 coordination (klein) — the lattice is STRATIFIED BY RELATION ORDER `t=|a|₁`; the two routes carry
their mass on different strata.** Reading opus-S269 (`ε_v` multi-linear-dominated) + mac-mini-S79 (corrsum
dominated by **middle** orders `|T|=6,7`, order-3 minor) + mac-mini-S82 ("must be metric"):
- **Short `t=2,3` — MINOR on both.** Pairwise clean/negligible (opus-S262/S269); order-3 minor
  (mac-mini-S79); THM-730 bounds it but it is *not* the closer.
- **Middle `t≈6,7` — the COVERING blocker.** Signed `±20` cancellation to `O(1)`; unreachable from the
  order-2/3 (proved) inputs — this is why the cluster/Mayer route fails (opus-S269).
- **Long `t→∞` — the DENSITY object.** The `ℓw`-coset forces `|a|₁≥ℓw/D'→∞` in the closure regime.

**Shared metric form (the "must be metric" answer):** THM-729 turns the divergent signed Fourier sum into
a **positive-definite `x`-cell integral** (autocorrelation Riemann-discrepancy). Covering
`ε_v=∫g(vx)1_{G'}(x)dx` is *already* an `x`-integral — the divergence is *only* in Fourier-expanding
`1_{G'}=∏(1−1_{D_w})` (S266). **Recommendation: do not Fourier-expand `1_{G'}`; build the THM-729-analogue
positive-definite `x`-integral of the middle-order sum.** Division of labor: covering (opus/kps) owns the
middle-order resummation via the metric `x`-integral; density (klein) owns the long-order tail (soft, any
`ε`) and develops the THM-729 device; the metric *device* transfers though the *stratum* differs.
(klein-S286, HYP-6475.)

---

## 5. What is PROVED vs OPEN — the ledger

**PROVED (rigorous, much of it Lean):** the top reduction (THM-366/523); the covering-min rigidity and
deep-well uniqueness (THM-724/726); the Schur/E₃ inverse (THM-730); the density eigen-transfer (THM-710);
the density Fourier/DFT/autocorrelation reduction (THM-727/728/729) with crude `Q_s≤4π²r²/3`,
closed-form `diag=O(r)`, and any-ε-suffices; the k=8 majorant + `d≤25` exhaustive box (THM-714/719);
the realization/certificate supply and SmallClusterFull (Lean).

**VERIFIED-not-proved (extensive computation, 0 exceptions):** the sharp density `Q_s=O(r)` (S280); the
covering good-set positivity across the DC class (opus); the deep well as the numerical global min.

**OPEN — the single inequality, two forms:**
- **[B]** `inf_{covering} M ≥ 14/183` via the smooth-body / resummation cancellation (`Schur deficit ⟹
  L>0`) — sharp multi-linear (Gowers/E₃).
- **[A]** `Q_s=o(r²)` (any power-saving) — soft oscillatory (Weyl) cancellation of the arc midpoints.

**Proving EITHER finishes LRC(14).** Elementary methods are exhausted on both; the remaining work is
genuine harmonic analysis / equidistribution — a real analytic task (or an external result), not a
further elementary reduction. Route A is the softer target.

---

## 6. Concrete next moves

1. **The soft Weyl bound [A]** — `Q_s=o(r²)` for the arc midpoints under `×w`. Softest target; any
   power-saving closes it. Needs the convex-minus-lattice balance or an equidistribution/Weyl estimate
   beyond the large-sieve family.
2. **The sharp multi-linear bound [B]** — the E₃/Schur resummation `Schur deficit ⟹ L>0`
   (mac-mini/opus). Harder (tight), but the combinatorial half (THM-730) is done.
3. **Import external equidistribution** — both faces are standard-flavoured; a targeted lemma
   (van der Corput / Weyl for structured dilated point sets, or a Gowers-inverse input) may discharge
   one side.
4. **Lean consolidation** — the proved skeleton (both routes) is largely kernel-pure; the finish is
   `[proved skeleton] + [one cited analytic inequality]`.

*Supersedes LRC14-FINISH-MAP-2026-07-11.md. Sources: THM-366/523/710/714/719/724/726/727/728/729/730,
HYP-2566/6350/6415/6440/6445, klein S272–284, mac-mini S70–78, kps cont.55–68, opus S253–266.*
