# LRC(14): the ρ\* density floor (OPEN-Q-108 / THM-527 part G) is **false** — inf ρ\* = 0

**Author:** opus session, 2026-06-29 (human prompt: "attack the real gap").
**Status:** VERIFIED (exact rational arithmetic + independent 200k grid). This is a
**negative result about a proof route**, not a statement about LRC(14)'s truth.
**Artifacts:** `04-computation/lrc14_rho_star_floor_refutation_opus_20260629.py`,
`05-knowledge/results/lrc14_rho_star_floor_refutation_opus_20260629.out`.

## The object (THM-527)
For `S = P ∪ L` (small part `P ⊆ {1..13}`, cluster `L`, co-offsets `E = Vmax − L`):
```
ρ*(P,E) = meas{ x∈[0,1) : ‖p x‖ ≥ 1/14 ∀p∈P  AND  maxgap{frac(e x):e∈E} > 2/7 }.
```
THM-527: `ρ*(P,E) > 0 ⟹ M(S) ≥ 1/14` (sufficient). The route (OPEN-Q-108 / part G
item 1) hoped to prove `inf ρ* ≥ c₀ > 0` over bounded-spread shapes, which would
give LRC(14).

## What is verified
1. **Canon arithmetic reproduced exactly:** `μ_3=1, μ_4=19/21, μ_5=9/14,
   μ_13=829/4620`; consecutive-cluster floor `= 1/84` at `k=9, P={1,2,3,12}`.
   The mechanism of `1/84` is a single window-0 sliver pinned between the `p=1`
   edge `1/14=12/168` and a `p=12` sub-danger `1/12−1/168=13/168` (width `1/168`
   per side, ×2).
2. **Unconditional partial floor (proved):** `ρ* ≥ μ(E) − |P|/7` (union bound,
   each `Bad_p` has measure `1/7`). Since `|P| = 13−k`, this gives `ρ* > 0`
   for **k ≥ 12** with no conjecture (k=12: `≥ 389/6468 ≈ 0.060`; k=13: `μ_13`).
3. **REFUTATION — `inf ρ* = 0`:** there are **admissible primitive covering
   13-sets** (a multiple of every `q∈{2..14}`) with `ρ* = 0`. Explicitly:
   - `S = {1,2,3,12,18,20,21,22,23,24,25,26,28}` (P={1,2,3,12}, E=[0,2,3,4,5,6,7,8,10]);
   - `S = {1,2,3,6,18,20,21,22,23,24,26,27,28}`;
   - `S = {1,2,3,18,20,21,22,23,24,25,26,27,28}`.
   Mechanism: at **spread ≥ 10**, window-0 has half-width `5/(7·spread) ≤ 1/14`,
   so it sinks entirely into `Bad₁`; then `{2,3}` cover the `½,⅓,⅔` windows.
4. **Not a counterexample to LRC(14):** all three sets are **lonely** —
   `M = 1/9, 5/41, 1/9` (all `> 1/14`), witnesses at `t = 17/27, 10/41, 11/45`.
   The witness denominators are **unrelated to Vmax=28** (consistent with THM-566:
   no uniform bounded-denominator witness). So `ρ*=0` is a structural **blind
   spot** of the criterion, not a failure of loneliness.

## Consequence for the proof program
The "lower-bound the density `ρ*`" strategy **cannot** prove LRC(14): its target
`inf ρ* ≥ c₀ > 0` is false on the very class (primitive covering sets) it was meant
to handle. The route is salvageable only by changing the certificate — e.g. taking
`max` of the good-period density **over admissible denominators `D`** (not just the
`Vmax` ruler), since loneliness of these sets lives at `D ∈ {27, 41, 45}`, not at
`Vmax`. The "attacker" structure (only speeds with `gcd(p,6)>1` damage the
low-denominator windows; coprime-to-6 speeds are gentle) explains both the small
consecutive floor and why the bounded attacker budget in `{1..13}` is the real
lever — and it is the same "small-prime" phenomenon flagged honestly in HYP-++2912.

## Recommended canonical integration (not done here, per "no silent canon override")
- Update THM-527 part G item 1 / OPEN-Q-108 status: **refuted as stated** (`inf ρ*=0`).
- Record `ρ* ≥ μ(E) − |P|/7 ⟹ ρ*>0 for k≥12` as a genuine unconditional lemma.
- Open the successor question: is `max_D ρ*_D(S) > 0` uniformly? (the multi-denominator criterion).

---

# Follow-up: the PRIME-DENSITY criterion (live successor to the refuted route)

**Artifacts:** `04-computation/lrc14_prime_density_criterion_opus_20260629.py`,
`05-knowledge/results/lrc14_prime_density_criterion_opus_20260629.out`.

**Criterion.** `N(S,p) = #{a∈Z/p : 14·min(r,p−r) ≥ p ∀s, r=(sa) mod p}`;
`N(S,p) ≥ 1 ⟹ M(S) ≥ 1/14` (sound: `t=a/p` is a witness).

**Verified (exact):**
1. **No blind spot.** Every `ρ*=0` set above is certified by a small prime witness
   (`N(S,19)=2`, `N(S,37)=2`, …). The defect that killed the Vmax-`ρ*` route is gone.
2. **Sidesteps THM-566.** For the divisor-loaded `S_B={1..11,13,84·lcm(1..B)}`, the
   density `N(S_B,p)/p` at a *fixed* large prime is **stable ≈ 0.010 for all B up to 100
   (a 43-digit speed)**. The unbounded *witness denominator* of THM-566 is irrelevant:
   at a prime coprime to the speeds the loaded runner is a generic unit, so the density
   tracks the fixed low-speed core.
3. **Exact, not a weakening.** `N(S,p)/p → lonely-measure(S) := meas{x:‖sx‖≥1/14 ∀s}`
   as `p→∞` over coprime primes, and `lonely-measure(S) > 0 ⟺ M(S) > 1/14`. So the
   criterion *equals* the target; it just makes it arithmetic over `Z/p`.
4. **Empirical floor.** `lonely-measure` over a 300-set covering sample: min ≈ **0.064**;
   the `ρ*=0` sets sit at 0.041–0.060; the resonant core `{1..11,13}` at 0.012. All `>0`.

**The clean remaining core (honest).** `N(S,p) = p − Σ|U_s| + Σ|U_s∩U_{s'}| − …`
(inclusion–exclusion) `= p·lonely-measure(S) + ` exponential-sum error. The union/first-
moment term is vacuous (`13·(1/7)=13/7>1`), so **positivity requires the second-order
cancellation** — bound the Gauss/Kloosterman-type error below the main term, uniformly
over covering 13-sets. This is the same analytic core as HYP-+2878 step 3 / "Node-3".

**Connection (git signal).** mac-mini's `HYP-3529/THM-578` (R-tail, 2026-06-29) proves a
**BV-Fourier uniform bound** `|R(M)| ≤ V_tot/12` by integration-by-parts over the same
`1/7`-sector structure. That is exactly an effective-equidistribution tool for the
exponential-sum error above. The prime-density criterion and mac-mini's R-tail bound are
two halves of the same analytic engine — worth joining.

**Net.** The `ρ*` route is dead (inf `ρ*`=0); the prime-density criterion replaces it with
a sound, blind-spot-free, THM-566-proof reformulation that is *equivalent* to LRC(14) and
lands precisely on the exponential-sum lower bound — the genuine hard core, now cleanly posed.

---

# Option (i): the exponential-sum analytic core — mapped and bottlenecked

**Artifacts:** `04-computation/lrc14_expsum_core_opus_20260629.py` + `.out`.

**Exact identity (verified).** For prime `p`, band `U={r:min(r,p-r)≤H}`, `H=⌊p/14⌋`,
`δ=|U|/p`, Dirichlet kernel `D_H(t)=Σ_{j=-H}^{H}e_p(-tj)`:
```
N(S,p) = Σ_{m=0}^{13} (-1)^m (1-δ)^{13-m} p^{1-m} Σ_{|T|=m} R(T),
R(T)   = Σ_{t_i≠0, Σ t_i s_i ≡ 0 (p)} Π D_H(t_i).
```
`m=0`: `p(1-δ)¹³` (the `(6/7)¹³` heuristic). `m=1`: `0`. **`m=2`: the dilate-overlap**
`(1-δ)¹¹ Σ_{i<j}(|a_{ij}U∩U| − |U|²/p)`, `a_{ij}=s_i s_j^{-1}`. The full expansion
reproduces `N` exactly (checked).

**Structural discovery — covering ⟹ forced resonance ⟹ non-truncation.**
A covering set must contain multiples of `2,…,14`, i.e. divisibility relations
`s_j | s_i`, giving **small integer ratios** `a_{ij}` and hence **persistently large**
band-overlaps `|a_{ij}U∩U|` at *every* prime. Consequence (verified): for the resonant
core, `main≈134, m2≈+105, rem(m≥3)≈−228`, net `N=12` — the `m≥3` terms are as large as
the main term and nearly cancel it. **So no absolute error bound (Cauchy–Schwarz, or a
naive BV/mac-mini majorant) can beat the main term in the hard regime** — positivity is
irreducible *cancellation*, not a boundable error. This is a sharp negative result: it
rules out the "main term wins" strategy and explains why every prior route stalled here.

**Payoff side — over-determination (verified, the live target).** The density of
witness primes `{p : N(S,p)≥1}` is **0.93–0.95 for every covering set**, including a
66-digit divisor-loaded speed; "bad" primes (bands cover `Z/p`) are a sparse minority.
A counterexample would need bands covering `Z/p` at *essentially every* prime —
forbidden by equidistribution of `{s_i/s_j mod p}` (Pólya–Vinogradov/Weyl). Proving the
bad-prime density is `<1` uniformly **is** LRC(14).

**Honest endpoint.** All three routes (Vmax-`ρ*`, prime-density, exp-sum) collapse to the
SAME irreducible statement, equivalent to the conjecture:
```
   lonely-measure(S) = meas{x : ‖sx‖ ≥ 1/14 ∀s∈S}  ≥  c₀ > 0   uniformly over
   primitive covering 13-sets    (empirically c₀ ≈ 0.012, attained near the resonant core).
```
What option (i) contributes: (1) the exact, verified exp-sum identity; (2) the dilate-
overlap form of the leading correction; (3) the proof that the obstruction is *forced
cancellation* (so error-bound strategies cannot close it); (4) the over-determination
target (positive witness-prime density) as the cleanest provable-looking successor. The
mac-mini BV-Fourier bound is the right tool for the *generic* tail (`m≥3` on non-resonant
pairs) but cannot by itself handle the resonant core — the two must be combined with a
cancellation (not absolute-value) argument.

---

# Option (i) -> over-determination: the bounded-core reduction & a 3-framework convergence

**Artifacts:** `04-computation/lrc14_overdetermination_opus_20260629.py` + `.out`.

**Part A - discretization structure theorem (PROVED).** For any 13-set `S` and prime `p`,
`max_a min_s ||sa/p|| >= M(S) - max(S)/(2p)` (`a=round(p t*)`, `||.||` is 1-Lipschitz). So if
`M(S)>1/14` then **every prime `p >= max(S)/(2(M(S)-1/14))` is a witness**, and the **bad
primes (`N=0`) are confined to `p < P0(S):=max(S)/(2(M-1/14))`** - verified (core: all bad
primes <=73 < P0=588). This makes "find a witness prime" effective and isolates the hard
primes as the *small* ones (the finite-check piece).

**Part B - the hard case is the BOUNDED RESONANT core (corroborated here).** Replacing
small speeds by far/large ones strictly **raises** `M` (`0.125->0.176->0.220->0.289`): far
elements *help*. The min-`M` covering sets are bounded and resonant (consecutive /
shared-factor); divisor-loaded sets reduce to the *same* `M` as the bounded core. The
*measure* can be pushed lower (to ~0.006 at an "almost all multiples of 3" set) but `M`
stays >=~0.083 - so the right uniform target is the **margin `M-1/14`**, not the measure.

**Convergence - three independent frameworks agree on the architecture:**
- this work (prime-density / exp-sum): hard case = resonant bounded core (`m>=3` tail is `Theta(p)`); finish needs signed cancellation, absolute bounds fail.
- **HYP-3131** (mac-mini, Lee-Yang/Asano): hard case = bounded core; far pushes zeros out; finish = `rho_bounded>1 => correlation inequality`.
- **kps-S31as** (Erdos-Turan/Fejer): hard case = AP/comonotone worst case; finish = signed Fejer/Erdos-Turan equidistribution.

All three independently conclude: (1) reduce to the **bounded resonant core** (far speeds
only help); (2) **no absolute/algebraic certificate** (Lee-Yang, Asano, SOS, moment-LP,
Bonferroni, Cauchy-Schwarz, BV) can cross the worst case - only **equidistribution /
signed cancellation** rescues it. Option (i) supplies the exp-sum *reason* (`Theta(p)`
resonant tail); kps supplies the *tool* (Fejer kernel = Erdos-Turan majorant); Part A
supplies the *effective bad-prime confinement*.

**The single remaining lemma (= LRC(14), now isolated):** a signed Erdos-Turan/Fejer
bound giving `lonely-measure(S) > 0` (equiv. `N(S,p)>0` at some prime) for the **bounded
resonant core** - the far structure and the small primes are handled (Parts A,B). A
concrete analytic-number-theory target, not another reformulation. **LRC(14) remains open.**

---

# Option (a): REFUTATION of HYP-2607 ("AP maximizes L_y") at k=12, 13

**Artifacts:** `04-computation/lrc14_AP_extremality_refutation_opus_20260629.py` + `.out`.

Working on the signed-equidistribution finish led to mac-mini's moment-LP route
(THM-534/537): `meas(S7(E)) <= L_y(E)=U_t(E)` is PROVED (dual certificate), and the
OPEN "finishing conjecture" HYP-2607 is that the **consecutive/AP cluster MAXIMIZES
`L_y` over all shapes E** (so `L_y(E) <= L_y(AP) <= cap_k`). It was verified
*exhaustively only at k=8* (AP the unique maximizer, 1716 shapes).

**My LP is validated exactly** against their anchor: `U4(consec_8)=2633/7350=0.358231`
(and `U2=0.4964, U3=0.4125` match THM-537's table).

**FINDING (rigorous, validated): HYP-2607 is FALSE at k=12 and k=13.**
| k | order t | U_t(AP) | max_E U_t | maximizer E | #shapes > AP |
|---|---|---|---|---|---|
| 12 | 2 | 0.7142 | **0.7330** | {0,2,3,4,5,6,7,8,9,10,12,14} | 14 |
| 13 | 2 | 0.7489 | **0.7619** | {0,1,..,10,12,15} | 13 |

The violation holds at **every moment order t=2,3,4** (t=2 is the order THM-537 uses
for k>=11). The maximizers are *near-AP* shapes (one shifted/added endpoint) -- exactly
the kind of perturbation that also broke the consecutive `1/84` floor (rho*) and the
THM-527 part-F consecutive extremality. The "AP is extremal" intuition fails at the
n=13 frontier in three different functionals now.

**Impact (honest, limited -- LRC(14) NOT threatened).** Like the rho* refutation, this
kills a clean *formulation*, not the conjecture. The binding case is small k (k=8,
exhaustively verified, AP unique max); at k=12,13 the cap margins are large, so
`L_y <= cap_k` very likely still holds -- but now via a **non-AP maximizer**, so the
moment-LP route's reduction `max_E L_y(E) = L_y(AP)` is invalid at k=12,13 and those
cases need re-verification (`max_E L_y(E) <= cap_k` directly, with the correct maximizer).
The proposed PROOF of HYP-2607 -- an "AP-orbit majorises the empty-sector moments" /
three-distance rearrangement (THM-537 D, kps-S31as Fejer finish) -- **cannot succeed as a
universal statement** and must be restricted to the binding small-k regime.

**Essential follow-up:** compute `cap_k` for k=12,13 and verify the non-AP maximizer's
`L_y` (0.733, 0.762) is still `<= cap_k`. If yes, LRC(14) survives at those k via a
patched argument; if no, the moment-LP loneliness bound itself fails there. **LRC(14)
remains open.**
