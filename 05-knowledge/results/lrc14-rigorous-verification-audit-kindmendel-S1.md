# Rigorous verification audit of the LRC(14) "witness-floor" proof

**Author:** kind-mendel-2026-06-22-S1
**Task:** rigorously verify the 14-runner Lonely-Runner proof the project has converged on.
**Scope:** the current sector / witness-floor route (THM-523/525/526/527/530/534/563,
the `LRCFourteenSkeleton` Lean DAG, kps-S32's "clean inequality", the legC synthesis).

---

## VERDICT (one line)

**LRC(14) is NOT proved.** The project has a sound, well-structured *conditional* reduction
plus overwhelming computational evidence, but at least **three analytic nodes remain open**.
This audit agrees with every relevant THM-file status line ("LRC(14) NOT proved") and finds
the recent message-level framing "VERIFIED-CLOSED / the clean inequality is closed"
(kps-S32, SESSION-LOG) to **overstate** the rigorous status. The conjecture itself is in no
computational doubt (true witness floors are ~5–7× the claimed threshold; see §5), but a
*rigorous uniform* floor is genuinely missing.

The foundation is legitimate: LRC for ≤12 moving runners is now a real theorem
(Sungkawichai–Trakulthongchai, arXiv:2604.23906, "Eleven, twelve, and thirteen lonely
runners"; building on Rosenfeld–Trakulthongchai for k≤9). So **LRC(14) = 13 moving runners
is genuinely the first open case**, and 14 = 2·7 being composite is exactly why the
polynomial method of those papers does not reach it. (Verified by web search this session.)

---

## 1. The logical skeleton (what the proof actually is)

Notation (THM-527/530/534): a counterexample must be a primitive covering 13-set `S = P ∪ L`,
`P ⊆ {1,…,13}` ("small part"), `L` the large cluster, cluster co-offsets `E` with `0∈E`,
`k=|E|`, `|P|=13−k`. `G_P = {x: ‖px‖≥1/14 ∀p∈P}` (safe set), `cap_k = min_{|P|=13−k} meas(G_P)`,
`p0(E) = meas(S7(E))` = measure of `x` where all seven 1/7-sectors are hit by `{frac(e_i x)}`,
`GOOD = {x: circular max-gap of {frac(e_i x)} > 1/7}`, `witnessG2 = meas(GOOD ∩ G_P)`,
`M(S) = max_t min_{v∈S}‖vt‖` (LRC(14) ⟺ `M(S) ≥ 1/14` for all primitive `S`).

The Lean assembly (`Verify.lean:942`, `lrc14_..._from_p0_positive_wide_bound_split_nodes`)
concludes `LRC14Statement` **from these hypotheses**, none of which it proves:

| node | statement | status |
|------|-----------|--------|
| reductions | counterexample ⊆ covering sets; `S=P∪L`; pigeonhole closes `k≤7` | **PROVED** (THM-523/525), elementary |
| `hbonf` | `nu+measGP−1 ≤ witnessG2` (Bonferroni) | **PROVED** (`LRCBonferroniMeasure`, sorry-free) |
| event side | `coverSet^c ⊆ goodSet`, measurability, `slowμ` prob. measure | **PROVED** (sorry-free) |
| **hp0cap** | `p0(E) ≤ cap_k − δ`, `δ>0`, for **all** E | **OPEN** (THM-534) — Gap 2 |
| **hmeasGP** | `cap_k ≤ meas(G_P)` | PROVED by *definition* of `cap_k` as the min (but see Gap 3) |
| **hpartA** | `0 < witnessG2 ⟹ M ≥ 1/14` | **OPEN** (THM-527 Part A) — Gap 1 |
| **herror** | `|finiteRho − witnessG2| ≤ arcCount/vmax` (finite-V_max ruler) | **OPEN / assumed** — Gap 1 |

So the Lean development proves an **implication**, not LRC(14). `rhoStar`/`shapeOf` are
`opaque`; `thm527_partA_density_pos_implies_reach` is an explicit `Prop` *obligation*. There is
**no** unconditional `theorem _ : LRC14Statement` anywhere (grep-confirmed this session).

---

## 2. What is rigorously PROVED (and I re-verified)

- **q-witness covering reduction (THM-523).** Elementary and correct: if `S` omits a multiple
  of some `q∈{2,…,14}` then `M(S) ≥ 1/q ≥ 1/14` at `τ=1/q`. A counterexample must be covering.
- **Easy-dominates-hard reduction (THM-525)** to a uniform `meas(G_C)` floor (OPEN-Q-108),
  using the genuinely-proved LRC(≤13) as the easy core. Non-circular.
- **Pigeonhole for small clusters:** `k≤6` unconditional, `k=7` a.e. (`maxgap ≥ 1/k > 1/7`).
- **Bonferroni** `meas(A∩B) ≥ measA+measB−1` (sorry-free) and the **event-side concreteness**.
- **THM-534 identity:** `meas(S7(E)) = p0 = Σ(−1)^{|A|}J(A,E)` and the dual `p0 ≤ L_y(E)`.
- **THM-563 periodicity:** `w·Δ_w` is exactly periodic (period `7·lcm(B)`) — genuinely proved.
- **The published exact constants.** I re-derived from definitions, independently
  (`04-computation/lrc14_independent_audit*_kindmendel.py`), and reproduced **exactly**:
  the per-size minima of `meas(G_P)` (6/7, 66/91, 55/91, 1979/4004, 2243/5880, 3029/10780 for
  |P|=1..6) and `p0(consec_8) = 481/1470`. The computational engine and the canon values are
  trustworthy.

## 3. What is only computationally VERIFIED (not proved)

- **`p0(E) ≤ cap_k` (the cap bound).** I re-ran the THM-534 search independently: over all
  **11,432** primitive bounded-spread clusters at k=8 (spread ≤16; their "11440"), all of k=9,
  k=10, and k=8 spread ≤20 (77,400 sets), the maximizer of `p0` is **consec** and `p0 < cap_k`
  with **zero** sets over cap. This corroborates THM-534's central claim — *over bounded spread*.
- **The legC three-piece structure** (k=9..12 genuine-wide doublets): careful, but rests on
  far-count monotonicity (HYP-2803, verified not proved), an **empirical** period-max bound
  (≤1.74; the rigorous Tornheim fallback only gives `M*≤48`), and finite checks (some still
  running at write time). It does **not** cover k=8 (binding) or k=13.

## 4. The OPEN nodes (the actual gaps)

**Gap 1 — Part A / slow-fast finite-V_max (hpartA + herror).** `witnessG2 > 0 ⟹ M ≥ 1/14` is
proved only in the `V_max → ∞` *limit*. The finite-V_max correction (the `arcCount/vmax`
ruler error, assumed in Lean) is unproved and is *most delicate exactly at the
boundary-achieving hard core* `{t,2t,…,12t,V}` and CRT relatives, where (kps-S4 reflection)
"the limit object achieves the boundary exactly" and "no compactness argument can reach" the
needed ε. This is the deepest gap.

**Gap 2 — uniform cap bound `p0(E) ≤ cap_k − δ` for ALL E (hp0cap).** THM-534 reduces it to
the scalar extremality **"consec maximizes `L_y(E)`"**, which is OPEN (verified only on
bounded-spread windows). The unbounded-spread / large-cluster region is handled only by the
legC computation with its own residuals (Gap 2 ⊇ legC's open items). k=8 and k=13 are not
closed by legC at all.

**Gap 3 — a genuinely uniform witness floor.** `hmeasGP` (`cap_k ≤ meas(G_P)`) is true by the
definition of `cap_k`, but the *clean* analytic lower bound on `witnessG2` that the proof can
actually invoke is Bonferroni, `witnessG2 ≥ cap_k − p0`, which is **lossy** (see §5). Getting
the true (large) floor uniformly needs the unproved resonant-nbhd-width lemma (kps-S32) or an
equivalent sharper-than-Bonferroni measure bound — itself ≈ as hard as Gap 2.

---

## 5. Reconciling the court case, its withdrawal, and kps-S32 (new, independently computed)

`02-court/active/CASE-p0-route-insufficiency.md` (WITHDRAWN) and kps-S32 appear to conflict.
I resolved this by computing the **true** `witnessG2` for the binding case `consec_8` with the
`cap_8`-achieving `P={1,5,7,8,9}` (`04-computation/lrc14_independent_audit4_kindmendel.py`,
exact rationals + grid):

```
meas(G_P)               = 2243/5880  = 0.381463   (= cap_8)
meas(coverSet ∩ G_P)    = 3079/35280 = 0.087273   <- the TRUE intersection
meas(coverSet^c ∩ G_P)  = 10379/35280= 0.294189   <- RIGOROUS lower bound on true G2
true G2 (grid, maxgap>1/7) ≈ 0.381433             <- nearly all of G_P is good
Bonferroni floor cap-p0 = 319/5880  = 0.054252    <- the clean bound (< m_P)
m_P                     = 14249/252252 = 0.056487
```

So **all three parties are partially right**:
- The court case is correct that the *clean Bonferroni* floor `cap−p0 = 0.0543 < m_P` at k=8.
- Its withdrawal is correct that only **positivity** is logically needed (and holds: I confirm
  `cap_k − p0(consec_k) > 0` for all k=8..13; `≥ m_P` fails *only* at k=8).
- kps-S32's *numerical* claim "true `G2 ≥ m_P` (several×)" is **correct**: the true floor is
  ≈0.29–0.38, i.e. **5–7× m_P**, even at the binding case.

The mechanism: Bonferroni bounds `meas(coverSet ∩ G_P)` by the *whole* `p0 = 0.327`, but the
true intersection is only `0.087` (most of `coverSet` lies outside `G_P`). That ~0.24 of slack
is the entire "k=8 failure." **The conjecture is not in doubt; the loss is in the only clean
bound we can currently prove.** kps-S32's resonant-nbhd argument would recover the true floor —
but its uniform "nbhd-width lemma" is exactly the unproved Gap 2/3 (kps-S32 itself flags
"REMAINING = the analytic RIGOR of the resonant nbhd-width lemma").

---

## 6. Bottom line / handoff

LRC(14) status: **first open case (real), strong conditional reduction, three open analytic
nodes, computationally overwhelming, NOT proved.** Highest-leverage targets, in order:

1. **Gap 1 (Part A finite-V_max).** Make `witnessG2 > 0 ⟹ M ≥ 1/14` rigorous at finite V_max,
   especially on the boundary-achieving core. This is the conceptual crux and the least
   addressed; everything else feeds a bound whose payoff is gated here.
2. **Gap 2 ("consec maximizes `L_y`").** A single scalar moment-functional extremality
   (cf. THM-531 AP-orbit invariance, THM-533 W-monotonicity). Closing it gives `p0 ≤ cap` for
   all E unconditionally and subsumes most of legC.
3. **Messaging:** stop labeling the floor "VERIFIED-CLOSED." It is "computationally verified;
   rigorous uniform lemma open." The canon THM files already say this; the SESSION-LOG should
   match.

Files: `04-computation/lrc14_independent_audit{,2,3,4}_kindmendel.py`,
`05-knowledge/results/lrc14_independent_audit_kindmendel_check{12,34,4,5}.out`.
Refs: THM-523/525/527/530/534/563, OPEN-Q-108, HYP-2847, kps-S32, legC synthesis,
CASE-p0-route-insufficiency, arXiv:2604.23906.
