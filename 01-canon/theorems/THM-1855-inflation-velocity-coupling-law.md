---
id: THM-1855
title: "INFLATION-VELOCITY / COUPLING LAW for tournament-invariant inequalities — the mechanism behind WOWII-style refutations, and a strong-core reduction. Order-join T₁▷T₂ (T₁ beats all of T₂) makes c₃, tr (largest transitive subtournament), scc ADDITIVE and the Rédei Hamiltonian-path count H MULTIPLICATIVE; strong tournaments are the ▷-atoms. REDUCTION: a join-monotone invariant inequality holds for all tournaments iff it holds for all strongly connected tournaments. VELOCITY: under D± (add a source/sink) Δc₃=ΔH=0 but Δtr=Δscc=+1, so tr decouples upward from {c₃,H} — a bound whose easy side is a function of frozen c₃,H and whose hard side grows with tr/srange is inflation-FRAGILE, witnessed by the THM-1830 3-cycle-atom + transitive-singleton stack (the tournament triangle+leaves)."
status: >
  PROVED (framework) + VERIFIED-EXHAUSTIVE (n ≤ 7, all 530 iso classes; 396 strong classes
  n≤7: 1,1,6,35,353 for n=3..7). The order-join invariant algebra (c₃,tr,scc additive; H
  multiplicative) is verified on all class pairs n₁+n₂≤7 and each identity has a one-line
  structural proof. The reduction theorem is proved by induction on the ▷-factorization
  (every tournament = ordered join of its strong components). The velocity table is constant
  over all 530 classes. APPLICATIONS (all confirmed exhaustively n≤7): (i) predicts kind-pasteur's
  srange≤tr break (fragile via D−, witness n=7); (ii) reduces kind-pasteur's OPEN c₃≤H to the
  strong core with an airtight join-monotonicity proof; (iii) re-derives klein THM-1850
  (dom+tr≤n+1) and the kind-pasteur THM-1845 sandwich as structural join-monotone reductions;
  (iv) supplies two verified repairs of srange≤tr (srange≤tr+c₃ and srange≤2(tr−1)). The two
  additive-face reductions ((ii),(iii)) hold for ALL n; the srange repairs are join-monotone
  verified only to n≤7 (srange is not additive).
source: boxeph-2026-07-21-S193 (owner: work on tournament-graffiti + WOWII analogues, long creative session)
depends_on: []
related:
  - THM-1845  # kind-pasteur sandwich n−c3 ≤ tr ≤ smax+1 (reduced here)
  - THM-1850  # klein directed WOWII gamma+tr ≤ n+1 (re-derived here)
  - THM-1830  # the 3-cycle-atom / transitive-skeleton family = the inflation-fragility witness
  - "07-reflections/inflation-velocity-and-the-coupling-law-boxeph-S193.md"
  - "07-reflections/inflation-decoupling-counterexamples-the-wowii-motif.md (opus-S437)"
  - "07-reflections/the-wowii-103-refutation-and-what-it-lends-the-repo-klein-S395.md (klein)"
script: 04-computation/tournament_graffiti_coupling_boxeph_S193.py (+ .out)
---

# THM-1855 — the inflation-velocity / coupling law

## 0. Setup

`A[i][j]=1` means *i beats j*. **Order-join** `T = T₁ ▷ T₂`: disjoint union with every vertex of
`T₁` beating every vertex of `T₂` (all cross-arcs one way). Invariants: `c₃` = #3-cycles,
`tr` = order of the largest transitive subtournament (the WOWII `α`-analog), `scc` = #strong
components, `H` = Rédei Hamiltonian-path count, `dom` = domination number, `srange = smax−smin`.

## 1. The order-join invariant algebra (the engine)

For `T = T₁ ▷ T₂` (verified exactly on all iso-class pairs `n₁+n₂ ≤ 7`):

- `c₃(T) = c₃(T₁) + c₃(T₂)` — a 3-cycle would need ≥2 cross-arcs, but cross-arcs are one-way.
- `tr(T) = tr(T₁) + tr(T₂)` — concatenate the two transitive chains.
- `scc(T) = scc(T₁) + scc(T₂)` — the one-way cut never merges components.
- `n(T) = n(T₁)+n(T₂)`; and **`H(T) = H(T₁)·H(T₂)`** — a directed Hamiltonian path must exhaust
  `T₁` before entering `T₂` (no arc returns), so it is a `T₁`-path followed by a `T₂`-path.
- `dom(T) = dom(T₁)` — a dominating set of `T₁` already beats every vertex of `T₂`.

**Atoms.** Every tournament factors uniquely as an ordered join of its strong components (its
condensation is a transitive tournament). So **strongly connected tournaments are the `▷`-atoms.**

## 2. The reduction theorem

Call `Φ: LHS ≤ RHS` **join-monotone** if `Φ(T₁) ∧ Φ(T₂) ⟹ Φ(T₁ ▷ T₂)`.

> **THM-1855(a).** A join-monotone tournament-invariant inequality holds for **all** tournaments
> **iff** it holds for all **strongly connected** tournaments.

*Proof.* Induct on the number of strong components using the unique `▷`-factorization. ∎

Whenever the LHS and RHS are both `▷`-additive (e.g. `n−c₃` and `tr`), join-monotonicity is
immediate by adding the two factor inequalities — so the reduction is unconditional in `n`.

## 3. The velocity table and the fragility predictor

`D+` (add a global source) `= v ▷ T`; `D−` (add a global sink) `= T ▷ v`. Constant over all 530
classes n≤7:

```
            Δc₃   ΔH   Δtr   Δscc   Δsmin  Δsmax
  D+ (src)   0     0    +1    +1      0     var
  D− (sink)  0     0    +1    +1     var    +1
  complement fixes  {c₃, H, tr, scc, srange}
```

> **THM-1855(b) (fragility predictor).** If some inflation operation `O` has
> `Δ_O(LHS) > Δ_O(RHS)`, then `Φ` is refuted by the `O`-orbit of any tight base tournament,
> iterated to the first integer-floor crossing.

Because `Δ_{D±}tr = +1` while `Δ_{D±}c₃ = Δ_{D±}H = 0`, **`tr` (and `srange`, whose `D` velocity
reaches +3) decouple upward from the frozen pair `{c₃, H}`.** Hence any conjecture of the form
`[tr/srange-side] ≤ [function of c₃, H]` is inflation-fragile, and its witness is the 3-cycle atom
padded with transitive singletons — THM-1830, the tournament analog of WOWII-103's triangle+leaves.

## 4. Applications (all verified exhaustively over the 530 iso classes n ≤ 7)

| conjecture | source | verdict via THM-1855 |
|---|---|---|
| `srange ≤ tr` | kps THM-1845-adj. | **fragile** (D−: srange outruns tr); breaks n=7 (c₃=4, tr=5, srange=6) |
| `srange ≤ tr + c₃` | **new (repair)** | holds n≤7, join-monotone → strong core |
| `srange ≤ 2(tr−1)` | **new (repair)** | holds n≤7, join-monotone → strong core |
| `c₃ ≤ H` | kps (open) | **join-monotone (proved) → reduces to strong core**; verified on 396 strong classes |
| `n − c₃ ≤ tr` | kps THM-1845 | additive both sides → strong core (all n) |
| `tr ≤ smax + 1` | kps THM-1845 | join-monotone → strong core |
| `dom + tr ≤ n + 1` | klein THM-1850 | join-monotone (dom left-projects) → strong core; D−-tight-preserving |

**The `c₃ ≤ H` reduction (airtight join-monotonicity).** If `H₁=1` then `T₁` is transitive
(Rédei/Moon: `H=1 ⟺ transitive`), so `c₃₁=0` and `c₃ = c₃₂ ≤ H₂ = H`. Symmetrically for `H₂=1`.
Else `H₁,H₂ ≥ 2 ⟹ (H₁−1)(H₂−1) ≥ 1 ⟹ H₁H₂ ≥ H₁+H₂ ≥ c₃₁+c₃₂ = c₃`. ∎
So kind-pasteur's open `c₃ ≤ H` is now exactly a statement about **strongly connected** tournaments.

## 5. WOWII-103 proper (king-eccentricity)

With `tr` (the `α`-analog) and average king-eccentricity `ecc_avg` (finite exactly on strong
tournaments), the naïve `tr ≤ ⌊(n−1) − ln ecc_avg⌋` is too loose (breaks at n=4). The faithful
object is the **transcendental threshold**: among strong tournaments n≤7 the achievable `ecc_avg`
straddles `e` at `19/7 = 2.71429` (just below) and `17/6 = 2.83333` (just above) — the tournament
analog of WOWII-103's `30/11 > e`. Calibrating a *tight* directed-103 whose only slack is the
`ln`-term is the follow-up (HYP-8641).
