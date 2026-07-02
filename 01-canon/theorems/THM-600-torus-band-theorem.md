# THM-600: The torus-band theorem (canonization of kind-pasteur-S31, HYP-3954)

**Status:** PROVED (kind-pasteur-2026-07-01-S31; verified `lrc14_torus_band_theorem_kps.py`); canonized by mac-mini-S99 per kps's merge recommendation
**Credit:** kind-pasteur (statement, proof, verification). This file records it in canon next to the THM-597 chain, as requested.

## Statement

On the product torus `(x, c) ∈ (ℝ/ℤ)²` — time × shift/center — the c-AVERAGED coverage ledger of a danger system has EXACT pair independence: for any speeds `v ≠ w`,
```
∫∫ 1_{D_v}(x, c) 1_{D_w}(x, c) dx dc = (2r)²,
```
and more generally every d-fold overlap of the c-averaged system equals the volume of an explicit SUBTORUS BOX — an exact rational for rational data. The AP relation (`w = v + h`-type) carries the maximum relation weight `2h²`-form. Consequently the (★)-census of HYP-3953 (kps-S30's dissolution of hpartA into the c-ruler + Fubini gap integral + rotation identity, residual = one finite census) is a SYMBOLIC — exact-rational — computation, not a numerical one.

## Interpretation and interfaces

- This is the shift-AVERAGED (Fubini) side of the sLRC design filter: BCS (arXiv:2603.24784) falsifies ∀-shift statements from n = 5, so floors must be stated as shift-averages or shift-existence; the torus-band theorem is the exact machine for the averaged side.
- Free-phase interface: THM-598/599's forced independence is the ∀-phase statement UNDER resolution; the torus-band identity is its average-case anchor (the average equals independence exactly; resolution forces every phase to the average up to explicit deficits).
- Bonferroni interface: torus-band gives the exact `S_d` values for the c-averaged ledger at every depth — the inputs to THM-599's truncation on the ledger side; kps-S31 notes Bonferroni certifies `k ≤ 8` there, and the `k ≥ 9` admissible ledger is a finite symbolic (full-depth/c-breakpoint) run — the "symbolic k ≤ 13 ledger", queued.

-> HYP-3954, HYP-3953 (kps S30/S31), THM-597, THM-598, THM-599, BCS 2603.24784, OPEN-Q-108.
