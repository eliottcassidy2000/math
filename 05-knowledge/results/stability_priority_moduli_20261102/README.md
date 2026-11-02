# Shared matching budgets, exact arithmetic precision, and guarded path transfers

2026-11-02. Current scope: **PROVED within this packet**, with exact
reproduction; no separate agent/human audit or literature-priority claim.

* **Weighted tournament stability is sharp for every n>=276.** The formula
  is K_n=(4/R)[1-(3/4)^R-1_(n even)(1/4)^R], R=floor(n/2), with the same
  complete fair-pair equality classification. Combined with the inherited
  n=3..9 proofs, the remaining interval is **10..275 OPEN**.
  [Proof](STABILITY_N276.md) · [new sharp graph lemma and intermediate proof](STABILITY_N500.md).
  The stronger discarded-degree inequality is W>=2D^2/n. The final proof
  explicitly incorporates incoming Session139's per-edge geometric-tail
  mechanism and keeps the binary baseline in both matching estimates.
* **Every finite FRACTRAN observation has an exact minimal residue modulus.**
  [P1-P6](FRACTRAN_PRECISION.md) retains priority and reachable-word guards,
  gives every-prime precision formulas, and classifies the reciprocal family
  (a/b,b/a) for every modulus and clock. With terminal residue unobserved,
  its *smallest arbitrary clocked observer* has h+2 states, despite needing
  modulus a*b^h as a residue chart. These are different minimization problems.
* **THM-4852 is corrected.** Prime support alone does not determine a rule
  from one primorial-base digit. States 2 and 8 in (1/4,3/2) are the small
  hostile. PRIMEGAME's actual squarefree case and phase law remain valid.
* **Joint alignment has a literal local cost.** [The pair reserve](PAIR_RESERVE.md)
  bounds cyclic probability by both endpoint variances and squared disagreement
  at a common third vertex. Equal variance summaries can lose this information.
* **A source-preserving H map can still lose the target.** Ordinal fair-pair
  padding sends the normalized stability defect to H-1. It therefore cannot
  lift a hypothetical small-order counterexample merely by preserving H.
  [Synthesis, recovered sources, and next tests](SYNTHESIS.md) connects this
  to the cited quasi-Riemann replication mechanism, incoming five-vertex
  zeta blindness, actual 5/11 context projectors, and square-sum/FRACTRAN guards.

[Audit and exact universes](AUDIT.md) · [normal replay](normal.json) ·
[optimized replay](optimized.json) · [source pins](SOURCE_PINS.json)

Reproduce from the repository root:

```powershell
python -X utf8 -B 05-knowledge/results/stability_priority_moduli_20261102/verify.py --out 05-knowledge/results/stability_priority_moduli_20261102/normal.json
python -O -X utf8 -B 05-knowledge/results/stability_priority_moduli_20261102/verify.py --out 05-knowledge/results/stability_priority_moduli_20261102/optimized.json
```

No extension of A000568, unrestricted Hamilton-path spectrum, all-size
strong-atom complement closure, LRC(14), Collatz, or RH is claimed. Square-sum
paths at orders 15,16,17,23,25 are independently counted here as finite controls;
the all-size theorem is credited, not rediscovered. The twelve-vertex 81/27
pair retains its inherited joint-response obstruction and is not classified
by these scalar or modular observations.
