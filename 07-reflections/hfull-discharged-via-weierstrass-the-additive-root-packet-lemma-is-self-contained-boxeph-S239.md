# `hfull` discharged via Weierstrass: the additive root-packet lemma is self-contained

*boxeph-2026-07-22-S239. Owner: keep working on `hfull` (the full-root Lagrange sum), think
Weierstrass. Builds on my S238 (`root_packet_eq_zero`, the additive-route core taking `hfull` as a
hypothesis), codex's `GMC2FullRootLagrange.sum_pow_div_derivative_nodal_eq_zero` (Lagrange for the
monic nodal). New Lean `GMC2FullRootPhi.lean`, kernel-pure.*

## The gap, and the Weierstrass idea

`root_packet_eq_zero` (S238) needed the full-root sum `∑_{α∈roots} αᵏ/Φ'(α) = 0` (`hfull`) as a
hypothesis. codex proved this for the **monic nodal** `∏(X−root)`, but the actual `Φ = X^M − tR` is
**non-monic** (leading coeff `−t·r_d`). The Weierstrass factorization bridges them:

> `Φ.map = C(leadingCoeff) · ∏(X − root)` (the polynomial *is* its leading coefficient times the
> Weierstrass/nodal product of its roots), so differentiating and evaluating at a root gives
> `Φ'(α) = leadingCoeff · nodal'(α)`, and the `Φ`-sum is `(1/leadingCoeff)` times the nodal sum `= 0`.

## Delivered: `hfull` discharged (kernel-pure)

`GMC2FullRootPhi` (`#print axioms = [propext, Classical.choice, Quot.sound]`):
- **`phi_map_eq`** — the Weierstrass product form: `Φ.map = C(lc) · Lagrange.nodal (aroots).toFinset id`
  (from `Splits.eq_prod_roots`, with `Multiset.dedup_eq_self` for the distinct roots).
- **`aeval_deriv_eq`** — the derivative-at-root proportionality `aeval α Φ' = lc · nodal'(α)`
  (`aeval_def`/`eval_map`/`derivative_map` to move to `Φ.map`, then `derivative_C_mul`).
- **`full_root_sum_eq_zero`** — the full-root Lagrange sum for the **non-monic** `Φ` is `0`
  (`k+1 < deg Φ`): reindex the root-subtype sum to `(aroots).toFinset` (`Finset.sum_subtype` + a
  `mem_aroots ↔ mem_rootSet` bridge), rewrite `Φ'(α)` by `aeval_deriv_eq`, factor out `1/lc`
  (`Finset.sum_div`), and close with codex's nodal Lagrange sum + `card_rootSet_eq_natDegree`.
- **`root_packet_eq_zero_selfcontained`** — the **self-contained** additive root-packet lemma:
  `b_k(S) ∈ F ⟹ b_k(S) = 0`, with `hfull` now discharged internally (`root_packet_eq_zero` ∘
  `full_root_sum_eq_zero`).

## State of the additive route

The additive-route **algebraic core is now fully self-contained and kernel-pure**: for `Φ` irreducible
over a char-0 field, splitting with distinct roots, the barycentric packet sum in `F` is `0` — with
**no `hfull` hypothesis, no THM-1550, no small-root product, no Hensel factorization, no Vieta**. The
whole additive skeleton — `GMC2GalRootAction` (transitive action) + codex's `GMC2LaurentShiftCheckA`
(additive orbit machinery + nodal Lagrange) + `GMC2RootPacketConcrete` (the concrete root-packet
lemma) + `GMC2FullRootPhi` (Weierstrass `hfull` discharge) — is kernel-pure.

The **only remaining input** is the `b = 1` wrapper: `∑_{α∈S_+} α^{M−1}/Φ'(α) = F(t) = ∑ D_m tᵐ` for
the positive-valuation packet, giving `b = 1 ∈ F` under `D_m = 0 ∀m`; the self-contained root-packet
lemma then gives `b = 0`, a contradiction, proving DvdK. This is the one shared valuation/Newton-polygon
input (kind-pasteur's verdict), but in additive/sum form — the partial-fraction residue identity
(THM-2101 §9's t-adic proof, or §8's transcendental one; cf. mac-mini-S163's root-free
`F(t) = [x⁰] xᴹ/(xᴹ−tR)`).

## Scope

Honest: `hfull` is **discharged** kernel-pure via the Weierstrass product form — the additive
root-packet lemma is now self-contained, removing the last hypothesis-level gap in the additive
algebraic core. Not full GMC(2): the `b = 1` residue-sum wrapper remains (the shared valuation input,
additive form). One checkpoint pushed.

Links: HYP-8985, HYP-8980 (S238 core), THM-2101 (codex), codex `GMC2LaurentShiftCheckA`, mac-mini-S163,
[[the-concrete-root-packet-lemma-the-additive-route-core-boxeph-S238]].
