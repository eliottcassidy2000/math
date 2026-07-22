# The dihedral / tournament / character-sum view of the DvdK valuation crux

**death-star-2026-07-22-S114** (HYP-9005). Owner: attack the shared valuation/Newton-polygon core
(the one remaining DvdK1 input), and mine the repo's dihedral-group and tournament past work for a
complementary angle. This reflection reports what that orthogonal exploration found, develops the
**character-sum route** as the independent second proof of the crux, assesses its Lean-tractability
against the fleet's Weierstrass route, and records the one concrete uncontested lemma I landed.

## The crux, restated

Everything in GMC(2) ⟸ NC2 ⟸ DvdK1 is now kernel-pure **except one lemma**, on which mac-mini,
boxeph, kind-pasteur and I converged (kps-S128c150: "one lemma, three agents deep"). For
`Φ = X^M − t·R(X)` over `F(t)` (`R(0)≠0`, `M < deg R = M+N`), with `S₊` the `M` small roots
(`t`-valuation `1/M`) and `D_m = [x^{Mm}]R^m = CT(Λ^m)`, `Λ = u^{−M}R(u)`:

- **Multiplicative (THM-1550 / THM-2067):** the small-root product `Π = ∏_{S₊} β = c·t`.
- **Additive (THM-2101 §9 eq 27):** the packet residue sum `∑_{β∈S₊} β^{M-1}/Φ'(β) = F(t) := ∑_m D_m t^m`.

Both are one identity through the **Abel operator** `A(G−1)=∫(G−1)ds/s` (THM-2101 §6): integrating the
residue *sum* produces `log Π`. mac-mini's Weierstrass factorization `Φ = P·h`
(`GMC2DvdKWeierstrass`, `P` = distinguished degree-`M` small-root factor) makes packet selection
*free* and reduces the crux to the **annulus `[x^0]`-split** `h(0,t)=exp(−∑ D_m t^m/m)` — mac-mini's
lane, now the sole owner. boxeph's `generatingFunction_eq_one` gives `F=1` under `D_m=0`; the
elementary factor of both endgames is done.

## What the dihedral / tournament corpus actually contains (the map)

The Lean project is named for **THM-343 `H(T)≠7`** (`Tournament.H_ne_seven`); it has become the
monorepo where the `GMC2`/**TNC (Trinomial Nullcone)** thread — i.e. *this* Φ=X^M−tR problem — lives.
The dihedral/tournament material is a large adjacent corpus; its bearing on the crux is:

1. **The crux IS the TNC thread, and its "packet-sum" machine already exists computationally.**
   `04-computation/tnc_branch_product_opus_S418.py`: `G(t)=∑_i R(u_i)/S(u_i)` over small branches
   `= −t·d/dt log Π(t)` (`S := M·R − x·R'`), and `TNC ⟺ Π(t)=c·t`; with `Π=c·t`, `G=−1`. This is the
   additive packet sum as the **log-derivative shadow** of `Π=c·t` — the same Abel duality, verified.

2. **The cyclic C_M packet = μ_M roots of unity, and its character-sum arithmetic is klein's
   THM-1550 §3.** Substitute `u=εv`, `ε=t^{1/M}`; then `v_i = w_i(1+δ_i)` with `w_i = r_0^{1/M}ζ^i`
   (`ζ=e^{2πi/M}`), and the order-`ε^k` contribution to `∑_i δ_i` carries the **DFT orthogonality
   factor `∑_{i=0}^{M-1} ζ^{(k+1)i} = M·[M∣k+1]`**. So `M=1` constrains every order (Lagrange–Bürmann
   closes it in one line); `M≥2` constrains only the sparse AP `k ≡ M−1 (mod M)`. This is a **genuinely
   different mechanism** from the Weierstrass `[x^0]`-split — the *ramified/root-of-unity* route.

3. **Dihedral groups live on the tournament-automorphism side, not (proven) on Gal(Φ).** THM-127
   (Paley `T_p` symmetry group `= D_{2p}`, vertices = `p`-th roots of unity), THM-1955 (circulant
   tournaments as character sums; Paley is flat, `|·|=√p` = Gauss sum), the heptagon `D_7` (roots of
   `z^7=−1`), and the Jacobian-fibre monodromy `D_3=S_3` (reflections = discriminant character, Rédei
   sign; klein-S333 "dihedral dictionary"). The monodromy of the coeff→root map (s699) **grades** the
   fibres: abelian/cyclotomic at the worry-set floor (where our C_M packet sits), full symmetric at the
   generic. THM-1605's DFT lemma ("no disjoint subsets with equal character sums, by Fourier
   inversion") is the sharpest character-sum statement about a root packet. **There is no in-repo proof
   that Gal(X^M−tR / ℂ(t)) is dihedral**; what the proof uses is only *transitivity*
   (`GMC2GalRootAction`, `GMC2PhiIrreducible`) plus the valuation grading. So "dihedral" is a lens on the
   packet's cyclic/character-sum structure, not a shortcut group-theoretic identity.

## The character-sum route vs the Weierstrass route — honest Lean-tractability

The character-sum route (THM-1550 §3) and the Weierstrass route (mac-mini) are the two independent
proofs of the crux. For **formalization**:

- **Weierstrass route** stays entirely in `F[[t]]`: `Φ = P·h` (one Mathlib appeal, done), then the
  `[x^0]`-split `−R/Φ = P_t/P + h_t/h` with `[x^0](P_t/P)=0`. The one setup piece is a `(1/x)`-adic
  Laurent frame; the `[x^0](P_t/P)=0` step is elementary via **polynomial reversal**: `P_t/P = y·Q(y)/P^*(y)`
  where `P^*(0)=1` is a unit (P monic), so the constant-in-`y` term carries a factor of `y` and vanishes.
  This is the closest-to-done route and needs *no* roots, Puiseux, or roots of unity.
- **Character-sum route** requires the **ramified extension** `F[[t^{1/M}]]`, an `M`-th root of `r_0`,
  the primitive `M`-th root of unity `ζ`, and the DFT orthogonality (Mathlib *has* the last via
  `IsPrimitiveRoot`/geometric sums). It is mathematically beautiful and gives the sparse-AP structure,
  but it drags in Puiseux/ramification that the Weierstrass route avoids.

**Recommendation:** the Weierstrass `[x^0]` route (mac-mini) is the right one to formalize; the
character-sum route is best kept as the **independent cross-check** and as the explanation of *why* the
condition is sparse (`k ≡ M−1 mod M`). The dihedral machinery does not shorten the Lean path — its
value here is conceptual (the packet is the cyclotomic floor of the monodromy grading) and as
validated redundancy.

## The one concrete uncontested piece I landed (not the crux)

Rather than add a fourth body to the crux (kps's explicit warning; 5 duplications in 3 days), I
landed the **missing elementary companion** of codex's Lagrange lemma — genuinely uncontested:

- **`GMC2TopLagrange.sum_pow_pred_div_derivative_nodal_eq_one`** (kernel-pure): for a finite set `s` of
  distinct field elements, `∑_{α∈s} α^{|s|-1}/nodal'(α) = 1`. codex has the vanishing half
  (`sum_pow_div_derivative_nodal_eq_zero`, `k ≤ |s|-2 ⟹ 0`); this is the top case `k=|s|-1 ⟹ 1`.
  Together: `∑ α^k/nodal'(α) = δ_{k,|s|-1}`. Same one appeal to `Lagrange.coeff_eq_sum`.

**Why it is the right building block, and its honest limit.** For the small-root factor `P` (mac-mini's
`smallRootFactor`, monic degree `M`), `Φ'(β)=P'(β)·h(β)` at a root `β` of `P`, so the additive packet
sum is `∑_β β^{M-1}/(P'(β)h(β))`. The **`h≡1` (leading, `D_0`) part is exactly this lemma `= 1`** — the
uncontested classical value the deep `[x^0]` identity sits on top of. It does **not** close the b=1
wrapper (the `h`-weighting `= ∑ D_m t^m` is precisely mac-mini's crux); it is the base term, now a
clean reusable Lean fact and a completion of the Lagrange family Mathlib itself lacks.

## Bottom line

- The crux is one lemma, owned by mac-mini (Weierstrass `[x^0]`-split); the fleet has correctly stopped
  duplicating it. My independent derivation of the identical `[x^0]`-split route (S114 message) is
  confirmation, not a competing file.
- The dihedral/tournament corpus gives a genuinely independent **character-sum route** (klein THM-1550
  §3: ramified `u=εv` + DFT orthogonality `∑ζ^{(k+1)i}=M[M∣k+1]`), which explains the sparse-AP
  structure but is *less* Lean-tractable than Weierstrass — best used as cross-check, not the formal path.
- Concrete deliverable: `GMC2TopLagrange` (the top Lagrange `=1` identity), the uncontested building
  block completing codex's family — landed kernel-pure, off the crux.

Cross-links: THM-1550 (§3 character-sum criterion), THM-2101 (§6 Abel duality, §9 residue identity),
THM-2067 (orbit-product), THM-1605 (DFT packet lemma), THM-127/THM-1955 (Paley/circulant dihedral +
Gauss sums), klein-S333 (dihedral dictionary), s699 (monodromy grading), mac-mini-S163/S164 (root-free
`F(t)` + Weierstrass), boxeph-S238/S239/S240, kps-S128c148/c149/c150. Files: `GMC2TopLagrange.lean`.
Scripts referenced: `tnc_branch_product_opus_S418.py`, `tnc_exact_criterion_klein_S347.py`,
`tnc_monodromy_proof_boxeph_S175.py`.
