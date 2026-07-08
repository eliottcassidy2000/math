---
source: opus-2026-07-08-S153
status: the E_j matched-tuple energies identified via the full-set relation lattice (E_2 = R2
  exact, E_j = support-2j relations); the resummed relation-lattice mean formula E[W] =
  sum_{n in L} prod psihat(n_i) verified. HONEST: the variance E_j-expansion is NON-PERTURBATIVE
  (mac-mini LEM-007: support terms grow, truncation overshoots 12x), so exact E_j do NOT give a
  uniform Var<=c*R2 bound -- reconciles/subsumes into the fleet's discrepancy route (LEM-005) +
  per-shape reframe (klein-S183). Concurrent with kps LEM-008 (per-subset overlap-mass kernels).
tags:
  - lrc14
  - covering-floor
  - variance
  - relation-lattice
  - additive-energy
  - non-perturbative
---

# The E_j energies are the relation lattice; E_2 = R2; but the expansion is non-perturbative

**opus-2026-07-08-S153 (HYP-5427).** Owner: derive the `E_j` triple/quad matched-tuple
energies exact. They are exactly the additive-relation-lattice counts — and `E_2 = R2` cleanly
— but the fleet's concurrent work (mac-mini LEM-007, klein-S183) shows the `E_j`-expansion of
the variance is **non-perturbative**, so this derivation is an exact structural identity, not a
route to the uniform bound. This note states both honestly.

## The full-set relation-lattice formula (resummed kps LEM-008)

kps-S83 LEM-008 gives the per-subset overlap mass `E[L_S] = sum_{m in Lambda_S} prod c_{m_a}`,
`Lambda_S = {m in Z^S : sum m = 0, sum m_a e_a = 0}` the balanced additive-relation lattice of
`S`. Summing the inclusion-exclusion `W = sum_S (-1)^|S| L_S` **resums** into a single full-set
lattice sum (the `psihat(0) = 1-theta` spectator weights absorb the `(-1)^|S|`):

> **`E[W] = sum_{n in L} prod_{i=1}^k psihat(n_i)`,  `L = { n in Z^k : sum n_i = 0, sum n_i e_i = 0 }`,**

`psihat(0) = 1-theta`, `psihat(m) = -c_m` (`c_m` the arc Fourier coeff). Verified to `1e-5` on
`{0,1,2}, {0,1,3}, {0,1,2,3}`. `L` is the FULL-set relation lattice (kps's `Lambda_S` are its
support-`S` slices); for a Sidon set `L = {0}` and `E[W] = (1-theta)^k = (6/7)^k` (the iid
value). The same resummation gives `Var(W) = sum_{n in L} K(n)` with `K(0)` the Poisson diagonal
and `K(n)` (`n != 0`) the resonance of relation `n`.

## E_2 = R2, and the E_j as relation supports

A relation `n in L` of the form `(1,-1,-1,1)` on `{i,j,k,l}` means `e_i + e_l = e_j + e_k` — an
**additive quadruple**. The pair (support-2-vector) resonance of `Var(W)` is carried exactly by
these support-4 relations, and their count is the additive energy: the full difference energy
`sum_{all d} r_d^2` splits as `k^2` (the `d=0` Poisson diagonal) `+ R2` (the reduced energy,
`d != 0`), so

> **`E_2 = R2`** (the pair matched-tuple energy = the reduced additive energy = the support-4
> additive-quadruple relations, once the `d=0` diagonal is placed in the Poisson term).

Higher up: `E_3` is carried by the support-6 relations (additive sextuples, `e_a+e_b+e_c =
e_d+e_e+e_f` and the like), `E_4` by support-8, etc. — each an exact relation-lattice count, the
`|S|>=3` siblings of `R2` (kps LEM-008's triangle/quad lattices are the per-subset generators).
So the `E_j` ARE exact: `E_2 = R2`, `E_j = ` the support-`2j` (and mixed) relation energies.

## Why exact E_j do not close the bound (the non-perturbative wall)

The catch, established rigorously by mac-mini LEM-007 (and I reconcile my S152 assembly to it):
the support-decomposition of `Var(W)` does **not truncate**. The support-`r` variance
contributions GROW (`|W_2|^2, |W_3|^2, |W_4|^2 = 0.077, 0.226, 0.932` at `k=11`), and the true
small `Var = 0.047` emerges from massive cancellation *between* support levels — the cross-terms
`W_r conj(W_s)` are large and negative. So:

- My S152 `Var ~ (1-theta)^{2(k-2)}(R2/2)c_2 * 1.28` was a heuristic that matched numerically but
  is **not** the rigorous support-2 term (`|W_2|^2 = 0.077 != 0.037`); the honest picture is
  mac-mini's non-truncating expansion.
- Therefore **exact `E_j` are exact building blocks but not a convergent series** — knowing
  `E_2 = R2`, `E_3`, `E_4`, ... exactly does NOT yield `Var(W) <= c*R2` uniformly, because the
  series is non-perturbative (kps-S82's "the resonance sign is non-perturbative"). This is why
  the uniform-`Var<=c*R2` target is "doubly dead" (klein-S183: mac-mini truncation fails + the
  per-family `c` is non-uniform).

## Where this leaves the k=11 leg (fleet consensus)

- The **overlap kernels are fully exact**: masses `E[L_S]` (kps LEM-008), variance kernels `c_j`
  (opus-S152), the relation-lattice `E[W]`/`Var(W)` (this note), `E_2 = R2`.
- The **uniform moment/`E_j` route is dead** (mac-mini non-perturbative; klein non-uniform `c`).
- The honest routes forward are the fleet's: **discrepancy** (LEM-005, the equidistribution of the
  phases controlling `Var` directly) and the **per-shape** reframe (klein-S183: `Var <= 2.02 E[W]^2`
  localized to high `R2`, binding = block+outlier, needing extended exhaustive or
  outlier-decorrelation). The exact kernels feed those, but do not replace them.

## Ledger

- VERIFIED/DERIVED: the full-set relation-lattice `E[W] = sum_{n in L} prod psihat(n_i)` (resummed
  kps LEM-008); `E_2 = R2` (pair energy = reduced additive energy = support-4 additive quadruples);
  `E_j` = the support-`2j` relation energies (exact structural identities).
- RECONCILED (honest): the `E_j`-expansion of `Var(W)` is non-perturbative (mac-mini LEM-007);
  exact `E_j` are building blocks, NOT a uniform bound; my S152 `1.28` heuristic corrected to the
  non-truncating picture.
- ROUTE: discrepancy (LEM-005) + per-shape (klein-S183), fed by the exact kernels.
- Files: (verification inline in session); builds on kps-S83 LEM-008, mac-mini-S57 LEM-007,
  klein-S183 (LEM-008 collision -> LEM-008 covering-variance), opus-S152 (c_j) / S151 (c_pair).
