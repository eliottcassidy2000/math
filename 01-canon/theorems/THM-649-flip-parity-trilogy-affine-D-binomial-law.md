---
id: THM-649
title: The flip-parity trilogy — (A) D(t) = c3(flip t) − c3(t) is AFFINE in the tile bits for every tiling, with coefficient +1 exactly on the two staircase legs, +2 on the apex, 0 on all interior tiles; (B) hence on grid-symmetric tilings D is a shifted 2·Binomial(n−2) (the σ-paired legs + fixed apex are the n−2 independent carriers, all interior orbits inert), centered by the flip involution; (C) corollaries: at odd n the support misses 0, so no grid-symmetric tiling has T(flip t) ≅ T(t) — blue self-loops do not exist at odd n; at even n the D=0 pool is the central binomial column and the self-loop count is 1, 2, 4 at n = 4, 6, 8 (2^{n/2−2}, verified), with ×2-mod-(n+1) among the witnesses (a permutation iff n is even)
status: (A), (B) and the odd-n corollary of (C) PROVED (derivation below; machine-verified EXACT over every tiling at n=4..7 for (A), over all gridsym tilings at n=4..9 for (B), 16/16 and 512/512 direct c3-parity checks at n=5,7). Even-n self-loop count: the conjectured 2^{n/2−2} law is **REFUTED at n=10** (klein-S173 full census: 48 loop tilings = **24 lines**, not 8; affine-D D=0 prefilter 286720 of 2^20 gridsym, exact backtracking-iso verification of every loop AND exact-iso completeness over all 14350 invariant survivors). The true sequence is 1, 2, 4, 24 at n = 4, 6, 8, 10. The multiplicative twist family ⟨2,−1⟩ mod n+1 accounts for all loops at n=4, HALF at n=6, NONE at n=8; the count jump has NO primitive-root explanation (MISTAKE-124: 2 is primitive mod 5, 9, 11 alike, non-primitive only mod 7) and matches neither 2^{n/2−2} (fails n=10) nor (n/2−1)! (fails n=8); the witness classification needs a non-affine mechanism (HYP-4961's 2-adic hierarchy insufficient as stated).
source: klein-2026-07-07-S168/S169/S170/S171 (HYP-4921, HYP-4931, HYP-4941, HYP-4951)
depends_on:
  - THM-643   # blue/black parity structure (T1-T5 family)
related:
  - THM-644   # gridsym-fiber law / anti-Rédei family
  - THM-647   # anti-Rédei proved (monad)
  - LEM-003   # orbit-stabilizer freeness
external: none (elementary).
---

# THM-649 — the flip-parity trilogy

## (A) The affine-D law (proved)

For any tiling `t` of the staircase (base path `n → … → 1`, tiles `(x,y)` with `x−y ≥ 2`),
let `D(t) := c₃(T(flip t)) − c₃(T(t))`, where `flip` complements every tile bit. Then

> `D(t) = D(0) + Σ_tiles c_{(x,y)}·bit_{(x,y)}` with
> `c_{(x,y)} = 2(s^b_x − s^b_y) + (tiledeg_x − tiledeg_y)`, which equals
> **+1 on the leg tiles** `(x,1)` (x = 3..n−1) **and** `(n,y)` (y = 2..n−2),
> **+2 on the apex** `(n,1)`, and **0 on every interior tile**.

*Proof.* `c₃ = C(n,3) − Σ_v C(s_v,2)` gives `D = −Σ_v [δ_v s_v + C(δ_v,2)]` with
`δ_v = tiledeg_v − 2s^t_v`. Expanding with `s_v = s^b_v + s^t_v`, the quadratic terms
`∓2(s^t_v)²` cancel identically, leaving `D = −K(n) + Σ_v w_v s^t_v`,
`w_v = 2s^b_v + tiledeg_v − 1` — linear in the tile-scores, hence affine in the bits with
per-tile coefficient `w_x − w_y`. Since `s^b_v = 1` for `v ≥ 2` (0 at v=1) and
`tiledeg_v = n−2` at `v ∈ {1,n}`, `n−3` otherwise, the coefficient vanishes unless an
endpoint of the tile is a path-endpoint, giving the stated values. ∎
(Machine-verified exactly over every tiling, n = 4..7.)

## (B) The binomial law on grid-symmetric tilings (proved)

The grid reflection `σ` maps `(n,y) ↔ (n+1−y, 1)`, pairing the two legs into `n−3`
orbits and fixing the apex. On gridsym tilings each leg-orbit contributes `0 or 2` to D,
the apex contributes `0 or 2`, and every interior orbit is inert:

> over the `2^{(m+f)/2}` gridsym tilings, `D` is distributed EXACTLY as a shifted
> `2·Binomial(n−2)` (i.e. `#{D = −(n−2)+2k} = C(n−2,k)·2^{(m+f)/2−(n−2)}`),
> centered at 0 because the flip involution negates D and is fixed-point-free.

(Machine-verified exactly at n = 4..9, through 65,536 gridsym tilings at n = 9.)

## (C) Corollaries and the even-n count

- **Odd n (proved):** the support `{−(n−2), …, n−2}` (step 2) misses 0, so
  `c₃(flip t) ≠ c₃(t)` for every gridsym `t`: `T(flip t) ≇ T(t)` — **blue self-loops do
  not exist at odd n** (the S168 theorem, now the no-zero-column corollary; equivalently:
  the diagonal `Σ_reps tiledeg = C(n−1,2) − (n−3)/2 = 2k²−2k+1` is odd at `n = 2k+1`).
- **Even n (counts census-exact; law open):** the `D = 0` pool is the central column
  (`C(n−2,(n−2)/2)·2^{…}`; n=8: 1280, census-exact; n=10: 286720), and the self-loop
  line count is **1, 2, 4, 24 at n = 4, 6, 8, 10** — the conjectured `2^{n/2−2}` is
  REFUTED at n=10 (klein-S173, doubly-verified exact census). The witness
  isomorphisms include the **doubling map `v ↦ 2v mod (n+1)`** — a permutation of the
  runners iff 2 is invertible mod n+1, iff **n is even** (the structural face of the
  parity dichotomy; the joint constraint system [gridsym = ×(−1)] + [twist = ×2] lives on
  the orbits of `⟨2,−1⟩ ≤ (Z/(n+1))^×` acting on vertex pairs, each twist an explicit
  GF(2) affine system). HONEST STATUS: the multiplicative family accounts for all loops
  at n=4 (×2 and ×2⁻¹, each a unique solution), half at n=6, none at n=8 — the remaining
  witnesses are non-affine; the full witness classification (and hence the Burnside proof
  of `2^{n/2−2}`) is OPEN.

## Files
`lrc_pairing_signflip_klein_S168.out`, `lrc_c3diff_law_selfloops_klein_S169.out`,
`lrc_binomial_proof_burnside_klein_S170.out`, `lrc_pairspace_dimension_klein_S171.out`
(all in 05-knowledge/results; scripts inline). The D-carriers being the triangle's
legs + apex ties the flip-parity structure to the project's geometric foundation
(CLAUDE.md, "everything is the triangle").
