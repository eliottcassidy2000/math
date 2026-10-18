---
id: THM-4695
title: "First-order universality of the quantum transfer: for every Virasoro minimal model Vir_(p,q) with N = (p-1)(q-1)/2 >= 3, the first-order transfer T_1 from Koszul 1-cycles of the derived arc space of Spec C[x]/(x^N) to the classical defect ker pi = grI/J equals (c - c_(2,2N+1)) T*_N modulo J, with T*_N depending only on N. The length-(N-1) part of the singular vector is alpha L_-4 L_-2^(N-2) + beta L_-3^2 L_-2^(N-3) with alpha = -2N(N-1)/3 - Nc/12, beta = 11N(N-1)/48 + 5Nc/96. Hence all non-free minimal models with the same N have the same first-order image, which contains THM-4671's relation b_N, and T_1 == 0 mod J exactly when p = 2"
status: >
  PROVED, given one CITED input: classical freeness of Vir_(2,2N+1).
  - Source: van Ekeren-Heluani, arXiv:1804.00017 = CMP 386 (2021), Thm 15.13. The PBW-filtration form goes back to
    Feigin-Frenkel 1993, section 3.
  - The term "classically free" is defined in Andrews-van Ekeren-Heluani (IMRN 2023), section 1.1, and in van
    Ekeren-Heluani arXiv:2103.06322, section 1.6.
  - Audited (s44_audit_20261018/REPORT.md): proofs (a)-(d) re-derived line by line; alpha, beta confirmed exactly over
    Q on 16 models; universality, lift-independence and the vanishing at c_(2,2N+1) reproduced with an independent
    engine to weight 24-30.
  - Parts (a) and (b) are elementary.
  - Part (c) is: (b) + an affine-dependence argument + the cited input at c_(2,2N+1).
  - THM-4671 (kind-pasteur), the closed form A = N(2N+1)(c - c_(2,2N+1))/8 on the first class, is the special case
    "T_1 on the Euler syzygy".
  FINITE-EXACT controls (mod 2^31 - 1, quantum_transfer_s44/):
  - (b) against exact singular vectors for 11 models (exact over Q; first_order_universality.json). The N = 4
    collinearity test uses only two points, so it is vacuous. The audit checked 16 models exactly, including
    N = 3, 7, 12;
  - (c) T_1 images equal as subspaces mod J, and T_1 proportional with ratio (c - c_0)/(c' - c_0) on a common cycle
    basis: (3,7) vs (4,5) through weight 22, (3,10) vs (4,7) through weight 25;
  - first-order dimension totals equal within each N: N = 6, 143 through weight 28; N = 9, 45 through weight 30.
  Builds on kind-pasteur-2026-10-18-S1 (singular_support_arc_spaces_classical_freeness_20261018.md, section 1.3 and
  open problem 6; engine vir_ss.py reused unchanged).
  Novelty (checked by the S44 audit): the general-c statement is not in the literature as stated. It has close
  precedents:
  - Andrews-van Ekeren-Heluani, section 2.4: their Ising lift is a first-order transfer of the N = 3 Euler syzygy,
    plus an L_-1^3 v term with symbol in J.
  - Feigin, Funkts. Anal. Prilozh. 45:4 (2011): the Euler relation is the associated graded of the level-3 singular
    vector of M(c_(2,2p+1), 2p).
  The factor is a Kac factor:
  - c - c_(2,2N+1) = (3/(2N+1)) Kac_3(c, 2N), with Kac_3(c,h) = (1/3)(3h^2 - (7-c)h + 2 + c) the level-3 Kac factor
    of M(c,h).
  - So THM-4671's coefficient is A = (3N/8) Kac_3(c, 2N).
  The spectral-sequence language is standard for filtered maps.
session: opus-2026-10-18-S44
source: 05-knowledge/results/li_spectral_sequence_quantum_transfer_20261018.md
scripts:
  - 04-computation/quantum_transfer_s44/vir_quantum_transfer_20261018.py (T_1, T_2 per bidegree; boundaries -> J; T_1 in grI)
  - 04-computation/quantum_transfer_s44/first_order_universality_20261018.py (+ first_order_universality.json)
  - 04-computation/quantum_transfer_s44/vir_transfer_orders_20261018.py (all orders; orders_*.out)
related: [THM-4671, HYP-9294, HYP-9315, HYP-9293, THM-4710, van Ekeren-Heluani 2021 (Thm 15.13), Andrews-van Ekeren-Heluani 2020, Feigin 2011, Feigin-Frenkel 1993, Arakawa-Linshaw 2019]
---

# THM-4695 — first-order universality of the quantum transfer

## Setting

- `p, q >= 2` are coprime, `N = (p-1)(q-1)/2 >= 3`, and `c = c_(p,q) = 1 - 6(p-q)^2/(pq)`.
- `V^c` is the universal Virasoro vertex algebra, and `v ∈ V^c` is the singular vector of weight `2N`.
- `I = U(Vir) v` is the maximal ideal: a single singular vector generates it (Feigin-Fuchs; Feigin-Frenkel 1993,
  section 3).
- `v` is normalized by `coeff(L_(-2)^N) = 1`. This coefficient is nonzero because Zhu's C_2-algebra is
  `R_L = C[x]/(x^N)`, i.e. by C_2-cofiniteness.
- Li filtration: `gr_F V^c = C[x_0, x_1, ...]` with `x_j <-> L_(-j-2)`, `wt x_j = j + 2`, and `L_(-1)` acting as
  `d x_j = (j+1) x_(j+1)`. The symbol of a vector is its top-length (number of modes) component.
- `J = (x_0^N)_∞` is the arc ideal, and `ker π = grI / J` is the classical defect (kind-pasteur S1).
- `K` is the Koszul complex of `f_n = d^n x_0^N` (`n >= 0`) on `C[x]`: the derived arc space of `Spec C[x]/(x^N)`
  (kind-pasteur S1, section 1.3).

**The transfer.** Put `u_n = L_(-1)^n v`, so that `σ(u_n) = f_n`.
- For a Koszul 1-cycle `z = Σ_n g_n e_n` of weight `D` and length `r`, i.e. `Σ g_n f_n = 0`, put
  `Φ(z) = Σ_n PBW(g_n) u_n ∈ I_D`, with `PBW` the ordered lift `x_(j_1)...x_(j_k) -> L_(-j_1-2)...L_(-j_k-2)`.
- The length-`r` component of `Φ(z)` is `Σ g_n f_n = 0`.
- **`T_1(z)` := the length-`(r-1)` component of `Φ(z)`.**

## Theorem

**(a) Well-definedness and linearity.**
- `T_1(z)` lies in `grI_(D,r-1)` and is well defined modulo `J`: other lifts or orderings change it by elements of `J`.
- `T_1` is `C[x][d]`-linear modulo `J`: `T_1(x_j z) ≡ x_j T_1(z)` and `T_1(dz) ≡ d T_1(z)`.

**(b) Sub-leading part of the singular vector.** The length-`(N-1)` part of `v` is
`α L_(-4) L_(-2)^(N-2) + β L_(-3)^2 L_(-2)^(N-3)`, with
```
α = -2N(N-1)/3 - Nc/12,          β = 11N(N-1)/48 + 5Nc/96.
```

**(c) Universality.** `T_1 = (c - c_(2,2N+1)) T*_N (mod J)`, where `T*_N` depends only on `N`.

**(d) Consequences.**
- All non-free minimal models with the same `N` have the same first-order image `Im T_1 = Im T*_N` (mod `J`). This
  image contains `b_N = x_0^(N-2) x_1 x_2 + ((2N-7)/12) x_0^(N-3) x_1^3`, the transfer of the Euler syzygy
  (THM-4671).
- `T_1 ≡ 0 (mod J)` iff `p = 2`.

## Proof

**(a)**
- *Lifts.* Changing the lift of `g_n` by an operator `δ` of length `<= deg g_n - 1` changes `Φ(z)` by `δ u_n`. Its
  length-`(r-1)` component is `σ(δ) f_n ∈ J`. Reordering PBW monomials produces commutator terms of lower length, and
  the same argument applies.
- *Membership in grI.* If `T_1(z) ≠ 0`, it is the top symbol of `Φ(z) ∈ I`.
- *Linearity in `x_j`.* `PBW(x_j g) - L_(-j-2) PBW(g)` has lower length.
- *Linearity in `d`.* `L_(-1)` acts on symbols as `d`, and `L_(-1) Φ(z) - Φ(dz)` is a sum of lifts of the same kind.
  Hence the correction terms lie in `J`.

**(b)** Write `v = L_(-2)^N + α L_(-4)L_(-2)^(N-2) + β L_(-3)^2 L_(-2)^(N-3) + (length <= N-2)`. Project `L_1 v = 0`
and `L_2 v = 0` to their length-`(N-1)` parts. Neither operator, applied to a vector of length `<= N-2`, produces
length `N-1`.
- **`L_1 v`, monomial `L_(-3)L_(-2)^(N-2)`.**
  - `L_1 L_(-2)^N|0> = Σ_i 3 L_(-2)^i L_(-1) L_(-2)^(N-1-i)|0>`, and `L_(-1) L_(-2)^k|0> = k L_(-3) L_(-2)^(k-1) + (lower)`.
    This contributes `3N(N-1)/2`.
  - `[L_1, L_(-4)] = 5 L_(-3)` contributes `5α`.
  - `[L_1, L_(-3)] = 4 L_(-2)`, applied to `L_(-3)^2`, contributes `8β`.
  - So `3N(N-1)/2 + 5α + 8β = 0`.
- **`L_2 v`, monomial `L_(-2)^(N-1)`.**
  - `[L_2, L_(-2)] = 4L_0 + c/2` gives `Σ_i (8(N-1-i) + c/2) = 4N(N-1) + Nc/2`.
  - `[L_2, L_(-4)] = 6L_(-2)` gives `6α`.
  - The `β` term only produces length `<= N-2`, because `[L_2, L_(-3)] = 5 L_(-1)` lowers the length.
  - So `4N(N-1) + Nc/2 + 6α = 0`.

Solving gives the formulas. Check: `(3,5)` (`N = 4`, `c = -3/5`) gives `α = -39/5`, `β = 21/8`, as in the exact singular
vector.

**(c)** `Φ(z) = Σ PBW(g_n) u_n`. The length-`(r-1)` part of `PBW(g_n) u_n` has two sources:
- `PBW(g_n)` applied to the length-`N` part of `u_n`, keeping the commutator terms of length `r - 1`;
- `g_n` times the length-`(N-1)` part of `u_n`.

`u_n = L_(-1)^n v`. Its length-`N` part is `f_n`, and its length-`(N-1)` part is `R_n + d^n σ_1(v)`:
- `R_n` comes from `L_(-1)`-reordering applied to `L_(-2)^N` and is independent of `c`;
- `σ_1(v) = α x_2 x_0^(N-2) + β x_1^2 x_0^(N-3)`.

All other ingredients (the PBW lifts, `J`, `K`) depend only on `N`. So `T_1 = A_N + α B_N + β C_N` (mod `J`), with
universal `A_N, B_N, C_N`. By (b) this is affine in `c`: `T_1 = A'_N + c B'_N`.

Now `Vir_(2,2N+1)` has the same `N`, and it is classically free (van Ekeren-Heluani), so `ker π = 0` there. Since
`Im T_1 ⊆ ker π`, we get `A'_N + c_(2,2N+1) B'_N ≡ 0 (mod J)`. Hence `T_1 = (c - c_(2,2N+1)) B'_N`. Put `T*_N = B'_N`.

**(d)**
- The images coincide because the scalar `c - c_(2,2N+1)` is nonzero for `p >= 3`: `c_(p,q) = c_(2,2N+1)` forces
  `p/q = 2/(2N+1)`.
- `T*_N ≠ 0`, because THM-4671's coefficient on the Euler syzygy is `N(2N+1)(c - c_(2,2N+1))/8`. ∎

## Remarks

- **The spectral sequence.** `Φ` is a filtered surjection onto `I` whose associated graded is the Koszul differential
  `K_1 -> C[x]` (image `J`). The spectral sequence of this two-term filtered complex has
  `E_1 = (Z_1(K), C[x]/J)` and converges to `(gr ker Φ, gr L)`.
  - Its differentials `d_k` are the higher transfers.
  - `ker π = grI/J` is filtered by transfer order, with order-1 piece `Im T_1`.
  - `E_1` is trivially model-independent. THM-4695 says that `d_1 = (c - c_(2,2N+1)) T*_N` is universal up to a
    scalar, so `E_2` is universal for `p >= 3`. The model enters only from `d_2` on.
  - The S44 data on transfer orders are HYP-9315.
- **Relation to THM-4671.** THM-4671 computes `T_1` on one class. THM-4695 shows that the whole of `T_1` carries the
  same factor `c - c_(2,2N+1)`. Three inputs then give `Im T_1 = C[x][d] b_N` modulo `J`:
  - `H_1` cyclic, generated by the Euler syzygy (HYP-9315(a)). This is KNOWN: a theorem for N = 2 (Kanade; Bai-Gorsky-
    Kivinen; Gorsky-Oblomkov-Rasmussen), and a consequence of Feigin 2011 Prop. 3.6 for all N, which omits one step.
  - `T_1(B_1) ⊆ J` (HYP-9315(d)), now proved; see HYP-9315.
  - Part (a) above.
  It was observed directly for Ising, `(3,5)` and `(4,5)` through weight 22.
- **Alternative proof of the vanishing at `c_(2,2N+1)`, without van Ekeren-Heluani** (audit m10). At `c_0`, the
  level-3 singular vector of `M(c_0, 2N)` lies in `ker Φ` with symbol `s_1`, so `T_1(s_1) ∈ J`. Then use HYP-9315(a)
  (Feigin), (d) and part (a).
