---
id: THM-4671
title: "First quantum correction to the arc space of a Virasoro minimal model: for coprime 2 <= p < q, N = (p-1)(q-1)/2 and V = Vir_{p,q}, the arc-ideal syzygy x_0 d(x_0^N) = N x_1 x_0^N lifts to u = L_{-2}L_{-1}v - N L_{-3}v in the maximal ideal, whose symbol is A b_N + B d^3(x_0^N) with the universal relation b_N = x_0^{N-2} x_1 x_2 + ((2N-7)/12) x_0^{N-3} x_1^3, A = N(2N+1)(c - c_{2,2N+1})/8 = 3(p-1)(p-2)(q-1)(q-2)(p^2 q^2 - (p-q)^2)/(16pq) and B = -(6N+c-10)/24; hence ker(J_infty R_V -> gr_F V) vanishes below weight 2N+3 and, for p >= 3, is spanned at weight 2N+3 by b_N"
status: >
  PROVED (elementary, from standard inputs).
  CITED inputs: Li's filtration of a Virasoro vertex algebra is the PBW-length filtration and gr_F V^c = C[x_0, x_1, ...]
  (Li 2005; Arakawa 2012; Andrews-van Ekeren-Heluani 2020 sec. 2); the maximal ideal of V^c at c = c_{p,q} is generated
  by one singular vector v of weight (p-1)(q-1) whose L_{-2}^N coefficient is nonzero, so R_V = C[x]/(x^N)
  (Feigin-Fuchs; Wang 1993; van Ekeren-Heluani 2021); the Rocha-Caridi character formula.
  The proof itself is the two-line degree count in Lemmas 1-2 below plus algebra.
  FINITE-EXACT confirmation: exact rational computation of u and its symbol for 12 models
  (3,4) (3,5) (4,5) (3,7) (3,8) (4,7) (5,6) (5,7) (3,10) and the p = 2 controls (2,7) (2,9) (2,11), all matching A and B
  (first_correction_check.py); the mod-p singular-support engine finds a minimal non-arc generator at weight 2N+3 in
  every p >= 3 model scanned and none for p = 2 (vir_ss.py, tables in the source note).
  Scope: this is the FIRST layer of ker(pi) only. It does not describe higher generators, finite generation, or the
  full singular support (see HYP-9294). Novelty: for N = 3 (Ising) the relation is Andrews-van Ekeren-Heluani's b
  (2020) up to the arc ideal; the closed form of A(p,q), the universality of b_N in N, and the lift mechanism were
  not found in van Ekeren-Heluani (2021), AEH (2020), Salazar (2024/26) or Arakawa-Linshaw (2019) (abstract/intro
  level search, 2026-10-18). Novelty is unverified beyond that search.
session: kind-pasteur-2026-10-18-S1
source: 05-knowledge/results/singular_support_arc_spaces_classical_freeness_20261018.md
scripts:
  - 04-computation/voa_singular_support/first_correction_check.py: exact check of A, B for 12 models
  - 04-computation/voa_singular_support/vir_ss.py: singular-support engine (grI, arc ideal J, minimal generators)
related: [HYP-9293, HYP-9294, Andrews-van Ekeren-Heluani 2020 (arXiv:2005.10769), van Ekeren-Heluani 2021 (CMP), Arakawa-Linshaw 2019]
---

# THM-4671 — the first quantum correction to the arc space of `Vir_{p,q}`

## Setting

- `p < q` are coprime integers `>= 2`, `c = c_{p,q} = 1 - 6(p-q)^2/(pq)`, `N = (p-1)(q-1)/2`.
- `V^c` is the universal Virasoro vertex algebra with PBW basis `L_{-k_1} ... L_{-k_r}|0>`, `k_1 >= ... >= k_r >= 2`.
  `G_r` is the span of monomials with at most `r` factors.
- **Li filtration = PBW filtration (CITED).** The Li degree of a monomial of weight `D` with `r` factors is `D - 2r`. So
  `gr_F V^c = ⊕ G_r/G_(r-1) = C[x_0, x_1, ...]`, with `x_j <-> L_(-j-2)` of weight `j + 2`. The translation `L_(-1)`
  becomes the derivation `d x_j = (j+1) x_(j+1)`.
- `v` is the singular vector of weight `2N`, normalized so that `L_{-2}^N` has coefficient 1. `I = U(Vir) v` is the
  maximal ideal and `V = V^c / I = Vir_{p,q}`.
- `R_V = C[x]/(x^N)`, `J = (x_0^N)_d` (the differential ideal), `C[x]/J = J_infty R_V` is the arc algebra, and
  `gr_F V = C[x]/grI`, where `grI` is the span of the top-length components (symbols) of elements of `I`.
- `pi : J_infty R_V -> gr_F V` is the canonical surjection, with `ker pi = grI / J`. `V` is **classically free** when
  `ker pi = 0` (van Ekeren-Heluani).
- Write `f = x_0^N`, `m_1 = x_0^(N-1) x_3`, `m_2 = x_0^(N-2) x_1 x_2` and `m_3 = x_0^(N-3) x_1^3`. These are the
  only monomials of weight `2N+3` and degree `N`.

## Statement

**Theorem.** Assume `N >= 3`, which holds whenever `p >= 3`. Put
```
b_N = x_0^(N-2) x_1 x_2 + ((2N-7)/12) x_0^(N-3) x_1^3,
A   = N(2N+1)(c - c_(2,2N+1))/8 = 3(p-1)(p-2)(q-1)(q-2)(p^2 q^2 - (p-q)^2) / (16 p q),
B   = -(6N + c - 10)/24.
```
1. `u := L_(-2) L_(-1) v - N L_(-3) v` lies in `I ∩ G_N`, and its symbol in `G_N/G_(N-1)` is
   `sigma_N(u) = A b_N + B d^3 f`.
2. `ker pi` is zero in every weight `< 2N+3`.
3. In weight `2N+3`, `ker pi` has dimension `1` if `p >= 3` and `0` if `p = 2`. For `p >= 3` it is spanned by the class
   of `b_N`.

So `b_N` is a relation of `gr_F Vir_{p,q}` that is not in the arc ideal, and its normal form depends only on `N`. The
model enters only through the scalar `A`, and `A` vanishes exactly on the boundary series `p = 2` (or `q = 2`).

## Proof

**Lemma 1 (subleading part of `v`).** Modulo `G_(N-2)`,
`v == L_(-2)^N + alpha L_(-4) L_(-2)^(N-2) + beta L_(-3)^2 L_(-2)^(N-3)`, with
```
alpha = -N(8(N-1) + c)/12,     5 alpha + 8 beta = -3N(N-1)/2.
```
*Proof.* The partitions of `2N` into `N-1` parts `>= 2` are `(4, 2^(N-2))` and `(3^2, 2^(N-3))`. Take the
`G_(N-1)/G_(N-2)` component of `L_1 v = 0`:
- `[L_1, L_(-2)] = 3 L_(-1)`, and `L_(-1) L_(-2)^m |0> == m L_(-3) L_(-2)^(m-1)|0>`. So `L_1 L_(-2)^N|0>` contributes
  `(3N(N-1)/2) L_(-3) L_(-2)^(N-2)`.
- `[L_1, L_(-4)] = 5 L_(-3)` and `[L_1, L_(-3)] = 4 L_(-2)` (twice). These contribute `5 alpha + 8 beta` on the same
  monomial.
- Terms of `v` in `G_(N-2)` cannot reach `G_(N-1)`, since `L_1` never raises length.

For `L_2 v = 0`, use `[L_2, L_(-2)] = 4L_0 + c/2`. Then `L_2 L_(-2)^N|0>` contributes `(4N(N-1) + Nc/2) L_(-2)^(N-1)`,
and `[L_2, L_(-4)] = 6 L_(-2)` contributes `6 alpha`. The `beta` term only produces `L_(-1)` insertions, which lower
the length. `□`

**Lemma 2 (the lift).** `u` has no length-`(N+1)` part, and
`sigma_N(u) = (N(N+1)/2 + 3 alpha) m_1 + (4 beta - 2 alpha) m_2 - 3 beta m_3`.

*Proof.* `L_(-1) L_(-2)^N|0> == N L_(-3) L_(-2)^(N-1) + (N(N-1)/2) L_(-5) L_(-2)^(N-2)` modulo `G_(N-2)`. This uses
`L_(-2) L_(-3) = L_(-3) L_(-2) + L_(-5)`. Left multiplication by `L_(-2)` then gives
`N L_(-3) L_(-2)^N + (N(N+1)/2) L_(-5) L_(-2)^(N-1)` modulo `G_(N-1)`.

The term `N L_(-3) L_(-2)^N` cancels against `-N L_(-3) v`. In the symbol `x_0 d f - N x_1 f = 0`, this cancellation is
exactly the arc syzygy.

The `alpha` and `beta` parts contribute their naive symbols:
- `alpha (x_0 d(x_2 x_0^(N-2)) - N x_1 x_2 x_0^(N-2))`;
- `beta (x_0 d(x_1^2 x_0^(N-3)) - N x_1^3 x_0^(N-3))`.

Expanding with `d x_j = (j+1) x_(j+1)` gives the stated coefficients. `□`

**Proof of 1.**
`d^3 f = 6N m_1 + 6N(N-1) m_2 + N(N-1)(N-2) m_3`. Substituting Lemma 1 into Lemma 2:
- the `m_1` coefficient fixes `B`;
- the `m_2` coefficient fixes `A`;
- the `m_3` coefficient is then satisfied identically in `(N, c)`. This is a polynomial identity, checked symbolically
  (residual 0).

Finally,
```
c_(2,2N+1) = 1 - 3(2N-1)^2/(2N+1) = -(12N^2 - 14N + 2)/(2N+1),
```
which gives the two forms of `A`.

**Proof of 2-3.** `I` is the image of the Verma module `M(c, 2N)`. By the Rocha-Caridi formula,
`char I = (q^(2N) - q^(2N+2p-1) - q^(2N+2q-1) + ...)/(q;q)_infty`, so `I_(2N+k)` has dimension `p(k)` for
`k < 2p-1`. Hence:
- `I_(2N) = <v>`.
- `I_(2N+1) = <L_(-1) v>`.
- `I_(2N+2) = <L_(-1)^2 v, L_(-2) v>`.
- `I_(2N+3) = <L_(-1)^3 v, L_(-3) v, u>` when `p >= 3`.

Their symbols are `f`, `d f`, `d^2 f`, `x_0 f`, `d^3 f`, `x_1 f` and `sigma_N(u)`.

The arc ideal `J` in these weights is spanned by `f`, `d f`, `{d^2 f, x_0 f}` and `{d^3 f, x_1 f}`, because
`x_0 d f = N x_1 f`.

So `grI = J` below weight `2N+3`, and in weight `2N+3`, `grI/J` is spanned by `sigma_N(u) == A b_N` modulo `J`.
- For `p >= 3`, `A != 0`, and `b_N` is not a multiple of `d^3 f` (its `m_1` coefficient is 0). So `b_N` spans the
  one-dimensional `ker pi_(2N+3)`.
- For `p = 2`, `A = 0`. Then `u` is the level-3 singular vector of the submodule (`2p-1 = 3`), and `ker pi_(2N+3) = 0`,
  consistent with classical freeness (van Ekeren-Heluani). `□`

## Checks, boundary, consequences

- **Ising (`N = 3`).** Andrews-van Ekeren-Heluani's generator `b = x_3 x_0^2/6 + x_2 x_1 x_0` equals
  `(1/108) d^3 f + (2/3) b_3`. So it is the same class.
- **Universality.** The 2-plane `grI_(2N+3, N) = span(d^3 f, b_N)` is the same for every `c` with `A(N, c) != 0`. It
  depends only on the `C_2`-algebra `C[x]/(x^N)`. The minimal model is recorded only by the scalar `A`.
- **Mechanism.** The relation is the quantum shadow of the first syzygy of the arc ideal: an `H_1` class of the Koszul
  complex of `(d^n f)_n`. In derived terms, a class in `pi_1` of the derived arc space of the fat point. Classical
  freeness at weight `2N+3` is exactly the vanishing of its transfer coefficient `A`.
- **Not claimed.** Generators of `ker pi` beyond weight `2N+3`, finite generation, and the bigraded character of
  `gr_F Vir_{p,q}` for `p >= 3`. The engine's data for these are in the source note and HYP-9294.
