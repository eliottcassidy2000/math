---
id: THM-4848
title: "Multiplier descent for run observers of circulant tournaments. Let T(S) be a circulant tournament on Z_q, q an odd prime, with multiplier stabiliser H = {m : mS = S}. Every run observer M_k(T) splits over the characters of Z_q into blocks B_j. (i) det M_k = det B_0 * N(det B_1). (ii) Each multiplier m in H conjugates B_1 to its Galois image by a monomial matrix, so det B_1 lies exactly in the fixed field K_H (no unit ambiguity) and det M_k = det B_0 * N_{K_H/Q}(det B_1)^|H|. (iii) Away from q, the nontrivial part of coker M_k is a Z[zeta_q]-module whose elementary divisors are constant on H-orbits of primes, so every Z/p^e occurs with multiplicity divisible by |H<p>|. For Paley P_q, K_H = Q(sqrt(-q)): odd-exponent primes of the quadratic norm split or ramify, and primes p = 1 (mod q) come in blocks of (q-1)/2. In P11, 89 and 62987 occur as (Z/p)^5, while in the controls R7 and R11 (H = 1) such primes occur once. Checked for P7 (k = 2, 3), P11, P19 (k = 2) and the <3>-cyclotomic class on Z_13. The local P7 three-run modules of CLOCK_MINIMAL_PAIRING C5-C6 are the case q = 7, k = 3."
status: >
  PRE-AUDIT checkpoint: independent audit pending.
  PROVED (elementary; mac-mini-2026-10-29-cycle101): (i)-(iii) and the Paley specialisation. The method (circulant
  block diagonalisation, then Galois descent through a group of automorphisms) is standard; the statement for run
  observers of every circulant tournament, every k, is assembled here.
  FINITE-EXACT (04-computation/experiments/cycle101_tournament_collatz_20261029):
  - exact block determinants, the factorisation in (i) and the invariance in (ii) for P7 (k = 2, 3), P11 (k = 2) and
    the four labellings of the <3>-cyclotomic class on Z_13 (k = 2);
  - the invariance in (ii) and the ninth-power law for P19 (k = 2), with a partial factorisation of the 310-digit
    quadratic norm;
  - P7 with k = 4 (840 x 840): beta, its content, and the local Smith data at the split primes 71, 79, 151, 2843 and
    7529 (PARI mathnfmodid; galois_local.py), with non-trivial instances of (iii) at 71 and 2843;
  - full Smith forms of M_k and B_0 and the multiplicity law (iii) for all of these, with the controls R7 and R11
    (H = 1). Large elementary divisors are tested through a gcd-free basis, primes below 10^7 split off.
  INHERITED, priority elsewhere: CLOCK_MINIMAL_PAIRING C5-C6 (overnight_research_20261029) certified the p-primary
  pieces of M_3(P7) as C7-clock modules. At 11, 79, 5861 and 392809 the source occupies one cubic Frobenius orbit and
  none of its dual, the period field is Q(sqrt(-7)), and the 17-part is trivial. Those are this theorem's statements
  at q = 7, k = 3, and they used THREE_RUN_GLUE's Smith certificate.
  NEW HERE: the general descent and multiplicity law, the exactness of det B_1 in K_H, the norm formula, (v) the
  converse conjugation (the occupied side is a chirality signal), and the P11, P19, Z_13 and control checks.
session: mac-mini-2026-10-29-cycle101
source: 05-knowledge/results/cycle101_tournament_collatz_dictionary_20261029.md
depends_on:
  - 01-canon/theorems/THM-4846-run-observers-of-tournaments-the-nauru-orientation-the-permanental-compound-t5-stars-and-the-ghost-lemma-behind-its-four-spurious-z6-factors.md
  - 05-knowledge/results/overnight_research_20261029/CLOCK_MINIMAL_PAIRING.md
  - 05-knowledge/results/overnight_research_20261029/THREE_RUN_GLUE.md
related:
  - 01-canon/theorems/THM-4849-circulant-tournaments-are-rational-collatz-cycles-shift-decimation-reversal-two-base-points-and-the-locking-theorem-behind-35-and-175.md
  - 05-knowledge/hypotheses/HYP-9446-circulant-collatz-dictionary-and-multiplier-descent-open-questions.md
scripts:
  - 04-computation/experiments/cycle101_tournament_collatz_20261029/p7_fourier.py (P7, k = 3)
  - 04-computation/experiments/cycle101_tournament_collatz_20261029/paley_fourier.py (P7, P11, k = 2)
  - 04-computation/experiments/cycle101_tournament_collatz_20261029/paley_fourier2.py, p19_norm.py (P19, k = 2)
  - 04-computation/experiments/cycle101_tournament_collatz_20261029/galois_smith.py (Smith forms, multiplicity law, controls)
  - 04-computation/experiments/cycle101_tournament_collatz_20261029/chirality.py, chirality_more.py, converse_conjugate.py (sides, (v))
  - 04-computation/experiments/cycle101_tournament_collatz_20261029/galois_local.py (P7, k = 4, local Smith data)
  - helpers from 04-computation/experiments/growth_gf_runs_fermat_20261029 (runs.py, fast_runs.py)
outputs:
  - 04-computation/experiments/cycle101_tournament_collatz_20261029/out/ (p7_fourier.json, paley_fourier.out/.json, paley_fourier2.json, p19_norm.out, galois_smith.out/.json)
---

# THM-4848 — multiplier descent for run observers

## 0. Setting

- *The observer.* `M_k(T)` is THM-4846's run observer. Its rows and columns are ordered `k`-tuples of distinct
  vertices, and `M_k[s, t]` counts the ordered covers by `k` directed paths, the `i`-th from `s_i` to `t_i`.
- *Automorphisms.* Every automorphism `g` of `T` gives `M_k[gs, gt] = M_k[s, t]`.
- *The circulant.* `T = T(S)` is circulant on `Z_q` (`q` an odd prime), `τ(x) = x + 1`, and
  `H = {m ∈ Z_q^* : mS = S}`. `H` has odd order and `−1 ∉ H`. Each `m ∈ H` gives the automorphism `μ_m(x) = mx`,
  and `μ_m τ μ_m⁻¹ = τ^m`.
- *The blocks.* `τ` acts freely on ordered tuples. Fix representatives `R` of its orbits, `ζ = e^(2πi/q)`, and put

      B_j[r, r'] = Σ_(s ∈ Z_q) M_k[r, τ^s r'] ζ^(js)    (|R| × |R|, entries in Z[ζ]),

  so `B_j = σ_j(B_1)` for `j ≠ 0`, with `σ_j : ζ ↦ ζ^j`.

## 1. Statements

**(i) Block factorisation.** `M_k` is similar over `Q(ζ)` to `⊕_j B_j`, and `det M_k = det B_0 · N_(Q(ζ)/Q)(det B_1)`.

**(ii) Descent.** For `m ∈ H` there are a permutation matrix `P` and a diagonal matrix `Δ` of `q`-th roots of unity
with `B_m = Δ P B_1 P⁻¹ Δ⁻¹`. Hence:
- `σ_m(det B_1) = det B_1` *exactly*, so `β := det B_1` lies in `O_(K_H)`, `K_H = Q(ζ)^H`;
- `det M_k = det B_0 · N_(K_H/Q)(β)^|H|`.

**(iii) Multiplicity law.**
- *Splitting.* For a prime `p ≠ q`, `coker(M_k) ⊗ Z_(p) ≅ coker(B_0) ⊗ Z_(p) ⊕ X_p`, where `X_p` is the cokernel of
  `B_1` over `Z_(p)[ζ]`.
- *Local pieces.* As a group, `X_p` is `⊕_𝔓 ⊕_i O_𝔓/𝔓^(e_i(𝔓))` over the primes `𝔓 | p` of `Q(ζ)`, each piece
  `≅ (Z/p^(e_i))^f`, where `f = ord_q(p)`.
- *The law.* The list `(e_i(𝔓))` is constant on `H`-orbits of primes. Consequently, for every `e ≥ 1`, the number
  of `Z/p^e` summands of `X_p` is divisible by `|H⟨p⟩|`, the order of the subgroup of `Z_q^*` generated by `H`
  and `p`.
- *When it says more than the residue degree.* The law is stronger than divisibility by `f` exactly when
  `H ⊄ ⟨p⟩`.

**(iv) Paley.** For `q ≡ 3 (mod 4)` and `S = QR_q`: `H = QR_q`, `K_H = Q(√−q)`, and

    det M_k(P_q) = det B_0 · N_(Q(√−q))(β)^((q−1)/2).

- *Odd exponents.* A prime with odd exponent in `N(β)` splits or ramifies in `Q(√−q)`, i.e. `(−q/p) ≠ −1`. An
  inert prime has norm `p²`.
- *Block sizes.* Each `Z/p^e` of `X_p` has multiplicity divisible by `(q−1)/2` when `(p/q) = 1`, and by `q − 1`
  otherwise.
- *Where this is non-trivial.* It beats `f` exactly at the primes whose order mod `q` is neither `(q−1)/2` nor
  `q − 1`; for `q = 7` and `q = 11` these are `p ≡ ±1 (mod q)`. The converse symmetry `x ↦ −x` only relates `M_k`
  to `M_k^T`, so `𝔭` and `𝔭̄` remain independent: chirality is allowed.

**(v) The converse conjugates (answers the question left open by the 79 sweep).**
- *The identity.* For every circulant `T(S)` and every `k`, `det B_1(M_k(T(−S))) = σ_(−1)(det B_1(M_k(T(S))))`.
- *For Paley.* `β(P_q^op) = β̄(P_q)`, so passing to the converse labelling swaps `𝔭 ↔ 𝔭̄`.
- *Orientation.* The labellings `T(QR)` and `T(NQR)` differ by `x ↦ −x`, which replaces `τ` by `τ⁻¹`. No
  automorphism of `P_q` does this, since `−1 ∉ QR`. So which prime of `Q(√−q)` carries `X_p` is an invariant of `P_q`
  with an orientation of its regular cyclic subgroup, i.e. a chirality signal, and the converse flips it.
- *The P7 example.* For the three-run observer, the occupied cubic of CLOCK_MINIMAL_PAIRING C8 at 79,
  `X³+13X²+12X−1`, becomes its dual `X³+67X²+66X−1` under the converse labelling. The translation's eigenvalues
  are inverted.
- *Checked* (`converse_conjugate.out`): P7 (`k = 2, 3`) and P11 (`k = 2`). For example
  `β(QR) = 11904 − 23232η` and `β(NQR) = 35136 + 23232η = β̄` for P7, `k = 2`.

## 2. Proofs

**(i).** Functions on tuples form a free `Z[⟨τ⟩]`-module with basis `R`, and `M_k` is linear over it, given by the
matrix `M̂[r, r'] = Σ_s M_k[r, τ^s r'] τ^s`. Then `Q[⟨τ⟩] ≅ Q × Q(ζ)` (`τ ↦ (1, ζ)`), and the Galois conjugates of
the second factor give the other characters. ∎

**(ii).**
- *Where `μ_m` sends the representatives.* Write `μ_m(r) = τ^(c(r)) π(r)`, with `π(r) ∈ R`; `π` is a permutation of
  `R` because `μ_m` permutes `τ`-orbits.
- *The entries.* By invariance and `μ_m τ^s = τ^(ms) μ_m`,

      M_k[r, τ^s r'] = M_k[π r, τ^(ms + c(r') − c(r)) π r'].

  Substituting `u = ms + c(r') − c(r)` in the definition of `B_1` gives

      B_1[r, r'] = ζ^(m⁻¹c(r)) · B_(m⁻¹)[πr, πr'] · ζ^(−m⁻¹c(r')).

  This is the stated monomial conjugation, for `m⁻¹ ∈ H`. So `det B_1 = det B_(m⁻¹) = σ_(m⁻¹)(det B_1)`.
- *The norm formula.* Transitivity of the norm gives
  `N_(Q(ζ)/Q)(β) = N_(K_H/Q)(N_(Q(ζ)/K_H)(β)) = N_(K_H/Q)(β^|H|)`. ∎

**(iii).**
- *Splitting.* `q` is a unit in `Z_(p)`, so `Z_(p)[⟨τ⟩] ≅ Z_(p) × Z_(p)[ζ]` as rings, because `x − 1` and `Φ_q` are
  comaximal (`Φ_q(1) = q`). Hence `coker M̂` splits as stated.
- *Local pieces.* After `p`-adic completion, `Z_p[ζ] = Π_𝔓 O_𝔓` is unramified with uniformiser `p`, so the
  elementary divisors at `𝔓` describe `X_p`.
- *Galois transport.* `σ_m` carries the divisors of `B_1` at `𝔓` to those of `σ_m(B_1) = B_m` at `σ_m𝔓`. By (ii)
  these are the divisors of `B_1` at `σ_m𝔓`.
- *Orbit sizes.* The decomposition group of `𝔓` is `⟨p⟩`, and the Galois group is abelian. So every `H`-orbit of
  primes has size `|H|/|H ∩ ⟨p⟩| = |H⟨p⟩|/f`. The count of `Z/p^e` is then
  `f · Σ_orbits |orbit| · n_e = |H⟨p⟩| · Σ_orbits n_e`.
- *Paley.* For `q ≡ 3 (mod 4)`, `(p/q) = (−q/p)` by reciprocity, and `QR_q⟨p⟩ = QR_q` or `Z_q^*`.
- *The converse.* `T^op = T(−S) ≅ T` via `x ↦ −x`, and covers of `T^op` are reversed covers of `T`. Hence
  `Π M_k Π⁻¹ = M_k^T`. In blocks this says `B_(−1) ~ B_(−1)^T`, which leaves the Smith data at `𝔭` and `𝔭̄`
  unrelated. ∎

**(v).** Covers of `T^op` from `s` to `t` are the reversed covers of `T` from `t` to `s`, so
`M_k(T^op)[s, t] = M_k(T)[t, s]`. With the same representatives `R`,

    B_j(T^op)[r, r'] = Σ_s M_k(T)[τ^s r', r] ζ^(js) = Σ_u M_k(T)[r', τ^u r] ζ^(−ju) = B_(−j)(T)[r', r].

Hence `det B_1(T^op) = det B_(−1)(T) = σ_(−1)(det B_1(T))`, and `T^op = T(−S)`. ∎

## 3. Data (FINITE-EXACT)

| case | size | det `B_0` | quadratic (or `K_H`) norm of `β` | non-trivial instances of (iii) |
|---|---:|---|---|---|
| P7, k = 2 | 42 | `3^3·7·13^3` | `2^14·3^2·7·1451` (`1451 ≡ 2`, split) | none possible at the primes present |
| P7, k = 3 | 210 | `−2^14·3^18·17^4` | `2^22·5^4·11^2·79·5861^2·392809` | none (C5–C6's modules) |
| P7, k = 4 | 840 | (local data only) | `β = a + bη` with `gcd(a, b) = 2^41·5^3·7^3·71^2`; split odd primes 71, 79, 151, 2843, 7529 | **71, 2843 ≡ 1 (mod 7): `X_p = (Z/71)^15`, `(Z/2843)^9`**, i.e. blocks of 3 as forced; also `(Z/79²)^3`, `(Z/151)^9`, `(Z/7529)^6` |
| P11, k = 2 | 110 | `3·11^2·19^3·23·43·241^2·22291^2` | `2^12·5^2·11^2·31·53·89·62987·665251·35096521·2871099599·4780984157` | **89, 62987 ≡ 1 (mod 11): `X_p = (Z/p)^5`** |
| P19, k = 2 | 342 | `−3^12·5^2·17^2·23^2·41^2·43·1097^2·8941·9661·14843·19421·104659^2·c51^2` | `2^8·5^2·11·47·61·263·c299` | (Smith form not computed) |
| `Z_13`, `H = ⟨3⟩`, k = 2 | 156 | (4 labellings, identical data) | quartic `K_H`; nontrivial primes 3, 1907, 2011 and a 170-digit coprime block | none (all have `3 | f`) |
| R7 control, `H = 1` | 42 | — | — | 29, 71, 337, 303073, 2277647 ≡ 1 (mod 7) each occur **once** |
| R11 control, `H = 1` | 110 | — | — | nine primes ≡ 1 (mod 11) each occur **once** |

Notes on the table:
- **P7, k = 3.** `β = 3257463411200 − 12008276134400 η` with `η = ζ + ζ² + ζ⁴ = (−1 + √−7)/2`. Its norm
  `a² − ab + 2b²` is the second column's entry, so `N(det B_1) = 2^66 5^12 11^6 79^3 5861^6 392809^3`.
- **The 2-part of P7, k = 3.** The nontrivial 2-part is `(Z/2)^21 ⊕ (Z/8)^3 ⊕ (Z/16)^9`, with multiplicities
  divisible by 3.
- **P11, 2-part.** 2 is inert in `Q(√−11)` and has even exponent. Its nontrivial 2-part `(Z/2)^30 ⊕ (Z/8)^10` has
  multiplicities divisible by `f = 10`.
- **P11, odd-exponent primes.** All are `≡ 1, 3, 4, 5, 9 (mod 11)`, i.e. split.
- **P19.** `σ_4(β) = β`, where 4 generates `QR_19`. The absolute norm is exactly the ninth power of the quadratic
  norm. The primes 11, 47, 61 and 263 have `(−19/p) = 1`, and 2 is inert with exponent 8. `c299` is a composite
  299-digit cofactor of the quadratic norm, and `c51` a composite 51-digit cofactor of `det B_0` (gp `isprime` = 0).
- **The controls.** They show that the repetition in P11 comes from `H`. Without multipliers, primes `≡ 1 (mod q)`
  enter `X_p` singly. That every odd prime in the controls' `X_p` is `≡ 1 (mod q)` is the generic behaviour of
  norms from `Q(ζ_q)`: a degree-`f` prime costs `p^f`. It is not a theorem here.
