---
id: HYP-9315
title: "The Li spectral sequence of a Virasoro minimal model. (a) [KNOWN] H_1 of the derived arc space of Spec C[x]/(x^N) is generated over C[x][d] by the Euler syzygy x_0 e_1 - N x_1 e_0 (modulo boundaries). (b) In degree N the minimal generators of grI sit at Li degrees {0, 3, 5, ..., 2p-3} (kind-pasteur's HYP-9294(a)), and the one at Li degree 2k+1 is a transfer of order exactly k. (c) Generators of degree N+i have transfer order <= p-2+i. (d) [PROVED] Koszul boundaries transfer into J. (e)-(g) Universality: the transfer-order tables T_(p,q)(L,i) agree for all models at relative level m = L+2i < min(2p-1, 2p'-1); for fixed p they agree for m < 2 min(q,q') - 1 (q-stability); the per-weight totals are forced by the characters, so the content is the refinement by length and order"
status: >
  Mixed: (a) KNOWN/CITED; (d) PROVED; (b), (c), (f), (g) OPEN with FINITE-EXACT evidence; (e) is a character identity.
  Independently audited: 04-computation/experiments/s44_audit_20261018/REPORT.md. Every number reproduced
  bidegree-for-bidegree by an independent engine at two primes, with 0 differences.
  - (a) KNOWN, missed by the session (MISTAKE-640):
    - N = 2 is a theorem: Kanade, thesis 2015 Thm 4.3.1 (JPAA 222 (2018)); Bai-Gorsky-Kivinen, Ramanujan J. 52 (2020)
      Thm 2.2; Gorsky-Oblomkov-Rasmussen, Exp. Math. 22 (2013) Cor. 3.20.
    - All N: a consequence of Feigin, Funkts. Anal. Prilozh. 45:4 (2011), Prop. 3.6, which omits one step.
    - Gorsky-Lewark, arXiv:1404.0623, Conj. 4, is the truncated all-N form.
    - The finite checks are now controls: S44 N = 2, 3, 4 to weight 24; audit N = 2..7.
  - (b) The real content is "order >= k" for k >= 3:
    - "order <= k" is forced by weight counting;
    - order >= 2 is forced by (a) + (d);
    - tested at k = 3 in (5,6), (5,7), (5,8), (5,9), (6,7) and at k = 4 in (6,7), including pre-registered runs that
      PASSED.
    - Given the positions, the generators themselves are forced, independent of c (THM-4710).
  - (c) order <= p-2+i holds in all data. Equality is attained for (p,i) = (3,1), (4,1), (5,1), (3,2) and is open for
    (6,1).
  - (d) PROVED; the proof is below (audit note N2).
  - (e) Character identity.
  - (f), (g) Universality, with 0 disagreements in all S44 tables:
    - 16 models; 388 + 452 bidegree comparisons and 2401 order counts (evaluate_corrected_windows.out);
    - the audit's 45 pairs agree.
    - The pre-registered windows P4/P6 (W(N) = {L <= N+i-1}) were WRONG: refuted by the audit's (5,6) run to
      weight 32, which S44 reproduced (MISTAKE-640).
session: opus-2026-10-18-S44
source: 05-knowledge/results/li_spectral_sequence_quantum_transfer_20261018.md
scripts:
  - 04-computation/quantum_transfer_s44/vir_transfer_orders_20261018.py (+ orders_*.out, vir_transfer_orders_*.json)
  - 04-computation/quantum_transfer_s44/koszul_h1_cyclicity_20261018.py (+ koszul_h1_cyclicity_N*_D24.json)
  - 04-computation/quantum_transfer_s44/order1_cyclic_check_20261018.py (+ order1_cyclic_*.json)
  - 04-computation/quantum_transfer_s44/vir_quantum_transfer_20261018.py
  - 04-computation/quantum_transfer_s44/evaluate_stability_prereg_20261018.py, evaluate_corrected_windows_20261018.py
  - 04-computation/experiments/s44_audit_20261018/ (independent engine, gens_orders_audit.py, stability_check.py, verma_relation_symbol.py)
related: [THM-4695, THM-4710, THM-4671, HYP-9294, HYP-9297, HYP-9293, Feigin 2011, Feigin-Frenkel 1993, Kanade 2018, Bai-Gorsky-Kivinen 2020, Gorsky-Oblomkov-Rasmussen 2013, Gorsky-Lewark 2014, Andrews-van Ekeren-Heluani 2023]
---

# HYP-9315 — the Li spectral sequence of a Virasoro minimal model

## Setting

See THM-4695. The filtered surjection `Φ: ⊕_n U(Vir_(<=-2)) e_n -> I`, `e_n -> L_(-1)^n v`, has associated graded
the Koszul differential of the derived arc space of `Spec C[x]/(x^N)`.
- Its spectral sequence has `E_1 = (Z_1(K), C[x]/J)` and converges to `gr L`.
- So `ker π = grI/J` carries the filtration by transfer order. Its `k`-th piece is `Im d_k`.
- **Precedent.** Feigin 2011 compares the associated graded of the BGG resolution of the `(2, 2N+1)` modules with
  the Koszul complex of jets. His complex `K(p)` is the `p = 2` instance of the `E_1` comparison here.
- The relative level of a bidegree (weight `D`, length `s = N+i`, Li degree `L = D - 2s`) is
  `m = D - 2N = L + 2i`.

**Computation.** `vir_transfer_orders_20261018.py` processes the Koszul 1-cycles of each weight in increasing length
and reduces `Φ(z)` top-down against the tops already found. The step at which a new symbol first appears is its
transfer order.
- The audit (N3) shows that the algorithm computes the canonical filtration
  `F_j = σ_s(A_(s+j) ∩ V^(<= s))`, with `A_r = Φ(chains of length <= r)`.
- So "order" is well defined.
- The order of a generator is read modulo decomposables.

## Statements

**(a) Cyclicity of `H_1` (KNOWN).**
- `H_1(K) = C[x][d] s_1 + B_1`, with `s_1 = x_0 e_1 - N x_1 e_0` of weight `2N + 3`.
- References: Kanade (N = 2); Bai-Gorsky-Kivinen (N = 2); Gorsky-Oblomkov-Rasmussen (N = 2, via Feigin); Feigin 2011
  Prop. 3.6 (all N, one step omitted); Gorsky-Lewark Conj. 4 (truncated, all N).
- With THM-4695(a) and (d) below, `Im T_1 = C[x][d] b_N` mod `J`.

**(b) Orders of the degree-`N` generators.**
- The minimal generators of `grI` in degree `N` are at Li degrees `0, 3, 5, ..., 2p-3` (HYP-9294(a)). The generator at
  Li degree `2k+1` is a transfer of order exactly `k`.
- What is forced:
  - A Koszul 1-cycle of length `N+j` has weight `>= 2N+2j+1`, and the unique minimal one is `x_0^(j-1) s_1`. So an
    order-`j` class has Li degree `>= 2j+1`, which forces "order `<= k`".
  - By (a) and (d), every generator other than `b_N` has order `>= 2`.
  - The genuine content is therefore "order `>= k`" for `k >= 3`.
- **Why the list stops at `2p-3`** (audit N5). At `m = 2p-1` the Verma map `M(c, 2N) -> I` acquires its first relation
  (BGG). There `ker Φ` is one-dimensional and its top symbol is exactly the Koszul cycle `x_0^(p-2) s_1`, checked for
  8 models. That cycle never transfers. For `m < 2p-1`, `Φ` is injective, so every cycle transfers.
- **The generators are forced** (THM-4710). Given the positions, the generator at Li `2k+1` is the unique lowest-weight
  vector admissible under the c-independent positive Witt action. So it depends only on `N`.

**(c) Higher degrees.**
- Generators of degree `N + i` have transfer order at most `p - 2 + i`.
- Equality is attained for `(p,i) = (3,1), (4,1), (5,1), (3,2)`, and open for `(6,1)` (maximum 4 through weight 44).
- It does not hold for every `i`: (3,5) has no degree-6 generator through weight 38.

**(d) Factorization (PROVED).** `T_1(B_1) ⊆ J`. Proof (audit N2):
- `B_1 = C[x][d]·{bd(e_0 ∧ e_n) : n >= 1}`, since `d(e_a ∧ e_b) = e_(a+1) ∧ e_b + e_a ∧ e_(b+1)` and `bd` commutes
  with `d`. By THM-4695(a) it suffices to treat `b = f_a e_b - f_b e_a`.
- (i) `L_k u_n ∈ C u_(n-k)` for `k >= 0`, since `v` is singular and `[L_k, L_(-1)] = (k+1) L_(k-1)`.
- (ii) `(u_a)_(-1)` is a sum of normally ordered monomials. So the length-`(2N-1)` part of `(u_a)_(-1) u_b` equals
  that of `PBW(f_a) u_b` mod `J`. Every other term is (at most `N-1` creation modes) times some `u_(b')`, and its
  length-`(2N-1)` part is (symbol)`·f_(b')`.
- (iii) Skew-symmetry:
  `(u_a)_(-1) u_b - (u_b)_(-1) u_a = Σ_(j>=1) (-1)^(j+1)/j! L_(-1)^j ((u_a)_(j-1) u_b)`,
  with `(u_a)_(m) = (-1)^a m^(falling a) v_(m-a)` and `v_(i) u_b = Σ_k C(b,k) i^(falling k) L_(-1)^(b-k) v_(i-k) v`.
- (iv) For `i >= 0`, `v_(i) v` has weight `4N - i - 1 < 4N`. So its terms have at most `N-1` creation modes from the
  length-`N` part of `v`, or come from lower-length parts of `v`. Either way, the length-`(2N-1)` part lies in
  `C[x]·{f_e} ⊆ J`, and `L_(-1)^j` preserves this.
- Hence `T_1(b) ∈ J`. The audit numerics agree, on all boundary generators for N = 3, 4, 6.

**(e) Forced totals (character identity).**
```
dim ker π(2N+m) = χ_(2,2N+1)(2N+m) - χ_(p,q)(2N+m) = p(m-3) - p(m-2p+1) - p(m-2q+1)
```
- This holds for `m < min(2(p+q), 4N+1)`. It comes from the BGG/Rocha-Caridi expansions: the Verma terms of `I` sit at
  relative levels `0, 2p-1, 2q-1, ...`, and for `(2,2N+1)` at `0, 3, 4N+1`.
- So equal totals of two models below `2 min(q,q') - 1` carry no information.

**(f) Universality and q-stability of the refined tables.** Write `T_(p,q)(L,i)` for the multiset of transfer orders at
length `N+i` and Li degree `L`.
- **(U)** `T_(p,q)(L,i) = T_(p',q')(L,i)` whenever `m = L + 2i < min(2p-1, 2p'-1)`.
- **(Q)** For fixed `p`, `T_(p,q)(L,i) = T_(p,q')(L,i)` whenever `m < 2 min(q,q') - 1`.
- Both windows are sharp. Every same-`p` pair first disagrees exactly at `m = 2 q_min - 1`:
  - (3,4) vs the rest at 7; (3,5) at 9; (3,7)/(3,8) at 13;
  - (4,5) at 9; (5,6)/(5,7) at 11, at `(L,i) = (9,1)`: `{3:2}` vs `{3:3}`.

**(g) Saturation in `p`.** For `k <= min(p,p') - 2`, the order-`k` counts agree whenever `m < 2 min(q,q') - 1`. The order
`p - 1` part is a proper part of the saturated count. For example, (9,1) has `{3:1}` for `p = 4` against `{3:3}` for
`p >= 5` with `q >= 7`.

**Sharper form (audit N5, OPEN).** The refined data at relative level `m` depend only on `{2p-1, 2q-1} ∩ [0,m]`, for
`m < 4N+1`.

## Evidence

The orders of all of `ker π`, through the stated weight:

| model | `N` | weight | order 1 | order 2 | order 3 | order 4 |
|---|---|---|---|---|---|---|
| (3,4) Ising | 3 | 30 | 511 | – | – | – |
| (3,5) | 4 | 30 | 516 | 34 | – | – |
| (3,7) | 6 | 28 | 143 | 28 | – | – |
| (3,8) | 7 | 30 | 143 | 33 | – | – |
| (3,10) | 9 | 30 | 45 | 7 | – | – |
| (3,11), (3,13) | 10, 12 | 32, 36 | 45 | 7 | – | – |
| (4,5) | 6 | 28 | 143 | 84 | 4 | – |
| (4,7), (4,9) | 9, 12 | 30, 36 | 45 | 30 | 3 | – |
| (4,11) | 15 | 42 | 45 | 30 | 3 | – |
| (5,6) | 10 | 32 | 45 | 30 | 13 | – |
| (5,7), (5,8) | 12, 14 | 36, 40 | 45 | 30 | 15 | – |
| (5,9) | 16 | 42 | 24 | 14 | 5 | – |
| (6,7) | 15 | 41 | 33 | 21 | 9 | 3 |

- **Pre-registration 1** (`PREREGISTRATION_orders_67_57.md`, commit 94fbe3b30b): PASSED. Vir_(6,7) orders 1-4 first
  appear at weights 33/35/37/39 at length 15, and Vir_(5,7) orders 1-3 at 27/29/31. Per the audit, only k = 3, 4 were
  informative.
- **Pre-registration 2** (`PREREGISTRATION_stability_20261018.md`, commit dcbb7ffa8a):
  - P5, P7, P8: PASSED.
  - P4, P6 as stated: REFUTED. The window `W(N) = {L <= N+i-1}` contains the forced break at `m = 2q-1`. S44's (5,6)
    data stopped at `m = 10`, so "all passed" was an artifact of the range.
  - On the corrected windows (f), (g): 0 disagreements in all data. This is post hoc and not pre-registered.
- **Generators.** kind-pasteur's degree-`N` generators have orders 1, 2, ..., `p-2` at Li 3, 5, ..., `2p-3` in every
  model; the audit checked this invariantly (gens_orders_audit.py). New generators found by the audit:
  - (5,6): (33,11) x2 {3,3}, (35,11) x2 {3,4};
  - (5,7): (37,13) x2 {3,3}, (38,13) {4}, (39,13) x3 {3,4,4};
  - (6,7): (41,16) {2}, (42,16) {3}, (43,16) x2 {3,3}, (44,16) {4}.

## Reading

- **Classical freeness and higher differentials.** Classical freeness is the vanishing of every `d_k`.
  - By THM-4695 and THM-4671, `d_1 = 0` iff `p = 2`.
  - For `p >= 3`, (b) predicts nonzero `d_1, ..., d_(p-2)` in degree `N`. The audit's structural reading (N5) explains
    why `d_(p-1)` is not needed there.
- **Universal and model-specific parts.** `E_2` is universal for each `N`. The model-specific content of the singular
  support lies in `d_2, d_3, ...`. Below `m = 2p-1` even these are universal (f), and in degree `N` they are rigid
  (THM-4710).
- **Literature.** Feigin-Frenkel 1993, Remark 2, notes that for general minimal models the symbols `S_n` "do not
  generate the ideal", and poses finding the other generators as a problem. (b), (f) and THM-4710 address it.

## Falsifiers

- A degree-`N` generator at Li degree `2k+1` whose transfer order is not `k`.
- A transfer of order `> p - 2 + i` in degree `N + i`.
- Two models whose refined tables differ at some `m < min(2p-1, 2p'-1)`. For the same `p`: at some
  `m < 2 min(q,q') - 1`.
