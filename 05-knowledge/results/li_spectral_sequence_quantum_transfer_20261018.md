# The Li spectral sequence: quantum transfer from derived arc spaces to the singular support

**Session:** opus-2026-10-18-S44. **Status:**
- PROVED: THM-4695 (one CITED input), THM-4710 (finite range), HYP-9315(d).
- KNOWN: HYP-9315(a), missed and repaired (MISTAKE-640).
- FINITE-EXACT / OPEN: HYP-9315(b), (c), (f), (g); HYP-9316.
- VERIFIED: section 6(i)-(ii).
- **Independently audited:** 04-computation/experiments/s44_audit_20261018/REPORT.md. THM-4695 is correct and every
  number reproduced, at two primes with a separate engine. Three MAJOR findings have been applied: M1 prior art for
  HYP-9315(a); M2 a false claim in section 6(iii); M3 a wrong pre-registered window. THM-4710 and section 8 were added
  after the audit.

**Owner prompt.**
> continue with another similar session, building on the ideas of other agents creatively

**What this builds on.**
- **kind-pasteur-2026-10-18-S1**,
  [singular_support_arc_spaces_classical_freeness_20261018.md](singular_support_arc_spaces_classical_freeness_20261018.md):
  - classical freeness;
  - the derived arc space of the fat point (section 1.3);
  - THM-4671, the first quantum correction;
  - the generator tables (HYP-9294);
  - open problem 6 ("test whether ker π is always the image of H_1 under a quantum transfer map; a spectral-sequence
    formulation of classical freeness");
  - the cusp "curiosity" (section 5a).
  Their engine `vir_ss.py` is reused unchanged.
- **mac-mini-2026-10-18-sixcolour**, [merge_18_19_11_22_collatz_hn_20261018.md](merge_18_19_11_22_collatz_hn_20261018.md),
  sections 2.2-2.2b: the κ table, the Hoffman margin 4.0000439 at `q = 31`, and the Delsarte/Schrijver LP.
- **S43**, [factorization_vertex_supersymmetric_bridge_20261018.md](factorization_vertex_supersymmetric_bridge_20261018.md):
  the derived critical locus and factorization-algebra framing.
- **Prior art located by the audit.**
  - Feigin, Funkts. Anal. Prilozh. 45:4 (2011). He compares the associated graded of the BGG resolution of the
    `(2,2N+1)` modules with the Koszul complex of jets. This is the closest precedent of sections 1 and 4.
  - Feigin-Frenkel 1993, Remark 2, poses the problem of finding the extra generators.
  - Andrews-van Ekeren-Heluani (Ising).
  - Kanade, Bai-Gorsky-Kivinen, Gorsky-Oblomkov-Rasmussen, Gorsky-Lewark (Koszul homology of jets of `x^N`).
  - Docampo-Miller-Overton-Walker 2026 (derived jets and arcs).
- Companion merge note (Pythagorean form, conic graphs, RP trick):
  [pythagorean_form_conic_graphs_rp_margin_20261018.md](pythagorean_form_conic_graphs_rp_margin_20261018.md).

## 1. The quantum transfer and its spectral sequence

Let `V^c` be the universal Virasoro vertex algebra at `c = c_(p,q)`. Let `v` be its singular vector of weight
`2N = (p-1)(q-1)` and `I` the maximal ideal. In the Li filtration, `gr V^c = C[x_0, x_1, ...]`, the arc ideal is
`J = (x_0^N)_∞`, and the classical defect is `ker π = grI/J`.

**The filtered surjection.**
```
Φ: ⊕_(n>=0) U(Vir_(<=-2)) e_n  ->  I,      Q e_n -> Q L_(-1)^n v,
```
filtered by naive length.
- Its associated graded is the Koszul differential `e_n -> f_n = d^n x_0^N` of the derived arc space of
  `Spec C[x]/(x^N)`.
- The spectral sequence has `E_1 = (Z_1(K), C[x]/J)` and converges to `(gr ker Φ, gr L)`.
- Feigin's `K(p)` (2011) is the `p = 2` instance of this comparison.

**Transfers.** A Koszul 1-cycle `z` lifts to `Φ(z) ∈ I`, whose naive top cancels.
- The first surviving component, one length lower, is `T_1(z)`.
- If it lies in `J` it can be corrected, and the next component gives `T_2`, and so on.
- `ker π` is filtered by transfer order. The audit (N3) shows this filtration is canonical.

**Classical freeness** means that every `d_k` vanishes.

## 2. THM-4695: `d_1` is universal up to a scalar

See [THM-4695](../../01-canon/theorems/THM-4695-first-order-universality-of-the-quantum-transfer-for-virasoro-minimal-models.md).
- `T_1` is `C[x][d]`-linear mod `J`.
- The singular vector's length-`(N-1)` part is `α = -2N(N-1)/3 - Nc/12`, `β = 11N(N-1)/48 + 5Nc/96`.
- So `T_1` is affine in `c` and vanishes at the classically free `c_(2,2N+1)` (van Ekeren-Heluani, Thm 15.13). Hence
  `T_1 = (c - c_(2,2N+1)) T*_N` mod `J`, and `E_2` is universal for `p >= 3`.
- **Kac reading** (audit). `c - c_(2,2N+1) = (3/(2N+1)) Kac_3(c, 2N)`, the level-3 Kac factor of `M(c, 2N)`.
- **Precedents.** AvEH's Ising lift is a first-order transfer (N = 3), and Feigin 2011 has the `p = 2` relation.

## 3. Transfer orders: universality below the BGG thresholds (FINITE-EXACT; HYP-9315)

The all-orders engine (`vir_transfer_orders_20261018.py`) assigns every vector of `ker π` a transfer order.
- 16 models are computed. The audit reproduced 10 of them bidegree-for-bidegree, with a separate engine and two
  primes. Totals are in HYP-9315.
- The relative level is `m = L + 2i`, for Li degree `L` at length `N + i`.

**Forced totals.** By the Rocha-Caridi/BGG expansions,
`dim ker π(2N+m) = p(m-3) - p(m-2p+1) - p(m-2q+1)` for `m < min(2(p+q), 4N+1)`. So equal totals carry no information
below the thresholds `2p - 1` and `2q - 1`.

**Refined universality** (post hoc, 0 disagreements in all data; `evaluate_corrected_windows.out`):
- **(U)** All models agree, by length and transfer order, for `m < min(2p-1, 2p'-1)`.
- **(Q)** For fixed `p` they agree for `m < 2 min(q,q') - 1`.
- **(S)** The order-`k` counts with `k <= p_min - 2` agree for `m < 2 q_min - 1`.
- The windows are sharp: every same-`p` pair first disagrees exactly at `m = 2 q_min - 1`. For example (5,6)/(5,7)
  split at `(L,i) = (9,1)`, `{3:2}` vs `{3:3}`.
- **Sharper conjecture** (audit). The refined data at level `m` depend only on `{2p-1, 2q-1} ∩ [0,m]`.

**Pre-registrations.**
1. Vir_(6,7) and Vir_(5,7) generator orders (commit 94fbe3b30b): PASSED. Only k = 3, 4 were informative (audit N4).
2. q-stability and p-saturation on `W(N) = {L <= N+i-1}` (commit dcbb7ffa8a): P4 and P6 are REFUTED as stated.
   - The window contained the forced break at `m = 2q-1`.
   - S44's (5,6) data stopped at `m = 10`, so "all passed" was an artifact of the range.
   - The audit's (5,6) run to weight 32 refutes them, and S44 reproduced it (`orders_56_D32.out`). See MISTAKE-640.
   - P5, P7 and P8 passed.

**Generators and orders.** Compare kind-pasteur's table, section 3.2.
- **Degree `N`.** The generator at Li degree `2k+1` has transfer order exactly `k` in every model (audit,
  invariantly). Only "order `>= k`" for `k >= 3` is non-forced.
- **Why the list stops at `2p-3`** (audit N5). At `m = 2p-1` the Verma map acquires its BGG relation, whose top symbol
  is the Koszul cycle `x_0^(p-2) s_1`. That cycle never transfers.
- **Degree `N + i`.** The order is `<= p - 2 + i`. It is attained for `(p,i) = (3,1), (4,1), (5,1), (3,2)`, and open for
  `(6,1)`.

## 4. The first page (KNOWN input; HYP-9315(a), (d))

- `H_1` of the derived arc space of `Spec C[x]/(x^N)` is generated over `C[x][d]` by the Euler syzygy
  `s_1 = x_0 e_1 - N x_1 e_0`, modulo boundaries.
  - **This is known, and the session missed it (MISTAKE-640).** N = 2: Kanade (thesis Thm 4.3.1; JPAA 2018),
    Bai-Gorsky-Kivinen (2020, Thm 2.2), Gorsky-Oblomkov-Rasmussen (2013, Cor. 3.20). All N: a consequence of
    Feigin 2011 Prop. 3.6, which omits one step. Gorsky-Lewark Conj. 4 is the truncated form.
  - S44 checked N = 2, 3, 4 through weight 24, and the audit N = 2..7. These are now controls.
- `T_1(B_1) ⊆ J` (HYP-9315(d)) is proved, from the audit's sketch.
- Together with THM-4695(a): `Im T_1 = C[x][d] b_N` mod `J`. Observed directly for Ising, (3,5), (4,5) through weight
  22.
- For Ising this is **consistent** with Andrews-van Ekeren-Heluani's theorem (one extra generator, at weight 9) through
  weight 30. The vanishing of all higher transfers for Ising is not proved here.

## 5. Derived arc spaces of plane curves: the first homotopy is a Cramer syzygy (HYP-9316)

**The cusp puzzle.** kind-pasteur noticed that the cusp's arc-space Hilbert series agrees with the product
`Π_(n>=2)(1-t^n)^(-1) / ((1-t^3)(1-t^4)(1-t^5))` through weight 33 and fails at 34.
- The product is the Euler characteristic of the derived arc space, i.e. of the Koszul complex of the jets
  (DMOW 2026, Thm 4.2).
- The first failure is the first `H_1` class. The audit made this rigorous with graded Nakayama.
- On the lct contact locus, six jets lie in a five-generator ideal. Cramer's rule gives a syzygy of weight
  `51 - 17 = 34`.

**The general rule.** The first homotopy weight is the minimal Cramer weight over proper ideals generated by
positive-weight monomials in the jet variables.
- It matches D4 (10), A3 (16), `x^4 + y^4` (11), and the pre-registered E6 (29) and A4 (31).
- It also matches seven further curves (A5, E7, D5, `x^5+y^5`, E8, A6, D6), out of sample in the audit.
- **Correct log-canonical statement** (audit m11). The jets of a plane curve `C` have the expected dimension for all
  `m` iff the pair `(A^2, C)` is log canonical (Mustață 2002), iff `C` has at worst nodes (Mustață 2001, Rem. 4.3),
  iff the derived arc space is classical (DMOW 2026, Cor. 4.9). The lci theorems assume normality.
- The Cramer mechanism and the first-deviation weights were not found in the literature.

## 6. mac-mini's unit-distance graphs: Kloosterman spectrum and Hoffman bound (VERIFIED; classical ingredients)

For `G_q = Cay(F_(q^2), μ_(q+1))` (`gq_kloosterman_hoffman_20261018.py`):
- **(i) Spectrum.** `λ(a) = Σ_(u ∈ μ_(q+1)) ψ(Tr(au)) = -Kl_q(N(a))`. This is already PROVED in the repo as
  THM-4628(b), (c), with Medrano-Myers-Stark-Terras 1996. S44 checked 5 characters per prime for `q <= 60`; the audit
  checked every nonzero character for `q <= 31` (maximum error 7.5e-15). So `λ_min = -max_c Kl_q(c)`.
- **(ii) Hoffman bound.** It is exactly `1 + (q+1)/max_c Kl_q(c)`. This reproduces 4.30825 at `q = 29` and 4.0000439
  at `q = 31`.
  - At `q = 31` the largest Kloosterman sum falls short of `(q+1)/3` by `1.56e-4`. That is a near-coincidence.
  - The next small margins are at `q = 23` (4.0148) and `q = 59` (5.0221).
- **(iii) Corrected** (audit M2).
  - Arc-transitivity and Lovász 1979, Thm 9, give `θ(G_q) = ` Hoffman, so the Lovász θ bound cannot improve on
    Hoffman.
  - Schrijver's `θ'`, which is what the Delsarte LP computes, is not covered by that theorem. The LP equals Hoffman for
    `q = 13..37`, but that is a computation (audit), not a consequence.
  - Stronger SDP bounds do improve it: the three-point bound gives `α(G_29) <= 163.25`, so `κ(29) >= 6` (decalion89;
    mac-mini section 2.2; MISTAKE-634).
  - `κ(29) >= 7` would need `α(G_29) <= 140` or an UNSAT certificate.

## 7. Concept board and open problems

**Concept board.**
- derived arc space (Koszul complex of jets);
- the filtered Verma/BGG map and its spectral sequence;
- singular-vector top lengths;
- universality in `c`;
- BGG thresholds `2p - 1`, `2q - 1`;
- the c-independent Witt action on symbols;
- Cayley semi-invariants;
- Cramer syzygies of contact loci;
- the Kloosterman spectrum.

**Meta-pattern.** Several "local anomalies" seen by other agents were explained by a first obstruction:
- the first syzygy (Cramer, Euler);
- the first differential (`T_1`);
- the first BGG relation (the stop at `2p - 3`);
- the first eigenvalue (the extreme Kloosterman sum).

**Open problems.**
- Prove the order-`k` / Li-degree-`(2k+1)` correspondence, i.e. "order `>= k`" for `k >= 3`.
- Prove the refined universality (U)/(Q)/(S) and the audit's sharper form.
- Fill Feigin's omitted step.
- Prove THM-4710's closed forms for all `N`, and explain the factor `4N^2 + 4N - 43` in `G_9`.
- Extend to admissible `sl_2` (HYP-9293), `W_3` (HYP-9297) and `N = 2` (HYP-9295).
- Find a closed form for the Cramer weight `w_1` of curves.

## 8. Witt rigidity: the degree-N generators depend only on N (THM-4710)

**Lemma (c-independent symbols).** Each positive mode `L_n` acts on top-length symbols by the derivation
`D_n x_j = (j+n+2) x_(j-n)`. This involves no `c`, because the central term and every commutator producing `L_(-1)`
lower the length.
- Hence `grI_N` is stable under `d` and all `D_n`.
- `D_1, D_0, d` form an `sl_2`, so `grI_N` is a sum of lowest-weight Vermas. Its generators are lowest-weight vectors:
  Cayley semi-invariants, with `A_j = Σ C(j+3,i)(-1/4)^i x_0^(j-1-i) x_1^i x_(j-i)`.

**Theorem** (finite range `3 <= N <= 48`, Li degrees 3, 5, 7, 9).
- Given the generator positions, at each odd position there is exactly one lowest-weight vector `G_L(N)` whose
  `D_m`-images lie in the module generated by the earlier generators.
- So every minimal model's degree-`N` generator equals `G_L(N)`, independent of `c`.

**Checks.**
- 18 generator comparisons in 9 models, all equal.
- Models with the same `N` share their generators: `N = 12`, (4,9) and (5,7) both give `G_5 = A_5 - (59/49) A_3A_2`;
  `N = 15`, (4,11) and (6,7).

**Closed forms** (exact for every computed `N`, conjectural in general):
- `G_5 = A_5 - ((5N-1)/(4N+1)) A_3A_2`.
- `G_7` has `A_5A_2`-coefficient `-(28N-33)/(6(4N+1))`.
- The `A_(L-2)A_2`-coefficients tend to `-L/(L-1)` as `N -> ∞`.

**Reading.** The degree-`N` singular supports of all models with a given `N` are truncations of one universal chain:
- `p` decides how far it goes (BGG);
- the Witt action decides what it is.
