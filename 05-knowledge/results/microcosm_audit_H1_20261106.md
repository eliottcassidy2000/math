# Audit H1 of the microcosm packet: THM-4900, HYP-9490 and the gate/doubling note

2026-11-06. Independent audit H1 of the checkpoint `ee4f5cfe49` (branch `mac-mini/session-20261106-microcosm`,
session mac-mini-2026-11-06-microcosm). The audit tried to break the claims. The packet itself was not edited.

**Overall verdict: PASS WITH CORRECTIONS.**
- Every PROVED statement of THM-4900 (1)–(6) and HYP-9490 Lemmas A–C re-derives on paper and checks exactly in
  independent code. All FINITE-EXACT tables reproduce to the printed digits.
- One genuine mathematical error in a (conditionally) proved statement: the HYP-9490 **Consequence** is stated with
  a multiplicative error term, `1 + 2Δc3·ρ·(1 + o(1))`. That form is false whenever c3(T) = c3(T′) but the cores
  differ in a higher sw_j. The correct form is additive, `1 + 2Δc3·ρ + o(ρ)`. The error is repeated in the HYP title,
  note §2 and the INDEX line. Small factual misstatements elsewhere are listed in §6.
- Two significant priority misses:
  - φ(T) = F_T(1) is THM-1950's total inverse-response s(T) (klein-2026-07-21-S400). Its bounds, the regularity
    characterisation, the velocity-addition law and the strong minima for n ≤ 8 are already there.
  - sw_1 = 2c3 at every order (indeed the whole size-weighted odd-cycle form of P_T at every n), the gate's
    cover operator, and the one-vertex tower counts are already in Session159/Session150 (2026-11-02/03).
- Smaller credit misses: Session212's class histograms (24, 228); the 2026-10-06 F21 notes and THM-4556's remark;
  THM-4644(3); HYP-9162/THM-4557; Session35 §6.
- Three small errors in note §3: a missing p ≡ 3 (mod 4) hypothesis, a garbled square-sum sentence, and "linearises
  the join".

## 1. Method

All checks use fresh code in `04-computation/experiments/microcosm_gate_resolvent_20261106/audit_H1/`. No packet
module is imported. `sh run_audit.sh` regenerates every `.out` in about 25 s.

| Script | Checks |
|---|---|
| `audit_tournament_lib.py` | Helpers. Explicit Γ, E, D (dominance [[M, M+I], [M−I, −M]]) and ⇒, with random relabelling of the shell vertices. Exact rational solves. Subset-DP path counts. Cover polynomial by q-fold subset convolution of h (the definition). |
| `audit_marked_pair.py` | THM-4900 (1)–(5), exactly. N is computed by the determinant lemma, N = det(zI − M + J) − det(zI − M), not by the adjugate or the packet's Schur route. |
| `audit_towers_numeric.py` | THM-4900 (2)–(3) numerically: spectra for k ≤ 14, gate-tower limits and contraction rates, D- and E-towers, F_T(z) > 0. |
| `hp_count.c` | Brute-force H and HC (rooted cycle DP), up to order 21. |
| `audit_tower_ocf.py` | THM-4900 (6) and HYP-9490 Lemmas A–B. P_T by subset convolution is compared with the odd-cycle set-partition (Grinberg–Stanley) form at odd and even n. The tower formula is checked against brute-force H of explicit towers. U_(d,r) is computed three ways: L on coefficient lists, the exponential walk W_r (the transpose route), and Taylor series of M^r e^y. |
| `audit_hyp9490_tables.py` | All HYP-9490 tables: N up to 1201, all 528 cores of orders 4–7, the rate table. |
| `audit_consequence_kappa0.py` | The Δc3 = 0 counterexample to the multiplicative Consequence. |
| `audit_data.py` | φ minima for n = 5–9 (gentourng -c). HC = 7 and 21 class counts and gates for n ≤ 8. Gated spectra for core orders ≤ 7. Note §4: face formula and its determination by L P_T. |
| `audit_mod_p_groups.py` | Note §3 group theory. |
| `audit_rate_and_k1_shape.py` | d = 3 against d = 4 normalised rate. A saddle-point normalisation. The shape of Γ^k(K1). |
| `audit_priority_hits.sh` | Prints file:line of every earlier in-repo appearance cited below. |

The packet's own `run_all.sh` was rerun in a scratch copy. All 13 committed `.out` files reproduce byte-identically.

## 2. Verdicts

| Item | Verdict | Main reason |
|---|---|---|
| THM-4900 (1) gate transfer, eigen-coordinates, P_Γ = P_E | **PASS** | Conventions match (Γ: s → T → t, t → s, which is Session210's a → b → T → a with a = t, b = s; E: s → t). Re-derived; exact on all 75 labelled cores n ≤ 4 and all 68 classes n = 5, 6. |
| THM-4900 (2) nested gates, limit, one-vertex spectrum | **PASS** (cosmetic) | Closed form exact on 60 (core, k) pairs. One-vertex formula exact for k ≤ 8. Spectrum to 1.1e-14 for k ≤ 14. The R(T) ≠ 0 argument is sound. |
| THM-4900 (3) z = 1, φ, bounds, φ(DT) = φ(T), towers | **PASS WITH CORRECTIONS** | Mathematics correct. φ is THM-1950's s(T), and the bounds with equality iff regular are KNOWN there (credit). The E-tower proof omits C_T(z) ≠ 0 (cosmetic). |
| THM-4900 (4) series law, golden map, Fibonacci orbit | **PASS WITH CORRECTIONS** | Correct. At z = 1 the series law is THM-1950's velocity addition (credit). "144/34" is unreduced, and the title's "convergents to phi^3" is imprecise (cosmetic). |
| THM-4900 (5) D and Γ | **PASS** | Exact for all classes n ≤ 5. THM-447(4) is cited correctly. Session211 is credited. |
| THM-4900 (6) sw_1 = 2c3, tower formula | **PASS WITH CORRECTIONS** | The proof is correct, line by line. But sw_1 = 2c3 at every n, and the identification of every sw_j with odd-cycle data at even n, are KNOWN (Session159 cycle_profiles (1)). As written, the file has no justification for "the core enters only through its odd-cycle statistics" at even n beyond sw_1 (THM-078 is odd n only). The d = 0 weight needs the convention Γ(∅) = TT_2. |
| THM-4900 front matter ("New here", depends_on) | **PASS WITH CORRECTIONS** | Credit: THM-1950, Session150/159, Session35 §6. |
| HYP-9490 Lemma A | **PASS** (credit, cosmetic) | Re-derived; three-way exact check for r ≤ 20, d ≤ 12. Session159 (2) is an earlier differential-operator form of the same operator. "The tower is a diffusion" is loose wording. |
| HYP-9490 Lemma B | **PASS** | Proofs correct (B.1 only needs, and only proves, monotonicity within a parity class). 0 violations for e ≤ 40, r ≤ 200. |
| HYP-9490 Corollary C | **PASS** | Correct, including d = 2: U_(0,r) = U_(2,r−1) exactly. |
| HYP-9490 Consequence | **PASS WITH CORRECTIONS (error)** | The multiplicative form is false at Δc3 = 0 (explicit order-5 example). The additive form follows from (a) and B.1, but the proof is not written. ρ_N is undefined. |
| HYP-9490 evidence tables | **PASS** (cosmetic) | Every entry reproduced by an independent route. "d = 4 agrees to within 1%" is false for N ≤ 44. Two values are rounded down. |
| HYP-9490 rate (c) | **PASS** as an OPEN conjecture | Honestly flagged. But the INDEX states "rho ~ (log N)^2/N^2" as a fact (overclaim). |
| Note §0, §1, §2 | **PASS WITH CORRECTIONS** | Credit lists incomplete. §2 repeats the Consequence error. The S84 sentence is now outdated. |
| Note §3 ({3,10}, {15,25}, {7,21}, Fano) | **PASS WITH CORRECTIONS** | Prime list lacks p ≡ 3 (mod 4). Square-sum sentence wrong. "linearises the join" is false. The class counts and the F21 reading are KNOWN in-repo. The typing should be KNOWN/ANALOGY. |
| Note §4 face formula | **PASS** | Re-derived from junction counts. Matches brute force on 17 random cores n ≤ 6. A_Γ(T) is recovered from L P_T alone, which confirms "KNOWN in content". |
| Note §5 data | **PASS WITH CORRECTIONS** | All numbers reproduce. The n ≤ 8 φ-minima are KNOWN (klein S400, THM-4051). "Second value" holds only for n ≥ 3. |
| Note §6 ledger | **PASS WITH CORRECTIONS** | Rows need the credits and typing below. |
| INDEX line | **PASS WITH CORRECTIONS** | Consequence form, the rate overclaim, credits. |

## 3. Findings in detail

### THM-4900 (1)–(5)

- **Gate transfer.** Ordering (s, core, t) gives zI − M_Γ = [[z, −uᵀ, 1], [u, B, −u], [−1, uᵀ, z]]. This is the
  stated orientation. For the gate, the Schur complement of B = zI − M is [[z+F, 1−F], [−1−F, z+F]]. For E it is
  [[z+F, −1−F], [1−F, z+F]]. Both have determinant z² + 1 + 2zF, so P_Γ = P_E.
- **Response.** The boundary equations za − b = 2 and zb + a = 2(1 − b)F give the stated F_Γ and the matrix G(z).
  det G = (z² + 3)² and tr G = 2z² + 10, so the eigenvalues are z² + 5 ± 2s = (s ± 1)². The eigenvectors (z, −2 ± s)
  check symbolically modulo s² = z² + 4.
- **Limit.** F_T(z) = z·uᵀ(z²I + MᵀM)⁻¹u > 0 for real z > 0, and P_T(z) > 0 there.
  - So R(T) ≠ 0. If Ψ−(T) = 0, F is already at the limit.
  - The observed per-gate contraction equals ((s−1)/(s+1))² to within a relative 4.4e-4 at k = 8–10, for
    z ∈ {1/2, 1, 2, 3, 4}.
  - The metallic-mean algebra is correct.
- **One vertex.** The bracket [(s+1)^(2k−1) + (s−1)^(2k−1)]/(2s) is a polynomial in s² of degree k − 1.
  - The root analysis (s = −i cot(θ/2), exclusion of s = 0, the pairing j ↔ 2k−2−j) is correct.
  - Γ^k(K1) is TT_(2k+1) with exactly the k concentric chords (s_i, t_i) reversed (checked for k ≤ 8).
- **φ.** For all 7412 classes of orders ≤ 8:
  - uᵀ(I − M)⁻¹u = uᵀ(I + MᵀM)⁻¹u;
  - 0 < φ ≤ n, with φ = n exactly for the 6 regular classes;
  - φ equals THM-1950's 𝟙ᵀ(I + K)⁻¹𝟙 (532 classes of orders ≤ 7; trivially, since uᵀXu = uᵀXᵀu).
- **Doubling.** Exact for 32 random cores of orders ≤ 8: φ(DT) = φ(T) and P_D(1) = 2ⁿP(1)².
  - Session35's D formulas were re-derived from D² = diag(2M² − I, 2M² − I) and checked exactly for n ≤ 5.
  - The D-tower telescoping checks for k ≤ 4. The E-tower tends to 1.
  - The double used, [[M, M+I], [M−I, −M]], is THM-447's. It is also the one in Session211 and HYP-9162.
- **Series law.** P_(A⇒B) = P_A P_B + N_A N_B and N_(A⇒B) = N_A P_B + P_A N_B, exactly, on 25 random pairs. This is
  THM-3322 (18) after z ↦ (1+z)/z.
- **Golden map and Fibonacci orbit.**
  - φ(ΓT) = (5φ+1)/(φ+1), φ(ET) = 1, K1 absorbing, φ(C3 ⇒ C3) = 3/5.
  - φ(Γ^k K1) for k = 0..10 is 1, 3, 4, 21/5, 55/13, 72/17, 377/89, 987/233, 1292/305, 6765/1597, 17711/4181
    = F_(2k+2)/F_(2k−1).
  - The multiplier at 2 + √5 is g′ = 4/(3+√5)² = φ_gold⁻⁴.
- **Statement 5.** Checked for all classes n ≤ 5:
  - φ(DΓT) = φ(ΓDT) = g(φ(T)) ≠ φ(Γ²DT);
  - D(ΓT) has no vertex of out-degree 1;
  - P_(DΓT) = P_(DET) (for n ≤ 4).

### THM-4900 (6)

The coefficient computation is right:
- e_2 of (f, f−1, …, f−n+1);
- Σ_π f = (n−1)n!/2;
- each transitive triple carries one directed 2-path and each cyclic triple three, giving C(n,3) + 2c3 in total;
- each 2-path carries weight 2(n−2)(n−3)!, and the non-adjacent position pairs contribute amounts independent of T;
- a(n) = 0 from P_(TT_n) = qⁿ.

Checks:
- sw_1 = 2c3 on all classes n ≤ 7 and on random tournaments n = 8, 9, 10.
- The full set-partition form P_T(q) = Σ over partitions into singletons (weight q) and odd blocks (weight
  2·HC(T[B])·q) equals P_T on 652 tournaments, including **even** n ≤ 8.
- The tower formula matches brute-force H on 280 explicit towers, and U_(d,r) = H(Γ^r TT_d) on 61 values (total order
  ≤ 19). The d = 0 values use Γ(∅) = TT_2.

**Priority.** Session159 `cycle_profiles/NOTE.md` (2026-11-03) proves, for every n via OCF,

    F_T(r) = Σ_K 2^|K| r^(n − Σ(|C|−1)),

and says explicitly that "the r^(n−2) coefficient is twice the number of directed triangles; the r^(n−4) coefficient
is twice the number of directed 5-cycles plus four times the number of disjoint pairs of directed triangles". So
"sw_1 = 2c3 for every order" is KNOWN, and the HYP's sw_2 test formula is Session159's. The same note gives the gate
operator as (2). Session159 `profile_operator/NOTE.md` gives the operator as (1), a positive second-order
differential form (2), and the iterated family G^m(K1) with h = 1, 3, 15, 123, 1479, 23955, 496479. Session150
`transfer/TRANSFER.md` (2026-11-02) has the triangle-wrapper transfer. Session210, which the packet credits for the
operator, came later (2026-11-06 17:52).

### THM-4900 (3)–(4): φ is THM-1950's s

THM-1950 (klein-2026-07-21-S400; file first committed 2026-07-21) defines the total inverse-response
s(T) = 𝟙ᵀ(I + K)⁻¹𝟙 with K = A − Aᵀ. Since uᵀ(I + M)⁻¹u = uᵀ(I − M)⁻¹u, this is φ. THM-1950 proves:
- s = ‖(I+K)⁻¹𝟙‖² ≥ 0 and s ≤ n;
- s = n iff T is regular;
- s(T₁ ⇒ T₂) = (s₁ + s₂)/(1 + s₁s₂), the "relativistic velocity-addition" law;
- s(K1) = 1.

Related sources:
- klein's reflection `07-reflections/the-skew-determinant-composes-by-velocity-addition-klein-S400.md` lists the
  strong minima of s: 3, 3, 2, 1, 0.667, 0.556 for n = 3..8.
- THM-4051 gives the complete strong order-7 range [2/3, 7].
- THM-3729 identifies det(I+K)·s, which is N_T(1), with an odd rooted-Pfaffian square sum.

What is new about φ, as far as the audit found:
- its reading as the value at the self-similar point z = 1 of Session35's marked response;
- φ(DT) = φ(T);
- the golden gate map and the Fibonacci orbit;
- the order-9 minimum and the minimiser counts and scores.

### HYP-9490

- **Lemma A.** The Leibniz computation with sinh^(i)(0) = [i odd] matches L q^e = q^(e+2) + Σ_k 2C(e,2k+1)q^(e−2k)
  term by term. ⟨P, e^y⟩ = P(1), and M = e^(−2cosh y)∂e^(2cosh y)∂. Three independent computations of U_(d,r) agree
  for r ≤ 20, d ≤ 12. The r = 1 identity H(ΓT) = P(1) + P(2) − P(0) agrees with brute force.
- **Lemma B.** Correct as proved.
- **Corollary C.** Correct as proved.
- **Consequence.** Under (a), the argument is short and should be written:
  - H(Γ^r T)/U_(n,r) = 1 + 2c3(T)ρ_(n,r) + Σ_(j≥2) sw_j(T)·U_(n−2j,r)/U_(n,r).
  - For j ≥ 2, Lemma B.1 gives U_(n−2j,r) ≤ U_(n−4,r).
  - U_(n−4,r)/U_(n,r) = ρ_(n−2,r)·ρ_(n,r) = o(ρ_(n,r)), by (a) at d = n − 2.

  This yields **1 + 2Δc3·ρ + o(ρ)**, not 1 + 2Δc3·ρ·(1 + o(1)). The two differ exactly when Δc3 = 0.
- **Counterexample to the multiplicative form** (`audit_consequence_kappa0.out`). The three order-5 classes with
  c3 = 4 have P = q⁵ + 8q³ + cq with c = 2, 4, 6 (H = 11, 13, 15). The multiplicative form predicts a ratio of
  exactly 1. The actual H(Γ^r T)/H(Γ^r T′) − 1 between c = 2 and c = 6 is −2.24e-4 at r = 20 and −7.9e-6 at r = 60.
  It is never 0.
- **Evidence at N ≈ 401.** The κ = 0 classes with dev ≠ 0 number 2 of 3, 11 of 12 and 60 of 61 at n = 5, 6, 7. The
  packet's own residual column measures the additive form, so the evidence supports the corrected statement.
- **Rate (c).**
  - d = 3 against d = 4 normalised values differ by 3.0% at N = 23/24, 1.2% at N = 43/44, then by at most 0.44%.
  - An optional upgrade, typed HEURISTIC: the walk heuristic of the HYP gives the height K by K(log K + 2) = N. The
    normalisation N²ρ/(log K + 2)² then increases monotonically, 0.782, 0.787, 0.802, 0.811, 0.817, 0.822, 0.826 for
    N = 23..243. That is consistent with Θ((log N)²/N²), whereas N²ρ/(log N)² is still decreasing. Neither settles (c).

### Note §3–§5

- **Group theory** (`audit_mod_p_groups.out`).
  - ⟨x+2, 2x⟩ mod 7 has order 21, with multipliers {1, 2, 4} = QR7. It equals Aut(QR7), computed by brute force over
    S7 (|Aut| = 21).
  - The sets x + {1, 2, 4} form the Fano plane. F21 has index 8 in PSL(3,2) (classical). |Aut QR23| = 253 (affine
    automorphisms counted).
  - The primes p < 200 with ord_p(2) = (p−1)/2 are 7, 17, 23, 41, 47, 71, 79, 97, 103, 137, 167, 191, 193, 199. Only
    those with p ≡ 3 (mod 4) (7, 23, 47, 71, 79, 103, 167, 191, 199) give a Paley tournament. So the note's condition
    "every prime p where 2 generates exactly the quadratic residues" is false as stated (17, 41, 97, …). It needs
    p ≡ 3 (mod 4), exactly as in THM-4644(3).
  - The identification F21 = {x ↦ 2^i x + b mod 7} = "doubling plus translation" is already in:
    - THM-4556's remark (typed ANALOGY there);
    - `mersenne_switch_parity_f21_compression_20261006.md` (Prop. 3 and its abstract);
    - `seven_twentyone_mersenne_openai_math_20261006.md` (§0, §2).

    HYP-9162/THM-4557 gives the genuinely tournament-level link between ×2 and F21: the skew-Hadamard doubling tower
    over C3 passes through QR7 and keeps Aut = F21 at every level.
- **Square sums.** Session211 has p(Q15) = p(Q16) = p(Q17) = 1, p(Q23) = 3, p(Q24) = 0, p(Q25) = 10. The note's
  "the first two square-sum path counts after the onsets 23 and 25 are exactly 3 and 10" is wrong. The counts *at*
  the onsets 23 and 25 are 3 and 10, and they are not the first counts.
- **"linearises the join".** In q = (H − 1)/2 the join is q + q′ + 2qq′. That is not linear (1 + 2q is
  multiplicative). The parallel with y = (φ − 1)/2 is an ANALOGY, not STRUCTURAL.
- **HC classes** (`audit_data.out`).
  - HC = 7: 24 classes at order 7 and 118 at order 8. HC = 21: 228 classes at order 8.
  - None of them is gated, and neither value occurs at smaller orders.
  - These counts are already in Session212 `cycle-levels/n7_hist.tsv` (7 → 24) and `n8_hist.tsv` (7 → 118,
    21 → 228). That data release is a complete isomorphism-class census through order 8.
- **Gated spectra.** Values for core orders 1–7 reproduce: the minimum is 2ⁿ + 1, the second value 3 + 3·2^(n−1)
  for n ≥ 3, and the values ≤ 40 are {3, 5, 9, 15, 17, 27, 33, 37}.
- **φ minima.** For n = 5..9 the minima are 2, 1, 2/3, 5/9, 1/2, attained by 1, 1, 1, 1, 2 strong classes. The
  minimiser scores are as stated. At n = 9 the two minimisers have scores [2,2,2,2,3,5,6,7,7] and
  [1,1,2,3,5,6,6,6,6]. The maximum is n at odd n, and 5 and 7 at n = 6 and 8.
- **§4 face formula.** The junction exponents ([A≠∅] + [B=∅] + [C≠∅] for s…t, and 2[B≠∅] for t…s) re-derive from
  the arc directions. The formula matches brute force. A_Γ(T) is reproduced from L P_T by inverting Session210 (7).

### Other priority searches (negative)

- No earlier appearance was found of: the gate's 2 × 2 transfer on (N, P); P_Γ = P_E; the metallic limit; the
  one-vertex cot² spectrum; φ(DT) = φ(T); the golden map φ ↦ (5φ+1)/(φ+1); the adjoint diffusion
  M = ∂² + 2 sinh(y)∂; the tower odd-cycle formula as an explicit identity; HYP-9490's first-order forgetting law.
- "gate tower" in THM-1286 and the July INDEX is an unrelated LRC notion.
- THM-4900 and HYP-9490 IDs are unique on the worktree and on `origin/main`.
- S84 has since pushed THM-4890 (one-apex collision balance) and THM-4891 (face-local coherence relaxations). Neither
  overlaps.
- A web search did not find 1, 3, 15, 123, 1479, 23955 in OEIS. This is not an OEIS API check, so the "not in OEIS"
  claim remains the authors'.

## 4. Required corrections

Severity: **error**, **overclaim**, **typing**, **credit**, **cosmetic**. "Old" is the exact current text.

### HYP-9490 (`05-knowledge/hypotheses/HYP-9490-gate-towers-forget-the-core-up-to-cyclic-triangles.md`)

1. **error** — title.
   - Old: `H(Gamma^r T)/H(Gamma^r T') = 1 + 2(c3(T) - c3(T')) rho_(n,r) (1 + o(1)) as r -> infinity`
   - New: `H(Gamma^r T)/H(Gamma^r T') = 1 + 2(c3(T) - c3(T')) rho_(n,r) + o(rho_(n,r)) as r -> infinity`
2. **error** — Consequence, first display.
   - Old: `    H(Γ^r T)/H(Γ^r T′) = 1 + 2(c3(T) − c3(T′))·ρ_(n,r)·(1 + o(1)).`
   - New: `    H(Γ^r T)/H(Γ^r T′) = 1 + 2(c3(T) − c3(T′))·ρ_(n,r) + o(ρ_(n,r)).`
   - Insert after the display:
     > *Proof (from (a)).* By THM-4900(6), H(Γ^r T)/U_(n,r) = 1 + 2c3(T)ρ_(n,r) + Σ_(j≥2) sw_j(T)U_(n−2j,r)/U_(n,r).
     > For j ≥ 2, Lemma B.1 gives U_(n−2j,r) ≤ U_(n−4,r), and U_(n−4,r)/U_(n,r) = ρ_(n−2,r)ρ_(n,r) = o(ρ_(n,r)) by
     > (a) at d = n − 2 (n ≥ 4; for n ≤ 3 there is no such term). Since ρ_(n,r) → 0 by (a) at d = n, dividing the
     > two expansions gives the display. ∎
     >
     > The error term is o(ρ), not 2Δc3·ρ·o(1). When c3(T) = c3(T′) the ratio is 1 + O(ρ_(n−2,r)ρ_(n,r)), and in
     > general it is not 1. For example, the order-5 classes with P = q⁵ + 8q³ + 2q and q⁵ + 8q³ + 6q both have
     > c3 = 4, and their towers differ by −2.2·10⁻⁴ at r = 20.
3. **error / cosmetic** — Consequence, second display.
   - Old:
     ```
     - With (b), which follows from (a) by Corollary C, the same holds against the one-vertex tower of the same total
       order:

           H(Γ^r T)/H(Γ^((N−1)/2) K1) = 1 + 2(c3(T) − ((n−1)/2)²)·ρ_N·(1 + o(1))    (n odd),
     ```
   - New:
     ```
     - Taking T′ = Γ^((n−1)/2)(K1), which has order n, compares T with the one-vertex tower of the same total order:

           H(Γ^r T)/H(Γ^((N−1)/2) K1) = 1 + 2(c3(T) − ((n−1)/2)²)·ρ_(n,r) + o(ρ_(n,r))    (n odd),

       and by (b) (Corollary C) ρ_(n,r) may be replaced by ρ_(d,r′) for any fixed d with d + 2r′ = N.
     ```
     (ρ_N is never defined, and (b) is not needed for the display itself.)
4. **cosmetic** — rate table bullet.
   - Old: `d = 4 agrees to within 1%.`
   - New: `d = 4 agrees to within 3.0% at N = 24, 1.2% at N = 44, and 0.5% from N = 84 on.`
5. **cosmetic** — next bullet.
   - Old: `If (a) holds, this is the same quantity` / `asymptotically (Corollary C).`
   - New: `It equals ρ_(3,r)/(1 + 2ρ_(3,r)) exactly (P_(C3) = q³ + 2q), so it is ρ_(3,r) to relative accuracy 2ρ < 4·10⁻⁴; no appeal to (a) is needed.`
   - In the table, `0.870` → `0.871` and `0.863` → `0.864` (0.87057 and 0.86357 rounded).
6. **cosmetic / credit** — Lemma A label.
   - Old: `**Lemma A (the tower is a diffusion; PROVED).**`
   - New: `**Lemma A (the tower is the r-th power of a diffusion generator; PROVED).**`
   - Add a sentence: `Session159 profile_operator/NOTE.md (2) had already written the gate operator as the positive second-order operator P_G = (z+z²)P″ + (3z+2z²)P′ + (z+z²)P on unordered-cover generating functions; the adjoint form below is different.`
7. **credit** — `related:`.
   - Old: `related: [THM-4900, Session210 graded-spectrum-gates (operator, Theorem 4/5), THM-078 (odd-cycle form), ...`
   - New: add `Session159 cycle_profiles (1) and profile_operator (1)-(2) (operator, all-order odd-cycle form of P_T, sw_2 = 2c5 + 4 pairs), Session150 transfer`.

### THM-4900 (`01-canon/theorems/THM-4900-gate-and-doubling-in-the-marked-resolvent-golden-towers-and-tower-ocf.md`)

8. **credit** — "New here", third item.
   - Old: `  - the self-similar point z = 1, the invariant phi and its golden gate map;`
   - New: `  - the self-similar point z = 1: phi = F_T(1) is THM-1950's total inverse-response s(T) = 1^T(I+K)^(-1)1 (its bounds, s = n iff regular, the velocity-addition law under =>, and the strong minima for n <= 8 are KNOWN: THM-1950, klein-2026-07-21-S400; order 7: THM-4051); new here are its reading as the self-similar point of D, phi(D T) = phi(T), and the golden gate map;`
9. **credit** — "New here", fourth item.
   - Old: `  - the doubling and padding towers;`
   - New: `  - the doubling and padding tower limits (one-line corollaries of Session35's F_D formula and its all-depth law C_(E^k) = ((z+1)/(z-1))^(2k) C);`
10. **credit** — "New here", last item.
    - Old: `  - the tower odd-cycle formula and sw_1 = 2 c3 for every order.`
    - New: `  - the tower odd-cycle formula (sw_1 = 2 c3 and the whole size-weighted odd-cycle form of P_T at every order are KNOWN: Session159 cycle_profiles/NOTE.md (1), 2026-11-03; re-proved here).`
11. **credit** — KNOWN ingredients and `depends_on`.
    - Add two bullets:
      - `phi = THM-1950's s (klein-2026-07-21-S400) with THM-4051 and THM-3729 (det(I+K)·s = N_T(1) = odd rooted-Pfaffian square sum)`;
      - `the gate's cover operator first in Session150 transfer/TRANSFER.md and Session159 profile_operator/NOTE.md (1)-(2) (2026-11-02/03), with the all-order cycle-packing form and the one-vertex family values 3, 15, 123, ..., 496479`.
    - In `depends_on:` add `THM-1950` and `session159 cycle_profiles + profile_operator NOTE.md`.
12. **credit** — statement 3.
    - Old: `Put φ(T) := F_T(1) = uᵀ(I − M)⁻¹u = uᵀ(I + MᵀM)⁻¹u.`
    - New: append ` This is THM-1950's total inverse-response s(T) = 𝟙ᵀ(I + K)⁻¹𝟙 with K = M (a scalar equals its transpose).`
    - Old: `- 0 < φ(T) ≤ n, with equality iff T is regular.`
    - New: `- 0 < φ(T) ≤ n, with equality iff T is regular (KNOWN: THM-1950).`
13. **credit** — statement 4.
    - Old: `In particular φ(A ⇒ B) = (φ_A + φ_B)/(1 + φ_A φ_B).`
    - New: `In particular φ(A ⇒ B) = (φ_A + φ_B)/(1 + φ_A φ_B), which is THM-1950's velocity-addition law for s.`
14. **credit / justification** — statement 6.
    - Old: `THM-078 identifies sw_j for odd n as the size-weighted counts of independent odd-cycle sets).`
    - New: `Session159 cycle_profiles/NOTE.md (1) identifies sw_j at every order as Σ 2^|K| over collections K of vertex-disjoint directed odd cycles with Σ_(C∈K)(|C|−1) = 2j; THM-078 is the odd-n case).`
    - Old: `- sw_0 = 1, and sw_1(T) = 2c3(T) for every n.`
    - New: `- sw_0 = 1, and sw_1(T) = 2c3(T) for every n (KNOWN: Session159 (1); an independent elementary proof is below).`
15. **cosmetic** — statement 6.
    - Old: `the tower enters only through universal weights U_(d,r) := H(Γ^r TT_d).`
    - New: `the tower enters only through universal weights U_(d,r) := H(Γ^r TT_d) (for d = 0 with Session210's convention Γ(∅) = TT_2).`
16. **cosmetic** — statement 4.
    - Old: `(1, 3, 4, 21/5, 55/13, 144/34, …)`
    - New: `(1, 3, 4, 21/5, 55/13, 72/17, …)`
17. **cosmetic** — title.
    - Old: `runs through Fibonacci convergents to phi^3`
    - New: `runs through the Fibonacci ratios F_(2k+2)/F_(2k-1) -> phi^3 (in y = (phi-1)/2, the even convergents of the golden ratio)`
    - Among the values 1, 3, 4, 21/5, 55/13, … only 4, 72/17 and 1292/305 (k ≡ 2 mod 3) are convergents of
      2 + √5 = [4; 4, 4, …].
18. **cosmetic** — proof of statement 3, E-tower.
    - Old: `So C → ∞ and`
    - New: `Since F_T(z) > 0 (proof of 2), C_T(z) ≠ 0, so C → ∞ and`
19. **cosmetic** — notation, throughout THM-4900, the note and the INDEX. φ denotes both the invariant and the golden
    ratio (as in "φ³ = 2 + √5" next to "φ(T)"). Suggest THM-1950's name s(T) for the invariant.

### Note (`05-knowledge/results/microcosm_gate_doubling_resolvent_20261106.md`)

20. **error** — §2.
    - Old: `    H(Γ^r T)/H(Γ^r T′) = 1 + 2Δc3·ρ_(n,r)·(1 + o(1)).`
    - New: `    H(Γ^r T)/H(Γ^r T′) = 1 + 2Δc3·ρ_(n,r) + o(ρ_(n,r)).`
21. **credit** — §0 "Older" list. After `THM-1370 (7 and 21 never occur).` add:
    > - THM-1950 (klein-2026-07-21-S400; with THM-4051, THM-3729): the total inverse-response s(T) = 𝟙ᵀ(I+K)⁻¹𝟙,
    >   which is φ, with its bounds, the regularity characterisation and the velocity-addition law; strong minima to
    >   n = 8.
    > - Session150 `transfer` and Session159 `profile_operator` / `cycle_profiles` (2026-11-02/03): the gate
    >   (triangle-wrapper) operator, its positive differential form, the all-order cycle-packing form of P_T
    >   (sw_1 = 2c3, sw_2 = 2c5 + 4·pairs), and the one-vertex family 3, 15, 123, …, 496479.
    > - THM-4556 (remark), the 2026-10-06 notes `mersenne_switch_parity_f21_compression` and
    >   `seven_twentyone_mersenne_openai_math` (F21 = doubling plus translation mod 7), THM-4644(3) (Paley primes),
    >   HYP-9162/THM-4557 (the doubling tower over C3 keeps F21).
    > - Session211 `arithmetic-cover-bridge` (the gate operator beyond tournaments).
22. **credit** — §0 Session212 bullet.
    - Old: `` - `cycle-levels`: complete fixed-order cycle supports through order 9. HC = 7 first occurs at order 7 and HC = 21 at `` / `order 8.`
    - New: append ` Its class histograms already give 24 classes with HC = 7 at order 7, and 118 and 228 classes with HC = 7 and HC = 21 at order 8.`
23. **cosmetic** — §0 S84 bullet.
    - Old: `The reservation has no content yet, so no comparison is possible.`
    - New: `S84 has since pushed THM-4890 (one-apex collision balance) and THM-4891 (face-local coherence relaxations of the permutohedron); neither overlaps this packet, and THM-4891 belongs in §4's "Across months" list.`
24. **credit** — §2.
    - Old: `**The condensation multiplies** (THM-3322's join law, re-proved in THM-4900(4)).`
    - New: `**The condensation multiplies** (THM-3322's join law; at z = 1 THM-1950's velocity-addition law; re-proved in THM-4900(4)).`
25. **credit** — §2.
    - Old: ``- The one-vertex tower's path counts 1, 3, 15, 123, 1479, 23955, 496479, 12725547, 393325815 have c3 = k². They are``
    - New: prefix `Session159 profile_operator/NOTE.md already computed this family through 496479 (with v_3(h) = 1 and its mod-9 law), and cycle_profiles/NOTE.md its first polynomials. `
26. **overclaim / typing** — §3.
    - Old: `- **The same affine normalisation** (STRUCTURAL).`
    - New: `- **The same affine normalisation** (ANALOGY).`
    - Old: `for path counts it linearises` / `the join, for the resolvent it conjugates the gate to F².`
    - New: `for path counts it writes the join as q + q′ + 2qq′ (1 + 2q is multiplicative), for the resolvent it conjugates the gate to the square of y ↦ 1 + 1/y.`
27. **error (wording)** — §3.
    - Old: `The first two square-sum path counts after the onsets 23 and 25 are exactly 3` / `and 10.`
    - New: `The square-sum path counts at the onsets 23 and 25 are exactly p(Q23) = 3 and p(Q25) = 10 (after p(Q15) = p(Q16) = p(Q17) = 1).`
28. **credit** — §3.
    - Old: ``  - This session adds only the class counts at those least orders (FINITE-EXACT, `hc_gate_free.out`, gentourng):``
    - New: ``  - The class counts at those least orders are already in Session212 `cycle-levels` (n7_hist.tsv, n8_hist.tsv); this session reproduces them with gentourng (`hc_gate_free.out`) and checks the (immediate) gate-freeness:``
29. **typing / credit** — §3.
    - Old: `- **Fano and the duality mod 7** (STRUCTURAL; NUMEROLOGY with respect to the gaps).`
    - New: `- **Fano and the duality mod 7** (KNOWN, classical; in-repo: THM-4556's remark and the 2026-10-06 F21 notes; ANALOGY as a reading of the tournament recursions; NUMEROLOGY with respect to the gaps).`
30. **error** — §3.
    - Old: `The same happens for every prime p where 2 generates exactly the quadratic residues:` / `p = 23 (|Aut QR23| = 253), 47, 71, 79, …`
    - New: `The same happens for every prime p ≡ 3 (mod 4) where 2 generates exactly the quadratic residues (THM-4644(3): then the doubling-invariant circulant tournaments on Z_p are exactly QR_p and its converse): p = 23 (|Aut QR23| = 253), 47, 71, 79, 103, …; for p ≡ 1 (mod 4) (17, 41, 97, …) QR_p is a graph, not a tournament.`
    - Optionally add: `HYP-9162/THM-4557: the doubling tower over C3 reaches QR7 and keeps Aut = F21 at every level.`
31. **credit** — §5.
    - Old: ``- **Smallest φ over strong tournaments** (`phi_strong_minimum.out`, nauty gentourng):``
    - New: ``- **Smallest φ over strong tournaments** (`phi_strong_minimum.out`, nauty gentourng; the values for n ≤ 8 are KNOWN from klein's S400 reflection and THM-4051, and n = 9 and the class counts are new):``
32. **cosmetic** — §5.
    - Old: `and the second value is 3 + 3·2^(n−1).`
    - New: `and, for n ≥ 3, the second value is 3 + 3·2^(n−1).`
33. **typing / credit** — §6 ledger.
    - Old: `| z = 1 doubling invariance of φ, φ ∈ (0, n], golden gate map, Fibonacci orbit | PROVED |`
    - New, two rows:
      - `| z = 1 doubling invariance of φ, golden gate map, Fibonacci orbit | PROVED |`
      - `| φ = THM-1950's s; φ ∈ (0, n], = n iff regular; velocity addition at z = 1; strong minima n ≤ 8 | KNOWN (THM-1950, klein S400, THM-4051); re-proved |`
    - Old: `| Relativistic addition over the condensation | PROVED here; KNOWN in substance (THM-3322) |`
    - New: `| Relativistic addition over the condensation | PROVED here; KNOWN in substance (THM-3322 at every z; THM-1950 at z = 1) |`
    - Old: `| Tower odd-cycle formula; sw_1 = 2c3 at every order | PROVED (from Session210 + THM-078) |`
    - New, two rows:
      - `| Tower odd-cycle formula | PROVED (Session210 Theorem 4 + Session159 (1)) |`
      - `| sw_1 = 2c3 (and the whole size-weighted form) at every order | KNOWN (Session159 cycle_profiles (1)); re-proved |`
    - Old: `| Least orders 7, 8 for HC = 7, 21 | KNOWN (THM-4889, Session212); class counts 24, 228 FINITE-EXACT here |`
    - New: `| Least orders 7, 8 for HC = 7, 21, with class counts 24, 228 | KNOWN (THM-4889; Session212 cycle-levels histograms); reproduced here |`
    - Old: `| x+2 / x·2 mod 7 = F21 = Aut(QR7), the cyclic Fano symmetry | STRUCTURAL |`
    - New: `| x+2 / x·2 mod 7 = F21 = Aut(QR7), the cyclic Fano symmetry | KNOWN (classical; THM-4556 remark, 2026-10-06 F21 notes, THM-4644(3)); ANALOGY as a reading |`

### INDEX line (`05-knowledge/hypotheses/INDEX.md`, line 2)

34. **error / overclaim**
    - Old: `H(Gamma^r T)/H(Gamma^r T') = 1 + 2 Delta c3 rho (1 + o(1)), rho ~ (log N)^2/N^2;`
    - New: `H(Gamma^r T)/H(Gamma^r T') = 1 + 2 Delta c3 rho + o(rho), rho conjecturally Theta((log N)^2/N^2) (normalised data still drifting);`
35. **credit**
    - Old: `is an exact doubling invariant, series composition is relativistic addition (THM-3322 in substance; K1 absorbing)`
    - New: `is an exact doubling invariant (phi = THM-1950's inverse-response s, whose bounds and velocity-addition law are KNOWN there), series composition is relativistic addition (THM-3322/THM-1950 in substance; K1 absorbing)`
36. **credit**
    - Old: `with sw_1 = 2 c3 at every order.`
    - New: `with sw_1 = 2 c3 at every order (KNOWN: Session159's cycle-packing form; the operator is first in Session150/159).`
37. **credit**
    - Old: `the least HC orders 7, 8 are KNOWN (Session210/211/212, THM-4889)`
    - New: `the least HC orders 7, 8 and their class counts 24, 228 are KNOWN (Session210/211/212, THM-4889)`
    - Update `audit PENDING` once the corrections are applied.

## 5. Not checked

- The literature priority of the new items (only the one web search above), and the "not in OEIS" claims.
- Session211 `score-permutohedron`'s literal HC = 7/21 switches, and THM-4889's (H, HC) data.
- The open parts: conjecture (a) and the rate (c). Only the finite data and the implications were audited.

## 6. Overall verdict

**PASS WITH CORRECTIONS.**

Genuine mathematical errors (candidates for a MISTAKE entry):
1. **HYP-9490 Consequence** (and the title, note §2, INDEX). The error term `2Δc3·ρ·(1 + o(1))` is false whenever
   Δc3 = 0 but a higher sw_j differs. Example: the order-5 classes with P = q⁵ + 8q³ + 2q and q⁵ + 8q³ + 6q. The
   correct conditional statement is `1 + 2Δc3·ρ + o(ρ)`. Its proof (Lemma B.1 plus (a) at d = n and d = n − 2) was
   also missing.
2. **Note §3, Paley primes.** "Every prime p where 2 generates exactly the quadratic residues" needs p ≡ 3 (mod 4).
   It fails for 17, 41, 97, ….
3. **Note §3, square sums.** The sentence misstates the counts: 3 and 10 are p(Q23) and p(Q25), not the "first two
   counts after the onsets".
4. **Note §3, "linearises the join".** q ↦ q + q′ + 2qq′ is not linear.
5. **HYP-9490 rate table** (minor numerical misstatement). "d = 4 agrees to within 1%" is false at N = 24 (3.0%)
   and N = 44 (1.2%).

The priority misses would also merit a MISTAKE line: φ = THM-1950's s, Session159's all-order sw_1 = 2c3 and gate
operator, and Session212's class counts. They do not affect any proof.
