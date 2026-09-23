# AMM 12592 uniform frontier: a Pólya-capacity gap C* ≥ 11/8, the golden constant as a natural boundary, and super-blocks that beat it at every tested scale

**Status.**
* **AUDIT (orchestrator, 2026-09-23) and PROMOTION.** Theorem A was re-derived by hand, with no gap found:
  * continuation by `|W_m| <= S^(d_m)`;
  * gluing on the connected set `Omega_0 ∩ Omega_1`, with Lemma A1's star-shapedness, vertical-interval connectivity and van Kampen;
  * the Catalan quotient and removability at `w = 1/4`;
  * Pólya, via capacity `1/R(W,0)` of the complement in `1/w`;
  * the Fatou–Gauss primitivity contradiction.

  Independent code, [`..._orchestrator_audit.py`](../../04-computation/experiments/amm12592_procgen_20260923_orchestrator_audit.py), re-verifies the `gamma = 7/20` two-disk certificate. Containment holds with a sampled maximum of `0.99833 < 1`, plus a Lipschitz margin. The closed form reproduces `log(8 sqrt3/9)` exactly at `c = 0` and gives `0.0059353920627 > 0`. The stored `gamma = 3/8` certificate re-verifies with `--verify-only`.

  Promoted as [THM-4467](../../01-canon/theorems/THM-4467-uniform-polya-capacity-gap-amm12592.md).
* **PROVED** (this note; the classical inputs are CITED):
  * **Theorem A (uniform gap).** Put `Λ(γ) := log R(W_γ, 0)` (defined in §3). If `Λ(γ) > 0`, then no exactly fair deterministic extractor has `T(n) ≤ (1+γ)n + D` for any constant `D`. The certificate at `γ = 3/8` (an explicit polynomial, §3.4) gives **`C* ≥ 11/8 = 1.375`**. A closed-form two-disk certificate at `γ = 7/20` gives `C* ≥ 27/20`.
    * This is the first positive gap that is uniform over all extractors.
    * It supersedes the corpus status "no uniform gap known" (ledger row; THM-3342 §1; MISTAKE-368), and it strictly generalizes THM-3342.
  * **Theorem A0 (one-point variant).** Pólya's theorem applied on `U_γ` itself gives `C* ≥ 1 + γ_0`, with `γ_0 ≈ 0.0955`.
  * **Theorem B (natural boundary).** Consider an extractor whose large dyadic blocks are separately balanced at horizon `2N`.
    * Its symmetrized 0-spine deviation is `poly + Σ_j ε_j w^(2^j)`, with `w = pq` and `ε_j = ±1`.
    * The continuation domain `W_γ` lies in the unit disk iff `γ ≥ γ* = log_5(φ²)`.
    * Hence the golden bound `C ≥ C_*` holds in that class. This is a third proof, after THM-3009 and Long §11.
    * **Corollary:** every extractor with `limsup T/L < C_*` has `φ` analytic at the golden point `w = −1` (that is, `p = −1/φ` and `pq = −1`).
  * **Lemma P.** A complement-symmetric extractor is fair iff `φ(w) := 2F(p(w)) − p(w) ∈ Z[[w]]` continues analytically across `w = 1/4`. Moreover `φ ≡ Σ_j w^(2^j) (mod 2)`.
  * **Lemma S (super-block reduction).** Complement-symmetric ratio-`B` super-blocks are equivalent to palindromic class polynomials of degree `(B−2)N`. Equivalently, they are the signed handoff states `S_N ∈ Z[w]` with `deg S_N ≤ (B−2)N/2` and `S_N(0) = ±1`. Long's fold then runs with the target `F_0(u) = ((1+u)/2)^(N−1) S_N((1−u²)/4)`.
* **FINITE-EXACT:**
  * **Long's verifier.** The verifier (SHA-256 `07801d62…ced61`) was re-run and reproduced byte for byte, including the extended `N = 1024` run and a re-exported certificate identical to the distributed one.
  * **Independent audit of Long's construction** for `N ≤ 512`, using our own identity, profile and deadline code. The lacunary identity `ε_N = +1` holds for `N ≤ 128`.
  * **Exact CP-SAT super-block table.** It includes `ρ_4(8) ≤ 28/19 < 3/2 < 14/9 = ρ_2(8)`.
  * **Three exactly verified integer super-blocks.**
    * `[64,256)` has `max T(L)/L = 53/34 = 1.55882`.
    * `[128,512)` has `max T(L)/L = 83/53 = 1.56604`.
    * `[256,1024)` has `max T(L)/L = 157/100 = 1.57`.
    * At those ratios, no separately balanced block `[N,2N)` exists. This is certified by Long's evaluation inequality at a single rational `θ`, in exact arithmetic.
  * **Exact fold thresholds** for real centers, `N = 64 … 2048`.
* **NUMERICAL:**
  * **`γ_1 = 0.3775246`.** So `C* ≥ 1.3775246` is exactly the reach of Theorem A. Two independent Laplace solvers agree to `10^−8` and match closed forms at `γ = 0`.
  * **Super-block asymptotics.** Optimal-contour (steepest-descent) bounds put the asymptotic fold threshold of the golden-zero family at about `1.570` (for `κ ≈ 0.09`). For the golden family itself they give `≈ C_* = 1.598`.
  * **Necessary thresholds of the same family** are about `1.5525`.
* **CITED:**
  * G. Pólya, Math. Ann. 99 (1928) 687–706 (rationality under transfinite diameter `< 1`).
  * P. Fatou, Acta Math. 30 (1906) (rational integer series `= A/B` with `B(0)=1`, `A`, `B` coprime in `Z[x]`).
  * The Hadamard–Fabry gap theorem.
  * C. D. Long, draft of 2026-09-20. Its uniform proof for `N ≥ 128` was read, not re-derived line by line. None of these were re-read on the web in this session.
* **OPEN:**
  * **Whether `C* < C_*`.** The finite-exact and numerical evidence says yes, at about `1.57`, via cross-annulus cancellation.
  * **Whether `C* < 3/2`.** Theorem A does not exclude it, and no evidence is available either way.
  * **The exact value of `C*`.** It lies in `[11/8, C_*]` (proved), narrowed numerically to `[1.3775, ≈1.570]`.
* **REFUTED:** no canonical theorem.
  * **Superseded:** the open status "no positive uniform lower gap is available" (ledger row, THM-3342 §1, MISTAKE-368), by Theorem A.
  * **Contradicted at block level, not refuted asymptotically:** the working expectation `C* = C_arch = C_*` (06-writeups remark 6; THM-3009 §4). Super-blocks beat every separately balanced block at scales 8, 64, 128 and 256 (FINITE-EXACT). The asymptotic claim `C* < C_*` is only NUMERICAL.

---

## 0. The question and the short answer

A deterministic procedure reads flips of a coin with unknown bias `p`. It must output an exactly fair bit, for every `p`, by flip `T(L)`, where `L` is the length of the initial constant run. The owner asks whether `C* = inf_extractors limsup T(L)/L` can be below `3/2`, or even `1 + ε`.

* **`1 + ε`: impossible for every `ε ≤ 3/8`.** This is PROVED with a certificate, so `C* ≥ 1.375`. It is also impossible for every `ε < 0.3775246`, but that part is NUMERICAL: it is the exact reach of the proved criterion, evaluated by a validated Laplace solver.
* **Below `3/2`: open.** The proved window is `11/8 ≤ C* ≤ C_* = 1.59799`.
  * The upper bound is Long's construction, reproduced here.
  * Separately balanced dyadic blocks cannot go below `C_*`.
  * Super-blocks, in which imbalance crosses the inner dyadic boundaries, beat the separately balanced class at every scale we tested. This is exact for integer blocks at `N = 8, 64, 128, 256`, and holds at the real-center level through `N = 2048`. Asymptotically they appear to reach `≈ 1.570` (numerical, not proved). This is still above `3/2`.

---

## 1. Setting

Bits are i.i.d. with `P(0) = p`, `q = 1 − p`. `L` is the initial run and `T(L) ≥ L + 1` the pathwise deadline. By the spine normal form (THM-2966), put `d_m = T(m) − m − 1`. A `T`-deadline fair extractor is equivalent to integers `0 ≤ w_{m,k}, v_{m,k} ≤ C(d_m,k)` with

```text
F(p) + G(1-p) = 1/2    on (0,1),
F(p) = sum_m p^m q W_m(p),        W_m(p) = sum_k w_{m,k} p^{d_m-k} q^k,
G(u) = sum_m u^m (1-u) V~_m(u),   V~_m(u) = sum_k v_{m,k} u^{d_m-k} (1-u)^k.
```

Both `F ∈ Z[[p]]` and `G ∈ Z[[u]]` have zero constant term.

---

## 2. Reconciliation with C. D. Long's draft (2026-09-20)

Source: GitHub `long-mathematics/deterministic-von-neumann-fair-extractor`, commit `6c5b433` (2026-09-22). The fetched files are in `scratch/procgen_amm/long_draft/`, with these SHA-256 digests:

* `.tex`: `4edfc10b…dcc2fef2`
* verifier: `07801d62aedf91b87e122bd96f53a9bd450d64dd4e265caaff364f3da0dced61`. This equals the digest printed in the manuscript and README.
* certificate JSON: `1be015fd…84cb2a9`

The verification runs and their results:

* `python3 scripts/verify_glazer_critical.py`: all exact checks pass.
* `--check-certificate`: all six blocks pass.
* `supplementary_checks.py`: all pass. This includes word-level realization for `N ≤ 16`.
* `--extra 128 256 512 1024`: 471,442 sites at `N = 1024`.
  * The output is identical to the recorded `verification_output.txt`.
  * The re-exported certificate is byte-identical to the distributed one.
  * Peak RSS was 235 MB.
* Our independent audit (`…_long_audit.py`) passes:
  * capacity and parity;
  * the packet identity, by an `O(N²)` Horner recursion (a different algorithm from his radix evaluation);
  * the deadlines, with `floor(βm)` computed by 80-digit mpmath and a separation guard;
  * `T ≤ ⌈C_*L⌉` and `T ≤ C_*L − log_5 L + 6`;
  * consistency with his own lower bound.
  * All pass for every dyadic `N ≤ 512`.

| Long's item | Corpus counterpart | Verdict |
|---|---|---|
| Thm 1.1: a fair extractor with `T(L) ≤ min(2N, ⌈C_*L⌉ − s_N)` for all `L`; `C_* = 1 + log_5 φ² = 1.5979874…` | THM-3029/3302/3329/3330 attain the golden floor only for `n ≤ 2047`, with slack `D0*(R) = 0, 1, 5, 15, …`. The all-`R` statement is OPEN in the corpus, and THM-3330 says plain rule A needs *linear* slack. | **NEW; closes the corpus's all-`R` attainment item.** Fold plus boundary-preserving lattice rounding reaches slack `D0 = −s_N ≤ 0`, where rule A failed. Reproduced by us for `N ≤ 1024`, and independently audited for `N ≤ 512`. |
| Log saving `T ≤ C_*L − log_5 L + 6`; Cor. 11.4 `1 ≤ Λ_log ≤ 2` | none | **NEW** |
| Thm 11.1: evaluation obstruction at `x = −θ*`, `θ* = φ^−2`, without complement symmetry. Cor. 11.3: `C_*` is sharp for separately balanced horizon-`2N` rules | THM-3009 (ARCH, the Taylor expansion at `u = −1`; its Stirling transfer was the declared residual debt) and THM-3027 | **Same constant, better proof.** It is one evaluation point. It discharges THM-3009's asymptotic debt for the ratio-2 class. In `u` coordinates the point is `u = √5`; in `w = pq` it is `w = −1` (§4). |
| §12: `o(N)` extra horizon still forces `C ≥ C_*` | THM-3024 (demoted, MISTAKE-361) | **NEW**, and consistent. Long explicitly excludes cross-block cancellation, which is exactly where §6 finds improvement. |
| §13: fairness for exchangeable sources; type-conditional fairness | none | **NEW.** It follows immediately from type balance. |
| §14: at horizon `2N`, dyadic `N` is necessary and sufficient; no two-bit output | THM-3007 (the two-parameter `[N,N+l)` version), THM-3343, THM-2160 §6.1 | **Overlap:** the special case `l = N` of THM-3007. The two-bit corollary is new but easy. |
| Remark: "not globally optimal when blocks can cancel" | ledger "uniform frontier open", HYP-9061 | **Agrees.** §6 of this note supplies finite-exact and numerical evidence that global schemes beat `C_*`. |
| History and conventions: `n_SE = L+1`; Wang's `7/4` refinement invalid; Wagon–Winkler `2L−1` | `06-writeups/amm12592-solution.tex` uses `L` | Agrees. |

**Disagreements: none of substance.**

* THM-3330's "plain rule A needs linear slack" remains true for rule A. It must not be read as "the golden floor profile needs slack": Long's rule closes every block with nonpositive slack.
* THM-3009's scope sentence "balanced block schemes" should be read as ratio-2 blocks (its own §11.2 says so). Ratio-4 super-blocks are strictly better at every scale tested (§6).

---

## 3. Theorem A: a uniform capacity gap

### 3.1 Domains

Write `S(p) = |p| + |1−p| ≥ 1`. Define

```text
Ω_0(γ) = { |p| S(p)^γ < 1 },            Ω_1(γ) = 1 − Ω_0(γ),
U_γ    = Ω_0 ∪ Ω_1 = { min(|p|,|1-p|) · S(p)^γ < 1 },
W_γ    = { p(1-p) : p ∈ U_γ },          Λ(γ) = log R(W_γ, 0),
```

where `R(W,0)` is the conformal radius at `0`.

**Lemma A1 (geometry).**
* `Ω_0` is star-shaped with respect to `0`: `|p| S^γ` is strictly increasing along rays, since `d/dt |1 − te^{iθ}| ≥ −1`.
* `Ω_0 ∩ Ω_1 = {max(|p|,|1−p|) S^γ < 1}` meets each vertical line `Re p = x` in an interval centred on the real axis, and that interval is nonempty iff `0 < x < 1`. So the intersection is connected.
* `U_γ` is therefore simply connected (van Kampen). It is invariant under `s(p) = 1 − p` and under conjugation.
* `w = p(1−p)` maps `U_γ` onto `W_γ` as a proper 2:1 map, branched only at `p = 1/2`. `W_γ` is simply connected and bounded.

### 3.2 Statement and proof

**Theorem A.** If `Λ(γ) > 0`, then no exactly fair extractor satisfies `T(n) ≤ (1+γ)n + D` for all `n`.

*Proof.*

1. **Continuation.** For complex `p`, `|W_m(p)| ≤ Σ_k C(d,k)|p|^{d−k}|q|^k = S(p)^{d_m}`. Since `d_m ≤ γm + D'`, `|p^m q W_m(p)| ≤ |q| S^{D'} (|p| S^γ)^m`. So the series for `F` converges locally uniformly on `Ω_0(γ)`, and its Taylor series at `0` is the formal `F`. The same holds for `G` on `Ω_0(γ)` in its own variable.
2. **Gluing.** `F(p) + G(1−p) = 1/2` holds on `(0,1) ⊂ Ω_0 ∩ Ω_1`. By Lemma A1 and the identity theorem, it holds on `Ω_0 ∩ Ω_1`. So `F` extends to `U_γ` by `F := 1/2 − G(1−·)` on `Ω_1`. Likewise `G` extends, and `F(z) + G(1−z) = 1/2` on `U_γ`.
3. **Symmetrization.** Put `Δ := F − G` and `Σ := F + G − 1/2`. Then `Δ∘s = Δ` and `Σ∘s = −Σ` on `U_γ`. Also `Δ ∈ Z[[z]]` and `Σ ∈ −1/2 + zZ[[z]]`.
4. **Quotient.** `z(w) = (1 − √(1−4w))/2 = Σ_{n≥1} Cat(n−1) w^n ∈ wZ[[w]]` inverts `w = z − z²`. The even functions descend:
   * `φ(w) := Δ(z) ∈ Z[[w]]`;
   * `ψ(w) := Σ(z)/(1−2z) ∈ −1/2 + Z[[w]]`. This uses `1/(1−2z) ∈ Z[[z]]`; the singularity at `z = 1/2` is removable because `Σ(1/2) = 0`.
   * Both `φ` and `ψ` are analytic on `W_γ`.
5. **Pólya's theorem** (1928): an integer power series that is analytic on a domain containing `0` with conformal radius `> 1` is rational. Equivalently, in the variable `1/w`, it is regular off a compact set of transfinite diameter `1/R < 1`. So `φ` and `ψ + 1/2` are rational. Hence `F = (φ(z(1−z)) + (1−2z)ψ(z(1−z)) + 1/2)/2` is rational, and so is `G(u) = 1/2 − F(1−u)`.
6. **Fatou–Gauss endgame.**
   * By Fatou, `F = A/B` and `G = C/D`, with `A,B` coprime in `Z[p]`, `C,D` coprime in `Z[u]`, and `B(0) = D(0) = 1`.
   * Then `C(1−p)/D(1−p) = (B−2A)/(2B)`, and both sides are reduced.
   * So `D(1−p) = 2λB` for some `λ ∈ Q^×`.
   * Both `D(1−p)` (whose value at `p = 1` is 1) and `B` (with `B(0) = 1`) are primitive, so `λ = ±1/2`.
   * Then `C(1−p) = ±(B−2A)/2 ∈ Z[p]` forces `B ∈ 2Z[p]`, which contradicts `B(0) = 1`. ∎

**Lemma A2 (quotient Robin constant).** `Λ(γ) = V(U_γ,0) + g_{U_γ}(0,1)`, where `V` is `log R` and `g` is the Green function.

*Proof.* For the proper 2:1 map `w`, `g_W(w(z),0) = g_U(z,0) + g_U(z,1)`. Expand at `z → 0`. ∎

**Theorem A0.** Applying step 5 to `F` on `U_γ` gives the same conclusion whenever `R(U_γ,0) > 1`. This is the one-point version, with `γ_0 = 0.0954903`. THM-3342 is its `γ → 0` shadow: `d_m = o(m)` implies `d_m ≤ γm + D` for every `γ > 0`.

**Why the symmetrization matters.** Integrality holds at both `p = 0` and `p = 1`. The 2:1 quotient turns this two-point integrality into one-point integrality on `W_γ`. The Green interaction `g_U(0,1) ≈ 0.20–0.35` is exactly what the quotient adds.

### 3.3 Numbers

Computed by a lightning Laplace solver: polynomial terms plus poles clustered at the reentrant corners `1/2 ± i·(2^{−2γ/(1+γ)} − 1/4)^{1/2}`. Two routes were used, on `U` (giving `V` and `G`) and directly on `W`; they agree to `10^−8`.

Validation at `γ = 0`, where `U_0 = D(0,1) ∪ D(1,1)`, against closed forms:
* `R(U_0,0) = 4√6/9`;
* `g(0,1) = ½ ln 2`;
* `R(W_0,0) = 8√3/9`.

| γ | V = log R(U,0) | g_U(0,1) | Λ(γ) | C = 1+γ |
|---|---|---|---|---|
| 0 | +0.084950 | 0.346574 | +0.431523 | 1 |
| 0.0955 | ≈0 | 0.298830 | +0.298822 | 1.0955 (one-point limit) |
| 0.2 | −0.082048 | 0.256534 | +0.174486 | 1.2 |
| 0.3 | −0.152089 | 0.223458 | +0.071369 | 1.3 |
| 0.375 | −0.200166 | 0.202387 | +0.002221 | 1.375 |
| 0.3775 | −0.201710 | 0.201732 | +0.000022 | 1.3775 |
| **γ_1 = 0.3775246** | | | **0** | **1.3775246 (quotient limit)** |
| 0.5 | −0.273299 | 0.172880 | −0.100420 | 1.5 |
| γ* | −0.325500 | 0.153677 | −0.171823 | 1.598 |
| 1 | −0.505197 | 0.098924 | −0.406273 | 2 |

**Consistency checks.**
* `Λ < 0` at `γ*` and at `γ = 1`, where constructions exist (Long's; the classical `2n`). So the obstruction does not fire there.
* The thresholds are stable to `10^−8` under three discretizations.

### 3.4 Certificates (rigorous)

* **`γ = 7/20` (closed form).** Take `V = D(c,ρ) ∪ D(1−c,ρ)` with `c = 13/200` and `ρ = 39/50`.
  * Containment `D(c,ρ) ⊂ Ω_0`: `Ω_0` is star-shaped and `0 ∈ D`, so it suffices to check the boundary circle. This was done by adaptive interval arithmetic in 172 pieces; the upper bound is `0.99990 < 1`.
  * Möbius map plus power map give, exactly,

    ```text
    log R(w(V),0) = log[ tan(A)(1/4+s²)/(κ s) ],
    s = (ρ²−(1/2−c)²)^{1/2},  β₀ = atan(s/(1/2+ρ−c)),
    κ = π/(2π−4β₀),  A = 2κ(atan 2s − β₀).
    ```

    At `c = 0`, `ρ = 1` this returns `log(8√3/9)`.
  * The interval evaluation gives `[0.0059354, 0.0059354] > 0`. Hence **`C* ≥ 27/20`**.
* **`γ = 3/8` (subordination).** Schwarz's lemma applies to `Φ_W∘P`. If `P` is analytic with `P(0) = 0` and `P(∂D) ⊂ W_γ` (so `P(D) ⊂ W_γ` by the winding argument, since `W_γ` is simply connected), then `R(W_γ,0) ≥ |P'(0)|`.
  * The certificate `P` is the truncated Riemann map at radius `0.9985`. It has 2047 real float64 coefficients, stored as hex in `amm12592_procgen_20260923_subordination_certificate.json`.
  * `b_1 = P'(0) = 1.00071974 > 1`, and `Σ k|b_k| = 1.91152`.
  * Membership `w ∈ W_γ ⟺ min(|p|,|1−p|)S^γ < 1` for `p = (1−√(1−4w))/2`. This was checked on `2^21` FFT nodes with explicit bounds:
    * node spacing and evaluation allowance: `δ = 2.865·10^−6`;
    * root displacement: `|Δp| ≤ 2δ/|1−2p|`, with `min|1−2p| = 1.319`;
    * Lipschitz constant: `S_max^γ + 2γ m_max`.
  * The maximum certified bound is `0.99990624 < 1`. Hence **`C* ≥ 11/8`**.
  * `--verify-only` re-checks the stored certificate on its own.
  * Rigor caveat. The `7/20` certificate is pure interval arithmetic.
    * The `3/8` certificate adds a floating-point FFT evaluation. The allowance `10^−9` used there dominates the standard a-priori FFT bound of about `10^−10` (roughly `log2(M) · u · √M · ‖b‖₂`).
    * An interval re-evaluation of the 2047-term polynomial would make it fully machine-checked.
    * Either way, the conclusion `C* > 1.35` rests only on interval arithmetic.

---

## 4. Structure: parity, lacunarity, and the golden point as a natural boundary

**Lemma P.** For complement-symmetric rules, `G(u) = u − F(u)`, so fairness reads `Δ(p) = Δ(1−p)` with `Δ = 2F − p`. So fairness holds iff `φ(w) = Δ(z(w))` has no branch point at `w = 1/4`. Moreover:
* `F_dev := Σ_m p^m q E_m` with `E_m = 2W_m − 1 ≡ 1 (mod 2)` in the Bernstein `Z`-basis, so `2F − p ≡ p (mod 2)`.
* `z(w) = Σ Cat(n−1)w^n` and `Cat(n−1)` is odd iff `n` is a power of 2 (checked for `n ≤ 4096`).
* Hence `φ ≡ Σ_j w^(2^j) (mod 2)`. Its reduction satisfies `X² + X = w` over `F_2(w)`, an Artin–Schreier equation with no rational root. This gives an independent reason why `φ` is never rational; it is the arithmetic input of Theorem A.

**Lacunary identity.** Long's Thm 11.1 gives `P_+ = ±1`, so every separately balanced horizon-`2N` block has 0-spine deviation exactly `ε_N (pq)^N`. For his construction `ε_N = +1`, verified as an exact polynomial identity for `N ≤ 128`. So `φ(w) = Σ_j ε_j w^(2^j)` up to a polynomial. This is a Hadamard-gap series, and the circle `|w| = 1` is its natural boundary.

**Theorem B.**
* `max_{closure W_γ} |w| = r_π(1+r_π)`, where `r_π(1+2r_π)^γ = 1`.
  * On the arc, `st = s^{1−1/γ} − s²` decreases in `s = |p|`.
  * `r_π = min_θ r_γ(θ)`.
* `r_π(1+r_π) ≤ 1` iff `r_π ≤ 1/φ` iff `√5^γ ≥ φ` (using `1 + 2/φ = √5`) iff `γ ≥ γ*`.

The lacunary `φ` must be analytic on `W_γ`, so `W_γ ⊂ D`, and therefore `γ ≥ γ*` for separately balanced schemes. At `γ = γ*`, `∂W_γ` touches the unit circle exactly at `w = −1`; this was checked numerically on a `γ` grid. So the golden constant is precisely the moment the continuation domain pokes out of the natural boundary.

**Corollary (golden-point analyticity).** If `limsup T/L < C_*`, then `W_γ ∋ −1` for some `γ < γ*`, so `φ` is analytic in a neighbourhood of `w = −1`. Any sub-golden extractor must therefore break the unit-circle natural boundary at the golden point. This rules out all lacunary schemes, and every scheme whose block contributions do not decay near `w = −1`.

---

## 5. Obstruction tools: failure anatomy

1. **2-adic structure (Kummer and Lucas carries).**
   * Inside a dyadic shell, parity cancels identically (THM-3009 Reduction A).
   * Globally it yields exactly one fact, `φ ≡ Σ w^(2^j) (mod 2)` (§4). That fact is irrationality, which is sufficient input for Pólya but carries no rate.
   * No 2-adic overconvergence is available. The 0-spine converges 2-adically only on `|p|_2 < 1` and the 1-spine only on `|1−p|_2 < 1`; these are disjoint. The 2-adic radius in `w` is exactly 1. So a Borel–Dwork product formula gains nothing at the prime 2.
2. **Capacity.**
   * The single-point Pólya bound gives `0.0955`.
   * The two-point quotient gives `0.3775`. This is the Cantor-capacity-type gain from integrality at both `0` and `1`.
   * The method is sharp for general integer series. Going beyond needs the Bernstein-box and positivity structure of the spine polynomials, not just continuation.
   * ARCH and evaluation capacity bounds are strictly block-local. Cross-shell Hall cuts fail (MISTAKE-361: deeper shells absorb any deficit with exponential room).
3. **Entropy/rate certificates** (the decoded artanh gate, HYP-9061): single-bias and rate arguments cannot see integrality, and Theorem A shows integrality is the whole obstruction.
4. **Spine normal form:** it is the starting point of Theorem A, so it did not fail. Its only limitation is that the analytic domain comes from the crude box bound `S(p)^{d_m}`.
5. **What would push the floor up.** A lower bound above `1.3775` must exclude the super-block states of §6. For the golden-zero family, the necessary conditions — natural boundary plus single-block evaluation at every complex point — are satisfied down to `C ≈ 1.5525`. Any proof that `C* ≥ 3/2` must therefore use more than analyticity and pointwise box bounds.

---

## 6. Construction side: super-blocks and the handoff state

**Lemma S.** A ratio-`B` super-block is `[N, BN)`, decided by flip `BN` and balanced as a whole. With Long's complement pairing, it is equivalent to the class polynomial `P_+(x) = Σ_i x^i (1+x)^{a_i} E_i(x)` being palindromic of degree `(B−2)N`, where `|[x^r]E_i| ≤ C(R_i,r)`, `[x^r]E_i ≡ C(R_i,r) (mod 2)`, and `P_+(0) = ±1`.
* Equivalently, the block's 0-spine deviation is `w^N S_N(w)`. The **handoff state** is `S_N ∈ Z[w]`, with `deg S_N ≤ (B−2)N/2`, `S_N(0) = ±1` and `S_N ≡ 1 + w^N (mod 2)` for `B = 4`.
* For `B = 2`, `S_N = ±1` is forced: the lacunary, golden case.
* In Long's `u` variable the fold target is `F_0(u) = ((1+u)/2)^(N−1) S_N((1−u²)/4)`. Here `1 + w = (5−u²)/4` vanishes at the golden point `u = √5`.

This is the smallest signed state space for the handoff. It is one integer polynomial per super-block, and it replaces the single corner sign `ε_N`.

### 6.1 Exact small super-blocks (CP-SAT; every solution re-verified exactly)

| B | N | best max T(L)/L on [N, BN) | status |
|---|---|---|---|
| 2 | 2, 4, 8, 16 | 3/2, 3/2, 14/9, 25/16 | exact; matches THM-3008 |
| 4 | 2, 4 | 3/2, 3/2 | exact (ratio-2 blocks give 14/9 on [4,16)) |
| 4 | 8 | ≤ 28/19 = 1.4737; 29/20 proved infeasible | 25/17 timed out, so the value lies in {16/11, 19/13, 22/15, 25/17, 28/19} (ratio-2 blocks give 25/16 on [8,32)) |
| 8 | 2, 4 | 3/2, 3/2 | exact (ratio-2 blocks give 25/16 on [4,32)) |

For `B = 4` the optimal handoff states all vanish at the golden point:
* `N = 2`: `S = w² − 1`;
* `N = 4`: `S = −(1+w)²(1+w²)`;
* `N = 8`: `S = −(1+w)²(1+w²)(1+w⁴)`, i.e. `S = −(1+w)(1−w^N)/(1−w)`.

For `B = 8` (`N = 2, 4`) the optimum `3/2` equals the `B = 4` value. The solver's states, for example `−(1 + w² + 2w³ + 2w⁴ + 2w⁵ + w⁶)` at `N = 2`, have `S(−1) = −1`: at these tiny scales the golden point does not bind.

Chaining this family gives `φ = −((1+w)/(1−w)) Σ_k (−1)^k w^(2^k)`, which still has the natural boundary `|w| = 1`: a rational factor cannot remove a Hadamard natural boundary, because the gap series admits no meromorphic continuation across any arc. So this family is sub-golden only at finite scale; by the argument of Theorem B it cannot beat `C_*` asymptotically.

### 6.2 Golden-zero families

The golden-zero families are

```text
S_N = (1+w)^(κN) (1−w^N)/(1−w^(κN)),
```

with `κN` a power of 2. They satisfy the parity condition and have a zero of order `κN` at `w = −1`.

**(a) Necessary conditions**, asymptotic, from a grid computation. The chained series must converge on `W_γ` (natural boundary). Long's evaluation inequality must also hold at every complex point, using both preimages:

| κ | natural boundary: C ≥ | single-block evaluation: C ≥ |
|---|---|---|
| 0 (lacunary) | 1.5980 | 1.5975 (grid) |
| 1/16 | 1.5634 | 1.5632 |
| 1/8 | 1.5528 | 1.5523 |
| 0.15 | 1.5552 | 1.5552 |
| 1/4 | 1.5803 | 1.5801 |

The golden obstruction is gone. (The full power `(1+w)^N`, κ = 1, is bad: its lemniscate `|w(1+w)| < 1` forces `C ≥ 1.959`.)

**(b) Exact fold thresholds.** These are real centers from Long's recursion with the target above, computed with exact integers. The criterion is: `ℓ1` masses `≤ 1` on levels `≥ 1`, and an exact Krawtchouk box at level 0. Profile `⌈C·L⌉`, `B = 4`.

| N | S=1 (golden family) | κ=1/32 | κ=1/16 | κ=1/8 |
|---|---|---|---|---|
| 64 | 1.5758 | 1.5649 | 1.5515 | 1.5539 |
| 128 | 1.5834 | 1.5683 | 1.5654 | 1.5598 |
| 256 | 1.5902 | 1.5739 | 1.5662 | 1.5674 |
| 512 | 1.5932 | 1.5753 | 1.5681 | 1.5712 |
| 1024 | 1.5952 | 1.5769 | 1.5691 | 1.5736 |
| 2048* | 1.5964 | 1.5775 | 1.5698 | 1.5748 |

* The `S = 1` column tends to `C_* = 1.5980`. A geometric extrapolation of its increments gives `1.5983`.
* The golden-zero columns saturate lower. The best golden-zero entry trails the golden column by `0.024` at N = 64, 128 and 256, then by `0.025`, `0.026` and `0.027` at N = 512, 1024 and 2048. The gap is not closing.
* (*) The `N = 2048` row comes from a scratch run of the same algorithm with bisection bracket `[1.54, 1.62]`. The default run stops at `N = 1024`; use `--fold-nmax 2048` for the full table (about 15 min).

**(c) Asymptotics by optimal contours.**
* The fold coefficients are the contour integrals

  ```text
  c^{(i)}_s = (−2)^i (2πi)^{−1} ∮ F_0(u) u^{−s−1} (u−1)^{−i} du
  ```

  over any contour enclosing `0` and `1`.
* The sharp exponent is the mountain-pass level of `Ψ = log|F_0|/N − v log|(1−u)/2| − σ log|u|` between `{0,1}` and `∞`, computed by grid percolation.
* The resulting asymptotic sufficient thresholds:
  * `S = 1`: **1.5967** (≈ `C_*`, within the grid error of about `10^−3`);
  * `κ = 1/32`: **1.5794**;
  * `κ = 1/16`: **1.5720**;
  * `κ = 0.09`: **1.5700**;
  * `κ = 1/8`: **1.5750**.
* The finite-`N` columns of (b) approach these values from below, as the golden column approaches `C_*`.
* Circles are not enough. For every `κ` tested, the circle-only Cauchy bound gives no better than about `1.598` (`1.5979–1.63`, exploratory scan), because the relevant saddles are complex and the target has oscillating sign. Deformed contours are essential.

**(d) Integer super-blocks (exact).**
* We took the fold center, applied Long's lattice rounding (his function, imported unmodified), and re-verified the result independently: box, parity, palindromic class polynomial, and equality with the target. This gives:
  * `[64, 256)`, `κ = 1/16`: `max T(L)/L = 53/34 = 1.55882`;
  * `[128, 512)`, `κ = 1/8`: `max T(L)/L = 83/53 = 1.56604`;
  * `[256, 1024)`, `κ = 1/16`: `max T(L)/L = 157/100 = 1.57`, with 171,545 packet sites.
* No separately balanced block reaches these ratios, and this is proved exactly. Long's inequality `Σ_i θ^i(1−θ)^{a_i}(1+θ)^{R_i} ≥ 1` fails for the ratio-2 profile `min(2N, ⌊C L⌋)`:
  * sum `= 0.6258` at `θ = 42/125` for `N = 64`;
  * sum `= 0.3134` at `θ = 181/500` for `N = 128` at `C = 83/53`;
  * sum `= 0.0429` at `θ = 93/250` for `N = 256` at `C = 157/100`.
  * THM-3009's ARCH bounds `ρ_2(64) > 1.5753`, `ρ_2(128) > 1.5828` and `ρ_2(256) > 1.5887` agree.
* So cross-annulus cancellation strictly beats separate balance at scales 64, 128 and 256. The margins are `0.017`, `0.017` and `0.019`, measured against ARCH. This is FINITE-EXACT.

**Honest extrapolation.** Three independent computations agree that ratio-4 super-blocks with a golden-point zero of order about `N/11 … N/16` have real centers at `C ≈ 1.570–1.572`:
* exact finite folds up to `N = 2048`;
* optimal-contour asymptotics;
* integer rounding that succeeds wherever it was tried (`N = 64, 128, 256`).

That is below `C_* = 1.598` by `0.026–0.028`. It is not a proof. The missing pieces are:
1. a rigorous, uniform-in-`N` version of the contour bound, i.e. explicit contours and interval certificates;
2. Long's rounding lemma for palindromic targets. His proof uses only the kernel structure and the profile shape; the computations suggest it transfers verbatim, but this was not written out.

Nothing seen so far goes below `3/2` asymptotically. The best asymptotic family value is about `1.570`, and its necessary-condition floor is about `1.5525`.

---

## 7. Best current answer to the owner

| Claim | Label |
|---|---|
| `C* ≥ 11/8 = 1.375`: no extractor has `limsup T/L < 11/8`, so a constant of `1 + ε` is impossible for `ε ≤ 3/8` | **PROVED** (Theorem A plus certificate; Pólya 1928 and Fatou 1906 CITED) |
| `C* ≥ 1.3775246` | **NUMERICAL** (the exact reach of Theorem A; the certificate covers `3/8`) |
| `C* ≤ C_* = 1 + log_5 φ² = 1.5979874…`, with `T(L) ≤ C_*L − log_5 L + 6` | **CITED** (Long) and **FINITE-EXACT** (reproduced to `N = 1024`; independently audited to `N = 512`) |
| `C_*` is optimal for separately balanced dyadic blocks | **PROVED** (THM-3009, Long §11, Theorem B) |
| Cross-annulus super-blocks beat every separately balanced block at scales 8, 64, 128 and 256 | **FINITE-EXACT** |
| `C* ≈ 1.57 < C_*` via ratio-4 super-blocks | **NUMERICAL** (not proved) |
| `C* < 3/2`? | **OPEN.** Not excluded by any known obstruction; no construction below `3/2` beyond tiny scales. |

So: **1 + ε is impossible for ε ≤ 3/8 (and, numerically, for ε < 0.3775), and 3/2 is open.** The provable window is `[1.375, 1.598]`. The evidence points to about `1.57` from above.

---

## 8. Hypothesis candidates (not filed; for the orchestrator)

* **H1 (super-blocks beat golden).**
  * Claim: for `κ = 1/16`, the states `S_N = (1+w)^{N/16}(1−w^N)/(1−w^{N/16})` on super-blocks `[N,4N)`, `N = 4^k`, yield integer fair extractors with `limsup T/L ≤ 1.573`. Hence `C* < C_*`.
  * Proof plan: (i) explicit contours for the fold integral, from the §6(c) mountain-pass level sets, checked with interval arithmetic uniformly in `v = i/N`; (ii) Long's five-unit lattice rounding lemma (§6 of his draft) for palindromic targets; (iii) finite certificates for `N < N_0`.
* **H2 (optimal handoff).** The infimum over all ratio-`B` handoff states (zero measures `μ` with mass `≤ (B−2)/2` and `S(0) = ±1`) of the asymptotic fold threshold decreases in `B`. Its limit `C_∞` lies in `[1.3775, 1.570]`. Deciding whether `C_∞ < 3/2` is a weighted-potential optimisation over `μ`: minimise the mountain-pass threshold.
* **H3 (Pólya sharpness).** Is `C* = 1 + γ_1 = 1.3775…`? This holds iff integer series with capacity exactly `1` on `W_{γ_1}` can be realised as spine functions. It is probably false (§5.5), but it is the cleanest target.
* **H4 (refined obstruction).** Combine Theorem A with the per-block evaluation inequality. Use Hankel determinants of `φ` restricted to block polynomials with `S_N(0) = ±1`. The goal is to exclude `C < 1.55` for all ratio-4 super-blocks.

## 9. Files and reproduction

* `04-computation/experiments/amm12592_procgen_20260923_long_audit.py` — the audit of Long's draft (§2, and the lacunary identity in §4).
* `04-computation/experiments/amm12592_procgen_20260923_polya_capacity.py` — Theorems A and A0, the tables, both certificates, Theorem B numerics, the parity check (§3, §4). The option `--verify-only` re-checks the stored certificate.
* `04-computation/experiments/amm12592_procgen_20260923_superblocks.py` — §6(a)–(d): CP-SAT, states, necessary thresholds, exact folds (`N ≤ 1024` by default), integer super-blocks, and contour asymptotics.
* `04-computation/experiments/amm12592_procgen_20260923_run.sh` — fetches Long's files if absent (generic user agent) and runs everything.
* `05-knowledge/results/amm12592_procgen_20260923.out` — the full log, with script SHA-256 digests.
* `05-knowledge/results/amm12592_procgen_20260923_subordination_certificate.json` — the `γ = 3/8` certificate polynomial.
* The `N = 2048` fold row was computed with the same code in a separate run (`scratch/procgen_amm/ilp/fold_fast.py`; the log is `scratch/procgen_amm/ilp/fold_fast_big.out`), not in the default run.
* Runtime of `run.sh`: about 10 minutes. Each step stays below 300 MB resident; the largest is Long's `N = 1024` verifier run at 235 MB.
* SHA-256 digests (as recorded in the `.out` header):
  * `long_audit.py` `b6faff92…a5e0be15`
  * `polya_capacity.py` `f01a4fd7…71d1fd3d`
  * `superblocks.py` `2daf65ac…bc98ced83`
  * `run.sh` `42062b69…b6149eb8`
  * certificate JSON `5eccbc0c…cf5a5fc5`
  * output `f3c4528a…17a92048`
