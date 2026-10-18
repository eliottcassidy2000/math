# Half-period twin classes and the CM harmonic log law

**Session:** opus-2026-10-18-S43. **Status:** PROVED (THM-4685, THM-4686) / CITED (KLMSY 2016) / OPEN (HYP-9310(K)
and (c), HYP-9311) / FINITE-EXACT. Follows THM-4669 and THM-4670 in the
[factorization-vertex-SUSY note](factorization_vertex_supersymmetric_bridge_20261018.md), section 3. Audited twice:
s43_audit_20261018/REPORT.md and s43b_audit_20261018/REPORT.md, MISTAKE-629 and MISTAKE-632.

**Correction after the S43b audit (MAJOR M1).** Several of the session's later findings are already in
Kibelbek-Long-Moss-Sheller-Yuan, *Supercongruences and complex multiplication* (J. Number Theory 164 (2016) 166-178,
arXiv:1210.4489):
- the CM harmonic law `(L_p)` of section 3 is their **Corollary 3**;
- the census twin `(H)`, i.e. classes 3 and 7, is their **Corollary 4**;
- the Lucas-block equivalence is their **Lemma 10** and Section 4.

The session cited this paper for its Theorem 2 but had not read it beyond the abstract. As a result these items were
filed as open conjectures and "pre-registered" predictions; they are now credited. What remains new here:
- the census itself (section 1);
- the three 4C classes (THM-4685(i));
- the law at `m = 64` (outside KLMSY) and the symmetry proof at `m = -64`;
- an alternative proof of Corollary 3 from KLMSY Theorem 2 across the 2-isogeny, with the equivalent
  `L(m) ⟺ F_m ≡ (-m/p) F_(4096/m) (mod p^2)` (THM-4686);
- the shifted-square twin `T_m` (no prior statement found);
- the four-family extension and the detector (HYP-9310(c), HYP-9311);
- the mod-`p^3` block statement at ordinary primes (HYP-9311(O)).

Several "pre-registered" claims below cannot be ordered by git: script and output were committed together. Only
HYP-9311(O) has a verifiable prediction-before-output commit order.

## 1. The census

Two integer sequences `u` and `v` are **half-period twins** if, at `n = (p-1)/2`,
`u(n) ≡ ε v(n) (mod p^2)` for every prime `p`, with `ε ∈ {1, (-1)^n}`, although `u ≠ v`.
- THM-4669 is the twin `A ~ (-1)^n s_7`.
- THM-4670 is the twin `b ~ Σ C(n+k,k)^3`.

[`halfperiod_twin_classes_20261018.py`](../../04-computation/experiments/halfperiod_twin_classes_20261018.py) generates
every product of 3 or 4 factors from 10 binomials in `(n, k)`. That gives 935 sums and 648 distinct nonzero sequences.
The script then groups the sequences by their residues at the half-period, both signs, and confirms each group on all
primes `5 <= p <= 200`. There are exactly **seven** twin classes:

| class | members (`Σ_k` of) | value mod `p^2` | governing form |
|---|---|---|---|
| 1 | `C(n,k)^2 C(n+k,k)` = Apéry `b`; `C(n+k,k)^3` | `F` | 4C `η(4τ)^6` |
| 2 | `C(n+k,k) C(2k,k) C(2n-2k,n-k)`; `C(n,n/2)^2 [n even]` | `4^n F` | 4C |
| 3 | `C(2k,k)^3`; `C(2n-2k,n-k) C(n+k,2k)^2` | `a(p)` | 7AB `η(τ)^3 η(7τ)^3` |
| 4 | `C(n,k)^3 C(2n,n+k)`; `C(n,k) C(n+k,k)^2 C(2n,n+k)` | `16^n F` | 4C |
| 5 | Apéry `A`; `(-1)^n s_7` | `a(p)` | 4A `η(2τ)^4 η(4τ)^4` |
| 6 | `C(n+k,k)^3 C(n,2k)`; `(-1)^n C(n+k,k) C(2k,k) C(2k,n) C(2n-2k,n-k)` | `±a(p)` mod `p` only | 4A |
| 7 | `C(2n-2k,n-k)^2 C(n+k,2k) C(2n,n+k)`; `C(n+k,2k)^3 C(2n,n+k)` | `a(p)` mod `p` only | 7AB |

Here `F = Σ_(k<=n) C(2k,k)^3/64^k ≡ a_(4C)(p)` (CITED). The identifications are exact to `p <= 300`
([`halfperiod_twin_class_forms_20261018.py`](../../04-computation/experiments/halfperiod_twin_class_forms_20261018.py),
library of all 21 M24 cycle-shape eta-products plus four CM forms, with characters and `2^(in) 3^(jn)` scales).

## 2. Mechanisms

At `n = (p-1)/2` every library binomial has a `p`-adic shadow (THM-4685, Lemma 1):
- `C(n,k) ≡ (-1)^k a_k (1 - pO_k)`;
- `C(n+k,k) ≡ a_k (1 + pO_k)`;
- `C(n+k,2k) ≡ (-1)^k a_k 4^(-k)`;
- `C(2n,j) ≡ (-1)^j (1 - pH_j)`;
- `C(2n-j,j) ≡ (-1)^j C(2j,j)(1 - p(H_(2j) - H_j))`.

A twin relation is therefore a harmonic congruence modulo `p`:
- **Classes 1, 2, 4:** the Dixon lemma `Σ a_k^3 O_k ≡ 0` from THM-4670. **PROVED** in THM-4685, together with the
  values `F`, `4^n F` and `16^n F` with explicit signs.
- **Class 5:** Clausen's identity (THM-4669).
- **Classes 3 and 7:** both reduce to `(H): Σ_(k<=n) C(2k,k)^3 (H_(2k) - H_k) ≡ 0 (mod p)`. This is KLMSY
  Corollary 4 ("for all primes `p > 3`"), so both classes are PROVED. It was checked here for all `p < 20000`.
- **Class 6:** reduces to a `4F3`-type congruence `(K)`. Holds for all `p < 20000`.

## 3. The CM harmonic law (= KLMSY 2016 Corollary 3)

**Reflection.** Substituting `k -> n-k` (THM-4685, Lemma 3) relates the sums at `z` and `1/z`, where `z = 64/m` and
`m` is the denominator in `C(2k,k)^3 m^(-k)`. At `m = 1`, `(H)` says that the harmonic sum vanishes. Reflected to
`m = 4096`, it says that the harmonic sum is `2 q_p(2)` times `F_4096`. Both are the cases `m = 1, 4096` of the guess
```
(L_p)   6 Σ_(k<=n) C(2k,k)^3 m^(-k) (H_(2k) - H_k) ≡ q_p(m) Σ_(k<=n) C(2k,k)^3 m^(-k)   (mod p).
```

**Prediction.** It was written in the script docstring but committed together with the output (8056a055ae), so the
order is not verifiable. It is also a published theorem: KLMSY Corollary 3, for all CM `λ ≠ 1`.

**Outcome: holds**
([`cm_harmonic_log_law_20261018.py`](../../04-computation/experiments/cm_harmonic_log_law_20261018.py)).
- `m = 1, -8, 16, -64, 64, 256, -512, 4096`: 0 failures in 301 primes each.
- `m = 3, 5, 7, 10, 100, -3, 2, -1, 32, 128`: 287-296 failures each.
- Extended run: the three independent open cases `m = 1, -8, 16` hold for all 2260 primes up to `20000`.

**Structure** (THM-4685; the non-self-dual points are KLMSY Corollary 3, re-derived in THM-4686):
- The law at `m` is equivalent to the law at `4096/m`. This pairs `{1, 4096}`, `{-8, -512}` and `{16, 256}`, which
  are exactly the 2-isogenous CM pairs: discriminants `-7/-28`, `-4/-16` and `-3/-12`.
- At the self-dual points `m = 64` (`j = 1728`) and `m = -64` (`j = 8000`) the law follows from the reflection,
  plus THM-4670 for `m = 64` when `n` is odd.
- Twin form: `T_m = Σ_j C(2j,j) C(p-1-j,j)^2 m^(-j) ≡ (1 - p q_p(m)/3) F_m (mod p^2)`. This is checked at all eight
  points for `p <= 200`. At `m = ±64` it reads `4^(p-1) T_(±64) ≡ F_(±64)`.

**Governing forms per CM point** (`j = (m-16)^3/m`;
[`cm_points_forms_legendre_20261018.py`](../../04-computation/experiments/cm_points_forms_legendre_20261018.py)):

| `m` | form | discriminant |
|---|---|---|
| `1` | 7AB | `-7` |
| `4096` | `(-1/p)` 7AB | `-28` |
| `-8`, `64` | 4C | `-4` |
| `-512` | `(2/p)` 4C | `-16` |
| `-64` | `(-1/p)` 8A | `-8` |
| `16` | `η(2τ)^3 η(6τ)^3` | `-3` |
| `256` | `(-1/p) η(2τ)^3 η(6τ)^3` | `-12` |

**Correction.** `η(2τ)^3 η(6τ)^3` is not an M24 cycle shape. The S43 remark that only M24 classes govern was a
property of the integer-sum library: the `m`-weighted CM family needs the `Q(sqrt-3)` form outside M24.

**Legendre squares** (CITED).
- The identity `P_N(x)^2 = Σ_k C(N+k,2k) C(2k,k)^2 ((x^2-1)/4)^k` is Z.-H. Sun's (1.7) (arXiv:1012.3898).
- The congruence `F_m ≡ P_n(x)^2 (mod p^2)`, with `x^2 = 1 - 64/m`, is in KLMSY Section 3, the route they use with
  Coster-van Hamme.
- Checked here for `N <= 25` and at all eight points.

**Archimedean analogue: REFUTED.** `Σ_k C(2k,k)^3 m^(-k) (H_(2k) - H_k) / Σ_k C(2k,k)^3 m^(-k)` is
`0.0010, 0.0173, -0.0075` at `m = 4096, 256, -512`, against `log|m|/6` of about `1`. The law is genuinely `p`-adic.

**Lucas-block form** (CITED: KLMSY Lemma 10 and Section 4;
[`cm_law_lucas_block_form_20261018.py`](../../04-computation/experiments/cm_law_lucas_block_form_20261018.py)).
`C(2(ap+j), ap+j) ≡ C(2a,a) C(2j,j)(1 + 2ap(H_(2j) - H_j)) (mod p^2)`, and the harmonic weight `6(H_(2j) - H_j)` is
exactly the `p`-shift of `C(2j,j)^3`. So the law is equivalent to the second-order Lucas factorization of the first
block:
```
Σ_(k=p)^(p+n) C(2k,k)^3 m^(-k) ≡ (8/m) Σ_(k=0)^(n) C(2k,k)^3 m^(-k)   (mod p^2).
```
This forces `F_(p^2) ≡ F_p^2 (mod p^2)`. Both hold at all CM points and fail generically. The converse holds only at
ordinary primes (S43b m9).

In the Beukers-Vlasenko picture this is the `s = 2` shadow of a mod-`p^(2s)` Dwork congruence at fixed points of the
excellent lift (their Conjecture 7.5). Their hypercubic `n = 3` family is exactly `3F2(1/2,1/2,1/2;1,1;64t^2)`, and
their Corollary 6.4 has a `log(t^σ/t^p)` correction matching `q_p(m)`. "CM points = fixed points" for `n = 3` is this
repository's heuristic (S43b m5).

**A `p`-adic CM detector** (OBSERVATION; script and output committed together;
[`cm_detector_lucas_block_scan_20261018.py`](../../04-computation/experiments/cm_detector_lucas_block_scan_20261018.py)).
The test is the block congruence `Σ_(k=p)^(2p-1) t_k ≡ t_1 Σ_(k<p) t_k (mod p^2)`, with `t_k = A(k) m^(-k)`, over
all `0 < |m| <= 10^5` and four families:
- `A = C(2k,k)^3`: the survivors are exactly Sun's eight CM values.
- `C(2k,k)^2 C(3k,k)`: 10 survivors.
- `C(2k,k)^2 C(4k,2k)`: 10 survivors.
- `C(2k,k) C(3k,k) C(6k,3k)`: the survivors `{-32768, -3375, 1728, 8000, 54000}`, exactly the rational CM
  `j`-invariants in range (this is the `j`-line family `3F2(1/6,1/2,5/6;1,1;1728/j)`).

All 33 survivors are supersingular at density `0.50-0.54`, measured modulo `p` independently of the test; the controls
lie at `0.01-0.09`.
- The S43b referee reproduced the lists by another method, with confirmation primes up to 600.
- The lists are complete in range relative to the CM signature.
- Every survivor is CM: `τ = j^(-1)(j)` is imaginary quadratic.

This is HYP-9310(c). In every family the detector congruence is equivalent to the harmonic law (S43b m3).

## 4. Literature

**KLMSY 2016** (arXiv:1210.4489v3; read in full by the S43b referee):
- Theorem 2: the mod-`p^2` CM supercongruences, with character `((1-λ)/p)`; v3 misprints it as `((λ-1)/p)`.
- Corollary 3: the CM harmonic law.
- Corollary 4: `(H)`.
- Lemma 10 and Section 4: the block form.
- Conjecture 9: mod `p^(2s)`. The case `s = 1`, every odd `m`, follows from Lemma 10 + Corollary 3 + Theorem 2, and
  also from THM-4686; KLMSY do not state it, and no later proof was found.

**Long-Ramakrishna** (Adv. Math. 2016, arXiv:1403.5232), Theorem 3: the `λ = 1` case mod `p^3`.

**Z.-W. Sun**, *Series with summands involving harmonic numbers* (arXiv:2210.07238), was read in this session.
- **Conjecture 29(i)** is a template for the Ramanujan-type, `k`-weighted sums
  `Σ C(2k,k)^3 m^(-k) [6(ak+b)(H_(2k) - H_k) + a]`, with `log|m|` on the archimedean side and `b(m^(p-1) - 1)` on the
  `p`-adic side.
- **Conjecture 36** (`m = 4096`, `42k + 5`) and **Conjecture 17(iii)** (`m = 1`, `21k + 8`, sum from `k = 1`) are
  instances. 17(iii) was re-checked modulo `p` for `p <= 43`.
- `(L_p)` is the **constant-term** analogue at the CM points, including the divergent point `m = 1`. It is KLMSY
  Corollary 3; Sun does not cite KLMSY.

Other sources:
- **The CM values** `Σ C(2k,k)^3 m^(-k) mod p^2`: Sun's Conjecture 5.2, proved in all cases by KLMSY Theorem 2;
  Z.-H. Sun covers some cases.
- **CDE:** Chowla-Dwork-Evans 1986.
- **The 4C value of `F`:** Ishikawa, van Hamme, Ahlgren.

## 5. Open

- **The class-6 congruence `(K)`** (HYP-9310). Verified for `p <= 10^5`; no reflection forces it; no literature found.
- **The "only if" direction of the CM detector** (HYP-9310(c)).
- **The laws for `C(2k,k)^2 C(3k,k)`, `C(2k,k)^2 C(4k,2k)` and `C(2k,k) C(3k,k) C(6k,3k)`** (HYP-9311). They use
  Sun's weights and are equivalent to the detector congruence. No mod-`p` reflection exists for these families
  (`family_reflection_probe_20261018.out`), so THM-4686's route does not apply. Candidate inputs: the `d = 3, 4, 6`
  analogues of KLMSY Theorem 2 (Mathematics 2013), or Beukers-Vlasenko.
- **The mod-`p^3` block statement at ordinary primes** (HYP-9311(O)), and KLMSY Conjecture 9 at `s >= 2`.
