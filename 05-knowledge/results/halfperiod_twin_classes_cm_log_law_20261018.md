# Half-period twin classes and the CM harmonic log law

**Session:** opus-2026-10-18-S43. **Status:** PROVED (THM-4685) / OPEN, strongly supported (HYP-9310) /
FINITE-EXACT / CITED. Follows THM-4669 and THM-4670 in the
[factorization-vertex-SUSY note](factorization_vertex_supersymmetric_bridge_20261018.md), section 3.

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
- **Classes 3 and 7:** both reduce to `(H): Σ_(k<=n) C(2k,k)^3 (H_(2k) - H_k) ≡ 0 (mod p)`. Holds for all
  `p < 20000`.
- **Class 6:** reduces to a `4F3`-type congruence `(K)`. Holds for all `p < 20000`.

## 3. The CM harmonic log law

**Reflection.** Substituting `k -> n-k` (THM-4685, Lemma 3) relates the sums at `z` and `1/z`, where `z = 64/m` and
`m` is the denominator in `C(2k,k)^3 m^(-k)`. At `m = 1`, `(H)` says that the harmonic sum vanishes. Reflected to
`m = 4096`, it says that the harmonic sum is `2 q_p(2)` times `F_4096`. Both are the cases `m = 1, 4096` of the guess
```
(L_p)   6 Σ_(k<=n) C(2k,k)^3 m^(-k) (H_(2k) - H_k) ≡ q_p(m) Σ_(k<=n) C(2k,k)^3 m^(-k)   (mod p).
```

**Pre-registration.** Written in the script docstring before the run: `(L_p)` holds at all eight of Sun's CM values
and fails at generic `m`.

**Outcome: PASSED**
([`cm_harmonic_log_law_20261018.py`](../../04-computation/experiments/cm_harmonic_log_law_20261018.py)).
- `m = 1, -8, 16, -64, 64, 256, -512, 4096`: 0 failures in 301 primes each.
- `m = 3, 5, 7, 10, 100, -3, 2, -1, 32, 128`: 287-296 failures each.
- Extended run: the three independent open cases `m = 1, -8, 16` hold for all 2260 primes up to `20000`.

**Structure** (THM-4685, all PROVED except the open cases):
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

**Legendre squares.** The exact identity `P_N(x)^2 = Σ_k C(N+k,2k) C(2k,k)^2 ((x^2-1)/4)^k` (checked for `N <= 25`)
and the shadow of `C(n+k,2k)` give `F_m ≡ P_n(x)^2 (mod p^2)` with `x^2 = 1 - 64/m`. This was checked at all eight
points. It is the standard bridge to the Coster-van Hamme CM theory, and is a candidate route to HYP-9310.

**Archimedean analogue: REFUTED.** `Σ_k C(2k,k)^3 m^(-k) (H_(2k) - H_k) / Σ_k C(2k,k)^3 m^(-k)` is
`0.0010, 0.0173, -0.0075` at `m = 4096, 256, -512`, against `log|m|/6` of about `1`. The law is genuinely `p`-adic.

**Lucas-block form** (PROVED equivalent;
[`cm_law_lucas_block_form_20261018.py`](../../04-computation/experiments/cm_law_lucas_block_form_20261018.py)).
`C(2(ap+j), ap+j) ≡ C(2a,a) C(2j,j)(1 + 2ap(H_(2j) - H_j)) (mod p^2)`, and the harmonic weight `6(H_(2j) - H_j)` is
exactly the `p`-shift of `C(2j,j)^3`. So the law is equivalent to the second-order Lucas factorization of the first
block:
```
Σ_(k=p)^(p+n) C(2k,k)^3 m^(-k) ≡ (8/m) Σ_(k=0)^(n) C(2k,k)^3 m^(-k)   (mod p^2).
```
This forces `F_(p^2) ≡ F_p^2 (mod p^2)`. Both hold at all CM points and fail generically. This is the `s = 2` shadow of
Beukers-Vlasenko's mod-`p^(2s)` Dwork congruence at fixed points of the excellent Frobenius lift (their Conjecture 7.5;
the CM points are the fixed points), whose hypercubic `n = 3` family is exactly `3F2(1/2,1/2,1/2;1,1;64t^2)`. Their
Corollary 6.4 has a `log(t^σ/t^p)` correction, which matches the `q_p(m)` here.

**A `p`-adic CM detector** (pre-registered, PASSED;
[`cm_detector_lucas_block_scan_20261018.py`](../../04-computation/experiments/cm_detector_lucas_block_scan_20261018.py)).
The test is the block congruence `Σ_(k=p)^(2p-1) t_k ≡ t_1 Σ_(k<p) t_k (mod p^2)`, with `t_k = A(k) m^(-k)`, over
all `0 < |m| <= 10^5` and four families:
- `A = C(2k,k)^3`: the survivors are exactly Sun's eight CM values.
- `C(2k,k)^2 C(3k,k)`: 10 survivors.
- `C(2k,k)^2 C(4k,2k)`: 10 survivors.
- `C(2k,k) C(3k,k) C(6k,3k)`: the survivors `{-32768, -3375, 1728, 8000, 54000}`, exactly the rational CM
  `j`-invariants in range (this is the `j`-line family `3F2(1/6,1/2,5/6;1,1;1728/j)`).

All 33 survivors are supersingular at density `0.50-0.54`, measured modulo `p` independently of the test; the controls
lie at `0.01-0.09`. The mod-`p^2` Lucas factorization therefore detects CM. This is HYP-9310(c).

## 4. Literature

Z.-W. Sun, *Series with summands involving harmonic numbers* (arXiv:2210.07238), was read in this session.
- **Conjecture 29(i)** is a template for the Ramanujan-type, `k`-weighted sums
  `Σ C(2k,k)^3 m^(-k) [6(ak+b)(H_(2k) - H_k) + a]`, with `log|m|` on the archimedean side and `b(m^(p-1) - 1)` on the
  `p`-adic side.
- **Conjecture 36** (`m = 4096`, `42k + 5`) and **Conjecture 17(iii)** (`m = 1`, `21k + 8`, sum from `k = 1`) are
  instances. 17(iii) was re-checked modulo `p` for `p <= 43`.
- `(L_p)` is the **constant-term** analogue at the CM points, including the divergent point `m = 1`. No statement of
  `(H)` or `(L_p)` was found there.

Other sources:
- **The CM values** `Σ C(2k,k)^3 m^(-k) mod p^2`: Sun's conjectures, proved by Z.-H. Sun and by
  Kibelbek-Long-Moss-Sheller-Yuan (arXiv:1210.4489). Attribution is from a secondary source.
- **CDE:** Chowla-Dwork-Evans 1986.
- **The 4C value of `F`:** Ishikawa, van Hamme, Ahlgren.

## 5. Open

- **HYP-9310:** prove `(L_p)` at `m = 1`, `-8` and `16`, and prove `(K)`.
  - Candidate routes: (a) a Legendre form of `T_m` plus Coster-van Hamme; (b) a WZ pair with a parameter, perhaps
    Guillera-Zudilin's divergent `Σ (21k+8) C(2k,k)^3`; (c) Gauss/Jacobi sums at the CM points.
- **Analogous laws** for `C(2k,k)^2 C(3k,k)`, `C(2k,k)^2 C(4k,2k)` and `C(2k,k) C(3k,k) C(6k,3k)` at their CM points.
  Sun's 29(ii)-(iv) suggest the harmonic weights `3H_(3k) + 2H_(2k) - 5H_k`, `4(H_(4k) - H_k)` and
  `3(2H_(6k) - H_(3k) - H_k)`.
