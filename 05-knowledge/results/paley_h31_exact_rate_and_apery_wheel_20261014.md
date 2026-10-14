# H(P_31) exactly, the Paley rate c₁ = 3/2, and where the Riemann hypothesis enters Paley tournaments

> **2026-10-14 extension and scope update.** The [wheel packet](wheel_apery_extension_20261014h.md) completes the all-odd-prime identity using Ahlgren--Ono and the collision correction. The [run packet](wheel_run_interface_20261014h.md) proves A4/A6/A8, identifies the exact degree-eight modular term, and supplies an independent upper bound making the reported H31 CRT reconstruction unique. The full c1=3/2 asymptotic remains conditional on a uniform run-tail bound; its first-order derivation and finite numerical support below are retained. The q31 DP itself is inherited, not rerun by this extension.

**Contents.**
- `H(P_31)` exactly: the June pre-registered prediction HYP-2371 is confirmed.
- The rate of THM-438's `R(q) → e`: `c₁ = 3/2`, derived and confirmed.
- Exact run polynomials, and pattern sums that vanish for `q ≡ 3 (mod 4)` (THM-4652).
- The wheel `W_4` is the Apéry modular form `η(2z)⁴η(4z)⁴` (HYP-9274).
- The least non-residue as a transitive run, the GRH channel.
- The Giuga/Korselt half of the brief (cross-referenced).

2026-10-14, session opus-2026-10-14-S34 (worktree `/tmp/math-wt-giuga-korselt-plusminus`).

Owner's prompt (verbatim): "think how guiga numbers and korselts criterion relate with guigas conjecture, which revolves around a single value which is both, see how this duality relates deeply with twin primes and even vs odd goldbach. again think in terms of linked lists and a number sometimes embodying the qualities of one greater or less than it. think about prior work on the multiplication/division even/odd half/whole and rational/irrational dualities (as in dedekind cuts) and think about the real/transcendental duality and prior work regarding irrationality exponents and consider elegant solutions to the continuum hypothesis, seeing how syncing together with prior work and bott periodicity and maybe an assumed partial or extended riemann hypothesis we can prove collatz or hadwiger nelson or the planar jacobian by seeing how number theory permeates through the structure of everything and discrete graphs and geometry/topology/analysis/abstract algebras just embody those structures through their lenses in a predictable way, think prior equidistribution, equidecomposability etc related work connects get creative, be free in your reasoning"

**Status.**
- **FINITE-EXACT:**
  - `H(P_27)` and `H(P_31)` (new values);
  - the exact run polynomials `A_4`, `A_6` and the two-cherry ratio;
  - the vanishing of the `K_4` and `K_5` sums for `q ≡ 3 (mod 4)`, and the `K_4` CM values for `q ≡ 1 (mod 4)`;
  - the `W_4` = Apéry identity at 27 primes;
  - the transitive-run table.
- **PROVED:** THM-4652 (pattern-sum vanishing and series-parallel polynomiality); the single-arc orbit reduction (S33).
- **DERIVED (first order):** `c₁ = 3/2`.
- **CONFIRMED:** HYP-2371 (pre-registered June 2026).
- **OPEN:** HYP-9274 (A)–(D).
- **CITED:** THM-438 and its addenda; Ahlgren–Ono (2000); Beukers (1987); Deligne (1974); Ankeny (1952); Bach (1990); Burgess (1957); Graham–Ringrose (1990); Erdős–Moser (1964).
- **Not yet audited** (section 8).

## 0. Answers in brief

1. **Giuga, Korselt, twin primes, Goldbach, CH, Bott: already covered today.** Mac-mini's [synthesis](d19_parity_shift_synthesis_20261014.md) (THM-4650/4651) and codex's packets cover this half of the brief:
   - [signed Giuga defects](giuga_korselt_defects_20261014f.md), [block lattice](giuga_block_lattice_20261014f.md), [successor guards](collatz_successor_guards_20261014f.md) and [finite-witness absoluteness](continuum_finite_witness_boundary_20261014f.md);
   - the "single value which is both" (`n ≡ p (mod p²(p−1))`) and the even-Giuga/odd-Carmichael parity divide;
   - the arithmetic-derivative linked list (`n' = n ± 1`);
   - twin primes flanking primary pseudoperfect numbers (`6 ↦ 30, 42`; `42 ↦ 1722, 1806`; `47058 ↦ 2214408306, 2214502422`, Sondow–Tsukerman);
   - the Goldbach parity barrier;
   - CH cannot decide finite witnesses;
   - Bott via `χ_8` and Milgram.

   This note does not repeat them. It takes up the one clause they leave least developed: **"an assumed partial or extended Riemann hypothesis … number theory permeates … discrete graphs … in a predictable way".** It answers it with the one fully proved channel in our corpus.
2. **Where an RH that is a theorem acts on tournaments.**
   - The Paley tournament's Hamiltonian-path count is controlled by quadratic character sums. THM-438 proved `R(q) = H(P_q)·2^(q−1)/q! → e` from them.
   - Here: the first exact `H(P_31)`, the exact first-order rate `c₁ = 3/2`, and the place where a genuine Hasse–Weil quantity (Deligne's Riemann hypothesis) first enters.
   - The low-order pattern sums are exact polynomials (Gauss and Jacobi sums only), and the complete patterns vanish because a tournament needs `χ(−1) = −1`.
   - The first surviving non-series-parallel pattern, the wheel `W_4`, carries exactly the coefficients of `η(2z)⁴η(4z)⁴`. That is the weight-4 form of Apéry's irrationality proof of `ζ(3)` (Beukers; Ahlgren–Ono).
   - This is literally "irrationality (`ζ(3)`) ↔ RH (Deligne) ↔ tournaments (Paley)".
3. **The GRH channel and the "linked list".** In `P_p` the consecutive residues `0, 1, …, n(p)−1`, up to the least non-residue `n(p)`, form a transitive subtournament.
   - Under GRH this run has length `< 2 log² p` (Ankeny; Bach).
   - Graham–Ringrose: `n(p) ≫ log p · log log log p` infinitely often. Whether their construction allows the residue class `p ≡ 3 (mod 4)` has not been checked.
   - Hence `trans(P_p)/log p` would be unbounded along those primes, and Paley tournaments are not Erdős–Moser-extremal there.
   - Section 5 has the table to `p = 283`.
4. **Collatz, Hadwiger–Nelson, JC2.** No assumed RH is known to imply any of them; mac-mini's §6 states this and we agree. What this note supplies is the precise mechanism in the one place where RH does act on our discrete objects, together with its order of magnitude.

## 1. `H(P_31)` exactly (FINITE-EXACT)

**Method.**
- `H(P_q)/|Aff| = N_1`, the number of Hamiltonian paths starting with the arc `0 → 1` (S33, section 5c).
- `N_1` is computed by a layered colex-rank subset DP: the S33 audit referee's layout, two layers in memory, about 9 GB.
- The DP runs in `uint32` residues modulo the primes `2147483647`, `2147483629` and `2147483587`, and the exact value is recovered by CRT. The modulus product is about `2^93`, and `N_1 ≈ 2^75`.

**Validation chain** (the compute-design doc's acceptance list):
- the same code reproduces `H(P_7), H(P_11), H(P_19), H(P_23)` exactly (`H(19) = 1172695746915`, `H(23) = 15760206976379349`);
- at `q = 27, 31` it matches the independent S33 and audit residues mod 256 (1 and 41);
- for `P_27` it also matches `H ≡ 3159 (mod 8424)` and `H/1053 ≡ 3 (mod 8)` (S33).

| `q` | `H(P_q)` | `R(q) = H 2^(q−1)/q!` | `q(1 − R/e)` |
|---|---|---|---|
| 7 | 189 | 2.400000000 | 0.8196 |
| 11 | 95095 | 2.439506173 | 1.1281 |
| 19 | 1172695746915 | 2.527147849 | 1.3360 |
| 23 | 15760206976379349 | 2.556980364 | 1.3648 |
| **27** | **418425945688463432031** | 2.578788368 | 1.3856 |
| **31** | **19876729650849384218015865** | **2.595511945** | 1.4001 |

**HYP-2371 (kind-pasteur, 2026-06-10, pre-registered): CONFIRMED.**
- The prediction was `R(31) ∈ [2.58949, 2.60249]`, central 2.59599. The exact value is `R(31) = 2.5955119`, which is `−0.00048` from the centre (half-width `0.0065`).
- Acceptance tests: `H` odd; `465 | H`; `H/465 = N_1 = 42745655163116955307561` odd (HYP-1713); `H ≡ 465 (mod 930)`. All pass.
- `H(P_31)` is a lower bound for OEIS A038375(31), which is known exactly only to `n = 11`.

## 2. The rate: `c₁ = 3/2` (DERIVED + NUMERICAL)

**Fits** (`paley_R_asymptotic_fit_20261014.out`):
- With `1 − R/e = Σ c_j/q^j`, order-2 and order-3 fits on the largest primes give `c₁ = 1.50156` and `1.50121`.
- The order-2 fit through `q = 19, 23` predicted `R(31)` with error `−5·10⁻⁶`.
- The prime power `q = 27`, never fitted, lies on the prime curve to `5·10⁻⁵` (universality).

**Derivation.** HYP-9274 (A) gives the details.
- `R(q) = E_σ Π_k (1 + χ(d_k))` over random orderings `σ`.
- Steps that are isolated, or odd in total, contribute 0 (THM-4652).
- A cherry weighs exactly `1/(q−2)`.
- Placing `j` separated cherries gives `C(q−2j, j)/(q−2)^j = (1/j!)(1 − (5/2) j(j−1)/q + …)`.
- Two cherries repel by the exact factor `(q−2)(q−9)/((q−4)(q−5)) = 1 − 2/q + …`.
- A quartic run weighs `2/((q−2)(q−4))`, for a total of `2/q` over all positions.
- Hence `R = e(1 − 5/(2q) − 1/q + 2/q) + O(q^(−2)) = e(1 − 3/(2q)) + O(q^(−2))`.

**With `c₁ = 3/2` fixed.**
- `q²(1 − R/e − 3/(2q))` is nearly constant: `−3.117, −3.109, −3.090, −3.097` at `q = 19, 23, 27, 31`. So `c₂ ≈ −3.08`.
- The holdout error at 31 is `1·10⁻⁵`.
- THM-438's forward note expected "`C ≈ 1.4`". That is the effective one-term value at `q = 31` (`1.4001`); the asymptotic coefficient is `3/2`.

**Pre-registered** (HYP-9274 (A)):
- `R(43) ∈ [2.62798, 2.62803]` and `R(47) ∈ [2.63531, 2.63535]` with `c₁ = 3/2` fixed;
- `[2.62729, 2.62868]` and `[2.63468, 2.63594]` model-agnostic.

## 3. Exact pattern sums, and why tournaments stay polynomial so long (THM-4652, FINITE-EXACT)

**Run sums** over `L + 1` distinct points:
- `A_2 = q(q−1)`;
- `A_4 = 2q(q−1)(q−3)` (11 primes);
- `A_6 = q(q−1)(q−3)(5q−29)` (6 primes; exact under leave-one-out interpolation);
- odd `A_L = 0`;
- leading coefficients `1, 2, 5` (Catalan, THM-438), now as exact polynomials.

**Why** (THM-4652):
- An isolated edge, an odd edge count, or an odd-reversal relabelling forces `0`.
- Series-parallel patterns are polynomials in `q`.
- The complete patterns are `±χ(Vandermonde)` and vanish when `χ(−1) = −1`, which is exactly when Paley tournaments exist.

**Contrast with Paley graphs.**
- For `q ≡ 1 (mod 4)` the same `K_4` sum equals `4a² − 2q` (`q = a² + b²`, `a` odd) for all 21 primes below 200. That is the CM-by-`Q(i)` arithmetic behind the classical 4-clique counts of Paley graphs.
- For `q ≡ 3 (mod 4)` it is 0 (24 primes); the `K_5` sum is 0 too (19 primes).
- Paley tournaments live exactly where this arithmetic is switched off.

## 4. The wheel `W_4` is the Apéry form (HYP-9274 (B), FINITE-EXACT at 27 primes)

**Why `W_4`.** The first 5-vertex pattern that survives THM-4652 and is not series-parallel is `K_5` minus two disjoint edges, the wheel `W_4`. Its sum is not a polynomial: cubic interpolation fails by about `10^5` at `q = 107`.

**The identity.** At every prime `5 <= q <= 109` (both classes mod 4),

`Σ*(W_4) = −q(q−1)(a_q + 4)`, where `f = η(2z)⁴η(4z)⁴ = Σ a_n q^n` (the `a_q` are computed in the script from the eta product).

**Context.**
- `f` is the weight-4 newform of level 8 in Beukers' supercongruence for Apéry's numbers. Apéry's numbers prove `ζ(3) ∉ Q`, and Ahlgren–Ono proved the supercongruence via Gaussian `₄F₃(1)`.
- Deligne's theorem gives `|a_q| <= 2q^(3/2)`; the observed ratio is at most 0.88.
- Expected proof: the scaled, non-distinct sum is the `χ`-twisted second moment `Σ_b χ(b)·a(E_b)²` of the Legendre family. That is the Ahlgren–Ono object; the distinctness corrections are polynomial.

**What it means for `R(q)`** (HYP-9274 (D), heuristic).
- A single run cannot glue into `W_4`: a trail has at most two odd-degree vertices, and `W_4` has four. Two runs can.
- The `W_4` channel therefore contributes about `a_q(f)/q^6` to `1 − R/e`, at most `2q^(−4.5)`.
- **This is the precise sense in which "an RH" permeates the Hamiltonian-path count:** a Deligne-bounded Fourier coefficient of the `ζ(3)` form, entering at relative order `q^(−4.5)` or lower.
- Other Eulerian non-series-parallel channels (the octahedron, from single runs of 12 steps or more) are not estimated.

## 5. The least non-residue as a transitive run (CITED + FINITE-EXACT)

For a prime `p ≡ 3 (mod 4)`, the set `{0, 1, …, n(p)−1}` induces a transitive subtournament of `P_p`. Only the initial run of residues works: transitivity needs all pairwise differences to be residues.

| statement | status |
|---|---|
| `trans(P_p) >= max(n(p), ⌊log₂ p⌋ + 1)` | PROVED (elementary + Erdős–Moser) |
| GRH ⟹ `n(p) < 2 log² p` | CITED (Ankeny; Bach) |
| `n(p) ≫ log p · log log log p` infinitely often, hence `trans(P_p)/log p` unbounded | CITED (Graham–Ringrose; restriction to `p ≡ 3 (mod 4)` not checked) |
| discrepancy bound alone gives `trans(P_p) <= 2√p + O(1)` | PROVED (elementary) |

**Data** (`paley_transitive_runs_20261014.out`, `p <= 283`, exact branch and bound, brute force for `p <= 19`):
- `trans(P_p)` grows from 2 to 14 (e.g. `trans(P_283) = 14` against `2 log₂ 283 = 16.3`), while `n(p) <= 7`;
- every row satisfies Erdős–Moser.

**Process note.** A first version of the search reused the clique-search trick of discarding explored vertices. That is invalid for ordered chains, and it gave `trans(P_19) = 4 < 5`. The built-in Erdős–Moser check caught it.

## 6. The other clauses of the brief, typed

**Dedekind cuts, rational/irrational, real/transcendental, irrationality exponents.**
- S31 (THM-4632) and mac-mini's THM-4633/4634 carry the irrationality-measure dictionary for Collatz.
- Here the new rhyme is that the Paley `W_4` count is governed by the modular form of Apéry's `ζ(3)` irrationality proof. That is an exact identity (FINITE-EXACT), not an analogy.

**Equidistribution and equidecomposability.**
- Laczkovich's circle squaring needs discrepancy bounds that are Diophantine in nature (mac-mini §6).
- The Paley analogue proved here is that the first-order discrepancy of path statistics is the rational `c₁ = 3/2`, with no arithmetic term until the `W_4` channel.

**Bott periodicity (ANALOGY).** The mod-8 phenomena in this thread are:
- HYP-9273 (B), the orbit counts `≡ 1 (mod 8)`;
- `χ_8`;
- Milgram's Gauss-sum signature.

The octahedron, the first single-run candidate in HYP-9274 (D), is the cross-polytope of dimension 3. A real link to Clifford periodicity is not claimed.

**CH and the Jacobian conjecture.**
- Every statement in this note is arithmetic: finite computations and identities between integers. By codex's finite-witness note, none of them can depend on CH.
- The Jacobian conjecture is untouched.

## 7. Open threads

1. Prove HYP-9274 (B) through the Legendre-family twisted second moment (Ahlgren–Ono `₄F₃`). It should be a short derivation from known identities, with the distinctness corrections as polynomials.
2. Derive `c₂ ≈ −3.08` exactly at second order. The ingredients are three-cherry correlations, cherry–quartic cross terms, the `A_6` leading term `5/q³`, and second-order placement.
3. Settle HYP-9274 (D): find the first order at which `R(q)` stops being rational, through the `W_4` channel or the octahedral channel.
4. `H(P_43)`: test the pre-registered intervals. Subset DP is out of reach (`2^41` subsets), and Karp inclusion–exclusion needs about `4·10^15` operations. A new algorithm is needed.
5. `A_8` exactly. It is `O(q^8)` per prime, so feasible to `q ≈ 31`. Test the scoped polynomiality conjecture.

## 8. Audit

Pending: an independent adversarial audit will be recorded here.

## 9. Reproduction and references

**Reproduction:**
- `paley_H31_exact_20261014.py` (about 3 minutes, 9 GB);
- `paley_R_asymptotic_fit_20261014.py`;
- `paley_rate_c1_ingredients_20261014.py` (about 10 minutes);
- `paley_A6_polynomiality_20261014.py` (about 5 minutes);
- `paley_W4_charsum_20261014.py` and `paley_W4_apery_form_20261014.py` (about 10 minutes);
- `paley_K4_charsum_cm_20261014.py`, `paley_K5_charsum_20261014.py`, `paley_transitive_runs_20261014.py 300`.

All are in `04-computation/experiments/`, with `.out` files beside this note.

**References.**
- **This repository:**
  - THM-438 (monad-explorer; `R(p) → e`, Catalan cluster integrals);
  - HYP-2371 and the [compute design](paley_H31_compute_design_kpc1.md) (kind-pasteur);
  - S33 [note](number_tournament_bridges_20261014.md) (`N_1` reduction, HYP-9273);
  - THM-582;
  - mac-mini's [D19/Giuga synthesis](d19_parity_shift_synthesis_20261014.md);
  - codex's Giuga/continuum packets (linked in section 0).
- **Literature:**
  - S. Ahlgren and K. Ono, A Gaussian hypergeometric series evaluation and Apéry number congruences, *J. reine angew. Math.* 518 (2000) 187–212;
  - F. Beukers, Another congruence for the Apéry numbers, *J. Number Theory* 25 (1987) 201–210;
  - P. Deligne, La conjecture de Weil I, *Publ. Math. IHÉS* 43 (1974);
  - N. C. Ankeny (1952); E. Bach, *Math. Comp.* 55 (1990);
  - D. A. Burgess (1957);
  - S. W. Graham and C. J. Ringrose (1990), lower bounds for least quadratic non-residues;
  - P. Erdős and L. Moser (1964), transitive subtournaments;
  - OEIS A038375.
