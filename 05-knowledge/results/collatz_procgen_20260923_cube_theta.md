# The cube-swap number (HYP-9127): where it sits, a formal-Padé reduction, and why every known mechanism stops

**Status.**
* **OPEN:** HYP-9127. Is `X = sum_(k>=0) rho^(k^3)`, with `rho = 2^10/3^9`, irrational in `Q_2`? Equivalently, does no rational have an eventually-`Y3` parity vector under `3x+1`? Not proved. Progress here is **partial**: a proved quadratic tier, proved obstructions, two proved conditional reductions and finite certificates. The real analogue `sum (2^10/3^9)^(k^3)` in `R` is also OPEN; no result on it was found.
* **AUDIT (orchestrator, 2026-09-23).** Re-derived by hand, with no gap found:
  * **Theorem Q.** Substituting `q = rho^(-2 alpha)` and `z = eps rho^(alpha+beta)` gives exactly Tschakaloff's series. Lemma P's bookkeeping rescales by `a/2`, which leaves the same bracket and the same `phi` threshold.
  * **Theorem U.** Set `D_K = b A_K - a 3^(9K^3)`. Then `v_2(D_K) = 10(K+1)^3` and `|D_K| <= (K+2) C H`, so the S-free part is `<= H^(1-1/mu_bar+o(1)) = H^(0.299)`. No subsum vanishes, by the 2-adic dominance of the least omitted term. `U(1/3)` itself is unproved; Theorem U is a reduction, not a proof.
* **PROVED** (this note; elementary given the cited inputs; exact checks in the `.out`):
  * **Theorem Q (quadratic tier).** Let `P` be any integer-valued quadratic with positive leading coefficient, `eps = +-1`, and `rho = 2^L/M_0` (`M_0` odd) with `mu_bar < phi`. Then `sum_(k>=k0) eps^k rho^(P(k))` is irrational in `Q_2`. So swap words whose swap positions are the values of one quadratic polynomial are settled under every `3x+r`, and under `5x+1` or `7x+1` whenever `mu_bar < phi`. Examples: triangular numbers, either branch of the pentagonal numbers. This proves HC6 of the hard-class note for single quadratic families and extends Theorem Y.
  * **Proposition Lv (Liouville exactness).** For positions given by a polynomial of degree `>= 2`, partial sums certify irrationality iff `mu_bar = 1`. For the maps `x/2, (mx+r)/2` this is exactly the subcritical case, which Lemma L already excludes. So Liouville settles no supercritical polynomial swap word.
  * **Proposition H.** Cubic positions admit no linear q-difference equation of any order, in one variable, for any rational power `q` of a formal `rho`. The proof goes through Garoufalidis's quadratic-degree theorem, and strengthens the hard-class note's "no first-order equation".
  * **Proposition M.** The functional equation of `G` is a Mahler equation with a unipotent matrix. By Kubota (as recorded by Adamczewski–Faverjon), no point is admissible, so Mahler's method in several variables does not apply.
  * **Proposition F and Corollary F3 (formal reduction).**
    * Every approximation argument that treats `rho` as a formal variable produces integer pairs `(A, B)` for the 0/1 cube series `Theta_3(x) = sum x^(k^3)`.
    * A Newton-polygon lemma certifies non-vanishing whenever the leading remainder coefficient is not divisible by `2^L`.
    * Consequently HYP-9127 follows from a purely combinatorial **2-adic zero estimate** (HC-CT1) for the lattice of such pairs. Siegel's lemma already supplies the gap with polynomial heights.
  * **Proposition S and Theorem U.**
    * No shift nearly preserves the cubes, so the fixed-dimension "return" compression of Erdős 1062(ii) is unavailable, and a Subspace argument needs about `(log H)^(1/3)` S-unit coordinates.
    * The uniform S-unit gap hypothesis `U(1/3)` implies HYP-9127.
    * The Browkin–Brzeziński n-conjecture does not.
* **FINITE-EXACT:**
  * The `Y3` block identity is re-verified with independent code modulo `2^20000` (Bernstein closed form plus parity iteration).
  * **No rational of small height.** `X` is not a rational `a/b` with `|a|, |b| <= 2^4999999`: PARI's `bestappr` returns no fraction at `10^7` bits, and Wang reconstruction in Python returns none at `10^5` bits.
  * **No integer relation of small norm** (Gram–Schmidt certificates at 40000 bits). None has Euclidean norm below:
    * `2^5714` among `1, X, ..., X^6`;
    * `2^6666` among `1, X, Y_1, ..., Y_4`, where `Y_j = sum k^j rho^(k^3)`;
    * `2^9999` among `1, X, Theta_2(rho), X Theta_2(rho)`;
    * `2^7999` among `1, X` and the values of `G` at the three non-tail points.

    Achieved exponents match the Dirichlet baseline to within 0.03%, and the algebraic control is found exactly.
  * **Formal Padé certificates beat the Liouville threshold at every tested size.** The diagonal `[n/n]` forms of `Theta_3` at `rho` have margin 5297 bits at `n = 1000`, and Siegel forms of height `<= 129` (mostly 1) exist for all tested `E <= 300`. In every case the Newton polygon certifies non-vanishing.
  * **Hankel census.** The Hankel determinants of the cube indicator grow like those of random 0/1 sequences with the same counting function, about `0.045 n log2 n` bits, so there is no Padé miracle for cubes. Those of the square indicator are about 4.4 times smaller than random.
  * **Zudilin analogue.** The cubic analogue of Zudilin's construction is worse than the partial sums for `n <= 6`: its exponents fall from 1.40 to 0.33.
  * **Theorem Q's forms** are exact for `n <= 12` on five quadratic families under four maps; the 2-adic valuation equals its prediction in every case.
  * **Repetition profiles.** Those of the square, `floor(k^(3/2))` and cube swap words all decay to 0, like `j^(-0.50)`, `j^(-0.48)` and `j^(-0.33)` in range.
* **CITED.**
  * `[P]` (primary text read): Ghidelli 2019; Bugeaud–Evertse; Zudilin 2005; Adamczewski–Faverjon 2018; the Erdős 1062(ii) exposition.
  * `[P, abstract]`: Bailey–Borwein–Crandall–Pomerance 2004; Garoufalidis; Corvaja–Zannier 2002; Krattenthaler–Rochev–Väänänen–Zudilin 2009.
  * `[R]` (named secondary source): Ridout, Bradshaw, Nesterenko, Bertrand, Duverney, Evertse, Tschakaloff, Bundschuh.
  * Metadata only: Väänänen–Wallisser, Bézivin, Barré-Sirieix–Diaz–Gramain–Philibert, Masser, Browkin–Brzeziński.
* **REFUTED** (as a proposed cheaper instance): swap positions `floor(k^(3/2))`. They share every obstruction of `Y3`: `Dio = 1`, no q-difference equation, and a nearly generic Hankel structure. They also need `N^(2/3)` S-unit terms instead of `N^(1/3)`.

Session `collatz-procgen-20260922`, cube-theta lane (machine `mac-mini`, 2026-09-23). It builds on the [hard-class note](collatz_procgen_20260922_hard_class.md) (section 2.4: Theorem Y, Lemma P, the `Y3` identity) and the [transversality foundry](collatz_procgen_20260922_transversality_foundry.md) (section 4.3, Theorem R and Lemma G). No HYP or THM file was created; candidates are in section 9.

* **Scripts:**
  * [`..._cube_theta_core.py`](../../04-computation/experiments/collatz_procgen_20260923_cube_theta_core.py) (exact 2-adic series, words, Bernstein formula);
  * [`..._cube_theta_ledger.py`](../../04-computation/experiments/collatz_procgen_20260923_cube_theta_ledger.py) (sections L0–L7);
  * [`..._cube_theta_lattice.py`](../../04-computation/experiments/collatz_procgen_20260923_cube_theta_lattice.py) (T1–T3);
  * [`..._cube_theta_pade.py`](../../04-computation/experiments/collatz_procgen_20260923_cube_theta_pade.py) (P1–P4);
  * [`..._cube_theta_run.sh`](../../04-computation/experiments/collatz_procgen_20260923_cube_theta_run.sh).
* **Output:** [`collatz_procgen_20260923_cube_theta.out`](collatz_procgen_20260923_cube_theta.out).
* **Web policy.** Every request used the generic user agent. A search-engine anomaly page, a Springer JavaScript challenge and a Project Euclid bot wall were met and **not** bypassed; the affected items are marked metadata-only or `[R]`. Downloaded texts were kept outside the repository.

## 0. Answers in brief

**Task 1 (literature).** Neither HYP-9127 nor its real analogue is known. Both lie outside every method on record.
* **Unit case, real** (`q = 1/b`). `sum b^(-k^3)` is trivially irrational. Bailey–Borwein–Crandall–Pomerance's density bound gives degree `>= 3`, and Bradshaw's nested gaps give degree `>= l` for `l`-th powers. Ghidelli (2019) proved `deg >= 4` for cubes, using Waring-type gaps. Transcendence is open.
* **Non-unit case, real** (`q = a/b`, `a >= 2`). No irrationality result was found. Liouville, Ridout and nested gaps all need the numerator 1.
* **p-adic.** The unit-like case (`M_0 < 2^L`) is irrational by Liouville (Proposition Lv). The p-adic literature found either covers first-order q-difference (Tschakaloff-type) series or Liouville-type super-fast gaps. The q-difference items are Bézivin 1988/1990 and Väänänen–Wallisser 1991 (metadata only); the fast-gap item is Çalışkan 2022 (abstract). None of them treats polynomial gaps of degree `>= 3`.
* **Mahler's method.** It needs matrices in the class M: no eigenvalue a root of unity, spectral radius `> 1`. `G`'s matrix is unipotent, so the method does not apply.
* **Subspace.**
  * Ridout gives transcendence of `sum b^(-lambda_n)` once `limsup lambda_(n+1)/lambda_n > 1`.
  * Bugeaud–Evertse's quantitative Subspace theorem handles about `(log N)^(3/2)` nonzero digits up to `N`; cubes need `N^(1/3)`.
  * Corvaja–Zannier 2002 treats lacunary values at S-unit points of Hadamard type.
* **Difficulty.** The real non-unit analogue is at least as far out of reach as Tschakaloff's problem would be without its q-difference equation.

**Task 2 (attempts).**
* **(a) Two-variable Padé / Hermite–Padé.**
  * Every construction that is formal in `rho` reduces to integer pairs `(A, B)` for `Theta_3` (Proposition F). This includes two-variable constructions evaluated at the tail points `(rho^(3s^2), rho^(3s))`.
  * Such pairs beat the Liouville threshold at every tested size: the natural height is `2^(14.26 E)`, against a 2-adic valuation of `10 W` with `W ~ 2E`. Non-vanishing is certified each time.
  * No family with provable asymptotics exists: diagonal Padé heights grow like `n log n` (P1), and the Siegel forms lack a 2-adic zero estimate.
  * The explicit two-variable analogue of Zudilin's construction is worse than Liouville (Proposition C).
* **(b) Subspace / S-units.**
  * **The growing dimension is intrinsic** (Proposition S). Every nonzero shift changes at least `2 N^(1/3) - O(1)` cube indicators below `N`.
  * **The bulk integer is costly.** In the 2-adic problem the non-S-unit "bulk" integer must be paid for at the archimedean place. This is the mirror image of Erdős 1062(ii).
  * **The n-conjecture is too weak.**
  * **Only finite partial results exist:**
    * no rational of height `<= 2^4999999`;
    * no algebraic relation of degree `<= 6` below norm `2^5714`;
    * the formal certificates.

    No special-denominator result was found; section 5.3 explains why.
* **(c) Conditional.** `U(1/3)` implies HYP-9127 (Theorem U). So does the 2-adic zero estimate HC-CT1 (Corollary F3). No implication from p-adic Schanuel is claimed: `X` is not an exponential-type value.
* **(d) Lattice hunt.** The exponent ratio `tau/m` is 1.0001 to 1.0003 for every family, against a Dirichlet baseline of 1. The Hadamard-lacunary control dips to 2/3 and the algebraic control is found exactly. So there is no hidden structure at 20000 and 40000 bits.
* **(e) Degree hierarchy (table in section 7).**
  * One quadratic family is PROVED for `mu_bar < phi` (Theorem Q). Zudilin's one-parameter family is optimal at `phi`, but diagonal Padé forms certify beyond `phi` at every tested `n` (up to `mu_bar = 1.86`); no proof.
  * Two quadratic families, such as the bilateral pentagonal numbers, are OPEN here.
  * Cubic and higher positions are OPEN at every supercritical `mu_bar`.
* **(f) Cheaper open instances.**
  * Square swaps with `mu_bar` in `[phi, 2)`. These are theta values with sub-random Hankel structure; the closest to `phi` is `5x+1` with 23 ones in 33 letters (`mu_bar = 1.61831`).
  * Bilateral quadratic families (two Tschakaloff values).
  * The near-critical cube swap `sum (2^19/3^12)^(k^3)` (`mu_bar = 1.001`), the cheapest cubic instance.
  * `floor(k^(3/2))` is not cheaper.

**Task 3 (context).**
* HYP-9127 is the 2-adic cousin of the **irrationality** problem for non-unit cubic theta values. That problem is open in `R`, with no literature found.
* It is not the cousin of the unit transcendence problem for `sum b^(-k^3)`. The 2-adic analogue of that problem sits at `mu_bar = 1`, where irrationality is elementary.
* Squares crossed the corresponding gap in 1919: Tschakaloff used the first-order q-difference equation. For cubes that equation provably does not exist (Proposition H), so HYP-9127 is the first place where the Tschakaloff mechanism is structurally absent.
* Progress would be one of the items of section 8.

## 1. Setting

* **The number.** `rho = 2^L/M_0` with `M_0 >= 1` odd, and `mu_bar = max(1, log2(M_0)/L)`. For `Y3` under `3x+1`, `L = 10`, `M_0 = 3^9` and `mu_bar = 0.9 log2 3 = 1.426466`.
* **The identity.** `Phi_2(Y3) = -1 - 512/(3^9-2^10) - (256/3^9)(X-1)` (hard-class note, section 2.4). It is re-checked in `.out` L1 with an independent Bernstein-formula code at `N = 20000`, and parity iteration reproduces 20000 letters.
* **Real values.** `X_R = 1.052024589801160144536635448506505408046` and `Phi_R(Y3) = -1.028116480848713921941372840810469987048`.
* **The swap word in general.** For blocks `B != B'` of length `L` with equal multiplier `M_0`, and swap set `S`, the hard-class block identity gives

  `Phi_T(Y_S) = -(1/M_0) [R_B/(1-rho) + (R_B' - R_B) Theta_S(rho)]`, with `R_B != R_B'`.

  Hence `Phi_T(Y_S)` is rational iff `Theta_S(rho) = sum_(s in S) rho^s` is rational.
* **Liouville ledger** (`.out` L2).
  * The partial sums `S_K = A_K/3^(9K^3)` satisfy `v_2(X - S_K) = 10(K+1)^3` exactly, against heights `2^(14.26 K^3 + O(1))`.
  * The gain is positive only for `K <= 7`, with a maximum of 376.8 bits at `K = 5`, and negative from `K = 8`.
  * The ratio `v_2/log2 H` tends to `1/mu_bar = 0.701`.

## 2. Literature (task 1)

| object | place | what is known | source |
|---|---|---|---|
| `sum q^(k^2)`, `q` algebraic, `0<|q|<1` | `C` | transcendental | Nesterenko 1996, via Bertrand 1997 Thm 4 and Duverney–Nishioka–Nishioka–Shiokawa 1996 [R via Ghidelli §2] |
| `T_q(z) = sum z^n q^(-n(n-1)/2)`, `q = q1/q2`, `z` rational | `R` | irrational if `log|q2|/log|q1| < (3-sqrt5)/2`, with an irrationality measure | Tschakaloff 1919, Bundschuh [R]; Zudilin 2005 [P] |
| `T_q(z)` and `q`-exponential values, integral `q` | `R` | not quadratic | Bézivin 1998 [R via Ghidelli]; Krattenthaler–Rochev–Väänänen–Zudilin 2009 [P, abstract] |
| `T_q(z)`, `q` in `Q` with `0<|q|_p<1` | `Q_p` | a linear independence measure; hypotheses UNVERIFIED | Väänänen–Wallisser, JNT 39 (1991) 225–236 [metadata] |
| p-adic q-difference solutions | `Q_p` | linear independence of values | Bézivin, Manuscripta 61 (1988) 103–129 and Acta Arith. 55 (1990) 233–240 [metadata] |
| `sum rho^(k^2)`, `rho = 2^L/M_0` | `Q_2` | irrational for `mu_bar < phi` | hard-class Lemma P / Theorem Y; Theorem Q here |
| same, transcendence | `Q_2` | no p-adic Nesterenko found; `j(rho)` is transcendental (Mahler–Manin), which does not bear on `theta(rho)` alone | Barré-Sirieix–Diaz–Gramain–Philibert, Invent. 124 (1996) 1–9 [metadata; statement R] |
| `sum b^(-k^3)`, `b >= 2` integer | `R` | irrational (trivial); not quadratic, by the density bound `#{1-bits <= N} > C N^(1/D)` for degree `D`; `deg >= l` for `l`-th powers (Bradshaw 2013); **`deg >= 4`** for cubes and `>= 5` for fourth powers; transcendence OPEN | BBCP, JTNB 16 (2004) 487–518 [P, abstract]; Bradshaw [R via Ghidelli]; Ghidelli, arXiv 1910.05076 Thm 1.1 [P] |
| `sum b^(-lambda_n)`, `limsup lambda_(n+1)/lambda_n > 1` | `R` | transcendental | Ridout 1957, via Schneider Satz 7 [R via Bugeaud–Evertse] |
| gap series with few digit changes | `R` | an algebraic irrational has `>= c (log n)^(3/2) (log log n)^(-1/2)` digit changes; `sum b^(-2^[j^eta])` is transcendental for `eta > 2/3` | Bugeaud–Evertse, arXiv 0709.1560 Thm 3.1, Cor. 3.2 [P] |
| lacunary series at S-unit points | `C_v` | Subspace-theorem transcendence, including Mahler's classical cases | Corvaja–Zannier, Compositio 131 (2002) 319–340 [P, abstract] |
| lacunary series at `U_m`-numbers | `Q_p` | Liouville-type, needing super-fast gaps | Çalışkan, C. R. Acad. Bulg. Sci. 75 (2022) 477–485 [P, abstract] |
| **`sum (a/b)^(k^3)`, `a >= 2`** | `R` | **no result found**; Liouville, Ridout and nested gaps all need `a = 1` | this note, section 5 |
| **`sum rho^(k^3)`, `mu_bar > 1`** (HYP-9127) | `Q_2` | **OPEN**; irrational for `mu_bar = 1` only (Proposition Lv) | this note |
| Mahler's method in several variables | `C` | admissible iff `T` is in class M, `T^k alpha -> 0` and `alpha` is `T`-independent; a root-of-unity eigenvalue excludes every point (Kubota) | Adamczewski–Faverjon, arXiv 1809.04823, Def. 1.2, Def. 3.1, Thm 3.4, Lemma 3.7 [P]; Masser, Invent. 67 (1982) [metadata] |
| q-holonomic sequences | formal | `deg_q a_n` is eventually a quadratic quasi-polynomial | Garoufalidis, arXiv 1005.4580 [P, abstract] |
| rotation-coded density sums | `R` | irrational via fixed-dimension (`<= 404`) three-place Subspace covering | Erdős 1062(ii) exposition (scratchpad copy) [P] |

**Ghidelli 2019 [P].**
* **Result.** For `l` in `{3, 4}` and an integer `q >= 2`, `theta_l(q) = sum q^(-n^l)` is either transcendental or of degree `>= l + 1` (Thm 1.1).
* **Method.** The proof abstracts Bradshaw's **nested gaps** into a linear-independence criterion (Thm 3.3). If `alpha f(1/q) + beta g(1/q) = 0`, then at a gap point the truncated sum is a rational with denominator `q^(n-1)` smaller than the tail, hence zero.
* **Why the unit numerator matters.** It makes a gap of *bounded* length sufficient.
* **Transfer to `Q_2`.** The 2-adic transfer would need gaps of length `> (mu_bar - 1) n` at position `n`. Sums of `j` cubes have relative gaps tending to 0, so the technique gives nothing for `mu_bar > 1` (section 4.1).
* **"Cubic theta."** The paper's "cubic theta series" means the 0/1 series of cubes. It is not the Borweins' modular cubic theta functions, which are what a keyword search returns.

**Verdict.** HYP-9127 is not known, and neither is its real analogue. The unit real case is known to degree `>= 4` and is open for transcendence. Every irrationality method on record uses either a unit numerator (Liouville, Ridout, nested gaps) or a functional equation (Tschakaloff, Mahler, Nesterenko); HYP-9127 has neither.

## 3. The formal framework (task 2a)

### 3.1 Proposition F (formal certificates; PROVED)

Let `S` be an infinite subset of `Z_(>=0)`, `Theta_S(x) = sum_(s in S) x^s` in `Z[[x]]`, `rho = 2^L/M_0` (`L >= 1`, `M_0` odd), and `xi = Theta_S(rho)` in `Z_2`. Let `A, B` in `Z[x]` have degree `<= E`, with `A != 0`, and write

`R := A Theta_S - B = sum_(j >= W') r_j x^j`, with `r_(W') != 0` and `v_2(r_(W')) < L`.

Put `P = M_0^E A(rho)` and `Q = M_0^E B(rho)`; both are integers. Then:
1. `v_2(P xi - Q) = L W' + v_2(r_(W'))`. In particular `R(rho) != 0`.
2. If `xi = a/b` in lowest terms, then `L W' + v_2(r_(W')) <= log2(|a| + |b|) + log2 max(|P|, |Q|)`.
3. `|P| <= ||A||_1 max(2^L, M_0)^E`, and likewise `|Q|` with `||B||_1 <= ||A||_1 #(S cap [0, E])`.

So the **margin** `L W' + v_2(r_(W')) - log2 max(|P|, |Q|)` bounds `log2(|a|+|b|)` from below for any rational equal to `xi`, and unbounded margins prove `xi` irrational.

*Proof.*
1. `R(rho) = A(rho) xi - B(rho)` (Cauchy product in `Q_2`).
2. For `j > W'`, `v_2(r_j rho^j) = v_2(r_j) + L j >= L(W'+1) > L W' + v_2(r_(W'))`, since `v_2(r_(W')) < L`. So the first term alone determines `v_2(R(rho))`.
3. `P xi - Q = M_0^E R(rho)`, and `M_0` is odd. This gives (1).
4. For (2): `b` is odd, because `xi` is in `Z_2`. Then `N = Pa - Qb = b(P xi - Q)` is a nonzero integer with `v_2(N)` equal to the value in (1), and `2^(v_2(N)) <= |N| <= (|a|+|b|) max(|P|, |Q|)`.
5. For (3): `|M_0^E rho^i| = 2^(Li) M_0^(E-i)` for `0 <= i <= E`. ∎

**Every formal construction is of this form.**
* The partial sums are `A = 1`, `B = S_K(x)`, with `E = K^3` and `W' = (K+1)^3`.
* Zudilin's forms (Lemma P) become polynomials after multiplication by `rho^(e_off)`, and their leading remainder coefficient is `R_n(q^(-n-1)) = prod_(i=1..n)(1 - rho^(2 alpha i))`. This equals `1 + O(x)`, which is odd. So Lemma P is an instance of Proposition F.
* A two-variable construction from `G(u,v) = 1 + uv rho G(uv^2 rho^3, v rho^3)` evaluated at the tail points `(rho^(3s^2), rho^(3s))` yields linear forms in `1` and `X` whose coefficients are Laurent polynomials in `rho`, so it is again a pair `(A, B)`.
* At other points `(rho^a, rho^b)`, `G` takes new values: T3 finds no relation between them and `X` below norm `2^7999`. Forms there say nothing about `X`.
* Type-I Hermite–Padé forms in `(1, Theta_3, Theta_3^(1))` behave the same way. A p-adic Nesterenko-type count would give `dim >= 3/mu_bar = 2.10`, so full independence would follow, but only with the same two missing ingredients: small heights and non-vanishing. This is heuristic: the p-adic criterion (Nesterenko, Manuscripta 139 (2012)) was seen as metadata only.

### 3.2 Corollary F3 (reduction of HYP-9127; PROVED)

**Statement.** Suppose that for some `w > mu_bar = 1.4265` and infinitely many `E`, some `A` in `Z[x]` of degree `<= E` and height `2^(o(E))` satisfies `A Theta_3 = B mod x^(ceil(wE))` with `deg B <= E`, and the leading coefficient of its remainder is not divisible by `2^10`. Then HYP-9127 holds.

**Siegel's lemma supplies the gap.** For every `w < 2` and every `E`, a nonzero `A` with the gap exists with `max |a_i| <= (E+1)^((w-1)/(2-w) + o(1))`.

*Proof.*
1. The margins are at least `10 ceil(wE) - 14.265 E - O(log E) - o(E)`, which tends to infinity. Apply Proposition F.
2. The gap conditions are `ceil(wE) - E - 1` homogeneous equations in `E + 1` unknowns with coefficients in `{0, 1}`. Siegel's lemma gives a nonzero solution bounded by `(E+1)^(M/(N-M))`. ∎

**The missing ingredient** is therefore exactly a **2-adic zero estimate**: some short vector of the gap lattice `Lambda_(E,W)` must have its first remainder coefficient not divisible by `2^10`.
* By the Newton polygon, `R(rho) = 0` would force a segment of slope `-L`, hence `v_2(r_(W')) >= L`. Any first coefficient of absolute value below 1024 already rules this out.
* Heuristically such a vector exists with probability `1 - 2^(-10)` per trial. No argument forces it: the functional `A -> R_A(rho)` could in principle vanish on the whole lattice, or on its short vectors. This is the classical impasse of Thue–Siegel methods without a functional equation.

**Two-form version.** If `X = a/b` and two forms `A_1, A_2` both vanish at `rho`, then `delta(rho) = 0`, where `A_1 B_2 - A_2 B_1 = x^W delta(x)` with `deg delta <= 2E - W`. The pair therefore certifies as soon as `v_2(A_2(0) r_W(A_1) - A_1(0) r_W(A_2)) < L`. This is the same parity condition in another form.

### 3.3 Data (`.out` P2, P4)

**Diagonal Padé `[n/n]` of `Theta_3`, evaluated at `rho = 2^10/3^9`.** The budget is `(2L - log2 M_0) n = 5.735 n`.

| `n` | `log2 h` | `W'` | `v_2(r_W')` | NP non-vanishing | `v_2(PX - Q)` (= prediction) | margin (bits) |
|---|---|---|---|---|---|---|
| 100 | 0.0 | 201 | 0 | yes | 2010 | 583 |
| 200 | 12.3 | 401 | 0 | yes | 4010 | 1145 |
| 400 | 67.2 | 801 | 2 | yes | 8012 | 2240 |
| 600 | 155.7 | 1201 | 0 | yes | 12010 | 3297 |
| 800 | 273.0 | 1601 | 0 | yes | 16010 | 4325 |
| 1000 | 451.4 | 2001 | 3 | yes | 20013 | 5297 |

* Each row is a FINITE certificate, much weaker than reconstruction (T1).
* The heights track the Hankel determinants: `log2 h(1000) = 451`, and `bits(H_1000) = 451` as well.

**Siegel forms (P4).** For every tested `(E, w)`, `E <= 300`, `w` in `{1.5, 1.75, 1.9}`:
* the LLL-shortest kernel vector has height 1, except height 129 at `E = 300`, `w = 1.9`;
* the first remainder coefficient is odd;
* the margins are at least `(10w - 14.265) E - log2 h`, and exceed `(10w - 14.265) E` whenever `h = 1`.

At small `E`, many of these are cube-gap coincidences: `A = 1 - x^(d_k)` with `d_k = (k+1)^3 - k^3` cancels `(k+1)^3`. For example, `E = 155` and `W' = 307` at `k = 5`. The relative gap of this family, `6/k`, tends to 0.

A proof needs growing `E`. There Siegel's lemma controls heights (`2^(O(log E))`) but not the leading coefficient's 2-adic valuation: at `E = 600`, `w = 1.75` the shortest vector already has `h = 655200` and `v_2(r_W') = 7`.

### 3.4 Why diagonal Padé cannot be made asymptotic: the Hankel census (`.out` P1)

Let `H_n = det(c_(i+j))_(i,j<n)` for 0/1 sequences `c`. The controls are random sequences with one element uniform in each `[k^d, (k+1)^d)`, so they have the same counting function.

| `n` | cubes | random ~ cubes (2 seeds) | squares | random ~ squares (2 seeds) | `floor(k^1.5)` | random ~ `floor(k^1.5)` |
|---|---|---|---|---|---|---|
| 400 | 68 | 15, 72 | 116 | 472, 458 | 589 | 768 |
| 800 | 275 | 324, 319 | 262 | 1075, 1120 | 1434 | 1801 |
| 1000 | 451 | 474, 434 | 361 | 1411, 1466 | 1930 | 2375 |
| 1400 | 665 | 729, 743 | 505 | 2221, 2193 | 3044 | 3509 |
| 2000 | 1063 | 1172, 1035 | – | – | – | – |

Entries are bits of `|H_n|`.

**Cubes.** Their determinants are statistically indistinguishable from the random controls: `bits / (n log2 n)` is 0.045 against 0.050 at `n = 1400`, and 0.0485 against 0.0534 and 0.0472 at `n = 2000`. There is **no Padé miracle**, which is FINITE-EXACT evidence for genericity.

**The growth is superlinear**, `kappa n log2 n`, while the certificate budget is linear (`5.735 n`). So diagonal Padé forms cannot prove HYP-9127, even though they would keep certifying up to about `log2 n ~ 5.7/kappa ~ 120`.

**Squares** are about 4.4 times smaller than random: `0.035 n log2 n` against `0.15 n log2 n`. This is the formal shadow of the theta structure, and it matters for the cheaper instances in section 7.

**Pade-table normality.** The cube Padé table is non-normal at 88 of the first 300 indices, against 138–163 for the random controls. Among `n <= 300`, 143 values of `H_n` are odd.

### 3.5 Proposition C (the cubic analogue of Zudilin's construction; FINITE-EXACT, with a heuristic reason)

**The construction.** Choose `c_0, ..., c_n` with `R(t) := sum_s c_s rho^(3ts(t+s)) = 0` for `t = 0..n-1`: the kernel of an `n x (n+1)` matrix. Then

`sum_s c_s rho^(-s^3)(X - S_(s-1)) = sum_(t>=n) R(t) rho^(t^3)`.

**Results** (exact, `.out` L4):

| `n` | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| exponent `v_2/log2 H` | 1.40 | 0.77 | 0.51 | 0.42 | 0.38 | 0.33 |
| partial sums `K = n` | 5.58 | 2.36 | 1.66 | 1.37 | 1.21 | 1.11 |

The construction is worse at every `n`.

**Reason.** `f(t,s) = ts(t+s)` is supermodular: `d^2 f/(dt ds) = 2(t+s) > 0`. So the 2-adic sizes of the minors come from anti-sorted matchings, about `n^4/6`, and their 3-adic heights from sorted ones, about `n^4/2`. There is no product formula like the Vandermonde's `f = 2ts` to absorb a common factor. The q-binomial cancellation that makes Zudilin's forms small has no cubic counterpart.

## 4. The structural obstructions (PROVED)

### 4.1 Proposition Lv (Liouville exactness)

Let `P` be in `Z[k]` with degree `d >= 2` and positive leading coefficient, and `S = {P(k) : k >= k0}`. The partial sums of `Theta_S(rho)` certify irrationality iff `mu_bar = 1`, i.e. iff `M_0 < 2^L`.

*Proof.*
1. `S_K = A_K/M_0^(P(K))` with `|A_K| <= (K+1) max(2^L, M_0)^(P(K))`.
2. `v_2(Theta_S(rho) - S_K) = L P(K+1)`.
3. So the gain is `L[P(K+1) - mu_bar P(K)] + O(log K)`.
4. It tends to `+infinity` if `mu_bar = 1` (because `P(K+1) - P(K) -> infinity`), and to `-infinity` if `mu_bar > 1` (because `P(K+1)/P(K) -> 1`). ∎

For `x/2, (mx+r)/2` with `M_0 = m^A`: `mu_bar = 1` iff `A/L < log_m 2`, the subcritical case. Lemma L (Monks–Yazinski) already excludes that case for every non-EP word. **So Liouville adds nothing for polynomial positions of degree `>= 2`.**

The near-critical cube swap `A = 12`, `L = 19` (`mu_bar = 1.001029`) has partial-sum gains up to `7.19*10^7` bits (at `K = 1944`), negative beyond `K ~ 2900` (`.out` L6). That is a huge finite certificate, and still no proof.

**Nested gaps (Bradshaw–Ghidelli) in `Q_2` (sketch).** Transcribe Ghidelli's Theorem 3.3 to a combination `sum_j R(j) rho^j` of powers of `X`, where `R(j)` is a combination of representation numbers by sums of cubes.
* At a gap point `n`, the truncation is a rational with denominator `M_0^(n-1)` and numerator of size about `max(2^L, M_0)^n`.
* It is forced to vanish only if the tail valuation `L(n + K)` exceeds `(n - 1) log2 M_0 + O(1)`, i.e. only if the gap length satisfies `K > (mu_bar - 1) n`.
* The supports of `X^j` are relatively dense: every `[n, (1+eps) n]` contains cubes, and also `k^3 + (j-1)` for suitable `k`. So the available gaps are `o(n)`, unless the combination cancels on whole intervals, which would be an identity among representation numbers.

For `mu_bar = 1`, gaps of bounded length suffice, which is why the unit results exist.

### 4.2 Proposition H (no linear q-difference equation; PROVED via Garoufalidis)

**Statement.** Let `y` be an indeterminate, `N >= 1`, `h != 0`, `q = y^h`, and let `P` be a polynomial of degree `d >= 3` with positive leading coefficient. Then `F(z) = sum_k y^(N P(k)) z^k` satisfies no nonzero equation `sum_(i=0..r) c_i(z) F(q^i z) = 0` with `c_i` in `Q(y)[z]`. With `y = rho^(1/N)` this covers every q-difference structure in which `q` is a rational power of `rho`.

*Proof.*
1. Clearing denominators, write `c_i(z) = sum_j c_ij z^j` with `c_ij` in `Q[y^(+-1)]`.
2. The coefficient of `z^n` gives `sum_j (sum_i c_ij y^(h i (n-j))) a_(n-j) = 0`, where `a_k = y^(N P(k))`. This is a linear recurrence whose coefficients are Laurent polynomials in `y` and `y^n`.
3. It is nontrivial: take `j0` minimal with some `c_(i j0) != 0`. The coefficient of `a_(n-j0)` is a nonzero Laurent polynomial in `Y = y^(hn)`, which vanishes for at most finitely many `n`.
4. So `(a_k)` is q-holonomic in Garoufalidis's sense. By his theorem [P, abstract], `deg_y a_k = N P(k)` is eventually a quadratic quasi-polynomial, which is false for `d >= 3`. ∎

The same holds for the positions `floor(k^(3/2))`, whose degree sequence is not a quasi-polynomial. So Zudilin's and Bézivin's q-difference methods have nothing to act on.

For squares, fixing `w = rho` leaves the one-variable linear map `u -> rho^2 u`. That is Tschakaloff's first-order equation. For cubes, the minimal closed system is two-dimensional, `(u, v) -> (rho^3 u v^2, rho^3 v)`, and it is nonlinear in `(u, v)`.

### 4.3 Proposition M (Mahler's method is inadmissible; PROVED via Adamczewski–Faverjon)

The homogenised series `G(u,v,w) = sum u^k v^(k^2) w^(k^3)` satisfies `G = 1 + uvw G(T(u,v,w))` with `T(u,v,w) = (u v^2 w^3, v w^3, w)`. Its exponent matrix is upper unitriangular, with columns `(1,0,0), (2,1,0), (3,3,1)`, so every eigenvalue is 1.

By Adamczewski–Faverjon, Lemma 3.7 (Kubota's observation), a matrix with a root-of-unity eigenvalue admits no point satisfying the non-vanishing condition (c). Moreover `T^k(1,1,rho) = (rho^(3k^2), rho^(3k), rho)` does not tend to 0, so Theorem 3.4 fails as well. The two-variable form with `w = rho` fixed is not a monomial map at all. ∎

### 4.4 Proposition S (the Subspace obstruction, quantified; PROVED) and Theorem U (conditional; PROVED)

**(a) No return compression.** Let `s != 0` and `C` be the set of cubes. The equation `x^3 - y^3 = s` has at most `O(|s|^(1/2))` solutions: `x - y` divides `s`, and `y <= sqrt(|s|/3)`. Hence `1_C(k)` and `1_C(k - s)` differ at `>= 2 floor(N^(1/3)) - O_s(1)` values of `k <= N`.

In Erdős 1062(ii) the analogous count along return times is bounded (`m <= 404`). That is what makes the three-place covering theorem apply. Here every shift identity `(rho^(-s) - 1) X = (finite) + sum(...)` carries a growing number of S-unit terms.

**(b) The bulk integer is costly 2-adically.** In the natural fixed-dimension setup, take the coordinates `y = (N, N S_K, N rho^((K+1)^3))`, with `N = 3^(9(K+1)^3)`: the second is an integer that is not an S-unit, the third is `2^(10(K+1)^3)`. Use the places `{inf, 2, 3}`, with the small form `X y_1 - y_2 - y_3` at 2 and coordinates elsewhere. The product of form values is `3^(9(K+1)^3 - 27K^2 + O(K)) 2^(-10(K+2)^3) = 2^((mu_bar - 1) L K^3 + O(K^2)) > 1`, because the non-S-unit integer is paid in full at `inf`.

In Erdős 1062(ii) the small form is archimedean and the integer `P` is free at the finite places, which is the mirror image. The only way to make the product small is to split `A_K` into its `K + 1` S-unit terms. With all coordinates S-units, the product is `2^(-10(K+1)^3) = H^(-0.70)`, far below `H^(-eps)`, but in dimension `K + 2`, which grows. The Subspace theorem says nothing about one point per dimension.

Known growing-dimension results (Bugeaud–Evertse, via Evertse–Schlickewei) handle about `(log N)^(3/2)` terms among `N` digits, i.e. about `(log log H)^(3/2)`. Cubes need about `(log H)^(1/3)`.

**Hypothesis `U(theta)`** (uniform S-unit gap, sparse regime; `S = {2, 3}`). For every `eps > 0` and `C >= 1` there is an `H_0` with the following property. Let `n <= (log H)^theta`, let `u_1, ..., u_n` be integer S-units with `H = max |u_i| >= H_0`, and let `c_i` be integers with `1 <= |c_i| <= C`. If `Sigma = sum c_i u_i` has no vanishing subsum, then the `{2,3}`-free part of `Sigma` is `>= H^(1-eps)`.

For fixed `n` this is the S-unit sum theorem (Evertse; van der Poorten–Schlickewei) [R]. Pigeonhole constructions reach S-free parts `H^(1 - n/log2 H)`, so they do not contradict `U(theta)` for `theta < 1`.

**Theorem U (PROVED).** `U(1/3)` implies HYP-9127. More generally, `U(1/d)` implies that `sum_k rho^(P(k))` is irrational for every `deg P = d` and every `mu_bar > 1`, with `S` the primes of `2 M_0`.

*Proof.*
1. Let `X = a/b`, with `b > 0` odd and `C = max(|a|, b)`.
2. Put `D_K = b A_K - a 3^(9K^3) = sum_(k=0..K) b u_k - a u_(K+1)`, where `u_k = 2^(10k^3) 3^(9(K^3-k^3))` and `u_(K+1) = 3^(9K^3) = H`. It has `n = K + 2 <= (log H)^(1/3)` terms.
3. `D_K != 0`, because `X != S_K`.
4. No subsum vanishes.
   * A subsum without `-a u_(K+1)` is a positive multiple of a sum of positive terms.
   * `b sum_(k in I) u_k = a u_(K+1)` would give `sum_(k in I) rho^(k^3) = X`, i.e. `sum_(k not in I) rho^(k^3) = 0` over an infinite set of `k`. That is impossible 2-adically: the term with the least `k` dominates.
5. But `v_2(D_K) = 10(K+1)^3` and `|D_K| <= (K+2) C H`. So the S-free part is `<= (K+2) C H^(1 - 10(K+1)^3/(14.265 K^3)) = H^(1 - 1/mu_bar + o(1))`, contradicting `U` with `eps < 1/mu_bar = 0.701`. ∎

**The n-conjecture does not suffice.** Browkin–Brzeziński 1994 (metadata verified; statement [R]) predicts `max |a_i| <= C rad^(2n-5+eps)` for `n`-term relations. Rationality of `X` only forces `rad <= 6|ab| H^(0.30+o(1))`. Then `rad^(2K+1)` is far above `H`, so there is no contradiction.

Effective S-unit bounds (Győry–Yu type) concern exact equations, not 2-adic smallness. No implication from p-adic Schanuel is claimed.

### 4.5 Remarks on what fails in `R`

In the real analogue the same partial sums have error `(2^10/3^9)^((K+1)^3)`. The exponent ratio tends to `1 - log a/log b = 0.299`, further from Liouville than the 2-adic `0.701`.

A two-place argument in the style of Proposition T needs `14.26(K+1)^3 > 28.5 K^3`, which fails. `X` does not converge 3-adically (`|rho|_3 > 1`). No place helps.

## 5. Partial results that were sought (task 2b)

1. **Heights.** The best unconditional statement is finite: **no rational `a/b` with `|a|, |b| <= 2^4999999` equals `X`** (T1).
   * PARI `bestappr` returned `[]` at `N = 10^7` bits, in 57 s; its default bound is `sqrt(2^N/2)`.
   * It was validated on planted rationals: one just below the bound is recovered, and one just above is not.
   * Python Wang reconstruction agrees at `10^5` bits.

   Earlier lanes had `2^4999`.
2. **Algebraic relations.**
   * `X` is not a root of an integer polynomial of degree `<= 6` whose coefficient vector has Euclidean norm `< 2^5714`.
   * There is no relation `c_0 + c_1 X + sum_(j<=4) c_(j+1) Y_j = 0` with `|c|_2 < 2^6666`.

   These are Gram–Schmidt certificates of the LLL basis at 40000 bits (T3); every exact relation lies in the lattice at every precision.
3. **Special denominators.** No proof that `X` is not in `Z[1/3]` was found.
   * The relation `3^e X = a` translates to `3^e A_K = a 3^(9K^3) mod 2^(10(K+1)^3)`.
   * Every residue class contains integers, and nothing constrains `a`'s size: the archimedean and 2-adic data are decoupled. This is the same growing-dimension problem.
4. **Formal certificates** (section 3.3) reach 5297 bits: weaker than item 1, but structured.

## 6. Experimental hunt for hidden structure (task 2d; `.out` T2, T3, P1)

* **T2: 2-adic approximation profile.** For each `N` on a geometric grid in `[64, 20000]`, take the Lagrange-reduced shortest `(a, b)` with `a = bX mod 2^N`, and compute the ratio `log2 max(|a|,|b|)/(N/2)`. A generic number stays near 1.

  | number | min ratio, `N >= 2000` | mean | `#N` with ratio `< 0.95` |
  |---|---|---|---|
  | `X` | 0.9985 | 0.9997 | 0 |
  | `Theta_2(rho)` | 0.9980 | 0.9996 | 0 |
  | `Y_1` | 0.9980 | 0.9997 | 0 |
  | random (2 seeds) | 0.9977, 0.9957 | 0.9997, 0.9996 | 0, 0 |
  | control `sum 2^(3^k)` (Hadamard gaps) | **0.6743** | 0.8784 | 46 |
  | control `sum 2^(k^3)` (unit cubic) | 0.9972 | 0.9996 | 0 |

  The small-`N` dips of `X` (0.368 at `N = 78`) are the partial sums `S_1, ..., S_3`. There is **no super-Dirichlet approximation**: the irrationality-exponent proxy is 2 at every scale.
* **T3: LLL relation hunts.** The lattice is `{c : sum c_i xi_i = 0 mod 2^N}` with `xi_0 = 1`, and the Dirichlet baseline for the shortest vector is `2^(N/m)`.

  | family (`m` numbers) | `N = 20000`: `log2 |c|` / baseline | `N = 40000`: ratio | certificate at 40000 bits |
  |---|---|---|---|
  | `1, X, X^2, X^3` | 4998.6 / 5000 | 0.9999 | no relation below `2^9999` |
  | `1, X, ..., X^6` | 2856.4 / 2857.1 | 0.9999 | below `2^5714` |
  | `1, X, Y_1, Y_2` | 4999.2 / 5000 | 0.9999 | below `2^9999` |
  | `1, X, Y_1, ..., Y_4` | 3332.3 / 3333.3 | 0.9998 | below `2^6666` |
  | `1, X, Theta_2, X Theta_2` | 4999.5 / 5000 | 0.9999 | below `2^9999` |
  | `1, X, G(rho,1), G(1,rho), G(rho,rho)` | 3999.1 / 4000 | 0.9999 | below `2^7999` |
  | `1, X, X^2, Y_1, X Y_1, Y_1^2` | 3332.7 / 3333.3 | 0.9999 | below `2^6666` |
  | control random | 4999.7 / 5000 | 0.9999 | – |
  | control `1, c, c^2, c^3`, `c^3 = 17` in `Z_2` | **4.1** / 5000 | 0.0004 | the relation, found EXACTLY |

  The achieved exponents `tau/m` lie between 1.0001 and 1.0003 for every `X` family. That is the Dirichlet/pigeonhole baseline, and there is no systematic excess, so **no Padé miracle** is visible. Together with P1 (cubes Hankel = random), this supports genericity of `X`.

## 7. The degree hierarchy (task 2e)

### 7.1 Theorem Q (quadratic tier; PROVED)

**Statement.** Let `P(k) = alpha k^2 + beta k + gamma` be integer-valued with `alpha > 0`, so that `a := 2 alpha` and `s := alpha + beta` are integers. Let `eps = +-1`, and let `rho = 2^L/M_0` (`M_0` odd, `L >= 1`) with `mu_bar < phi`. Then `sum_(k>=k0) eps^k rho^(P(k))` is irrational in `Q_2`, for every `k0`.

**Corollary (swap words).** Let `T` be an affine 2-adic shift map and `B != B'` blocks of length `L` with equal multiplier `M_0`, and let `P(k) >= 0` be injective for `k >= k0`. The word with `B'` at block indices `P(k)` and `B` elsewhere has an irrational Bernstein number if `mu_bar = max(1, log2(M_0)/L) < phi`, and so does every word that eventually agrees with it.
* Under `3x+r`, `mu_bar <= log2 3 < phi`, so every such word is settled. Examples: triangular numbers, `k(3k-1)/2`, `k(3k+1)/2`, `2k^2 + k`.
* With a third block `B''` satisfying `R_(B'') = 2R_B - R_(B')`, the sign-twisted words `eps = -1` are settled too.

*Proof.*
1. **Reduction.** WLOG `gamma = 0`. Then `P(l) = a l(l-1)/2 + s l`. Put `q = rho^(-a)` (`|q|_2 = 2^(aL) > 1`) and `z = eps rho^s`. Then `sum_(l>=0) eps^l rho^(P(l)) = T_q(z) = sum_l z^l q^(-l(l-1)/2)`, and the finitely many terms with `l < k0` are rational.
2. **The forms.** Use Zudilin's forms exactly as in Lemma P: `R_n(T) = prod_(j=1..n)(1 - q^j T) = sum_k C_k T^k` and

   `I_n = sum_(t>=1) R_n(q^(-t)) z^(t+m) q^(-(t+m)(t+m-1)/2) = A_n T_q(z) - B_n`.
3. **Size of `I_n`.**
   * `R_n(q^(-t)) = 0` for `1 <= t <= n`.
   * For `t >= n+1`, every factor `1 - q^(j-t)` is a 2-adic unit.
   * The `t`-th term has absolute value `2^(-L P(t+m))`.
   * `P(l+1) - P(l) = a l + s > 0` for large `l`, so for large `n` the term `t = n + 1` dominates, and `|I_n|_2 = 2^(-L P(n+m+1))`. In particular `I_n != 0`.
4. **Exponent bookkeeping.** `C_k` has `q`-degrees in `[k(k+1)/2, k(k+1)/2 + k(n-k)]`.
   * The monomials of `A_n` are `+-rho^e` with `e` in `[-k(a(n+m) + s), -(a(k^2 + km) + s k)]`.
   * The monomials of `B_n` are the same, plus `P(l)` for `0 <= l <= k+m`, whose maximum over `k` is `P(m)`.
   * Hence every monomial lies in `[-e_off, e_+]`, with `e_off = n(a(n+m) + s) + O(1)` and `e_+ = max(0, P(m)) + O(1)`.
   * `D_n = 2^(L e_off) M_0^(e_+)` clears all denominators, and `|D_n A_n|, |D_n B_n| <= 2^n (n+m+1) max(2^L, M_0)^(e_off + e_+)`.
5. **The contradiction.** If `T_q(z) = a'/b'`, then `G_n = a' D_n A_n - b' D_n B_n = b' D_n I_n` is a nonzero integer with `v_2(G_n) = L(e_off + P(n+m+1))`. So

   `L(e_off + P(n+m+1)) <= log2(|a'|+|b'|) + n + log2(n+m+1) + mu_bar L (e_off + e_+)`.
6. **Asymptotics.** Take `m = floor(xn)` and divide by `(a/2) L n^2`. The inequality becomes `(x^2 + 4x + 3) - mu_bar (x^2 + 2x + 2) <= O(1/n)`. At `x = 1/phi` the left side equals `(x^2 + 2x + 2)(phi - mu_bar) > 0`, since `(x^2+4x+3)/(x^2+2x+2) = phi` there. This is a contradiction for large `n`.
7. **Sign twist.** Replacing `z` by `-z` changes neither 2-adic sizes nor heights. ∎

**Checks** (`.out` L3). For five families (squares, triangular, both pentagonal branches, `2k^2+k`) and `n <= 12`, `v_2` equals the prediction in every case.

| map (`mu_bar`) | margins |
|---|---|
| `3x+1` (1.4265), `5x+1` `1^6 0^4` (1.393) | all positive and growing (squares under `3x+1`: `m = 0` from +27 to +411, `m = n/phi` from +54 to +1305, for `n <= 12`) |
| `5x+1` `1^7 0^3` (1.6253 `> phi`) | `m = 0`: negative from `n = 7` (-186 at `n = 12`); best `m`: still +634 at `n = 40` (L3b), but the leading coefficient is `-0.0265 alpha L n^2`, so eventually negative |
| `5x+1` `1^8 0^2` (1.8575) | `m = 0` negative from `n = 3`, `m = n/phi` from `n = 4` (-923 at `n = 12`) |

### 7.2 Where `phi` bites, and whether it can be raised

| map | supercritical iff | Theorem Q covers |
|---|---|---|
| `3x+r` | `A/L > 0.6309` | all (`A/L < 1`) |
| `5x+1` | `A/L > 0.4307` | `A/L < 0.6969` |
| `7x+1` | `A/L > 0.3562` | `A/L < 0.5764` |
| `9x+1` | `A/L > 0.3155` | `A/L < 0.5104` |

**Nearest open square instances** (`L <= 40`):
* `5x+1`, 23 ones in 33 letters: `mu_bar = 1.61831` (`phi + 0.00028`);
* `7x+1`, 15 in 26: `1.61963`;
* `5x+1`, `1^7 0^3`: `1.62535`.

Zudilin's one-parameter family is optimal at `x = 1/phi`: the value `phi` is the maximum of `(x^2+4x+3)/(x^2+2x+2)`. No improvement of `phi` for 2-adic theta values was found in the literature. For general `T_q(z)`, Zudilin (2005) reproves Tschakaloff–Bundschuh with the same threshold `gamma_0 = (3-sqrt5)/2`, and no later improvement was found; whether one exists is UNVERIFIED.

**Diagonal Padé forms (P3) do beat `phi` at every tested `n`.** At `n = 800` the margins are +2747 at `mu_bar = 1.6253`, +889 at `1.8575` and +7677 for `7x+1` (15 in 26). This holds because the square indicator's Hankel data sit 4.4 times below random. Their heights (`0.34 n` bits at `n = 800`) still grow like `n log n`, so this is FINITE-EXACT evidence, not a proof. A structural bound `log2 h <= c n` with `c < L(2 - mu_bar)` would settle square swaps up to `mu_bar < 2 - c/L` (HC-CT4).

### 7.3 The table

| positions `S` | functional equation | method | `3x+r` | `5x+1` |
|---|---|---|---|---|
| arithmetic progression | rational generating function | geometric series | EP (rational) | EP |
| Hadamard, `s_(k+1)/s_k >= lambda > mu_bar` | – | Liouville / Theorem D | PROVED | PROVED iff `lambda > mu_bar` |
| one quadratic family (squares, triangular, one pentagonal branch, `alpha k^2 + beta k + gamma`), optional sign `(-1)^k` | first-order q-difference (Tschakaloff) | **Theorem Q** | **PROVED** | PROVED iff `mu_bar < phi`; OPEN above (cheaper instance) |
| two quadratic families (bilateral pentagonal `T_q(rho) + T_q(rho^2)`, squares together with triangular numbers) | two Tschakaloff values, `z_1/z_2` not in `q^Z` | would need Tschakaloff–Bundschuh linear independence, transcribed 2-adically (not done) | OPEN here | OPEN |
| signed bilateral pentagonal `prod (1 - rho^n)` | `F(z) = (1-z) F(rho z)` (q-exponential) | KRVZ-type Padé, not transcribed | OPEN here | OPEN |
| cubes, any cubic, degree `>= 3` | **none** (Proposition H); Mahler inadmissible (Proposition M) | Liouville only at `mu_bar = 1` (subcritical) | **OPEN** (HYP-9127 at 1.4265) | OPEN |
| `floor(k^(3/2))`, other non-polynomial growth | none | – | OPEN, and not cheaper (section 8) | OPEN |

## 8. Cheaper open instances (task 2f) and what counts as progress (task 3)

**Ranking of open instances beyond Theorem D** (`Dio = 1` for all, `.out` L7). They are ordered by the mechanisms that still act on them.
1. **Square swaps at `mu_bar` in `[phi, 2)`.**
   * These have a first-order q-difference equation, a modular theta structure (the value is `(theta_3(rho) + 1)/2` on the Tate curve), and a sub-random Hankel table.
   * Diagonal Padé certificates exist at every tested `n`.
   * The closest case to the proved range is `5x+1` with 23 ones in 33 letters (`mu_bar = 1.61831`).
2. **Two quadratic families** (bilateral pentagonal). These should be within reach of a 2-adic transcription of Tschakaloff–Bundschuh linear independence for small `mu_bar`.
3. **The near-critical cube swap**, `sum (2^19/3^12)^(k^3)` in `Q_2` (`A = 12`, `L = 19`, `mu_bar = 1.001029`).
   * This is the cheapest cubic instance.
   * Any proof must beat Liouville by a fixed factor with no functional equation, which is exactly the missing mechanism.
   * Partial sums already exclude rationals of height `<= 2^(7*10^7)`.
4. **`Y3` itself** (HYP-9127, `mu_bar = 1.4265`).
5. **`floor(k^(3/2))`: not cheaper** (REFUTED as a cheaper instance).
   * Proposition H applies.
   * L7 gives a repetition profile `LPF(j)/j` decaying like `j^(-0.48)` in range: 0.061 at `2^12`, 0.009 at `2^18`, 0.003 at `2^21`. So `Dio = 1`. The argument: matched gap runs have length at most `(8/3) sqrt(k)` gaps near index `k`, because `(j+i)^(3/2) - (k+i)^(3/2)` drifts by `(3/4)(j-k)/sqrt(k)` per step. This gives `LPF(j) = O(j^(2/3))`. It is a sketch, not audited.
   * Its Hankel determinants are 0.87 of random.
   * A Subspace argument needs `N^(2/3)` terms.

**What would count as progress on HYP-9127:**
* **(P1)** An irrationality proof for `sum rho^(k^3)` at *any* `rho` with `mu_bar > 1`, for instance at `mu_bar = 1.001`, or for the real `sum (a/b)^(k^3)` with any `a >= 2`. This would be the first cubic result beyond Liouville, in either world.
* **(P2)** The 2-adic zero estimate HC-CT1, a statement about the 0/1 cube series alone (Corollary F3).
* **(P3)** Any uniform-in-`n` S-unit gap in the sparse regime `n <= (log H)^theta`, `theta >= 1/3` (Theorem U).
* **(P4)** A sub-`n log n` bound on the Padé heights of `Theta_3`. This would contradict the random-like data of P1, so it is not expected.
* **Not progress on HYP-9127, though useful:** raising `phi` for 2-adic theta values, and 2-adic linear independence of two Tschakaloff values. These settle cheaper instances 1 and 2.

## 9. Hypothesis candidates (for `INDEX.md` numbering by the integrator; no HYP files created)

* **HC-CT1 (2-adic zero estimate for cube Padé lattices; OPEN).** For some `w` in `(1.4265, 2)` and infinitely many `E`, the lattice `{A in Z^(E+1) : A Theta_3 = B mod x^(ceil(wE)), deg B <= E}` contains a vector of height `2^(o(E))` whose first remainder coefficient is not divisible by `2^10`.
  * It implies HYP-9127 (Corollary F3).
  * Evidence: true for every tested `E <= 1000` (diagonal forms, `v_2(r_W') <= 3`) and `E <= 300` (Siegel forms, `r_W'` odd).
* **HC-CT2 (`U(1/3)`; OPEN).** Uniform S-unit gap for `n <= (log H)^(1/3)` terms (section 4.4). It implies HYP-9127 and every higher-degree analogue.
* **HC-CT3 (no Padé miracle for cubes; OPEN, FINITE-EXACT evidence).** `log2 |H_n(cubes)| = (kappa + o(1)) n log2 n`, with the same `kappa` as random sequences with the cube counting function (about 0.045–0.05 at `n <= 2000`).
* **HC-CT4 (square Padé heights; OPEN).**
  * The diagonal Padé heights of `sum x^(k^2)` satisfy `log2 h_n <= c n` for some constant `c`.
  * Proposition F would then settle square swaps for `mu_bar < 2 - c/L`. With the measured `c ~ 0.34` (`n = 800`) and `L = 10`, that means `mu_bar < 1.966`, i.e. every `5x+1` block with `A/L < 0.84`.
  * Data: the heights grow like about `0.035 n log2 n`, so the linear bound is probably false. A proof either way would be informative.
* **HC-CT5 (near-critical cube; OPEN).** `sum (2^19/3^12)^(k^3)` is irrational in `Q_2`. It is the smallest cubic instance of PC.
* **HC-CT6 (bilateral pentagonal; OPEN).** `T_q(rho) + T_q(rho^2)`, with `q = rho^(-3)`, is irrational in `Q_2` for `mu_bar < mu_2` for some explicit `mu_2 > 1`.
* **HC-CT7 (square swap beyond `phi`; OPEN).** `Phi_(5x+1)` of the square-swap word with blocks of 23 ones in 33 letters (`mu_bar = 1.61831`) is irrational.

## 10. Reproduction

```bash
cd <worktree>
bash 04-computation/experiments/collatz_procgen_20260923_cube_theta_run.sh          # full run
#   writes 05-knowledge/results/collatz_procgen_20260923_cube_theta.out
#   ~6 minutes (ledger 11 s, lattice 71 s, Pade 244 s), one process at a time, peak RSS ~380 MB; deterministic (fixed seeds).
#   Needs python3 + gmpy2 + python-flint; PARI/GP (gp) for the 10^6 / 10^7-bit bestappr rows.
bash 04-computation/experiments/collatz_procgen_20260923_cube_theta_run.sh --quick  # ~1.5 minutes
#   writes scratch/procgen_cube/cube_theta_quick.out
```

**Sections of the `.out`:**
* **Ledger.**
  * L0: constants.
  * L1: `Y3` identity, independent check, and real values.
  * L2: Liouville ledger.
  * L3 and L3b: Theorem Q forms.
  * L4: cubic Zudilin analogue.
  * L5: threshold map.
  * L6: near-critical cube.
  * L7: repetition profiles.
* **Lattice.**
  * T1: reconstruction certificates.
  * T2: approximation profile.
  * T3: LLL relation hunts with Gram–Schmidt certificates.
* **Padé.**
  * P1: Hankel census.
  * P2: diagonal Padé of `Theta_3` at `rho`.
  * P3: diagonal Padé of `Theta_2` beyond `phi`.
  * P4: Siegel forms.

**Certificate methods.** All 2-adic quantities are exact integers modulo `2^N`. Valuations of linear forms are compared with the Newton-polygon or Zudilin prediction, and matched in every row. PARI's `bestappr` bound was validated on planted rationals. Relation certificates use exact rational Gram–Schmidt on the LLL basis.

**Scratch exploration** (not needed for the `.out`): `scratch/procgen_cube/{hankel_probe,pade_probe,siegel_explore,dio_probe,t_speed}.py`. Literature helper scripts are in `scratch/procgen_cube/lit/`; downloaded texts were moved out of the worktree.

## Sources

`[P]` means the primary text was read in this lane, `[P, abstract]` that only the abstract was read, `[R]` that a named secondary source was used, and "metadata" that only the bibliographic record (Crossref) was seen. UNVERIFIED marks everything else.

* **Cubic and lacunary series.**
  * L. Ghidelli, *Arithmetic properties of cubic and biquadratic theta series*, arXiv 1910.05076v1 (2019). [P: abstract, sections 1–2, Theorem 1.1, Remark 2.1, Definitions 3.1–3.2, Theorem 3.3, references.]
  * R. Bradshaw, *Arithmetic properties of values of lacunary series*, MSc thesis, Ottawa (2013). [R via Ghidelli: `deg theta_l(q) >= l`.]
  * D. H. Bailey, J. M. Borwein, R. E. Crandall, C. Pomerance, *On the binary expansions of algebraic numbers*, J. Théor. Nombres Bordeaux 16 (2004) 487–518. [P, abstract: `#(1-bits <= N) > C N^(1/D)`.]
  * Y. Bugeaud, J.-H. Evertse, *On two notions of complexity of algebraic numbers*, arXiv 0709.1560v1. [P: section 1 (1.1)–(1.3), Theorem 3.1, Corollary 3.2, and the Ridout / Schneider Satz 7 remark.] The journal version is UNVERIFIED.
  * D. Ridout (1957); Th. Schneider, *Einführung in die transzendenten Zahlen* (1957), Satz 7. [R via Bugeaud–Evertse.]
  * P. Corvaja, U. Zannier, *Some new applications of the subspace theorem*, Compositio Math. 131 (2002) 319–340. [P, abstract via Crossref; the exact lacunary statements are UNVERIFIED.]
  * G. Çalışkan, *Some lacunary power series and Mahler's U_m-numbers in p-adic domain*, C. R. Acad. Bulgare Sci. 75 (2022) 477–485. [P, abstract.]
* **q-series and theta values.**
  * W. Zudilin, *An elementary proof of the irrationality of Tschakaloff series*, arXiv math/0506086 (J. Math. Sci. 146 (2007)). [P: Theorem with `gamma_0 = (3 - sqrt5)/2`, construction (3)–(4).]
  * L. Tschakaloff (1919); P. Bundschuh (Satz 2); O. Szász; F. Bernstein–O. Szász. [R via Zudilin.]
  * C. Krattenthaler, I. Rochev, K. Väänänen, W. Zudilin, *On the non-quadraticity of values of the q-exponential function and related q-series*, Acta Arith. 136 (2009) 243–269, arXiv 0812.2921. [P, abstract.]
  * A. B. Dixit, V. Kumar, S. S. Pathak, arXiv 2211.03030. [P, abstract.]
  * K. Väänänen, R. Wallisser, JNT 39 (1991) 225–236. [Metadata; hypotheses UNVERIFIED.]
  * J.-P. Bézivin, Manuscripta Math. 61 (1988) 103–129 and Acta Arith. 55 (1990) 233–240. [Metadata.]
  * J.-P. Bézivin, Math. Nachr. 190 (1998). [R via Ghidelli.]
  * J.-P. Bézivin, P. Robba, Ann. of Math. 129 (1989). [Metadata.]
  * Yu. V. Nesterenko (1996); D. Bertrand, Ramanujan J. 1 (1997), Theorem 4; D. Duverney, K. Nishioka, K. Nishioka, I. Shiokawa, Proc. Japan Acad. 72 (1996) 202–203; D. Duverney (1993, 1995). [R via Ghidelli. The DNNS and Springer pages were behind bot or JavaScript walls and were not bypassed.]
  * K. Barré-Sirieix, G. Diaz, F. Gramain, G. Philibert, *Une preuve de la conjecture de Mahler–Manin*, Invent. Math. 124 (1996) 1–9. [Metadata; statement R.]
  * Yu. V. Nesterenko, *On a criterion of linear independence of p-adic numbers*, Manuscripta Math. 139 (2012). [Metadata.]
* **Functional equations.**
  * B. Adamczewski, C. Faverjon, *Mahler's method in several variables I*, arXiv 1809.04823v1. [P: Definition 1.2, Remark 1.3, Definition 3.1, Remark 3.2, Definition 3.3, Theorem 3.4, Lemma 3.7.]
  * D. Masser, *A vanishing theorem for power series*, Invent. Math. 67 (1982) 275–296. [Metadata; R via Adamczewski–Faverjon.]
  * S. Garoufalidis, *The degree of a q-holonomic sequence is a quadratic quasi-polynomial*, arXiv 1005.4580. [P, abstract.]
* **S-units.**
  * J.-H. Evertse, S-unit sums (Compositio 1984); A. J. van der Poorten, H. P. Schlickewei (1982). [R, standard; not re-read.]
  * J. Browkin, J. Brzeziński, *Some remarks on the abc-conjecture*, Math. Comp. 62 (1994) 931. [Metadata; the n-conjecture statement is R.]
  * J.-H. Evertse, H. P. Schlickewei (2002). [R via Bugeaud–Evertse.]
* **Erdős 1062(ii).** L. Kruer, J. Kohlmeyer, exposition of the accepted Lean proof (22 September 2026), session scratchpad copy. [P: sections 3–5, Propositions 5.1–5.2, the dimension bound `m <= 404`.]
* **Repository.**
  * [hard-class note](collatz_procgen_20260922_hard_class.md) (Theorem Y, Lemma P, `Y3` identity, HC1, HC6);
  * [transversality foundry](collatz_procgen_20260922_transversality_foundry.md) (Theorem R, Lemma G, Lemma L);
  * `05-knowledge/hypotheses/HYP-9127-cube-swap-cubic-theta.md`.
