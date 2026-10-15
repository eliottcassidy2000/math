# Least resistance: why the Paley run sums are polynomials, where the first resistance is an elliptic curve, and how rigidity circuits climb a ladder of modular forms

**Contents.**
- Series-parallel ("least resistance") reduction makes the Paley run sums polynomials.
- The first resistance is an elliptic curve (a Y-junction is a Legendre curve).
- Through twelve steps, single runs see a single modular form: the wheel and the octahedron are both the Beukers-supercongruence form.
- The `K_6 − P_6` terminal is the symmetric square of a weight-3 level-32 newform. That is mac-mini's HYP-9276, confirmed here independently in a second normalisation, with the 2-adic `2^7` law traced to a congruence of the newform.
- Wheels are power sums of Jacobi sums (PROVED). `W_6` is the weight-6 level-8 newform.
- A dictionary: 2D rigidity circuits (`|E| = 2|V| − 2`) are exactly the graphs whose double covers are Calabi–Yau. They carry new motives of weight `|V| − 2`.
- Planar Fourier duality (PROVED) fixes this Calabi–Yau line and swaps the Fano and general-type sides. The Fano cube equals `q` times its dual, the octahedron, which refutes this session's own pre-registered "Laman-sparse ⟹ polynomial".
- The owner's triangular-Fibonacci, `F = U + S`, mod-8, Mahler and two-point-Chowla themes, typed.
- Companion note on the owner's "128" request, the 1/3–2/3 and gold partition conjectures, Mahler's middle third and Borsuk: [one_twenty_eight_atlas_balance_gold_20261014.md](one_twenty_eight_atlas_balance_gold_20261014.md).

2026-10-14, session opus-2026-10-14-S35 (worktree `/tmp/math-wt-fractal-recursion-23`).

Owner's prompt (verbatim), quoting `2q(q−1)(q−3)(q−7)(7q−55)` from S34: "consider how these values, and similar formulas you aim to discover, relate deeply to natural fractal modes of recursion least resistance compression through the natural numbers. think of triangular-fibbonacci numbers and F=U+S {10,2,1} and 10=7+3 and as many other related themes to these as we have investigated, and see how they are the hidden part behind mod 8 bott periodicity doubling/squaring exponentiation as the noncommutative unbalanced synthesis of multiplication and addition / try to see the geometry deeply behind processes like the annotated and importantly, other similar ones you are able to identify. merge in mahler's 3/2 related work and other famous problems hinging on these first two primes and see what underlies them. also consider the attached paper for possible inspiration. feel free to leverage maximal creativity in you theory building about possible ways the natural numbers all sit together in harmony for path of least resistance recursion, think of functions as directed graphs freely and see structures between filds of math". The attached paper is the 87-page manuscript *Ordinary two-point correlations of multiplicative functions* (OpenAI, 2026-09-24), read in S25 as M3.

**Status.**
- **PROVED:**
  - the partition-lattice identity for the run sums;
  - THM-4652's reduction calculus and vanishing lemmas, which apply to every core met here;
  - `U(W_4) = q(q−1)·M(q)` by scaling;
  - `M(q) = −(q + a_q)` (Ahlgren–Ono; codex; the S34 audit);
  - the Y-junction is a Legendre curve (classical).
- **FINITE-EXACT:**
  - `A_2 … A_8` reproduced;
  - `A_10` and `A_12` exactly, at 14 primes plus the forced zeros;
  - the octahedron identity at all odd primes `5 … 107`;
  - the `K_6 − P_6` table at all 124 primes `3 … 691` (exact `O(q³)`), and its square law `Y_q = |b(q)|² − q² + 2χ₋₄(q)q`, with `b` the level-32 weight-3 newform. The identification is mac-mini's HYP-9276 (to 30,000), and this session confirms it independently;
  - the wheel laws `T_4`, `T_6` at all 39 primes `q ≡ 3 (mod 4)` in `[7, 400]`;
  - the circuit-reduction tests and the cube identity `U(Q_3) = q·U(K_(2,2,2))` of sections 4c–4d.
- **PROVED (new here):**
  - wheels are power sums of Jacobi sums, `Tr((DS)^n) = Σ_ψ (χ(−1)J(ψ, χ))^n`;
  - planar Fourier duality, `U(G) = χ(−1)^{|E|/2} q^{(|V| − |V*|)/2} U(G*)`.
- **REFUTED (pre-registered here):** "Laman-sparse pure cores are polynomial". The cube is the counterexample.
- **OPEN:** HYP-9275 (A); HYP-9278 (the rigidity/Calabi–Yau dictionary).
- **CITED:** Duffin (1965); Ahlgren–Ono (2000); Beukers (1987); Deligne (1974); Luo (1989); Evans–Pulham–Sheehan (1981); Mahler (1968); Flatto–Lagarias–Pollington (1995).
- **ANALOGY / NUMEROLOGY:** typed in place.
- **Not yet audited** (section 9).

## 0. Answers in brief

1. **The geometry behind the annotated `A_8 = 2q(q−1)(q−3)(q−7)(7q−55)`.** `A_L` is the signed self-avoiding walk sum of the Paley tournament. It equals a Möbius-weighted sum, over set partitions of the walk's `L + 1` points, of the character sums of the quotient multigraphs:

   `A_L = Σ_π μ(π)·U(path/π)`.

   - Only partitions without consecutive points in a block count, and there are exactly `Bell(L)` of them: 4140 at `L = 8`, 4,213,597 at `L = 12`.
   - **Series composition** of edges is convolution and **parallel composition** is pointwise product, in the three-element algebra `{χ, 1, δ}` (THM-4652). The rules are `χ*χ = 1 − qδ`, `χ·χ = 1 − δ`, and so on.
   - These are literally the resistor rules of an electrical network. Series-parallel graphs, i.e. those with no `K_4` minor (Duffin 1965), reduce completely, and give **polynomials** in `q`.
   - The factors `(q − 3)(q − 7)` of `A_8` are "too few points" zeros: a walk of 8 steps needs 9 distinct points, which `F_3` and `F_7` lack.
2. **The first resistance.** The classical irreducible network is the **Wheatstone bridge `K_4`**.
   - For Paley tournaments it vanishes, because `χ(−1) = −1` makes complete patterns `±χ(Vandermonde)`.
   - The irreducible step that survives is the **Y-junction**: summing a degree-3 vertex gives `Σ_x χ((x−a)(x−b)(x−c)) = −a(E)`, the trace of the Legendre elliptic curve through `a, b, c`. Hasse's Riemann hypothesis bounds it by `2√q`.
   - Where an electrical Y–Δ transform is algebraic, here it is an elliptic curve.
3. **The first surviving core is the wheel `W_4`**, and it is the second moment of the Y-junctions.
   - Scaling one rim vertex gives `U(W_4) = q(q−1)·M(q)`, with `M(q) = Σ_λ χ(λ)a(E_λ)²` the twisted second moment of the Legendre family.
   - Ahlgren–Ono give `M(q) = −(q + a_q)`, where `a_q` are the coefficients of `η(2z)⁴η(4z)⁴`, the level-8 form of Beukers' supercongruence for the Apéry numbers (`ζ(3)`).
   - Deligne bounds them by `2q^(3/2)`.
4. **The similar formulas.** The compiler gives exactly:
   - `A_10 = 2q(q−1)(21q⁴ − 639q³ + 7028q² − 32098q + 48662) + 34·q(q−1)·a_q` (S34 audit, codex, tiling-modular-atoms and this compiler agree);
   - `A_12 = P_12(q) + 28q(q−1)(17q − 251)·a_q`, with `P_12` given in section 3.
   - At 12 steps seven cores appear: `K_4`, `K_5`, three six-vertex cores that vanish by odd-reversal symmetry, the wheel, and the **octahedron**. The octahedron is again the same form:

     `U(K_{2,2,2}) = q(q−1)·((5 − q)·a_q + 3q)`, at all odd primes `5 … 107`.

   **So through twelve steps, the Paley run sums live in `Z[q] ⊕ Z[q]·a_q`: one modular form.** The forced zeros `A_L(3) = A_L(7) = A_L(11) = 0` then *determine* `a_3 = −4`, `a_7 = 24` and `a_11 = −44` from the polynomial parts.
5. **A second arithmetic object: a symmetric square.** Tiling-modular-atoms' terminal `K_6 − P_6` enters the weighted `R(q)` expansion and is not in `Z[q] ⊕ Z[q]a_q`. Section 4 gives the details.
   - Its values `Y_q = U/(q(q−1))` satisfy `Y_q + q² − 2χ₋₄(q)q = |b(q)|²` at all 124 primes `3 … 691`.
   - Here `b = q + 4iq³ + 2q⁵ − 8iq⁷ − 7q⁹ − …` is the non-CM newform of `S_3(Γ_0(32), χ₋₄)`, with coefficient field `Q(i)` (McCarthy–Papanikolas' `f_2`). This session built it from level-32 eta quotients.
   - mac-mini found the same identity first (HYP-9276, THM-4653; 3,243 primes to 30,000), in the injective normalisation `W_q = Y_q + 7 + 14χ − 2χq`. The two pipelines agree.
   - The `≡ 1 (mod 128)` signature is `Y_q − 1 = (|b(q)| − q − 1)(|b(q)| + q + 1)` together with the congruence `|b(q)| ≡ ±(q + 1) (mod 16)`.
   - `Y_q = 1` exactly when `|b(q)| = q + 1` (`q = 3, 7, 31, 199, 223, …`).
   - `b(127) = 0` at the Mersenne prime 127.
   - `Y_191 = 1 + 2^16` and `Y_383 = 1 − 2^17`, because both factors are powers of 2.
5b. **Wheels and the rigidity ladder** (sections 4b–4c).
   - `U(W_n) = q·Σ_ψ (χ(−1)J(ψ, χ))^n` (PROVED).
   - `U(W_6) = q(q−1)(q² − q·a_q − b6_q)`, with `b6` the weight-6 level-8 newform (FINITE-EXACT, 39 primes).
   - Every new motive met so far sits at a 2D rigidity circuit, a graph with `|E| = 2|V| − 2` whose double cover is Calabi–Yau, and has weight `|V| − 2`.
   - Planar Fourier duality maps the sum of a plane graph to that of its dual (PROVED). The Calabi–Yau count `|F| = |V|` is self-dual.
   - The cube, the first pure core in the Fano regime, is not a polynomial: `U(Q_3) = q·U(octahedron)`.
6. **The owner's numbers** (section 5): the exact facts are listed there, and the rest is typed NUMEROLOGY.
   - The affine group `|Aff(P_q)| = T_{q−1}` is a triangular number. It is Fibonacci exactly at `q = 3, 7, 11` (3, 21, 55; Luo's theorem). Those are exactly the orders where the maximum Hamiltonian path count is known and attained by Paley.
   - `H(P_q) ≡ T_{q−1} (mod 8)` (HYP-9273 B), with period 16 in `q`.
   - `10 = 7 + 3`: the two roots of `A_8` sum to 10, and their product is `21 = |Aut(P_7)|`.
7. **Mahler 3/2 and the 2–3 problems** (section 6): cross-referenced to tiling-modular-atoms' carry packets. The common obstruction is the Diophantine nature of `log₂ 3`, whose continued fraction has the partial quotient 55 at index 14.
8. **The attached paper** (section 6): the integer two-point Chowla preprint is the `Z` shadow of what is exact over `F_q`.
   - The two-point finite-field correlation is a Jacobi sum.
   - Higher correlations split into polynomial main terms plus Deligne-bounded modular terms.
   - The first such term is the Beukers form at the wheel.

## 1. The compiler (PROVED identity + FINITE-EXACT evaluation)

`paley_run_partition_compiler_20261014.py`:
- **Möbius inversion** over the partition lattice, `μ(π) = Π_B (−1)^(|B|−1)(|B|−1)!`.
- **Exact reduction** in `{χ, 1, δ}`: loops, parallel products, pendant vertices, series convolutions, and edge splitting (`χ`, delete, or contract).
- **Cores** (pure `χ`, minimum degree `>= 3`) are canonicalised and recorded separately.
- **Exactness.** Arithmetic is exact at 14 primes.

**Validation.**
- It reproduces `A_2, A_4, A_6, A_8` exactly at all 14 primes; these are the closed forms of HYP-9274 (C).
- `A_10(7) = 0` and `A_12(7) = A_12(11) = 0` come out without being imposed.
- Brute force gives `A_10(11) = 478720`.

## 2. Least resistance, the Wheatstone bridge, the Y-junction and the wheel (PROVED + CITED + FINITE-EXACT)

| step | electrical network | Paley character algebra | arithmetic |
|---|---|---|---|
| series | resistances add | convolution: `χ*χ = 1 − qδ`, `χ*1 = 0`, `1*1 = q` | Jacobi sums: polynomials |
| parallel | conductances add | product: `χ·χ = 1 − δ` | polynomials |
| series-parallel graph | reduces to one resistor | reduces to a polynomial (THM-4652 part 5) | `Z[q]` |
| Wheatstone bridge `K_4` | first irreducible network | `±χ(Vandermonde)`, vanishes (THM-4652 part 4) | 0 for tournaments; `4a² − 2q` for Paley graphs (CM by `Q(i)`) |
| Y-junction (degree 3) | Y–Δ is algebraic | `Σ_x χ((x−a)(x−b)(x−c))` | `−a(E_λ)`, Legendre curve, Hasse |
| wheel `W_4` | — | `q(q−1)·Σ_λ χ(λ)a(E_λ)²` | `−q(q−1)(q + a_q)`, weight-4 level-8 Beukers form, Deligne |
| octahedron `K_{2,2,2}` | — | squares of quartic traces | `q(q−1)((5−q)a_q + 3q)`: the same form |
| `K_6 − P_6` (tiling's `V_q`) | rigidity circuit | Legendre bilinear form (THM-4653) | `|b(q)|² − q² + 2χ₋₄(q)q`: symmetric square of a weight-3 level-32 newform (HYP-9276; section 4) |
| wheel `W_6` | rigidity circuit | sixth power sum of Jacobi sums | `(q−1)(q² − q a_q − b6_q)`, weight-6 level-8 newform (section 4b) |

**The six-vertex cores at 12 steps.**
- `K_6` minus (`P_3` plus an edge), `K_(2,4)` plus a matching, and `K_(3,3)` plus an edge each have an automorphism reversing an odd number of edges. So they vanish (THM-4652 part 4; `paley_core_symmetry_20261014.out`).
- The octahedron, the wheel and `K_6 − P_6` have no such automorphism, and their sums are nonzero.

## 3. `A_10` and `A_12` exactly (FINITE-EXACT)

- `A_10 = 2q(q−1)(21q⁴ − 639q³ + 7028q² − 32098q + 48662) + 34q(q−1)a_q`.
- `A_12 = P_12(q) + 28q(q−1)(17q − 251)·a_q`, where `P_12(q) = 4q(q−1)(33q⁵ − 1540q⁴ + 28243q³ − 249550q² + 1027794q − 1488872)`. The `W_4` coefficient in `A_12` is `6688 − 408q` and the octahedron coefficient is `−68` (`paley_A12_assembly_20261014.out`).
- In both, the polynomial part alone does not vanish at the too-few-points primes; the modular term supplies the balance:
  - `A_10(3) = 816 − 816`;
  - `A_10(7) = −34272 + 34272`;
  - `A_12(7) = P_12(7) + c_12(7)·a_7 = 3725568 − 3725568` (by cores: reducible `4912152`, wheel `−4989264`, octahedron `77112`);
  - `A_12(11) = −8673280 + 8673280` (reducible `−5764440`, wheel `7986000`, octahedron `−2221560`);
  - `A_12(3) = −134400 + 134400`, a check outside the interpolation set (`q = 3` is not among the 14 primes).
- **Catalan leading terms (FINITE-EXACT for `k <= 6`; mechanism ANALOGY).** The leading coefficients of `A_2, A_4, …, A_12` are `1, 2, 5, 14, 42, 132`, the Catalan numbers `C_1 … C_6`. So `A_{2k} = C_k q^(k+1) + O(q^k)` through `L = 12`.
  - Heuristic: `S² = J − qI` makes `S/√q` a conference (orthogonal-like) matrix. The top-degree quotients are closed walks folded onto trees with every edge used once each way. Those are the Catalan structures of Wigner's semicircle moments, and of the leading Weingarten asymptotics of Haar orthogonal matrices (Collins–Śniady).
  - The individual Möbius weights `μ(π)` are not all 1 (a tree vertex visited `j` times carries `(j−1)!`), so this is not a proof that the leading coefficient is `C_k` for all `k`.
  - The arithmetic enters only below this free-probability leading order: the wheel at relative order `q^(−5/2)` (Deligne) inside `A_10`.

## 4. The `K_6 − P_6` terminal is a symmetric square (mac-mini's HYP-9276, confirmed independently; 2-adic mechanism)

**Credit.** The identification is mac-mini-2026-10-14-wilson's, pushed earlier today:
- THM-4653 (PROVED reduction to a Legendre bilinear form, with an `O(q log q)` FFT evaluator);
- HYP-9276 (FINITE-EXACT at all 3,243 primes `5 <= q <= 30,000` and 12 prime powers; Sato–Tate; conditional independence from the wheel; literature).

This session found the same law independently, before reading HYP-9276, by a different route.
1. **Exact table.** `Y_q = U(K_6 − P_6)/(q(q−1))` was computed exactly at all 124 primes `3 <= q <= 691` with an `O(q³)` evaluation. Fixing `x_0 = 0` and `x_1 = 1` leaves a sum of Legendre-curve traces times a two-variable sum (`paley_k6p6_fast_table_20261014.py`). It agrees with the brute-force values of the earlier table at all 33 overlapping primes.
2. **Square law.** `Y_q + q² − 2χ₋₄(q)q` is a perfect square `c_q²` at all 124 primes, with `0 <= c_q <= 2q`. The values `c_q/(2q)` have the non-CM Sato–Tate moments: the mean of `(c/2q)²` is 0.231 against 0.25, with standard error 0.022.
3. **The newform from eta quotients.** The weight-3, level-32 eta quotients with character `χ₋₄` (holomorphic, cuspidal, `|r_δ| <= 10`) span a 4-dimensional Hecke-stable space, on which `T_3` has characteristic polynomial `x²(x² + 16)`.
   - The kernel is old: `η(4z)⁶` and `η(8z)⁶`, CM by `Q(i)`.
   - The factor `x² + 16` is the newform pair `b, b̄`: `b = q + 4iq³ + 2q⁵ − 8iq⁷ − 7q⁹ − 4iq¹¹ − 14q¹³ + 8iq¹⁵ + 18q¹⁷ − 12iq¹⁹ + …`.
   - Sanity checks: `b(9) = b(3)² + 9` and `b(15) = b(3)b(5)`.
   - `b(q)` is real for `q ≡ 1 (mod 4)` and purely imaginary for `q ≡ 3 (mod 4)`. These are HYP-9276's coefficients (gp `mf`; McCarthy–Papanikolas' `f_2`).
4. **The identity.** `Y_q = |b(q)|² − q² + 2χ₋₄(q)q` holds at all 124 primes (`paley_k6p6_level32_form_20261014.py`). Equivalently, `Y_q = χ₋₄(q)·(tr Sym²ρ_b(Frob_q) + 2q)`.
5. **Normalisations reconciled.** HYP-9276's `W_q` is the injective sum, over pairwise distinct vertices. This session's `Y_q` is the unrestricted sum `U`. Brute force at `q = 7, 11, 13, 19, 23, 29, 31` gives `W_q − Y_q = 7 + 14χ₋₄(q) − 2χ₋₄(q)q` exactly, and reproduces HYP-9276's `W_7 = 8`, `W_11 = −112`, `W_19 = −224`, `W_23 = 1064`, `W_31 = 56`. At Paley primes `Y_q − 1 = M_q`, in HYP-9276's notation.

**The 2-adic signature, explained down to a congruence of `b`.**
- At Paley primes, `Y_q − 1 = (|b(q)| − (q+1))(|b(q)| + (q+1))`. This factorisation is also in HYP-9276 (d).
- Both `|b(q)|` and `q + 1` are `≡ 0 (mod 4)`, and the two factors differ by `2(q+1)`. So `v_2(Y_q − 1) = v_2(c − q − 1) + v_2(c + q + 1)`, where the smaller valuation is `1 + v_2(q+1)` whenever the larger exceeds it.
- The data (124 primes, plus the extended run below) give the congruence

  `|b(q)| ≡ ±(q + 1) (mod 2^m)`, with `m = 4` for `q ≡ 11 (mod 16)`, `m = 5` for `q ≡ 3 (mod 16)`, and `m = 6` for `q ≡ 7 (mod 8)`.

  Hence `v_2(Y_q − 1) >= 7, 8, 10, 11` for `q ≡ 11, 3, 7, 15 (mod 16)`. The value 7 is attained at every `q ≡ 11 (mod 16)` tested.
- This refines HYP-9276 (e), which records `2^7` for `q ≡ 3 (mod 8)` and `2^10` for `q ≡ 7 (mod 8)` and calls them unexplained. For `q ≡ 1 (mod 4)`, `c_q ≡ ±(q+1) (mod 8)` only, and `v_2(Y_q − 1) = 2` exactly.
- **Proposed mechanism (ANALOGY; a proof route).** McCarthy–Papanikolas compute `b(p)` from Hijikata's trace formula and Schoof's count of Legendre curves with `Z/2 × Z/8` torsion. Such curves have `16 | #E(F_p)`, i.e. `a_p ≡ p + 1 (mod 16)`. An Eisenstein-type congruence `b ≡ (1 + q)·(unit) (mod 16)` is what that structure would produce.
- **Extended range (FINITE-EXACT).** `paley_k6p6_two_adic_mersenne_20261014.py` re-implements THM-4653's FFT evaluator with an exact `int64` final sum, and checks it against the `O(q³)` table. It then tests the square law and the depth table at all primes `q ≡ 3 (mod 4)` in `[30,000, 200,000]`, beyond HYP-9276's range: [EXT200K]

**Exact values and Mersenne/Fermat coincidences (FINITE-EXACT; the coincidences are NUMEROLOGY with an explained mechanism).**
- `Y_q = 1` iff `|b(q)| = q + 1`: `q = 3, 7, 31, 199, 223` below 700. HYP-9276 starts at `q = 5` and lists `7, 31, 199, 223, 1151, …`; `q = 3` belongs to the list too.
- `b(127) = 0`, at the Mersenne prime `127 = 2^7 − 1`. The congruence forces `|b(127)| ≡ ±128 (mod 2^7)`, i.e. `|b(127)| ∈ {0, 128}`, and the value is 0. So `Y_127 = 1 − 128²`.
- `|b(191)| = 320 = 192 + 128` gives `Y_191 = 1 + 128·512 = 1 + 2^16 = F_4`. `|b(383)| = 128` gives `Y_383 = 1 − 256·512 = 1 − 2^17 = −M_17`.
- Mersenne primes `2^p − 1` for `p = 13, 17, 19`: [MERSENNE]

## 4b. The wheel ladder: wheels are power sums of Jacobi sums (PROVED), `W_6` is a weight-6 newform (FINITE-EXACT)

Let `D = diag(χ(x))` and `S_{xy} = χ(y − x)`. Then `U(W_n) = q·T_n` with `T_n = Tr((DS)^n)`: the hub is at 0, the spokes give `D`, and the rim gives `S`.

**Proposition (PROVED).** `T_n = Σ_ψ (χ(−1)J(ψ, χ))^n`. The sum runs over all `q − 1` multiplicative characters, and `J(ψ, χ) = Σ_t ψ(t)χ(1 − t)`.

*Proof.* Terms with some `x_i = 0` vanish. Put `t_i = x_{i+1}/x_i`, so that `Π t_i = 1`. Then

`χ(x_i)χ(x_{i+1} − x_i) = χ(x_i)²χ(t_i − 1) = χ(t_i − 1)`,

and `(x_1, t_1, …, t_{n−1})` parametrises the tuples bijectively. Hence `T_n = (q−1)·φ^{*n}(1)` for `φ(t) = χ(t − 1)` on the cyclic group `F_q^×`. Fourier inversion gives `T_n = Σ_ψ φ̂(ψ)^n`, with `φ̂(ψ) = Σ_t χ(t − 1)ψ(t) = χ(−1)J(ψ, χ)`. ∎

So the nonzero eigenvalues of `DS` are exactly the Jacobi sums `χ(−1)J(ψ, χ)`. For `q ≡ 3 (mod 4)` the odd `T_n` vanish, because `S` is skew.

**Identifications (FINITE-EXACT).** These were found by exact overdetermined fits: `T_6` needed 3 of 8 basis functions at 13 primes (`paley_wheel_ladder_fit_20261014.py`). They were then confirmed at all 39 primes `q ≡ 3 (mod 4)` in `[7, 400]` through the Jacobi-sum formula (`paley_wheel_jacobi_20261014.py`).
- `T_4 = −(q − 1)(q + a_q)`, with `a = η(2z)⁴η(4z)⁴` (weight 4, level 8). This is the S34 wheel.
- `T_6 = (q − 1)(q² − q·a_q − b6_q)`, where `b6 = q + 20q³ − 74q⁵ − 24q⁷ + 157q⁹ + 124q¹¹ + 478q¹³ − 1480q¹⁵ + …` is the newform of `S_6(Γ_0(8))`, the `T_3 = 20` eigenvector. This space is built from `a(z)` times the weight-2 forms `E_2(z) − dE_2(dz)`, `d = 2, 4, 8`, and `T_3` has characteristic polynomial `(x − 20)(x + 12)²`.
- `T_8` is **not** in the span of `q^j` and the Tate twists of all level-2, 4 and 8 eigenforms of weights `<= 8`. `S_8(Γ_0(8))` was computed in full, with dimension 5: newforms `a_3 = 44` and `a_3 = −84`, and the old level-2 form `a_3 = 12`. OPEN, with the prediction that weight-8 forms of level 16 or 32 enter.

The pattern `W_{2m} ↔ weight 2m` is what Birch–Ihara-type moment formulas lead one to expect (ANALOGY). Here it comes from the power sums of one Jacobi-sum family.

## 4c. The rigidity–Calabi–Yau dictionary (HYP-9278; ANALOGY with exact support)

**Why `|E| = 2|V| − 2`.** `U(G) = Σ_x Π_{uv ∈ E} χ(x_v − x_u)` counts points, up to trivial terms, on the double cover of `P^{|V|−2}` (coordinates modulo translation and scaling) branched along the `|E|` hyperplanes `x_u = x_v`. A double cover of `P^n` branched in degree `2d` has `K = π^*O(d − n − 1)`. So it is Calabi–Yau exactly when `|E| = 2(|V| − 1)`, Fano-like below that, and of general type above.
- In rigidity theory, `|E| = 2|V| − 3` with Laman sparsity is minimal rigidity in the plane. `|E| = 2|V| − 2` with all proper subgraphs sparse is a **2D rigidity circuit**.
- Two-trees, the maximal series-parallel graphs, are Laman, built by Henneberg type-I moves. These are precisely THM-4652's polynomial regime.
- Odd `|E|` always gives `U = 0`, since scaling by a non-square flips the sign.

| graph | `|V|`, `|E|` | rigidity type | double cover | Paley sum | status |
|---|---|---|---|---|---|
| series-parallel (2-trees and below) | `|E| <= 2|V| − 3` | Laman | Fano or below | polynomial | PROVED (THM-4652) |
| `K_(3,3) − e`, prism `− e` (reducible by THM-4652 once the degree-2 vertices are series-reduced) | 6, 8 | Laman-sparse | Fano | `−q(q−3)`, `q` (both classes mod 4) | FINITE-EXACT (holdout primes) |
| cube `Q_3`: the first pure core in the Fano regime | 8, 12 | Laman-sparse | Fano | `q·((5 − q)a_q + 3q)`, **not** polynomial: `q` times its dual, the octahedron | FINITE-EXACT (32 primes); section 4d |
| Y-junction (3 neighbours) | — | — | elliptic curve (4 branch points) | Legendre trace, weight 2 | classical |
| `K_4 = W_3` | 4, 6 | circuit | K3 (double sextic) | 0 for tournaments; CM `η(4z)⁶`, weight 3 | CITED |
| `W_4` | 5, 8 | circuit | CY3 (double octic) | `η(2z)⁴η(4z)⁴`, weight 4, level 8 | PROVED (S34) |
| `K_6 − P_6` | 6, 10 | circuit | CY4 | `Sym²` of the weight-3 level-32 newform (GL(3), motivic weight 4) | FINITE-EXACT (HYP-9276; here) |
| `K_(2,4) + 2K_2` (Paley graphs) | 6, 10 | circuit | CY4 | `q(q−1)b_q²`: square of the weight-3 CM form | FINITE-EXACT (HYP-9275 B) |
| `W_6` | 7, 12 | circuit | CY5 | weight-6 level-8 newform `+ q·`(weight 4) | FINITE-EXACT (4b) |
| `K_4` plus a degree-2 vertex | 5, 8 | contains the circuit `K_4` | — | `−b_q` (the `K_4` form), 0 for tournaments | FINITE-EXACT |
| octahedron minus 2 edges at one vertex | 6, 10 | contains the circuit `W_4` | — | `q + a_q` (the `W_4` form) | FINITE-EXACT |
| octahedron `K_(2,2,2)` | 6, 12 | over-braced | general type | `(5 − q)a_q + 3q`: the `W_4` form, twisted | FINITE-EXACT |

**Reading (HYP-9278, after the cube).**
- New non-Tate motives appear at rigidity circuits, with motivic weight `|V| − 2`.
- Graphs containing a smaller circuit reduce to that circuit's motive, with Tate twists.
- For plane graphs, a graph also carries the motives of the circuits of its dual (section 4d).
- The first formulation, "Laman-sparse pure cores are Tate", was registered as a prediction before the cube was computed. **It was refuted.** The cube is Laman-sparse, yet it carries the `W_4` form, which it inherits from its dual, the octahedron.
- "Least resistance" is the reducible regime, the series-parallel graphs. It is not the whole Fano side.

The dictionary is not proved: the arrangement singularities need crepant resolutions, and the weight assignment is only checked on the rows above.

## 4d. Planar Fourier duality: the Calabi–Yau line is self-dual (PROVED; FINITE-EXACT checks)

**Proposition (PROVED).** Let `G` be a connected plane graph with `|E|` even, and give the dual edge of `u → v` the orientation crossing it from the face on its left to the face on its right. Then

`U(G) = χ(−1)^{|E|/2}·q^{(|V| − |V*|)/2}·U(G*)`.

*Proof.*
1. Write `χ(d) = q^{−1}Σ_y χ̂(y)ψ(yd)`. The vertex sum then becomes `q^{|V| − |E|}Σ_{flows y}Π_e χ̂(y_e)`.
2. On a plane graph, flows are face potentials modulo constants, which contributes a factor `1/q`.
3. `χ̂ = χ(−1)G·χ` with `G² = χ(−1)q`, and Euler's formula gives `|V*| = |E| − |V| + 2`. ∎

**Checks** (`paley_planar_fourier_duality_20261014.py`):
- `U(Q_3) = q·U(K_(2,2,2))` at all 32 primes of the cube run;
- `K_6 − P_6` (planar and self-dual) at `q = 5, 7, 11, 13`;
- `W_4` and `W_5` at several primes each;
- three random plane graphs.

**Consequences.**
- Duality preserves the Calabi–Yau count `|E| = 2|V| − 2`, which for plane graphs is `|F| = |V|`.
- It exchanges the Fano side (`|E| < 2|V| − 2`) with general type. This is a mirror-like involution on the graph side.
- Every circuit found so far to carry a new motive is self-dual or has a self-dual count: `K_4`, `W_4`, `W_6` and `K_6 − P_6` are all self-dual plane graphs.
- The quadratic character is a Fourier eigenfunction, `χ̂ = χ(−1)G·χ`, and the Gauss sum's phase is `1` or `i`. This is the precise "squaring" (`G² = χ(−1)q`) and the precise mod-4/mod-8 phase in the owner's "doubling/squaring and mod 8" prompt. Its Bott-periodicity reading (Weil representation, Maslov index mod 8) is ANALOGY.



## 5. The owner's numbers, typed

| item | status |
|---|---|
| `|Aff(P_q)| = q(q−1)/2 = T_{q−1}` (triangular) | EXACT |
| `T_{q−1}` is Fibonacci iff `q − 1 ∈ {1, 2, 6, 10}` (Luo 1989: the triangular Fibonacci numbers are 1, 3, 21, 55); with `q ≡ 3 (mod 4)` a prime power: `q = 3, 7, 11` | EXACT |
| `q = 3, 7, 11` are exactly the orders where the maximum `H` is known (A038375 to `n = 11`) and attained by Paley | EXACT; the link to Fibonacci is NUMEROLOGY |
| `H(P_q) ≡ T_{q−1} (mod 8)` at all 7 known `q`, i.e. 3, 5, 7 at `q = 3, 7, 11` and 1 at `q = 31`; period 16 in `q` | EXACT given HYP-9273 (B) |
| `A_8` roots 3, 7 sum to `10 = 7 + 3` (the repository's `F = U + S` for `p²qr`) and multiply to `21 = |Aut(P_7)| = T_6 = F_8` | roots EXACT (too few points); the rest NUMEROLOGY |
| `A_6, A_8` non-integer roots `29/5 = 6 − 1/5` and `55/7 = 8 − 1/7` (`55 = T_10 = F_10`) | EXACT roots; the pattern `2k − 1/(2k−1)` stops being testable at `A_10`, which is not a polynomial |
| the `a_q` coefficient 34 in `A_10` (`= F_9`); `log₂ 3 = [1; 1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55, …]` | EXACT values; NUMEROLOGY links |
| "mod 8 Bott periodicity" | EXACT mod-8 facts: `H ≡ T_{q−1}`; `χ_8 = (2/q)` decides whether doubling `x ↦ 2x` is an automorphism (`q ≡ 7 (mod 8)`) or an anti-automorphism (`q ≡ 3`); the `Y_q ≡ 1 (mod 128)` signature. The Clifford/K-theory periodicity itself is ANALOGY only |
| "noncommutative unbalanced synthesis of multiplication and addition" | EXACT reading: `AGL(1,q) = F_q ⋊ F_q^×` acts on the Hamiltonian paths of `P_q`. The tournament splits it into automorphisms (square multipliers) and anti-automorphisms (non-squares), and `H = |AGL|·F + (|AGL|/2)·P` (S33 split) |

## 6. Mahler 3/2, the 2–3 problems, the attached paper, functions as graphs

**Mahler and the 2–3 problems (CITED + cross-reference).**
- Tiling-modular-atoms' same-day packets carry the Mahler side: [the carry interface](mahler_carry_interface_20261014i.md) and [the triangular quotient of Mahler's ceiling map](recursion_carry_bridges_20261014i.md). Both build on THM-2228.
- Related 2–3 problems: Collatz cycles; Mahler's Z-numbers (Mahler 1968; Flatto–Lagarias–Pollington 1995, `limsup − liminf >= 1/p` for `ξ(p/q)^n`); Waring's `g(k)` through `{(3/2)^k}`; Erdős' ternary digits of `2^n`; Furstenberg's `×2 ×3`; Catalan (8, 9); Pillai–Tijdeman.
- All of them need control of `k log₂ 3 mod 1` beyond what linear forms in logarithms give.
- The continued fraction of `log₂ 3` has the large partial quotient 55 at index 14 (between `301994/190537` and `16785921/10590737`), the source of very close approaches `2^A ≈ 3^l` at that scale. That it is `T_10 = F_10` is NUMEROLOGY.

**The attached paper (M3; PREPRINT CLAIM; read in S25).**
- It claims ordinary two-point Chowla for Liouville with a log saving.
- Over `F_q`, the quadratic character is the analogue of Liouville. The two-point correlation is an exact Jacobi sum, and every higher pattern correlation splits by the compiler of section 1 into a polynomial main term plus core terms. The first core term is a Deligne-bounded modular coefficient, the wheel's Beukers form.
- M3's apparatus (prohibited paths, closed words, centred traces) organises integer correlations by the same kind of words and closed walks. The finite-field version of this note shows what its "exceptional words" would be if a modular structure existed (ANALOGY).

**Functions as directed graphs.** Every quantity in this note is a count of graph homomorphisms of a path (or a union of paths) into the Paley Cayley digraph, decomposed over its quotients. That is the "functions as graphs" view made exact. The Collatz and aliquot functional graphs of S32/S33 differ: they are out-degree-1 graphs of a single map, and their quotients are cycles, where the arithmetic sits.

## 7. Hostile controls

1. **The electrical analogy** is exact for the combinatorics: the reduction steps and Duffin's characterisation. It is not a statement that some resistance equals a character sum.
2. **"One modular form through twelve steps"** is FINITE-EXACT. The core list is complete for these reduction rules, and the octahedron identity is checked at 26 primes, not proved. At larger `L` new cores will appear; `K_6 − P_6` already enters the weighted `R(q)` expansion.
3. **The `K_6 − P_6` identification** is mac-mini's (HYP-9276), and this session's agreement is a second pipeline, not a proof. The 2-adic congruence `|b(q)| ≡ ±(q+1) (mod 2^m)` is FINITE-EXACT; the `Z/2 × Z/8` mechanism is a proposal. The Mersenne and Fermat values are coincidences of a known shape: both factors being powers of 2.
4. **The rigidity–Calabi–Yau dictionary** is an organising hypothesis, not a theorem, and it has already lost one clause. "Laman-sparse pure cores are polynomial" was refuted by the cube, then repaired by planar duality. "Calabi–Yau" ignores the singularities of the arrangement double covers.
5. **The `T_6` identification** came from an exact overdetermined fit (3 nonzero coefficients among 8 candidates, at 13 primes) and was then confirmed at 39 primes by an independent evaluation route. `T_8` failed in the same basis, which is reported as such.
6. **The Catalan leading terms** are FINITE-EXACT for `k <= 6`, and the free-probability mechanism is heuristic.
7. **The owner's numerology** (Fibonacci, 10 = 7 + 3, 55) is recorded with its exact sources. No mechanism links it to the run sums beyond the too-few-points roots.
8. **Nothing here bears on Collatz, Hadwiger–Nelson or JC2 directly.** It locates precisely where the Riemann hypothesis for varieties acts on a family of discrete structures: at rigidity circuits.

## 8. Open threads

1. Prove the octahedron identity, presumably through the quartic-family second moment and Ahlgren–Ono.
2. Prove HYP-9276. Explain the congruence `|b(q)| ≡ ±(q+1) (mod 2^m)`, perhaps through McCarthy–Papanikolas' `Z/2 × Z/8` Legendre curves.
3. Identify `T_8`, the wheel `W_8` (weight 8, level 16 or 32?), and test the dictionary on 7-vertex circuits other than wheels.
4. `A_14`: Bell(14) = 190,899,322 partitions, out of reach for the Python compiler. A compiled version is needed to see whether `K_6 − P_6` enters single runs.
5. Derive the coefficients 34, `6688 − 408q` and `−68` combinatorially, as Möbius-weighted counts of core embeddings.

## 9. Audit

Pending.

## 10. Reproduction and references

**Reproduction.** Scripts in `04-computation/experiments/`, with `.out` files beside this note:
- `paley_run_partition_compiler_20261014.py` (`L = 2..10` in under a minute; `L = 12` in about 10 minutes);
- `paley_A12_assembly_20261014.py`;
- `paley_six_vertex_cores_20261014.py`;
- `paley_core_symmetry_20261014.py`;
- `legendre_twisted_moment_apery_20261014.py`;
- `paley_k6p6_fast_table_20261014.py` (`Y_q` exactly at all primes to 691, 70 s);
- `paley_k6p6_level32_form_20261014.py` (the level-32 newform from eta quotients; the square law at 124 primes);
- `paley_k6p6_two_adic_mersenne_20261014.py` (FFT evaluator after THM-4653; square law and 2-adic depths to 200,000; Mersenne primes);
- `paley_wheel_ladder_20261014.py`, `paley_wheel_ladder_fit_20261014.py`, `paley_wheel_jacobi_20261014.py` (wheel traces, fits, Jacobi-sum confirmation to 400);
- `paley_cube_fano_core_20261014.py` (the cube, `O(q⁴)`, cross-checked by brute force; 32 primes);
- `paley_planar_fourier_duality_20261014.py` (duality checks);
- `paley_k6p6_two_adic_valuation_law_20261014.py` (`|b(q)|` at every prime `q ≡ 3 (mod 4)` below 200,000; the valuation law).
- `paley_rigidity_dictionary_tests_20261014.py` (the Fano-regime and circuit-reduction rows of section 4c).

**References.**
- **Literature:**
  - R. J. Duffin, Topology of series-parallel networks, *J. Math. Anal. Appl.* 10 (1965);
  - S. Ahlgren and K. Ono, *J. reine angew. Math.* 518 (2000);
  - F. Beukers, *J. Number Theory* 25 (1987) (supercongruence), and *Astérisque* 147–148 (1987) (the modular irrationality proof, level 6);
  - P. Deligne (1974);
  - M. Luo, On triangular Fibonacci numbers, *Fibonacci Quart.* 27 (1989);
  - R. J. Evans, J. R. Pulham and J. Sheehan (1981);
  - K. Mahler, *J. Austral. Math. Soc.* 8 (1968);
  - D. McCarthy and M. A. Papanikolas, A finite field hypergeometric function associated to eigenvalues of a Siegel eigenform, *Int. J. Number Theory* (2015), arXiv:1205.1006 (the newform `f_2`, via HYP-9276);
  - G. Laman, On graphs and rigidity of plane skeletal structures, *J. Engrg. Math.* 4 (1970); L. Henneberg (1911);
  - C. Meyer, *Modular Calabi–Yau threefolds* (Fields Inst. Monographs 22, 2005), double octics;
  - L. Flatto, J. C. Lagarias and A. D. Pollington, *Acta Arith.* 70 (1995);
  - OEIS A038375.
- **This repository:**
  - THM-4652;
  - mac-mini's THM-4653 and HYP-9276 (the `K_6 − P_6` identification);
  - HYP-9274, HYP-9275, HYP-9278;
  - S34 [note](paley_h31_exact_rate_and_apery_wheel_20261014.md);
  - codex's [wheel](wheel_apery_extension_20261014h.md) and [run](wheel_run_interface_20261014h.md) packets;
  - tiling-modular-atoms' [next cores](paley_next_cores_20261014i.md), [Mahler carries](mahler_carry_interface_20261014i.md) and [Mahler triangular quotient](recursion_carry_bridges_20261014i.md);
  - THM-2228;
  - S31, S33, S25 notes.
