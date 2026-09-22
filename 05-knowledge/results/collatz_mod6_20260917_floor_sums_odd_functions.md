# The odd-function floor law: why the three floor sums see squarefreeness, primality, and (for even powers) class numbers

**Status: PROVED elementary identities (S1-S4, S6a) with FINITE-EXACT controls; CITED Dirichlet class-number formula with independent finite verification (S5); one OPEN table (quartic deviation D_4); one SCOPE verdict (no object-level map from floor sums to Redei "odd cycles" or to 3n+k cycle counts). No novelty claim for the squarefree / prime characterizations, the fixed-layer formula, or the cube-root rectangle reciprocity: all of these are already PROVED in [arithmetic_braids2_20260917_floor_reciprocity.md](arithmetic_braids2_20260917_floor_reciprocity.md) (its equations (3)-(5), (8)-(10), (12)); what this lane adds is the named odd-function law with the fixed layer `Z_f(n)`, the Pillai gcd-sum form of the bilinear sum with the multiplicative inequality `P(n) >= 2n-1`, the corollary `4N3 - (N1+N2) = 2(P-2n+1) - Z_3`, and the even-power class-number statement. Nothing here bears on Collatz convergence or cycle uniqueness.** Session collatz-mod6-20260917 (mac-mini), wave two, lane `floor_sums_odd_functions`; script recovered 2026-09-21 and re-run unchanged (all checks pass).

## Inheritance and concept board

The three sums and their finite truth sets (`n<400`, `n<120`, `n<200`), the relation `2n*sum floor(k^e/n) = 2 sum k^e - n(n-1) + n Z_e(n)` (odd `e`, `n<150`), `Z_7(n) = n/rad(n)-1` (`n<200`) and the even-exponent failure at `(p,e) = (3,2)` are the session lead's firsthand fact F2; the `3n+k` cycle counts `4,9,5,7,13` for `k = 1,5,7,11,13` are the lead's F3 (re-censused independently in a scratch script on a bounded seed range, not part of the lane output; the counts agree).

**Closest proved mechanism.** [arithmetic_braids2_20260917_floor_reciprocity.md](arithmetic_braids2_20260917_floor_reciprocity.md), Section 2, equation (6): for a finite set with an involution and `h(iota x) = -h(x) (mod n)`, `sum floor(h/n) = (sum h)/n - (|X|-Z)/2`. Section 1 below is that lemma with `X = {1,..,n-1}`, `iota(k) = n-k`, `h = f` an odd polynomial, and the fixed layer named `Z_f(n)`; the braids2 note already remarks that (6) covers odd integer polynomials. The same note proves (its (4), (12)) `R_d(n) = prod p^{a - ceil(a/d)}` with `R_d = n/rad(n)` iff `d >= max a`, which is Section 2 here; its (7)-(9) is the rectangle reciprocity of Section 3 for every degree; its (3), (5) is Section 4's `N3 = (n-2)(n-1)^2/4 + Z2(n)/2` with `Z2(n) = sum_{i=1}^{n-1}(gcd(i,n)-1)`, and it proves "(3) iff prime" by a zero-divisor construction. Behind both notes stands the retained-fibre principle of [THM-2422](../../01-canon/theorems/THM-2422-operation-fibres-summand-closure-and-twin-center-ancestry.md): the exceptional (fixed / boundary) fibre must be kept, never averaged away.

**Canonical hostile.** `n = 4`: half-integer baselines with correction `1/2` in both the cubic and bilinear sums; and `(f, n) = (x^2, 3)`, the minimal failure of the odd law for an even function (floor sum `1` versus law value `2/3`).

**Corrected near misses.** (i) "the correction for even powers is `h(-p)/2`": it is exactly `h(-p)` for `p = 3 (mod 4)`, `p > 3` (Section 5). (ii) "`Z_7(n) = n/rad(n) - 1` for all `n`": true iff every prime exponent of `n` is at most `7`, first failure `n = 256` where `Z_7 = 63` and `n/rad(n) - 1 = 255` (braids2 (12), re-verified). (iii) The braids2 near miss "cubes permute the nonzero residues mod `p`" is not used anywhere below: only the antipodal symmetry `r_k + r_{n-k} in {0, n}` is.

**Least-used sidecar.** The even-modulus fixed point `k = n/2`, whose residue lies in `{0, n/2}` and therefore contributes `(n/2)[r != 0]` exactly, so even `n` needs no separate case; and the multiplicativity of `P(n)/n = sum_{d|n} phi(d)/d`, which gives the inequality `P(n) >= 2n-1` with equality iff `n` is prime.

**Signed Collatz inheritance.** Negation conjugacy `C_+(-n) = -C_-(n)` and the three positive `3n-1` cycles (the fixed point `1`, the 2-cycle `(5,7)`, and a 7-cycle) are [arithmetic_braids_20260917_summand.md](arithmetic_braids_20260917_summand.md) Section 7 and [arithmetic_braids_20260917_collatz.md](arithmetic_braids_20260917_collatz.md) (B7 and the four `3n+1` cycles on `Z` including `(-5,-7)`, FINITE-EXACT bounded-word). The nine cycles of `3n-5` in [arithmetic_braids2_20260917_signed_cycles.md](arithmetic_braids2_20260917_signed_cycles.md) are, under negation, the `9` cycles of `3n+5` in F3. The repo meaning of "odd cycles" is Redei parity: [THM-001](../../01-canon/theorems/THM-001-redei.md), [LEM-020](../../01-canon/theorems/LEM-020-redei-involution-parity-layers.md), [THM-1745](../../01-canon/theorems/THM-1745-leaf-graded-arborescence-filtration-721-shadow.md) (h-spectrum = odd numbers minus `{7,21}`). SCOPE: [THM-3341](../../01-canon/theorems/THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors.md) (Gaussian squaring of triples, Pell hypotenuses) and THM-4139 / THM-4146 (the `x^2 - 29/16` three-cycle) were read for this lane; no map found from any floor sum here to squares of Gaussian integers, Pell orbits, or quadratic three-cycles.

| Lane | Object | Question | Hostile / retained coordinate |
|---|---|---|---|
| Anchor | odd `f`, residues `r_k` of `f(k)` mod `n` | exact value of `sum floor(f(k)/n)` | fixed point `k = n/2`; fixed layer `Z_f(n)` |
| Niche | `x^e`, `n_e = prod p^{ceil(a/e)}` | when is the fixed layer empty? | `e = 1` degenerate; `n = 256` for `e = 7` |
| Bridge | bilinear `ij` | which defect sees primality? | `n = 6`: `Z_3 = 0` but `P - 2n + 1 = 4` |
| Wildcard | even `e` mod prime `p` | what replaces the pairing? | `p = 3` (`w = 6`); `p = 5 (mod 8)` quartic table |

Script: `04-computation/experiments/collatz_mod6_20260917_floor_sums_odd_functions.py`; regenerated output `05-knowledge/results/collatz_mod6_20260917_floor_sums_odd_functions.out`. All arithmetic is exact (integers / `Fraction`); cube roots are monotone integer walks; class numbers are counted from reduced forms; every check is an explicit `raise`, active under `python3 -O`.

## 1. The odd-function floor law (PROVED; specialization of braids2 (6))

**Theorem S1.** Let `f in Z[x]` with `f(-x) = -f(x)`, `n >= 2`, and `Z_f(n) = #{1 <= k <= n-1 : n | f(k)}`. Then

    sum_{k=1}^{n-1} floor(f(k)/n) = (1/n) sum_{k=1}^{n-1} f(k) - (n-1)/2 + Z_f(n)/2.

*Proof.* Write `f(k) = n q_k + r_k`, `0 <= r_k < n`, so `sum floor(f(k)/n) = (sum f(k) - sum r_k)/n`. Since `f` is an odd polynomial, `f(n-k) = f(-k) = -f(k) (mod n)`, so `r_k + r_{n-k} = 0 (mod n)`, i.e. `r_k + r_{n-k} in {0, n}`, and the value `0` occurs iff `r_k = 0` iff `r_{n-k} = 0`. If `n` is even the involution `k -> n-k` has the single fixed point `k = n/2`, where `2 r_{n/2} = 0 (mod n)` forces `r_{n/2} in {0, n/2}`, i.e. `r_{n/2} = (n/2)[r_{n/2} != 0]`. Hence every `k` with nonzero residue contributes exactly `n/2` on average over its orbit, and

    sum_k r_k = (n/2) * #{k : r_k != 0} = (n/2)(n - 1 - Z_f(n)).

Divide by `n`. QED. Only `f(-x) = -f(x) (mod n)` was used, so the law holds for every odd function `Z/n -> Z/n`; this is braids2 equation (6) with `|X| = n-1`.

FINITE-EXACT: verified for `f = x, x^3, x^5, x^3+x, 3x, x^7-5x^3+x` and all `2 <= n <= 400`; residue pairs and fixed-point residues printed for `n = 4, 6, 10, 12` (e.g. `n = 6`, `f = x^3`: pairs `(1,1,5), (2,2,4)`, fixed residue `3 in {0,3}`; `n = 12`: fixed residue `0 in {0,6}`).

**REFUTED for even `f`** (minimal witness): `(x^2, n = 3)`, floor sum `1`, law value `2/3`. Also `x^4` (`5` vs `14/3`), `x^2 + x` (`2` vs `13/6`) and the constant `2` (`0` vs `1/3`) all first fail at `n = 3`; `n = 2` holds for every `f` because the only term `k = 1` is the fixed point. This is the lead's `(p,e) = (3,2)` and braids2's "even-degree antipodal pairing is false".

## 2. The fixed layer `Z_{x^e}(n) = n/n_e - 1` and the squarefree truth set (PROVED; no novelty claim; density CITED)

**Theorem S2.** For `e >= 1` put `n_e = prod_{p^a || n} p^{ceil(a/e)}`. Then `Z_{x^e}(n) = n/n_e - 1`.
*Proof.* `n | k^e` iff `a <= e v_p(k)` for every `p^a || n` iff `v_p(k) >= ceil(a/e)` iff `n_e | k`; `n_e | n`, so the multiples of `n_e` in `[1, n-1]` number `n/n_e - 1`. QED. (This is `R_e(n) - 1` with braids2 (4).)

**Corollaries.** (a) `Z_e(n) = 0` iff `e = 1` or `n` squarefree. (b) Since `(sum k^3)/n - (n-1)/2 = (n-1)(n-2)(n+1)/4`,

    sum_{k=1}^{n-1} floor(k^3/n) = (n-2)(n-1)(n+1)/4 + Z_3(n)/2,

so the user's identity (1) holds **iff `n` is squarefree** (braids2 (1)). (c) For every odd `e >= 3`, `sum floor(k^e/n) = (sum k^e)/n - (n-1)/2` iff `n` is squarefree; for `e = 1` it holds for all `n` (degenerate member; braids2 (10)). (d) `Z_7(n) = n/rad(n) - 1` iff every exponent `a <= 7`; first failure `n = 256`, `Z_7(256) = 63` versus `255` (braids2 (12)).

CITED (Gegenbauer 1885; Hardy-Wright, Theorem 333): the squarefree integers have natural density `1/zeta(2) = 6/pi^2`. That is the entire content of the user's `6/pi^2` remark: the truth set of (1), and of the odd-`e` identity for each fixed odd `e >= 3`, has density `6/pi^2`.

Typed analogy: source = truth set of (1) (and of every odd-`e` identity); target = squarefree integers; map = `n -> Z_e(n) = n/n_e - 1` (identity holds iff `Z_e(n) = 0`); preserved = "squarefree", both directions (PROVED); lost = nothing (the deviation `Z_e/2` is itself the exact local invariant); sidecar = the exponent profile `(a_p)` through `n_e`; decisive test `n = 4, 9, 12` (deviations `1/2, 1, 1/2` for `e = 3`). No Collatz input enters this map.

FINITE-EXACT: `Z` formula for `1 <= e <= 9`, `2 <= n <= 2000`; truth set of (1) on `[2, 2000]` = the `1214` squarefree numbers exactly; `e = 3, 5, 7, 9` on `[2, 1200]` iff squarefree; `e = 1` everywhere on `[2, 1200]`; `Z_7` first failure `256`; `#squarefree <= 10^6 = 607926` against `6/pi^2 * 10^6 = 607927.1` (ratio `0.607926`).

| n | floor sum (e=3) | RHS(1) | Z_3(n) | deviation = Z_3/2 |
|---|---|---|---|---|
| 4 | 8 | 15/2 | 1 | 1/2 |
| 8 | 96 | 189/2 | 3 | 3/2 |
| 9 | 141 | 140 | 2 | 1 |
| 12 | 358 | 715/2 | 1 | 1/2 |
| 36 | 11010 | 22015/2 | 5 | 5/2 |

## 3. Cube-root lattice reciprocity (PROVED; no novelty claim: braids2 (7)-(9))

**Theorem S3.** Let `Y = floor((n-1)^3/n)`, `N1 = sum_{x=1}^{n-1} floor(x^3/n)`, `N2 = sum_{y=1}^{Y} floor((ny)^{1/3})`. Then (i) `Y = (n-1)(n-2)`; (ii) `N1 + N2 = (n-1)^2 (n-2) + Z_3(n)`; (iii) `N2 = (3n-5)(n-2)(n-1)/4 + Z_3(n)/2`, so identity (2) holds iff `n` is squarefree.

*Proof.* (i) `(n-1)^3 = n(n^2 - 3n + 3) - 1`. (ii) In the box `[1, n-1] x [1, Y]`, `N1` counts the points with `ny <= x^3` (all inside the box since `floor(x^3/n) <= Y`) and `N2` counts the points with `x^3 <= ny` (inside since `ny <= (n-1)^3` gives `floor((ny)^{1/3}) <= n-1`). Every box point lies in at least one set; the overlap is the curve `x^3 = ny`, whose box points are exactly the `x in [1, n-1]` with `n | x^3` (then `y = x^3/n in [1, Y]` automatically), i.e. `Z_3(n)` of them. So `N1 + N2 = (n-1)Y + Z_3(n)`. (iii) Subtract Section 2(b): `N2 = (n-1)(n-2)(4(n-1) - (n+1))/4 + Z_3/2`. QED.

Cube roots are exact: `(ny)^{1/3}` is nondecreasing in `y`, so walk `x` upward while `(x+1)^3 <= ny`. FINITE-EXACT: (ii) and (iii) for all `2 <= n <= 220` and for `n in {225, 243, 250, 256, 289, 300, 331, 343, 360, 361, 385, 399, 400}`; truth set of (2) on that universe = squarefree exactly (covers the lead's `n < 120`).

| n | N1 | N2 | N1+N2 | (n-1)^2(n-2) | Z_3 | N2 - RHS(2) |
|---|---|---|---|---|---|---|
| 7 | 60 | 120 | 180 | 180 | 0 | 0 |
| 8 | 96 | 201 | 297 | 294 | 3 | 3/2 |
| 9 | 141 | 309 | 450 | 448 | 2 | 1 |
| 12 | 358 | 853 | 1211 | 1210 | 1 | 1/2 |
| 31 | 6960 | 19140 | 26100 | 26100 | 0 | 0 |
| 32 | 7676 | 21161 | 28837 | 28830 | 7 | 7/2 |

## 4. The bilinear sum via Pillai's gcd-sum (PROVED; the Pillai form and the inequality are this lane's, "(3) iff prime" is braids2 (5))

**Theorem S4.** Let `P(n) = sum_{d|n} d phi(n/d) = sum_{k=1}^{n} gcd(k, n)` (Pillai's function, A018804). Then

    S(n)  := sum_{i,j=1}^{n-1} (ij mod n) = n (n^2 - P(n))/2,
    N3(n) := sum_{i,j=1}^{n-1} floor(ij/n) = ((n-1)^2 n^2/4 - S(n))/n = (n-1)^2 (n-2)/4 + (P(n) - 2n + 1)/2,

and `P(n) >= 2n - 1` with equality iff `n` is prime; hence identity (3) holds **iff `n` is prime**.

*Proof.* (i) Extend `i, j` to `[0, n-1]` (added terms vanish). Fix `i`, put `d = gcd(i, n)`, `i = d i'` with `gcd(i', n/d) = 1`. Then `ij mod n = d (i' j mod n/d)`, and as `j` runs over `[0, n-1]` the value `i' j mod n/d` takes each residue in `[0, n/d - 1]` exactly `d` times, so `sum_j (ij mod n) = d^2 (n/d)(n/d - 1)/2 = n(n-d)/2`. There are `phi(n/d)` values of `i` with `gcd(i, n) = d`, so `S(n) = (n/2) sum_{d|n} phi(n/d)(n - d) = (n/2)(n^2 - P(n))` using `sum_{d|n} phi(n/d) = n`. (ii) `sum_{i,j} ij = (n(n-1)/2)^2`; subtract `S` and divide by `n`; simplify with `(n-1)^2(n-2) = n^3 - 4n^2 + 5n - 2`. (iii) `P(n)/n = sum_{d|n} phi(d)/d` is multiplicative with value `1 + a(1 - 1/p)` at `p^a`. For `n = p`: `P = 2p - 1`. For `n = p^a`, `a >= 2`: `P/n = 1 + a - a/p >= 1 + a/2 >= 2 > 2 - 1/n`. For `n` with two distinct prime factors: `P/n >= (2 - 1/p)(2 - 1/q) >= (3/2)(5/3) = 5/2 > 2`. So `P(n) > 2n - 1` unless `n` is prime, and identity (3), i.e. `P(n) = 2n - 1`, holds iff `n` is prime. QED. For prime `n`, (i) is the classical permutation argument `j -> ij mod n`.

Relation to braids2 (5): `Z2(n) = sum_{i=1}^{n-1}(gcd(i, n) - 1) = P(n) - 2n + 1`, so `N3 = (n-2)(n-1)^2/4 + Z2(n)/2` there and here agree; the new content is naming `P` and the multiplicative proof of the inequality.

**Corollary S4b (PROVED).** `4 N3(n) - (N1(n) + N2(n)) = 2(P(n) - 2n + 1) - Z_3(n)`; at primes both corrections vanish and `4 N3 = N1 + N2 = (n-1)^2 (n-2)`. The two corrections are independent: `n = 6` has `Z_3 = 0`, `P - 2n + 1 = 4`; `n = 8` has `Z_3 = 3`, `P - 2n + 1 = 5`. The bilinear sum sees the whole divisor lattice (through `P`); the cubic sums see only the exponent profile (through `n_3`). Braids2 phrases the same separation as "reduced ring versus field" with separator `n = 6`.

FINITE-EXACT: (i), (ii), (iii) and the gcd-sum identity for all `2 <= n <= 260` (truth set of (3) = primes exactly, covering the lead's `n < 200`); S4b for all `2 <= n < 150`.

| n | P(n) | P-2n+1 | N3 | (n-1)^2(n-2)/4 | deviation |
|---|---|---|---|---|---|
| 4 | 8 | 1 | 5 | 9/2 | 1/2 |
| 5 | 9 | 0 | 12 | 12 | 0 |
| 6 | 15 | 4 | 27 | 25 | 2 |
| 7 | 13 | 0 | 45 | 45 | 0 |
| 8 | 20 | 5 | 76 | 147/2 | 5/2 |
| 9 | 21 | 4 | 114 | 112 | 2 |
| 12 | 40 | 17 | 311 | 605/2 | 17/2 |
| 15 | 45 | 16 | 645 | 637 | 8 |

## 5. Even powers see class numbers (CITED Dirichlet + PROVED reduction + FINITE-EXACT p < 500)

**Setup (PROVED).** For a prime `p` and `e >= 1` let `g = gcd(e, p-1)`. The map `k -> k^e` sends `(Z/p)^*` onto the unique subgroup `H_g` of index `g`, `g`-to-one, so `sum_{k=1}^{p-1} (k^e mod p) = g sum_{r in H_g} r`. If `-1 in H_g` (iff `(p-1)/g` is even) then `H_g` is symmetric under `r -> p-r`, `sum_{H_g} r = p(p-1)/(2g)`, and the odd-law value `(sum k^e)/p - (p-1)/2` is exact; this recovers Section 1 for odd `e` (then `g` is odd and `(p-1)/g` even). If `(p-1)/g` is odd (forcing `g` even) there is a genuine deviation

    sum_{k=1}^{p-1} floor(k^e/p) = (sum k^e)/p - (p-1)/2 + D_e(p),   D_e(p) = (p-1)/2 - g sum_{H_g} r / p.

**Quadratic case `g = 2`, `p = 3 (mod 4)`.** Let `R`, `N` be the sums of the quadratic residues / non-residues in `[1, p-1]`; `R + N = p(p-1)/2`. CITED (Dirichlet's class number formula for `Q(sqrt(-p))`; Davenport, *Multiplicative Number Theory*, ch. 6): `R - N = sum_a (a/p) a = -p h(-p)` for `p > 3` (unit count `w = 2`). Hence `R = p(p-1)/4 - p h(-p)/2` and

    sum_{k=1}^{p-1} floor(k^2/p) = (p-1)(2p-1)/6 - (p-1)/2 + h(-p)      (p = 3 (mod 4), p > 3).

The deviation from the odd law is **exactly `h(-p)`**, not `h(-p)/2`. **REFUTED at `p = 3`** (`w = 6`): the displayed formula would give `5/3` while the actual sum is `1` (formula error `-2/3`); the odd-law deviation there is `1/3` with `h(-3) = 1`, consistent with the general Dirichlet form `R - N = -(2p/w) h(-p)` at `w = 6` (that `w`-form is the same CITED source; only the `p = 3` value `1/3` is verified in the lane output). For `p = 1 (mod 4)` the deviation is `0`: the even power sees nothing. For general even `e` (PROVED from the setup): `D_e(p) = h(-p)` iff `gcd(e, p-1) = 2` and `p = 3 (mod 4)`; `D_e(p) = 0` whenever `(p-1)/g` is even.

`h(-p)` is computed independently by counting reduced primitive forms `(a, b, c)` with `b^2 - 4ac = -p`, so the Dirichlet recollection is itself checked on the finite universe. FINITE-EXACT: for all `50` primes `p = 3 (mod 4)` below `500` with `p > 3`, the reduced-form count equals `-(1/p) sum (a/p) a`, `R - N = -p h`, and the floor-sum deviation equals `h(-p)`; all primes `p = 1 (mod 4)` below `500` have deviation `0`; for `e in {2, 4, 6, 8, 10, 12}` and all primes `3 <= p < 500` both cases of the general statement hold.

| p | sum floor(k^2/p) | odd-law value | deviation | h(-p) |
|---|---|---|---|---|
| 3 | 1 | 2/3 | 1/3 | 1 |
| 7 | 11 | 10 | 1 | 1 |
| 11 | 31 | 30 | 1 | 1 |
| 19 | 103 | 102 | 1 | 1 |
| 23 | 157 | 154 | 3 | 3 |
| 31 | 293 | 290 | 3 | 3 |
| 43 | 575 | 574 | 1 | 1 |
| 47 | 695 | 690 | 5 | 5 |
| 59 | 1105 | 1102 | 3 | 3 |
| 67 | 1431 | 1430 | 1 | 1 |
| 71 | 1617 | 1610 | 7 | 7 |
| 79 | 2007 | 2002 | 5 | 5 |
| 83 | 2217 | 2214 | 3 | 3 |
| 103 | 3439 | 3434 | 5 | 5 |
| 107 | 3713 | 3710 | 3 | 3 |
| 127 | 5255 | 5250 | 5 | 5 |

**Precise statement.** Odd powers see only squarefreeness (Section 2); even powers see class numbers. **OPEN:** the quartic deviation `D_4(p)` for `p = 5 (mod 8)` (`g = 4`, `(p-1)/4` odd) is tabulated only: `6/5, 2, -2, 2, -2, 2, -6, 10, 6, 2, -6, 14, -2` for `p = 5, 13, 29, 37, 53, 61, 101, 109, 149, 157, 173, 181, 197`; it equals neither `h(-4p)` (`2, 2, 6, 2, 6, 6, 14, 6, 14, 6, 14, 10, 10`) nor `2h(-4p)` on all `p = 5 (mod 8)` below `500` (both verdicts `False`). No formula is claimed.

## 6. Typing "odd cycles ~ odd functions ~ odd powers"

| slogan word | object | status | map / verdict |
|---|---|---|---|
| odd functions | negation conjugacy `T_{-k}(-n) = -T_k(n)`, `T_k(n) = (3n+k)/2^{v_2(3n+k)}` | PROVED (`v_2(-m) = v_2(m)`, `-(3n+k) = 3(-n) + (-k)`); inherited from summand Section 7 (`C_+(-n) = -C_-(n)`) | source: S1 pairing `f(n-k) = -f(k) (mod n)`; target: `T_k` on odd integers; map: both are "the map commutes with `x -> -x`" (mod `n` in S1, on `Z` for Collatz); preserved: orbit structure up to global sign; lost: nothing; sidecar: the sign of `k`; test: the `3n-1` cycle `(5,7)` is the `3n+1` cycle `(-5,-7)` (`T_1(-5) = -7`, `T_1(-7) = -5`). FINITE-EXACT for `k in {+-1, +-5, +-7, +-13}` and odd `|n| <= 10001`. |
| odd powers | `x^e` is odd iff `e` is odd | PROVED (S1, S2) | literally "the exponent is odd, so `x^e` is an odd function"; even `e` destroys the pairing (`r_{n-k} = r_k`) and the residue mass becomes a character sum, hence `h(-p)` (Section 5). Nothing deeper. |
| odd cycles | Redei THM-001 (odd number of Hamiltonian paths), LEM-020 (involution with fixed layer, parity law), THM-1745 (`{7,21}` hole) | SCOPE (no object-level map found) | source: S1 residue multiset `{r_k}`; target: Redei witness sets; map: none between the objects. Shared **proof template** only: an involution (`k <-> n-k`; `tau <-> 1-tau`; path reversal in Route A of THM-001) whose non-fixed orbits cancel and whose fixed layer carries the residue. S1's cancellation is exact and integer-valued (`sum r_k = (n/2)(n-1-Z)`), strictly stronger than a parity law. Decisive test: `Z_3(8) = 3`, a fixed layer of size `3`, whereas every Redei fixed layer has size `1`; the `{7,21}` hole has no counterpart (the truth sets here are squarefree / prime, cofinite in nothing). |
| "three cycles" of `3n-1` | F3 cycle counts `4, 9, 5, 7, 13` for `k = 1, 5, 7, 11, 13` | no map found | `n = 4` and `n = 9` are both non-squarefree yet are the cycle counts for `k = 1` and `k = 5`; no predicate of Sections 2-4 survives. The `9` at `k = 5` is the negation image of the nine `b = -5` cycles of braids2 `signed_cycles`; the `4` at `k = 1` is the braids collatz note's bounded-word census of `3n+1` cycles on `Z`. |

Honest verdict: beyond the word "odd" (= commutes with negation), "odd cycles" shares only the involution proof template with "odd functions / odd powers".

## Reproduction

    python3 04-computation/experiments/collatz_mod6_20260917_floor_sums_odd_functions.py > 05-knowledge/results/collatz_mod6_20260917_floor_sums_odd_functions.out

About `9 s` (the output stamps `t=8.9s`), far below the RAM cap, ends with `ALL CHECKS PASSED`; `python3 -O` gives identical output modulo the timing stamps. The script recovered on 2026-09-21 passed unchanged; the regenerated output is identical (modulo timing) to the previously stored one.

## Stopping boundary / next question

Everything the user's three sums can say has been said, and most of it was already said in braids2 `floor_reciprocity`: one antipodal involution, one fixed-layer invariant `n/n_e - 1`, one gcd-sum `P(n)`, and for even exponents Dirichlet's `h(-p)`. This lane's additions are the named odd-function law, the Pillai form with `P(n) >= 2n-1` iff-prime, the corollary `4N3 - (N1+N2) = 2(P-2n+1) - Z_3`, and the class-number contrast. The single live question produced is the **quartic deviation** `D_4(p)`, `p = 5 (mod 8)`: it is a biquadratic-character sum, so the cheapest next test is to compare the table with the class number of `Q(i, sqrt(-p))` or with `h(-4p)` weighted by the representation `p = a^2 + b^2`. A second cheap extension is the composite-modulus even-power sum (small `n = pq` with both primes `3 (mod 4)`), where the deviation should decompose over `h(-d)` for `d | n`, `d = 3 (mod 4)`. Neither touches Collatz; the only rigorous Collatz bridge from this lane is negation conjugacy, which is inherited, and the braids2 dynamical hostiles (first loss and first gain of squarefreeness under one accelerated step) already show the floor defect is not monotone under the Collatz operation.
