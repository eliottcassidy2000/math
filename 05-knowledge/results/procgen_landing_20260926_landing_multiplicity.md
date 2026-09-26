# The landing multiplicity of the thin-divergence recursion: the worst case is exactly `ceil((k - D)/log_2 3)`, not `k`; the recursion is saturated at `a*(mu) = mu lambda*/h* - 3/2`; averaged bounds would suffice, but no orbit-blind count gives `mu < 1`

Lane `landing`, session `collatz-procgen-20260922`, 2026-09-26.
Scripts: `04-computation/experiments/procgen_landing_20260926_{lib,run}.py`, `procgen_landing_20260926_scan.c`.
Output: [procgen_landing_20260926.out](procgen_landing_20260926.out). It comes from one runner. Every `ok:` line is a `check(...)` that raises on failure, and the file ends with `ALL CHECKS PASSED`.

## Status

**Bottom line.** THM-4499's exponent `a* = lambda*/h* - 3/2 = -0.9862` is **not improved**. The landing multiplicity is now understood exactly in the worst case. The averaged version that the recursion could use is **OPEN** and out of reach of every orbit-blind argument. Nothing here bears on the existence of divergent orbits. Novelty is not claimed.

**PROVED (full hand proofs in sections 1-3).**

- **Theorem 1 (what the multiplicity costs, both directions).**
  - (i) Suppose the recursion's *averaged* multiplicity satisfies `#dippers <= C L^mu N(X 2^(-D)) + O(L^2)` for all large `X` and all depths `D` in the bootstrap range, with `L = log_2 X` and `0 <= mu <= 1`. Then `N(X) <= K X^(h*) L^a` for every `a > a*(mu) = mu lambda*/h* - 3/2`. For `mu = 0` this also holds with `a = -3/2` itself.
  - (ii) Saturation. `psi(X) = X^(h*) L^(a*(mu))` satisfies the recursion inequality `(R_mu)` at **every** depth `D`, with the no-dip count in THM-4499's form. So the recursion with multiplicity `L^mu` cannot give any exponent below `a*(mu)`.
  - THM-4499's remark is confirmed: multiplicity `L^(1/2)` gives `lambda*/(2h*) - 3/2 = -1.2431`, and `mu = 0` gives the ballot floor `-3/2`.
- **Theorem 2 (the worst case, exactly).** Setting: a positive segment of distinct integers, scale `X`, window `k = floor(log_2 X)`, depth `D = theta log_2 X >= 1`.
  - **Lemma S (shell).** Every dipper of a landing index `j` lies in the single dyadic shell `(2^D y_j, 2^(D+1) y_j]`, and the step into `j` is a halving.
  - **Lemma O.** Two dippers of `j` are separated by at least one odd letter.
  - **Bound.** Hence `m(j) <= ceil((k - D)/log_2 3)` for `b > 0`. For `b < 0` the bound is `ceil((k - D + log_2(3/2))/log_2 3)` when `y_j >= k|b|`.
  - **Consequence.** The pigeonhole constant `k` of THM-4476/4499 becomes `0.631 (k - D)`. This is a constant factor, so `a*` is unchanged.
- **Proposition C (climbs and W-bit hovers).**
  - An odd run contributes at most 2 dippers to any landing point.
  - A hover in a `W`-bit band followed by a crash spreads its dippers over at most `ceil(W) + 1` landing points.
- **Lemma A (depth averaging).** A pair `(i, j)` serves at most one integer depth, so `sum_D m_D(j) <= k - 1`. The best of `W` consecutive depths has average multiplicity at most `k/W`, but its no-dip set is paid at the deepest depth. Measured at that depth, the effective multiplicity drops by at most the factor `0.966` (best at `W = 2`), so the effect on `a*` is nil.
- **Proposition H (orbit-blind splits are linear).** For `b > 0`, `D >= 3` and every `M >= 2`, some residue class modulo `2^l` with `l <= 2 log_2 3 (M-1) + D + 2` consists of heavy dippers (at least `M` dippers on one landing point), apart from its `O(l)` smallest members.
  - Hence `#H_M(X) >= X 2^(-(2 log_2 3 (M-1) + D + 5))` for `X >= b l 2^(l+8)`.
  - Any split of the form "light landing points `<= M` each, heavy dippers counted as integers `<= X`" closes only if `M >= ((1-h*)L - D - a log_2 L - O(1))/(2 log_2 3) = Theta(L)`, which again gives `mu = 1`.
- **Corollary (answer to the brief).** Every approach the brief lists falls into one of three cases. A worst-case bound is sharp at `Theta(L)`. An orbit-blind or 2-adic average is linear by Proposition H and section 3.3. An average over depths is neutral by Lemma A. The approaches:
  - deterministic time;
  - hover-then-crash with strip and ballot counts;
  - two-place information;
  - depth averaging.

  The exponent `a*` is exactly what the method produces from its inputs. **No improved thin-divergence theorem results.** The only change is cosmetic: THM-4499's recursion holds with `ceil((k - D)/log_2 3)` in place of `k`.

**FINITE-EXACT.**

- **Exhaustive worst case.** Every `y <= 2^(k+1) - 1` was examined as a first dipper, with exact 128-bit comparisons. This covers 144 cells for `b = 1` (`k = 10..25`, nine half-integer depths `D` in `[1, 6]`), 144 cells for `b = -1`, and 78 cells each for `b = 5` and `b = -5` (`k <= 22`).
  - For `b = 1, 5, -1` the maximum equals `ceil((k-D)/log_2 3)` in **every** cell.
  - For `b = -5` every cell is within the `b < 0` bound.
  - There are zero violations of Lemmas S and O.
- **Constructed hostile integers.** Sturmian hover-then-halve integers with `k = 30, ..., 400` (`X` up to `2^401`), for `b = +-1` and three depths each, reach the bound or the bound minus 1 in all 36 cases (18 exactly).
- **Records.** All 90 delay records and all 61 path records up to `10^12` (OEIS A006877, A006884; 144 distinct starts) were checked at `L = 20, ..., 60` and three depths. Every landing point satisfies Lemmas S and O and the bound.
  - Actual Collatz orbits reach the worst case at small scales. The orbit of 13255 has a landing point with 12 dippers, which is exactly `ceil((20-2)/log_2 3)`, at `L = 20`, `D = 2`.
  - At the bootstrap's depth `theta_X` the maximum is 21 of 22 at `L = 40` (orbit of 4578853915) and 24 of 34 at `L = 60`. So the maximum is `0.40-0.53 L` on actual orbits.
  - The mean multiplicity over landing points is only 2.5-4.3, and the averaged ratio `#dippers/N(X 2^(-D))` is 0.13-0.60.
- **2-adic hostile word.** Twenty hover-crash-climb cycles attain the bound at every cycle (`23` at `k = 40`, `48` at `k = 80`). The integer realising the word (967 bits) shows no landing point inside the word at its own scale.

**CORRECTIONS** to the parallel note `collatz_landing_20260926_multiplicity_reassessment.md` and to HYP-9161. These are on `origin/main` from session opus S7 and were not visible in this worktree at launch.

1. "Residue classes realise multiplicity `k`" should read `ceil((k-D)/log_2 3) ~ 0.631 (k - D)` (Theorem 2 and the exhaustive scan).
2. **Climb-then-drop.** The claim that "every index of the climb lands at the same point" is false.
   - A climb step multiplies by more than `3/2`, so at most 2 climb points share a dyadic shell. Hence they share a landing point at most in pairs.
   - In the constructed example, 40 climb points spread over 24 landing points.
   - The probe maxima `0.25-0.35 L` do not come from the initial climbs. On the starts `2^L - 1` the climb lies above `X`, and `2^(L-1) + 1` has no climb. The maximal landing points take 0 or 1 dipper from the initial odd run; their dippers are short hovers later in the orbit (section 5.4).
3. **HYP-9161's hostile example.** A 4-bit hover of `0.2 L` steps does not give multiplicity `0.2 L` at one landing point. Its dippers split over at most 5 landing points (Proposition C), so the average is at least `0.04 L`.
   - This is still linear, so HYP-9161's substance (only the average can be polylogarithmic, and it is orbit-coupled) stands.

**OPEN.** HYP-9161 (averaged multiplicity `O(L^beta)` with `beta < 1`) is unchanged.
- In the shell language of Lemma S it asks for local times: for how many steps can one orbit revisit the dyadic shell lying `D` bits above its later landing points?
- It cannot follow from any orbit-blind count (Proposition H) or from 2-adic data (section 3.3).

**Not claimed.** Any Collatz statement; any bound on divergent orbits beyond THM-4499; novelty.

## 0. Setting

**The map and the constants.**

- `T_b(x) = x/2` for `x` even and `(3x+b)/2` for `x` odd, with `b` odd.
- `alpha = log_2 3`, `rho* = log_3 2`, `h* = h(rho*) = 0.9499555`.
- `lambda* = log_2(rho*/(1-rho*))/alpha = 0.488077` (THM-4487, THM-4499).

**The segment.** As in THM-4476, section 1.1, fix a positive one-signed segment `y_0, y_1, ..., y_M` (`M <= infinity`). Its terms are pairwise distinct positive integers with `y_(i+1) = T_b(y_i)`.

**Scale, window, depth.**

- The scale is `X >= 2`, with `L = log_2 X` and window `k = floor(L)`.
- The depth is `D = theta L` for real `D > 0`. THM-4476 writes the dip threshold as `y_i X^(-theta)`; this note writes it as `y_i 2^(-D)`.
- In the scripts `D` is rounded to a multiple of `1/64`, and all comparisons are exact.

**Dippers and landing points.**

- An index `i <= M - k` with `y_i <= X` is a **dipper** if `y_(i+s) < y_i 2^(-D)` for some `1 <= s <= k`.
- Its **landing index** is `j = i + s` for the least such `s`.
- The **landing multiplicity** `m(j) = m_(X,D)(j)` is the number of dippers with landing index `j`.
- `#Dip(X, D) = sum_j m(j)`, and `Lambda(X, D)` is the set of landing indices. Every landing index has `y_j < X 2^(-D)`.

**The recursion** (THM-4476 section 1.5, THM-4499 section 3). The indices with `y_i <= X` are:

- at most `k` end indices `(E)`;
- the dippers `(D)`;
- the no-dip indices `(ND)`, which inject into `F_b(X, theta)`.

So

```text
N(X) <= k + #Dip(X, D) + #F_b(X, theta),       #Dip(X, D) = sum_(j in Lambda) m(j),       (0.1)
#F_b(X, theta) <= 2|b| X^(0.585+theta) + 4 D_s X^(h*) L^(-3/2) 2^(lambda* D) e^(s D)      (THM-4499, Lemma 1.4c).
```

THM-4476 bounds `#Dip <= k N(X 2^(-D))`: each landing index is an orbit index below `X 2^(-D)`, and the pair (landing index, dip time) determines the dipper. **The landing multiplicity is the constant in this one inequality.**

## 1. What the multiplicity costs: Theorem 1

### 1.1 What the recursion needs

The recursion uses the multiplicity only through the sum in (0.1). So the relevant quantity is the **averaged** multiplicity `#Dip(X, D)/N(X 2^(-D))`. Three bounds of decreasing strength all feed it the same way:

- **worst case:** `m(j) <= M` for every `j`;
- **average over landing points:** `#Dip <= M |Lambda|` (the form of HYP-9161);
- **average against the counting function:** `#Dip <= M N(X 2^(-D))`.

Each implies the next, because `|Lambda| <= N(X 2^(-D))`. Moreover, the bound may depend on the orbit, on `X` and on the chosen depth. The depth may even be chosen after looking at the orbit, since (0.1) holds for every `D`.

**Hypothesis A(mu)** (for one orbit, `0 <= mu <= 1`). There are `C` and `X_0` such that for all `X >= X_0` and all `1 <= D <= theta_1 L`:

```text
#Dip(X, D) <= C L^mu N(X 2^(-D)) + C L^2.
```

The additive `C L^2` absorbs, for example, landing points with small values (Theorem 2(b)).

### 1.2 Theorem 1

**Theorem 1.** Let `0 <= mu <= 1` and `a*(mu) = mu lambda*/h* - 3/2`.

- **(i) Sufficiency.** If an orbit satisfies A(mu), then for every `a > a*(mu)` there is `K` with `N(X) <= K X^(h*) L^a` for all `X >= 2`. If `mu = 0`, this also holds for `a = -3/2`.
- **(ii) Saturation.** Let `C > 0` and `E = C^(-lambda*/h*)`. The function `psi(X) = X^(h*) L^(a*(mu))` satisfies, for every `X >= 4` and every real `D` in `[1, L/2]`,

  ```text
  psi(X) <= C L^mu psi(X 2^(-D)) + E X^(h*) L^(-3/2) 2^(lambda* D).        (R_mu)
  ```

  So the family of inequalities `(R_mu)` does not imply `N(X) = o(X^(h*) L^(a*(mu)))`. This holds for any choice of depths, even orbit-dependent ones, and for any bootstrap.

**Values.**

| `mu` | `a*(mu)` | reading |
|---|---|---|
| 1 | `-0.986211` | THM-4499 |
| 1/2 | `lambda*/(2h*) - 3/2 = -1.243105` | THM-4499's remark, confirmed |
| 0 | `-3/2` | the ballot floor of THM-4495 |

The whole prize of the multiplicity is therefore the factor `(log X)^(lambda*/h*) = (log X)^(0.5138)`. The parallel S7 note states the same formula and prize; Theorem 1 adds the converse (ii) and the endpoint `mu = 0`.

### 1.3 Proof of (i)

The proof is THM-4499's section 3 with `k` replaced by `C L^mu`; only the exponent bookkeeping changes.

**Case `mu > 0`.**

- Put `eta = (a - a*(mu))/8`, `c_1 = (mu + eta)/h*`, `theta_X = c_1 (log_2 L)/L` (so `D = c_1 log_2 L`) and `s = min(0.1, eta ln 2/c_1)`.
- The no-dip term of (0.1) is at most `4 D_s X^(h*) L^(-3/2 + lambda* c_1 + eta)`. Its exponent is `a*(mu) + eta (lambda*/h* + 1) <= a - 6.48 eta <= a - 2 eta`.
- For the dipper term put `Y = X 2^(-D) = X L^(-c_1)`. By the induction hypothesis,

  ```text
  C L^mu N(Y) <= C L^mu K Y^(h*) (log_2 Y)^a <= C 2^(|a|) K X^(h*) L^(mu - c_1 h* + a) = C 2^(|a|) K X^(h*) L^(a - eta),
  ```

  using `log_2 Y in [L/2, L]`.
- For `X >= X_3` with `C 2^(|a|) L^(-eta) <= 1/4` and `L^(-2 eta) <= 1/4` (and THM-4499's `X_2` conditions), the induction step reads `N(X) <= X^(h*) L^a [(1 + 4 D_s)/4 + K/4] <= K X^(h*) L^a` for `K >= 8 D_s + 1`.
- The base case, the absorption of `k + 2|b| X^(0.585+theta) + C L^2` into `X^(h*) L^(a - 2 eta)`, and the induction on `floor(X)` are as in THM-4499.

**Case `mu = 0`, `a = -3/2`.** Fix a depth `D_0` with `C 2^(3/2) 2^(-D_0 h*) <= 1/4`. For `L >= 2 D_0`:

- `(L - D_0)^(-3/2) <= 2^(3/2) L^(-3/2)`, so `C N(X 2^(-D_0)) <= (K/4) X^(h*) L^(-3/2)`.
- The no-dip term is `4 D_s 2^(lambda* D_0) e^(s D_0) X^(h*) L^(-3/2)`, a constant times `X^(h*) L^(-3/2)`.
- Lemma 1.4c applies because `D_0/L <= theta_1` for large `L`.
- The induction closes for `K >= 2(4 D_s 2^(lambda* D_0) e^(s D_0) + 1)`. ∎

### 1.4 Proof of (ii)

Since `a = a*(mu) < 0` and `D <= L/2`, we have `((L - D)/L)^a >= 1`, so

```text
C L^mu psi(X 2^(-D)) >= C L^mu 2^(-D h*) psi(X).
```

Put `t = 2^(D h*)/(C L^mu)`.

- If `t <= 1`, the first term is at least `psi(X)/t >= psi(X)`.
- If `t > 1`, then `2^(lambda* D) = (t C L^mu)^(lambda*/h*) >= C^(lambda*/h*) L^(mu lambda*/h*)`. The second term is then at least `E C^(lambda*/h*) X^(h*) L^(-3/2 + mu lambda*/h*) = psi(X)`. ∎

**Remarks.**

- THM-4499's true no-dip count has an extra factor `(D + 1)`. Lemma M's control shows `M_k(y) ~ c (y + 1) 2^(hk) k^(-3/2) 2^(lambda* y)`. That factor only enlarges the second term, so (ii) holds a fortiori against the sharp count.
- The control (section A of the `.out`) evaluates both sides for `mu in {0, 1/4, 1/2, 3/4, 1}`, `L in {2^4, ..., 2^20}` and every integer `D in [1, L/2)` (capped at 400), with `K = C = 1`. The smallest value of `log_2(RHS/psi)` is `0.92`.
- **Reading.** Given the no-dip count (sharp as a count of residue classes, THM-4499) and a multiplicity `L^mu`, the exponent `a*(mu)` is exactly what (0.1) yields. **Every improvement of `a*` must lower `mu` or bound the no-dip part of the orbit, not the no-dip set.** Sections 2-3 test the first option.

## 2. The worst case: Theorem 2

Throughout, `j` is a landing index with dippers `i_1 < i_2 < ... < i_m` (`m = m(j) >= 1`) and `D >= 1`.

### 2.1 Lemma S (one dyadic shell)

**Lemma S.** Suppose `b > 0`, or `b < 0` and `y_j >= |b|`. Then:

- the step into `j` is a halving: `y_(j-1) = 2 y_j`;
- no dipper of `j` is `j - 1`;
- every dipper satisfies `2^D y_j < y_(i_r) <= 2^(D+1) y_j`.

*Proof.*

1. Let `i` be a dipper of `j`, so `y_j < 2^(-D) y_i`.
   - If `i <= j - 2`, then `j - 1` lies strictly between `i` and `j`, so it is not a dip time for `i`: `y_(j-1) >= 2^(-D) y_i > y_j`.
   - If `i = j - 1`, then `y_j < 2^(-D) y_(j-1) < y_(j-1)`.
   - Either way `y_(j-1) > y_j`.
2. An odd step satisfies `(3y + b)/2 > y` iff `y > -b`. This always holds for `b > 0`. For `b < 0` it holds because `y_(j-1) > y_j >= |b|`. So the step `j - 1 -> j` is even: `y_(j-1) = 2 y_j`.
3. If `i = j - 1` were a dipper, then `y_j < 2^(-D) 2 y_j`, which forces `D < 1`. So `i <= j - 2`.
4. The non-dip condition at `j - 1` then gives `2^(-D) y_i <= y_(j-1) = 2 y_j`. The dip condition gives the lower edge. ∎

### 2.2 Lemma O (one dipper per odd letter)

**Lemma O.** Under the hypotheses of Lemma S, every stretch `[i_r, i_(r+1))` between consecutive dippers contains an odd step.

*Proof.* Otherwise all steps in the stretch halve, so `y_(i_(r+1)) = 2^(-(i_(r+1) - i_r)) y_(i_r) <= y_(i_r)/2 <= 2^D y_j`. This contradicts the lower edge of Lemma S for `i_(r+1)`. ∎

### 2.3 Theorem 2 and its proof

**Theorem 2.**

- **(a)** If `b > 0`, then `m(j) <= ceil((k - D)/alpha)` for every landing index `j` and every `D >= 1`.
- **(b)** If `b < 0`, then `m(j) <= ceil((k - D + log_2(3/2))/alpha)` for every landing index with `y_j >= k|b|`. The landing indices with `y_j < k|b|` have distinct values, so there are fewer than `k|b|` of them. They carry fewer than `k^2 |b|` dippers in total, an additive term absorbed by Hypothesis A's `C L^2`.

*Proof.*

1. Let `n = j - i_1 <= k`, and let `O` be the number of odd steps among the positions `i_1, ..., j - 1`. By Lemma O, `O >= m - 1`.
2. Composing the steps, where an odd step is `y -> (3/2) y (1 + b/(3y))`:

   ```text
   y_j = y_(i_1) (3^O / 2^n) prod_(odd t in [i_1, j)) (1 + b/(3 y_t)).
   ```

3. **(a)** For `b > 0` the product is at least 1. With the dip condition `y_(i_1) > 2^D y_j` this gives `y_j > 2^D y_j 3^O/2^n`, i.e. `O alpha < n - D <= k - D`. Hence `m - 1 <= O < (k - D)/alpha`, and since `m` is an integer, `m <= ceil((k - D)/alpha)`.
4. **(b)** Every `t in [i_1, j)` has `y_t > y_j`. The dippers exceed `2^D y_j`, and every other index in the stretch is at least `2^(-D) y_(i_1) > y_j`.
   - So each factor is at least `1 - |b|/(3 y_j) >= 1 - 1/(3k)`.
   - There are `O <= n <= k` factors, so by Weierstrass's inequality the product is at least `2/3`.
   - Then `3^O < (3/2) 2^(n - D)`, i.e. `O alpha < k - D + log_2(3/2)`, and the bound follows as in (a). ∎

**Remarks.**

- The deterministic time constraint (dip time `s > D`) is implicit: `O >= 0` forces `n > D`. It gives nothing beyond (a).
- **Two places.** The proof uses exactly two pieces of information: the real place (the size ratio `y_j/y_(i_1)`) and the 2-adic place (the letters, Lemma O). Sharpness (2.4) shows that no further two-place information improves the worst case. Every parity word is a residue class (Terras), and for large `y` the carries are negligible.
- **Cost for THM-4499.** Replacing `k` by `ceil((k - D)/alpha) <= 0.631 k + 1` in the recursion lowers the constant `K` only: `mu = 1` in Theorem 1.
- **Injective invariant sets** (THM-4476, Corollary 4). The same bound holds verbatim. `T_b` is injective on the set, so the dippers of a landing point lie on its unique backward chain inside the set, which is a segment of distinct integers.
- **Parallel with THM-4478.** THM-4478 bounds the integer ancestors of one endpoint per (time, odd count) by `floor(k/3) + 1`. That bound concerns the inverse tree. Inside one orbit there is one ancestor per time, and Lemma O plus the ratio bound play the corresponding role: one dipper per odd count, over at most `(k - D)/alpha` odd counts. THM-4478 is not needed for a single orbit, or for an injective invariant set, where the same holds.

### 2.4 Sharpness

**Exhaustive (FINITE-EXACT, section B of the `.out`).**

- Every landing configuration whose dippers are `<= X` has its first dipper in `[1, X]`.
- The C scan therefore takes each `y <= X = 2^(k+1) - 1` as a first dipper. It computes `y`'s landing index and counts, with exact 128-bit comparisons, the dippers of that index among `y`'s successors. It checks Lemmas S and O, the halving into `j`, and Theorem 2 in every configuration.
- The worst case over all integers equals `ceil((k - D)/alpha)` in all 144 cells for `b = 1` (`k = 10, ..., 25`; `D in {1, 1.5, ..., 4, 5, 6}`), in all 144 cells for `b = -1`, and in all 78 cells for `b = 5`.
- For `b = -5` (`k <= 22`) all 78 cells are within bound (b):
  - 74 cells are equal to `ceil((k - D)/alpha)`;
  - one cell is 1 above it (allowed by the `log_2(3/2)` term);
  - three cells are below it (`k = 10, 11, 12` with `D = 6`), because of the `y_j >= k|b|` restriction.
- There are zero violations. Excerpt (`b = 1`, maximum / bound):

```text
   k |  D=1.0   D=2.0   D=3.0   D=4.0   D=5.0   D=6.0
  15 |  9/9     9/9     8/8     7/7     7/7     6/6
  20 | 12/12   12/12   11/11   11/11   10/10    9/9
  25 | 16/16   15/15   14/14   14/14   13/13   12/12
```

**Construction (PROVED lower bound + FINITE-EXACT, section C).** Work with the pure walk `S_t = o_t alpha - t`, where `o_t` is the number of odd letters among the first `t`.

- A Sturmian word visits the unit band `[-a, 1-a)` at every odd count `o = 0, ..., O`, at the times `n(o) = floor(o alpha + a)`.
- Between visits it uses the letters `1` or `10`. It is followed by `floor(D) + 1` halvings.
- The offset is chosen so that `frac(O alpha + a) = 1 - frac(D) - eps`. Then the landing shell of Lemma S is the visit band shifted by `eps`, and in the pure walk every visit is a dipper.
- Terras realises the word by a residue class. Taking its largest member below `X 2^(a-1)/1.01` gives actual integers.
- For `b = +-1`, `k in {30, 40, 64, 100, 200, 400}` and `D in {1.05 log_2 k, 0.1k, 0.2k}`, the actual multiplicity is `ceil((k-D)/alpha)` in 18 of the 36 cases and one less in the other 18. At `k = 400`, `D = 9.08` it is 247 of 247 (`m/k = 0.618`).
- The rigorous asymptotic statement with a uniform margin is Proposition H below. It is weaker by a factor 2 because it keeps only visits with margin `1/8`. Keeping all visits needs an irrationality measure of `log_2 3`; by Baker's theorem the loss is then `O(1)` (CITED, not written out). So **the worst case is `(k - D)/alpha + O(1)`**. The upper bound is proved in full. The matching lower bound `(k - D)/alpha - O(1)` rests on Baker (CITED). Equality holds wherever it was computed.

### 2.5 Climbs and W-bit hovers (Proposition C)

**Proposition C.** Let `b > 0`.

- **(a) Odd runs.** At most 2 indices of an odd run (consecutive odd steps) are dippers of the same landing index.
- **(b) W-bit hovers.** Take a stretch `[p, q)` whose values stay within `W` bits of each other (`max/min <= 2^W`), followed by halvings down to below `2^(-D-1)` times its minimum. Suppose every dipper among the stretch's indices lands in the crash. Then those dippers use at most `ceil(W) + 1` landing indices. Hence their average multiplicity is at least `(number of dippers in [p, q))/(ceil(W) + 1)`.

*Proof.*

- (a) An odd step multiplies by more than `3/2`, and `(3/2)^2 > 2`. So a half-open dyadic shell contains at most 2 points of an odd run. Apply Lemma S.
- (b) By Lemma S, a landing index `j` in the crash collects only dippers in `(2^D y_j, 2^(D+1) y_j]`. The crash values are `y_j = 2^(-r) y_q` for consecutive `r`, so these shells are consecutive dyadic shells. A range of width `W` bits meets at most `ceil(W) + 1` of them. ∎

**Controls (section G).**

- `1^40 0^35`, realised by `y = 72464114801907393888255` at `X = 2^100` with `D = 6.97`: the 40 climb dippers spread over 24 landing points, at most 2 each.
- A random 4-bit hover of 60 steps followed by 13 halvings: 60 dippers on 4 landing points, with multiplicities `7, 10, 15, 28`.
- The probe segments of the parallel note: the maximal landing points take 0 or 1 dipper from the initial odd run (section 5.4).

## 3. Averaged multiplicities

By section 1.1, any averaged bound feeds the recursion exactly as a worst-case bound does, even an orbit-dependent bound with an orbit-dependent depth. The question is whether some averaged bound with `mu < 1` can be *proved*. This section tests the three natural averages:

- over depths (Lemma A);
- over integers, by counting residue classes (Proposition H);
- over the parity word, i.e. the 2-adic data (section 3.3).

### 3.1 Lemma A (depth averaging) and why it is neutral

**Lemma A.** Let `b > 0` (or `b < 0` with `y_j >= |b|`), let `j` be an index and `i < j`.

- The set of real depths `D >= 1` for which `i` is a dipper landing at `j` lies in a half-open interval of length at most 1. So it contains at most one integer.
- Consequently `sum_(D = 1, 2, ...) m_D(j) <= k - 1`.
- For every `D_0` and `W >= 1`, some integer `D in [D_0, D_0 + W)` has `#Dip(X, D) <= (k - 1) N(X 2^(-D_0))/W`.

*Proof.*

1. The condition is `y_j < 2^(-D) y_i <= mu_ij`, where `mu_ij = min_(i<t<j) y_t` (and `mu_ij = +infinity` if `i = j - 1`). Equivalently `D in [log_2(y_i/mu_ij), log_2(y_i/y_j))`.
2. For `i <= j - 2` the length is `log_2(mu_ij/y_j) <= log_2(y_(j-1)/y_j)`.
   - If `y_(j-1) < y_j` the interval is empty.
   - Otherwise the step into `j` halves (Lemma S, step 2), and the length is at most 1.
3. For `i = j - 1` the set is `{D : D < log_2(y_(j-1)/y_j)}`, which contains no integer `D >= 1`.
4. Hence each of the at most `k - 1` indices `i in [j - k, j - 2]` serves at most one integer depth.
5. Every landing index at depth `D >= D_0` lies below `X 2^(-D_0)`. Summing over them gives `sum_(D_0 <= D < D_0 + W) #Dip(X, D) <= (k - 1) N(X 2^(-D_0))`. A minimum is at most the average. ∎

**Neutrality.** Take the best of the `W` depths in (0.1) and write `D' = D_0 + W - 1` for the deepest one, where the no-dip set must be paid. In the scale of `D'` the dipper term is

```text
(k/W) N(X 2^(-D_0)) = (k 2^((W-1)h*)/W) * [2^(-(W-1)h*) N(X 2^(-D_0))],
```

and `psi` satisfies the corresponding inequality with the factor `k 2^((W-1)h*)/W` in place of `C L^mu`.

- The factor `2^((W-1)h*)/W` is smallest at `W = 2`, where it equals `0.965907` (`.out`, section A).
- Replacing `psi(X 2^(-D_0))` by `2^((W-1)h*) psi(X 2^(-D'))` costs a factor `((L - D_0)/(L - D'))^(a*) = 1 - o(1)`. So Theorem 1(ii) applies at depth `D'` with `mu = 1` and `C = 0.965907 (1 - o(1)) k/L`, and it holds for every `C > 0`. Depth averaging therefore changes the multiplicity by at most the factor `0.966`, and `mu = 1` stays.

**Weighted depths.** A weighted choice over depths is no better. The dipper side can only use `max_D w_D` per pair `(i, j)`, and the no-dip side pays `sum_D w_D 2^(lambda* D)`. By the same convexity the optimum sits at a single depth.

**Control (section E).** 18929 windows (`k = 40`) of the record orbits, 757160 pairs `(i, j)`: every pair serves at most one integer depth.

### 3.2 Proposition H (orbit-blind splits need a linear threshold)

A natural averaged argument splits the dippers at a threshold `M`:

- a landing index with `m(j) <= M` costs at most `M`;
- in a heavier landing index, every dipper except the last `M` is *heavy*.

A heavy dipper `z` has at least `M` further dippers at its landing index, so it lies in

```text
H_(M+1)(X) = { z <= X : z is a dipper of its own forward segment, and its landing index receives >= M+1 dippers from that segment }.
```

Therefore `#Dip <= M |Lambda| + #(orbit ∩ H_(M+1)(X))`. An **orbit-blind** version bounds the last term by `#Z` for some set `Z ⊇ H_(M+1)(X)` that does not depend on the orbit. Examples:

- the integers whose `k`-word has at least `M+1` returns to a unit band before a crash;
- a strip set counted by THM-4495 or by the Robin bounds.

**Proposition H.** Let `b > 0`, `D >= 3` (dyadic) and `O >= 1`. Build the word:

- put `eps = 10^(-6)`, `a = frac(1 - frac(D) - eps - O alpha)`, and `n(o) = floor(o alpha + a)` for `0 <= o <= O`;
- `w` is the Sturmian hover with visits at the times `n(o)` (letters `1` or `10`), followed by `c = floor(D) + 1` halvings;
- its length is `l = n(O) + c <= O alpha + D + 2`.

Call a visit `o` *safe* if `o = O` or `||(O - o) alpha + D|| >= 1/8`, where `||x||` is the distance from `x` to the nearest integer. Then:

- **(i)** there are at least `floor(O/2) + 1` safe visits;
- **(ii)** for every `X` with `l <= floor(log_2 X)` and every `y ≡ r_w (mod 2^l)` with `16 b l < y <= X/4`, every safe visit of the segment of `y` is a dipper (scale `X`, depth `D`) landing at index `l`;
- **(iii)** hence, taking `O = 2M - 2`, `#H_M(X) >= X 2^(-(2 alpha (M - 1) + D + 5))` for every `M >= 2` and every `X >= b l 2^(l + 8)`.

*Proof.*

- **(i)** If `r` and `r + 1` are both unsafe, then `||alpha|| <= ||r alpha + D|| + ||(r+1) alpha + D|| < 1/4`. But `||alpha|| = 2 - alpha = 0.415`. So among `r = O - o in [1, O]` no two consecutive values are unsafe, which gives at least `floor(O/2)` safe values, plus `r = 0`.
- **(ii), pure walk.** Visits sit at `V_o = frac(o alpha + a) - a in [-a, 1 - a)`. The landing index `j = l` sits at `S_j = V_O - c`, and `frac(O alpha + a) = 1 - frac(D) - eps` by the choice of `a`. So

  ```text
  V_o - D - S_j = frac(o alpha + a) + eps ≡ -(r alpha + D)  (mod 1),     r = O - o.
  ```

  Lemma S's conditions for `o` (dip: `> 0`; non-dip at `j - 1`: `<= 1`) hold with margin `||r alpha + D|| >= 1/8` for safe visits. For `r = 0` they hold exactly, since `y_(n(O)) = 2^c y_j` and `D < c <= D + 1`.
- **(ii), carries.** Write `T^t(y) = 2^(S_t) (y + h_t)` with `0 <= h_t = sum_(odd p < t) b/(3 * 2^(S_p))` (THM-4478, section 2). On the hover `S_p > -1`, so `h_t < 2bl/3` and `0 <= log_2(1 + h_t/y) < bl/y < 1/16`. This shift is nondecreasing in `t`, so it moves every difference above by less than `1/16`, which is below the margin.
- **(ii), first dip.**
  - Hover values after a visit lie at most one bit below it, and `D >= 3`.
  - Crash values before `j` are at least `S_j + 1`, which sits at least the margin above the threshold.
  - All values are at most `4y <= X`, and `l <= k`.
- **(iii)** Let `o_0` be the first safe visit. The map `y -> T^(n(o_0))(y)` is injective on the class, and every image lies in `H_M(X)` once `O = 2M - 2` (so there are at least `M` safe visits). The class meets `(16 b l, X/4]` in at least `X 2^(-l-2) - 16 b l 2^(-l) - 1 >= X 2^(-l-3)` integers when `X >= b l 2^(l+8)`, and `l + 3 <= 2 alpha (M-1) + D + 5`. ∎

**Consequence.** Suppose an orbit-blind split closes the bootstrap at exponent `a`. Then it needs `#Z <= K' X^(h*) L^a`, while `#Z >= #H_(M+1)(X) >= X 2^(-(2 alpha M + D + 5))`. So

```text
M >= ((1 - h*) L - D - a log_2 L - 5 - log_2 K')/(2 alpha) = Theta(L),
```

and the light part costs `C L^1`: `mu = 1` again.

- The rigorous slope is `(1-h*)/(2 alpha) = 0.0158`.
- The true heavy counts decay by only `0.7`-`1.0` bits per unit `M` (section B, `k = 25`). The true threshold is therefore about `(1-h*) L/0.7 ~ 0.07 L`.
- The Robin strip bounds and the ballot counts of THM-4495 can only *upper*-bound `#H_M`. They sharpen the slope, not the linearity.

**Controls (section D).** Five classes (`O = 10, ..., 120`; `D = 3, ..., 9.25`; `l = 19, ..., 200`), 40 random members each. In every member every visit was a dipper (`m = O + 1`), not only the safe ones.

### 3.3 Two places: the parity word alone permits hostility at every scale; one integer fixes one window

**2-adic data.** Terras makes every parity word the word of a 2-adic integer. The word of section H consists of 20 cycles, each a hover (aligned as in 2.4), then a crash of `floor(D) + 1 = 5` halvings, then a climb back to the starting level. In the pure-walk model, with positions `o_t alpha - t` at the design scale:

| design `k` | multiplicity at every cycle's landing point | mean over all landing points | `#Dip/N(X 2^(-D))` |
|---|---|---|---|
| 40 | 23 = `ceil((40 - 4)/alpha)` | `0.42 k` | 6.6 |
| 80 | 48 = `ceil((80 - 4)/alpha)` | `0.45 k` | 14.6 |

The averaged ratio grows linearly with the design window. So **no argument that uses only the parity word can bound the averaged multiplicity.** This includes strip and ballot counts, the Sturmian structure, and anything else 2-adic.

**The real place.** An integer `y <= X` prescribes only `floor(log_2 X)` letters. The integer realising the 967-letter word has 967 bits. At its own scale (`k = 969`, `D = 10.42`) the 5-bit crashes are not dips, and the word contains **no landing point at all**. Hostility at every cycle is a 2-adic phenomenon. For an actual orbit, whether hostile windows recur after the first one is decided by the orbit's arithmetic, not by the word.

The method's only real-place input is the count of integers `<= X` per residue class. That is the setting of Proposition H, and it permits `Theta(L)`.

The two-place carrier of the crossroads note (section 3: normalised centres in `[M, M + O(M^gamma)]` along a tail above `M`) says the carries along a divergent tail are relatively `O(M^(gamma-1))`. So a divergent orbit follows the pure-walk model faithfully. It does not constrain the word.

### 3.4 What HYP-9161 needs, in the language of Lemma S

By Lemmas S and O, `m(j)` is a **local time**: the number of visits, at distinct odd counts, that the orbit pays to the one dyadic shell `(2^D y_j, 2^(D+1) y_j]` during the `k` steps before the final crash into `j`. HYP-9161 asks that these local times be `O(L^beta)` on average over landing points. Section 5.3 shows that on actual record orbits they are:

- 2.5-4.3 on average;
  - at the bootstrap's depth `theta_X` that is `0.55-0.7 D`;
  - at fixed `theta` the mean grows much more slowly than `D` (`0.24-0.67 D` at `theta = 0.2`);
- up to `0.6 L` at the worst landing point (`L = 20`), and `0.40-0.53 L` at `theta_X`.

By Proposition H and section 3.3 a proof must use how the orbit continues beyond one window. This is an oscillation statement, as the parallel S7 note also concludes.

## 4. Answer to the brief: every approach, its type, and what it gives

The recursion accepts every row's bound, worst-case or averaged, orbit-dependent or not (section 1.1).

| approach | worst case or averaged | best bound on the multiplicity | status | effect on `a*` |
|---|---|---|---|---|
| THM-4476 pigeonhole (one dipper per dip time) | worst case | `k` | PROVED | `a*(1) = -0.9862` |
| deterministic time (a dip needs `s > D`) | worst case | `k - D` | PROVED | none |
| shell + odd separation, real and 2-adic places together (Theorem 2) | worst case | `ceil((k - D)/alpha) ~ 0.631 (k - D)` | PROVED, sharp (exhaustive to `k = 25`, constructed to `k = 400`, records reach it) | none (constant only) |
| hover-then-crash counted by residue classes, strip and ballot sets (Proposition H) | averaged over integers, orbit-blind | the threshold must be `Theta(L)` | PROVED (obstruction) | none |
| 2-adic information alone (section 3.3) | averaged over the word | no bound: `#Dip/N` grows linearly | PROVED (construction) | none |
| depth averaging (Lemma A) | averaged over depths | effective factor `>= 0.966 k` | PROVED | none |
| average over the landing points of one orbit (HYP-9161) | averaged, orbit-coupled | `O(L^beta)` conjectured; observed means 2.5-4.3 | OPEN | would give `a*(beta)`; `-3/2` if `beta = 0` |

**Improved thin-divergence theorem: none.**

- THM-4499 stands as stated. The only change is cosmetic: its recursion (R'') holds with `ceil((k - D)/alpha)` in place of `k`, which changes `K`.
- By Theorem 1(ii), `a* = lambda*/h* - 3/2` is the exact output of the one-window method with its present inputs.
- Going below `a*` requires one of two things:
  - an orbit-coupled bound on the averaged local times (HYP-9161);
  - a bound on the no-dip *part of the orbit* that beats the no-dip *set*, which is THM-4487's "one free window per element" wall.

## 5. Controls (FINITE-EXACT; all from the `.out`)

### 5.1 Exhaustive worst case

This is section B of the `.out`, summarised in 2.4.

- **Cells.** `b = 1, -1`: `k = 10..25`, `D in {1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6}`. `b = 5, -5`: `k = 10..22`, `D in {1, 1.5, 2, 3, 4, 6}`.
- **Matches with `ceil((k - D)/alpha)`.** 144/144, 144/144 and 78/78 cells for `b = 1, -1, 5`. For `b = -5`, 74/78 (all 78 within the `b < 0` bound).
- **Heavy counts** (`b = 1`, `k = 25`, `D = 2`, `X = 2^26 - 1`). The share `log_2(H_M/X)` for `M = 1, ..., 15` is

  ```text
  -0.17 -0.95 -1.66 -2.35 -3.03 -3.73 -4.48 -5.31 -6.21 -7.22 -8.38 -9.79 -11.58 -13.90 -17.21
  ```

  This is about 0.7 bits per unit `M` until the window boundary. Exactly 443 integers (`2^(-17.21) X`) attain the maximum 15.

### 5.2 Constructed hostile integers (section C)

| `b` | `k` | `D` | `m(j)` | `ceil((k-D)/alpha)` | `m/k` |
|---|---|---|---|---|---|
| 1 | 30 | 5.156 | 16 | 16 | 0.533 |
| 1 | 64 | 12.797 | 33 | 33 | 0.516 |
| 1 | 100 | 6.969 | 59 | 59 | 0.590 |
| 1 | 200 | 8.031 | 121 | 122 | 0.605 |
| 1 | 400 | 9.078 | 247 | 247 | 0.618 |
| 1 | 400 | 80.000 | 202 | 202 | 0.505 |
| -1 | 400 | 40.000 | 227 | 228 | 0.568 |

Across all 36 cases (`b = +-1`, `k in {30, 40, 64, 100, 200, 400}`, three depths each), the gap to the bound is 0 in 18 cases and 1 in 18. Each segment was checked to consist of distinct integers, with all dippers `<= X` and Lemmas S and O holding.

### 5.3 Records up to `10^12` (section F)

This covers 144 distinct starts (90 delay records, 61 path records). `D = theta_X L` uses `theta_X = 1.05 log_2 L/L`, the bootstrap's choice. The table lists landing points, dippers, the mean multiplicity, the maximum (and its record start), the maximum's ratio to `ceil((k-D)/alpha)`, and `#Dip/N(X 2^(-D))`:

```text
   L  depth       D   landing  dippers  mean-mult  max-mult (record)       max/bound   #D/N(X2^-D)
  20  theta_X  4.531     1929     4884     2.532      9 (  59152641055)      0.900         0.434
  20  0.1L     2.000     3108     8415     2.708     12 (        13255)      1.000         0.603
  30  theta_X  5.156     3178     9106     2.865     15 (        13255)      0.938         0.440
  40  theta_X  5.594     3973    13582     3.419     21 (   4578853915)      0.955         0.410
  40  0.2L     8.000     2776     7806     2.812     16 (   4578853915)      0.762         0.263
  50  theta_X  5.922     4353    17061     3.919     22 (     80049391)      0.786         0.385
  60  theta_X  6.203     4454    18835     4.229     24 (      3542887)      0.706         0.394
  60  0.2L    12.000     2123     6170     2.906     15 (    210964383)      0.484         0.132
```

All 15 rows are in the `.out`. The pattern:

- The maximum is a large fraction of the worst case (`0.48-1.0`), so a worst-case bound `o(L)` is false on actual Collatz orbits.
- The mean over landing points is 2.5-4.3, growing slowly with `L` and roughly proportional to `D`.
- The recursion's actual averaged ratio `#Dip/N(X 2^(-D))` is below 1 throughout.

This is the behaviour HYP-9161 predicts. These are record orbits that reach 1 (eventually periodic), so they illustrate but do not test the hypothesis for divergent orbits.

### 5.4 The probe segments of the parallel note (section G)

`D = theta_X L`. The last column is the largest number of dippers that any landing point takes from the initial odd run.

```text
   L  start        D      max-mult  dippers(index range)  initial odd run  max from the run
  20  2^L-1       4.531      4       64..70                   20              0
  40  2^L-1       5.594     11      149..177                  40              0
  60  2^(L-1)+1   6.203     21       54..105                   1              1
  80  2^L-1       6.641     19      416..461                  80              0
  80  2^(L-1)+1   6.641     11      263..292                   1              1
```

The maximal landing points collect dippers from stretches well after the start. On `2^L - 1` the initial run climbs above `X` at once, and `2^(L-1) + 1` has no run. So the "climb" attribution in the parallel note does not hold (Proposition C(a): at most 2 in any case).

### 5.5 Proposition H classes and the 2-adic word (sections D and H)

**Proposition H classes.** Five classes: `(O, D, l)` = `(10, 3, 19)`, `(20, 4, 36)`, `(40, 5.5, 69)`, `(80, 7, 134)`, `(120, 9.25, 200)`. Forty random members each, at `X = 2^(k+1) - 1` with `k = 40, ..., 260`.

- In every member all `O + 1` visits are dippers of index `l`.
- The safe visits alone number 9, 16, 30, 62 and 91.

**The 2-adic word** is tabulated in 3.3. Its 967-bit integer realisation follows the whole word (checked). At the integer's own scale it has no landing point inside the word.

## 6. Failures, caveats and what was not done

**Where a first attempt was wrong.**

- The first 2-adic control used an unaligned hover offset, and its designed landing points received only 1-6 dippers. The aligned construction (visit band = landing shell) fixed this. The last cycle of a finite word has a truncated window and is excluded (one extra cycle is built).
- An early library version compared depths as fractions with huge denominators and was unusably slow. Depths are now dyadic with denominator at most 64, compared exactly after a float fast path with a `1e-7` guard.

**Scope of the worst case.** The worst case is taken over **segments of distinct positive integers** (all integer starts), not over non-eventually-periodic orbits, which are unknown. The hostile integers and the record orbits reach 1. They show that no bound derived from windows can beat `ceil((k-D)/alpha)`. They do not show that a divergent orbit contains such windows.

**Sharp lower bound.**

- The asymptotic lower bound proved in full is Proposition H's `floor(O/2) + 1` visits, with a margin, at every large member of a class.
- The `(k - D)/alpha - O(1)` bound for every large `k` needs a polynomial lower bound for `||q log_2 3||`. That follows from Baker's theorem on linear forms in two logarithms (CITED; not written out).
- The exact equality is FINITE-EXACT: exhaustive for `k <= 25`, constructed within 1 for `k <= 400`.

**The case `b < 0`.** The bound for `b < 0` carries `+ log_2(3/2)` and needs `y_j >= k|b|`. In one `b = -5` cell the worst case is 1 above `ceil((k-D)/alpha)`, which the `b < 0` bound allows. The recursion needs both signs, because negative segments are conjugated to `T_(-b)`. The additive term for small landing values is `O(k^2 |b|)`.

**Theorem 1(ii) is about the recursion, not about orbits.** It shows the inequalities cannot give more. It says nothing about the true size of `N(X)` for a divergent orbit, which is expected to be `O(log X)`.

**The 2-adic model has no carries and no injectivity.** It models word-level information only, which is the point of section 3.3.

**Records.** The record lists are CITED from OEIS b-files (A006877, A006884). Record status was not re-verified; the orbits themselves are exact.

**Not attempted.**

- An orbit-coupled proof of HYP-9161.
- The Lean formalisation.
- Changing THM-4499's file. The cosmetic replacement of `k` is left to its owner.

## 7. Reproduction

```text
python3 -u 04-computation/experiments/procgen_landing_20260926_run.py
```

The runner compiles `procgen_landing_20260926_scan.c` with `cc -O2` into `scratch/procgen_landing/`. It runs the four scans at most two at a time and writes `05-knowledge/results/procgen_landing_20260926.out`.

The OEIS b-files are taken from `scratch/procgen_landing/oeis/`. If absent, they are copied from the family27 lane's scratch cache or fetched with `curl` using the generic User-Agent `Mozilla/5.0 (research; math-repo)`; no personal data is sent.

- **b-file hashes.** `b006877.txt` sha256 `4bdd9c264ffb33f6252adad6716f8f0ad55c84687d4300a9bfbe646a1a81327f`; `b006884.txt` sha256 `657ec2fbc0ade279ce052005ae11234ea795c47b0fdaf9add6fff5d4f214de86`.
- **Run on this machine:** Apple M-series, Python 3.10.0, Apple clang 17.
- **Hashes** (raw LF bytes):


```text
a4a5b54094a27678ba6531b787e5035b7a66c926779c07ba7e17267ddd75b491  04-computation/experiments/procgen_landing_20260926_lib.py
4d41acf48658a531e52a797d788aa28e5ac0ca0b4790ef316f9c100af864d458  04-computation/experiments/procgen_landing_20260926_run.py
59ba0b059794d5300b029e62bd5bdbe6f45833af9cf726981020e8f1e0b8042d  04-computation/experiments/procgen_landing_20260926_scan.c
31bb595b698ab29bcb4f504d139c7ef48bef9b09ad89bf7043c5adaea7962723  05-knowledge/results/procgen_landing_20260926.out
```

- **Wall time and memory.** 184.4 s wall time. Peak RSS: runner 14.4 MB, largest child 47.1 MB (the C compiler; each scan process is under 1 MB).
- **Checks.** 130 `ok:` lines, ending with `ALL CHECKS PASSED`.
- **Temporary files.** Scan tables, the compiled binary and the b-file cache live in `scratch/procgen_landing/` and are not committed.
