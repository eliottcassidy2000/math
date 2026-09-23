# The hard class: the Diophantine exponent decides Theorem S, a q-series method goes past it, and C2 is a coupled Z-number problem

**Status.**
* **PROVED** (this note; elementary given the cited inputs; independent exact checks in the `.out`):
  * **Theorem D.** `Phi_T(w)` is irrational whenever `Dio(w) > eta(w)`. Here `Dio` is the Adamczewski–Bugeaud Diophantine exponent and `eta` the height exponent (`= max(1, beta log2 m)` for `x/2, (mx+r)/2`).
  * **Corollary D1.** Every eventually Sturmian word, and every word `u S(s)` with `s` Sturmian and `S(01) != S(10)` (this includes every recurrent quasi-Sturmian word), under every affine 2-adic shift map with `eta < 5/3 + 4 sqrt10/15 = 2.50994`. This covers every slope for `3x+r` and for Mahler's map, and **every slope for `5x+1`** (the parent note had `alpha < 0.804`). For `7x+1` it covers `beta < 0.894`. The threshold is sharp for the method.
  * **Corollary D2.** Characteristic Sturmian words are covered for `eta < 1 + limsup q_(n+1)/q_n`, which is at least `2.618`.
  * **Lemma J'.** Codings of rotations satisfy `Dio >= mu*(k0, L)`.
  * **Lemma F.** For fixed points of primitive substitutions, `Dio` is the supremum of weighted prefix-pattern ratios, so every certificate is finite.
  * **Proposition Y.** The square-swap word `Y` is explicit, supercritical (`beta = 9/10`), of zero entropy and of discrepancy width 1, and has `Dio(Y) = 1`: it lies beyond Theorem D.
  * **Lemma P and Theorem Y.** A 2-adic Tschakaloff–Padé argument shows that no rational has an eventually-`Y` parity vector. The same holds for every square-swap block word under `3x+r`.
  * **Proposition M.** C2 for integers is equivalent to a coupled Z-number statement. The decoupled statement MZ implies C2.
  * **Propositions T and Rnd** and two calibrations.
* **FINITE-EXACT:**
  * The Dio engine reproduces the Bugeaud–Kim equality case (`rep = 1.66229`, `Dio = 2.50991`) on `2^20` letters.
  * 31 isometry-checked families of Theorem-R certificates (every one passes).
  * Censuses: 9,716 random supercritical primitive substitutions, each with a finite certificate under `3x+1`, and 160 golden multi-arc rotation codings.
  * Twelve C2 strip words (random, greedy, Sturmian-perturbed, `Y`, `Y3`): best Theorem-R gains are 5–32 bits for the positive-entropy ones and 78 bits for `Y`. For `Y3` the gains are positive only below `|UV| ~ 2.7*10^4` (best 1533 bits) and negative beyond. For each tested word, no rational of height `<= 2^4999` has that parity vector (for each C2 strip word it was checked by reconstruction at `N = 10^4` bits).
  * The `Y` identity holds modulo `2^30000` and the Tschakaloff linear forms are exact for `n <= 14`.
* **CITED:**
  * Bugeaud–Kim 2019; Adamczewski–Bugeaud 2007 (ETDS and Annals) and 2011; Berthé–Holton–Zamboni 2006.
  * Cassaigne 1998, via Damanik–Lenz; Durand and Durand–Host–Skau; Adamczewski 2003.
  * Zudilin 2005, Tschakaloff and Bundschuh; Väänänen–Wallisser 1991 (metadata only); Ghidelli 2019; Dubickas 2006; FLP 1995.
* **OPEN:**
  * Collatz (T1), the Periodicity Conjecture (PC), and C2 at every supercritical slope and every width `>= 1`.
  * The cube-swap word `Y3` (HC1). This is the smallest open instance found.
  * Whether some linearly recurrent word under `3x+r` has `Dio <= eta` (HC2).
  * HC3–HC6.
* **REFUTED:**
  * The description of HARD as "positive entropy". `Y` and `Y3` are zero-entropy words that no proved periodic-approximant mechanism reaches.
  * The extrapolation "periodic approximants settle every Sturmian word under every affine map". Bugeaud–Kim's word `s'` under `19x+1` has `Dio = 2.50994 < eta = 2.602`, so Theorem D does not apply, and every periodic approximant tested has negative gain.
  * "The Dio barrier bounds every approximation argument". Theorem Y beats it for `Y` with Padé approximants of anomalously small height.
  * "Word-free Mahler statements can imply C2 for wide strips". MZ is vacuous whenever every `K_q` has length `>= 1`.

Session `collatz-procgen-20260922`, hard-class lane, 2026-09-23. It builds on the [transversality foundry](collatz_procgen_20260922_transversality_foundry.md) (Theorem R, Lemma H, Theorem S, Proposition B, C1–C6), the [in-house discrepancy theorem](collatz_guards_20260921_discrepancy.md) and [foundry v4](collatz_procgen_20260922_foundry.md), section 3b. Collatz, PC and C2 stay **OPEN**. No priority is claimed: experts may know Theorem D, Lemma F and Lemma P in other guises.

* **Scripts:**
  * [`..._hard.py`](../../04-computation/experiments/collatz_procgen_20260922_hard.py) (sections H0–H6, C1–C4);
  * [`..._hard_words.py`](../../04-computation/experiments/collatz_procgen_20260922_hard_words.py) (exact generators, quadratic irrationals via `isqrt`);
  * [`..._hard_verify.py`](../../04-computation/experiments/collatz_procgen_20260922_hard_verify.py) (independent exact verifier);
  * [`..._hard_lpf.c`](../../04-computation/experiments/collatz_procgen_20260922_hard_lpf.c) (suffix array plus Longest-Previous-Factor array);
  * [`..._hard_run.sh`](../../04-computation/experiments/collatz_procgen_20260922_hard_run.sh).
* **Output:** [`collatz_procgen_20260922_hard.out`](collatz_procgen_20260922_hard.out).

## 0. Answers in brief

1. **C3, the provable side.** Theorem R's reach is governed by a single number, the Diophantine exponent `Dio(w)` of Adamczewski–Bugeaud.
   * **Theorem D.** `Dio(w) > eta(w)` implies `Phi_T(w)` is not rational.
   * In the requested shape: infinitely many prefixes `U V^e` with `e >= mu + (mu-1)|U|/|V| + delta(1 + |U|/|V|)` and `|UV| -> infinity`. Theorem S's options A and B are the cases `U = empty` and `U = s[0, j0]` (section 1.4).
   * Bugeaud–Kim's theorem `rep <= sqrt10 - 3/2` for every Sturmian **and quasi-Sturmian** word gives `Dio >= 2.50994`. This settles (a) for every map with `eta < 2.50994`, including every slope of `5x+1`, and the threshold is sharp.
   * (b) Rotation codings get the chain bound `Dio >= mu*(k0, L)` of Lemma J', and (c) substitution fixed points get finite certificates (Lemma F). In both classes the guaranteed exponent drops below `mu` for some parameters; those cases are **marked**. Every tested member under `3x+1` was nevertheless certified.
   * The barrier is real outside these classes. `Y` (zero entropy, width 1) has `Dio(Y) = 1`. Yet `Phi_2(Y)` is a 2-adic theta value, and a 2-adic Tschakaloff–Padé argument proves it irrational (Theorem Y).
   * The cube-swap word `Y3` has `Dio = 1` and a cubic-theta Bernstein number, so it is beyond both mechanisms: OPEN.
2. **C2, the open side.** No proof.
   * Capacity dies on exponential growth, and Theorem R dies because strip words have `Dio = 1` (almost surely for random ones, Proposition Rnd). Monks–Yazinski and Tao-type density see nothing.
   * Proposition M makes the (D1) reduction exact. With the corrected sign `E_s = -Phi_R(sigma^s w) > 0`, C2 for integers is **equivalent** to a *coupled* Z-number statement: no `L != 0` and no strip word `w` with `frac(L 3^(a_s)/2^s) = frac(E_s(w))` for all `s`.
   * The decoupled Mahler relaxation MZ implies C2. MZ is vacuous for wide strips, and non-vacuous and heuristically true for narrow ones (`P_MZ < 0` for slope `9/10`, `W <= 3/2`).
   * FLP, Dubickas and Koksma do not apply.
   * Proved partial results:
     * C2 on every strip word with `Dio > mu`, in particular linear complexity `c < 1/(mu-1)`;
     * C2 on square-swap words (Theorem Y);
     * two DEFECT-blind calibrations.
   * For positive integers alone, nothing new.
3. **The missing mechanism.** It is not about entropy. Periodic approximants with natural heights reach exactly `Dio > eta`, and `Y` and `Y3` are zero-entropy words below that line. So C2 is not the smallest instance: single explicit words are. `Y` fell to a genuinely different mechanism: Padé approximants from a q-difference equation, whose rationals have anomalously small heights. The first word we could not settle is `Y3` (HC1). Section 4 gives the paragraph.

## 1. Theorem D: the Diophantine exponent of the word against the height exponent of the map

### 1.1 Setting

We use the notation of the parent note, section 4.3.
* `T(x) = (m_e x + r_e)/2` for `x = e (mod 2)`, with `m_e` odd and `r_e = e (mod 2)`.
* `T^|z|(x) = (M_z x + R_z)/2^|z|` on the cylinder of `z`.
* `Phi_T(u v^inf) = P/Q`, where `Q = M_u D_v`, `P = 2^|u| R_v - R_u D_v` and `D_v = 2^|v| - M_v`.
* The gain is `g(u,v) = lambda(w, u v^inf) - log2(|P| + |Q|)`.
* Theorem R: if `g` is unbounded, then `Phi_T(w)` is not rational.

**Height exponent.** `eta(w) := max(1, lim_n max_(|z|=n, z factor of w) (1/n) log2 M_z)`. The limit exists by subadditivity.
* For `m0 = 1` (the maps `3x+r`, `5x+1`, ...), `eta = max(1, beta+ log2 m1)`, where `beta+` is the upper Banach density of 1s.
* If `w` has a uniform frequency `beta` (for instance bounded discrepancy), then `eta = mu(beta) = max(1, beta log2 m1)`.

**Diophantine exponent** (Adamczewski–Bugeaud, ETDS 27 (2007), section 2 [P]). `Dio(w)` is the supremum of the `rho` for which there are prefixes `U_n V_n^(w_n)` with `|U_n V_n^(w_n)| >= rho |U_n V_n|` and `|V_n^(w_n)|` strictly increasing.

The engine computes it as `Dio(w) = 1 + limsup_j LPF(j)/j`, where `LPF(j) = max_(a<j) LCE(a, j)`. For `U = w[0,a)` and `V = w[a,j)`, the common prefix of `w` and `U V^inf` has length `j + LCE(a,j)`, and the maximising `a` gives the best pattern with `|UV| = j`. Bugeaud–Kim's function is `r(n) = n + min{j : LPF(j) >= n}`.

### 1.2 Lemma H' (heights)

For every `eps > 0` there is a `C_eps` with
`log2(|P| + |Q|) <= (eta + eps)(t + p) + log2(t + p + 1) + C_eps`
for all `u = w[0,t)` and `v = w[t, t+p)`. If `log2 M_z <= eta|z| + C` for all factors `z` (bounded discrepancy), then `eps = 0` is allowed.

*Proof.*
1. By definition of `eta`, `M_z <= C_eps 2^((eta+eps)|z|)` for every factor `z`.
2. From `R_(ze) = m_e R_z + r_e 2^|z|` we get `R_z = sum_(i<|z|) r_(z_i) 2^i M_(z(i,|z|))`, where `z(i,|z|)` is the suffix after position `i`.
3. Each term is at most `max|r| C_eps 2^((eta+eps)|z|)`, because `eta + eps >= 1`.
4. So `|R_z| <= C' |z| 2^((eta+eps)|z|)`, and the bound follows from `|D_v| <= 2^p + M_v`, `M_u M_v = M_(uv)` and `|P| <= 2^t|R_v| + |R_u||D_v|`. ∎

(For words with bounded discrepancy this is the parent note's Lemma H.)

### 1.3 The theorem

**Theorem D.** Let `w` be a word that is not eventually periodic (EP). If `Dio(w) > eta(w)`, then `Phi_T(u' w)` is irrational for every finite `u'`.

*Proof.*
1. Take `eta < rho < Dio(w)` and `eps = (rho - eta)/2`. There are prefixes `U_n V_n^(e_n)` with `|U_n V_n^(e_n)| >= rho |U_n V_n|` and `|V_n^(e_n)|` strictly increasing.
2. `|U_n V_n| -> infinity`. Otherwise infinitely many `n` share one pair `(U, V)` while `e_n -> infinity`, and then `w = U V^inf` is EP.
3. With `u = U_n` and `v = V_n`, the word `U_n V_n^(e_n)` is a prefix of both `w` and `u v^inf`, so `lambda >= rho |uv|`.
4. By Lemma H', `g(u, v) >= eps |uv| - log2(|uv| + 1) - C_eps`, which tends to infinity.
5. Theorem R now shows that `Phi_T(w)` is not rational.
6. For `u'w`: `T^|u'|` maps `Phi_T(u'w)` to `Phi_T(w)` and preserves rationality. ∎

### 1.4 The exact repetition condition, with options A and B

Write a prefix pattern as `U V^e`, with `|U V^e| = |U| + e|V|`. Then

`|UV^e| >= rho |UV|` is equivalent to `e >= rho + (rho - 1)|U|/|V|`.

So Theorem D asks for infinitely many prefixes `U V^e` with `e >= mu + (mu - 1)|U|/|V| + delta (1 + |U|/|V|)` for a fixed `delta > 0`. For bounded-discrepancy words the logarithmic form is enough: `sup (|UV^e| - mu|UV| - log2(|UV|+1)) = +infinity`. It matters only at the threshold. The options of Theorem S (parent note section 4.3, `q = q_n`, with `j0 < j1` the first defects) are the following.
* **Option A:** `U = empty`, `V = s[0, q)`. It needs `e_A >= (j0 - 1 + q)/q > mu + delta`, i.e. an initial power of exponent `> mu`; `ice(w) > mu` suffices.
* **Option B:** `U = s[0, j0]`, `V = s[j0+1, j0+q]`. Here `e_B >= (j1 - j0 - 2 + q)/q >= (q_(n+1) + q - 2)/q`, and it needs `e_B > mu + (mu - 1)(j0 + 1)/q + delta(...)`.
* If A fails at scale `n`, then `j0 <= (mu - 1) q + 1`, and B holds once `q_(n+1) > mu(mu - 1) q_n + O(1)`. Theorem S's hypothesis `mu(mu - 1) < limsup q_(n+1)/q_n` is exactly "A or B at infinitely many scales". In `Dio` language it gives `Dio(s) >= mu_S(L)`, the root of `mu(mu - 1) = L`, which is at least `1.8668`.
* Lemma J' (section 2.2) extends this to a chain A, B1, ..., B_(k0).

### 1.5 What periodic approximants cannot do

1. **Proposition L (only periodic approximants, unless PC fails).** Let `r != Phi_T(w)` be a rational with odd denominator. Then `v_2(Phi_T(w) - r) = lambda(w, w(r))`, where `w(r)` is the parity vector of `r` (isometry). Under PC, `w(r)` is EP. So every rational that a Liouville argument can use is `Phi_T` of an EP word agreeing with `w` for `lambda` letters.

   **Caution.** Its height may be far below the *natural* height `2^(eta|uv|)`. The Dio barrier is exact only for certificates whose heights obey Lemma H' with the matching lower bound, i.e. without large gcd reductions.
   * In every certificate we computed, the canonical (minimal) representation had `gcd(P, Q) <= 2^8` (`.out`, "canon.gcd").
   * Theorem Y (section 2.4) is the counterexample to extrapolating this: its Padé rationals have heights about `2^(0.95 lambda)`, against natural heights of about `2^(1.43 lambda)`.
2. **The Subspace refinement has the same threshold.** Take `S = {inf, 2, 3}` and `X = (M_u 2^p, M_u M_v, P)`. Use the forms `X1, X2, X3` at `inf` and `3`, and `X1, X2, xi(X1 - X2) - X3` at `2`. Every S-unit cancels, the product is `|P| |P|_3 2^(-lambda)`, and the Subspace Theorem needs `lambda >= (1 + eps) log2 |P|`, which is again `Dio > eta`. Hence C1 (transcendence) has the same reach as Theorem D. Adamczewski–Bugeaud's p-adic Theorem 6 (Annals 2007 [P]) concerns Hensel digit expansions, whose natural height exponent is 1. It does not apply to `Phi_T`.
3. **Proposition T (two places; PROVED).** Let `w` have bounded discrepancy around `alpha` with `alpha log2 m > 1`, let `w` be non-EP, and let `Dio(w) > 2`. Then `Phi_2(w)` and the real value `Phi_R(w)` of the same series are not both rational.

   *Proof.*
   1. A periodic approximant `r = Phi(u v^inf)` equals both its 2-adic and its real sum (section 3.2, (iv)).
   2. `|Phi_2(w) - r|_2 = 2^(-lambda)` and `|Phi_R(w) - r| <= C 2^(-(eta - 1) lambda)`.
   3. Suppose `x = Phi_2(w)` and `y = Phi_R(w)` are both rational. Then `x != y` (section 3.2, (iv)). For `N = (x - r)(y - r) != 0`, the product formula gives `|N|_inf |N|_2 >= 1/(odd part of den N) >= 1/(C Q^2)`.
   4. The left side is `<= C' 2^(-eta lambda)`. So `lambda <= 2|uv| + O(log|uv|)`, which contradicts `Dio > 2`. ∎

   For `3x+1` this is weaker than Theorem D. It says something new only for maps with `eta > 2`. The Subspace version would need both values algebraic (HC4).

## 2. The classes (C3)

### 2.1 Sturmian and quasi-Sturmian words

**Corollary D1 (PROVED).** Let `w = u S(s)`, where `s` is Sturmian, `S` a morphism with `S(01) != S(10)`, and `u` finite. This includes every eventually Sturmian word and every recurrent quasi-Sturmian word (Cassaigne, DLT'97 [R, statement via Damanik–Lenz, arXiv math-ph/0105034, Prop. 2.2, read]). Then `Phi_T(w)` is irrational for every affine 2-adic shift map with `eta(w) < 5/3 + 4 sqrt10/15 = 2.5099407`.

*Proof.*
1. Bugeaud–Kim [P]: Thm 3.4 gives `rep(s) <= sqrt10 - 3/2` for every Sturmian word, and Thm 3.6 gives the same for every quasi-Sturmian word (`p(n) = n + k` for `n >= n0`).
2. Lemma 10.3 [P]: `rep = 1` iff `Dio = infinity`, `rep = infinity` iff `Dio = 1`, and otherwise `rep = Dio/(Dio - 1)`. Hence `Dio(w) >= (sqrt10 - 3/2)/(sqrt10 - 5/2) = 5/3 + 4 sqrt10/15`.
3. `S` is non-erasing, since `S(01) != S(10)`. So factors of `S(s)` have bounded discrepancy about `beta = (|S(0)|_1 (1-theta) + |S(1)|_1 theta)/(|S(0)|(1-theta) + |S(1)| theta)`, where `theta` is the slope of `s`, and `eta = mu(beta)`.
4. Apply Theorem D to `S(s)`; the prefix `u` is handled by `T^|u|`. ∎

A route that avoids Cassaigne's theorem and Thm 3.6: `S` maps every prefix pattern `U V^e` of `s` to the prefix pattern `S(U) S(V)^(e')` of `S(s)`. Since `|S(z)| = c|z| + O(1)` on factors of `s`, the ratios are preserved asymptotically, so `Dio(S(s)) >= Dio(s) >= 2.50994` by Thm 3.4 alone.

**Consequences.**
* `3x+r` (`eta <= log2 3`) and Mahler's map (`eta = log2 3`): every slope and every intercept.
* `5x+1` (`eta <= log2 5 = 2.3219`): **every slope**. The parent note had `alpha < 0.804` (Theorem S) and `alpha < 0.8614` (ADQZ route).
* `7x+1`: `beta < 0.8941`. `9x+1`: `beta < 0.7918`. `19x+1`: `beta < 0.5909`.
* Quasi-Sturmian words are a genuinely new class. Example (`.out` H2): `q = 00 S(s)` with `S(0) = 1101`, `S(1) = 111` and `s` the golden Sturmian word. It has `p(n) - n = 5` for `11 <= n <= 40` and `beta = 0.887`, and it is PROVED under `3x+1` (`mu = 1.406`), `5x+1` (`mu = 2.060`) and `7x+1` (`mu = 2.490`).

**Corollary D2 (PROVED, per word).**
* (a) If the slope has unbounded partial quotients, then `Dio = infinity` (Bugeaud–Kim Thm 3.3 [P]; Adamczewski–Bugeaud 2011, Prop. 11.1 [P]), and every map is covered.
* (b) For the characteristic word `c_alpha`, and for `0 c_alpha` (the lower mechanical word with intercept 0, checked on `>= 1.5*10^5` letters), `Dio >= ice(c_alpha) = 1 + limsup q_(n+1)/q_n >= 1 + phi = 2.618`. This is Berthé–Holton–Zamboni 2006, Thm 1.2 and its proof [P], with `ice <= Dio`.
  * Example: slope `[0;1,8,(1)]` under `7x+1` has `mu = 2.5155`, beyond every uniform bound, and is PROVED.
  * For the intercept `1/3` at the same slope, `Dio_est = 3.49`, but that is FINITE-EXACT only.
* (c) Per slope, `mu(mu - 1) < limsup q_(n+1)/q_n` (Theorem S's options). Example: `[0;1,(10)]` under `7x+1`.

**Sharpness (PROVED + FINITE-EXACT).** Bugeaud–Kim's `s'` has slope `(sqrt10-2)/3 = [0;(2,1,1)]` and intercept `1/3`. Its `rep` equals `sqrt10 - 3/2` (their (3.1)), so `Dio(s') = 2.50994` exactly.
* The engine gives `rep = 1.66229` and `Dio = 2.50991` on `2^20` letters.
* Under `19x+1`, the complement of `s'` (`beta = 0.61257`) has `eta = 2.6022 > Dio`. Theorem D gives nothing, and the exact gains of the best approximants are all negative (best `-8.3` bits).
* So `2.50994` is the exact uniform threshold for Sturmian words in periodic-approximant methods.
* Whether `Phi_(19x+1)` of that word is irrational is OPEN (HC5).

### 2.2 Codings of a rotation by finitely many arcs (Lemma J')

**Setting.**
* `alpha` is irrational; `f: R/Z -> {0,1}` is right-continuous with a finite discontinuity set `Delta` of `k` points; `x_i = f(rho + i alpha)`.
* Split `Delta` into `k0` classes modulo `Z alpha + Z`. Choose each class's base point so that its members are `e + m alpha` with offsets `0 <= m <= M`.
* `L = limsup q_(n+1)/q_n`, where `p_n/q_n` are the convergents of `alpha`.
* `mu*(k0, L)` is the root `> 1` of `mu(mu^k0 - 1) = L`, and `mu* = infinity` if `L = infinity`.

**Lemma J' (PROVED).** `Dio(x) >= max(1 + 1/k, mu*(k0, L))`.

*Proof.*
1. **The complexity bound.** Factors of length `n` correspond to the arcs cut out by the at most `kn` points `Delta - j alpha` (`0 <= j < n`), so `p(n) <= kn`. Then `Dio >= 1 + 1/k` by the pigeonhole argument in the proof of Theorem 3 of Adamczewski–Bugeaud, ETDS 2007 [P].
2. **Defect arcs.** Fix `n`, `q = q_n` and `delta = q alpha - p_n`. If `x_(i+q) != x_i`, a discontinuity `e` lies in the half-open arc between `y_i = rho + i alpha` and `y_i + delta`, i.e. `y_i` lies in an arc `A_e` of length `|delta|`. Let `J_e = {i in Z : y_i in A_e}`.
   * (J1) Distinct elements of `J_e` differ by at least `q_(n+1)` (best approximation).
   * (J2) If `e' = e + m alpha`, then `J_(e') = J_e + m`.
   * (J3) Every `q_n + q_(n+1)` consecutive integers meet `J_e` (three-gap theorem, as in Lemma J).
   * So the defects of a class lie in blocks `[a, a + M]` with `a in J_(e_C)`.
3. **Components.** Let `[c_t, c'_t]` be the components of the union of the blocks that meet `[0, infinity)`. Once `q_(n+1) > k0(M+1)`, a component has length at most `M' = k0(M+1)`, by (J1).
4. **Pigeonhole.** Among the first `k0 + 1` components, two anchors belong to the same class, so `c_(k0+1) >= c_1 + q_(n+1) - M'`.
5. **The options.** Put `U = x[0, c'_t + 1)` and `V = x[c'_t + 1, c'_t + 1 + q)`, plus the option with `U` empty. They give ratios `r_0 = (max(c_1, 0) + q)/q` and `r_t >= (c_(t+1) + q)/(c_t + M' + 1 + q)`.
6. **The chain.** If all `r_t <= rho` for `t <= k0`, induction gives `c_(k0+1) + q <= rho^(k0) (c_1 + q) + C`, hence `q_(n+1) <= rho(rho^(k0) - 1) q_n + C'`.
7. So if `L > rho(rho^(k0) - 1)`, infinitely many scales give a pattern with ratio `> rho` and `|V^e| >= q_n -> infinity`. ∎

For Sturmian words, `k0 = 1` (the class `{0, -alpha}`, `M = 1`), and Lemma J' is Theorem S's `mu(mu-1) < L`.

**Constants** (`.out` H0):

| `k0` | `mu*(k0, phi)` | `mu*(k0, 2)` | `mu*(k0, 1+sqrt2)` | `mu*(k0, 3)` |
|---|---|---|---|---|
| 2 | 1.4536 | 1.5214 | 1.5876 | 1.6717 |
| 3 | 1.3079 | 1.3532 | 1.3972 | 1.4526 |

**PROVED instances.**
* Every rotation coding whose slope has unbounded partial quotients, under every map (`L = infinity`).
* Every 2-class coding under `3x+r` whose partial quotients are not eventually in `{1, 2}` (then `L >= 3`, and `mu*(2, 3) = 1.6717 > log2 3`).
* For eventually-golden slopes with `k0 = 2`, `beta < 0.917`. Example: the golden coding with 1-set `[0, 9/10)` under `3x+1` (`mu = 1.4265 < 1.4536`).
* FINITE-EXACT, per convergent `q_n <= 17711`: the best option ratio is always `>= mu*(2, q_(n+1)/q_n)` (`.out` H3).

**MARKED.** With bounded partial quotients and `k0 >= 3` classes, the guarantee falls below `mu`. For example, `mu*(3, phi) = 1.3079 < mu(0.9) = 1.4265`, and the complexity bound `1 + 1/k` is even weaker. Here the needed exponent exceeds what the class guarantees.
* **Census** (FINITE-EXACT, not a proof). 160 golden codings, each with a 0-set made of 2–6 arcs `[c_i, c_i + {m alpha})` (bounded discrepancy) at rational `c_i`. Under `3x+1`, `Dio_est > mu` in all 160, with minimum margin `+0.313` and median `Dio_est = 2.364`. Under `5x+1`, 124 of 160.

### 2.3 Linear complexity, linearly recurrent words, primitive substitutions

* **Linear complexity (PROVED).** If `p(n) <= c n` for infinitely many `n`, then `Dio >= 1 + 1/c` (Adamczewski–Bugeaud, ETDS 2007, proof of Thm 3 [P]). Hence PC holds for such words whenever `c < 1/(eta - 1)`. Under `3x+1`: `c < 1/(beta log2 3 - 1)`, i.e. every `c` near the critical slope, `c < 2.345` at `beta = 0.9`, and `c < 1.71` at `beta = 1`.
* **Linear recurrence (PROVED, MARKED).** Take the definition of Durand [P, arXiv 0807.4430, Def. 2]: every return word to `u` has length at most `K|u|`. The prefix of length `n` returns within `Kn`, so `ice >= 1 + 1/K` and `Dio >= 1 + 1/K`. Primitive substitutive words are linearly recurrent (Durand–Host–Skau [R via Durand]).
  * The guarantee `1 + 1/K` is below `mu` once `K > 1/(mu - 1)`: **marked**.
* **Lemma F (PROVED).** Let `sigma` be primitive, with Perron–Frobenius eigenvalue `theta > 1`, positive left eigenvector `l` (`l M = theta l`, where `M_(bc) = |sigma(c)|_b`), and aperiodic fixed point `x = sigma^inf(a)`. Write `l(z) = sum_b |z|_b l_b`. Then

  `Dio(x) = sup { l(U V^e)/l(U V) : U V^e a prefix of x }`.

  *Proof.*
  1. (>=) If `U V^k V'` is a prefix (with `V'` a proper prefix of `V`), then so is `sigma^j(U) sigma^j(V)^k sigma^j(V')`, and `sigma^j(V')` is a prefix of `sigma^j(V)`.
  2. Perron–Frobenius gives `|sigma^j(z)| = theta^j (c l(z) + o(1))`, so the length ratios tend to `l(UV^e)/l(UV)`, while `|sigma^j(V)^...| -> infinity`.
  3. (<=) Primitive substitutive words have uniform frequencies, so along long patterns `l(z)/|z|` tends to a constant uniformly. ∎

  **Consequence.** Every prefix pattern with `l`-ratio `> eta` is a *finite certificate*: PC holds on that fixed point, so the method is semi-decidable per substitution. Certificates transfer to letter-to-letter codings, since `Dio(tau(x)) >= Dio(x)`, and to eventual versions.
* **Census (FINITE-EXACT; each certified word is PROVED).**
  * There were 9,716 random primitive binary substitutions with images of length `<= 16` and aperiodic supercritical fixed points; 11 periodic ones were excluded by Morse–Hedlund.
  * Under `3x+1`, every one has a certificate from its first 300 letters, with minimum margin `+0.295`. Of these, 3,944 have `|theta_2| < 1`, so they have bounded discrepancy (Adamczewski, TCS 307 (2003), Thm 13 [P]) and lie in C2 strips.
  * Example: `sigma(0) = 0110010100`, `sigma(1) = 11110111111`. It has `beta = 0.828`, `mu = 1.313`, and the certificate `1.6081` comes from the pattern `U = x[0,158)`, `V = x[158,169)`. It is PROVED. Its exact gains grow (2,819 bits at `|UV| = 8328`), and `|theta_2| = 5.17`, so its discrepancy grows like `N^0.69`.
* **The drift side.** Under `5x+1`, 99 of the 9,716 have no certificate; the same `sigma` has `mu = 1.9235 > 1.6081`, and its exact gains are negative. These are explicit linearly recurrent words for which the method is limited, on the side where PC is expected to fail.

### 2.4 Beyond the Dio barrier: the square-swap word Y (Theorem Y) and the cube-swap word Y3

**Definition.** `Y = b_0 b_1 b_2 ...` with `b_n = 1^8 0 1` if `n` is a nonzero square, and `b_n = 1^9 0` otherwise. `Y3` is the same with nonzero cubes.

**Proposition Y (PROVED).** Let `Z` be `Y` or `Y3`.
* (i) Every block has nine 1s, so `a_s - 0.9 s` lies in `[-0.1, 0.9]` (width 1). `Z` lies in the supercritical strip of slope `9/10` (under `3x+1`, `mu = 1.4265`).
* (ii) `p_Y(n) <= n^2/40 + 2n + 50` (zero entropy), and similarly `p_(Y3)(n) = O(n^(3/2))`.
* (iii) `Z` is not EP.
* (iv) `Dio(Z) = 1`, with `LPF(j) <= 13 sqrt(j) + 60` for `Y` and `LPF(j) = O(j^(2/3))` for `Y3`.

*Proof of (iv) for `Y`.* Let `U V^e` be a prefix with `a = |U|` and `p = |V|`, so that `Y[i] = Y[i+p]` on `I = [a, a + l)`. The zeros of `Y` sit at `10n + 9` (normal block) or `10n + 8` (swapped block).
* **Case `p != 0, +-1 (mod 10)`.** A zero in `I` must map to a zero, so `I` contains no zero and `l <= 11`.
* **Case `p = 1 (mod 10)`.** Every zero in `I` belongs to a swapped block. Between two swapped blocks there are normal zeros, so `I` contains at most one zero and `l <= 22`.
* **Case `p = 9 (mod 10)`.** Every zero in `I` is normal and maps to a swapped block. Two consecutive normal blocks would map to two squares differing by 1 or 2, which is impossible. So `l <= 22`.
* **Case `p = 10m`.** For blocks inside `I`, `b_n = b_(n+m)`. If `I` contained two consecutive swapped blocks `k^2` and `(k+1)^2`, then `k^2 + m = K^2` and `(k+1)^2 + m = K'^2` with `K > k`. Then `2k + 1 = K'^2 - K^2 >= 2K + 1 > 2k + 1`, a contradiction. So `I` contains at most one swapped block. With `k0` the least `k` such that `10 k^2 >= a`, this gives `l <= 40 k0 + 20 <= 13 sqrt(a) + 60`.
* **Conclusion.** Every pattern has ratio at most `1 + (13 sqrt(j) + 60)/j -> 1`.

For `Y3`, use `K'^3 - K^3 >= 3K^2 + 3K + 1 > 3k^2 + 3k + 1` and note that cubes differ by at least 7. For (ii), a window meeting two swapped blocks starts before `10(n/20 + 1)^2`. ∎

FINITE-EXACT check: the maximiser gives `LPF(j)/sqrt(j) -> 12.64`, against the proof's `40/sqrt10 = 12.65`, so the bound is sharp. The profile decays like `1 + 12.6/sqrt(j)`. The best exact Theorem-R gain for `Y` is `77.6` bits at `|UV| = 280`, and the gains are negative for all `|UV| >= 1030`. Theorem D does not apply.

**The block identity (PROVED; checked exactly modulo `2^30000`).**
1. `T^10` acts on the two blocks by `x -> (3^9 x + 3^9 - 2^9)/2^10` and `x -> (3^9 x + 3^9 - 2^8)/2^10`.
2. Inverting 2-adically, with `rho = 2^10/3^9` (`|rho|_2 = 2^-10`), gives `Phi_2(Y) = -(1/3^9) sum_n c_n rho^n`, where `c_n = 19171` if `n` is not a square and `19427` if it is:

   `Phi_2(Y) = -1 - 512/(3^9 - 2^10) - (256/3^9) sum_(k>=1) rho^(k^2)`.
3. So `Phi_2(Y)` is rational if and only if `theta(rho) = sum_(k>=0) rho^(k^2)` is rational in `Q_2`.
4. The same expression summed in `R` is `Phi_R(Y) = -1.0281165761...`. By Nesterenko's theorem, theta values at algebraic `q` with `0 < |q| < 1` are transcendental ([R]: Duverney–Nishioka–Nishioka–Shiokawa, via the search results), so `Phi_R(Y)` is transcendental. The p-adic analogue of Nesterenko's theorem is not available.

**Lemma P (2-adic Tschakaloff; PROVED).** Let `rho = 2^L/M_0` with `M_0 > 1` odd, `L >= 1`, and `mu_bar = max(1, log2(M_0)/L) < phi = 1.618...`. Then `theta(rho) = sum_(k>=0) rho^(k^2)` is irrational in `Q_2`. If `mu_bar < 3/2`, the unshifted construction (`m = 0`) already suffices.

*Proof.* This is Zudilin's construction for Tschakaloff's series [P, arXiv math/0506086], transferred to `Q_2`.
1. **The linear form.** Put `q = rho^-2` (so `|q|_2 = 2^(2L) > 1`) and `z = rho`. Then `theta = T_q(z) = sum_l z^l q^(-l(l-1)/2)`. Let `R_n(T) = prod_(j=1..n) (1 - q^j T) = sum_k C_k T^k`, with `C_k in Z[q]` and `sum |coefficients| = 2^n`. Put `m = round(n/phi)` or `m = 0`, and

   `I_n = sum_(t>=1) R_n(q^-t) z^(t+m) q^(-(t+m)(t+m-1)/2) = A_n theta - B_n`,

   where `A_n = sum_k z^-k C_k q^(k(k-1)/2 + km)` and `B_n = sum_k z^-k C_k q^(k(k-1)/2 + km) sum_(l <= k+m) z^l q^(-l(l-1)/2)`. (Reindex with `l = t + m + k`.)
2. **Its 2-adic size.** `R_n(q^-t) = 0` for `1 <= t <= n`, and `|R_n(q^-t)|_2 = 1` for `t > n`, because `|q^(j-t)|_2 < 1`. The `t`-th term has size `2^(-L(t+m)^2)`, so `|I_n|_2 = 2^(-L(n+m+1)^2)` exactly, and `I_n != 0`.
3. **Clearing denominators.** Every monomial is `+-rho^e` with `e in [-e_off, m^2]`, where `e_off = n + 2n(n+m)`. `N_n = 2^(L e_off) M_0^(m^2)` makes `N_n A_n` and `N_n B_n` integers of size `<= 2^n (n+m+1) max(2^L, M_0)^(e_off + m^2)`.
4. **The contradiction.** If `theta = a/b`, then `G = a N_n A_n - b N_n B_n = b N_n I_n` is a nonzero integer with `v_2(G) >= L(e_off + (n+m+1)^2)` and `|G| <= (|a|+|b|) H_n`. So
   `L n^2 [(3 + 4x + x^2) - mu_bar (2 + 2x + x^2)] + O(n) <= log2(|a|+|b|)`, with `x = m/n`.
5. The bracket is positive for `mu_bar < sup_x (x^2+4x+3)/(x^2+2x+2) = phi`, attained at `x = 1/phi`. At `x = 0` it is positive for `mu_bar < 3/2`. ∎

This is the p-adic transcription of Zudilin's condition `gamma < (3-sqrt5)/2`, with `1/mu_bar` in the role of `1 - gamma`. The p-adic Tschakaloff series has been studied before, for example by Väänänen and Wallisser, J. Number Theory 39 (1991) 225–236 ("linear independence measure for ... `sum q^(n(n-1)/2) z^n`, `q in Q`, `0 < |q|_p < 1`"). We read only the metadata and a search snippet, so whether their hypotheses cover `rho = 2^10/3^9` is UNVERIFIED. The proof above is self-contained.

FINITE-EXACT (`.out` H5): the linear forms were computed exactly for `n <= 14`.
* `v_2(N_n I_n)` equals the predicted `L(e_off + (n+m+1)^2)` in every case.
* The margins `v_2 - log2 H_n` grow from `+27` to `+518` (`m = 0`) and from `+54` to `+1798` (`m = round(n/phi)`).

**Theorem Y (PROVED).** Let `B != B'` be blocks of equal length `L` with the same number `A` of 1s, and let `T = x/2, (mx + r)/2` with `mu_bar = max(1, A log2 m / L) < phi`. Then no rational has an eventually `Y_(B,B')` parity vector, where `Y_(B,B')` uses `B'` at the nonzero squares and `B` elsewhere.

In particular this holds for `Y` under `3x+1` (`mu_bar = 1.4265`), and for **every** square-swap block word under `3x+r`, since `mu_bar <= log2 3 < phi`.

*Proof.* The block identity gives `Phi_T(Y_(B,B')) = -(1/m^A) sum_n c_n rho^n`, where `rho = 2^L/m^A`, `c_n in {R_B, R_(B')}` and `R_B != R_(B')` (distinct words of the same length have distinct cylinders). So rationality is equivalent to the rationality of `theta(rho)`, and Lemma P applies. `T` maps rationals to rationals and back, so the eventual version follows. ∎

Checked for `5x+1` (`.out` H5):
* `B = 1^6 0^4` has `mu_bar = 1.393 < phi`, and its margins are positive.
* `B = 1^8 0^2` has `mu_bar = 1.858 > phi`, and its margins are negative (`-591`): Lemma P does not apply.

**`Y3`: the smallest open instance found.** The same identity with cubes gives

`Phi_2(Y3) = -1 - 512/(3^9-2^10) - (256/3^9) sum_(k>=1) rho^(k^3)`

(exact modulo `2^30000`). Neither mechanism reaches it.
* `Dio(Y3) = 1` (Proposition Y). The exact gains are positive only up to `|UV| ~ 2.7*10^4`, as predicted by `LPF = O(j^(2/3))`, and negative beyond.
* The partial sums have heights `2^(14.26 K^3)` against 2-adic errors `2^(-10(K+1)^3)`, so Liouville fails; it would need `log2(3^9)/10 < 1`.
* `sum z^k rho^(k^3)` satisfies no first-order q-difference equation. Only a two-variable one exists: `G(u,v) = 1 + u v rho G(u v^2 rho^3, v rho^3)` for `G = sum u^k v^(k^2) rho^(k^3)`. So Lemma P's construction does not apply.
* Even in `R`, for `q = 1/(integer)`, cubic theta values are known only to have algebraic degree `>= 4` (Ghidelli, arXiv 1910.05076 [P, abstract]).

**Status: OPEN (HC1).**

### 2.5 Classification table

`beta` is the frequency of 1s and `mu = max(1, beta log2 m1)`. "Condition" is what Theorem D needs.

| word class | guaranteed `Dio` (source) | condition | `3x+r` | `5x+1` | `7x+1` |
|---|---|---|---|---|---|
| eventually Sturmian | `>= 2.50994`, sharp (BK 3.4 + 10.3) | `mu < Dio(s)` | PROVED, all slopes and intercepts (and Mahler) | **PROVED, all slopes** | PROVED for `beta < 0.894`; per word beyond |
| Sturmian, unbounded partial quotients | infinity (BK 3.3; AB 2011, Prop 11.1) | none | PROVED | PROVED | PROVED |
| Sturmian, characteristic (`c_alpha`, `0 c_alpha`) | `>= 1 + limsup q_(n+1)/q_n >= 2.618` (BHZ Thm 1.2) | `mu < 2.618` | PROVED | PROVED | PROVED for `beta < 0.9326` |
| Sturmian, per slope (options A/B) | root of `mu(mu-1) = L` | `mu(mu-1) < L` | – | – | e.g. `[0;1,(10)]` |
| eventually quasi-Sturmian `u S(s)` | `>= 2.50994` (BK 3.6 + 10.3) | `mu < 2.50994` | PROVED | PROVED | PROVED for `beta < 0.894` |
| `p(n) <= cn` infinitely often | `>= 1 + 1/c` (AB ETDS 2007) | `c < 1/(mu-1)` | PROVED if `c < 1/(beta log2 3 - 1)` | PROVED if `c < 1/(beta log2 5 - 1)` | same shape |
| rotation coding, `k` endpoints, `k0` classes | `>= max(1+1/k, mu*(k0,L))`; infinity if `L` infinite (Lemma J') | `mu <` that | PROVED under the condition; **MARKED** for bounded partial quotients and `k0 >= 3`; census 160/160 empirically `> mu` | MARKED; 124/160 | MARKED |
| linearly recurrent, constant `K` | `>= 1 + 1/K` (Durand) | `mu < 1 + 1/K` | **MARKED** for `K > 1/(mu-1)` | MARKED | MARKED |
| primitive substitution fixed point, primitive morphic, eventually | exact: sup of `l`-weighted prefix ratios (Lemma F) | a finite certificate `> mu` | PROVED per word; census 9,716/9,716 certified; class-wide OPEN (HC2) | 9,617/9,716; explicit method-limited members | per word |
| square-swap block words (zero entropy) | `Dio = 1` (Prop Y): Theorem D fails | – | **PROVED by Lemma P** (Theorem Y) | PROVED iff `mu_bar < phi` (Lemma P) | same |
| cube-swap `Y3` (zero entropy, width 1) | `Dio = 1` | – | **OPEN (HC1)** | OPEN | OPEN |
| positive-entropy strip words (C2) | `= 1` almost surely (Prop Rnd) | – | **OPEN (C2)** | OPEN (expected false somewhere: drift) | OPEN |

### 2.6 Numerics for part 1 (`.out` H1–H6)

* **Engine calibration.** Fibonacci: `Dio = 2.618028`, `rep = 1.61804`. Thue–Morse: `5/3` and `5/2`. Bugeaud–Kim `s'`: `2.509912` and `1.66229`. Random: `1.005`.
* **Exact certificates.** Every certificate has the Bernstein residue validated by parity iteration, and `v_2(Phi - P/Q) = lambda` exactly.
  * Sturmian `[0;1,7,(1)]` under `5x+1` and `7x+1`, where only Bugeaud–Kim applies: best gains `7781.7` and `1876.2` bits (intercept 0), `12232.0` and `8582.2` (intercept `1/3`).
  * `[0;1,8,(1)]` under `7x+1` (BHZ, characteristic): `1573.9` bits.
  * `s'` under `19x+1`: `-8.3` bits (sharpness).
  * Quasi-Sturmian `q` under `3x+1`, `5x+1` and `7x+1`: `17137.0`, `7894.7` and `1806.4` bits.
* **Reconstruction at `N = 10^4` bits.** Every tested word gives either no candidate, or a candidate that leaves the word at letter `N` or later.

## 3. C2, the open side

**C2.** No rational has a non-EP parity vector with `sup_s |a_s - alpha s| < infinity` for some `alpha in (log_3 2, 1)`.

### 3.1 (a) Why every existing mechanism stops

* **Capacity.** The in-house proof needs `n_j <= C j` (D2c). A C2 orbit has `n_s >= c lambda^s`, with `lambda = 3^alpha/2 > 1`: `1.0212` at `alpha = 0.65` and `1.3439` at `0.9`. So `#{s : n_s <= X} <= log X / log lambda + C`.
  * The orbit has natural density 0 and even logarithmic density 0.
  * The counting inequality `N <= #{admissible integers <= C lambda^N}` always holds.
  * The density-zero stopping set `E_eps` contains the orbit without contradiction.
  * Proposition B's source-density step has nothing to bite on.
* **Theorem R.**
  * **Proposition Rnd (PROVED).** Consider a rational-slope strip automaton with no forced cycle (true for all 16 strips in `.out` C1). Choose letters independently at free states, fairly or with bias `1 - eps` versus `eps`. Then almost surely `w` is non-EP and `Dio(w) = 1`.

    *Proof.* Free states recur within `G` steps. A fixed pair `(a, j)` has a repetition of length `l` with probability `<= c^(l/G - 1)`, where `c < 1`. The union over `a < j` at `l = C log j` is summable, so Borel–Cantelli gives `LPF(j) = O(log j)`. ∎
  * So Theorem D cannot apply, and the gains are bounded.
  * Measured (`.out` C2), best exact gains:
    * random words: `18.7`, `18.5`, `30.5`, `32.1` and `15.4` bits;
    * greedy low-repetition words: `4.7`, `6.6`, `11.6` and `13.2` bits;
    * the Sturmian-perturbed word: `21.8` bits;
    * `Y`: `77.6` bits.
  * The per-scale estimates are negative from `|UV| ~ 2^7` to `2^9` on (for `Y`, from `2^10`).
* **Monks–Yazinski.** Their liminf bound `a_s/s >= log_3 2` holds for every C2 word (`alpha > log_3 2`), so it says nothing.
* **Tao-type density.** These are log-density or measure statements. A C2 orbit has log density 0, and the set of 2-adic integers with strip parity vectors is Haar-null. These statements are DEFECT-blind.

### 3.2 (b) The reduction: C2 is a coupled Z-number problem

Take `T(x) = x/2`, `(3x+1)/2`. Let `d_l` be the position of the `l`-th 1, `a_s = #{i < s : w_i = 1}` and `lambda_s = 3^(a_s)/2^s`.

**Proposition M (PROVED).** Let `w` be strictly supercritical (`liminf a_s/s > log_3 2`).
* **(i)** `S_inf(w) = sum_l 2^(d_l) 3^(-l)` converges in `R`. Put `Phi_R(w) := -S_inf(w)`: this is the real sum of Bernstein's series, whose 2-adic sum is `Phi_2(w)`. Put `E_s := S_inf(sigma^s w) = -Phi_R(sigma^s w) > 0`. Then `E_s = 2 E_(s+1)` if `w_s = 0`, and `E_s = 1/3 + (2/3) E_(s+1)` if `w_s = 1`.
* **(ii)** For every real `x`, `f_(w_(s-1)) o ... o f_(w_0)(x) = lambda_s (x - Phi_R(w)) - E_s`, where `f_0(y) = y/2` and `f_1(y) = (3y+1)/2`.
* **(iii)** Let `x = c/b` be rational (`b` odd) and `L = x - Phi_R(w)`. Then `x` has parity vector `w` iff `b(lambda_s L - E_s) in Z` for all `s >= 0`. For `x in Z`, this says `frac(L lambda_s) = frac(E_s)` for all `s`.
* **(iv)** Suppose `w` has bounded discrepancy. If `x` has parity vector `w` and `w` is not EP, then `L != 0` and `T^s(x) -> sign(L) * infinity`, and `x > 0` implies `L > 0`. If instead `w = u v^inf` is EP and supercritical, then `Phi_R(w) = Phi_2(w)` as rationals, so `L = 0`.

*Proof.*
* (i) `d_l <= l/alpha'` for large `l`, with `alpha' > log_3 2`, so the terms are at most `(2^(1/alpha')/3)^l`. The recursion comes from reindexing the tail.
* (ii) Induction: `y/2` and `(3y+1)/2` act on `lambda_s L - E_s` exactly as in (i).
* (iii) Integrality of `y_s` and `y_(s+1) = f_(w_s)(y_s)` forces `b y_s` to be even (for `w_s = 0`) or odd (for `w_s = 1`), so the parities are automatic. (For integers this is the key point: the parity bit is not extra data.)
* (iv) If `L = 0`, then `y_s = -E_s` is bounded, so the orbit is finite in `(1/b)Z` and hence EP. The EP case: the real and 2-adic sums are the same geometric-series expression, and it converges in both places because `v` is supercritical. ∎

FINITE-EXACT (`.out` C3): `(1)^inf -> -1`, `(110)^inf -> -5` and four more EP words give real sums equal to `P/Q` to within `1e-45`. The tail recursion holds to `1e-20`.

**Signs.** The task's `E_s` is the **negative** of the real Bernstein value of the tail, and `E_s > 0`. With `L = x + S_inf = x - Phi_R(w)`, one has `frac(L lambda_s) = frac(E_s) = frac(-Phi_R(sigma^s w))`.

**Equivalence.** C2 for integers (slope `alpha`, width `W`) holds **iff** there are no strip word `w` and no `L != 0` with `frac(lambda_s(w) L) = frac(E_s(w))` for all `s >= 0`. Here `lambda_s = lambda^s 3^(delta_s)`, `lambda = 3^alpha/2` and `delta_s = a_s - alpha s` is bounded. For rationals with denominator `b`, replace this by `frac(b lambda_s L) = frac(b E_s)`. For a given word, the only possible `L` is `Phi_2(w) - Phi_R(w)`, and it is a candidate only when `Phi_2(w)` is an integer. The multiplier sequence and the target are **coupled** through the same word.

**Decoupling (MZ).**
* **Definition.** Let `K_q = {E(w') : w' a strip tail from state q}`. MZ(`alpha, W`) says there are no `L != 0` and no strip path `(q_s)` with `frac(L lambda_s) in frac(K_(q_s))` for all `s`.
* **Implication (PROVED).** MZ implies C2 for integers on that strip, since `E_s in K_(q_s)`.
* **When MZ is vacuous (PROVED).** Graph-directed IFS: `E(1w) = 1/3 + 2E(w)/3` and `E(0w) = 2E(w)`. If the one-step images of the hulls overlap at every state, then `K_q = [lo_q, hi_q]`. If moreover every hull has length `>= 1`, then every pair (`L`, path) satisfies the condition, so MZ is false and the implication MZ ⟹ C2 carries no information.
* **The table** (`.out` C1). Put `m_q = min(1, |K_q|)` and `P_MZ = log2 lambda + log2 rho(adjacency weighted by m_q)`.
  * All `2/3` strips: MZ is vacuous (`m_q = 1`).
  * `3/4` and `4/5` strips of width `>= 1.5`: `P_MZ` between `+0.35` and `+0.89`, so MZ is heuristically false.
  * Width-1 strips of slope `3/4`, `4/5`, `9/10`, and the strip `9/10`, `W = 3/2`: `P_MZ = -0.478, -1.362, -3.514, -0.099`. These are upper bounds where 'overlap' fails. There MZ is non-vacuous and heuristically true.
* **So:**
  * **no fixed-target Mahler statement can imply C2 on wide strips;**
  * on narrow strips a single new ingredient would suffice (HC3): a Z-number theorem for the jittered sequence `3^(a_s)/2^s` with a state-dependent target of small measure. `Y` lives in the `(9/10, 1)` strip.

**Relation to FLP, Dubickas and Koksma.**
* Flatto–Lagarias–Pollington 1995 ([R], via the parent note): the spread of `frac(xi (p/q)^n)` is `>= 1/p`.
* Dubickas, Bull. London Math. Soc. 38 (2006) 70–80 ([P], abstract): a spread inequality for `frac(xi alpha^n)` with `alpha > 1` algebraic.

Both concern **fixed** multipliers with no target coupling. For rational `alpha = A/B`, along `s = kB` the multiplier is `theta^k 3^(e_k)` with `theta = 3^A/2^B` and a bounded integer jitter `e_k`. Even with the jitter frozen, FLP would give only a spread `>= 3^(-A)` (`3^(-9)` at `alpha = 9/10`) for a *fixed* target set. In C2 the target values `frac(E_(kB))` are set by the same word, so a lower bound on the spread of `frac(L lambda_s)` cannot contradict them unless it is coupled to `w`.

Koksma's metric theorem says that for fixed `w`, almost every `L` makes `frac(L lambda_s)` equidistributed. But by (iv), at most one `L` per word is a candidate: DEFECT-blind.

**Verdict.**
* C2 is **equivalent** to the coupled Z-number statement.
* It is **implied by** the decoupled MZ, which is vacuous for wide strips.
* It **neither implies nor is implied by** Mahler's Z-number problem, which concerns another map (THM-2228).
* It is **not implied by** FLP, Dubickas or Koksma.

### 3.3 (c) Partial results and conditional statements

1. **C2 on the reachable part (PROVED).** C2 holds for every strip word with `Dio(w) > mu(alpha)` (Theorem D). This includes:
   * words with `p(n) <= c n` infinitely often, `c < 1/(alpha log2 3 - 1)` (e.g. `c < 2.345` at `alpha = 0.9`);
   * Sturmian and quasi-Sturmian strip words;
   * certified substitutive words (3,944 with bounded discrepancy in the census);
   * all square-swap block words (Theorem Y), which have zero entropy and `Dio = 1`.

   So C2 contains PROVED zero-entropy families on both sides of the Dio barrier.
2. **Width.**
   * For irrational `alpha`, `W < 1` gives only Sturmian words (Theorem S).
   * For rational `alpha`, `W < 1` gives only EP words, and the closed width-1 strips already have positive entropy (`.out` C1: `0.1` bits/letter at `9/10`).
   * Every strip that allows free choices at a positive density of times contains words with `Dio = 1`. This follows from Proposition Rnd for rational slopes; the same argument works for irrational slopes, whose free times have density `W - 1`.

   So "C2 for `W < W0`" needs a new mechanism for every `W0 > 1`.
3. **Calibrations (PROVED, DEFECT-blind).**
   * For every non-atomic measure on words, almost every word has irrational `Phi_2`, because `Q cap Z_2` is countable.
   * C2(`alpha, W`) fails for at most **countably many slopes** `alpha`: a rational counterexample determines `w` and hence `alpha = lim a_s/s`.
4. **Positive integers only.** Nothing new. By Proposition M, integrality is equivalent to the coupled Z-number condition, and integrality plus growth yields only `L > 0`. The capacity argument needs polynomial growth and has none.
5. **Conditional statements.**
   * MZ(`alpha, W`) implies C2(`alpha, W`) (HC3).
   * The two-place Subspace statement (HC4) would force `Phi_R(w)` to be transcendental for every rational counterexample with `Dio(w) > 1`.
   * "Repetition forcing": if rational parity vectors in the strip had `Dio > mu`, then C2 would follow.

### 3.4 (d) Numerics (`.out` C2, H5, H6)

| word (`3x+1`) | `mu` | entropy (bits/letter) | `Dio` profile at scales `2^8 .. 2^14` | best exact gain | reconstruction at `10^4` bits |
|---|---|---|---|---|---|
| strip `3/4`, `W=2`, random | 1.189 | 0.33 | 1.086 → 1.003 | 18.7 | none of height `<= 2^4999` |
| strip `3/4`, `W=2`, greedy | 1.189 | 0.33 | 1.074 → 1.003 | 4.7 | none |
| strip `4/5`, `W=1.5`, random / greedy | 1.268 | 0.43 / 0.42 | 1.134 → 1.004 / 1.080 → 1.004 | 18.5 / 6.6 | none / candidate leaves at 10000 |
| strip `9/10`, `W=1.5`, random / greedy | 1.427 | 0.32 / 0.33 | 1.231 → 1.007 / 1.125 → 1.007 | 30.5 / 11.6 | candidate leaves at 10000 / none |
| strip `9/10`, `W=2`, random / greedy | 1.427 | 0.28 / 0.28 | 1.219 → 1.007 / 1.128 → 1.007 | 32.1 / 13.2 | none / candidate leaves at 10000 |
| strip `sqrt2/2`, `W=3/2`, random / perturbed | 1.121 | 0.43 / 0.25 | 1.094 → 1.003 / 1.184 → 1.012 | 15.4 / 21.8 | none |
| `Y` (width 1) | 1.427 | 0 | 1.707 → 1.097 → 1.025 at `2^18` | 77.6 (at `|UV| = 280`) | none; **irrational by Theorem Y** |
| `Y3` (width 1) | 1.427 | 0 | 2.23 → 1.48 → 1.20 at `2^18` | 1533.5 (at `|UV| = 7840`) | candidate leaves at 10002; OPEN |

The real Bernstein values, 30 digits each, and `Phi_2 mod 2^64` are printed per word. For example, the random `3/4` strip word has `Phi_R = -2.88147259529235767727736647323`, and `Y` has `Phi_R = -1.02811657612401143029094617068`.

## 4. What this means for the unified missing mechanism of foundry v4 (task 3)

Foundry v4 located the missing mechanism in the right place but described it on the wrong axis. The axis is the Diophantine exponent, not entropy. Every proved every-word argument built on periodic approximants reaches exactly the words with `Dio(w) > eta(w)`: Theorem R, S and D, the Subspace refinement, and any Liouville argument whose approximants have their natural heights (Proposition L). Below that line there are zero-entropy words: `Y` and `Y3` (width 1, `Dio = 1`), Bugeaud–Kim's Sturmian `s'` under `19x+1`, and linearly recurrent substitutive words under `5x+1`. So C2 is **not** the smallest instance; single explicit words are, and testing them is productive. The first one, `Y`, fell to a genuinely different mechanism. Its Bernstein number is a 2-adic theta value, and the Padé approximants of the q-difference equation `F(x) = 1 + rho x F(rho^2 x)` are rationals of height about `2^(0.95 lambda)`, far below the natural `2^(1.43 lambda)` of their eventually periodic parity words (Theorem Y). The mechanism family therefore has at least two axes: repetitions (`Dio`) and functional equations satisfied by the Bernstein number. HARD should be redefined as the supercritical words with `Dio <= eta` and no such functional equation. The cube-swap word `Y3` is the first explicit member: zero entropy, inside C2 at slope `9/10` and width 1. Its question is exactly whether `sum (2^10/3^9)^(k^3)` is irrational in `Q_2` (HC1), the 2-adic, non-unit-numerator analogue of the open cubic-theta problem. C2 remains the right smallest *class*: it is the supercritical analogue of Proposition B, the place where capacity and transversality both stop, and its generic words have `Dio = 1` and no functional equation. `Y3` is the right smallest *test case*: a method for C2 must in particular settle HC1.

## 5. Hypothesis candidates (for `INDEX.md` numbering by the integrator; no HYP files created)

* **HC1 (cube swap; OPEN).** `Phi_(3x+1)(Y3)` is irrational, where `Y3` has `b_n = 1^8 0 1` iff `n` is a nonzero cube and `1^9 0` otherwise. Equivalently, `sum_(k>=1) (2^10/3^9)^(k^3)` is irrational in `Q_2`.
  * This is the smallest open instance found: zero entropy, discrepancy width 1, `Dio = 1`, and no first-order q-difference equation.
  * FINITE-EXACT: no rational of height `<= 2^4999`.
  * Checked-but-unproved: the heuristic of section 1.5 predicts no Liouville certificate from periodic approximants.
* **HC2 (linearly recurrent barrier; OPEN).** Does every linearly recurrent binary word with uniform frequency `beta > log_3 2` satisfy `Dio(w) > beta log2 3`? If so, PC holds on all linearly recurrent words under `3x+r`.
  * Evidence: minimum margins `+0.295` over 9,716 substitutions and `+0.313` over 160 rotation codings.
  * Counter-evidence: the prefix-pattern constraint alone does not force it. A DFS finds binary words of length 60 with every prefix ratio `< log2 3` (scratch).
* **HC3 (decoupled Mahler, narrow strips; OPEN).** MZ(`9/10, 3/2`) holds (heuristic `P_MZ <= -0.099`); it implies C2(`9/10, 3/2`). The same for the width-1 strips (`P_MZ <= -0.48`).
* **HC4 (two places; sketch, not audited).** For strictly supercritical, bounded-discrepancy, non-EP `w` with `Dio(w) > 1`, `Phi_2(w)` and `Phi_R(w)` are not both algebraic. The route is Schlickewei's subspace theorem with the forms of section 1.5 at `2` and `inf`: the product is `2^(eta(|uv| - lambda))`. Proposition T is the elementary version: `Dio > 2` rules out "both rational".
* **HC5 (Sturmian beyond 2.50994; OPEN).** `Phi_(19x+1)` of the complement of Bugeaud–Kim's `s'` is irrational. Theorem D provably does not apply (only a gcd accident or a new mechanism could certify it), and PC is not expected for `19x+1` in general.
* **HC6 (q-holonomic swap sets).** If the generating function of `S` satisfies a first-order q-difference equation over `Q(x)` (squares, triangular numbers, values of quadratic polynomials), then the Lemma P construction settles the `S`-swap block words under `3x+r`. `S` = cubes is the first case outside.

## 6. Reproduction

```bash
cd <worktree>
bash 04-computation/experiments/collatz_procgen_20260922_hard_run.sh           # full run
# builds scratch/procgen_hard/hard_lpf from ..._hard_lpf.c and writes 05-knowledge/results/collatz_procgen_20260922_hard.out
# full run: about 4 minutes, one process, peak RSS about 170 MB, deterministic (fixed seeds)
bash 04-computation/experiments/collatz_procgen_20260922_hard_run.sh --quick   # about 1 minute, writes scratch/procgen_hard/hard_quick.out
```

**Sections of the `.out`:**
* H0: constants.
* H1: engine calibration.
* H2: Sturmian and quasi-Sturmian certificates.
* H3: Lemma J' check and the rotation census.
* H4: Lemma F, the substitution census and the drift side.
* H5: `Y`, the block identity, the Tschakaloff certificate, and the `5x+1` square-swap words.
* H6: `Y3`.
* C1: strips (entropy, `K_q`, overlap, `P_MZ`, forced cycles).
* C2: strip words.
* C3: Proposition M checks.
* C4: capacity and Monks–Yazinski measurements.

**Certificate methods:**
* `Phi_T(w) mod 2^N` is computed by divide and conquer and validated by parity iteration.
* Approximants `P/Q` are exact integers. Every certificate is checked through `v_2(Phi - P/Q) = lambda`, and `lambda` is recomputed directly from the word.
* Canonical (minimal) representations are used for the gcd statistics.
* Reconstruction follows each candidate until it leaves the word.
* The Tschakaloff forms are exact integers.

**Scratch exploration** (not needed for the `.out`): `scratch/procgen_hard/{explore_*,search_*,dfs_prefix,tschakaloff_check}.py`.

**Warning.** Dio estimates from finite windows are lower-bound-type, since the largest scale is capped by the window. The rotation and substitution censuses are FINITE-EXACT evidence: only certified words (a finite pattern with ratio `> mu`) count as PROVED.

## Sources

`[P]` means the primary text was read in this lane, `[R]` that a named secondary source was used, and UNVERIFIED that neither was possible.

* **Diophantine exponents and repetitions.**
  * Y. Bugeaud, D. H. Kim, *A new complexity function, repetitions in Sturmian words, and irrationality exponents of Sturmian numbers*, Trans. Amer. Math. Soc. 371 (2019) 3281–3308. [P, arXiv 1510.00279v3 read: Def. 3.2, Thms 3.3, 3.4 (with (3.1) and the example `s'`), Def. 3.5, Thm 3.6, Thms 4.2, 4.3, Defs 10.1, 10.2, Lemma 10.3.] The pages are from a search result.
  * B. Adamczewski, Y. Bugeaud, *Dynamics for beta-shifts and Diophantine approximation*, Ergodic Theory Dynam. Systems 27 (2007) 1695–1711. [P, author PDF read: section 2, Condition `(*)_rho` and `Dio`; Thm 1, threshold `log M(beta)/log|beta|` (the exact analogue of `eta` here); the proof of Thm 3, `p(n) <= cn` gives `Dio >= 1 + 1/c`.]
  * B. Adamczewski, Y. Bugeaud, *On the complexity of algebraic numbers I*, Ann. of Math. 165 (2007) 547–565. [P, arXiv math/0511674: Condition `(*)_w`, Thm 5 (Subspace criterion), Thm 6 (p-adic, Hensel expansions).]
  * B. Adamczewski, Y. Bugeaud, *Nombres réels de complexité sous-linéaire : mesures d'irrationalité et de transcendance*, J. reine angew. Math. 658 (2011) 65–98. [P, preprint: (2.2) `mu(xi) >= Dio`, Thm 2.1, section 5 (primitive fixed points have finite `Dio`), Prop. 11.1, `Dio <= Ind`.]
  * V. Berthé, C. Holton, L. Q. Zamboni, *Initial powers of Sturmian sequences*, Acta Arith. 122 (2006) 315–347. [P: section 2.2 `ice`, Prop. 2.1, Thm 1.1, Thm 1.2 and its proof, Cor. 3.5; ADQZ 2001 [3] for initial squares, [R].]
* **Structure of the classes.**
  * J. Cassaigne, *Sequences with grouped factors*, DLT'97, Aristotle Univ. Thessaloniki (1998) 211–222. [R: via D. Damanik, D. Lenz, arXiv math-ph/0105034, Prop. 2.2 (read), and Bugeaud–Kim section 3.]
  * F. Durand, *Linearly recurrent subshifts have a finite number of non-periodic subshift factors*, arXiv 0807.4430. [P: Def. 2, Prop. 5 (from Durand–Host–Skau), Props 6, 10.] The journal version (ETDS 20 (2000)) is UNVERIFIED.
  * F. Durand, B. Host, C. Skau, ETDS 19 (1999) 953–993. [R]
  * B. Adamczewski, *Balances for fixed points of primitive substitutions*, Theoret. Comput. Sci. 307 (2003) 47–75. [P: Thm 13, Cor. 15.]
* **q-series.**
  * W. Zudilin, *An elementary proof of the irrationality of Tschakaloff series*, arXiv math/0506086. [P: Theorem and the construction (1)–(7).]
  * L. Tschakaloff 1919 and P. Bundschuh. [R, via Zudilin.]
  * K. Väänänen, R. Wallisser, *A linear independence measure for certain p-adic numbers*, J. Number Theory 39 (1991) 225–236. [Crossref metadata and search snippet only; hypotheses UNVERIFIED.]
  * L. Leinonen, arXiv 1312.3819. [P, abstract and introduction; archimedean only.]
  * L. Ghidelli, *Arithmetic properties of cubic and biquadratic theta series*, arXiv 1910.05076. [P, abstract.]
  * Nesterenko's theorem and the transcendence of `sum q^(n^2)` (Duverney–Nishioka–Nishioka–Shiokawa). [R, search results only.]
* **Fractional parts.**
  * A. Dubickas, *Arithmetical properties of powers of algebraic numbers*, Bull. London Math. Soc. 38 (2006) 70–80. [P, abstract.]
  * Flatto–Lagarias–Pollington 1995. [R, via the parent note.]
  * Koksma's metric theorem. [R]
* **Collatz.** Parent notes and their sources: Bernstein 1994, Lagarias 1985, Bernstein–Lagarias 1996, Monks–Yazinski 2004, Tao 2022. [P in the parent lanes.]
* **Repository.**
  * [transversality foundry](collatz_procgen_20260922_transversality_foundry.md);
  * [in-house discrepancy](collatz_guards_20260921_discrepancy.md);
  * [foundry v4](collatz_procgen_20260922_foundry.md);
  * THM-2228 (Mahler map).
