# Discrete–continuous bridges for Collatz: natural boundaries, Mahler/Cobham, harmonic trees and Hex

**Status.** PROVED here, with hand proofs below and finite ingredients checked by the scripts:
* Theorem A: a periodicity lemma for every map `n/2, (qn+r)/2` with `q, r` odd.
* Corollaries B and C: the natural-boundary dichotomy. Collatz holds iff `B` is rational, iff algebraic, iff D-finite, iff `B` continues across one arc.
* The three controls: SHEET (`3n-1`), DRIFT (`5n+1`) and DEFECT (the planted map `T1`, unconditionally).
* Proposition D: Cobham-type rigidity fails for Mahler equations with roots of unity.
* Proposition E: the Terras kernel identity.
* Proposition F: SHEET for generating functions.
* Theorem G: every live backward tree has branching number at least `1.2334` for `3n±1`, `1.1227` for `5n+1` and `1.0870` for `7n+1`. So simple random walk on every such tree is transient.
* Proposition H: Althöfer's game has an infinite play from every odd `n >= 3`.

**KNOWN (priority).** The coordinator's candidate theorem (lead 1, parts (a)–(c)) is already in print for `3x±1`: Bell–Lagarias, Acta Arith. 170 (2015) 101–120, Theorems 1.1–1.3. Their proof uses the same commutator trick. Their Theorem 1.3 proved the Berg–Opfer conjecture that the three `3n-1` basin series have natural boundaries. Theorem A extends their Theorem 1.1 to every odd multiplier `q`, by a routine argument that does not need Skolem–Mahler–Lech.

**FINITE-EXACT.**
* Berg–Opfer's basin counts are reproduced exactly to `10^7`.
* The `3n-1` basin series satisfy no pure 2- or 3-Mahler equation of small order and degree.
* Kernel counts, tree growth `4/3`, branching-number profiles, harmonic-measure dimensions and the convergence of Tutte embeddings are all computed.

**OPEN.**
* Collatz itself.
* Whether any nonempty proper `T`-invariant set is 2-automatic, or 2-Mahler.
* Whether `T` is ergodic on the profinite integers (§6).

**No lead is refuted.** Lead 1 is not new, the coordinator's sixth-root equation is correct, and the Hex analogy fails at its first step.

Session `collatz-procgen-20260922` (mac-mini), lane "continuous bridges", 2026-09-25. No HYP or THM file was created.
* Scripts: `04-computation/experiments/procgen_continuous_20260925_{natural_boundary,mahler,tree_harmonic,hex_game,run}.py`.
* Output: [`procgen_continuous_20260925.out`](procgen_continuous_20260925.out).
  * Parts A–D; every check raises on failure.
  * Total about 35 s, one process at a time, peak resident memory 510–540 MB (part C).

Notation: `T_{q,b}(n) = n/2` for even `n` and `(qn+b)/2` for odd `n`. The plus sheet is `T = T_{3,1}` and the minus sheet `T_{3,-1}`. `B(z) = sum_{n in basin(1)} z^n` for `T`.

## 0. Ranking of the four angles by proof potential

| rank | angle | best statement | SHEET | DRIFT | DEFECT | proof potential |
|---|---|---|---|---|---|---|
| 1 | natural boundary: Collatz iff `B` continues across one arc of `|z|=1` | PROVED; KNOWN for `3x±1` (Bell–Lagarias 2015); general `q` here | **aware**: the `3n-1` basin series have natural boundaries | **aware**: `5n+1` | **aware**: the `T1` series has a natural boundary | Low. The criterion is an exact equivalence, so it carries the whole difficulty, and no known mechanism continues a 0/1 series across an arc short of periodicity |
| 2 | Mahler/Cobham: two pure Mahler equations imply `B` rational, hence Collatz | conditional implication PROVED; obstructions PROVED | the same mixed equation has natural-boundary 0/1 solutions on the negative side | `5n+1` likewise | the `T1` series satisfies one pure (2-Mahler) equation | Very low. Nothing supplies a pure 3-Mahler equation, and roots of unity provably destroy the rigidity |
| 3 | random walk, branching number, harmonic measure, Tutte embedding of the backward trees | transience and branching-number bounds PROVED | **blind** (PROVED by transport) | threshold-blind: transient for `5n+1` and `7n+1` too (PROVED) | blind | None for Collatz; a clean bridge and a negative result |
| 4 | Hex ⟺ Brouwer (Gale) versus Althöfer's `3n±1` game | infinite plays exist (PROVED); P-positions exist | n/a | n/a | n/a | None: the Hex mechanism has no analogue |

The common thread is Corollary 7 of the [mod-192 note](collatz_procgen_20260924_inverse_tree_mod192.md), seen from four sides. Every continuous object built from residues, Haar measure, graph structure or sizes is invariant under `x -> -x`. Only the natural-boundary criterion avoids this, because it is an equivalence and so names the positive cone. The price is that it inherits the full difficulty of the problem.

## 1. The natural-boundary equivalence (lead 1)

### 1.1 The periodicity lemma for all `qn+r` maps

Let `q >= 1` and `r` be odd with `q + r >= 2`, so that `T = T_{q,r}` maps `Z_{>0}` into itself. Let `r'` be the largest divisor of `r` that is prime to `q`. For example `r' = 1` when `r = ±1`, and `r' = 5` for `3n+5`.

Call `S ⊆ Z_{>0}` *eventually invariant* if there is an `n_1` such that `n ∈ S ⟺ T(n) ∈ S` for all `n >= n_1`. This covers basins, unions of components and, with `n_1 = m+1`, backward orbits `O^-(m)`.

**Theorem A (PROVED; checked in NB2).** An eventually invariant `S` is eventually periodic iff, for all large `n`, it coincides with `{n : n mod r' ∈ X}` for some `X ⊆ Z/r'` closed under `x -> 2x` and `x -> qx`.

In particular, if `r' = 1` an eventually periodic, eventually invariant set is finite or cofinite. An exactly invariant nonempty set that is eventually periodic is all of `Z_{>0}`.

*Proof.* Let `s = 1_S`, with `s(n+P) = s(n)` for `n >= n_0 >= n_1`.

1. *The period can be taken odd.* If `P = 2P'`, then for `n >= n_0`,
   `s(n) = s(2n) = s(2n+P) = s(2(n+P')) = s(n+P')`,
   using `T(2m) = m`. Iterating leaves an odd period `M`.
2. *Two relations on `Z/M`.* Write `s(n) = σ(n mod M)` for `n >= n_0`. Every class mod `M` contains large even and large odd numbers, because `M` is odd. Applying `T(2m) = m` to large `n = x (mod M)` gives `σ(2x) = σ(x)`. For odd `n = x (mod M)`, `T(n) = (qn+r)/2`, so `σ(x) = σ((qx+r)·2^-1)`. Combined with the first relation, `σ(x) = σ(qx+r)`.
3. *The classes of the relation.* Let `~` be the equivalence on `Z/M` generated by `x ~ 2x` and `x ~ qx + r`. Write `M = M_1 M_2`, where the primes of `M_1` divide `q` and `M_2` is prime to `q`.
   * Mod `M_1`, `q` is nilpotent. So `C(x) = qx + r` satisfies `C^J(x) = u (mod M_1)` for every `x` once `J` is large, where `u = r/(1-q)` is the unique fixed point of `C` mod `M_1` (`1-q` is a unit mod `M_1`).
   * Hence every `x` is `~` to a point of `F = {x = u mod M_1}`.
   * On `F`, `(u,y) ~ (2u,2y) ~ C^J(2u,2y) = (u, C_2^J(2y)) ~ (u,2y)`, and `(u,y) ~ (u, C_2 y)`. So `~` restricted to `F` contains the relations generated by `A_2 : y -> 2y` and `C_2 : y -> qy + r` on `Z/M_2`, and both are permutations there.
   * Their commutator is a translation: `A_2 C_2 A_2^-1 C_2^-1 (y) = y + r`.
   * So the classes on `F` are unions of cosets of `dZ/M_2`, with `d = gcd(r, M_2) = gcd(r', M)`. Modulo `d`, `A_2` acts as `×2` and `C_2` as `×q`.
   * Conversely, reduction mod `d` is equivariant (`qx + r = qx mod d`), so distinct `⟨2,q⟩`-orbits mod `d` are never joined.
   * So the `~`-classes are exactly the preimages of the `⟨2,q⟩`-orbits on `Z/d`.
4. *Conclusion.* `σ` is a function of `n mod d`, invariant under `×2` and `×q`, and `d | r'`. If `r' = 1`, `σ` is constant, so `S` is finite or cofinite. If `S` is exactly invariant and nonempty, it contains `2^k s_0` for every `k`. Then `σ = 1`, and every `n` lies in `S` because `2^k n` does for large `k`. ∎

**The case `3 | M` for Collatz itself.** The coordinator asked for this case; no CRT is needed.
* `6x+1 ~ 2x ~ x ~ 3x+1 ~ 6x+2`.
* `12x+1 ~ 4x ~ x ~ 3x+1 ~ 12x+4`.

When `3 | M` (with `M` odd), `6x` and `12x` run over `3Z/M`. So `v ~ v+1` and `v ~ v+3` for every `v = 1 mod 3`. Hence the class `1 + 3Z/M` is a single `~`-class, the class `2 + 3Z/M` joins it through `v ~ v+1`, and `x ~ 3x+1` sends `3Z/M` into it. When `3 ∤ M`, `v ~ v+1` holds for all `v` directly. For `3n-1`, replace `+1` by `-1` throughout.

**Check (NB2).** For all 90,180 triples (odd `M <= 1001`, 9 values of `q`, 20 values of `r`), the number of classes equals the number of `⟨2,q⟩`-orbits on `Z/gcd(M, r')`. Examples: `3n+5` has 2 classes mod 945, and `15n+7` has 3.

### 1.2 The dichotomy

**Corollary B (PROVED).** Suppose `r' = 1`. For an eventually invariant `S`, the series `F_S(z) = sum_{n in S} z^n` is either a polynomial (`S` finite) or `z/(1-z)` minus a polynomial (`S` cofinite). Otherwise `|z| = 1` is its natural boundary, and then `F_S` is transcendental over `C(z)` and not D-finite.

*Proof.*
* An infinite 0/1 series has radius 1.
* By the Pólya–Carlson theorem it is rational or has the unit circle as natural boundary. The theorem was conjectured by Pólya and proved by Carlson, Math. Z. 9 (1921) 1–13; the statement used here is the one read in Bell–Lagarias, Thm 2.2.
* A rational series with coefficients in `{0,1}` has eventually periodic coefficients: the recurrence state takes finitely many values.
* Theorem A then applies.
* Algebraic and D-finite functions continue along every path that avoids finitely many points, so neither can have a natural boundary. ∎

**Corollary C (Collatz, PROVED).** For `T = T_{3,1}` and `B = F_{basin(1)}` the following are equivalent:
* Collatz;
* `B` is rational;
* `B` is algebraic;
* `B` is D-finite;
* `B` continues analytically to a neighbourhood of a single point of `|z| = 1`;
* `(1-z)B(z)` is bounded on `|z| < 1`;
* `B(1/2) = sum_{n in basin(1)} 2^-n = 1`.

*Proof of the boundedness item:* `H^∞ ⊂ H^2` forces `sum |s(n) - s(n-1)|^2 < ∞`, so `s` is eventually constant. The last item holds because `B(1/2) <= 1`, with equality iff the basin is everything.

If Collatz fails, `B` and every backward-orbit series `f_{1,m}` has natural boundary (Bell–Lagarias Thm 1.2).

The coordinator's sixth-root form checks out exactly; see §2.1.

### 1.3 The three controls (all PROVED)

* **SHEET (`3n-1`).**
  * An independent census to `10^7` finds exactly the cycle minima `1, 5, 17`. The basin counts at `X = 10, 10^2, ..., 10^7` agree **exactly** with Berg–Opfer's Table 1.3; at `10^7` they are `3,273,791`, `3,244,985` and `3,481,224` (densities `0.32738`, `0.32450`, `0.34812`).
  * Each basin is nonempty and proper, so by Corollary B each series has natural boundary. This is Bell–Lagarias Thm 1.3, which proved Berg–Opfer's conjecture.
  * NB6 is a sanity check: no period `p <= 3000` fits on `[10^7 - 10^6, 10^7]`; the minimal mismatch fraction is `0.23`.
  * **Real-number form (NB7).** The basins are the binary digits of
    * `κ_1 = 0.9578097492942029409434...`,
    * `κ_5 = 0.0421819948142079932369...`,
    * `κ_17 = 0.0000082558915890658195...`.

    Each is irrational by Theorem A. Their sum is 1 iff the `3n-1` conjecture holds, and `B(1/2) = 1` iff Collatz holds.
* **DRIFT (`5n+1`).**
  * The cycles `{1,3,8,4,2}`, `{13,33,83,208,104,52,26}` and `{17,43,108,54,27,68,34}` are checked, and `r' = 1`.
  * So every nonempty proper invariant set, in particular `basin(1)`, has natural boundary.
  * Every backward-orbit series of `5n+1` is irrational, since no backward orbit is cofinite. This is the `5n+1` analogue of Bell–Lagarias Thm 1.3; it is not in their paper, which treats `3x+k` only.
* **DEFECT (planted `T1`).**
  * `T1(3·2^m) = 3·2^(m+1)` for `m >= 1`, and `T1 = T` otherwise.
  * `T^-1(3·2^m) = {3·2^(m+1)}`, so `T`-orbits from outside `S1 = {3·2^m : m >= 1}` never enter `S1`, and `basin_T1(1) = basin_T(1) \ S1`. The implication `T(n) ∈ S1 ⟹ n ∈ S1` is checked on `[1, 10^7]`.
  * The series is irrational **unconditionally**. The steps of Theorem A only need, in each residue class, representatives outside the density-zero set `S1`. The resulting `σ` is constant, which contradicts `2^k ∈ basin_T1(1)` together with `S1 ∩ basin_T1(1) = ∅`.
  * Hence it has natural boundary. Its coefficients below `2^71` are those of `z/(1-z) - sum_(m>=1) z^(3·2^m)` (Barina, CITED): a rational function minus a Hadamard-lacunary series. It equals that series iff Collatz holds.
  * So the criterion *sees* density-zero planting.

**General maps.** For `r' = 1`, every map `n/2, (qn+r)/2` with at least two components gives natural boundaries for all nonempty proper invariant sets. For `r' > 1` this can fail. Under `3n+5` the multiples of 5 form an invariant, periodic set, whose dynamics is `5×(3n+1)`. The correct statement is then Theorem A's: a class union mod `r'` closed under `×2` and `×q`.

### 1.4 Literature and priority (sources read unless marked)

* **Bell–Lagarias**, "3x+1 inverse orbit generating functions almost always have natural boundaries", Acta Arith. 170 (2015) 101–120, doi 10.4064/aa170-2-1, arXiv:1408.6884. The arXiv v1 was read in full.
  * *Theorem 1.1.* For `3x+k` with `k = ±1 mod 6`, a finite union of backward orbits has a rational series iff it is eventually a union of classes mod `|k|` closed under `×2` and `×3`.
  * *Theorem 1.2.* `f_{1,m}` has natural boundary for every `m` except possibly `m = 1, 2, 4, 8`, and Collatz is equivalent to their rationality.
  * *Theorem 1.3.* Every `f_{-1,m}` has natural boundary.
  * *Theorem 1.4.* For general `k`, all but finitely many `m`.
  * *Method.* Skolem–Mahler–Lech plus Claims 1–3. Claim 3 is the commutator `S1 S3 S1^-1 S3^-1 (r) = r + k`. Their Claim 2 handles `3 | d`; our nilpotent-CRT step does the same job for every prime dividing `q`.
  * **So lead 1 is a rediscovery**, and the repository should cite Bell–Lagarias for it. Our additions are modest:
    * Theorem A for every odd `q`, with the exact classification by `r'`;
    * the `5n+1` corollary;
    * the unconditional DEFECT control;
    * the independent reproduction of Berg–Opfer's counts.
* **Berg–Meinardus**, Results Math. 25 (1994) 1–12 (doi 10.1007/BF03323136) and Rostock. Math. Kolloq. 48 (1995) 11–18. Only the zbMATH reviews (Zbl 0810.11013, 0861.11008) were read. According to the reviews, the papers show that Collatz holds iff every solution of `h(z^3) = h(z^6) + (1/3z) sum_ν λ^ν h(λ^ν z^2)` analytic in the disc is of the form `h_0 + h_1 z/(1-z)`, and the 1995 paper adds that entire solutions are constant. [R]
* **Opfer**, "An analytic approach to the Collatz 3n+1 problem", Hamburger Beiträge zur Angewandten Mathematik 2011-09. The preprint was read, including its cover note.
  * The paper reformulates Berg–Meinardus as the kernel `K` of two operators `U, V` and claims `K = Δ_2`, hence Collatz.
  * **The known gap, in the author's own words.** A note dated 17 June 2011 on the cover says that the reasoning on p. 11 is incomplete: the claim that the vertices of the "annihilation graph" contain every even number `>= 6` exactly once. The note withdraws the claim of proof, "at least temporarily".
  * Our reading: that graph is generated from 8 by the backward algorithm (Lemma 4.13, Algorithm 4.14), i.e. it *is* the backward tree. The incomplete step is the Collatz conjecture itself.
* **Berg–Opfer**, "An analytic approach to the Collatz 3n+1 problem for negative start values", Comput. Methods Funct. Theory 13 (2013) 225–236 (doi 10.1007/s40315-013-0017-z). The Hamburg preprint 2012-11 was read; the published version was not, because Springer's page is behind a bot check, which was not bypassed.
  * They reformulate the `3n-1` problem: the conjecture holds iff exactly three linearly independent 0/1 solutions exist (their Thm 2.3). Their Table 1.3 is reproduced exactly above.
  * They *conjecture* the natural boundaries (preprint Remark 2.5; "Conjecture 2.4" of the published version, as cited by Bell–Lagarias) and propose Fabry's gap theorem as the tool.
  * Fabry needs a density-zero support, whereas the basins have empirical densities near `1/3`. Bell–Lagarias's route through Pólya–Carlson is the one that works.
* **Other work.**
  * Neklyudov, arXiv:2106.11859 (abstract read): the Berg–Meinardus operator has no nontrivial fixed points in `H^2(D)`. `B` is not in `H^2`, so this does not bear on `B`.
  * Efrem, arXiv:2510.06736 (abstract read): functional equations for generalized maps; nothing on natural boundaries.
  * Siegel, arXiv:2111.07882: withdrawn (zbMATH note read).

## 2. The Mahler/Cobham structure (lead 2)

### 2.1 The exact functional equations (PROVED; checked in NB1)

For `T_{q,b}` and any `T`-invariant `s`, with `h(z) = sum_{n>=1} s(n) z^n`:

```text
(BM)  h(z^q) - h(z^{2q}) = z^{-b} (1/q)  sum_{j<q}  ω_q^{-j(q+b)/2}  h(ω_q^j z^2)
(W)   h(w^q) - h(w^{2q}) = w^{-b} (1/2q) sum_{j<2q} ζ_{2q}^{-(q+b)j} h(ζ_{2q}^j w)
(D)   h(z) + h(-z) = 2 h(z^2)
```

Here `ω_q = e^{2πi/q}` and `ζ_{2q} = e^{2πi/(2q)}`.
* For `q = 3, b = 1`, (BM) is Berg–Meinardus's equation, since `ω^{-2j} = ω^j`.
* For `q = 3, b = -1`, it is Berg–Opfer's (2.6).
* For `q = 3, b = 1`, (W) is **the coordinator's form** `B(w^3) - B(w^6) = w^{-1}(1/6) sum_j ω^{-4j} B(ω^j w)` with `ω = e^{2πi/6}`, verified.

Comparing coefficients, (BM) and (W) are each equivalent to the pair of relations `s(2n) = s(n)` and `s(2m+1) = s(qm + (q+b)/2)`, that is, to invariance.

Checks:
* exactly, on random sequences that are constant on the components of the truncated graph (`q ∈ {3,5,7}`, both signs; residual 0 on up to 19,929 coefficients);
* on 50-point perturbations, which the identities detect;
* numerically, at random `|z| <= 0.85` on the three `3n-1` basins (residuals below `5e-16`).

### 2.2 What would suffice (PROVED implications from CITED theorems)

* **(S1)** `B` satisfies a nontrivial linear `k`-Mahler equation *and* one for `l`, with `k, l` multiplicatively independent (for example 2 and 3). Then `B ∈ C(z)`, by Adamczewski–Bell, *A problem about Mahler functions*, Ann. Sc. Norm. Super. Pisa (5) 17 (2017) 1301–1355 (abstract read; it proves that a power series satisfying both a `k`- and an `l`-Mahler equation is rational, as Loxton and van der Poorten conjectured), and Schäfke–Singer, JEMS 21 (2019) 2751–2792 (abstract read). Hence Collatz, by Corollary C.
* **(S2)** `basin(1)` is both 2-automatic and 3-automatic. Then it is eventually periodic by Cobham, Math. Systems Theory 3 (1969) 186–192. That paper was not read; the statement is the standard one, and Schäfke–Singer's abstract says their results give a new proof of Cobham's theorem. Hence Collatz, by Theorem A.
* **(S3)** `B` is D-finite. Then Collatz holds, by Corollary B.

Mahler-ness by itself buys nothing more. By Bell–Coons–Rowland (JIS 16 (2013) 13.2.10; abstract read), which reproves Bézivin, Nishioka and Randé, a D-finite Mahler function is rational, and a Mahler function meromorphic in the disc is rational or has the unit circle as natural boundary. That is the same dichotomy that Corollary B already gives.

### 2.3 Why the Collatz equation falls outside those theorems (all PROVED)

* **(O1) Roots of unity kill the rigidity (Proposition D).** Take the 3-smooth numbers `S = {2^a 3^b}` and `F = F_S`. Then
  * `F(z) + F(-z) = 2F(z^2)`, and
  * `F(z) + F(ωz) + F(ω^2 z) = 3F(z^3)`,

  which is exactly `s(2n) = s(n)` and `s(3n) = s(n)` (MB3). Yet `S` is infinite of density 0. Its counting function is about `(ln N)^2/(2 ln2 ln3)`: 142 points to `10^6`. So `F` is irrational, and by the Fabry–Faber gap theorem (statement read in Bell–Lagarias Thm 2.1) it has a natural boundary. More generally, every `S_A = {2^a 3^b m : m ∈ A}`, with `A` any set of integers prime to 6, satisfies both equations. Only countably many series are rational, so uncountably many of these solutions are irrational.

  So **no analogue of Adamczewski–Bell or Schäfke–Singer exists for Mahler equations with roots of unity**. The Collatz equation (D)+(BM) is of exactly that type: the sieve by parity is a sum over roots of unity.
* **(O2) SHEET, for generating functions (Proposition F).** Let `S ⊆ Z \ {0}` be `T_+`-invariant.
  * `F_+(z) = sum_{n in S, n>0} z^n` solves `E_{3,+1}`, and `F_-(z) = sum_{n in S, n<0} z^{-n}` solves `E_{3,-1}`, because `T_+(-m) = -T_-(m)`. This was checked directly for all `m <= 10^5` in MB4. Berg–Opfer note the related change `z -> 1/z`.
  * So *the 3n+1 equation on all of `Z` is the pair `(E_{3,+1} in z, E_{3,-1} in 1/z)`*, and the second has three natural-boundary 0/1 solutions.
  * Hence any argument that forces rationality using only features shared by the two equations proves a falsehood. Those features include the Mahler-with-roots-of-unity shape, the 2–3 multiplicative structure, 0/1 coefficients, growth, and Hadamard-type operations.
  * What differs is precisely the monomial `z^{-b}` and the sieved class `n = -b mod 3` in `h(z^3) - h(z^6) = z^{-b}(1/3) sum_j ω^{bj} h(ω^j z^2)`, i.e. the sign. This is the analytic face of Corollary 7 of the mod-192 note.
* **(O3) DRIFT.** `E_{5,1}` has natural-boundary 0/1 solutions (§1.3).
* **(O4) DEFECT: one pure equation is not enough.**
  * `B_T1(z) = z/(1-z) - G(z)` with `G(z) = z^6 + G(z^2)`. It is 2-Mahler and irrational, while `T1` has a divergent orbit (the support `N \ S1` is 2-automatic).
  * By Adamczewski–Bell it is not 3-Mahler; consistently, the MB1 search finds a 2-Mahler equation for it and no 3-Mahler one.
  * So the *second* pure equation in (S1) is indispensable, and it is the one no structure supplies.

### 2.4 What can be proved about the class of `B`

* `B` satisfies (D) and (BM). These are linear equations in the operators `σ_2`, `σ_3`, `σ_6` and the rotations by sixth roots of unity. On 0/1 solutions of such invariance equations a Pólya–Carlson dichotomy holds (Corollary B): rational or natural boundary. No sign-blind subclass can separate `3n+1` from `3n-1` (O2).
* **Proposition E (Terras kernel identity, PROVED).** For any `T`-invariant `s`, `s(2^e n + r) = s(3^a n + T^e(r))`, where `a` is the number of odd steps among the first `e` steps of `r`. This follows from Terras's `T^e(2^e n + r) = 3^a n + T^e(r)`.

  So the 2-kernel of an invariant set consists of progression subsequences `n -> s(3^a n + t)`, and `s` is 2-automatic iff only finitely many of them are distinct.

  **FINITE-EXACT (MB2).** For each of the three `3n-1` basins, the cumulative 2-kernel counts for `e <= 16` coincide **exactly** with the number of distinct pairs `(a, T^e(r))`:

  ```text
  1, 2, 4, 7, 12, 21, 38, 69, 127, 235, 438, 819, 1535, 2883, 5425, 10218, 19275
  ```

  So no two pairs give the same subsequence, and the statistic is identical for all three basins. It is a type-level quantity, equal to 1 for the Collatz basin.

  The 3-kernel counts are essentially maximal: `(3^{e+1}-1)/2` up to `e = 8`, and `88570`–`88573` at `e = 10` against the maximum 88573. Thue–Morse gives 2 as a control.
* **FINITE-EXACT (MB1).** No pure 2- or 3-Mahler equation with `(order, degree)` in `{(1,40), (2,30), (3,24), (4,19)}` holds for any `3n-1` basin series on 3000 coefficients. This was decided by exact rank mod `2^31-1`; full rank mod `p` implies no equation over `Q`, hence none over `C`. The solver *finds* the equations of all five controls: Thue–Morse (2-Mahler), the base-3 Cantor set (3-Mahler), an eventually periodic set, the all-ones truncation, and `B_T1` (2-Mahler only).

**Verdict on lead 2.** The route "prove both pure equations, then apply Adamczewski–Bell and Theorem A" is sound, but its first step has no source. The one equation Collatz does supply is mixed and uses roots of unity. That class is provably non-rigid (O1), sheet-symmetric (O2) and drift-symmetric (O3). One pure equation does not suffice even with DEFECT considered (O4).

## 3. Random walk, harmonic measure and Tutte's embedding on the backward trees (lead 3)

### 3.1 The branching number

Call a positive integer `x` *live* when its odd predecessor branch is ever legal:
* `x ≢ 0 (mod 3)` for `3n±1`;
* `x ≢ 0 (mod 5)` for `5n+1`;
* `x mod 7 ∈ {1, 2, 4}` for `7n+1`.

Non-live points have only the bare doubling ray `x, 2x, 4x, ...`.

**Theorem G (PROVED; checked in TC1).** For every live `x`, the backward tree of `x` has branching number (Lyons 1990)

```text
br >= λ_3 = u_3^{-1/2} = 1.2334422187,   2u^3 + u^2 = 1,              for 3n+1 and 3n-1;
br >= λ_5 = u_5^{-1/4} = 1.1227075951,   2u^5 + u^4 + u^3 + u^2 = 1,  for 5n+1;
br >= λ_7 = 1.0870353564  (numerical root of the analogous equation),  for 7n+1.
```

Consequences:
* Simple random walk on the tree is transient. So is simple random walk on every component of the undirected `T`-graph that contains a live point, by Rayleigh monotonicity.
* Bernoulli bond percolation with `p > 1/λ_3 = 0.8107` has an infinite cluster, since `p_c = 1/br`.

*Proof.*
1. **The ladders.** By Proposition 11 of the mod-192 note and its analogues, the odd predecessors of a live `x` sit on its doubling chain at positions `h0, h0+o, h0+2o, ...`, with `o = ord_q 2`: `o = 2, 4, 3` for `q = 3, 5, 7`. Each rung is at tree distance `h0 + oj + 1` from `x`.
2. **The rung residues.** The rungs satisfy `p_{j+1} = 4p_j ± 1`, `16p_j + 3` or `8p_j + 1`, so `p_j mod q` advances by a constant unit.
   * For `q = 3` exactly one rung in three is a multiple of 3, hence a bare ray.
   * For `q = 5` one rung in five is dead.
   * For `q = 7` three rungs in seven are live.
3. **The flow recursion.** The max-flow strength with capacities `λ^{-|e|}` obeys `g(v) = λ^{-1} sum_{children c} min(1, g(c))`, with `g = 1` on the truncation level `D`.
4. **The induction.** Take `φ = λ - 1` and argue from level `D` upwards. Write `K = D - |v|` and `Σ = sum_{live rungs} λ^{-(h0+oj+1)}`. If every live child has `g >= φ`, walking up the doubling chain of a live `v` (every intermediate bound stays `<= 1` because `φ <= λ - 1`) gives
   `g(v) >= φ (Σ - λ^{-K}/(λ-1)) + λ^{-K} >= φ + λ^{-K}(1 - φ/(λ-1)) >= φ`.
   Here the term `λ^{-K}` is the leaf at the end of the chain, and the last step needs the worst ladder type to satisfy `Σ >= 1`.
   * The worst type is `h0 = o-1` with the first live rung as late as possible: the first rung is dead for `q = 3, 5`, and the first three for `q = 7`.
   * At `λ = λ_q` its sum equals exactly 1: `(u^2+u^3)/(1-u^3) = 1` with `u = λ^{-2}` for `q = 3`. TC1 checks this over all types, and sharpness at `1.001 λ_q`.
5. **The limit.** So `g_D(root) >= φ > 0` uniformly in `D`. Every cutset contains a finite one, so `inf_Π sum_Π λ^{-|e|} > 0`, and `br >= λ`. Transience follows from Lyons 1990 (Lyons–Peres, *Probability on Trees and Networks*, Thm 3.5, read), and percolation from Thm 1.8 there. ∎

**The actual trees (TC1).**
* The uniform bound holds on the truncated trees: `min g_D(v; λ_q) = 0.533` against `φ = 0.233` for 3n±1, `0.449` against `0.123` for `5n+1`, and `0.454` against `0.087` for `7n+1`.
* The nodes whose ladder contains the root, where the tree cuts the root cycle, are excluded; they keep positive flow.

**Where the branching number actually sits (FINITE-EXACT, TC2–TC3).**
* The level sizes grow by `1.3333` per level for every `3n±1` root, which is `E[N_d] = (4/3)^d`, and by `1.20` (`5n+1`) and `1.14` (`7n+1`), i.e. `1 + 1/q`.
* The flow strength `g_D(root; λ)` is flat in `D` up to `λ = 1.32`. It decays for `λ >= 1.34`; at `λ = 1.40`, `D = 30 -> 46` gives `0.111 -> 0.051` for the root 1. This profile is consistent with `br = gr = 4/3`, and `6/5` for `5n+1`.
* Only the lower bounds in Theorem G are PROVED.

### 3.2 Harmonic measure and Tutte's barycentric embedding

Simple random walk with a reflecting root, stopped at depth `D`, exits with the unit current flow. At each vertex this flow splits in proportion to the effective conductances of the child branches, which gives an exact recursion.

**The Tutte bridge (PROVED, elementary).** Pin the depth-`D` leaves on the unit circle at `exp(2πi · 0.m_1 m_2 ... m_D)`, where `m_k = 1` for an odd move, so that subtrees occupy dyadic arcs. Place every interior vertex at the average of its neighbours; this is Tutte's barycentric condition. The result is `h(v) = E_v[position of the exit leaf]`, the harmonic-measure transform of the boundary map.

So convergence of the Tutte embedding as `D -> ∞` is transience made visible: the walk converges to an end, a point of a Cantor set of move-words in `Z_2`. It was computed by a two-pass `α/β` recursion (TC5). The root position settles: its increments per 4 levels are `1.1e-3`, `5.5e-3`, `8.8e-4` and `1.4e-3` on the four `3n±1` trees, decreasing in every case.

**Dimension drop (FINITE-EXACT, TC4).** The entropy dimension of harmonic measure in the 2-adic word metric:

| tree | harmonic dimension | `log2(growth)` | E-frequency on harmonic rays |
|---|---|---|---|
| `3n+1`, root 1 | 0.355 | 0.415 | 0.179 |
| `3n-1`, roots 1 / 5 / 17 | 0.341 / 0.358 / 0.361 | 0.415 | 0.165 / 0.190 / 0.205 |
| Haar model (multitype GW), two samples | 0.349 / 0.330 | 0.415 | 0.201 / 0.169 |
| `5n+1`, root 1 | 0.243 | 0.262 | 0.125 |
| `7n+1`, root 1 | 0.141 | 0.190 | 0.067 |

The Collatz trees behave like their Haar model: supercritical, transient, and with the Lyons–Pemantle–Peres dimension drop, as for Galton–Watson trees (LPP 1995, abstract read). The scatter between roots matches the scatter between two Haar samples.

### 3.3 Typing

* **SHEET-blind (PROVED; checked TC6).** `x -> -x` maps the `T_+`-tree of `x` onto the `T_-`-tree of `-x`; levels 0–30 agree as sets. Every invariant of rooted trees is shared. The sheets differ only in *which* roots are positive.
* **Drift-threshold-blind (PROVED).** Transience holds for `3n+1` (drift `log(√3/2) < 0`) and equally for `5n+1` and `7n+1` (positive drift). The values `br ≈ 1 + 1/q` measure the 3-adic (`q`-adic) legality rate, not the archimedean drift. The drift appears only in size-weighted data: the averaged tree series `1/(1-g(s))` of the mod-192 note (Prop. 10), whose residue is `1/drift`, and which is itself sign-blind.
* **DEFECT-blind.** Planting `T1` changes a tree only by one bare ray.
* **Cheeger constant 0 (PROVED).** Bare rays `3m, 6m, 12m, ...` have isoperimetric ratio `2/L`.
* **Planarity and the Tutte polynomial carry nothing.** The Collatz graph is a forest plus one cycle, whose Tutte polynomial is a monomial.
* **Schreier graphs.** The finite graphs `x ~ 2x, x ~ 3x+1` on `Z/M` are connected by Theorem A. `x -> -x` makes them isomorphic to the `3x-1` graphs, so any spectral-gap statement about them is sign-blind as well.

## 4. Hex ⟺ Brouwer and Althöfer's game (lead 4; under 15% of the effort)

Gale's theorem (Amer. Math. Monthly 86 (1979) 818–827; the paper was not read, and the statement is as commonly cited, UNVERIFIED here) says that Hex cannot end in a draw, and that this is equivalent to the 2-dimensional Brouwer theorem. It is a statement about *every completely filled board*.

**Proposition H (PROVED; checked HG1).** In Althöfer's game every odd `n >= 3` admits an infinite play. Of `3n+1` and `3n-1` exactly one is `2 mod 4`, and that move gives the odd number `(3n±1)/2 > n`. Repeating it never reaches 1; 10^4 such moves from 3 reach a 1762-digit number.

So "every complete play has a winner", the Hex/Brouwer mechanism, is **false** here. The prize question concerns optimal play only: no position has infinite remoteness, i.e. the least and greatest fixed points of the retrograde operator coincide.

* **Strategy stealing has no purchase.** The mover loses at P-positions, which are plentiful. The capped retrograde (HG2) gives 43 of the 100 odd `n < 200`: `7, 13, 23, 25, 29, ...`. The repository's density below `2^32` is 0.48 ([game note](collatz_procgen_20260922_althofer_game.md)).
* **Determinacy is available but does not help.** Gale–Stewart determinacy of open games (Ann. Math. Studies 28 (1953) 245–266; bibliographic data only) yields for each start either a forced win for the mover or a strategy for the opponent that avoids losing. It does not exclude draws.
* **Verdict: no mechanism.**

## 5. The bridges in one table

| discrete object | continuous object | sees | blind to |
|---|---|---|---|
| basin indicator | power series on the disc: rational or natural boundary (Cor. B) | SHEET, DRIFT, DEFECT (an equivalence) | nothing, and it supplies no mechanism |
| invariance relations | Mahler operators with roots of unity (BM), (W), (D) | nothing beyond invariance | rigidity (O1); SHEET via `z <-> 1/z` (O2) |
| residue classes | clopen invariant sets of `Ẑ_odd` (Theorem A): trivial | — | density-zero defects are not clopen |
| basin as digits | the real `κ = B(1/2)` (NB7) | Collatz iff `κ = 1` | — |
| backward tree | electrical network, random walk, harmonic measure, Tutte embedding | branching rate `1+1/q`, transience | SHEET (transport), drift threshold, DEFECT |
| sizes | Dirichlet series `1/(1-g(s))`, residue `1/drift` (mod-192 note) | DRIFT | SHEET |
| game | Hex/Brouwer | — | not applicable: infinite plays exist |

## 6. Remarks and open questions (no claims)

* **A measurable version of Theorem A (OPEN; HEURISTIC).** Theorem A says the clopen `T`-invariant sets of the profinite integers are trivial. The measurable analogue asks whether `T` is ergodic for Haar measure on `Ẑ`.
  * `T` is non-singular there.
  * On `Z_2 × Z_3` a contraction argument suggests yes. The `Z_3`-coordinate of `T^n(x)` forgets `x_3` once `a_n` is large, so invariant sets depend only on `x_2`, where `T` is an exact Bernoulli endomorphism.
  * The primes `p >= 5` need a compact-group-extension argument.
  * If true, every Besicovitch-limit-periodic invariant set would have density 0 or 1. Consistently, the residue-equidistributed `3n-1` basins, of density about `1/3`, would not be limit-periodic.
  * Such a statement would be DEFECT- and SHEET-blind, so it has no proof potential for Collatz. The literature was not checked.
* **2-automatic invariant sets (OPEN).** By Proposition E, a 2-automatic nonempty proper invariant set would need `s(3^e m - 1) = s(3^{e'} m - 1)` for some `e < e'` and all `m`; these are the kernel elements for `r = 2^e - 1`. Whether this is impossible is open. It would not give Collatz; see (O4).

## 7. Sources

**Read.**
* Bell–Lagarias arXiv:1408.6884v1, in full.
* Opfer, Hamburg preprint 2011-09, in full, with its cover note.
* Berg–Opfer, Hamburg preprint 2012-11: §§1–2 and the references.
* zbMATH reviews of Berg–Meinardus 1994 and 1995.
* Lyons–Peres, *Probability on Trees and Networks* (author PDF), Ch. 1 and Thm 3.5.
* Abstracts: Adamczewski–Bell (arXiv:1303.2019), Schäfke–Singer (arXiv:1605.02616), Bell–Coons–Rowland (arXiv:1210.2070), Lyons 1990 (OpenAlex), Lyons–Pemantle–Peres 1995 (OpenAlex), Neklyudov (arXiv:2106.11859), Efrem (arXiv:2510.06736).
* Crossref records for all DOIs quoted.

**Not read (bibliographic data only).**
* Carlson 1921, Pólya 1916 and Szegő 1922. The Pólya–Carlson and Fabry–Faber statements were read in Bell–Lagarias §2.
* Cobham 1969.
* Furstenberg 1967. Topological `×2×3` is not used; it would only enter §6.
* Gale 1979 (UNVERIFIED statement) and Gale–Stewart 1953.
* Tutte 1963, "How to draw a graph", Proc. LMS 13, 743–767. Only the barycentric condition is used.
* The published CMFT version of Berg–Opfer.

**Web policy.** All retrieval used the generic User-Agent `Mozilla/5.0 (research; math-repo)`, through the arXiv, Crossref, OpenAlex and zbMATH APIs, the Hamburg preprint server and Lyons's page. No bot check was bypassed and no personal data was sent.

## 8. Reproduction

```bash
cd <worktree>
python3 04-computation/experiments/procgen_continuous_20260925_run.py    # writes the .out; about 35 s
# or the parts separately (stdout = sections, stderr = memory):
python3 04-computation/experiments/procgen_continuous_20260925_natural_boundary.py   # NB1-NB7; 3n-1 census to 10^7
python3 04-computation/experiments/procgen_continuous_20260925_mahler.py             # MB1-MB4
python3 04-computation/experiments/procgen_continuous_20260925_tree_harmonic.py      # TC1-TC7
python3 04-computation/experiments/procgen_continuous_20260925_hex_game.py           # HG1-HG2
```

* Peak resident memory is 510–540 MB (part C), and each part runs alone.
* Every script ends with `ALL CHECKS PASSED`, and every check is an explicit `raise`.
