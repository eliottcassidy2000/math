# Seven and twenty-one on the Mersenne line, the openai/math collection read for compression, and a closed hypothesis: odd shift distances, Mersenne plateaus, debt resolution, the negacyclic clock, and Aut(T_k) = F_21

**Session:** mac-mini-2026-10-06-mod1819 (worktree `math-wt-chessboard-20261006`), 2026-10-06, continuing the mod-18/19 note.

**Owner's prompt (verbatim, with two uploaded PDFs: *The joint Dickman law for consecutive integers* and *Snaky in 21 Maker moves*):**
"think of how ideas ananlgous to {7,21} relate deeply with the idea that every odd-exponent Mersenne switch has an odd shift distance D, and 37 of the 60 odd exponents in [3,121] produce collision. also explore deeply and intricately this repository of work and look to blend its ideas with ours in useful ways and extend concepts meaningfully toward proofs https://github.com/openai/math thinking in terms of information compression creatively especially https://github.com/openai/math/tree/main/preprints/Integer-multiplication-below-n-log-n-September-23-2026 explore around there very open mindenly, pursuing any possible connection between themes or similar looking equations, even if the topics seem like they could not be related"

**Canon:**
* [THM-4556](../../01-canon/theorems/THM-4556-the-mersenne-line-is-a-chain-of-debt-states-odd-shift-distance-2-adic-periodicity.md): the Mersenne line is a chain of debt states, with odd shift distance, 2-adic periodicity and a certified density.
* [THM-4557](../../01-canon/theorems/THM-4557-doubling-a-doubly-regular-tournament-keeps-its-automorphism-group-hyp-9162.md): doubling keeps the automorphism group; **HYP-9162 is proved** (the theorem is known: Hanaki 2020, Thm 3.4; the proof here is new).

**Hypotheses:**
* [HYP-9213](../hypotheses/HYP-9213-mersenne-collatz-trajectories-coalesce-plateaus-of-odd-step-time.md): Mersenne coalescence.
* [HYP-9214](../hypotheses/HYP-9214-reset-two-debt-resolves-with-probability-tending-to-one.md): debt resolution tends to 1.

**Resolved:** [HYP-9162](../hypotheses/HYP-9162-sierpinski-tournament-tower-frobenius-21.md).

**Scripts:** [Mersenne plateaus](../../04-computation/experiments/sevens_20261006_mersenne_plateaus.py), [DRT doubling](../../04-computation/experiments/sevens_20261006_drt_doubling_aut.py) and [debt-resolution trend](../../04-computation/experiments/sevens_20261006_debt_resolution_trend.py). The first two `.out` files end in `ALL CHECKS PASSED`, and the third records the HYP-9214 table.

**Concurrent work on the same prompt.** The session opus-2026-10-06-S18 answered the same prompt in parallel: [mersenne_switch_parity_f21_compression_20261006.md](mersenne_switch_parity_f21_compression_20261006.md). Its four additions:
* **Parity law (Theorem 2).** For *every* source `n = 2^(r+1) t − 1`, the trailing-ones lags come in pairs `{D, D+1}`, and the first member has a parity fixed by `r` and `t mod 4`. So the least lag is odd for every reset-2 source. This generalizes THM-4556 (iii) beyond the Mersenne case `t = 1`.
* **The append maps.** `7 = A_1^3(0)` and `21 = R^3(0)`, where `A_1(x) = 2x + 1` and `R(x) = 4x + 1`, whose fixed points are this note's `−1` and `−1/3`. The two maps generate `F_21` exactly at 7.
* **Mirror clocks.** They are CRT-independent at every finite depth.
* **Certified density.** It reaches 0.1556 at template total 27. The median template total of an actual least-lag merge is 247.

We cite it rather than repeat it.

**Sources read.** The openai/math collection: 722 manuscripts, read from their LaTeX sources, read-only. The integer-multiplication paper was read in depth. Six parallel readers covered clusters of related papers: Hadamard/tournaments, Artin/prime predecessors, Dickman/Snaky, self-similar measures/entropy, information/computation, and Diophantine/DFT. Every connection below is typed.

## 0. Answers in brief

1. **Why the shift distance is odd (PROVED, THM-4556 (iii)).**
   * The reset switch glues each Mersenne pair `{2^(2k−1) − 1, 2^(2k) − 1}` at equal time: in ternary repunits, `U(R_(2k−1)) = oddpart(R_(2k))`.
   * So an odd-exponent number that meets an odd-exponent partner also meets the even one above it.
   * Hence the least shift is odd, and the nearest partner is `2^(2k) − 1 = 3·(4^k − 1)/3`: a backward leaf whose cofactor maps to 1, the `63 = 3·21` pattern.
2. **What "37 of 60" is (FINITE-EXACT, HYP-9213).**
   * The odd-step stopping time `σ(2^a − 1)` (OEIS A390816) is constant on long plateaus: 104 values for `a <= 6000`.
   * "37 of 60" counts odd exponents that are not first on their plateau. That share rises to 98.0% (2450 of 2500) for odd `a` in `[10^3, 6·10^3]`.
   * Switches are periodic in `a` on the 2-adic clock of 3. A finite class search certifies lower density `>= 15719/131072 = 0.1199` among odd exponents.
3. **Seven and twenty-one.**
   * `7 = 111_2` and `21 = 10101_2` approximate two distinguished 2-adic points. `−1` is the repelling fixed point of `(3x+1)/2` (Mersenne/run direction). `−1/3` is where `v_2(3x+1)` is infinite (trunk direction).
   * `F_21 = {x ↦ 2^i x + b mod 7}`, and 7 is the only Mersenne prime with `⟨2⟩ = QR_p`.
   * **Underneath all of it is Catalan's equation `2^3 + 1 = 3^2`** (Mihăilescu: the only consecutive powers). It makes `2^6 − 1 = 63 = 7·9` Zsigmondy's exceptional case: 7 sits in the cyclic half `2^3 − 1`, and the negacyclic half `2^3 + 1 = 9` is a pure power of 3.
   * The same identity, as `9 = 8 + 1`, generates 73% of the degenerate word collisions in the census.
4. **The multiplication paper.**
   * Its key ring `Z[y]/(y^r + 1)` has a 3-adic twin: in `Z/3^(n+1)`, `2^(3^n) = −1`. Our clock tower splits by CRT into a cyclic half `2^(3^n) − 1` (primes 7, 73, 262657, …) and a negacyclic half `2^(3^n) + 1` (`3^(n+1)` times 19, 87211, 163, …). In particular `19·27 = 2^9 + 1`.
   * Its power saving is a subcritical recursion, `s/W < m` subcalls per scale. That is a similarity dimension below 1, like our survivor dimension `h`.
   * Its rank saving rests on Alon's orthogonal-representation construction.
5. **A closed hypothesis (THM-4557).** Reading the Hadamard papers, the session proved HYP-9162: `Aut(T_k) = F_21` for all `k >= 3`. More generally, doubling any doubly regular tournament on at least 7 vertices keeps its automorphism group.
   * The independent audit found the theorem in the literature: Hanaki 2020, Theorem 3.4, with a different proof.
   * So HYP-9162 was already known when posed. Ours is a new proof (MISTAKE-571).
6. **Toward proofs.**
   * The reset-2 "debt" now has a precise probability, and it **grows with size**: 0.19 at 16 bits, 0.82 at 2048, 0.91 at 8192 (HYP-9214).
   * A logarithmic-density proof of HYP-9214 would make the rewrite program's obstruction a density-zero set.
   * The joint Dickman law gives, conditionally, that `P^+(n)` and `P^+(U(n))` are independent over odd `n`.
   * Nothing in the collection supplies an input for the thin sets our program needs.
## 1. Seven and twenty-one on the Mersenne line (PROVED unless typed otherwise)

**Two distinguished 2-adic points.** `7 = 111_2` and `21 = 10101_2` approximate two special 2-adic points of the Collatz dynamics (`7 ≡ −1 mod 8`, `21 ≡ −1/3 mod 64`):

* `−1 = ...1111_2` is the fixed point of `f_1(x) = (3x+1)/2`. It is the limit of the Mersenne numbers `M_a = 2^a − 1`, and it is the point where THM-4555's switches live: the run `1^r` means "`n` is `r` digits close to `−1`".
* `−1/3 = ...010101_2` is the root of `3x + 1 = 0`, where `v_2(3x+1)` is infinite. It is the limit of the trunk `a_k = (4^k − 1)/3 = 1, 5, 21, 85, …`, the direct predecessors of 1 (`3a_k + 1 = 4^k`).

The two points are linked by `3·(−1/3) = −1`. At finite level this reads `2^(2k) − 1 = 3·a_k`, so `63 = 3·21`. Among Mersenne numbers, the even-exponent ones are exactly those divisible by 3, i.e. those whose third approximates `−1/3`.

**Ternary repunits.** Put `R_a = (3^a − 1)/2 = 3R_(a−1) + 1`. This is the orbit of 0 under the *undivided* map `x -> 3x + 1`, and `R_a = 111…1_3`. Then:

* `U^a(M_a) = oddpart(R_a)` for `a >= 2`: the binary repunit passes through the ternary repunit.
* Example: `7 -> 11 -> 17 -> 13 = 111_3`.
* Example: `63 -> … -> 91 = 7·13`.

**The reset pairing.** For even `a = 2k >= 4`, `R_(a−1)` is odd and `U(R_(a−1)) = oddpart(R_a)`. Hence `U^a(M_a) = U^a(M_(a−1))`: the pair `{M_(2k−1), M_(2k)}` merges at equal odd-step time. This is the reset switch (Ahmed 2016; `checked_switch_phase19` (3)). Its standard-step shadow is OEIS A193688's `a(2n) = 1 + a(2n−1)` (Herrera 2023).

For odd `a`, `R_a = 3·2^v·oddpart(R_(a−1)) + 1` with `v = 1 + v_2(a−1)`. This is precisely a debt state `(e, h) = (1, v)` of `checked_switch_phase19` (4). **The Mersenne line is a chain of debt states.**

**Odd shift distance (answers the owner's question).**

*Lemma.* Let `a` be odd, `b < a` odd, and suppose `M_a` and `M_b` merge at equal odd-step time. Then `M_a` also merges with `M_(b+1)`.

*Proof.* `M_(b+1)` and `M_b` merge at time `b + 1` by the pairing, and `b + 1 < a <=` the merge time. ∎

Hence the **least** shift `D = a − b` is always odd. The S18 note's parity law (Theorem 2) extends this to every source `2^(r+1) t − 1`: lags come in pairs `{D, D+1}`. The nearest partner `M_b`, with `b` even, is then a multiple of 3: `M_b = 3·a_(b/2)`, a leaf of the backward tree whose cofactor reaches 1 in one step.

So the least-shift switch from an odd exponent lands on a "63 = 3·21" number. The whole partner set is a union of reset pairs, so other switches need not: `M_31` merges with `M_30` and `M_29`, and later with `M_26` and `M_25`.

**Plateaus (FINITE-EXACT; A390816).** Let `σ` be the number of odd steps to 1.

* Equal-time merging of `M_a` and `M_b` is the same as `σ(M_a) = σ(M_b)` together with a first meeting before 1.
* `σ(M_a)` takes only **58 values for `a <= 1200` and 104 for `a <= 6000`**. The count of distinct values grows roughly like `A^0.37`.
* Its level sets are unions of reset pairs, with up to 59 pairs in one set.
* The "37 of 60" is the fraction of odd `a` that are not the first of their level set. That fraction rises to 2450/2500 = 0.980 for odd `a` in `[10^3, 6·10^3]`. Of these, 0.882 share `σ` with `a − 1`.
* For `a <= 6000`, every coincidence of `σ` is an equal-time merge (FINITE-EXACT; the converse direction is not proved).
* Near `a = 6000` the trajectories run in two interleaved streams (`σ = 27728` and `28762`).

**2-adic periodicity (PROVED) and a certified density (FINITE-EXACT).**

* The reduced 2-adic word of `M_a` up to total `K` is the word of `2·3^(a−1) − 1 mod 2^(K+1)`. So it depends only on `a mod 2^(K−2)`, the 2-adic clock of 3.
* Every collision switch therefore propagates along a residue class. For example, `(a0, D, K) = (95, 1, 9)` gives switches for all `a ≡ 95 (mod 128)`, each merging at step `(a − 1) + 3`.
* A rigorous class search certifies a lower density of **`15719/131072 = 0.1199`** of switching exponents among odd `a` (`K = 20`, classes mod `2^18`). Allowing even `D` certifies the same classes. The S18 note pushes the certified share to 0.1556 at `K = 27`.
* The mirror: backward-minimality of `M_a` runs on the 3-adic clock of 2 (`a mod 2·3^(k−1)`, THM-4554(v)). The forward switches run on the 2-adic clock of 3.

**Debt resolution grows with size (NUMERICAL, seeded).**

* A random first-reset-2 source of `B` bits merges at equal odd-step time with `(n−1)/2` with frequency 0.19, 0.33, 0.41, 0.50, 0.61, 0.70, 0.79, 0.82 for `B = 16, 32, …, 2048`. The share of merges within 10 steps stays near 12%, so the growth comes from late coalescence.
* Even debt heights resolve more often than odd ones (0.68 vs 0.57 at `B = 256`). This is consistent with the 239 seeds, 92% of which have odd debt height.
* THM-4555's 15% for `n < 2·10^4` is a small-number effect.

**The group behind `{7, 21}` (PROVED identification; ANALOGY in use; nothing in THM-4556 uses it).**

* `F_21 = {x -> 2^i x + b (mod 7)}`: doubling and translation modulo `M_3`. Its analogue on `Z/M_a` has order `a(2^a − 1)`; doubling mod `M_a` is the cyclic rotation of `a`-bit words.
* On the Sierpiński tower, `F_21` acts on the low three bits with orbits of size 7, and the multiples of 8 induce a smaller tower: `63 = 8·7 + 7`.
* On the Mersenne line, the reset involution `2k − 1 <-> 2k` acts on the low bit of the exponent with orbits of size 2, and the level sets of `σ` are unions of these orbits.
* In both, a low-bit symmetry's orbits are the units of a coarse self-similar structure. The odd-`D` lemma says every switch lands on the top of an orbit, as the HYP-9162 reduction lands on the unique source of the invariant set `D`.
## 2. Integer multiplication below n log n, read against the Collatz program

**What the paper does.** OpenAI, *Integer multiplication below n log n*, 2026-09-23: [openai/math, preprints/Integer-multiplication-below-n-log-n-September-23-2026](https://github.com/openai/math/tree/main/preprints/Integer-multiplication-below-n-log-n-September-23-2026). We read it from its LaTeX sources.

* **The bound.** It multiplies two `n`-bit integers in `O(n (lg n)^(1−κ))` time, with `κ = 2^−182`, on a fixed multitape Turing machine. This disproves the `n log n` optimality conjecture in that model.
* **Why `n log n` is the barrier.** In the Harvey–van der Hoeven pipeline (Gaussian resampling, Bluestein's chirp, then Nussbaumer's synthetic transforms over `Z[y]/(y^r + 1)`), the cost is dominated by scanning `Θ(n)` bits at each of `Θ(log n)` transform levels.
* **The fix: two tape procedures.**
  * An interchange of two address fields of `u` bits costs `O(V u^τ)`. Its intermediate data are XOR combinations.
  * A simultaneous Hadamard butterfly `H_0(u, v) = ((u+v)/2, (u−v)/2)` on `d` selected address bits costs `O(V d^λ')`. Its intermediate data are exact Gaussian dyadics.
* **Where the speed-up comes from.** Both procedures are recursive linear networks.
  * A network on `W` wires makes `s < Wm` smaller calls when the parameter is multiplied by `m`.
  * The recursion is `F(e) <= (s/W) F(e/m) + O((eK)^τ + 1)`, with `τ = 1 − 2^−50`.
  * The rank saving comes from triples of a 100-set. Over `F_2`, the identity is "gather then scatter" (incidence matrix `M`, rank `<= 100`) plus a sparse correction on the pairs meeting in one point: `I = M M^T + A_1`.
  * Every scratch wire is restored ("transparent computation": Buhrman–Cleve–Koucký–Loff–Speelman).
  * The subspace labels come from Alon's 1998 orthogonal-representation construction, with the form `x^T (I − J/9) y` and `⟨t_S, t_T⟩ = |S ∩ T| − 1`.

**What it shares with our objects.**

1. **The negacyclic ring and the 3-adic clock (PROVED identities; ANALOGY in use).**
   * The paper's key ring is `Z[y]/(y^r + 1)`. There `y^r = −1`, so multiplying by `y^j` is a signed shift.
   * In `Z/3^(n+1)`, the number 2 plays the role of `y`: `2^(3^n) ≡ −1` and `2` has order `2·3^n`. In the clock's "log frame" `a -> 2^a`, multiplication by 2 is a signed cyclic shift.
   * The tower splits by CRT exactly as `x^(2r) − 1 = (x^r − 1)(x^r + 1)` does: `2^(2·3^n) − 1 = (2^(3^n) − 1)(2^(3^n) + 1)`.
     * The cyclic half holds the primes whose 2-order divides `3^n` (7, 73, 262657, 2593, 71119, 97685839).
     * The negacyclic half is `3^(n+1)` times the primes with 2-order exactly `2·3^m` (19, 87211, 163, 135433, 272010961).
   * In particular **`19·27 = 513 = 2^9 + 1`**. This is the whole reason the mod-19 and mod-27 clocks coincide.
   * `63 = 7·9` straddles the halves at level 2: 7 is cyclic, and 9 is the 3-adic negacyclic part. `21 = 63/3` takes one factor of 3 from the negacyclic side.
2. **Index maps.** The paper places coefficients on distinct prime cyclic axes by CRT, and Bluestein turns a DFT into a convolution. Rader's algorithm turns a prime-length DFT into a cyclic convolution of length `p − 1`, indexed by powers of a primitive root.
   * Our `(Z/27)^* ≅ F_19^*` via `2 -> 2` is the coincidence of the Rader index maps of 19 and of the units mod 27.
   * Every prime `2·3^k + 1` (7, 19, 163, 487, 1459, …) has a 3-smooth Rader length `2·3^k`. Requiring 2 to be the primitive root is a choice of generator, under which Rader's index map is the Collatz clock of level `k + 1`.
   * That choice holds iff `k` is even and 2 is not a cube mod `p`. It happens for `k ∈ {0, 2, 4, 16, 320, 782, 1252, 1454, 5480}` (`k <= 9000`; FINITE-EXACT, Lucas certificates), correcting the list `3, 19, 163, 86093443`, which is complete only below 320.
   * The identification of `(Z/27)^*` with `F_19^*` is at the index level only: the kernels `ζ_19^(2^j)` and `ζ_27^(2^j)` have different spectra.
   * The collection's exact-DFT paper uses no primitive roots, so relative to it this is ANALOGY.
3. **Subcritical recursion = similarity dimension (ANALOGY).**
   * The paper's cost exponent `log_m(s/W) < τ < 1` solves a Moran equation `(s/W) m^(−τ) = 1`.
   * Our forward survivors have the analogous saving `1 − h = 0.050` (`W_k ~ 2^(hk)`): describing a non-descending parity word costs `h < 1` bits per step.
   * The backward sieve has no saving: its survivors have positive measure, `s_∞ ∈ [0.288, 0.299]`.
4. **Dyadic coefficients (ANALOGY).**
   * The paper's intermediate values live in `Z[i, 1/2]`. Our collision values `f_u(−1) = N/2^(A−1)` live in `Z[1/2]`.
   * The butterfly `H_0` is our Sierpiński generator `H_2 = [[1,1],[−1,1]]` up to sign and scale (Sylvester versus skew).
5. **Changing frames.** The paper's central device stores `D_v g` instead of `g`, choosing frames `D_v` so that each edge is cheap.
   * Collatz has no single frame: `×3` is a shift in base 3, `÷2` is a shift in base 2, and `+1` is a carry.
   * In the double-base number system (sums of `2^i 3^j`) both `×3` and `÷2` are free exponent shifts, and all the difficulty sits in the `+1`.
   * Our inverse numerators `B_w = Σ 3^(j−1−i) 2^(S_i)` are {2,3}-representations whose terms form an antichain under divisibility, of Erdős–Lewin / Blecksmith–McCallum–Selfridge type (A355915). They are *not* Dimitrov–Imbert–Mishra double-base chains, whose terms are totally ordered by divisibility.
   * For fixed length the representation is unique (2-adic peeling). THM-4555's collisions are coincidences of two such representations of different lengths, i.e. vanishing sums of signed {2,3}-units (section 5d).

**NUMEROLOGY (recorded, not used).**
* The form `I − J/9` has eigenvalue `1 − 100/9 = −91/9`, and `91 = 7·13 = U^6(63) = oddpart(3^6 − 1)`.
* `κ = 2^−182`, with `182 = 2·91`. In the paper, `κ` comes out of a chain of parameter choices starting at `τ = 1 − 2^−50`. No mechanism links it to Collatz.
## 3. A standing hypothesis closes: Aut(T_k) = F_21 for all k (THM-4557; known since Hanaki 2020, new proof)

Reading the collection's Hadamard papers (*The circulant Hadamard conjecture*; *Exact Fourier certificates for complex Hadamard matrices of order six*; *The maximum number of mutually unbiased bases in dimension six*; *A proof of Seymour's second neighbourhood conjecture*), the session's reader proved the second-copy lemma to which section 5 of the mod-18 note had reduced HYP-9162. The papers supplied the mood (kernels at maximum support rank; characteristic-2 contradictions via the all-ones `J`), not a technique.

* **Statement.** For every doubly regular tournament `T` on `4t + 3 >= 7` vertices, the skew-Hadamard doubling `D(T)` has `Aut(D(T)) = Aut(T)`. Hence `Aut(T_k) = F_21` for all `k >= 3`.
* **Prior art (found by the audit; MISTAKE-571).** This is Theorem 3.4 of A. Hanaki, arXiv:2011.06141 (2020), proved there with triple intersection numbers. The intransitivity goes back to Faradžev–Klin–Muzichuk (1994). So HYP-9162 was already a corollary of the literature. What is new here is the proof.
* **Proof idea.**
  * The out-neighbourhood of a second-copy vertex `x'` is `T` with the arcs inside `N^−(x)` reversed.
  * Double regularity would force the odd-order skew matrix of `T[N^−(x)]` to kill the indicator of `N^+(j) ∩ N^−(x)`.
  * Mod 2 that matrix is `J − I`, of rank `m − 1`, so its kernel is the constants. But the indicator has `t + 1` of `2t + 1` entries.
* **Checks.** We verified it ourselves (FINITE-EXACT, nauty): `T_3..T_8`, Paley `P_7..P_43`, and double doublings.
* **Where it fails.** At `t = 0` the argument fails, and there `|Aut|` jumps from 3 (`T_2`) to 21 (`T_3 = P_7`).

**Side results from that reading (PROVED or FINITE-EXACT as marked; not used elsewhere).**
* THM-4552(iii)'s transversality of knight directions is exactly mutual unbiasedness of Weyl–Heisenberg eigenbases: `det(v, w) ∈ Z_n^×` (PROVED, standard criterion). The four knight directions are pairwise MU iff `gcd(n, 30) = 1`.
* Barker-7 is the Legendre sequence of `P_7`. Its minus-set is a `(7,3,1)` Fano difference set whose multiplier group gives `F_21 = Z_7 ⋊ ⟨2⟩` (classical).
* Among Mersenne primes, `⟨2⟩ = QR_p` only for `p = 7` (PROVED: `k = 2^(k−1) − 1` forces `k = 3`). **This is a second, number-theoretic reason that 7 and 21 are special.** `F_21` is the affine group of `⟨2⟩` acting on `F_7`, and `⟨2⟩` is exactly the set of squares only at `p = 7`.
* Tao's `S_6` (order-6 complex Hadamard) realises the outer automorphism of `S_6` on row versus column permutations (FINITE-EXACT; probably known).

## 4. Artin primitive roots and prime predecessors

The papers read: *Primitive roots for every admissible integer base*; *Simultaneous primitive roots*; *The Poisson–Dirichlet law for prime predecessors*; *Prime predecessors with an even number of prime factors*; *Weighted dilation graphs, smooth shifted primes and totient fibers*.

All of them are unrefereed preprints. The primitive-root theorem rests on the collection's quasi-Riemann-hypothesis zero-free strip. It gives `>> x/(log x)^2` primes in `(x, 2x)` with `a` primitive.

* **Correction to our list (FINITE-EXACT, Lucas certificates).** The primes `2·3^k + 1` with 2 primitive are `k ∈ {0, 2, 4, 16, 320, 782, 1252, 1454, 5480}` for `k <= 9000`. Odd `k` never works, because then `p ≡ 7 (mod 8)`. The earlier list `3, 19, 163, 86093443` is complete only for `k < 320`.
* **Zsigmondy (PROVED, classical).**
  * The primes of 2-order exactly `2·3^k` are the primitive divisors of `Φ_(2·3^k)(2)`. The tower's level 2, `2^6 − 1 = 63`, is exactly Zsigmondy's exceptional case for base 2: no prime has order 6.
  * So the level-2 "new prime" 7 lives in the cyclic half (order 3), and the negacyclic half `2^3 + 1 = 9` is a pure 3-power.
  * This is the {7, 21} level of section 2's splitting: 7 is cyclic, and `21 = 7·3` borrows its 3 from the negacyclic side.
* **The right form of the Rader coincidence (PROVED, Chebotarev).** `(Z/3^(k+1))^*` is a 2-equivariant quotient of `F_p^*` iff `2·3^k | ord_p(2)`. These primes have density `(17/24)·3^(2−k)/8`, and the observed values match.
* **Nothing transfers.** The papers' primes have `p − 1` with a prime factor above `x^0.9`, the opposite of 3-smooth. Our objects (`2·3^k + 1`, `2^K − 3^d`, `Φ_(3^k)(2)`) are exponentially thin.
* **Cycle denominators (PROVED, elementary).** A fixed prime `p ∤ 6` divides `2^K − 3^d` with frequency `[F_p^* : ⟨2, 3⟩]/(p − 1)` over full periods. So primes at which 2 or 3 is primitive are the *least* likely divisors of cycle denominators.
* **NUMEROLOGY.**
  * Lenstra's GRH density `(4/5)A = 0.299165` of 2-primitive primes in `p ≡ 1 (mod 3^k)` sits `4.5·10^-5` above our upper bound for `s_∞`.
  * Hasse's `7/24 = 0.291667` (primes with `ord_p(2)` odd) sits inside `[0.28820, 0.29912]`.
  * So do a dozen other constants. There is no mechanism behind any of these.
## 5. The two uploaded papers: the joint Dickman law and Snaky in 21 Maker moves

**The joint Dickman law** (OpenAI, 2026-09-24).
* `log P^+(n)/log n` and `log P^+(n+1)/log n` are asymptotically independent Dickman in natural density. This resolves the Erdős–Pomerance conjecture, and `P^+(n) < P^+(n+1)` has density 1/2.
* **Method.** Bin characters `f_x`, `g_x` are invariant under fixed multipliers, and Fourier inversion reduces the law to a mixed decorrelation. A nonzero correlation produces a profile on `(0, ∞) × Ẑ`. An amplifier counts factorizations `n = am`, `n + 1 = cl` with rough squarefree `a, c`. Cauchy–Schwarz then yields a positive energy, which a lag-graph substitution and short-interval averages (Matomäki–Radziwiłł(–Tao)) destroy.

* **A Collatz corollary (CONDITIONAL on the paper; derivation by the session's reader; not independently audited).**
  * The paper states no progression version.
  * Its proof gives the joint law on every divisibility class `{gcd(n, q) = d}`. The amplifier's coefficients are rough, so `q | am` iff `q | m`, and the restricted energy is dominated by the paper's unrestricted one.
  * For odd `n`, the bin labels see `U(n)` exactly as `3n + 1`, because `f_x(2) = 1`.
  * Hence **over odd `n`, `P^+(n)` and `P^+(U(n))` are asymptotically independent Dickman, and `P^+(U(n)) > P^+(n)` has density 1/2.**
  * Numerically, to `10^7`, the joint-to-product ratio is 0.950 (it was 0.923 at `10^6`) and the ordering density is 0.5151. The bias fits `1/2 + 0.244/log X`.
  * Unit classes and the multi-step version need new work.
* **The trailing-ones switches seen multiplicatively (PROVED, one line).**
  * `m + 1 = (n+1)/2^D`, so a switch preserves `oddpart(n + 1) = t` and `P^+(n + 1)`.
  * THM-4555's switch families are exactly the fibres of `n -> oddpart(n + 1)`.
  * The hostile Mersenne family is the fibre `t = 1`, the only one in which `n + 1` has no odd prime.
  * Along the ones-run, `T(n) + 1 = (3/2)(n + 1)`, so `max(3, P^+(t))` is conserved.
  * Beyond this, the labels are blind to 2 and 3, and Collatz merging is decided by 2 and 3 alone. There is no further link.
* **Reference correction.** arXiv:1511.09141 (on consecutive integers with equal Collatz height) is Elia–Tucker. Its counterexamples are to conditions conjectured by Garner.

**Snaky in 21 Maker moves** (OpenAI, 2026-09-25).
* **Result.** Maker completes the Snaky hexomino within 21 moves on the empty infinite board. The certificate is 728 numbered conditions `(A, T, h)` closed under a combination rule: at a pivot `p`, `T = {p} ∪ ⋃T_i`, `A = (⋃A_i ∪ ⋂T_i) \ {p}`, `h = 1 + max h_i`. Its final card has `A = ∅`, a 251-cell envelope, and height 21.
* **Our re-evaluation (FINITE-EXACT, the reader's own parser).** 4089 references and 1620 combination nodes (898 inline). Cards per height 1..21: 6, 42, 66, 94, 100, 86, 65, 51, 56, 44, 32, 20, 26, 12, 12, 5, 3, 3, 1, 3, 1.
* **Comparison with a Collatz certificate (ANALOGY, with the disanalogies spelled out).**
  * Snaky's certificate is finite modulo an infinite symmetry group (translations and `D_4`). Collatz has no symmetry group. THM-4555's uniform switches play that role, but only along collision families.
  * Snaky's height is bounded (21). Collatz stopping times are unbounded, so a Collatz certificate needs induction.
  * "Breaker must miss one child envelope" is a covering lemma, made finite because one Maker claim owns `⋂T_i`. Collatz covers are pure "for all": the 2-adic point `−1` lies in every all-ones cylinder and can never be pre-empted.
  * Each uniform family accumulates at `−1` for one value of `t`. So finitely many families cannot cover a neighbourhood of `−1`.
* **NUMEROLOGY, dissolved.**
  * `21 = (2^6 − 1)/3` and `728 = 3^6 − 1` for a six-cell target: the 728 is an encoding artifact (the DAG has 1620 nodes; other certificates have 2016 and 616 cards).
  * `722 = 728 − 6` is the collection's manuscript count, a moving snapshot.
  * `κ = 2^−182` comes from `2κ = 2^−75·2^−56·2^−50`.
  * Among integers up to 1000, 14% have the form `(2^a ± 1)/d` or `(3^b ± 1)/d` (up to small adjustments).
## 5b. Self-similar measures and entropy

Papers read: *The entropy-rate dimension formula for self-similar measures on the line* (`dim μ = min{1, h_RW/χ}` with exact overlaps allowed, which settles Varjú's conjecture); *Arithmetic classification and non-Pisot singularity for Bernoulli convolutions*; *Sharp binary-information contraction on the discrete cube* (the Courtade–Kumar conjecture); and *Positive metric entropy for the standard map* (skimmed).

* **The Collatz IFS.** The maps `f_c(x) = (3x+1)/2^c` with Haar weights `2^−c` have letter entropy `H = 2` bits and Lyapunov exponent `χ = 2 − log_2 3 = 0.415` bits. Their stationary law is `X = Σ 3^k 2^(−A_(k+1))`.
  * **No exact overlaps (PROVED).** The slope fixes the length and the total, and the {2,3}-representation `B_w` then fixes the word (2-adic peeling).
  * **The dimension formula does not apply directly.** `f_1` expands and the alphabet is countable. Its prediction `dim = 1` is CONDITIONAL on extensions to systems that contract only on average.
* **The second Moran root is a tail index (PROVED).**
  * `E[(3/2^c)^s] = 3^s/(2^(s+1) − 1)` equals 1 at `s = 0` and `s = 1`.
  * Goldie's renewal theorem gives `P(X > t) ~ 1.911/t`, matched by Monte Carlo (1.887, 1.907, 1.912 at `t = 100, 300, 1000`).
  * Since `g(s) = E[M^(s−1)]` with `M ∈ {1/2, 3/2}`, the root `s = 2` of our Moran function `g` is the Kesten–Goldie tail index 1 of the stationary law.
* **Collisions versus exact overlaps (PROVED).**
  * For reduced `u`, `2^K f_u(−1)` has 2-adic valuation 1. So colliding words have equal totals and different lengths.
  * Each fibre of `w -> f_w(−1)` at length `n` has at most `n + 1` words. Hence `H(f_(W_n)(−1)) >= 2n − log_2(n + 1)`: collisions cost at most `log_2(n+1)` bits and **cannot change any dimension formula**.
  * In `Q_3`, the address `rev(v)1^∞` codes `f_v(−1)`. Collisions are exactly where this 3-adic coding fails to be injective, on the null set of addresses ending in `1^∞`. **This is the 3-adic analogue of a dyadic rational having two binary expansions: the reset switch `(a) ~ (2, a−2)` is Collatz's `0.0111… = 0.1000…`.**
  * The uniform families of THM-4555 are the shift along that tail.
  * These statements come from the session's reader and are not independently audited.
* **Collision counts.** The primitive sporadic counts continue as 3187, 6128, 11845, 23036, 44921, 87835 for totals 21..26 (FINITE-EXACT). The growth ratio rises from 1.90 to 1.955, passing `2^h = 1.932`.
  * "×1.9 per unit total" is therefore pre-asymptotic.
  * The share `spor_K/2^K` falls only from 0.00159 to 0.00131. Whether it tends to 0 is open.
* **min ρ is a 3-adic load (PROVED).**
  * On the critical line the entropy per odd step is `h log_2 3` bits. Hence `2^(h log_2 3)/3 = 3^(h−1) = min ρ`: THM-4554's constant is the number of typical critical words per letter, divided by the factor 3 by which a 3-adic cell shrinks.
  * The real-line dimension formula breaks down exactly on our line: there `χ = 0`, the two Kesten roots merge, and no stationary probability exists.
* **Similar equations with one mechanism.** `σ(2^a − 1)/a -> 2/(2 − log_2 3) = 4.819` (HEURISTIC; observed 4.79 at `a = 6000`). This is exactly `H/χ` for the Collatz IFS: both are the mean exponent divided by the drift.
* **The cube.**
  * The Terras parity bits are `v_k = x_k XOR c_k(x_(<k))`, with `c_k` a carry.
  * Only `v_0` and `v_1` attain the Courtade–Kumar bound under bit noise. Every later parity falls strictly below it (PROVED).
  * The information each parity carries about noisy bits decays at nearly the rate of a random Boolean function (avalanche; NUMERICAL).
* **A conditional route.** Combining Kaimanovich–Vershik, Brofferio's boundary theorem and Ledrappier's entropy–dimension relation at each place would give `dim_3` of the critical Syracuse measure `= h`. Two inputs are unchecked: an S-adic Ledrappier–Young formula and the boundary at `χ = 0`.
## 5c. Information, computation and the rewrite compiler

Papers read:
* *Conditional information under deterministic coordinate sweeps* and *Optimal-order mixing of the Thorp shuffle*: the Thorp shuffle mixes in `Θ(log n)` sweeps, and memory lasts exactly one sweep.
* *Prefix instructions and incompressible flows* and *Incompressible box transport and finite computation*: halting is embedded in smooth incompressible flows.
* *Simulating one-tape time in two-fifths power space*.
* *L = RL = BPL* (skimmed).

The facts below are as established by the session's reader. They are not independently audited unless marked.

* **Collatz is one deterministic sweep (PROVED; exhaustive at `k = 16`).**
  * The Terras map `n mod 2^k -> parity word` sets `v_j = x_j XOR g_j(lower bits)`. It is a binary-tree automorphism, an element of the Sylow 2-subgroup of `S_(2^k)`.
  * The transfer operator averages the two preimages of each class, so `K` steps kill every zero-mass density mod `2^K`. This is the Thorp one-sweep identity with zero noise.
  * Hence `H(n mod 2^K | first j parity bits) = K − j`. A single orbit reveals `n` one bit per step, with no loss and nothing to contract.
* **What is not trivial is merging (NUMERICAL).** Flipping the lowest bit of `n` changes parity bit `j` with probability about `(1 − M_j)/2`. Here `M_j` is the probability that the two orbits have merged by step `j` (0.61 at `j = 59`). These merges are exactly what the compiler's switches certify, and HYP-9214 measures them.
* **The product formula as an information ledger (PROVED).**
  * Each branch `(3x+1)/2^a` has `|·|_∞ |·|_2 |·|_3 = 1`.
  * Per odd step, `a` 2-adic bits are erased, `log_2 3` bits are written into the 3-adic record, and `a − log_2 3` bits of size are lost. On average this reads `2 = 1.585 + 0.415`.
* **Compressibility and its limit (PROVED).**
  * Any `n` with glide `>= log_2 n` has Kolmogorov complexity `<= h log_2 n + O(log log n)`. A divergent orbit has at least `log_2 x − O(1)` such points below `x`, while counting allows up to `x^0.95`. Closing that gap is OPEN.
  * More decisively, the parity sequences of the negative cycles `−1, −5, −17` are computable, come from integers, and expand at every prefix (margins 0.585, 0.170, 0.095 bits).
  * So **"low complexity plus an h-dimensional null set" can never exclude divergence: any proof must use positivity.** This is an information-theoretic form of the repo's sheet-blindness barrier.
* **Universality (ANALOGY).** Generalized Collatz maps are universal (Conway; Kurtz–Simon). The flow papers make particle reachability `Σ^0_1`-complete. Undecidability of a class says nothing about the instance `3x + 1`.
* **The compiler as a compressor.**
  * **PROVED:** the greedy back-reference parse is depth-optimal.
  * **PROVED:** "delete `D` trailing ones" is a pop of `D` digits on the 2-adic digit stack, so the compiler is a string-rewriting system whose total termination is Collatz (compare Yolcu–Aaronson–Heule).
  * **PROVED (dichotomy):** for rules guarded by residue classes, the uncovered set on each run length is finite or of positive density.
  * **PROVED (guard inequality):** a rule with descent factor `s` fires on a set of `(2,3)`-adic measure at most `s`, so each certified bit of descent costs at least one bit of guard.
  * **NUMERICAL, the archived compiler rerun unchanged:**

    | `N` | 10^4 | 2·10^4 | 4·10^4 | 8·10^4 |
    |---|---|---|---|---|
    | seeds | 239 | 463 | 893 | 1638 |

    The seed density per dyadic shell falls like about `0.65/log_2 n`. This is LZ78-like decay with a linearly growing dictionary, and shows no sign of a finite certificate.
  * **NUMERICAL, an independent check of THM-4554:** the backward-minimal fraction among non-multiples of 3 below `10^7` is 0.2970, inside the bracket `[0.28820, 0.29912]`.
* **A 3-adic Thorp problem (FINITE-EXACT + OPEN).**
  * The worst-case covering exponent `A_cov(j)` is the largest, over classes mod `3^j`, of the smallest word total landing in the class. It is 2, 6, 9, 11, 14, 16, 18, 20, 22, 25, 27, 29 for `j <= 12`.
  * Proving `A_cov(j) = O(j)` would be the 3-adic analogue of Thorp's `Θ(d)` mixing.
## 5d. Diophantine and transform papers

Papers read: *The irrationality exponent of π is 2*; *Short Egyptian fractions*; *Subset Sum in time 2^(0.49n)*; *The additive indecomposability of the primes*; *An explicit power saving for the exact discrete Fourier transform*; *Single-fold Diophantine representations* (skimmed). All are as read by the session's reader and not independently audited.

* **`μ(log_2 3)` (ANALOGY-level assessment; CONDITIONAL implications).**
  * The `μ(π) = 2` method (an interpolation determinant on `(e^z, z, …, z)`, with Taylor rows tested against a common function because `e^(2πij) = 1`) plausibly extends to `μ(log r) = 2` for rational `r`. That is the reader's reading, not the paper's claim.
  * It does not reach `log_2 3`, which needs the homogeneous form `K log 2 − d log 3`. Its proof is also ineffective.
  * Known: `μ(log_2 3) <= 5.117` (Wu–Wang 2014, from the zbMATH review). The figure 5.125 is Salikhov's bound for `μ(log 3)`.
  * If `μ(log_2 3) = 2`, then `|2^K − 3^d| >= (log 2/2) 3^d d^(−1−ε)`, a cycle's least element would satisfy `x_min < d^(2+ε)/(3 log 2)`, and gaps between 3-smooth numbers would be `>= n/(log n)^(1+ε)`. None of this would be explicit.
* **Collisions are vanishing S-unit sums (PROVED; FINITE-EXACT census).**
  * A collision is a vanishing sum of `ℓ + ℓ' + 2` signed {2,3}-units. For example `2^7 + 3^3 = 2^5 + 3·2^4 + 9·2^3 + 3` is our `(8, c) ~ (4, 1, 1, c+2)`.
  * Each projective class holds at most one collision, so by the S-unit theorem the non-degenerate collisions are finite for each number of terms. This is ineffective.
  * Census for total `<= 23`: 1, 5, 27 and at least 137 non-degenerate collisions with 6, 8, 10, 12 terms. The 8-term ones are 439, 599, 1063, 8189, 8207, stable since total 15.
  * About half of all collisions are degenerate, and their minimal vanishing subsums are mostly **Catalan's `9 = 8 + 1` (73%)**.
  * So the census grows because the number of terms grows with the total. Counts: 172072 values at total 27, with the ratio climbing to 1.959.
* **Representation counts (FINITE-EXACT).**
  * For fixed length the representation is unique (2-adic peeling, i.e. Terras).
  * The number of representations of odd `N <= 2·10^5` is 0, 2, 4 or 6, with frequencies 83793, 13948, 2224 and 36. The root twin doubles every count.
* **Egyptian fractions, subset sum, indecomposability of the primes (ANALOGY, no transfer).**
  * Membership in our representation family is decided in `O(ℓ)` by peeling: a 2-adically superincreasing knapsack.
  * Chain values are equidistributed modulo 5 to 17 and miss only `0 mod 3`.
  * The Egyptian-fraction descent needs density `> 1/2`, while the represented odd `N` have density about 16%.
* **Exact DFT.**
  * The power saving is again a subcritical recursion: `θ = log_m(m − Δ/2^71)` with `m = 10^6`, built on `C = (1/2)[[1+i, 1−i], [1−i, 1+i]]` with `C^2 = swap`.
  * It uses Good–Thomas and Bluestein, and no Rader step.
## 6. Information compression: one reading of all of the above

The owner asked us to think in terms of information compression. Seven compressions run through the session's objects. Each is typed.

1. **A Collatz trajectory is the evaluation of a {2,3}-representation (PROVED, definitional).**
   * `U^j(n) = (3^j n + B'_w)/2^(A_w)`, and the inverse numerator `B_w = Σ 3^(j−1−i) 2^(S_i)` is built by the Horner recursion `B <- 3B + 2^(S_i)`.
   * Its terms form a 3-smooth antichain (Erdős–Lewin type), not the divisibility-ordered double-base chains of fast scalar multiplication.
   * In a double-base representation both `×3` and `÷2` are free exponent shifts, so the Collatz map's whole difficulty is the `+1`, a carry. This is also the step on which exact multiplication ends (carry propagation after the transforms).
   * THM-4555's collisions are coincidences of two such representations: vanishing sums of signed {2,3}-units.
   * By the S-unit theorem (Evertse; van der Poorten–Schlickewei) the non-degenerate ones are finite for each number of terms. That is PROVED but ineffective; the census finds 1, 5, 27, at least 137 for 6, 8, 10, 12 terms.
   * The degenerate ones come mostly (73%) from Catalan's `9 = 8 + 1`.
   * The primitive sporadic collision values grow by a ratio drifting from 1.90 toward 2 (1.955 at total 26), against `2^(A−1)` words. Whether their share tends to 0 is open.
   * Either way, collisions cost at most `log_2(n+1)` bits of entropy at `−1` (section 5b).
2. **Certificates compress like back-references (ANALOGY with exact content).**
   * A switch `n => m` replaces `n`'s certificate by a pointer to a smaller certified `m`, plus a finite word pair. This is an LZ77-style back-reference.
   * Uniform switches are dictionary entries that serve infinitely many sources.
   * The compiler's 239 residual seeds below `10^4` are the literals.
   * On the Mersenne line the dictionary is periodic (THM-4556 (iv)). Its certified coverage is `>= 0.1199` at `K = 20`. The literals are the plateau leaders: 104 of the first 6000 exponents.
3. **Coalescence (HYP-9213, HYP-9214).**
   * The Collatz graph is a tree, so trajectories eventually merge. The finding here is that they merge *at equal odd-step time* far more often than independence would suggest, and more often as the numbers grow: 0.19 at 16 bits, 0.82 at 2048 bits.
   * The stopping-time function of `2^a − 1` has about `A^0.36` distinct values, and 98% of odd exponents in `[10^3, 6·10^3]` share their value with a smaller one.
   * In compression terms, the information in a long trajectory is mostly shared with a smaller one.
4. **Survivor entropy (PROVED, THM-4554).**
   * Non-descending parity words have entropy `h = 0.949956` bits per step: `W_k = Θ(2^(hk) k^(−3/2))`. So describing a forward survivor saves `1 − h = 0.050` bits per step.
   * The backward sieve saves nothing: its survivors have Haar measure `s_∞ ∈ [0.28820, 0.29912]`.
   * The Moran duality is the statement of *which side compresses*.
5. **Subcritical recursions (ANALOGY).**
   * The multiplication paper's saving is a recursion with `s/W < m` subcalls per scale, i.e. a similarity dimension `log_m(s/W) < 1`. The skeleton is the same as the Moran equation whose root is `h`.
   * Its rank saving routes the identity on 161700 triples through a rank-100 bottleneck plus a sparse correction (`I = MM^T + A_1` over `F_2`). This is Alon's orthogonal-representation construction, the linear-algebra family of Lovász's theta and Shannon capacity.
   * Collatz has no such low-rank channel. The `+1` couples all digits (carry), which is why the double-base frame does not trivialise it.
6. **Frames (ANALOGY).**
   * The paper stores `D_v g` in place of `g`, so that each operation is a shift in its own frame.
   * The 3-adic clock is such a frame for the backward map: in the index coordinate `a -> 2^a mod 3^(n+1)`, multiplication by 2 is a signed cyclic shift (`2^(3^n) ≡ −1`), exactly Nussbaumer's `y` in `Z[y]/(y^r + 1)`. In section 2's CRT splitting the tower's numbers fall into a cyclic half and a negacyclic half.
   * The forward map has the mirror frame, the 2-adic clock of 3: `a mod 2^(K−2)` decides the word of `2^a − 1` (THM-4556 (iv)).
   * No single frame serves both directions. That is the "two places couple only through size" barrier in another language.
7. **Merge certificates for consecutive integers (ANALOGY).** See section 5 on the joint Dickman law.

## 7. Toward proofs: what is now usable

* **Closed today.**
  * HYP-9162, by THM-4557. That is a new proof of Hanaki's 2020 theorem; the hypothesis was already a corollary of the literature.
  * The odd-shift question, by THM-4556 (iii).
  * The 3-multiple switch, negatively in the uniform sense (THM-4555).
* **Now precise.**
  * The reset-2 debt is the event that a debt state `3·2^h M + 1` meets `M`'s orbit at the right time. On the Mersenne line this is the plateau continuation (THM-4556 (ii)).
  * Its probability grows with size (HYP-9214).
  * **A proof of HYP-9214 in a logarithmic-density form would make the rewrite program's obstruction a density-zero set.** That is the natural "almost all" target, in the spirit of Tao's almost-all theorem but for *coalescence* rather than descent.
* **A new finite lever.**
  * THM-4556 (iv)–(v) turns Mersenne switching into a finite automaton per depth `K`.
  * The certified density `δ_K` (0.1199 at `K = 20`) can be pushed by tree search over uncertified classes only, since certified classes stay certified.
  * The increments `δ_K − δ_(K−1)` (about 0.006 and slowly decreasing) decide whether `δ_K -> 1`.
* **Not usable (typed negatives).**
  * The Artin and prime-predecessor theorems produce primes with a huge prime factor of `p − 1`. Our objects are exponentially thin.
  * The circulant Hadamard machinery needs a regular group, and `T_k` has none for `k >= 4`.
  * The joint Dickman law does not see 2-adic merging (section 5).
## 8. Verdicts

| claim | status |
|---|---|
| `U^a(2^a − 1) = oddpart((3^a − 1)/2)`; reset pairing in ternary repunits; Mersenne debt states | PROVED (THM-4556 (i)–(ii)) |
| `σ(M_a) = σ(M_(a−1))` iff lag 1 in stopping time; merge ⇐ meeting before 1; converse | PROVED; converse FINITE-EXACT (`a <= 6000`) |
| least shift odd; nearest partner `2^(2k) − 1 = 3(4^k − 1)/3`; partner sets are unions of reset pairs | PROVED (THM-4556 (iii)) |
| 2-adic periodicity mod `2^(K−2)`; propagation of collision switches | PROVED (THM-4556 (iv)) |
| certified density `>= 15719/131072` of switching odd exponents | FINITE-EXACT (THM-4556 (v)) |
| 104 values of `σ(2^a − 1)` for `a <= 6000` (`~A^0.36`); 0.980 share | FINITE-EXACT; conjecture HYP-9213 |
| debt resolution 0.19 to 0.91 from 16 to 8192 bits | NUMERICAL; conjecture HYP-9214 |
| `Aut(D(T)) = Aut(T)` for DRTs on `>= 7` vertices; `Aut(T_k) = F_21` (HYP-9162) | PROVED (THM-4557; known, Hanaki 2020; new proof) |
| CRT splitting of the clock tower into cyclic and negacyclic halves; `19·27 = 2^9 + 1`; Zsigmondy and Catalan at 63 | PROVED (classical identities) |
| 7 is the only Mersenne prime with `⟨2⟩ = QR_p` | PROVED |
| joint Dickman law for `(P^+(n), P^+(U(n)))` over odd `n` | CONDITIONAL (on the OpenAI paper; reader's derivation) |
| collisions = vanishing S-unit sums; non-degenerate ones finite per term count | PROVED (ineffective) |
| collisions cost `<= log_2(n+1)` bits of entropy; reset = "two expansions" | PROVED (reader; not independently audited) |
| second Moran root = Goldie tail index; `min ρ` = 3-adic load | PROVED |
| negative cycles block every information-only exclusion of divergence | PROVED (reader; not independently audited) |
| subcritical recursions in multiplication/DFT vs Moran dimension; frames; Rader index maps | ANALOGY |
| `κ = 2^−182`, Snaky's 21/728, the 722 manuscripts, `−91/9` | NUMEROLOGY (dissolved) |
| Collatz | OPEN |

## 9. Audit record

**THM-4556, HYP-9213, HYP-9214 (independent audit, own code).** All mathematics passes. The auditor recomputed `σ(M_a)` for `a <= 6000` independently, reproducing the OEIS A193688 b-file exactly. Corrections, all applied and logged as MISTAKE-572:

* (ii): only the lag-1 stopping-time equivalence is proved. The merge converse is checked for `a <= 6000`.
* (iii): "nearest partner"; partner sets are unions of reset pairs; the explicit size bound for "no merge before time `a`".
* (iv): the hypotheses `|u'| = |u| + D` and `3^(a−1+|u|) > 2^K`; "by step … exactly then if `u` is shortest".
* (v): the exact fraction `15719/131072` among odd exponents; even `D` certifies the same classes.
* (vi): 0.95 corrected to 0.980.
* Remarks: "poles" replaced; the false claim that every switch lands on a multiple of 3 restricted to least-shift switches; `F_21` typed as ANALOGY.
* HYP-9213: "shares `σ`", not "merges".
* HYP-9214: the median-lag row is noisy; "within 10 steps" is a fraction of sources; the sibling model has `k ∈ Z \ {0}`; a committed script regenerates the table.

**THM-4557 (independent audit, own code: dreadnaut plus a nauty-free enumerator).**
* Every step passes, and `D(T_k) = T_(k+1)` holds entrywise.
* Checked on `T_k` for `k <= 9` (`|Aut T_9| = 21`, 511 vertices), Paley primes to 83 and `P_27`.
* Checked on **438 further doubly regular tournaments** of orders 15–35, neither Paley nor tower.
* **Major correction: the theorem is known** (Hanaki 2020, Theorem 3.4; Faradžev–Klin–Muzichuk 1994 for intransitivity). Applied, and logged as MISTAKE-571.

**Not independently audited:** the readers' own derivations quoted in sections 4–5d (the Dickman divisibility-class extension, the entropy bound at `−1`, the information-theoretic propositions, and the S-unit census). They are typed as such.

## 10. Reproduction

```bash
python3 04-computation/experiments/sevens_20261006_mersenne_plateaus.py 6000   # ~90 s
python3 04-computation/experiments/sevens_20261006_drt_doubling_aut.py         # ~5 s, needs dreadnaut
python3 04-computation/experiments/sevens_20261006_debt_resolution_trend.py    # ~6 s
```

The readers' scripts (Snaky re-evaluation, Dickman sieve, collision census, Artin/Lucas certificates, Courtade–Kumar computations) are in the session scratchpad and are not part of the repo record.
