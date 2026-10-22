# Audit, second pass: the post-audit continuation of THM-4750

**Independent adversarial audit (second pass), 2026-10-21.**

**Scope.** This pass covers only the material that mac-mini-2026-10-21-fanopsp added after the first audit (`camion_function_audit_20261021.md`, whose notation is used here):

* §5b, Statement 5′: Teichmüller rigidity, `Λ_p(a) = σ_p(a) − ⟨a⟩^(−1/2)`, the Camion lift `κ`, and the quotient decomposition;
* the reworded 5(a) and 5(e);
* the upgrade of 6(b) and 6(c) to "PROVED in outline";
* the r = 4 hypergraph data in §6b;
* the base-3 and connected-graph Camion–Wieferich lists;
* `cw_flint_np.py`;
* the "p-adic avatars" paragraph of §7.

**Rules followed.** The auditor is an independent Claude subagent. It edited no canon, hypothesis or results file and ran no git write command. C code was compiled outside the repository. Scripts and outputs are in `04-computation/experiments/camion_wieferich_20261021/audit/`, files `c1_*` to `c9b_*`.

**Files audited.**

* THM-4750 = `01-canon/theorems/THM-4750-the-camion-function-half-speed-interpolation-and-the-camion-lift.md`
* HYP-9361 = `05-knowledge/hypotheses/HYP-9361-camion-function-open-laws-split-inert-equidistribution-tower-heights-uniformity-k-speed.md`
* Session25 = `05-knowledge/results/session25_spectrum_wiring_20261021/session25-congruence-review/ADDITIVE_PRIME_POWER_PROOF.md`

References of the form `THM-4750:381` give a file and line number. Line numbers refer to THM-4750 and HYP-9361 as they stood when this audit began (working tree at 20:49 MDT, uncommitted), before the session started applying this audit's findings.

## 0. Verdict

**The new mathematics of §5b is correct.** I re-derived (a)–(d) line by line and found no error. A few one-line steps are missing (findings 1 and 3), and the proof in fact gives more than it states (finding 1, OK+).

**Every new FINITE-EXACT number I tested reproduces with independent code.** That covers rigidity, `Λ_p`, `κ`, the height formula, the decomposition (d), the r = 4 table, the forced roots, and the base-3 and connected lists. The independent code uses a different recurrence (`T′ = S′·T²`), FLINT power series mod `p²`, a different hypergraph recurrence, and gcd root counting. Highlights:

* **Rigidity holds directly.** `s_a(p^k) ≡ s_(a′)(p^k) (mod p^(k+1))` for `a ≡ a′ (mod p)`, with no interpolation and up to mod `3^10`, `5^7`, `7^6`, `11^5` and `13^5`.
* **Designed heights.** The integer bases `κ(a) + 3p^h` have tower height exactly `h`, for every `h` tested.
* **All bases.** Over all 794 integer bases below `5^4` and `7^3`, the tower height equals `v_p(a − κ(a))` with no exception.
* **Session reruns.** Every committed output of the new session scripts reproduces byte for byte.

**One MAJOR finding.**

* **P2-M1. Statement 6(b) is false as written, under its new label "PROVED in outline"** (THM-4750:381).
  * The congruence `f(n + 3·2^j) ≡ f(n) (mod 2^(j+2))` has an exception set that grows linearly in `j`. It fails at `n = 5, 6` for `j = 11`, at `n ≤ 8` for `j = 14`, and at `n ≤ 10` for `j = 17`.
  * The text records exceptions only for `j ≤ 8` ("`n ≤ 3` for `j ≤ 7`, `n ≤ 4` for `j = 8`").
  * The outline proves exactly the statement with the range `n ≥ 2(j−1)/3`. Finding 13 has the details.

**MINOR findings.** Eleven findings below carry a MINOR tag: 1, 3, 6, 7, 8, 17, 20, 22, 23, 24 and 25. The most substantive is finding 20.

* **Finding 20.** The remark that connected graphs have "their own lift and therefore their own list" (THM-4750:364, HYP-9361:71) does not hold in the sense of §5b.
  * The circulant target `2^((n−1)/2) − 1` makes the tower relation drift, so `E(p^r) ≡ εE(p^(r−1)) + (1 − ε)`.
  * For `(2/p) = −1` the prime powers that pass are not an initial segment: at `p = 5`, 5 and 125 pass and 25 fails.
  * Both listed primes, 5 and 197, have `(2/p) = −1`.
* **The others.** They are a missing existence line in 5′(a), the `p`-adic-base conventions, an almost tautological cross-check, wording at `p = 7`, a mislabelled printout, a title inconsistency in HYP-9361 (d), and the 10^6 search being cited before it has finished.

**Recommended labels.**

* **§5b (a)–(d): PROVED.** Add the lines in findings 1 and 3. The finite-level form (finding 1) can replace the limit statement.
* **6(c): PROVED in outline is honest.** The range qualifier is present, and every numerical ingredient checks out.
* **6(b): PROVED in outline only after the range is added** (finding 13).
* **FINITE-EXACT, confirmed:** the r = 4 data; the base-3 list; the connected list, now extended to `p ≤ 10^4` with no new hits.
* **The 10^6 search: pending.** At audit time it has no hits on the sub-ranges already covered (finding 22).

## 1. Findings

*Numbering note.* In a draft of this report, read while it was being written, the findings from the MAJOR one onward were numbered one lower (the MAJOR finding repeated the number 12). The final numbering is: P2-M1 = finding 13, connected graphs = finding 20, the 10^6 search = finding 22. Citations of "second audit finding 12" or "finding 19" in THM-4750 or HYP-9361 mean findings 13 and 20.

### Statement 5′: rigidity, the limit, the lift (checklist 1)

**1. OK, with a MINOR gap and an OK+ sharpening. The rigidity proof of 5′(a) is correct** (THM-4750:321–326).

*The steps check out.*

* Since `a_(j,0) = C(C(0,2), j) = 0` for `j ≥ 1`, the term `U = T_a − A` is a combination of `p^j δ_d ∗ E_α` with `1 ≤ d ≤ 2j`.
* Products in the Hurwitz ring add the `d`'s, by `δ_d ∗ δ_e = C(d+e,d)δ_(d+e)` and `E_α ∗ E_β = E_(α+β)`. So every term of `Σ_(k≥1)(−A^(−1)U)^k A^(−1)` has `d ≥ 1`.
* The grade-0 term `A^(−1)` is the same for all lifts of `a mod p`. Indeed `G = θR_0` depends only on `p`, `ε` and the precision, and `A` is the group-algebra element whose image is `ω(a)^C(n,2)`.
* The representation of `A^(−1)` need not be unique, but its image, the Hurwitz inverse of `t_0`, is unique, and only that image is evaluated. In particular `s_(ω(a)) = e − A^(−1)` exactly.

*The bound is stronger than written.* At `n = p^k`, `v_p(C(p^k,d)) = k − v_p(d)`; the text's `k − v_p(d!)` is weaker but correct. Since `d ≤ 2j` and `p^j > 2j`, `v_p(d) ≤ j − 1`. So every `d ≥ 1` term has valuation at least `j + k − v_p(d) ≥ k + 1` for unit nodes, and more for nonunit nodes. Hence, for every `k ≥ 1`,

    a ≡ a′ (mod p)  ⟹  s_a(p^k) ≡ s_(a′)(p^k)   (mod p^(k+1)).

This finite-level rigidity is what the data test (finding 5).

*The double limit is harmless.*

* `ε^k s_a(p^k)` is a single sequence in `Z_p`.
* The tower relation for `s_a` itself is Session25 §7 at `n = p^(k−1) ≥ 2k − 1`: `ε^k s_a(p^k) ≡ ε^(k−1) s_a(p^(k−1)) (mod p^k)`. So the sequence is Cauchy, and `ε^k s_a(p^k) ≡ σ_p(a) (mod p^(k+1))`.
* At precision `p^r` no limit in `r` is needed.

*MINOR gap: the existence of the limit is not argued* (THM-4750:310 claims "exists").

* The proof shows only that the grade-`≥ 1` part dies.
* The tower relation above gives existence in one line.
* Alternatively, for `β = θb` with `b ∈ R_0^×`, `ε^k β^(p^k) = θ·ω(b)·⟨b⟩^(p^k)`, which equals `θω(b)` mod `p^r` once `k ≥ r − 1`. Nonunit nodes tend to 0. So the "Teichmüller projection" of `A^(−1)` converges.

**2. OK. Statement 5′(b) is correct** (THM-4750:327).

* `ε_a^k c_a(p^k) = ε_a^(2k) ⟨a⟩^((p^k−1)/2) → ⟨a⟩^(−1/2)`.
* `(p^k − 1)/(p − 1) ≡ k (mod 2)`.
* `x ↦ ⟨a⟩^x` is continuous on `Z_p`.

**3. OK, with MINOR omissions. Statement 5′(c) is correct** (THM-4750:328–331).

*The steps check out.*

* With `x = σ_p⟨a⟩^(1/2) ≡ 1`, `v(σ_p − ⟨a⟩^(−1/2)) = v(x − 1) = v(x² − 1) = v(σ_p²⟨a⟩ − 1)`. Since `x + 1 ≡ 2` is a unit, this equals `v(a − ω(a)σ_p^(−2))`.
* Uniqueness of `κ(a)` follows from `κ(a′) = κ(a)` for every `a′ ≡ a (mod p)`.
* *General base.* Session25 §7 gives the law for `S_a`, and `a^(φ(p^r)/2) ≡ (a/p) (mod p^r)` gives it for `c_a`. Together they give `D_a(p^r) ≡ ε_a D_a(p^(r−1)) (mod p^r)`.
* *Initial segment and freezing.* Telescoping gives `ε^k D_a(p^k) ≡ Λ_p(a) (mod p^(k+1))`, so the THM-4746 (9) lemma holds verbatim.
* Numerically: 168 tower relations, the designed heights and the all-bases test of finding 5.

*MINOR omissions.*

* **(i)** `σ_p ≡ 1 (mod p)` is used but not shown. It follows from `σ_p ≡ ε s_a(p) (mod p²)` and `s_a(p) ≡ ε (mod p)`.
* **(ii)** `h_a(p)` is never defined in §5b. Define it as `max{k : p^k | D_a(p^k)}`.
* **(iii)** The Notation (THM-4750:306) defines everything for an integer "base `a` prime to `p`", while 5′(a) and 5′(c) speak of `p`-adic lifts and "the unique `p`-adic base". The fix needs one sentence:
  * for `a ∈ Z_p^×`, `s_a` is defined by the same recurrence;
  * `s_a(n) mod p^r` depends only on `a mod p^r`, since it is an integer polynomial in `a`;
  * "`p^k` is a Camion pseudoprime to base `a`" means `p^k | D_a(p^k)` in `Z_p`.
  * Then `κ(a)` is a genuinely `p`-adic object. Integer bases `≡ κ(a) (mod p^K)` have height `≥ K`, and an integer base passes every `p^k` only if it equals `κ(a)`.
* **(iv)** `c_a(n)` for even `n` is not stated. It should be 0, as for base 2.

**4. OK, with an OK+ formula. Statement 5′(d) is correct** (THM-4750:332). The tower gives `εD(p) ≡ Λ_p (mod p²)` and `⟨a⟩^(−1/2) ≡ 1 + p q_a(p)/2`.

*OK+: an explicit formula for the Camion residue* (relevant to HYP-9361:128, "closed form?").

* As in `cw_fast.c`, `s_a(p) ≡ a^C(p,2) + p·u_p(a) (mod p²)`, where `u_p(a) = [x^p] 1/T_(a,<p)(x)` in `F_p[[x]]` and `T_(a,<p) = Σ_(k<p) a^C(k,2) x^k/k!`.
* Also `a^C(p,2) = (ε + pe′)^p ≡ ε (mod p²)`.
* Hence

      ρ_C(a) ≡ ε_a · [x^p] (Σ_(k<p) a^C(k,2) x^k/k!)^(−1)   (mod p).

  The right side visibly depends only on `a mod p`. This gives an elementary proof of the level-1 rigidity, and a concrete object for the Gauss-sum question.
* Checked at 301 primes (finding 9).

**5. OK. The finite-level rigidity holds in direct computation** (`c1_rigidity_direct.c`, `c1_run.sh`, `.out`; `c1b_rigidity_residue1.out`).

*Method.*

* `s_a(n) mod p^(k+1)` up to `n = p^kmax`, by the recurrence `s(m+1) = T(m+1) − Σ_(k<m) C(m,k) s(k+1) T2(m−k)`, with `T2 = T ∗ T`. This recurrence comes from `T′ = S′T²`, an identity different from the session's Wright form. It was self-tested against Wright for `n ≤ 60`.
* Each run uses 7 lifts of the residue: `a, a+p, a+2p, a+5p, a−p+p²`, the Teichmüller lift `ω(a)`, and a pseudo-random lift.

*Rigidity results.*

* 95 level checks, plus 20 more for the residue 1 (where `σ_p(1) = 1` exactly). All passed. The residues are 2, 3, 5 and `p − 1` at `p = 7, 11, 13`; 2, 3, 4 at `p = 5`; and 2 alone at `p = 3` and `p ≥ 17`.
* The highest levels reached: `p = 3` mod `3^10`; `p = 5` mod `5^7`; `p = 7` mod `7^6`; `p = 11` mod `11^5`; `p = 13` mod `13^5`; `p = 17, 19, 23` mod `p^4`; `29 ≤ p ≤ 43` mod `p³`.
* All 168 tower relations (for `s` and for `D`) hold.

*Height results* (`c2_lambda_kappa.out`, `c3_heights_general_residue.out`, `c3b_height_all_bases.out`).

* `Λ_p` (direct, from `D(p^k)`) equals `σ_p(ω(2)) − ⟨2⟩^(−1/2)` at every level, for 13 primes `3 ≤ p ≤ 43`. Here `⟨2⟩^(−1/2)` is computed as `⟨2⟩^((p^(K−1)−1)/2)`, not by the session's Newton root.
* `v(Λ_p) = v(2 − κ(2))`: 2 at `p = 3` and `13`, 1 otherwise.
* My `κ(2)` agrees with the session's at 7, 11 and 13.
* **Designed heights.** The base `κ(a) + 3p^h` gives `v(D(p^k)) = k+1` for `k < h` and `= h` from `k = h` on, for every `h` tested:
  * residue 2: `h = 1..6` at `p = 5`, `h = 1..5` at 7, `h = 1..4` at 11 and 13;
  * residues 3, 5 and `p − 1` at `p = 5, 7, 11, 13`.
* **All bases.** All 500 integer bases below `5^4` and all 294 below `7^3` give tower height (direct) `= v_p(a − κ(a))`, with 0 mismatches.
  * The histograms are exactly the proportions `(1 − 1/p)p^(1−h)` below the cap, and `p^(1−cap)` in the capped bin: `400/80/16/4` at `p = 5` and `252/36/6` at `p = 7`.

**6. MINOR. The "`Λ_p` cross-check" of THM-4750:338 is nearly tautological.**

* `camion_lift.py` interpolates `s` and `D = s − c` on the same nodes.
* The circulant class sequence is geometric, `2^((p−1)/2)(1+pz)^l`, and its interpolation at `z_∞ = −p/(p−1)` is `(1+pz)^(−1/(p−1)) = ⟨2⟩^(−1/2)` identically, to the precision used.
* So "MATCH" checks only the sign conventions and the square-root routine.
* The genuine test is finding 5: direct `D(p^k)` against `σ` computed from the Teichmüller base.
* The second sentence of that line, "equals `D(p²) mod p³`", is right, but its evidence is `tower_height_direct_checks.out`, not the §5b outputs.

**7. MINOR. The `κ` test at `p = 7` is mis-described** (THM-4750:339: "`v_p(D_κ(p^k)) = 2, 3, 4` … at `p = 7, 11, 13`").

* `camion_lift_check.out` line 1 shows `[2, 3, 3]` at `p = 7`, capped at `K = 3`. The interpolation there cannot go beyond `K = 3`, because it needs `p > 2K − 1`.
* The claim is true. Audit c2 gives `[2, 3, 4, 5, 6]` at `p = 7` with `K = 6`, by direct values. Cite it, or say "capped at 3 for `p = 7`".

**8. MINOR, cosmetic. The Lagrange-degree printout is mislabelled** (`camion_residue_symmetry.py:29`, cited at THM-4750:343).

* The script prints "min Lagrange degree: 9 (generic 10)".
* With the value at 0 free, the generic minimum is `p − 2`: the `x^(p−1)` coefficient `−Σ f` can always be killed.
* The observed degrees 9, 11, 15, 17 at `p = 11, 13, 17, 19` are exactly `p − 2`. So "generic" in the THM is right, but the printout suggests degree deficiency.

**9. OK. The decomposition (d) holds at all 301 primes `5 ≤ p ≤ 2000`, computed three ways** (`c4_decomposition_d.py`, `.out`).

* **(A) Direct:** `ε s_2(p)` and `ε s_ω(p)` mod `p²` by the `T′ = S′T²` recurrence.
* **(B) FLINT:** the formula of finding 4.
* **(C) Session:** `camion_residue_2000.out`, line by line in `(ρ_C, q_2, q_C)`.

Results:

* 0 failures in (A) against (B) and in the decomposition.
* 0 differences from the session.
* `ρ_C(2) = 0` exactly at 199, 809 and 929.
* The session's test with `s_ω` against `s_2` is itself a genuine level-1 rigidity test, not a tautology.

**10. OK. The session reruns are byte-identical.** All compiled or run in scratch:

* `camion_lift.py`, `camion_lift_check.py`, `camion_lift_map.py`;
* `camion_residue_symmetry.py`, `camion_residue_fish.py`;
* `camion_residue.c 5 2000`;
* `hyper_r_second_order.c 4 5 2200`;
* `cw_base.c` and `cw_connected.c`, which match apart from the timing lines;
* `hyper3_class2_check.py` and `two_adic_regulator.py`.

The newly committed `second_order_direct_3_59.out` records `direct-law checks=17158 failures=0`, which closes the first audit's finding 9.

### 5(a) and 5(e) wording (checklist 2)

**11. OK. Both match the first audit's recommendations** (THM-4750:286, 292–296, 299).

* 5(a) states `ε^a D(p^a) ≡ Λ_p (mod p^(a+1))` and proves it by the contraction argument.
* 5(e) states `p ≥ 2K − 1`, the data `D(p) mod p^K` and `D(p + j(p−1)) mod p^(K−1)` for `j ≤ 2K − 4`, and the `K = 3` formula.

*Numerical check* (`c2_lambda_kappa.out`). My own implementation of the full interpolation and of the cheap `K = 3` formula equals the direct `Λ_p`:

* mod `p³` for all 12 primes `5 ≤ p ≤ 43`;
* mod `p⁴` for `7 ≤ p ≤ 23`;
* including the boundary cases `p = 2K − 1`.

*OK+.* The `mod p^(a+1)` congruence also follows from the tower relation alone, by telescoping: `Λ − ε^aD(p^a) = Σ_(j>a)(ε^jD(p^j) − ε^(j−1)D(p^(j−1)))`, with each term `≡ 0 (mod p^j)`. That route needs no Camion function, and it holds for every base and for `p = 3`.

### 6(b) and 6(c) "PROVED in outline" (checklist 3)

**12. OK. The ingredients are stated correctly** (`c8_two_adic.py`, `.out`; arithmetic in `Z/2^160[ζ_3]`).

* **Newton polygon of `T`.** The segments are slope 0 (2 zeros), 3/2 (2), 4 (2), 5 (2), 8 (2) and 19/2 (2). So exactly two zeros lie in the closed unit disc, and the next pair sits on `|x| = 2^(3/2)`.
* **Valuations.** For the root `ρ = ζw` (Newton iteration to `2^160`): `v(w − 1) = 2`, `v(w³ − 1) = 2`, `v(C) = 0`, `v(r − 1) = 2` and `v(t − 1) = 3`, with `N(r) = N(t) = 1`.
* **Remainder.** It is better than "`(3/2)n − O(1)`": `v_2(r_n) = 3n/2 + 3, 9/2, 3, 7/2` for `n ≡ 0, 1, 2, 3 (mod 4)` (all `n ≤ 120`), so `v_2(r_n) ≥ 3n/2 + 3`.
* **Cross-check.** `f(n)` from the truncated series equals the exact `s(n)/n!` for `n ≤ 25`.

**13. MAJOR (P2-M1). Statement 6(b) is false as written, under the label "PROVED in outline"** (THM-4750:381–385).

*What the text says.* The statement is the unqualified congruence `f(n + 3·2^j) ≡ f(n) (mod 2^(j+2))`. The only exceptions it lists are "`n ≤ 3` for `j ≤ 7`, and `n ≤ 4` for `j = 8`", which are the first audit's data limit. A reader would take the congruence to hold for `n ≥ 5` and all `j`.

*Counterexamples* (`c8_two_adic.out`, block (5)).

| `j` | `n` | `v_2(f(n+3·2^j) − f(n))` | needed `j + 2` |
|---|---|---|---|
| 11 | 5, 6 | 12 | 13 |
| 13 | 7 | 14 | 15 |
| 14 | 8 | 15 | 16 |
| 17 | 9, 10 | 18 | 19 |

*The full failure sets for `j ≤ 18`, `n ≤ 60`* (block (4)).

* `j ≤ 4`: none.
* `j = 5, 6`: `n ≤ 2`.
* `j = 7`: `n ≤ 3`.
* `j = 8, 9, 10`: `n ≤ 4`.
* `j = 11, 12`: `n ≤ 6`.
* `j = 13`: `n ≤ 7`.
* `j = 14, 15, 16`: `n ≤ 8`.
* `j = 17, 18`: `n ≤ 10`.

*Consistency.* These agree, for `j ≤ 8`, with the first audit's exact `n ≤ 1500` data, which used a different method.

*Exact rule* (FINITE-EXACT for all `j ≤ 18` and `n ≤ 60`): the congruence holds **iff** `v_2(r_n) ≥ j + 2`. This is precisely what the outline gives:

    f(n + 3·2^j) − f(n) = Tr(Cρ^(−n)(w^(−3·2^j) − 1)) + (r_(n+3·2^j) − r_n),

where the first term has valuation `≥ j + 2`.

*Correct statement.* `f(n + 3·2^j) ≡ f(n) (mod 2^(j+2))` whenever `v_2(r_n) ≥ j + 2`. In particular it holds for all `n ≥ 2(j − 1)/3`, using `v_2(r_n) ≥ 3n/2 + 3`, which the outline must prove (it now says only `−O(1)`). Failures do occur below that range, for every tested `j` with `5 ≤ j ≤ 18`.

*Related sub-claim.* The `mod 2^(j+3)` claim on `n ≡ 1 (mod 3)` is correctly restricted to `j ≤ 8` (THM-4750:382). Its exception set also grows: `n = 7` from `j = 12`, and `n = 10` from `j = 16`.

*Provenance.* The first audit's finding 23 ("the same argument gives the quasi-period") omitted this range. The status line THM-4750:16 should not be read as covering the unqualified congruence.

**14. OK. The label of 6(c) is honest** (THM-4750:388–402).

* The range qualifier is present: "whenever `2 + v_2(n − n*)` is below the decay rate of `r_n`".
* All the numbers the outline needs are confirmed: `v(C) = 0`, `v(r − 1) = 2` and `Tr = Cζw^(−n)(1 − t r^(−n))`. I re-derived the last on `n ≡ 2 (mod 3)`, using `σ(ζ) = ζ²`.
* *One implicit step to add.* `n* = log t / log r` needs `t ∈ r^(Z_2)`, not just `log t ∈ Z_2·log r`.
  * This holds because `N^1 ∩ (1 + 4O)` is torsion-free of rank 1.
  * `r` generates it, since `v(r − 1) = 2` is minimal.
  * `t` lies in it, since `v(t − 1) = 3 ≥ 2`.
  * It also gives `v_2(n*) = 3 − 2 = 1`, as the first audit found.

### r = 4 hypergraphs (checklist 4)

**15. OK. Code review of `hyper_r_second_order.c`: no bugs.**

* *Exponent.* `C(n,r)` is computed exactly in u128 and reduced mod `p(p−1) = φ(p²)`, which is valid because 2 is a unit mod `p²`. The "+PH" is harmless.
* *Recurrence.* The Pascal row `C(n−1, ·)` is updated correctly, so the rooted recurrence is right.
* *Mahler coefficients.* Computed from `r + 2` points by forward differences.
* *Degree-r law test.* `p | a_1..a_r` and `p² | a_(r+1)`.
* *Root evaluation.* The Newton form `Σ c_i l^(i)/i!` with falling factorials and inverse factorials is correct for `i ≤ r < p`.
* *Class range.* `m ∈ [r+3, r+p+1]` covers each residue mod `p − 1` once, starting one above the proved threshold `r + 2`.
* *Prime set.* `gcd(ord_p 2, 24) = 1` is the same set as `gcd(ord_p 2, 6) = 1`, so r = 3 and r = 4 use the same 60 primes, with `Σ(p−1) = 60444`.

**16. OK. The r = 4 table re-derives exactly, with independent code** (`c7_hyper_r.c`, `c7_hyper_analyze.py`, `c7_hyper_r4_raw.out`, `c7_hyper_r3_raw.out`, `c7_hyper_analyze.out`).

*Method.*

* `2^C(n,r)` is built by the Pascal rule in the exponent, `P_j(n+1) = P_j(n)·P_(j−1)(n)`, so no binomial exponent appears.
* Connected counts come from `b′ = T′·T^(−1)` through the Hurwitz inverse.
* Root counts are `deg gcd(g, x^p − x)`, computed with FLINT.

*Table.*

* My parser of `hyper4_second_order_5_2200.out` and my recomputation both give:
  * class 1: 23 / 19 / 16 / 0 / 2;
  * class 2: 0 / 18 / 32 / 1 / 9;
  * class 3: 0 / 21 / 26 / 0 / 13;
  * sporadic: 21 / 19 / 13 / 0 / 3.
* There are 0 per-prime, per-class differences from the session.
* The quartic law has 0 failures, and 50 classes are not genuinely quartic.
* The r = 3 table also reproduces (22/25/1/12, 0/31/0/29, 21/21/0/14), with 58 classes not genuinely cubic.

*Forced roots.*

* Classes 2 and 3 vanish at `l = −1` (`n = 2`, `n = 3`) in 60 of 60 cases each.
* The low bases `m = 2` and `m = 3` already obey the quartic law mod `p²` at all 60 primes. This is the mechanism the THM states (THM-4750:441), now checked for r = 4.
* For r = 3, base `m = 2` obeys the law at all 60 primes and base `m = 3` fails at all 60, in agreement with audit b8.

*Special cases.*

* **The anomalous class 2 at `p = 71`** (3 roots) is `g = 65(x+1)²(x+12)(x+59)`. It is a double root at the forced point, a probability-`1/p` degeneration.
* **The only non-genuine zero classes** are `(23, class 1)`, where `g` is linear (r = 4), and `(191, class 64)` (r = 3).

*The `S_r` numbers are exact.*

* `S_4` fixed points: `3/8, 1/3, 1/4, 0, 1/24`.
* Forced root plus `S_3`: `0, 1/3, 1/2, 0, 1/6`.
* `S_3` fixed points: `1/3, 1/2, 0, 1/6`.
* The expected counts in THM-4750:431–439 follow from these.

**17. MINOR. The title of HYP-9361 (d) contradicts its body.**

* HYP-9361:3 and :86 speak of "forced rational roots at exact zeros **inside the analytic range**".
* HYP-9361:100 and THM-4750:441 say correctly that these zeros lie **below** the proved threshold.
* Suggested wording: "at exact zeros where the class polynomial already holds (below the proved threshold; FINITE-EXACT)".

### New Camion–Wieferich lists (checklist 5)

**18. OK. The rotation congruences are right.** The derivations:

* **Base 3.** `Z_p` acts on labelled semicomplete digraphs on `Z_p`, with orbits of size 1 or `p`.
  * The fixed points are the circulant ones: for each `d ≤ (p−1)/2`, one of 3 states is chosen for all pairs `{x, x+d}`. That gives `3^((p−1)/2)` of them.
  * All of them are strong, because any nonzero difference generates `Z_p`.
  * The strong set is rotation-invariant, so `s_3(p) ≡ 3^((p−1)/2) (mod p)`.
* **Connected graphs.** There are `2^((p−1)/2)` circulant graphs (`S = −S`), and one is connected iff `S ≠ ∅`. So `cg(p) ≡ 2^((p−1)/2) − 1 (mod p)`.
* Both hold with 0 failures for all primes up to `10^4` (c5).

**19. OK. The lists reproduce with an independent algorithm, and the connected list is extended** (`c5_cw_lists_flint.py`, `c5_connected_direct.c`, `c5_direct_crosscheck.py`, with outputs).

*Method.*

* FLINT power series mod `p²`:
  * `s_a(p) = (p−1)!·[x^(p−1)] T′/T²`;
  * `cg(p) = (p−1)!·[x^(p−1)] G′/G`.
* Only coefficients `x^k` with `k ≤ p − 1` enter, so every denominator is a unit.
* Self-tested against exact integers for `p ≤ 31`.

*Results for `p ≤ 10^4`.*

| structure | Camion–Wieferich primes | Wieferich primes (same base) |
|---|---|---|
| base 2 | {3, 13} | {1093, 3511} |
| base 3 | {17, 31} | {11} |
| connected graphs | {5, 197} | — |

* The connected search was extended from the session's 4000 to `10^4`, with no new hits.

*Cross-checks.*

* The mod-`p` identity `D_a(p)/p = u_p − e′` agrees with the mod-`p²` values at every prime.
* Value by value for `p ≤ 2000`, two direct recurrences written for this audit agree with FLINT, with 0 mismatches:
  * `T′ = S′T²` for `s_3`;
  * `b′ = T′T^(−1)` for `cg`.

*Coverage.*

* Base 3: `[5,3000] ∪ [3001,7500] ∪ [7501,10000]`, with 428 + 520 + 279 = 1227 primes. That is all odd primes `≠ 3` below `10^4`.
* Connected graphs: `[3,4000]`, 549 primes. That is all odd primes below 4000.

**20. MINOR, but substantive. "Each Hurwitz-closed structure has its own lift and therefore its own list" does not carry over to connected graphs** (THM-4750:364–366; HYP-9361:71 and :77; `c9_connected_tower.py`, `c9b_connected_rigidity.py`, with outputs).

*(i) The proof does not cover connected graphs.*

* 5′ is proved for `s_a = e − T_a^(−1)`. For connected graphs Session25 gives `Δb = ΔT ∗ T^(−1)`.
* The 5′ argument does not transfer verbatim: the `d ≥ 1` terms are evaluated at `n = p^k − 1`, where `C(p^k − 1, d)` is a unit.
* Numerically the rigidity still holds: `b_a(p^k) mod p^(k+1)` is lift-independent at `p = 3, 5, 7, 11` (`k ≤ 4, 3, 2, 2`). So a lift plausibly exists, but it is unproved.

*(ii) The tower drifts.*

* The circulant target `2^((n−1)/2) − 1` has a constant term, which is not `ε`-twisted.
* Session25's law for `b` (`n ≥ 2r`) therefore gives, for `E(n) = cg(n) − 2^((n−1)/2) + 1`,

      E(p^r) ≡ ε E(p^(r−1)) + (1 − ε)   (mod p^r)     (p^(r−1) ≥ 2r).

  This was verified with exact integers at `p = 3, 5, 7, 11, 13`.
* Take the natural connected analogue of the Camion test: `n` passes if `cg(n) ≡ 2^((n−1)/2) − 1 (mod n)`. For `ε = −1` and `p ≥ 5`, the drift forces `E(p²) ≡ 2 (mod p)`, so `p²` never passes. The prime powers that pass need not form an initial segment: at `p = 5`, `v_5(E(5)) = 2`, `E(25) ≡ 2 (mod 25)` and `v_5(E(125)) = 4`, so 5 and 125 pass and 25 fails.

*(iii) Both listed primes are in this case.* Both 5 and 197 have `(2/p) = −1`, so for them `p² | E(p)` does not make `p²` a pseudoprime.

*What the list is.* `{5, 197}` is a correct list of primes with `p² | E(p)`, a second-order rotation congruence. It is not a list of "height ≥ 2" primes. The base-3 list is unaffected: there is no constant term, and the general-base 5′ theory applies, so `h_3(p) ≥ 2` exactly at 17 and 31 for `p ≤ 10^4`.

### `cw_flint_np.py` (checklist 6)

**21. OK. The script is correct, and validated.**

*Review.*

* `prefix_prod`: a Hillis–Steele inclusive scan. The right-hand side is evaluated before the assignment, so there is no aliasing. With `ar[0] = 1` it yields `k!` for `k ≥ 0`.
* `vpow`: a correct elementwise square-and-multiply.
* Reducing the exponent `C(k,2)` mod `p − 1` is valid by Fermat.
* Overflow: residues are below `2^20`, so products are below `10^12 ≪ 2^63`.
* `e′ = ((2^((p−1)/2) − ε) mod p²)/p`, and `inverse_series_trunc(p+1)[p]` is `u_p`.

*Validation* (`c6_validate_cw_flint_np.py`, `.out`; `c6b_large_p_crosscheck.py`, `.out`).

* The session's per-prime code path, loaded verbatim, gives the same `q(p)` as the direct recurrence at all 429 primes `3 ≤ p ≤ 3000`. The comparison is of values, not just hits, and the hits are {3, 13}.
* At `p = 100003`, `499979` and `999983`:
  * the vectorised arrays `2^C(k,2)`, `k!` (all `k`) and `1/k!` equal plain loops;
  * `q(p)` equals the independent mod-`p²` formula `(p−1)![x^(p−1)]T′/T²` (94550, 377978 and 555418).

**22. MINOR. The 10^6 search is cited as a result before it has finished** (THM-4750:45 and :362; HYP-9361:14, "extended to 10^6").

* At 21:28 the four jobs (PIDs 85430–85433, started 20:05) had printed only progress lines, and **no hits**.
* Coverage so far: [100001, 418199], [507001, 643061], [711001, 812467] and [867001, 952573].
* Hits are printed immediately with `flush`, so there are no Camion–Wieferich primes in those sub-ranges.
* The remaining sub-ranges were still running (rate about 1 s per prime near `10^6`).
* Until the final "searched" lines appear, cite the search as in progress.

### "What this says" (checklist 7)

**23. MINOR, wording. The "p-adic avatars" paragraph is supported, but it is an interpretation** (THM-4750:453–456).

Each sub-claim is a proved fact:

* `Λ_p` is the class-1 Camion function at `n = 0` and gives the heights: 5(b), 5(c).
* `n ≡ 2 (mod 3)` is the only 2-adic zero class, with a simple zero: 6(a), 6(c) in outline.
* `D(3) = 0` is the only exact zero with `n ≥ 3`, and it forces a split zero when the root is simple and `p ≥ 5`: 4(d), 4(e).

What is interpretive:

* "Avatar" is a reading.
* The K2 identification relates classes mod 3 at `p = 2` to classes mod `p − 1` at odd `p` by analogy, through `F_4^×` replacing `F_p^×`.
* Suggest: "**Reading.** The universal residues … can be seen as …", with "`p ≥ 5`" added to the triangle sentence.

### Minor provenance

**24. MINOR. A comment slightly overstates the code's scope.** The header comment of `camion_lift.py` (lines 2–6) presents its claims as the derivation being tested. The script computes `σ_p` by interpolation, valid for `p ≥ 2K − 1 = 5`, and it was run only for `7 ≤ p ≤ 43`. That is fine, but the comment "Exactly one lift of each residue mod p² is Camion-Wieferich" is a consequence of 5′(c), not something the script checks.

**25. MINOR. The status line needs updating after this pass.** THM-4750:17–19 still calls §5b and the r = 4 data "not yet independently audited".

## 2. Independent checks run

| id | script (in `audit/`) | method | result |
|---|---|---|---|
| c1 | `c1_rigidity_direct.c`, `c1_run.sh` | direct `s_a(p^k) mod p^(k+1)` by `T′ = S′T²`; 7 lifts per residue | 95 + 20 level checks of rigidity, all pass (to mod `3^10`, `5^7`, `7^6`, `11^5`, `13^5`); 168 tower relations OK |
| c2 | `c2_lambda_kappa.py` | direct `Λ_p`, `σ_p(ω)`; `⟨2⟩^(−1/2)` by exponent; designed heights; own 5(e) | `Λ = σ − ⟨2⟩^(−1/2)` at all levels (13 primes); `κ` matches the session; designed heights `h = 1..6` exact; 5(e) full and cheap formulas match, mod `p³` and `p⁴` |
| c3 | `c3_heights_general_residue.py` | as c2, residues 3, 5, `p−1` | all designed-height patterns as predicted (11 pairs `(p, a)`) |
| c3b | `c3b_height_all_bases.py` | every base `a < 5^4`, `7^3` | height `= v_p(a − κ(a))` for all 794 bases |
| c4 | `c4_decomposition_d.py` | direct mod `p²` (two lifts) + FLINT `ρ_C = ε[x^p]1/T_(<p)` + session file | (d) at 301 primes, 0 failures; `ρ_C = 0` at 199, 809, 929 |
| c5 | `c5_cw_lists_flint.py`, `c5_connected_direct.c`, `c5_direct_crosscheck.py` | FLINT mod `p²` series; direct recurrences | base 3 {17, 31}, connected {5, 197} (`p ≤ 10^4`); base 2 {3, 13}; 0 mismatches against direct for `p ≤ 2000` |
| c6 | `c6_validate_cw_flint_np.py`, `c6b_large_p_crosscheck.py` | the session code path replayed; plain-loop arrays | 0 mismatches of `q(p)` at 429 primes; arrays exact near `10^6`; three large `q(p)` confirmed |
| c7 | `c7_hyper_r.c`, `c7_hyper_analyze.py` | Pascal-exponent `2^C(n,r)`, `b′ = T′T^(−1)`, gcd root counts | r = 3 and r = 4 tables exact, 0 per-class differences; forced roots 60/60; low bases obey the law; `p = 71` double root |
| c8 | `c8_two_adic.py` | `Z/2^160[ζ_3]`, Newton polygon, unit root, remainder, quasi-period scan | outline ingredients confirmed; `v_2(r_n) ≥ 3n/2 + 3`; **6(b) fails for `n < 2(j−1)/3`** (P2-M1); exact rule `⟺ v_2(r_n) ≥ j + 2` |
| c9 | `c9_connected_tower.py` | exact integers | connected-graph tower drift `(1 − ε)`; at `p = 5`, 5 and 125 pass and 25 fails |
| c9b | `c9b_connected_rigidity.py` | exact mod `p^(k+1)` | connected `b_a(p^k)` lift-independent (numerical; unproved) |

## 3. Required corrections (summary)

1. **P2-M1 (THM-4750:381–385).**
   * Restate 6(b) with its range: "`f(n + 3·2^j) ≡ f(n) (mod 2^(j+2))` whenever `v_2(r_n) ≥ j + 2`, in particular for all `n ≥ 2(j−1)/3`".
   * In the outline, replace "`(3/2)n − O(1)`" by the bound actually used (observed `v_2(r_n) ≥ 3n/2 + 3`).
   * Replace the exception data by the growing sets (finding 13), with `j = 11`, `n = 5` as the first counterexample beyond `n ≤ 4`.
2. **Finding 20 (THM-4750:364; HYP-9361:71, 77).**
   * Drop "own lift … therefore" for connected graphs, and describe 5, 197 as the primes with `p² | E(p)`.
   * Record the drifted tower `E(p^r) ≡ εE(p^(r−1)) + 1 − ε`, and the fact that for `(2/p) = −1` (both listed primes) `p²` never passes.
3. **Findings 1, 3 (THM-4750:306–331).**
   * Add the existence line, through the tower relation for `s_a`.
   * Optionally replace the limit statement by the finite-level rigidity `s_a(p^k) ≡ s_(a′)(p^k) (mod p^(k+1))`.
   * Add `σ_p ≡ 1 (mod p)`, the definition of `h_a(p)`, the convention for `p`-adic bases, and `c_a(n) = 0` for even `n`.
4. **Findings 6, 7 (THM-4750:338–339).** Replace the tautological `Λ` cross-check by audit c2, and fix "2, 3, 4 … at `p = 7`".
5. **Finding 14 (THM-4750:393).** Add one sentence: `t ∈ r^(Z_2)` because `v(t − 1) ≥ 2` and `r` generates `N^1 ∩ (1 + 4O)`.
6. **Finding 17 (HYP-9361:3, 86).** "inside the analytic range" should become "where the class polynomial already holds (below the proved threshold)".
7. **Finding 22 (THM-4750:45, 362; HYP-9361:14).** Mark the 10^6 search as in progress until the final lines are written.
8. **Findings 8, 23–25.** Fix the printout label; label the "avatars" paragraph a reading and add `p ≥ 5`; adjust the comment in `camion_lift.py`; update the status line after this pass.
9. **Optional OK+ (findings 4, 11).**
   * Record the explicit formula `ρ_C(a) ≡ ε_a [x^p](Σ_(k<p) a^C(k,2)x^k/k!)^(−1) (mod p)` in HYP-9361 (f).
   * Record the telescoping proof of the `mod p^(a+1)` congruence.
