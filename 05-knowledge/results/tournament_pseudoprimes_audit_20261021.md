# Audit: THM-4746, THM-4747, HYP-9357, HYP-9358 and the tournament-pseudoprime note (2026-10-21)

**Auditor.** Independent adversarial audit, 2026-10-21, worktree `math-wt-fanopsp`.
* The five files under audit were read in full and **not modified**. Nothing was committed or pushed.
* All verification code is new and lives in `04-computation/experiments/tournament_pseudoprimes_20261021/audit/`. No session function was reused; session outputs were only *compared* against.

**Files audited.**
* `01-canon/theorems/THM-4746-camion-pseudoprimes-…-the-triangle.md`
* `01-canon/theorems/THM-4747-hamiltonian-cycles-under-automorphisms-…-paley-tournaments.md`
* `05-knowledge/hypotheses/HYP-9357-paley-hamiltonian-path-maximality-…-3-7-11.md`
* `05-knowledge/hypotheses/HYP-9358-camion-pseudoprimes-open-laws-….md`
* `05-knowledge/results/tournament_pseudoprimes_camion_fermat_20261021.md` (the note)

## 1. Summary verdict

**THM-4746 statements 1–6 are correct.**
* The group-ring proof of Theorem 1 holds at every step. It was checked by hand and replayed numerically over `F_p` and `F_(p²)`.
* The rotation proof of Corollary 2, the mod-3 identity, the 2-part argument and the Korselt-type criterion are also correct.
* Every exact number in statements 1–6 was reproduced.

**THM-4746 contains three MAJOR errors.**
1. **Connected graphs are quasi-periodic.** The "hostile control" claim in statement 7 is false, and so is the thesis "the mechanism is the linear order of the strong components".
   * Connected labelled graphs satisfy `cg(m+p−1) ≡ (2/p)·cg(m) (mod p)` for **every m ≥ 2**.
   * The auditor proves this in §6, Lemma C. It is also the `k = 1` case of a published theorem (Mani–Stones 2016).
   * The session's test failed only because it demanded the law already at `m = 1`.
2. **The table of `2^k − 1`, `k ≤ 64`, is wrong.**
   * `2^6 − 1 = 63` and `2^12 − 1 = 4095` *are* Camion pseudoprimes (direct computation, three code paths).
   * `k = 59` is decidable: it fails at 179951.
   * The 14 non-squarefree `k` were never decided, although the file says every `k ≠ 59` was.
3. **`D` is quasi-periodic modulo `p²`.** "`D` is not periodic modulo `p²`" is false beyond `m ≤ 2`.
   * The auditor proves `s(m + p(p−1)) ≡ (2/p)·s(m) (mod p²)` for all `m ≥ 3` (§6, Lemma B).
   * This **proves HYP-9358(d)** (the mod-9 law).
   * It also proves the parenthetical equivalence "`p²` passes ⟺ `p² | D(p)`".
   * It extends the criterion to odd parts that are cube-free; this was checked against direct computation on 346 cases.

**The other four files.**
* **THM-4747** is correct. All data were reproduced exactly. Its one error is the class count: 262, not 264.
* **HYP-9357.** The data are correct, but:
  * the Kantor citation points to the wrong paper;
  * the "PROVED given Kantor" deduction has a fixable gap;
  * it misses Alon's and Wormald's prior conjecture that the cyclic (rotational) tournament is the asymptotic maximiser.
* **HYP-9358.**
  * (a): the evidence stands, and its "equivalently" is now proved.
  * (b) and (c) were reproduced (one list is incomplete).
  * (d) is provable.
  * The hostile control in (e) is wrong.
* **The note.** All census claims (§6) were reproduced with independent code. It shares errors 1–3 and has one further false heuristic ("climbs towards e").

## 2. Confirmed

Script names refer to `audit/`. "Direct" means `s(n) mod n` from the binomial recurrence reduced mod `n`, with no use of the criterion.

| claim | how it was recomputed | result |
|---|---|---|
| `s(1..16)` (A054946), `D(r)` for `r ≤ 16` with all factorisations (`D(11) = 2⁵·3²·5·11·2226142198097`, `13² ‖ D(13)`, `3511 \| D(15)`) | `a1`: exact integers via the Hurwitz inverse of `2^C(m,2)`; sympy; brute-force enumeration of all labelled tournaments for `n ≤ 6` (strong counts, circulant counts, Camion's theorem) | all equal; A054946 data equal (OEIS JSON) |
| `D(1) = D(2) = D(3) = 0`, `D(r) > 0` for `r ≥ 4` | proof checked (non-regular strong tournaments exist for `r ≥ 4`); exact `D(r)` for `r ≤ 200` | confirmed |
| Theorem 1, the proof | step by step (see the paragraph after this table) | correct, including `p = 3` and the `m = 0` term |
| Theorem 1, numerics | `a2 period`: every odd prime `p ≤ 400`, `1 ≤ m ≤ 4p` (55,540 checks, `s` and `D`); `a21`: `p = 43, 127, 2731, 8191`, `m ≤ 3p` | 0 failures; `s(p−1) ≢ 0` in general, so `m ≥ 1` is needed |
| the group-ring mechanism itself | `a3`: for every odd `p ≤ 31` (`K = F_p` or `F_(p²)`), solved `τ = ι(t)`, checked `aug t = 1`, `(z−1)^p = 0`, inverted `t`, and checked `s(m) = −Σ_(a≠0) d_a (g₀a)^m` for `1 ≤ m ≤ 4p+2` | all pass; `Σ_x d_x = 1` at `m = 0` with `d_0 ≠ 0` |
| Corollary 2, both proofs | `a19`: `t_k = 2^(k(kp−1)/2) ≡ ε^k 2^C(k,2)`, `F_k ≡ ε^k s(k)`, `s(kp) ≡ F_k (mod p)` for `p ≤ 13`, `k ≤ 8`; brute force on `Z_6` (32 invariant tournaments, 24 strong) | confirmed; strong components of fixed tournaments are unions of cosets (the condensation has no nontrivial automorphism) |
| statement 3 | identity re-derived; `a8`: `s(m) mod 3` and `mod 5` for `m ≤ 200` | confirmed |
| Theorem 4 (criterion, including the 2-part) | `a2 direct`: **every** composite `4 ≤ n ≤ 2500` (2,132 numbers); `a4`: criterion from exact `D(r)` | agrees on all 1,321 `n ≤ 2000` with squarefree odd part (508 odd, 813 even); no passing non-squarefree `n` violates the necessary condition; session counts 472 / 56 / 142 reproduced (`a4b`) |
| all composite `n ≤ 600` passing | `a4_korselt_check.out` | list printed |
| odd composite Camion pseudoprimes `≤ 1200` | `a4` | identical to `camion_psp_1200.out` (92 numbers) |
| `21, 1023, 561, 1105, 1729` pass; `341, 85, 25, 27` fail | direct | confirmed (`341`: `r_31 = 11`, `D(11) ≡ 22 (mod 31)`) |
| families | all 9 Carmichael numbers `≤ 20000` direct; every `2p` (167) and `3p` (119) `≤ 2000` direct; `3b(3b−2)` members `195, 399, 1023, 1443, 4623, 12099, 16383, 19599` direct; `67108863 = 3·2731·8191` via the proved criterion (`r = 1, 3, 3`) plus Theorem 1 numerics at 2731 and 8191 | confirmed |
| squarefree list to 20000 | `a13`: criterion with exact `D` (`r_p ≤ 137` always) | 2145 members, identical to `camion_squarefree_list.txt`; `2p`: 1228, `3p`: 857, exotic: 60 (37 odd, 23 even) |
| Camion–Wieferich primes `≤ 3000` | `a2 wief`: `s(p) mod p²` for all 429 odd primes | exactly `{3, 13}` |
| prime squares | direct | `9, 169` pass; `25, 27, 49, 81, 121, 125` fail |
| mod-9 pattern | `a8`: `s(m) mod 9`, odd `m ≤ 1401` | holds to 1401 (and is proved in §6) |
| zero classes and the Poisson table (HYP-9358(b)) | `a14`: EGF-inversion over `F_p` (a different algorithm), `p ≤ 400`, all `r` | odd classes identical to `camion_verify.out`; Poisson table reproduced exactly (31, 25, 13, 4, 2, 0; mean 0.9467) |
| statement 6 | `a9`: equivalence with `ord_n(2) \| (n−1)²/2` for every odd composite `n ≤ 20000` | 0 mismatches; `85, 133, 481, 511, 949` pass the count test and are not base-2 pseudoprimes; `645` is one; the 21/85/341/1023/561/1105/1729 table is correct |
| Burnside reading | `a9`: brute-force orbit counts for `m = 3, 5, 7` | `(4, 2), (208, 112), (299600, 239648)` = the formulas |
| statement 7, semicomplete part | `a10`: brute force `n ≤ 5` gives `1, 1, 15, 543, 51969` = the `3^C(m,2)` recurrence (linear condensation); `(3/p)` law for `5 ≤ p ≤ 61` | confirmed; for `p = 3` it holds trivially with `(3/3) = 0` |
| THM-4747 proofs (fixed-vertex lemma; orbit analysis; reflections are the only path-fixing anti-automorphisms; Burnside; parity) | by hand | correct |
| THM-4747 data | `a11` (subset DP; full-DP cross-check `n ≤ 17`) and `a12` (all permutations, **every** connection set `n ≤ 11`, 62 sets, 0 mismatches) | 262 multiplier classes for odd `3 ≤ n ≤ 23` (1, 1, 2, 4, 4, 6, 16, 16, 30, 88, 94); no congruence violation; naive test 0/4, 0/16, 0/88; multiset of `(H, C, ρ_0, μ, A)` identical to the session's |
| Paley data | `a11` | `P_7: (189, 24, 9, 9)`, `P_11: (95095, 5505, 1729, 185)`, `P_19: (1172695746915, 34358763933, 6857869865, 573057)`, `P_23: (15760206976379349, 374127973957716, 62293308207033, 63871533)`; `t = 1, 100, 200928444, 1478766695485`; `tc mod q = 1, 2, 2, 4, 14` |
| HYP-9357 table, `n ≤ 23` | `a18` | every row reproduced (`Z_9{2,3,4,8}`: 3357; rotational is the maximum at 13: 3711175, 15: 198464295, 17, 19, 21, 23); A038375 = `1,1,3,5,15,45,189,661,3357,15745,95095` (OEIS) |
| designs | `a20` (nauty `dreadnaut`) | `Aut(D(7)) = 168`, `Aut(D(11)) = 660`, both 2-transitive; not 2-transitive for `q = 19, 23, 27, 31, 43, 47`; multiplier groups: Singer `⟨2⟩` of order `m = 3, 5, 7` vs QR of order `(q−1)/2`; Fano point stabiliser has order 24 with the element-order profile of `S_4` |
| census `n ≤ 9` | `a15`/`a16`: own C code from `gentourng -q n`; `\|Aut\|` = number of Hamiltonian paths `Q` with `P_0[i] ↦ Q[i]` an automorphism; `Σ n!/\|Aut\| = 2^C(n,2)` for every `n` | all claims hold: `tc` misses exactly `{7, 21}` below 60 at `n = 7, 8, 9`; `H` misses exactly `7, 21` below the ceiling (611 at `n = 8`, 2883 at `n = 9`); `C` has no gaps at `n = 8`; at `n = 7` it misses 16 and 18–23, with maximum 24 attained only by `P_7`; at `n = 9`, 21 values in `[215, 248]`; strong holes 47, 55, 67 at `n = 8` and 77, 79, 83 at `n = 9`; the session's `census_n*.txt` multisets are identical |

**Theorem 1, step-by-step verdict.**
* `ι` is a unital ring homomorphism: `E_x ∗ E_y = E_(x+y)` and `E_0 = e`.
* `T mod p` satisfies `T(m+p−1) = ε T(m)` for all `m ≥ 0`, including `m = 0`.
* The `p−1` exponentials `E_(g₀a)` lie in the `(p−1)`-dimensional solution space, and their Vandermonde matrix on nodes `g₀a` is invertible, so they span it.
* `aug(t) = τ_0 = 1`.
* In characteristic `p`, `([x]−1)^p = [px]−1 = 0`, so the augmentation ideal is generated by nilpotents and is nilpotent. Hence `t` is a unit.
* Inverses in the commutative ring `H(K)` are unique, and `(e−σ)∗τ = e` is Wright's recurrence.
* The `x = 0` term matters only at `m = 0`. The audit found `d_0 ≠ 0`, and `s(p−1) mod p ≠ 0` (e.g. `s(4) ≡ 4 (mod 5)`), so the restriction `m ≥ 1` is genuinely needed.
* For `p = 3` the argument works over `F_9` with `g₀ = i`.

## 3. Errors

### MAJOR

**E1. Connected labelled graphs are quasi-periodic. The "hostile control" and the "linear order is the mechanism" thesis are false.**

*Where.*
* THM-4746 title (7), "connected graphs (unordered components) have none".
* THM-4746 status ("PROVED: statements 1-7"; "the controls of statement 7").
* THM-4746 §1, lines 77–81.
* Note §1 (line 67) and §2 (lines 88–90: "Linearity of the condensation is essential … connected graphs obey none").
* HYP-9358(e) (line 50: "Hostile control: connected graphs … have no quasi-period").

*Problem.* For every odd prime `p` and every `m ≥ 2`, `cg(m+p−1) ≡ (2/p)·cg(m) (mod p)`, where `cg` counts connected labelled graphs (A001187). This is the same period and the same multiplier as for strong tournaments.
* Proof: §6, Lemma C (shift operator; three lines).
* Numerics (`a10b`): every odd `5 ≤ p ≤ 101` and `1 ≤ m ≤ 8p`. The **only** failure is at `m = 1`, where `cg(p) ≡ (2/p) − 1` while `(2/p)·cg(1) = (2/p)`.
* `camion_controls.py` tested (anti)periods only *from `m = 1`*, so this single index hid the law.
* The finite sentence in §1 ("no period or antiperiod `L < 2p` … for any `p ≤ 23`") is literally true. The conclusions drawn from it are false.
* Published prior art: Mani–Stones (2016), Proposition 1.11 as quoted in Kotek–Makowsky (arXiv:2112.06581). For `p ≥ 3`, `b ≢ 1 (mod p)`, `n ≥ p^k` and `n > p`, `T(K_n;1,b)` is congruent mod `p^k` to `b^(φ(p^k)/2)` times the term at `n − φ(p^k)`.
  * `T(K_n;1,2)` is the number of connected labelled graphs. With `b = 2` and `k = 1` the proposition is exactly the law above.
  * The HTML rendering of that formula was partly garbled. The `k = 1` statement therefore rests on the auditor's own proof (Lemma C) and numerics.
  * The `mod p²` analogue was also seen in `a17`: it fails only at `m ≤ 2`.

*Fix.*
* Delete "connected graphs (unordered components) have none" from the title and status.
* Rewrite §1 and note §1–§2 as follows. The mechanism is that `T mod p` lies in the image `ι(K[G])` of a finite local group algebra. That image is closed under Hurwitz inversion (ordered decompositions) and under the shift `Δ` (rooted or unordered decompositions, via `Δ(cg) = ΔT ∗ T^(−1)`).
* Linearity of the condensation is **not** essential. It only fixes the threshold (`m ≥ 1` versus `m ≥ 2`).
* HYP-9358(e) needs a new hostile control. Rooted counts such as `m·s(m)` are a real negative for period `p−1`; they have quasi-period `p(p−1)`.

**E2. The `2^k − 1` (`k ≤ 64`) table is false.**

*Where.* THM-4746 status ("every 2^k − 1 with k ≤ 64 except k = 59 (undecided), the passes being exactly k = 4, 10, 14, 26, 34, 38, 62") and §5, line 139. The note's status row "Mersenne `k ≤ 64` | FINITE-EXACT" inherits the error.

*Problem.*
* `mersenne_camion.py` skips every `2^k − 1` whose odd part is not squarefree. Its own `.out` labels them `squareful`. These are `k = 6, 12, 18, 20, 21, 24, 30, 36, 40, 42, 48, 54, 60, 63`. They were never decided, yet the THM says every `k ≠ 59` was.
* Two of them pass. Direct `s(n) mod n` gives:
  * `63`: `s(63) mod 63 = 2 = c(63)` (exact integer `s(63)` and the C code agree);
  * `4095`: `s(4095) mod 4095 = 128 = c(4095)` (C code and an independent pure-Python forward recurrence agree).
* The other twelve fail. Their first failing primes are 19, 11, 127, 17, 151, 19, 11, 43, 17, 19, 11, 73, respectively. `2^18 − 1` also fails directly.
* `k = 59` is not undecided. `r_179951 = 29737`, and `D(29737) ≡ 46737 (mod 179951)` by two different algorithms: the binomial recurrence (`a2`) and EGF inversion over `F_p` (`a6`). So `2^59 − 1` fails at 179951. The session's script simply capped `r ≤ 3000`.

*Fix.* The composite `2^k − 1` with `k ≤ 64` that are Camion pseudoprimes are exactly `k ∈ {4, 6, 10, 12, 14, 26, 34, 38, 62}`. None is undecided (`a5_mersenne.out`). Note also:
* `63 = 4³ − 1` is the `q = 3` member of the `4^q − 1` family. It is excluded by `b ≥ 5` only because `9 | 63`.
* By Lemma B the 9-part never obstructs, so `63` passes for a structural reason.

**E3. "`D` is not periodic modulo `p²`" is false beyond `m ≤ 2`; HYP-9358(d) is provable.**

*Where.* THM-4746 §5, line 147; note §2, line 119; HYP-9358(d), line 44 ("not periodic modulo `p²` with any period dividing `2p(p−1)`. So a proof needs a mod-`p²` refinement…").

*Problem.* For every odd prime `p` and all `m ≥ 3`: `s(m+p(p−1)) ≡ (2/p)·s(m)` and `D(m+p(p−1)) ≡ (2/p)·D(m) (mod p²)`. Hence `D(m+2p(p−1)) ≡ D(m) (mod p²)` for `m ≥ 3`.
* Proof: §6, Lemma B. This *is* the mod-`p²` refinement HYP-9358(d) asks for.
* Numerics (`a17`, `a17b`): `p ≤ 31`, `m ≤ 3p(p−1)+10`. Failures occur only at `m ∈ {1, 2}`; the threshold `m ≥ 3` is sharp.
* No committed session script tests periodicity modulo `p²`. The claim is unsupported as well as false.

*Consequences (all checked).*
1. **HYP-9358(d) is a theorem.** For `p = 3`: `D(m+6) ≡ −D(m) (mod 9)` for `m ≥ 3`, with `D(3) = 0`, `D(5) = 540 ≡ 0`, and `D(7) ≡ 6 (mod 9)`. So for odd `m`, `9 | D(m)` iff `m ≢ 1 (mod 6)` or `m = 1`, exactly as conjectured. Also `9 | D(m)` for every `m ≥ 3` with `3 | m`, so the factor 9 never obstructs (`a22`: 170/170 cases `≤ 2500`).
2. **The Camion–Wieferich equivalence is proved.** Take `m = p`: `D(p²) ≡ (2/p)·D(p) (mod p²)`, so `p²` passes iff `p² | D(p)`. Directly verified for `p ≤ 97` (`a7`).
3. **The criterion extends to cube-free odd parts.** If `p² ‖ n`, require `p² | D(r)` with `r ≡ n (mod p(p−1))` and `3 ≤ r ≤ p(p−1)+2`. This agrees with direct computation on all 346 composite `n ≤ 2500` whose odd part is cube-free but not squarefree (`a23`).

*Fix.* Replace the sentence by Lemma B. Move HYP-9358(d) to PROVED, with Lemma B as the proof. Optionally record the observed mod-`p^k` analogue: `s(m+φ(p^k)) ≡ (2/p)·s(m) (mod p^k)` holds for `m ≥ 2k−1` in all tested cases (`k = 3`: `p ≤ 7`; `k = 4`: `p = 3`). It is conjectural for `k ≥ 3`; the same proof should work, with derivatives up to order `2(k−1)`.

### MINOR

* **E4. Class count.**
  * *Where:* THM-4747 title (5), §5 line 99, and note §3 line 146 ("264 circulant classes").
  * *Problem:* There are 262 multiplier classes (1, 1, 2, 4, 4, 6, 16, 16, 30, 88, 94). `check_congruences.py` reads `paley19.out` and `paley23.out` *in addition to* `circ_19.out` and `circ_23.out`, which already contain `P_19` and `P_23`. So 264 = 262 + 2 duplicates.
  * *Fix:* say 262.
* **E5. Wrong Kantor citation.**
  * *Where:* HYP-9357 status.
  * *Problem:* W. M. Kantor, *J. Combin. Theory A* 38 (1985) 66–74 is "Homogeneous designs and geometric lattices". It classifies 2-transitive designs with `λ = 1` and does not contain the `(11,5,2)` biplane.
  * *Fix:* cite W. M. Kantor, "Classification of 2-transitive symmetric designs", *Graphs Combin.* 1 (1985) 165–166, DOI 10.1007/BF02582940.
* **E6. Gap in the "PROVED (elementary, given Kantor)" deduction (HYP-9357).**
  * *Problem 1:* Different multiplier groups show only that the difference sets `QR_q` and the Singer set are *inequivalent*. 2-transitivity needs the *designs* to be non-isomorphic.
  * *Bridge:* for prime `q = 2^m − 1`, `q ‖ |PGL(m,2)|`. So any regular cyclic subgroup is a Sylow `q`-subgroup, and they are all conjugate. A design isomorphism can therefore be chosen to intertwine the translation groups, which makes it affine and hence an equivalence of difference sets.
  * *Problem 2:* Kantor's list is up to complementation. The complement case reduces to the same `PG(m−1,2)` case.
  * *Problem 3:* `2^m − 1` is never a proper prime power, so only prime `q` arise.
  * *Problem 4:* "the multiplier group of the Singer set is `⟨2⟩`, of order `m`" needs a citation.
  * *Auditor's check:* the conclusion is right. `D(31)` vs `PG(4,2)` and `D(127)` vs `PG(6,2)` have different block-triple intersection profiles, so they are non-isomorphic (`a20`). Nauty confirms that `D(q)` is not 2-transitive for `q ∈ {19, 23, 27, 31, 43, 47}`.
  * *Fix:* add the four points above.
* **E7. Unproved equivalence asserted.**
  * *Where:* THM-4746 §5, "(`p²` passes)", and HYP-9358(a), "equivalently".
  * *Problem:* both assert the equivalence without proof; the status calls it finite agreement for `p ≤ 1500`.
  * *Fix:* it is true; cite Lemma B (consequence 2).
* **E8. Pointer to a list that is not there.**
  * *Where:* THM-4746 §5, line 152: "The ones not of the form `3p` up to 20000 are listed in the results note".
  * *Problem:* they are not. The note lists only family members; the full list is only in `camion_squarefree_list.txt`.
  * *Fix:* point to the `.txt` file or add the list.
* **E9. Incomplete sporadic list.**
  * *Where:* HYP-9358(c), line 39: "products using sporadic classes (1239, 3765, 8401, 12019, 14317 below 20000)".
  * *Problem:* these are the *odd* ones. The even sporadic members `186, 530, 1898, 3290, 7202, 8282, 17066` are omitted (`a13b`).
* **E10. "All 142 even n ≤ 400".**
  * *Where:* THM-4746 status and note §2.
  * *Problem:* there are 163 even `n ≤ 400` with squarefree odd part. The 142 are those whose odd primes are `≤ 100`. The auditor checked all 813 even `n ≤ 2000`.
* **E11. "The moduli 3 and 5 test nothing" (THM-4746 §3 heading).**
  * *Problem:* modulo 5 the test is not empty when `4 | n`: then `r_5 = 4` and `D(4) = 24 ≢ 0`, so e.g. `20, 40, 60` fail at 5 (direct).
  * *Fix:* restrict the heading to odd `n` (and `n ≡ 2 mod 4`).
* **E12. Missing hypotheses.**
  * Title (5): "(all residues in the triangle class 3)". The factor 3 has residue 1.
  * Title (6) needs "odd `n`".
  * Title (7): the stated proof of the `(3/p)` law needs `p ≥ 5`; for `p = 3` it holds trivially with `(3/3) = 0`.
  * §4: "`Z_p` always contains `{1, 2, 3}`" needs `p ≥ 5`, since `Z_3 = {1, 2}`.
  * Note §2 table: "3p … 9, 15, …" lists `9`, which is not squarefree and is not one of the 857.
  * THM-4747 (3): "in particular it is odd" uses that `μ` is odd. This holds because `M_S` is a subgroup of the cyclic group `Z_p^*` not containing `−1`; say so.
* **E13. "H·2^(n−1)/n! climbs towards e" (note §4, line 164) is unsupported and likely wrong.**
  * `e` is only the Adler–Alon–Ross lower bound.
  * Wormald proves `h'(n) > 2.85588·n!/2^(n−1)` for infinitely many `n`, using blow-ups of the *cyclic (rotational)* tournament `R_d`, and conjectures the limit is about 2.855958.

## 4. Overclaims and wording

* **THM-4746 title (6) and §6; note line 123.** "The triangle is exactly what separates the two tests" is an overclaim.
  * Only the sets of **universal** residues differ by exactly `{3}`.
  * The tests also disagree for non-universal reasons. `85`, `341` and `645` pass the count test through count-test zero classes but fail the strong test. `1239`, `3765`, … pass the strong test through sporadic strong classes.
  * Suggested wording: "the universal residue sets of the two tests differ exactly by the triangle class `r = 3`".
* **THM-4746 statement 7 typed PROVED.** The connected-graph half is false (E1). The control results were FINITE-EXACT statements from which a false general conclusion was drawn.
* **"Theorem 1 reads it as Euler's criterion" (§2).** This is acceptable as a reading. Theorem 1's proof uses `2^((p−1)/2) ≡ (2/p)` as input, so no new proof of Euler's criterion is implied.
* **Title (5): "D(r) > 0 for r ≥ 4, so every Carmichael number, every 2p and every 3p…".** The "so" is a non sequitur. Those families follow from `D(1) = D(2) = D(3) = 0`; `D(r) > 0` only shows that no other residue is universal.
* **HYP-9357 title, "maximises … among circulants … iff q is 3, 7 or 11".** `P_27` and other prime-power Paley tournaments are not circulants. Restrict to prime `q`, or say Cayley tournaments.
* **HYP-9358(e).**
  * As stated it is essentially a theorem: the THM-4746 argument applies verbatim whenever `T mod p` is `(p−1)`-quasi-periodic, and Lemma C shows unordered decompositions work too.
  * The candidate "tournaments with a marked vertex" fails for period `p−1` immediately: `m·s(m)` has quasi-period `p(p−1)`.
* **ANALOGY items.** The ANALOGY and DICTIONARY items of note §7 (7/21/341 cyclotomic, `{2,3,11}`, Fano/tetrahedron, sunflowers) are correctly typed. Arithmetic checked:
  * `Φ_3(2) = 7`, `Φ_3(4) = 21 = 7·3`, `Φ_5(4) = 341`;
  * Cipolla residues;
  * Mersenne primes `≡ 7 (mod 8)` and Wagstaff primes `≡ 3 (mod 8)`;
  * `F_7 \ e = M(K_4)` with stabiliser `S_4`;
  * arXiv:2610.26686 is Bucić–He–Lee, "Proof of the Erdős–Szemerédi sunflower conjecture" (the bound `3(2 − 1/(16k³))^n` was not checked).
  * One detail: in `C_3^[k]` every line is a cyclic triangle, but not every cyclic triangle is a line (e.g. `(0,0), (1,0), (2,1)`). The note's wording is consistent with this, but should not be read as a bijection.

## 5. Missing credits and prior art

* **Specker–Blatter (1981; Specker 1988).**
  * The labelled count of any MSO-definable class of structures with one binary relation satisfies linear recurrences modulo every `m`, so it is ultimately periodic mod every `m`.
  * "Strong tournament" is MSO-definable: a tournament axiom plus "no nonempty proper vertex set is closed under out-arcs".
  * So ultimate periodicity of `s(n) mod m` for **every** `m`, including `p²`, was known in principle. Theorem 1 and Lemma B give the explicit periods and multipliers.
  * Statement and extensions are in Kotek–Makowsky, arXiv:2112.06581 (Thm 2.9).
* **A. P. Mani and R. J. Stones**, "The number of labeled connected graphs modulo prime powers", *SIAM J. Discrete Math.* (2016), DOI 10.1137/15M1024615.
  * As quoted in Kotek–Makowsky (Prop. 1.11, Conj. 1.12; the formula was partly garbled in the HTML), they give an explicit law for connected graphs: `C_n` is congruent mod `p^k` to `2^(φ(p^k)/2)` times `C_(n−φ(p^k))`, for `n > p`, `n ≥ p^k`.
  * Its `k = 1` case is re-proved here as Lemma C.
  * This is the closest prior art to Theorem 1, and it contradicts the hostile control (E1).
* **N. Alon**, "The maximum number of Hamiltonian paths in tournaments", *Combinatorica* 10 (1990). **N. C. Wormald**, "Tournaments with many Hamilton cycles" (preprint, users.monash.edu.au/~nwormald/papers/hamtourn.pdf).
  * Alon found it plausible that the cyclic tournament `R_n = Z_n{1..(n−1)/2}` (the "rotational circulant") maximises Hamilton cycles.
  * Wormald conjectures `H(R_n) = h(n)` for large odd `n`, analogously for paths. His blow-ups of `R_d` give the best known constant, about 2.8559.
  * HYP-9357 part 3, and the "only if" direction of part 1 for large `q`, are therefore special cases of an existing conjecture and should be credited.
  * The census data show that the naive all-`n` version fails at small `n`: `C(P_7) = 24 > C(R_7) = 17`, and at `n = 9`, `249 > 222`.
* **Adler–Alon–Ross**, *Random Structures Algorithms* 18 (2001): the lower bound `(e − o(1))·n!/2^(n−1)` (relevant to E13).
* **Wilson and Fermat by rotation (THM-4747 §4).**
  * The un-weighted argument is classical. Sagan, "Congruences via abelian groups", *J. Number Theory* 20 (1985), credits Peterson (1872).
  * See also Anderson–Benjamin–Rouse, *Amer. Math. Monthly* 112 (2005), and Golomb's necklace proof.
  * The tournament weighting cancels by Fermat, so §4 should say "remark" and cite.
* **Gessel**, congruences via group actions on labelled graphs, as credited in Sagan (1985). This is the general method behind Corollary 2's second proof.
* **What was not found.**
  * No prior statement of `s(m+p−1) ≡ (2/p)·s(m)` for strong tournaments.
  * No prior statement of the Korselt-type criterion or the "Camion pseudoprime" notion.
  * The odd Camion pseudoprimes, all Camion pseudoprimes, the exotic list and the count-test list `85, 133, 341, 481, 511, 561, 645, 949, 1105` all return no OEIS match (queried 2026-10-21).
  * A054946 records no congruence.
  * THM-4747's cycle congruences are routine Burnside counting, and no source states them.

## 6. Auditor's proofs supporting the fixes

**Notation.**
* `H(Λ)` is the Hurwitz ring (binomial convolution) over a commutative ring `Λ`, and `E_x = (x^m)_m`.
* `Δ` is the left shift `(a_m) ↦ (a_(m+1))`, so `ΔE_x = x·E_x`.
* `ε = (2/p)`, and `K`, `g₀`, `G = g₀F_p`, `ι`, `t` are as in THM-4746 Theorem 1.

### Lemma C (connected graphs)

**Statement.** For every odd prime `p` and every `m ≥ 2`: `cg(m+p−1) ≡ ε·cg(m) (mod p)`.

*Proof.*
1. Root the recurrence at the last vertex. Its component has `k` further vertices, so `T(m+1) = Σ_k C(m,k)·cg(k+1)·T(m−k)`. In the Hurwitz ring this reads `ΔT = (Δcg) ∗ T`.
2. Modulo `p`, `T = ι(t)` and `ΔT = ι(t̃)` with `t̃ = Σ c_x·x·[x]`. Also `T` is invertible in `ι(K[G])` (THM-4746).
3. Hence `Δcg ≡ ι(t̃·t^(−1)) = Σ_(x∈G) f_x E_x`, i.e. `cg(m+1) ≡ Σ_x f_x x^m`.
4. For `m ≥ 1` the term `x = 0` vanishes, and `x^(p−1) = ε` for `x ≠ 0`. ∎

At `m = 1` the law fails in general, because `cg(p) ≡ ε − 1`.

### Lemma B (modulo `p²`)

**Statement.** For every odd prime `p` and every `m ≥ 3`: `s(m+p(p−1)) ≡ ε·s(m)` and `D(m+p(p−1)) ≡ ε·D(m) (mod p²)`.

*Proof.*

**Setup.**
* Let `R = Z_p[g₀]`, where `g₀ = 1` if `ε = 1`. If `ε = −1`, `g₀` is the Teichmüller lift of an element of order `2(p−1)` of `F_(p²)^*`, so `g₀^(p−1) = ε` exactly.
* Let `Λ = R/p²` and `X = {g₀ζ : ζ ∈ μ_(p−1)} ∪ {0}`, so that `X mod p = G`.
* Put `∂E_x = (m x^(m−1))_m` and `∂²E_x = (m(m−1)x^(m−2))_m`. For `x = 0` these are `δ_1` and `2δ_2`.
* Let `A` be the `Λ`-span of `{E_x, p∂E_x, p∂²E_x : x ∈ X}`.

**(i) `A` is a ring.**
* `E_x∗E_y = E_(x+y)`, and `∂^jE_x ∗ E_y = ∂^jE_(x+y)`.
* Write `x + y = z + pw` with `z ∈ X` and `w ∈ R`. Then `E_(z+pw) ≡ E_z + pw·∂E_z` and `p∂^jE_(z+pw) ≡ p∂^jE_z (mod p²)`.
* Products of two multiples of `p` vanish.

**(ii) `T ∈ A`.**
* Write `2 = ω(1+pz')` with `ω` the Teichmüller lift of 2. Then `2^N ≡ ω^N(1 + Npz') (mod p²)`.
* The sequence `a(m) = ω^C(m,2)` satisfies `a(m+p−1) = ε·a(m)` exactly. The reason is that `ω^((p−1)/2) = ±1 ≡ 2^((p−1)/2)`.
* So `a = Σ_ζ c_ζ E_(g₀ζ)` with `c_ζ ∈ R`; the Vandermonde matrix on Teichmüller nodes has unit determinant.
* Then `C(m,2)·a(m) = Σ c_ζ·(x_ζ²/2)·∂²E_(x_ζ)(m)`, and `T ≡ a + pz'·C(m,2)·a ∈ A`.

**(iii) `T` is a unit of `A`.**
* Reduction mod `p` maps `A` onto `ι(K[G])` with kernel `N`, and `N² = 0`.
* Lift `ι(t^(−1))` to some `Y_0 ∈ A`. Then `T∗Y_0 = e − n` with `n ∈ N`.
* Put `Y = Y_0∗(e+n)`. Then `T∗Y = e`.
* By uniqueness of inverses and Wright's recurrence, `Y = e − S`.

**(iv) Every `f ∈ A` satisfies `f(m+L) ≡ ε·f(m) (mod p²)` for `m ≥ 3`, where `L = p(p−1)`.**
* For `x = g₀ζ` we have `x^L = ε` and `pL ≡ 0 (mod p²)`. The derivative terms keep this property.
* The `x = 0` generators are supported on `m ≤ 2`.

**(v) Conclusion.**
* Apply (iv) to `e − S`.
* For `c`: `c(m+L) = 2^(L/2)·c(m)` for odd `m`, and `2^(p(p−1)/2) ≡ ε^p = ε (mod p²)`. ∎

The failures at `m ∈ {1, 2}` observed in `a17` are exactly the `δ_1` and `δ_2` supports.

## 7. Reproduction

```bash
cd 04-computation/experiments/tournament_pseudoprimes_20261021/audit
sh run_audit.sh        # compiles a2, a6, a11, a12, a14, a15 and runs every script; about 5.5 min on the mac-mini (measured 5:20)
```

Individual steps, each script writing the matching `.out`:

| step | command |
|---|---|
| exact values, brute force `n ≤ 6`, factorisations | `python3 a1_exact_bruteforce.py` |
| Theorem 1 numerics | `./a2_smod period 400 4` |
| group-ring replay | `python3 a3_groupring_check.py` |
| direct `s(n) mod n`, all composite `n ≤ 2500` | `./a2_smod direct 4 2000; ./a2_smod direct 2001 2500` |
| criterion vs direct | `python3 a4_korselt_check.py; python3 a4b_session_counts.py` |
| Mersenne table | `./a2_smod single 262143; python3 a5_mersenne.py` (uses `a2_smod` and exact `D`); independent `k = 59` check: `./a6_egf_prime 29737 179951` |
| Camion–Wieferich | `./a2_smod wief 3000` |
| prime-square identity | `python3 a7_prime_square_identity.py` |
| statement 3 and mod 9 | `./a2_smod seq 1401 45; python3 a8_mod3_mod5_mod9.py` |
| statement 6 and Burnside | `python3 a9_count_test.py` |
| controls (E1) | `python3 a10_controls.py; python3 a10b_connected_extended.py` |
| circulants and THM-4747 | `for n in 3 5 … 23; do ./a11_circ $n all; done`; `./a12_circ_brute n`; `python3 a12_join.py; python3 a18_thm4747_checks.py`; Paley: `./a11_circ 19 paley; ./a11_circ 23 paley` (about 2 s, 250 MB) |
| lists and zero classes | `python3 a13_list20000.py; python3 a13b_exotic_classify.py; ./a14_zero_classes 400; python3 a14_zero_classes_check.py` |
| mod `p²` / `p^k` (E3) | `python3 a17_mod_p2_quasiperiod.py; python3 a17b_mod_pk_quasiperiod.py; python3 a22_nine_criterion.py; python3 a23_square_criterion.py` |
| Corollary 2 rotation proof | `python3 a19_rotation_proof.py` |
| HYP-9357 designs | `python3 a20_designs.py` (needs `dreadnaut`) |
| large primes, `2p`/`3p` direct | `python3 a21_large_primes_period.py` |
| census | `for k in 3 … 9; do gentourng -q $k \| ./a15_census $k > a15_census_n$k.txt; done; python3 a16_census_analyze.py` |

**Literature checks.**
* OEIS JSON (`curl` with a user agent): A054946, A038375, and the four pseudoprime sequences above.
* Web: the Mani–Stones abstract; Kotek–Makowsky arXiv:2112.06581 (HTML v5); Kantor's Springer record; Wormald's preprint PDF (read pages 1–9); the arXiv:2610.26686 abstract.
