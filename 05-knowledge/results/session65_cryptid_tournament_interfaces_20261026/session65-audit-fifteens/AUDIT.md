# Session 65 independent audit: Frobenius defects and the ternary fifteens

**Date:** 2026-10-26. **Verdict:** PASS for the elementary pair-orbit theorem, its involution specialization, and the stated arithmetic/ternary applications. One source-numbering correction was requested and fixed before this final pin. No substantive mathematical defect found.

This is an independent audit of the frozen sibling `session65-fifteens`, not a proof of Sloane's or Erdős's conjecture, and not a novelty determination. The repository checkout was read-only; all audit outputs are in this sibling directory. `audit-manifest.json` pins the final input bytes and the independent checker/receipts.

## 1. Final input and correction

- `REPORT.md`: SHA-256 `8af7ff2e90048025eeef63f7072409e43e87f0e9c9d891d688148112ffc5c397`.
- `GENERAL_PERMUTATION.md`: SHA-256 `989cada5edf052e062e75772b284980d0379c28b0aabec77f73fc853591e3bdd`.
- Both source checkers and their normal/optimized receipts were inspected and replayed; replay bytes match the frozen receipts exactly.
- **Corrected:** Saye's v2 PDF labels the lifting result **Lemma 3.1**, not Lemma 1. The report and source manifest now use the PDF numbering. This changed no proof or computation.

Closest mechanism: a signed cycle admits exactly those failure patterns whose parity equals its return twist. Canonical hostile: a swapped two-vertex pair. Corrected near miss: permutation sign detects parity but not minimum defect, witnessed by `(12)(34)`. Least-used indispensable coordinate: the exponent's parity coordinate in the ternary inverse limit.

## 2. Independent proof of the full polynomial

Let a pair-orbit have length `L`, and transport a chosen orientation around it. Write `d_i=1` if the tournament fails to agree with transport across the `i`th comparison and zero otherwise. The product of all comparison signs telescopes: the total number of failures is even when the return fixes both endpoints, and odd when the return swaps them.

Conversely, choose any failure subset with this required parity and either orientation of the initial edge. Propagate around the orbit. The parity condition is exactly what makes the final requirement consistent, so there are exactly two orientations realizing that failure subset. Consequently the factor is

`2 Σ_(j of required parity) C(L,j)z^j = (1+z)^L + ε(1−z)^L`.

Independent pair-orbits use disjoint edge variables. Multiplying establishes the stated formula for every finite labelled vertex set, including the empty set.

For a negative orbit, some power of the vertex permutation swaps its endpoints. The endpoints therefore lie on the same vertex cycle, and that cycle must be even. With cycle length `2a`, only antipodal pairs can be swapped; their unordered-pair orbit has length `a` and is unique in that vertex cycle. Conversely every even vertex cycle supplies this antipodal negative orbit. This proves both directions of the classification without relying on a census.

Thus, with `e` even vertex cycles and `q` unordered-pair orbits,

- the minimum defect is `e`;
- every defect has parity `e`;
- the number of minimizers is `2^q ∏_(even C)(|C|/2)`;
- an invariant tournament exists exactly when the permutation has odd order.

For an involution with `f` fixed vertices and `b` transpositions, the pair-orbits consist of `C(f,2)` positive singletons, `b` negative singletons, and `m=b(f+b−1)` positive doubletons. Their factors give

`2^[C(f,2)+b+m] z^b(1+z²)^m`.

This verifies the sharp minimum, parity, all intervening supported defects, and minimizer count. The conclusion is about the **space of all tournaments**, not an intrinsic arithmetic comparison relation. No ordinary/anti-automorphism identification is used.

## 3. Arithmetic specialization and source hypotheses

Aricheta's Theorems 1.1 and 3.1 assume `p∤N`; their count `Q_p(X)` is the number of non-prime-field supersingular points, not the number of conjugate pairs. At `N=1`, the formula gives `(s_p−r_p)/2=g(X_0^+(p))` with no prime excluded. In particular `p=2,3` are within the cited theorem; the report correctly restricts its separate class-number formulas to `p>3`. The paper also states the Ogg rationality/genus equivalence and lists the Monster prime factors. [Aricheta, v2 PDF, pp. 1–3 and 8](https://arxiv.org/pdf/1809.07421v2).

**Independent finite arithmetic route:** rather than reuse the lane's Legendre-parameter enumeration, the audit constructs a smooth Weierstrass representative for every one of the 1,369 elements `j∈F_37[t]/(t²−2)`, verifies its `j` formula, and evaluates the coefficient of `x^36` in `(x³+Ax+B)^18`. Its zero set is exactly

`{8, 3+10t, 3−10t}`.

Field exponentiation confirms the stated Frobenius transposition. Independent prime-field point counts give only `j=8` with 38 points. The defect distribution on this three-point set is `4z+4z³`, so minimum defect and cited genus both equal one. This supplies an explicit non-Ogg hostile rather than merely inferring one from a class-number table.

The zero-defect interpretation is consequently valid, but does not supply a canonical preferred tournament, a new classification of the Ogg primes, a proof of why there are fifteen, or a connection from that set to Sloane's exceptional exponent. The full unbounded Ogg classification remains **CITED**, not established by the finite census through 211.

## 4. Ternary statements and the inverse-limit trap

The source distinction is correct. Sloane's exceptional exponents are `0,1,2,3,4,15` for omission of digit zero; Erdős's are `0,2,8` for omission of digit two. Saye reports verification through exponent `2·3^45`; that computation was **not** rerun here. The exact primitive-root/lifting input appears as Lemma 3.1. [Saye, v2 PDF](https://arxiv.org/pdf/2202.13256v2). Sloane's own sequence still marks the no-further-exceptions assertion open. [OEIS A102483](https://oeis.org/A102483).

The independent direct census through exponent 1,000 recovers all three displayed digit-omission lists, including `1,3,9` for omission of one. Exact arithmetic confirms

`(222222222222222)_3 → 2^15=(1122221122)_3 → 64=(2101)_3 → 0`.

There are **six**, not five, digits equal to two in the middle power. The 1973 scanned page was separately viewed; the report's decision not to rely on its displayed intermediate is sound. Its persistence motivation and one-way consequence are correctly distinguished from equivalence. [Sloane, 1973, p. 98](https://neilsloane.com/doc/persistence.2.jpg).

For every `k≥1`, let `u_k=2·3^(k−1)` and `n_k=u_k/2−1`. Since 2 is a primitive root modulo `3^k`, `2^(u_k/2)=−1`. Dividing by two gives

`2^n_k ≡ (3^k−1)/2 (mod 3^k)`.

Also `n_(k+1)−n_k=u_k`. Thus the all-ones tail survives coherently at every depth, and adding positive multiples of `u_k` gives arbitrarily large exponents at each **fixed** depth. This does not give one integer exponent realizing the whole infinite branch. Under

`lim Z/(2·3^(k−1)) ≅ Z/2 × Z_3`,

the branch is `(0,−1)`. An ordinary integer with 3-adic coordinate `−1` would be the ordinary integer `−1`, whose parity is one. Hence no ordinary integer realizes the branch. This verifies the report's stronger exclusion of *all* ordinary integers, not merely nonnegative ones.

The signed-tournament use agrees with the directly read PROVED dependency `THM-4754-three-adic-sierpinski-tournament-collatz-arc-flip-steiner-midpoint.md`: multiplication by two reverses the first-nonzero-digit orientation, while multiplication by four preserves it. The `01→11→21` kernel cycle shows that the zero-avoidance predicate is not invariant under the kernel. A shared residue carrier therefore does not transport complete-digit avoidance as an invariant tournament property.

## 5. Machine citations and explicit exclusions

Stérin–Woods's Theorems 1 and 2 give the historical 15-state/2-symbol and 5-state/4-symbol Erdős-counterexample machines. These are **CITED** correctness results; this audit did not reprove their compiler or execute a downloaded simulator. [Stérin–Woods, v2 PDF](https://arxiv.org/pdf/2107.12475v2).

The current community registry lists a 14-state Erdős machine and a 6-state/4-symbol Sloane machine as not independently verified. The audited report preserves that qualification; it does not turn the registry's claim into a theorem. Its parser-only instruction count is correctly distinguished from semantic correctness. [Cryptids, revision 9619](https://wiki.bbchallenge.org/w/index.php?title=Cryptids&oldid=9619).

The report's other-fifteens table is treated as inherited or cited background, not as a fresh independent audit of Mazur, pentagonal tilings, or the complete Apéry-like catalogue. No equality of counts is promoted to an arithmetic map.

## 6. Executed verification and final stopping boundary

`independent_check.py` is new audit code using only the standard library. It imports no lane implementation. Its defect census uses bitstring pushforward and Hamming distance rather than the lane's edge-by-edge comparison loop. Its arithmetic path uses general Weierstrass Hasse coefficients rather than Legendre roots.

Executed successfully:

```sh
python3 independent_check.py > receipt.json
python3 -O independent_check.py > receipt-optimized.json
cmp receipt.json receipt-optimized.json
```

Finite universes:

- every vertex permutation and every tournament for `n=0,…,5`;
- every tournament for each of the 11 vertex-cycle types at `n=6`;
- **484,918** tournament/permutation pairs in the combined exhaustive orientation census;
- all **5,760** vertex permutations at `n=6,7` for signed-orbit structure, factorization, parity, minimum, minimizer count, and the involution specialization;
- every `j` in `F_37²` for the independent Hasse test and every `j∈F_37` for point counts;
- complete digit-omission census through exponent 1,000;
- every exponent residue and all three lifts at depths `1,…,10`;
- exact all-ones branch congruences through depth 100.

The two frozen lane checkers were also replayed in normal and optimized modes, with byte-identical outputs matching their submitted receipts. All checkers use explicit raising conditions rather than removable assertions. No subprocess, shell, network call, or external downloaded implementation occurs inside the independent checker.

**Final boundary:** all assigned proof, source-hypothesis, pinning, and bounded verification obligations are satisfied. There is no execution/access blocker. Sloane, Erdős, recent-machine correctness, conceptual explanation of the number fifteen, and any proposed map between the Ogg set and ternary exceptional exponents remain outside the proved claims.
