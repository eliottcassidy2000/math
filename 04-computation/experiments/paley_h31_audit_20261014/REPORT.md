# Audit of opus-2026-10-14-S34 (paley_h31_exact_rate_and_apery_wheel_20261014)

Independent adversarial referee: a general-purpose subagent, 2026-10-14, 100 tool uses, about 47 minutes. The harness blocked report files, so the author saved this from the referee's reply. The referee's scripts and outputs are in this directory; no other file was modified by the referee.

**Verdict.** The headline numbers are right.
- `H(P_27)` and `H(P_31)` are reproduced exactly by a genuinely different algorithm.
- THM-4652 holds.
- The `c₁ = 3/2` ingredients are right and now exactly provable.
- The `W_4` identity is true.

The errors are in HYP-9274 (C)/(D), the "where arithmetic enters R(q)" story, and in one false attribution. HYP-9274 (B) was mistyped as open: it is a short corollary of Ahlgren–Ono (2000).

## Errors

- **E1. `A_10` is not a polynomial, so the scoped conjecture in HYP-9274 (C) is refuted.**
  - The referee's exact evaluation (Möbius inversion over set partitions plus the THM-4652 class-`F` calculus) gives `A_10 = 2q(q−1)(21q⁴ − 639q³ + 7028q² − 32098q + 48662) + 34·q(q−1)·a_q`, with `a_q` the coefficients of `η(2z)⁴η(4z)⁴`.
  - An independent parallel DFS gives exactly this value at `q = 11, 19, 23` (478720, 226059264, 1142766592).
  - The polynomial part alone does not vanish at `q = 3, 7`; `A_10` does, because `a_3 = −4` and `a_7 = 24`.
  - `A_8 = 2q(q−1)(q−3)(q−7)(7q−55)` holds, and it matches the June values in THM-438's kind-pasteur addendum.
- **E2. "A single run cannot glue into `W_4`" is false.** The odd-degree argument applies to quotients only. A 10-step trail glues into `W_4` plus a length-2 path between non-adjacent rim vertices, a graph with two odd vertices. Summing out the middle vertex with `χ*χ = 1 − qδ` leaves the `W_4` core. The octahedron is not the first single-run candidate.
- **E3. The leading `W_4` term comes from 2-, 3- and 4-run configurations, with coefficient 7, not about 1.**
  - The `W_4`-core coefficients in `Σ*` (canonical orientation, with `Σ(W_4) = −q(q−1)(a_q + q)`) are `(2,6): −14`, `(4,4): −6`, `(2,2,4): +128`, `(2,2,2,2): −960`.
  - With placement weights, `R/e ∋ −7·a_q/q⁶` (two runs `+17`, three runs `−64`, four runs `+40`). So `1 − R/e ∋ +7a_q/q⁶`, bounded by `14q^(−4.5)`.
- **E4. `η(2z)⁴η(4z)⁴` is not the form of Apéry's irrationality proof.** Beukers' modular proof of `ζ(3) ∉ Q` (Astérisque 147–148, 1987) works at level 6. The level-8 form is that of Beukers' supercongruence `A((p−1)/2) ≡ a(p) (mod p²)` (J. Number Theory 25, 1987), proved by Ahlgren–Ono.

## Imprecise

1. **HYP-9274 (B) is a theorem.**
   - For distinct nonzero `a, b`, scaling gives `W = (q−1)[Σ_(λ≠0,1) χ(λ)g(λ)² + q − 3]`.
   - Ahlgren–Ono Lemma 2.2 (also Pal–Roy–Sadek, Funct. Approx. 69 (2023), Theorem 4.5, `k = 2`) gives `q³·₄F₃(1) = Σ_(λ≠0) χ(λ)g(λ)²`.
   - Ahlgren–Ono Theorem 6 gives `q³·₄F₃(1) = −a_q − q`.
   - So `W/(q−1) = −a_q − 4` for every odd prime. Only the Paley-pattern reading is new.
   - For prime powers the Frobenius trace replaces `a_q`: over `GF(27)`, `W/26 = −264`, against `−156` from `a_27`.
2. **Universality is overclaimed.**
   - The `q = 27` offset (`+4.8·10⁻⁵`) is ten times the `q = 31` holdout error, while fit orders 2 and 3 agree at 27 to `8·10⁻⁷`. It is a resolved deviation.
   - The arithmetic term `−7e·t_q/q⁶` is of the same size as the fit precision.
   - `c₁ … c₄` are universal structurally; `R(q)` is not.
3. **`A_4`, `A_6` and the two-cherry formula are PROVED** by the referee's symbolic engine (validated against brute force). Bonus: `A_(2,4) = 2q(q−1)(q−3)(q−7)(q−13)`.
4. **The pre-registered intervals were rounded inward.** Unrounded: `[2.6279778, 2.6280293]`, `[2.6353086, 2.6353544]`, and model-agnostic `R(47)` `[2.6346761, 2.6359440]`. They should be rounded outward.
5. **The holdout figure picked the better of two variants** (the `c₂`-only variant gives `3.5·10⁻⁵`). The `c₂` uncertainty had no stated procedure; the fits span `−3.04 … −3.15`.
6. **THM-438 addendum.** It needs the "first-order, not rigorous" caveat, and the label "ADDENDUM-9" collides with an existing addendum.
7. **The main note said `c₁ = 3/2` was "proved here"**; it is derived at first order.
8. **The June `A_8` values** in THM-438's kind-pasteur addendum were overlooked.
9. **THE-ZOO** still listed HYP-2371 as un-run.
10. **THM-4652 part 4, `K_5 − e`:** "reverses one edge" holds only for consecutive labels. The parity is still odd, and the Vandermonde sentence is the clean proof.
11. **The mod-256 cross-checks** use the same algorithm, so they are not independent (moot now).
12. **`H(P_43)`** takes about 2 days with full affine canonical forms (orbits of 903), not "a new algorithm is needed".
13. **Missing caveat and citation:** Graham–Ringrose in the transitive-run docstring; cite Evans–Pulham–Sheehan (1981) for `K_4 = 4a² − 2q`.

## Verified (methods)

- **A. The `H` values.** Karp inclusion–exclusion with `Z_31` rotation representatives (69,273,666 classes, 397 s), three 58-bit moduli with one held back as a check. GF(27) was built from a different model (`y³ + 2y² + 1`).
  - `H(P_31) = 19876729650849384218015865` and `H(P_27) = 418425945688463432031`.
  - Acceptance tests: `H` odd; `465 | H`; `N_1 = 42745655163116955307561`, which is odd, `≡ 1 (mod 8)` and `≡ 41 (mod 256)`; `H ≡ 465 (mod 930)`.
  - `R(31) = 2.595511945211` lies inside the interval, `−0.0004781` from the centre.
  - The method also reproduces `q = 7 … 23`.
- **B. THM-4652.** All proof steps check. 700 random multigraphs at `q = 7, 11, 19` gave 0 failures, and 269 series-parallel graphs matched brute force. `W_4` is the unique non-series-parallel survivor among the 34 graph classes on 5 vertices.
- **C. `c₁`.** All four ingredients are correct, and there is no other order-`1/q` term. The fits reproduce at 50 digits.
- **D. The `W_4` identity.** Brute force at all primes up to 67 and at 113, 127, 131; the `O(q²)` reduction for all primes below 1000.
  - `a_q` computed by the pentagonal-number method, cross-checked through Beukers' supercongruence for `17 <= p <= 1005`; Ahlgren–Ono normalisation checked for `p < 100`.
  - `f` is LMFDB 8.4.a.a.
  - No prior statement of the Paley/`W_4` link was found.
- **E. Transitive runs.** The full table (32 primes up to 283) reproduced by a different exact algorithm. `n(p)`, the initial-run property, `trans <= √(4p+1)` and the citations all check.
- **F. Novelty and documentation.** OEIS has no match for the `H(P_q)` or orbit-count sequences. A038375 is known exactly to `n = 11`. The `K_4` claims hold for 167 primes below 1000. The backlog and OPEN-QUESTIONS figures are accurate; the Karp check now serves as design B.

## Files

`audit_H_karp_ie.py` (+ `_small.out`, `_27_31.out`), `audit_thm4652.py`, `audit_run_polys_symbolic.py`, `audit_A10_cores.py`, `audit_A10_dfs_parallel.py`, `audit_W4_channel_multirun.py`, `audit_multirun_check.py`, `audit_W4_apery.py`, `audit_fits.py`, `audit_transitive.py` and `audit_K4_cm.py`, each with its `.out`.

## Disposition (author)

All errors and imprecisions were applied, merged with the concurrent codex and tiling-modular-atoms extensions to HYP-9274. Recorded as MISTAKE-609.
