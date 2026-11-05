# Audit D — THM-4876 statement 7 (end-count bound `e_a <= (h - 1)/2`), 2026-11-04

Independent adversarial audit of ONE statement: THM-4876 (7), its proof paragraph **7.** and the clause "(7) ..." of the
YAML title, in `01-canon/theorems/THM-4876-*.md` as of 87ac336954 plus the author's two self-corrections (the `sigma`
term and the `TT_(k+2)` relabelling), re-checked unchanged at HEAD ed565d1871 (audit C applied meanwhile). Nothing in
canon was edited. Scripts and
receipts: `04-computation/experiments/chanpak_spine_20261104/audit_D/` (written from scratch; the session's scripts were
read, not imported).

## Verdicts

| Object | Verdict |
|---|---|
| Statement 7 (bound, duality, equality family, consequences, order-9 census claim) | **PASS WITH CORRECTIONS** (all minor: wording, one convention, one receipt) |
| Proof 7 | **PASS** (correct as written; one step redundant, see D7) |
| Title clause (7) | **PASS** (optional wording, D6) |

The mathematics is right. No violation of `e_a <= (h - 1)/2`, of its dual, of the sole-exit law or of the floor
`(3h + 1)/2` was found on any tournament class of orders 3-9 (bound) and 3-8 (law), nor on any labelled tournament of
orders 3-6, nor by an independent permutation brute force at orders 3-7.

## Proof 7, step by step (all confirmed)

1. `O = N+(a)` nonempty; for `P = (x_1, ..., x_(n-1), a)` the least `i` with `x_i in O` exists and `i <= n - 2`
   automatically (`x_(n-1) -> a` since `P` ends `x_(n-1), a`, so `x_(n-1) notin O`).
2. `phi(P) = (x_1, ..., x_(i-1), a, x_i, ..., x_(n-1))` is a Hamiltonian path: `x_(i-1) -> a` because `x_(i-1) notin O`
   (`i >= 2`), `a -> x_i` because `x_i in O`; for `i = 1` it is `(a, x_1, ...)`, covered by the "if `i >= 2`" clause.
3. `phi(P)` does not end at `a` (`a` sits at position `i <= n - 2`) and its last vertex `x_(n-1)` beats `a`. Both are
   needed: they put `phi(P)` in `F`, and `F` is disjoint from the `e_a` paths ending at `a` and from the `e_O` paths
   ending in `O`, which is what makes `|F| = h - e_a - e_O` (equality, not only `<=`).
4. `phi` is injective (delete `a`, append it).
5. `e_O >= 1` when `e_a >= 1`: `a` is in the last strong component `X_s` (all Hamiltonian paths end there), `|X_s| >= 3`
   since `a` is not a sink, `a` has an out-neighbour `z in X_s`, Camion gives a Hamiltonian cycle of `X_s`, cut it
   before `z`, prefix Hamiltonian paths of the earlier components. Correct. It is also redundant: `2 e_a <= h - e_O <= h`
   and Rédei's parity (`h` odd, `2 e_a` even) already give `2 e_a <= h - 1`. What the paragraph buys is the sharper
   `2 e_a + e_O <= h`, which the statement does not record (D7); it holds with no exception through order 9.
6. Dual by the converse tournament: paths reverse, ends become starts, sinks become sources, `h` is unchanged. Correct.
7. Equality family: the only backward arc of `TT_(k+2)` with `(1, k + 2)` reversed is `k + 2 -> 1`, every
   non-transitive Hamiltonian path reads `(S increasing), k + 2, 1, (complement increasing)` with `S <= {2, ..., k + 1}`,
   every `S` works (`s -> k + 2` for `s >= 2`, `1 -> everything but k + 2`), `h = 2^k + 1`; `k + 1` is last iff
   `k + 1 notin S`: `2^(k-1)` paths. Correct for `k >= 1` (D3). In general `e_j = 2^(j-2)` for `2 <= j <= k + 1`,
   `e_1 = e_(k+2) = 1`, so `k + 1` is the unique attaining vertex for `k >= 2`, all three vertices for `k = 1` (`C3`).
8. Consequences. THM-4862 (a), quoted: "For `T + x` with `x -> a` and `v -> x` (`v != a`): `h(T + x) = h(T) - e_a(T)
   + [P A_a](T)`", `[P A_a](T) = sum_{S <= V - a} h(T[S]) s_a(T[V - S])`. The canon file of THM-4862 does not state the
   empty-prefix convention; its source packet does (`two_branch_rearming_sigma_sweep_20261104/PROOFS.md`, line 3:
   "`h(∅) = 1`"), and the `S = ∅` term `s_a(T)` is exactly what the paths starting at `a` need in the injection
   "split a Hamiltonian path of `X` at `a`" (D2). With it `[P A_a](X) >= h`, `sigma >= 0`, and `sigma =
   sum_S e_O(S) s_a(V - S)` exactly (paths of `X[S]` ending in `N+(a)`), i.e. the path-pair formula of statement 1
   with the cut `({a}, V - a)`. Arithmetic: `2h - (h - 1)/2 = (3h + 1)/2`, rise `>= (h + 1)/2`, and with statement 5's
   `k_min(h) = min {(H(X + v) - h)/2 > 0}` this is `k >= (h + 1)/4`. Correct. For a sink `a` the extension is
   `H`-neutral (`[P A_a] = h(X - a) = h`), so "non-sink" is the right hypothesis.
9. Relation to THM-4873 (i): sinks end `h`, non-sinks `<= (h - 1)/2 <= h - 3` for `h >= 5`, so no vertex ends `h - 1`
   or `h - 2`; at `h = 3`, `e_a = 1 = h - 2` is allowed and (census, every order `<= 9`) occurs only for
   `TT_j => C3` with `a` in the `C3`, as THM-4873 (i) states. Statement 7 also gives THM-4873 (i)'s first half
   (`e_a != h - 1`) for every `h >= 3`, since `(h - 1)/2 < h - 1` there.

## Corrections (quote, then repair)

- **D1 (statement 7, wording).** "the cases `h <= 3` of THM-4873 (i) are its `TT_j => C3` configurations": THM-4873
  (i) is stated for `H >= 3`, and at `h = 1` every non-sink ends `h - 1 = 0` paths, so there is no `h = 1` case.
  Repair: "the case `h = 3` of THM-4873 (i) (`e_a = h - 2 = 1`, only for `TT_j => C3` with `a` in the `C3`)".
- **D2 (statement 7, convention).** "`sigma >= 0` because splitting a Hamiltonian path of `X` at `a` is an injection
  into the pairs counted by `[P A_a]`": state `h(X[∅]) = 1` (THM-4862's canon text omits it; its packet has it), since
  the paths starting at `a` land in the `S = ∅` term.
- **D3 (statement and title).** "`TT_(k+2)` with the arc `(1, k + 2)` reversed (`h = 2^k + 1`...)": add `k >= 1`
  (`k = 0` is `TT_2`, `h = 1 != 2`). For `k = 1` the family member is `C3` and all three vertices attain the bound,
  which the proof's parenthetical ("for `m >= 4`") already says; the statement's "at the penultimate vertex" is fine.
- **D4 (statement 7, receipt).** "Checked exhaustively on all tournaments of orders `<= 9`: no violation; the maximum
  ... at order `9` is `64/129`": the repository holds no session receipt for order 9 (`endcount_ratio_n8.out` stops at
  order 8; `increment2_order9.py` computes no end counts and its `.out` has none; the commit message of b81b5ff2a7
  says the end-count output was carried there, but the file contains only the increment ladder). The claim is TRUE:
  audit C's `audit_C/del9.out` line 73 ("max e_a/H over non-sinks (H>=3): 64/129", committed in ed565d1871) and
  this audit's `audit_D/classes_n9.out` (all 191,536 classes, 1,716,944 non-sink pairs, 0 violations, max `64/129`
  at `R_9`, 24 s) both confirm it. Repair: cite one of these receipts (and audit C's `128/257` at order 10).
- **D5 (YAML `scripts:`).** The statement-7 scripts `endcount_ratio.py` (+ `_n8.out`) and
  `verify_stmt7_consequences.py` (+ `.out`) are not listed; add them. Note the session's own check of the sole-exit
  law covers labelled orders `<= 6` only; audit D extends it to every class of orders `<= 8` (and `9`, below).
- **D6 (title, optional).** "(dually for sources)" reads as if sources were bounded; say "(dually, a non-source
  starts at most `(h - 1)/2`)".
- **D7 (proof, optional).** The `e_O >= 1` paragraph is not needed for the stated bound (Rédei parity suffices,
  step 5); keep it and record what it proves: `2 e_a + e_O <= h`, i.e. `e_a <= (h - e_O)/2` with `e_O >= 1`, tight on
  the family (`e_O = e_(k+2) = 1`).

## Recomputed (audit D, own code: `endcount_audit_D.c`, subset DP; `crosscheck_perm_D.py`, permutations)

Classes from `gentourng` (whichever arc convention the ASCII format uses, the class list is closed under converse up
to isomorphism, so every class is processed once; both the bound and its dual were checked on every vertex).

| order | classes | non-sink `(X, a)` | max `e_a/h` (at) | equality `2e_a = h - 1`, `h >= 3` | sole-exit cases | `sigma = 0` | max `sigma` | `2H(X+v) = 3h + 1` | min `(H(X+v) - h)/h` |
|---|---|---|---|---|---|---|---|---|---|
| 3 | 2 | 5 | 1/3 (`C3`) | 3 | 5 | 3 | 3 | 3 | 2/3 |
| 4 | 4 | 14 | 2/5 (`R_4`) | 4 | 14 | 4 | 9 | 3 | 2/3 |
| 5 | 12 | 56 | 4/9 (`R_5`) | 5 | 56 | 7 | 27 | 3 | 2/3 |
| 6 | 56 | 324 | 8/17 (`R_6`) | 6 | 324 | 19 | 105 | 3 | 10/17 |
| 7 | 456 | 3,136 | 16/33 (`R_7`) | 7 | 3,136 | 68 | 387 | 3 | 18/33 |
| 8 | 6,880 | 54,584 | 32/65 (`R_8`) | 8 | 54,584 | 383 | 1,827 | 3 | 34/65 |
| 9 | 191,536 | 1,716,944 | 64/129 (`R_9`) | 9 | 1,716,944 | 3,659 | 8,187 | 3 | 66/129 |

`R_m` = `TT_m` with `(1, m)` reversed. Violations of the bound, of the dual, of `2e_a + e_O <= h`, of the law
`H(X + v) = 2h - e_a + sigma` (direct count of `H(X + v)` against `[P A_a]` from its definition), of `sigma >= 0`
and of `2H(X + v) >= 3h + 1`: none, at every order listed (order 9: all 191,536 classes, 1,716,944 non-sink
extensions, `classes_n9_soleexit.out`, about 8 minutes). `h >= 5` with `e_a in {h - 1, h - 2}`: none; `h = 3`, `e_a = 1`: always `TT_j => C3`, `a` in
the `C3`. Labelled orders 3-6: 18 + 224 + 4,800 + 190,464 = 195,506 non-sink cases, `sigma = 0` in 6 + 48 + 520 +
8,880 = 9,454, max `sigma` 105 (the session's three numbers reproduced); HYP-9469 (1)'s "53 of 420" is 33 non-sink
`sigma = 0` class pairs plus the 20 sinks (`sigma = 0` for a sink). The permutation brute force agrees with the DP on
every aggregate at orders 3-7 (classes) and 3-5 (labelled). Family: `m = 3..10`, `h = 2^(m-2) + 1`,
`e_(m-1) = 2^(m-3)`, unique attaining vertex for `m >= 4`, `sigma = 1` at `m - 1`, `H(X + v) = 3·2^(m-3) + 3`.

## Observations (not corrections; candidate sharpenings)

- Through order 9 the equality cases with `h >= 3` are exactly `TT_j => R_(k+2)` at vertex `k + 1` (`k >= 1`; the
  three `C3` vertices for `k = 1`): one class per `(j, k)`, `n - 2` per order, `h = 2^k + 1` only. The statement does
  not claim uniqueness; a proof would start from step 5 (`e_a = (h - 1)/2` forces `e_O = 1`, hence `M = 1` and
  `|O| = 1`).
- The floor `H(X + v) = (3h + 1)/2` is attained only at `h = 3` (`TT_j => C3`, `a` in the `C3`): the 3 "tight"
  pairs per order. On `R_(k+2)` `sigma = 1`, so the rise is `(h + 3)/2` and `k = (h + 3)/4` (`33/129` at order 9, audit
  C's least ratio). For `h >= 5` the census floor is `(3h + 3)/2`, `k >= (h + 3)/4`.
- `sigma = 0` is rare and shrinking: 60%, 29%, 12.5%, 5.9%, 2.2%, 0.70%, 0.21% of non-sink class pairs at orders 3-9.

## Priority

`git grep` over `01-canon` and `05-knowledge/hypotheses` (and `05-knowledge/results`) for `(h-1)/2`, `(H-1)/2`,
"ends at most", "end at most", "end count", "half of the paths": no statement of `e_a <= (h - 1)/2` or an equivalent
before 2026-11-04. The only neighbours are THM-4873 (i) (opus-S77, 2026-11-04: `e_a != H - 1`, `e_a = H - 2` iff
`TT_j => C3`; its lemma `e_v = M e_v(C)` combined with statement 7 on `C` gives `e_a <= (h - M)/2`), THM-4869's
"uniquely Hamiltonian end-count lemma" (a different object: end counts in uniquely Hamiltonian tournaments), and the
sole-exit law THM-4862 (a) (klein, 2026-11-04). Audit C's report (line 134) reached the same conclusion. The
`(H-1)/2` hits elsewhere (THM-1980, HYP-884/993/994/1041) are the quantity `(H - 1)/2 = ` number of odd-cycle
families, unrelated. Statement 7 is new within the repository.
