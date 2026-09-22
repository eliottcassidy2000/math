# Lean paste audit, wave 6: square-sum graph loops, the Q_14 merger, the Delta=4 ladder, and the vacuous sign field

**Status:** PROVED (core Lean 4.30.0, `decide`, no Mathlib, no `sorry`, no `native_decide`, axioms at most `propext, Quot.sound`; scratch file, not a build target): `Q_14` connected, `Q_13` and `Q_12` disconnected, component counts `3, 2, 1, 1, 1, 1` for `n = 12..17`, the pasted adjacency reflexive at value `2`, the ladder statement, `15` and `21` composite, the pasted certificate structure inhabited for every `n >= 1`. FINITE-EXACT: the session lead's `Q_n` probe (components `n <= 32`, degree-`<= 1` vertices `n = 14..32`, Hamiltonian paths `n <= 32` by exhaustive backtracking with witnesses). REFUTED (minimal witnesses): `square_sum_graph` as a `SimpleGraph` (value `2`, `2+2 = 4`), "merger of the 3 historical components at `N = 14`" (`2` components at `13`), the "`Delta = 4` prime ladder" (`15 = 3*5`, `21 = 3*7`, `17-11 = 6`), `is_positive_sheet : True` as a sign guard (a one-inhabitant field), "at `N = 15` the number `4` is the edge to avoid" (`4` is interior, `12,4,5`). CITED: the minimal descent-certificate type (compression_and_lean_audit S15); OEIS A090461. SCOPE: no map found from the summand-graph triple `{1,4,6}` (THM-2422) to the three components of `Q_12`. OPEN: unchanged (`GlobalTallyMargin`). Results note of session `collatz-mod6-20260917` (machine `mac-mini`), wave 2026-09-22, lane `lean_paste_audit_w6`; not a reserved canon ID.

Script [collatz_mod6_20260922_w6_lean_paste_audit_w6.py](../../04-computation/experiments/collatz_mod6_20260922_w6_lean_paste_audit_w6.py), output [collatz_mod6_20260922_w6_lean_paste_audit_w6.out](collatz_mod6_20260922_w6_lean_paste_audit_w6.out), Lean scratch file [collatz_mod6_20260922_w6_lean_paste_audit.lean](../../04-computation/lean/standalone/collatz_mod6_20260922_w6_lean_paste_audit.lean). Claims S1-S9 match the `.out` sections.

## Inheritance and concept board

The pasted Lean block is the third one this session has audited; the first two are typed in [compression_and_lean_audit](collatz_mod6_20260922_compression_and_lean_audit.md) (sections 6-7: `4` theorems, `4` sorries, nothing elaborating; the correct minimal `DescentCertificate`) and [paley_fano_octonion_design](collatz_mod6_20260922_paley_fano_octonion_design.md). The formal side inherits the opus [CollatzBlueprintAudit package](../../04-computation/lean/CollatzBlueprintAudit/README.md) (`descent_iff_affine_inequality`, `GlobalTallyMargin`) and the standalone files `collatz_mod6_20260922_descent_certificate.lean` and `collatz_mod6_20260922_certificate_iff_descent_audit.lean` under `04-computation/lean/standalone/`; the sheet-blind / sign-specific typing is the [counterexample portrait](collatz_mod6_20260922_counterexample_portrait.md). The summand graph and its `{1,4,6}` exclusion are canon THM-2422 (`01-canon/theorems/`), the reflection [summand-graph-fermat-zeckendorf](../../07-reflections/summand-graph-fermat-zeckendorf.md) and [arithmetic_braids_20260917_summand](arithmetic_braids_20260917_summand.md); nothing about them is re-derived here. **Closest proved mechanism:** `decide` on a Boolean reachability closure is exactly the CatalanEllipticAudit pattern (`12` axiom-free certificates, cited by path above); here it settles the square-sum connectivity questions in one `lake env lean` call (S8). **Canonical hostile:** the value `2` of `Fin 14` at index `1`: `2+2 = 4` makes the pasted `Adj` reflexive, so `loopless := sorry` can never be discharged and the `SimpleGraph` does not exist as pasted (S1). **Corrected near miss:** `square_sum_fourteen_connectivity` is a TRUE statement (S2, S8), but the comment attached to it is off by one merger: `Q_12` has `3` components `[1,3,6,8,10], [2,7,9], [4,5,11,12]`, `13` joins the first and third (`13+3 = 16`, `13+12 = 25`) leaving `2`, and `14` joins `[2,7,9]` (`14+2 = 16`, `14+11 = 25`) leaving `1`. **Least-used sidecar:** `trivial` DOES close the ladder statement in core Lean without `Finset`/`Set` (it delegates to `decide` on list membership, S6/S8), so the paste's one working tactic was working for a reason the paste did not state; the reason also shows why the statement is empty (it is a `4`-element table lookup).

Typed non-analogy (SCOPE, no map found). Source: the summand-graph exclusion set `{1,4,6}` of THM-2422 (three numbers that the distinct-summand closure of `{2,3}` never produces; the reflection's "three modes" unary / binary bridge / ternary). Target: the three components of `Q_12`. Map attempted: send each excluded number to the component containing it. Result: `1` and `6` lie in the SAME component `[1,3,6,8,10]` and `4` in `[4,5,11,12]`, while the third component `[2,7,9]` contains none of `{1,4,6}`; the three summand modes are not the three square-sum components. Preserved: nothing; lost: the triple itself; sidecar: `2` (the loop witness of S1) is the lowest vertex of the last component to merge. Test: S2 of the `.out`.

## 1. The pasted adjacency is not loopless (S1, REFUTED)

`Adj x y := exists k, (x.val+1)+(y.val+1) = k^2` on `Fin n` allows `x = y`. Self-loops with value `<= 32` sit at `(index, value) = (1,2), (7,8), (17,18), (31,32)` (values `2v` square: `4, 16, 36, 64`); among values `1..14` there are `2` of them, at `2` and `8`. Lean: `adjPaste_not_irreflexive : adjPaste 1 1 = true` and `pasted_loopless_fails : ¬ (∀ i, i < 14 → adjPaste i i = false)`, both by `decide`. The fix is `x != y`; with it `adj_irrefl_upto_32` and `adj_symm_upto_32` hold by `decide`. The session lead's remark that `18` has only the neighbour `7` until `31` appears is exactly the `(17,18)` loop being excluded (S3: `18`'s neighbours in `Q_30` are `[7]`, in `Q_31` `[7, 31]`).

## 2. Connectivity and the merger, formalized (S2, S8; PROVED)

Values `1..n`; `stepR` adds one Boolean matrix-vector round; `closure n (n-1) s` is reachability from `s`; `connected n := all v, reach n 1 v`; `numComponents n` counts the values not reached from any smaller value. All by `decide`:

| theorem | statement | tactic |
|---|---|---|
| `square_sum_fourteen_connectivity` | `connected 14 = true` | `decide` |
| `square_sum_thirteen_disconnected` | `connected 13 = false` | `decide` |
| `square_sum_twelve_disconnected` | `connected 12 = false` | `decide` |
| `components_twelve .. components_seventeen` | `numComponents 12..17 = 3, 2, 1, 1, 1, 1` | `decide` |

Independent Python census (S2): components `1,2,2` for `n = 1,2,3`, `3` for `4 <= n <= 12`, `2` at `13`, `1` for `14 <= n <= 32`; edge counts `0,0,1,1,2,...,46`. The session lead probe is confirmed exactly. The pasted comment "merger of the 3 historical components" is REFUTED: `Q_13` has `2` components, and the `14`-merger joins `[2,7,9]` to the rest.

## 3. Leaves and Hamiltonian paths (S3, S4; FINITE-EXACT)

Degree-`<= 1` vertices: `n = 14: [8,9,10]`, `15: [8,9]`, `16: [8,16]`, `17: [16,17]`, `18: [16,17,18]`, `19: [16,18]`, `20..30: [18]`, `31, 32: []`. Hamiltonian paths (exhaustive backtracking, witnesses printed in the `.out`) exist for `n = 1, 15, 16, 17, 23, 25, 26, 27, 28, 29, 30, 31, 32` and for no other `n <= 32`; this is the prefix of OEIS A090461 (CITED; the script's in-sandbox fetch failed with `URLError`, a manual `curl` with the research user agent returned the sequence `15, 16, 17, 23, 25, ...` and the name "Numbers k for which there exists a permutation of the numbers 1 to k such that the sum of adjacent numbers is a square"; its comment that every `k >= 25` qualifies is CITED from the session lead's report, not checked here). Both lead paths validate. The `n = 15` lead path uses only the squares `[9, 16, 25]` (the paste's "braiding `9,16,25`" is a correct description of that witness), the `n = 23` path uses `[4, 9, 16, 25, 36]`. In `Q_15`, `deg(4) = 2` with neighbours `[5, 12]` and `4` is interior to the lead path (`12, 4, 5`); the endpoints are the two leaves `8, 9` (degrees `1, 1`), so "the number `4` is the edge to avoid" is REFUTED and "the `N = 15` endpoints are `8` and `9`" is the forced statement. `n = 18` fails by the leaf count alone (`3` leaves); `n = 19..22, 24` fail only by the exhaustive search.

## 4. The Delta = 4 ladder (S6; TRUE statement, REFUTED reading)

`[3,7,11,17] + 4 = [7,11,15,21]`, so `delta_four_prime_step_invariance` is TRUE; in core Lean both `by decide` (list `all`/`contains` form) and `by trivial` (Prop form `∀ p, p ∈ [3,7,11,17] → p+4 ∈ [7,11,15,21]`) close it (S8, exit `0`). But `15 = 3*5` and `21 = 3*7` are composite (`fifteen_not_prime`, `twentyone_not_prime`), the differences are `[4, 4, 6]` (`ladder_not_ap : 17 - 11 ≠ 7 - 3`), and the only difference-`4` prime triple below `100` is `(3, 7, 11)`. The "prime progression regulating the carry" is a `4`-row lookup table with two composite outputs; it carries no Collatz content.

## 5. The vacuous sign field and the correct certificate (S7; REFUTED as a guard, CITED shape)

`is_positive_sheet : True` has exactly one inhabitant; `is_positive_sheet_unique` (two certificates at `n = 7` have `rfl`-equal sign fields). With `L, K_L, B_L` free the inequality `3^L n + B_L < 2^{K_L} n` is met for every `n >= 1` by `L = 0, K_L = 1, B_L = 0` (`n = 1..1000` checked; Lean `pastedCertificateAlwaysInhabited` for all `0 < n`), and no orbit is mentioned, so the structure is inhabited on both sheets for every start: it is neither a sign guard nor a descent certificate. The correct minimal type is the cited `DescentCertificate n := ∃ t K L B, 0 < t ∧ 2^K * iterate collatz t n = 3^L * n + B ∧ 3^L * n + B < 2^K * n` (compression_and_lean_audit S15, where it is also proved equivalent to `∃ t > 0, iterate collatz t n < n`). Witnesses restated here: plus sheet `n = 3, t = 6`, orbit `3,10,5,16,8,4,2`, `K = 4, L = 2, B = 5`, `32 = 32 < 48`; minus sheet (`3n-1`) `n = 3, t = 3`, orbit `3,8,4,2`, `K = 3, L = 1, B = 7`, `16 = 16 < 24` (`descent_three`, `descent_three_minus`, axiom-free). The certificate shape is SHEET-BLIND by construction in the portrait's sense: the sign enters only through which map `iterate` binds, never through the inequality. The content-bearing target remains the package's `GlobalTallyMargin` with the actual counters (OPEN).

## 6. Verdict table (S9)

| pasted item | verdict | witness |
|---|---|---|
| `square_sum_graph` as pasted | REFUTED as `SimpleGraph` | value `2`: `2+2 = 4`; fix `x != y` |
| `square_sum_fourteen_connectivity` | TRUE, PROVED by `decide` | `connected 14 = true` |
| "merger of the 3 historical components" | REFUTED (off by one) | components `3, 2, 1` at `12, 13, 14` |
| `Q_13`, `Q_12` disconnected | PROVED by `decide` | `connected 13 = false`, `connected 12 = false` |
| `delta_four_prime_step_invariance` | TRUE (`trivial` and `decide`) | `[7,11,15,21]` |
| "`Delta = 4` prime ladder / AP" | REFUTED | `15 = 3*5`, `21 = 3*7`, `17-11 = 6` |
| `is_positive_sheet : True` | REFUTED as sign guard (vacuous) | inhabited for all `n >= 1` with `L=0, K=1, B=0` |
| minimal correct certificate | CITED (compression_and_lean_audit S15) | sheet-blind; `n = 3` witnesses on both sheets |
| session lead `Q_n` probe | FINITE-EXACT confirmed | S2-S4 |
| "at `N = 15`, `4` is the edge to avoid" | REFUTED | `4` interior (`12,4,5`), `deg(4) = 2` |

## 7. Reproduction

```text
cd /tmp/math-wt-collatz-mod6-b/04-computation/experiments
python3 collatz_mod6_20260922_w6_lean_paste_audit_w6.py > ../../05-knowledge/results/collatz_mod6_20260922_w6_lean_paste_audit_w6.out
python3 -O collatz_mod6_20260922_w6_lean_paste_audit_w6.py   # identical except the two timing lines
# the Lean step the script runs (exit 0, 21 theorems, 0 sorry in code, 0 axiom declarations, 0 imports):
cd /tmp/math-wt-collatz-mod6-b/04-computation/lean/CollatzBlueprintAudit
lake env lean /tmp/math-wt-collatz-mod6-b/04-computation/lean/standalone/collatz_mod6_20260922_w6_lean_paste_audit.lean
```

Toolchain `leanprover/lean4:v4.30.0`; `#print axioms` reports `[propext, Quot.sound]` for the `decide` theorems and no axioms for the certificate witnesses. Runtime about `6` s, memory negligible.

## Stopping boundary / next question

Everything in the paste's Lean block is now either machine-checked (connectivity, disconnectedness, the ladder table, the loop) or refuted at its smallest instance (loops, the merger count, the composite ladder targets, the vacuous sign field); no residual formal object remains on the square-sum side, and the summand-graph triple `{1,4,6}` has no map to the square-sum components. The one open formal target is unchanged and already stated in the package (`GlobalTallyMargin`). The next question with content is the compression lane's: a residue class `n = r mod 2^k` (the guards note's `23 mod 256` first-reset cylinder is the candidate) for which the certificate time and counters are given by a closed formula, so that an infinite family `∀ j, DescentCertificate (r + 2^k j)` can be proved by arithmetic rather than `decide` on single starts. Not attempted here.
