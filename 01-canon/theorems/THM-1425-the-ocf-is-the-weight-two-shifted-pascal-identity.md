---
id: THM-1425
title: "THE OCF IS THE SHIFTED-PASCAL FIBONACCI IDENTITY WITH WEIGHT 2 INSTEAD OF 1 — Σ_k C(m−k,k)·1^k = F(m+1) counts independent sets in the path; Σ_k C(m−k,k)·2^k = J(m+1) = (2^{m+1}−(−1)^{m+1})/3 is the SAME sum with the OCF's weight, so whenever the odd-cycle intersection graph is a path the Hamiltonian-path count is a JACOBSTHAL number. The OCF itself is confirmed exhaustively: H(T) = Σ over sets of vertex-disjoint directed odd cycles of 2^{#cycles}, zero mismatches over all 8/64/1024/32768 labelled tournaments at n = 3..6. RÉDEI IS ITS MOD-2 SHADOW — every nonempty collection contributes an even 2^{|S|}, so only the empty collection survives mod 2 and H ≡ 1, a one-line derivation of oddness. AND THE 60 IS REAL BUT IT MOVES: Fibonacci's last digit has period 60 mod 10 and 560 mod 1001, while the weight-2 sequence has period 4 mod 10 and EXACTLY 60 mod 1001 = 7·11·13 — because the weight-2 period mod q is governed by ord_q(2) and ord_1001(2) = 60. Same constant, opposite modulus, swapped by the very weight-1 → weight-2 substitution that turns Fibonacci into the OCF."
status: >
  OCF form VERIFIED-EXACT by exhaustive enumeration of every labelled tournament at
  n = 3,4,5,6 (8/64/1024/32768; zero mismatches against the subset-DP Hamiltonian-path
  count).  The shifted-Pascal identities are PROVED (standard transfer matrix) and
  verified m = 0..11.  Rédei-as-mod-2-shadow is PROVED given the OCF form.  All period
  claims are COMPUTED (eventual period with explicit pre-period), not asserted.
  The "odd-cycle intersection graph is a path" hypothesis is a CONDITIONAL — no claim is
  made here about how often it holds; see Honest scope.
source: mac-mini-2026-07-20-S129 (owner: "think how fibonacci arises from summing pascals
  triangle shifted"; and "1001 = three sixties relates to the fact the final digit of the
  Fibonacci sequence in base 10 repeats every 60 terms")
related:
  - THM-070   # the OCF -> Claim A chain
  - THM-1415  # kind-pasteur §III/§IV: the "three sixties" deflation and the Faulhaber
              # diagonal law -- a DIFFERENT and complementary attack on the same prompt
  - THM-1420  # this session's companion: no linear invariants
script: 04-computation/ocf_is_weighted_fibonacci_macmini_S129.py (+ .out)
---

# THM-1425 — the OCF is Fibonacci with weight 2

**One line.** Fibonacci counts independent sets in a path. The OCF counts the same thing with
each chosen element worth `2` instead of `1`. That single substitution turns `F` into the
Jacobsthal numbers, turns the shifted-Pascal identity into the tournament formula, makes
Rédei's oddness a one-line remark — and moves the number 60 from modulus 10 to modulus 1001.

## (A) The OCF form, confirmed exhaustively

> `H(T) = Σ` over sets of **vertex-disjoint directed odd cycles** of `2^{#cycles}`.

| `n` | labelled tournaments | mismatches vs subset-DP `H` |
|---|---|---|
| 3 | 8 | **0** |
| 4 | 64 | **0** |
| 5 | 1 024 | **0** |
| 6 | 32 768 | **0** |

## (B) The shifted-Pascal bridge

The shallow diagonals of Pascal's triangle sum to Fibonacci because `C(m−k,k)` counts the
independent sets of size `k` in the path `P_m`. Give each chosen vertex a weight `w`:

| weight | `Σ_k C(m−k,k)·w^k` | sequence |
|---|---|---|
| `w = 1` | `F(m+1)` | 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144 |
| `w = 2` | `J(m+1) = (2^{m+1}−(−1)^{m+1})/3` | 1, 1, 3, 5, 11, 21, 43, 85, 171, 341, 683, 1365 |

Both verified `m = 0..11`. **The OCF is the `w = 2` sum on the odd-cycle intersection
graph** (vertices = directed odd cycles, edges = share a vertex; "independent" = vertex-
disjoint). Hence:

> **Corollary.** If a tournament's odd-cycle intersection graph is a path `P_m`, then
> `H(T) = J(m+1) = (2^{m+1} − (−1)^{m+1})/3` — a Jacobsthal number.

## (C) Rédei as the mod-2 shadow

Every **nonempty** collection contributes `2^{|S|} ≡ 0 (mod 2)`; the empty collection
contributes `1`. Therefore `H(T) ≡ 1 (mod 2)`:

> **Rédei's oddness is one line from the OCF.** (Verified: every `H` odd at `n = 3..6`.)

This also explains the divisibility `3 |` nothing and the oddness of every Jacobsthal number
`J(k)`, `k ≥ 1` — the same `2^{|S|}` mechanism. It further explains THM-1420(A): the tiling
fibre `t = H/|Aut|` is odd because `H` is, which forces `A000568(n)` even.

## (D) The 60 is real, and it MOVES

The prompt's observation — Fibonacci's last digit repeats every 60 — is correct, and the
weight-2 substitution does something better than destroy it. Computed eventual periods
(pre-period shown where nonzero):

| modulus | Fibonacci (`w=1`) | weight-2 (`w=2`) |
|---|---|---|
| 2 | 3 | 1 *(pre 1)* |
| 5 | 20 | 4 |
| **10** | **60** | 4 *(pre 1)* |
| **1001 = 7·11·13** | 560 | **60** |

> **The 60 does not vanish under the tournament weighting — it migrates from modulus 10 to
> modulus 1001.**

The reason is structural, not numerological. Fibonacci's period mod 10 is
`lcm(π(2), π(5)) = lcm(3, 20) = 60`. The weight-2 sequence has the closed form
`J(k) = (2^k − (−1)^k)/3`, so **its period modulo `q` is governed by `ord_q(2)`** — and
`ord_1001(2) = 60`, which is also exactly why base-1000 digit grouping detects divisibility
by 7, 11 and 13. Two genuinely different 60s, exchanged by the same substitution that turns
Fibonacci into the OCF.

This is a **complement to, not a contradiction of, THM-1415 §III**, which deflates the
"three sixties" coincidence quantitatively (60 is the second most frequent lcm of subsets of
`{1..12}`, at 4.90%, so bare agreement at 60 is cheap). That deflation stands. What is added
here is that *one* of the sixties — `ord_1001(2)` — is not a bare coincidence but the actual
period of the repo's own weighted recursion, reached from the Fibonacci side by the exact
substitution the OCF performs.

## Honest scope

- The OCF form is **verified**, not proved here; it is the repo's existing formula (THM-070
  and the GS-OCF bridge) and this is an independent exhaustive confirmation at `n ≤ 6`.
- **The path hypothesis in (B) is conditional.** Nothing here says how often a tournament's
  odd-cycle intersection graph is a path — for most tournaments it is dense, and at large `n`
  the corollary applies to few or no tournaments. The identity is exact where it applies;
  its *scope* is unmeasured. Measuring it is the obvious next step.
- The Jacobsthal numbers here are the **sequence** `(2^k−(−1)^k)/3`, unrelated to the
  **Jacobsthal function** `g(n)` that appears in the LRC threads. Same name, different object
  — do not cross-link them.
- Period claims are eventual periods with explicit pre-periods (the `w=2` sequence is not
  purely periodic mod 2 or mod 10; an earlier run of this script reported "None" for those
  because it tested pure periodicity only).
- (D) is an identity plus an observation about where a constant lands. It is not evidence
  that 1001 or 60 plays a structural role in tournaments.

*Artifacts:* `04-computation/ocf_is_weighted_fibonacci_macmini_S129.py` (+out).
*Credits:* THM-070 / the GS-OCF bridge for the OCF; THM-1415 §III–IV (kind-pasteur) for the
concurrent and complementary attack on the same prompt.
