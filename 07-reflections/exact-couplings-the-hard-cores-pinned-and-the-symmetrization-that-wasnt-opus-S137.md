---
source: opus-2026-07-07-S137
status: exact computations (engine coupling, hard cores, U-moments) + one theorem-let + one
  parity clarification; heavy same-hour fleet convergence documented honestly
tags:
  - lonely-runner
  - LRC14
  - hard-cores
  - intersected-ledger
  - symmetrization
  - palindrome-conjecture
  - exact-engine
---

# Exact couplings: the hard cores pinned, and the symmetrization that wasn't

**opus-2026-07-07-S137.** Owner worklist session (engine↔ledger coupling, the two hard cores,
palindromic-extremizer testing, symmetrized moments, μ tables, …). The fleet is now so dense
that four agents completed overlapping items within the same hours — monad-S4 (exhaustive
hard-core ledgers), mac-mini-S44 (cores + the same vacuity theorem-let, pushed minutes ahead),
kps-S63 (Schur refutation + palindrome refinement), boxeph (the 2-anchor discharge). Everything
below states what my runs add: **independent exact cross-validation, several exact firsts, and
two clarifications.**

## 1. The engine now computes the intersected ledger exactly

`ρ*_θ(P,E) = meas(G_P ∩ {maxgap_E > θ})` — G_P is an explicit rational interval list, so the
order-cell engine's per-subcell superlevel intervals intersect it exactly.
Validation: every THM-530 size-minimum reproduced (`6/7, 66/91, 55/91, 1979/4004, 2243/5880`,
`m_P` at the |P|=10 minimizer); **mac-mini-S42's probe minimum is now exact:
ρ*({1,5,9,11,13}, {0..7}) = 35771/105105 ≈ 0.340336** (their 80k grid: 0.3390).

## 2. The two hard cores are pinned (same-hour convergence with monad-S4)

Exact: `meas(G_{9..13}) = 40247/90090`, `meas(G_{10..13}) = 1577/3003`. My structured sweep +
exact descent found the minima **ρ*(CORE-8) = 246287/630630 ≈ 0.3905 (6.9× m_P)** and
**ρ*(CORE-9) = 1619/4290 ≈ 0.3774 (6.7× m_P)**, both at the **AP co-offset** — monad-S4's
independent exhaustive sweeps (3432/1287 shapes) give the *identical exact rationals*, so the
values are now double-engine-verified. Structural readout: the core minima sit at the
**compressed end** (diameter 7–8 — inside kps-S60's already-proved intersected-ledger zone),
while the spread/CRT regime that motivated the cores has `R ≈ 1.00` (even `R = 1.013 > 1` at
14-aligned spread — *positive* correlation) and 8–9× headroom; my exact descent could not push
below 0.446/0.524 there. The load-bearing residual of the window factoring is, in exact terms,
**compressed shapes already under the ledger + a spread tail where avoidance is trivially
strong.**

## 3. The symmetrization that wasn't (mac-mini-S44 first; independent proof + exact check)

The owner's "reversal-symmetrized moments" item dissolves: `config(E, 1−x)` is the reflection
of `config(E, x)` through 0 (`frac(−ex) = 1 − frac(ex)`), and reflection preserves the gap
multiset — so **U(E, 1−x) = U(E, x) pointwise**, and likewise every symmetric gap functional.
Reversal-symmetrization of U (or μ, or any gap moment) is *identically vacuous*; mirror
families have identical U functions, which is why mac-mini-S43's E[U]-min "mirror pair" is an
exact degeneracy. My exact cross-moments confirm: `E[U(x)U(1−x)]/E[U²] = 1.000` on all ten
families tested. (mac-mini-S44 pushed the same one-liner minutes earlier via THM-639-A — cede
priority, this is the independent confirmation.) The genuine symmetrization freedom left is in
**E-space arrangement** (which step word, given the multiset) and in **event families** (CE
cells under the dihedral action), not in x-space.

## 4. Exact firsts that stand

- **E[U] exact** (first exact values of opus-S131's uncovered length): `E[U](AP₁₃) =
  8599/67914`, monad record `14585/135828`, 2AP∪{13} `1847/16170`; the S41 mirror pair —
  identical exact value `170908753/1425854430`, as §3 forces.
- **G_P-restricted PZ certificates**: by Cauchy–Schwarz on the restricted measure,
  `ρ*(P,E) ≥ (∫_{G_P}U)²/∫_{G_P}U²` — exact rationals per family: 0.261–0.466 on the hard-core
  shapes = **4.6–8.2× m_P, rigorous per-family floors** (the moment computation is exact; the
  inequality is elementary). A uniform-over-E version of these two moments is now a clean,
  bounded target for mac-mini's PZ-on-U lane — on the cores, with G_P fixed, both moments are
  explicit finite objects.

## 5. The palindrome conjecture: parity clarification, then an exact refutation at even D

A palindromic 12-step word has even sum, so **at odd diameter no palindrome exists at all**:
odd-D minimizers are mirror pairs *by parity*, not variational preference — kps-S63's "D=15
min = mirror pair" is the parity triviality. The conjecture's content is even D only. And
there the exact per-D exhaustive sweep (all compositions up to reversal, gcd 1,
`lrc_palindromy_perD_exact_opus_S137.out`) settles it:

| D | classes | exact min μ | minimizer (steps) | palindromic? |
|---|---|---|---|---|
| 12 | 1 | 477/1078 ≈ .4425 | 1¹² (AP) | ✓ |
| 13 | 6 | .5589 | (1¹¹,2) | forced-no (parity) |
| 14 | 42 | 2539/5390 ≈ .4711 | (2,1¹⁰,2) | ✓ |
| 15 | 182 | .5604 | (1¹¹,4) | forced-no |
| 16 | 693 | 4166/8085 ≈ .5153 | (2,2,1⁸,2,2) | ✓ |
| 17 | 2184 | .5556 | (1¹¹,6) | forced-no |
| **18** | **6216** | **26821/48510 ≈ .5529** | **(1¹¹,7) = {1..12,19}** | **✗ REFUTES even-D palindromy** (best palindrome .5673) |
| 19 | 15912 | 272/495 ≈ .5495 | (1¹¹,8) | forced-no |

**The true per-D landscape:** two known structures trade the lead — the **palindromic
2-block interlacings** (D=14: (2,1¹⁰,2); D=16: (2,2,1⁸,2,2); D=20: (2⁴,1⁴,2⁴) = monad's
record; D=22: (2⁵,1²,2⁵) = the prim-sat family) and the **far-element ladder
{1,…,12,w}** (D=15..19), whose μ decreases in w toward the deep-well regime. Neither wins
uniformly; the crossovers are non-monotone. So the palindromic-extremizer conjecture is
**false per-D** (exact counterexample at D=18) and survives only as: *the global minimizer
(the AP) is palindromic*, plus the observation that all strong palindromic performers are
2-block interlacings — the structures the parity/bisection mechanism already explained.

## 6. Where the (A′) question now stands (fleet state, for the record)

boxeph's same-day fusion is the milestone: the skeleton needs only `μ ≥ T_k` (monad's weakened
bars), and the **2-anchor tail** `P(max(gap∋0, gap∋1/2) > 1/7)` clears every `T_k` with margins
0.15–0.31 — my S134 anchored-tail statistic, which *fails* at the exact AP bar, *suffices* at
the DAG bar. Full (A′) extremality (my S136 certified boxes, isolation gaps γ_k) remains true
and now decorative to the DAG: the load-bearing object is the 2-anchor joint law — observer +
antipode — where the 2-adic step structure (THM-580, anchor ½) is the natural proof tool.
