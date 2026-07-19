---
id: THM-1270
title: The D=9 rung CONFIRMED at N=30031 — M(F_9(30031)) = 9/270287 exactly (fourth consecutive out-of-sample tower confirmation, primorial 30030, window width 5.5·10⁻¹⁰, computed in ONE second) — with the COMPOSITE-SKIP LEMMA (two lines: at tower-N every composite binder's rung is dead because a prime factor of 2D−1 lies in the primorial and divides Q_D) and the TOWER-CLOSURE LEMMA (all lower prime rungs dead: N ≡ 1 mod p ⟹ p | Q_D), plus the GHOST-ENUMERATION evaluator (O(De) clearance per candidate; 142/142 validation rows exact).
status: >
  PROVED — the composite-skip and tower-closure lemmas (both from one congruence:
  Q_D ≡ 2D−1 mod ℓ when N ≡ 1 mod ℓ); the L1 floor M ≥ 9/270287 (THM-1271 Lemma 1,
  witness a = 31799, direct residue check over all 30031 speeds). PROOF-BACKED EXACT —
  M(F_9(30031)) = 9/270287 via THM-1271's e-channel reduction with the ghost-enumeration
  clearance (gated: 142/142 odd-N table rows vs the full evaluator, 0 mismatches;
  reproduces 6/1271 and 7/16183).
source: death-star-2026-07-19-S59e (HYP-7910; owner: run the D=9 rung at N=30031)
depends_on:
  - THM-1271  # Lemmas 1-3 + the e-channel reduction (the machinery)
  - THM-1257  # the tower being extended
related: [THM-1256, HYP-7905, kind-pasteur S128c86 (the tie side: F_2 Goddyn-Wong periodicity)]
scripts:
  - 04-computation/lrc_D9_rung_N30031_deathstar_S59e.py -> 05-knowledge/results/lrc_D9_rung_N30031_deathstar_S59e.out
---

# THM-1270 — the tower reaches primorial 30030

## 1. The result

```text
F_9(30031) = {1, …, 30029, 30031} ∪ {270270}        (30031 speeds; x = 9·30030)
M(F_9(30031)) = 9/270287  EXACTLY                    (Q = 30032·9 − 1 = x + 17)
rung ATTAINED, strictly inside W_30031 = (1/30032, 2/60063), width 5.544·10⁻¹⁰.
Witness: a = 31799 = 9·17⁻¹ mod Q (THM-1271 L1's closed form; min distance 9
over all 30031 speeds, checked directly).  Exact upper bound: the ghost
evaluator, 1 second.
```

Gate check: `30031 = 30030 + 1 ≡ 1 (mod L_9 = 2·3·5·7·11·13)` and
`30031 ≡ 9 ≢ 1 (mod 17)` — the prediction OPEN is confirmed. The tower's
confirmed rungs now ride the primorials **6, 30, 210, 2310, 30030**
(D = 3, 4, 6, 7, 9), with every confirmation after D=3 made out-of-sample.

## 2. The composite-skip and tower-closure lemmas (PROVED)

> **Lemma (one congruence, two consequences).** Let N ≡ 1 (mod ℓ) for an odd
> prime ℓ. Then `Q_D = (N+1)D − 1 ≡ 2D − 1 (mod ℓ)`, so `ℓ | Q_D ⟺ ℓ | 2D−1`.
> Consequently, at any tower-N (N ≡ 1 mod every odd prime below the level's
> binder):
> (a) **every lower PRIME binder p′ = 2D′−1 < p divides its own Q_{D′}** — and
> since `p′·a ≡ ±D′ (mod Q_{D′})` would force `p′ | D′` (impossible:
> `p′ = 2D′−1 > D′` and `gcd(p′, D′) | 1`), the D′-rung pair is unsolvable:
> **all lower prime rungs are DEAD**;
> (b) **every COMPOSITE binder 2D′−1 is dead the same way**: any prime factor
> `ℓ | 2D′−1` satisfies `ℓ ≤ (2D′−1)/3 <` the level binder, hence lies in the
> primorial, hence `ℓ | Q_{D′}` while `ℓ ∤ D′` — unsolvable again. ∎

Part (b) upgrades THM-1257's empirical "D=5 and D=8 never open" from
observation to theorem **at tower-N** (at generic N the composite skip remains
empirical). Verified at N=30031: gcd(2D−1, Q_D) = 5, 7, 3, 11, 13, 15 for
D = 3..8 — all > 1, all rungs dead; D=9 is the least live rung, and the
computation shows it is not merely live but attained.

## 3. The ghost-enumeration evaluator

THM-1271 reduced exact `M(F_D(N))` to the e-channel; the remaining cost was the
O(N) clearance sweep per candidate. The ghost enumeration replaces it: in the
channel at modulus S with `gcd(a, S) = g′`, the elements at distance r are the
representatives ≤ N of `u ≡ ±(r/g′)·(a/g′)⁻¹ (mod S/g′)` (empty unless
`g′ | r`; if `S/g′ ≤ N` the zero class kills the candidate outright). Since
`S/g′ > N`, each class has AT MOST ONE representative ≤ N, so

```text
c_eff = min( least r ≤ De with a base representative , De )
```

costs O(De) modular operations. Validation: 142/142 odd-N table rows
(D = 3..6, members AND non-members) match the full evaluator; F_6(211) and
F_7(2311) reproduced (0.1 s). N=30031 runs in ONE second — five orders of
magnitude below the naive cost. The next rung, D=10 (binder 19) at
`N = 510511` (≡ 1 mod 510510; 510511 mod 19 = 0 ≢ 1, gate predicted OPEN
→ M = 10/5105119), is within the same machinery's reach in minutes — filed as
a lead rather than run here.

## 4. Scope, typed

- PROVED: the two lemmas of §2; the floor M ≥ 9/270287 (L1 + direct check).
- PROOF-BACKED EXACT (lemmas L1–L3 + THM-1002 + gated computation): the upper
  bound, hence M = 9/270287 and window membership.
- The general-N gate law for D=9 (which other N open it) is untested beyond
  this point; the cascade predicts `{N ≡ 1 mod 30030, N ≢ 1 mod 17}` with
  possible extra branch conditions — same status as D=4's gate before the
  [8,100] sweep.
- Cross-link: kind-pasteur S128c86's Goddyn–Wong periodicity (F_2 ties the
  floor at every N ≡ 1 mod 6) is the m=2/tie face of this same
  binder-competition family; the tie side, the gap side (this tower), and the
  degrade side (F_3(31), F_5(211)) are one object with three regimes.
