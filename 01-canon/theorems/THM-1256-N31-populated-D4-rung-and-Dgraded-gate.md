---
id: THM-1256
title: The N=31 first gap is POPULATED — {1..29, 31, 120} attains the D=4 slack-1 rung 4/127 (refuting THM-1255 §4's band-law prediction), the single-far strata are completely classified for N = 14..31 (members exactly 3/59, 3/77, 4/127 at N = 19, 25, 31), and the mod-30 mediant gate is the D=3 slice of a D-GRADED gate family (binder 2D−1 at Q = (N+1)D−1) with an arithmetic see-saw — N ≡ 1 (mod 5) closes the D=3 gate and opens the D=4 gate at once.
status: >
  PROVED — M({1..29,31,120}) = 4/127 exactly (three independent methods: pair-sum
  evaluator, witness residue check at t = 55/127, full scan over all q ≤ 260 ⊇ all
  candidate denominators); the single-far classification N = 14..31 (absorption lemma +
  finite check, PER-INSTANCE base-floor certificates M(B) > 2/(2N+1) computed exactly —
  no open-LRC citation used; zero base-floor alarms, so incidentally no LRC(N+1)
  near-violation among the 18·N bases). VERIFIED (evidence-grade) — no additional N=31
  member across bordered (452,905), two-defect (342,432), targeted-multiple (10,440)
  and needle-repair species; species D (band descent) yielded 0 within budget — typed
  honestly as NO EVIDENCE from that instrument, not a negative.
source: death-star-2026-07-19-S59b (HYP-7890; the owner-directed test of THM-1255 §4's N=31 prediction)
depends_on:
  - THM-1255  # the absorption lemma + classification machinery (its §4 conjecture is corrected by this file)
  - THM-1002  # pair-sum denominator lemma (exactness of the evaluator)
related:
  - HYP-4516  # the mod-30 binder gate = the D=3 slice of the D-graded gate
  - HYP-7890  # this test
  - HYP-7840 / THM-1240  # opus bound-D at N=13: the D=4 rung 4/55 is the same object one gate down
  - mac-mini-S29  # M(F(31)) = 1/32 degrade (replicated here, P1)
scripts:
  - 04-computation/lrc_N31_bandlaw_test_deathstar_S59b.py -> 05-knowledge/results/lrc_N31_bandlaw_test_deathstar_S59b.out
---

# THM-1256 — N=31 is populated: the D=4 rung opens where the mediant dies

## 1. The discovery

```text
V = {1,...,29, 31, 120}   (31 speeds, primitive; = {1..31}\{30} ∪ {4·30})
M(V) = 4/127  ∈  (1/32, 2/63) = W_31        [witness t = 55/127; active pair (7,120)]
```

Verified three ways: the THM-1002 pair-sum evaluator; the direct residue check at
`t = 55/127` (all `55v mod 127` at distance ≥ 4, with the pair `(7,120)` at `±4`); and
a brute maximum over ALL `q ≤ 260 ⊇ {2·max = 240}`, which by THM-1002 §1 is a superset
of every candidate denominator. `127 = 32·4 − 1`: the family sits exactly on the
**slack-1 rung** `D/((N+1)D − 1)` at `D = 4`, order `k = q − N·D = 3`.

This **refutes** THM-1255 §4's band-law prediction (N=31 empty), which I posed
yesterday-equivalent and tested today at the owner's direction. The instrument that
refuted it is THM-1255's own absorption classification — the conjecture and its
counterexample came out of the same machinery, in the right order.

## 2. The complete single-far classification, N = 14..31 (PROVED)

For every `N ∈ [14,31]` and every defect `i`, the base `B = {1..N}∖{i}` has more than
12 speeds, so the settled-LRC floor is unavailable; instead each base carries a
**per-instance certificate**: `M(B)` computed exactly and checked `> θ = 2/(2N+1)`
(zero failures across all 18·N bases — incidentally, no LRC(N+1) near-violation among
them). Then THM-1255's absorption lemma gives finite `X₀` and the finite check decides
the stratum. Complete member list over all `N ∈ [14,31]`, all `i`, ALL `x` unbounded:

```text
N = 19:  {1..17, 19, 54}    M = 3/59     (D=3 rung; i=18, x=3·18, binder 5)
N = 25:  {1..23, 25, 72}    M = 3/77     (D=3 rung; i=24, x=3·24, binder 5; 77=7·11 composite)
N = 31:  {1..29, 31, 120}   M = 4/127    (D=4 rung; i=30, x=4·30, binder 7)
ALL OTHER (N, i, x): NONE.
```

The N=19 and N=25 rows are simultaneously the validation gates (predicted in advance by
the D=3 gate; found organically by the instrument) — the N=31 result is produced by a
detector proven able to find exactly this kind of member. `M(F(31)) = 1/32` (the
canonical D=3 family's degrade to the floor, mac-mini-S29) is replicated exactly.

## 3. The D-graded gate (the corrected mechanism)

The canonical-D single-far family is `F_D(N) = {1..N}∖{N−1} ∪ {D(N−1)}`. Its far
element binds `b = 2D−1` at `Q = D(N−1) + (2D−1) = (N+1)D − 1` — the slack-1
denominator — so the attained value, when every looser competitor branch dies, is the
slack-1 rung `D/((N+1)D−1)`. HYP-4516's mod-30 gate is the `D = 3` slice
(`b = 5`, `Q = 3N+2`). At `D = 4` the competitor at `b = 5` sits at `Q' = 4N+1`, and

```text
5 | 4N+1  ⟺  N ≡ 1 (mod 5)  ⟺  5 | 3N+2,
```

so **the congruence that closes the D=3 gate is the same one that removes the D=4
gate's b=5 competitor**: an arithmetic see-saw. (Parity kills b=2 at `Q' = 4N−2` for
free; `3 | 4N−1 ⟺ N ≡ 1 (mod 3)` kills b=3; N = 31 ≡ 1 (mod 30) satisfies everything
at once.) The see-saw is necessary-not-sufficient: N = 16 ≡ 1 (mod 15) has NO member
(P2), so the full D=4 gate carries further conditions — extracting its exact mod-lcm
form from the P2 data, and testing its predicted next opening (N = 61?), are the named
follow-ups.

## 4. What this feeds

- **opus THM-1240/HYP-7840 (the N=13 wall):** the D=4 slack-1 rung IS `4/55` at N=13 —
  the mediant of `(1/14, 3/41)`. This theorem shows D=4 rungs are REALIZED in nature
  when the gate arithmetic conspires, so `(1/14, 3/41)` emptiness cannot lean on any
  "D ≥ 4 is never attained" heuristic — it must be N=13-specific. Concretely:
  `F_4(13) = {1..11, 13, 48}` fails because its b=5 competitor at `Q' = 53` is ALIVE
  (`53 ≢ 0 mod 5`, N = 13 ≡ 3 mod 5) — and my S59 atlas proves nothing single-far
  enters `(1/14, 3/41)` at all. The remaining danger for opus's interval is
  non-single-far D=4 realizers, exactly where boxeph-S130's mod-17/19/23 CRT stack
  should be pointed.
- **The rung ladder reading (boxeph-S130):** attained first-gap values across N are
  slack-1 rungs `D/((N+1)D−1)` with D graded by the gate — the certificate-rung
  ladder and the gate family are one object seen from the two sides.
- **Scope, typed honestly (the boxeph-S130 rule):** single-far strata N=14..31 =
  PROVED CLASSIFIED (bounded check with per-instance certificates; bounds in the
  out-file). Non-single-far species at N=31 = bounded evidence (counts above).
  Non-single-far species at N=14..30 = NOT TESTED here at all. Full first-gap
  emptiness at any N ≥ 14 outside {19,25,31} = OPEN (only its single-far stratum is
  closed).
