---
id: THM-1065
title: The doubling family is tight EXACTLY when n ≡ 1 (mod 6) — for G(n) = {1,…,n−2, n, 2(n−1)} (remove n−1, add its double), M(G(n)) = 1/n if n is even, 2/(2n+1) if n ≡ 3,5 (mod 6), and 1/(n+1) — TIGHT — iff n ≡ 1 (mod 6). The mechanism is exact: G(n) misses precisely the residue class ±(n−1) modulo BOTH 2n and 2n+1, so a clearance-2 witness exists at q iff (n−1) is invertible mod q; and gcd(n−1,2n) = gcd(n−1,2) while gcd(n−1,2n+1) = gcd(n−1,3). Both witnesses die exactly when 2 | (n−1) AND 3 | (n−1), i.e. 6 | (n−1). THE "6" IS 2·3 FROM THE TWO NEIGHBOURING DENOMINATORS. This yields an infinite family of sporadic tight instances (n = 7,13,19,25,31,…), contains Goddyn–Wong (n=13) as its second member, and PROVABLY MISSES n=12.
status: The forward direction (6 ∤ (n−1) ⟹ NOT tight) is PROVED — explicit witnesses, structurally verified n=5..39. The converse (6 | (n−1) ⟹ tight) is VERIFIED-EXACT for n ≤ 60 (the full M formula matched with ZERO mismatches) but NOT proved in general: it needs "no witness other than q ∈ {2n, 2n+1, n+1} beats 1/(n+1)".
source: mac-mini-2026-07-18-S118 (owner: characterize n where remove n−1, add 2(n−1) is tight)
depends_on:
  - THM-1031  # the covering lemma and the Farey sup-companion
  - THM-1001  # the safe-interval element bound
external: LRC(n) SETTLED (for M ≥ 1/(n+1)).
related:
  - HYP-7460  # the S117 finding that both large sporadics are this one doubling
  - HYP-7440  # the rigid run n=8..12 / sporadic locus
  - THM-733/734  # the GW / near-AP tile
---

# THM-1065 — the doubling family, characterized mod 6

**One line.** Goddyn–Wong's doubling works at `n = 13` and fails at `n = 11` for a
completely explicit reason: the construction leaves exactly one residue class open,
`±(n−1)`, modulo each of the two neighbouring denominators `2n` and `2n+1`, and the
escape witness at that class exists precisely when `n−1` is invertible there.

## Statement

For `n ≥ 5` put `G(n) = \{1,2,…,n−2\} ∪ \{n\} ∪ \{2(n−1)\}` (remove `n−1`, add its
double). Then

```
                 | 1/n           if n is even                     (witness at q = 2n)
   M(G(n))   =   | 2/(2n+1)      if n ≡ 3 or 5 (mod 6)            (witness at q = 2n+1)
                 | 1/(n+1)       if n ≡ 1 (mod 6)      ← TIGHT
```

> **`G(n)` is tight ⟺ `n ≡ 1 (mod 6)` ⟺ `6 | (n−1)`.**

Tight at `n = 7, 13, 19, 25, 31, 37, 43, 49, 55, …` (verified exactly to `n = 60`).

## The mechanism (forward direction PROVED)

Write the residues of `G(n)` in signed form (`|·| ≤ q/2`):

- **mod `2n`:** `G(n) ≡ \{1,…,n−2\} ∪ \{n\} ∪ \{−2\}`, so `±G(n)` covers every class
  except **exactly `±(n−1)`**.
- **mod `2n+1`:** `G(n) ≡ \{1,…,n−2\} ∪ \{n\} ∪ \{−3\}` (since `2n−2 ≡ −3`), so again
  `±G(n)` misses **exactly `±(n−1)`**.

A time `t = m/q` has every clearance `≥ 2/q` iff `vm ≢ 0, ±1 (mod q)` for all `v ∈ G(n)`,
i.e. iff `u := m^{-1}` avoids `±G(n)`. By the two displays the only available class is
`u = ±(n−1)`, so **a clearance-2 witness at `q` exists iff `(n−1)` is invertible mod `q`.**
And

```
gcd(n−1, 2n)   = gcd(n−1, 2n − 2(n−1)) = gcd(n−1, 2),
gcd(n−1, 2n+1) = gcd(n−1, 2n+1 − 2(n−1)) = gcd(n−1, 3).
```

Hence: the `q = 2n` witness exists **iff `n` is even**, giving `M ≥ 2/(2n) = 1/n`; and the
`q = 2n+1` witness exists **iff `3 ∤ (n−1)`**, giving `M ≥ 2/(2n+1)`. Both values exceed
`1/(n+1)`, so in either case `G(n)` is **not tight**. Both witnesses fail simultaneously iff
`2 | (n−1)` and `3 | (n−1)`, i.e. **`6 | (n−1)`**. ∎ (forward direction)

*The `6` is `2·3`, contributed one prime each by the two neighbouring denominators `2n` and
`2n+1`.* Structural claims (missing class, gcd rule, realized witness with clearance ≥ 2)
verified exactly for `n = 5..39`.

## Converse (verified, not proved)

For `6 | (n−1)` the two escapes above are gone and the maximizer returns to `q = n+1` with
clearance 1, i.e. `M = 1/(n+1)`. Verified exactly for every `n ≤ 60` — the three-case `M`
formula matched with **zero** mismatches, and the observed maximizer denominator was `2n`,
`2n+1`, `n+1` in exactly the predicted cases. What is missing for a proof is the statement
that no *other* denominator beats `1/(n+1)` once these two are dead.

## Consequences

1. **An infinite family of sporadic tight instances**: `n = 7, 13, 19, 25, 31, …`, with
   Goddyn–Wong (`n = 13`) as the second member. `n = 7` (`{1,2,3,4,5,7,12}`) is the first.
2. **It provably misses `n = 12`** (`n−1 = 11`, `6 ∤ 11`) — so the most dangerous known
   sporadic family cannot populate the `n=12` branch. This is a genuine, if partial, piece
   of evidence for the `n = 12` rigidity conjecture.
3. It explains part of the sporadic locus `\{4,5,7,13\}`: `7` and `13` are the first two
   members of this family. `n = 4, 5` are **not** (`6 ∤ 3, 4`) and come from a different
   construction — so at least two distinct sporadic mechanisms exist.

## Honest scope

- **Proved:** the forward direction (`6 ∤ (n−1) ⟹ not tight`) with explicit witnesses.
- **Verified, not proved:** the converse, and the exact `M` formula, for `n ≤ 60`.
- This characterizes **one family**. It does not characterize tight instances in general,
  and it does not close the `n = 12` sporadic branch — it only shows this family misses it.

*Artifacts:* computed inline (three-case `M` formula `n=5..60`; structural verification
`n=5..39`). Credits: THM-1031 (covering lemma), THM-1001, HYP-7460 (which isolated this
family as the common form of both large sporadics).
