---
id: THM-1065
title: The doubling family is tight EXACTLY when n ≡ 1 (mod 6) — for G(n) = {1,…,n−2, n, 2(n−1)} (remove n−1, add its double), M(G(n)) = 1/n if n is even, 2/(2n+1) if n ≡ 3,5 (mod 6), and 1/(n+1) — TIGHT — iff n ≡ 1 (mod 6). The mechanism is exact: G(n) misses precisely the residue class ±(n−1) modulo BOTH 2n and 2n+1, so a clearance-2 witness exists at q iff (n−1) is invertible mod q; and gcd(n−1,2n) = gcd(n−1,2) while gcd(n−1,2n+1) = gcd(n−1,3). Both witnesses die exactly when 2 | (n−1) AND 3 | (n−1), i.e. 6 | (n−1). THE "6" IS 2·3 FROM THE TWO NEIGHBOURING DENOMINATORS. This yields an infinite family of sporadic tight instances (n = 7,13,19,25,31,…), contains Goddyn–Wong (n=13) as its second member, and PROVABLY MISSES n=12.
status: **REDISCOVERY — the characterization is the r=n−1, m=2 slice of the PUBLISHED Goddyn–Wong criterion (see the S120 attribution note below); NOT new.** The proof given here is independent and elementary, and the exact M-values for the non-tight cases ARE net-new. Both directions proved unconditionally (S119 supplies the converse). Forward: explicit escape witnesses at q = 2n, 2n+1. Converse: a Farey ring argument using only the classical adjacency fact (adjacent denominators in F_N sum to > N) — no appeal to LRC. Verified: the three-case M formula exact for n = 5..60; the converse's neighbour criterion matches 6 | (n−1) exactly for n = 5..49.
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

> ## ⚠ ATTRIBUTION CORRECTION (mac-mini-S120) — this result is NOT new
>
> Corpus mining (THM-709's addendum, klein-S253's literature merge of Perarnau–Serra
> §4) shows the **published Goddyn–Wong acceleration criterion**:
>
> > replacing `r ∈ [n]` by `mr` preserves tightness **iff `r` has a common factor with
> > every element of the interval `[(n+1−r), m(n+1−r)−1]`**.
>
> For this family (`r = n−1`, `m = 2`) that interval is
> `[n+1−(n−1), 2(n+1−(n−1))−1] = [2,3]` for **every** `n`, so the criterion reads
> *"`n−1` shares a factor with both 2 and 3"* — i.e. exactly **`6 | (n−1)`**. Verified
> identical to this theorem's prediction for all `n < 50`.
>
> **So the characterization below is the `r=n−1, m=2` slice of GW, not a new theorem,
> and the "2·3" I presented as a discovery is literally GW's interval `[2,3]`.** I
> claimed novelty for it in the S118/S119 letters; that was wrong and is retracted.
>
> **What remains net-new here:** (i) an independent, elementary, unconditional proof
> (escape-witness + Farey ring) that never invokes LRC; (ii) the *exact* `M` values in
> the non-tight cases (`1/n` for even `n`, `2/(2n+1)` for `n ≡ 3,5 mod 6`) — GW's
> criterion is a tight/not-tight predicate and does not supply these; (iii) the
> mechanism, i.e. *why* the interval is `[2,3]` (below). This mirrors THM-709, whose
> stated net-new content over GW was likewise "the exact escape values".
>
> Also from the same merge, for the height-bound thread: **Erdős/Jacobsthal-linked tight
> families with `v_max = 2n − Θ(log n)` exist**, and **Pomerance**: `n < v_max < 2n −
> c·log²n ⟹ NOT tight`. Both are consistent with — and sharpen — the `max(A) ≤ 2n`
> conjecture (HYP-7450): the true extremal growth is `2n − Θ(log n)`.

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

## Converse (PROVED — the Farey ring argument)

We must show: `6 | (n−1) ⟹` the danger sets of `G(n)` at level `L = 1/(n+1)` cover the circle.

**Step 1 (reduction to a ring).** The AP `\{1,…,n\}` already covers, *elementarily*: by the
classical Farey fact that adjacent denominators in `F_N` satisfy `i+j > N`, every gap remnant
`(1 − L(i+j))/(ij)` is `≤ 0` at `L = 1/(n+1)`. So for any `t` some `v ∈ \{1,…,n\}` has
`‖vt‖ ≤ L`; if `v ≠ n−1` then `v ∈ G(n)` and we are done. Hence `G(n)` can only fail at times
covered **solely** by `n−1`. Such `t` lie in an arc of `n−1`: write `t = a/(n−1) + ε` with
`(n−1)|ε| ≤ L`. (If `a/(n−1)` is *not* primitive, its reduced denominator `m' ≤ (n−1)/2 ≤ n−2`
is itself a speed of `G(n)` with the larger radius `L/m'`, which swallows the whole arc — so
only **primitive** `a/(n−1)` matter.)

**Step 2 (the doubled speed covers the inside).** `‖(2n−2)t‖ = ‖2(n−1)t‖ = 2(n−1)|ε|`, which is
`≤ L` exactly when `(n−1)|ε| ≤ L/2`. So `2n−2` covers the inner half, leaving the

> **dangerous ring** `L/(2(n−1)) < |ε| ≤ L/(n−1)`.

**Step 3 (when a Farey neighbour covers the ring).** Let `p/i` be a Farey neighbour of
`a/(n−1)` in `F_n`; by unimodularity it sits at distance `1/(i(n−1))`, and its arc has radius
`L/i` (legitimate: `i ≤ n−2` or `i = n`, both speeds of `G(n)`). Its arc reaches the ring iff

```
1/(i(n−1)) − L/i  ≤  L/(2(n−1)).
```

Multiplying by `i(n−1)` and using `1 − L(n−1) = 1 − (n−1)/(n+1) = 2/(n+1) = 2L`, this becomes
`2L ≤ iL/2`, i.e. **`i ≥ 4`**. (The `2n−2`-centres other than `a/(n−1)` itself sit at distance
`≥ 1/(2(n−1))` with radius `L/(2(n−1))`, hence reach only to `(1−L)/(2(n−1)) > L/(n−1)` — they
never enter the ring, so they cannot help.)

**Step 4 (small denominators are exactly 2 and 3).** So `G(n)` is tight iff every primitive
`a/(n−1)` has **both** `F_n`-neighbours of denominator `≥ 4`. Now:

- **denominator 1** — the neighbours of `0/1, 1/1` in `F_n` are `1/n, (n−1)/n`, never `a/(n−1)`;
- **denominator 2** — adjacency to `1/2` forces `|2a − (n−1)| = 1`, hence `n−1` **odd**;
- **denominator 3** — adjacency to `1/3` or `2/3` forces `|3a − (n−1)| = 1` or
  `|3a − 2(n−1)| = 1`, each impossible when **`3 | (n−1)`**.

(Both are genuine adjacencies when the congruence permits, since `2 + (n−1) > n` and
`3 + (n−1) > n`.) Therefore both small-denominator neighbours are excluded **iff** `2 | (n−1)`
and `3 | (n−1)`, i.e. **`6 | (n−1)`**, and then the ring is covered and `G(n)` is tight. ∎

Verified: the neighbour criterion matches `6 | (n−1)` exactly for `n = 5..49`, and the killing
neighbour is precisely denominator `3` when `2|(n−1)` only, denominator `2` when `3|(n−1)` only,
and both when neither divides.

## The 2·3 appears twice — and both are GW's interval [2,3] (S120)

The "coincidence" **dissolves completely**, and in the end all three views are the same
published fact:

- **Forward:** an escape at denominator `q` exists iff `gcd(n−1, q) = 1`.
- **Converse:** a Farey neighbour of `a/(n−1)` at denominator `i` exists iff the
  unimodular equation `|ia − (n−1)p| = 1` is solvable, i.e. iff `gcd(i, n−1) = 1`.

**These are the same condition** — `gcd(·, n−1) = 1`. And they are applied to the same
numbers, because `2n = 2(n−1)+2` and `2n+1 = 2(n−1)+3`: *the escape denominators reduce
mod `(n−1)` exactly onto the Farey denominators 2 and 3.* The range is forced too — a
clearance-2 escape needs `2/q > 1/(n+1)`, so `q ≤ 2n+1`, and `q = 2(n−1)+i` gives
`i ∈ {1,2,3}`; `i = 1` (`q = 2n−1`) is self-defeating because there `n ≡ −(n−1)`, so the
speed `n` itself fills the class `n−1` and no class is left open. Hence `i ∈ {2,3}`.

And `\{2,3\}` is precisely GW's interval `[(n+1−r), m(n+1−r)−1] = [2,3]`. So: one
arithmetic condition, three coordinate systems (escape denominators, Farey denominators,
GW's interval).

## (superseded framing) The 2·3 appears twice, from opposite directions

This is the striking part. The same factorisation `6 = 2·3` controls both halves, by
different mechanisms:

| | source of the 2 | source of the 3 |
|---|---|---|
| **forward** (escape exists) | `gcd(n−1, 2n) = gcd(n−1,2)` | `gcd(n−1, 2n+1) = gcd(n−1,3)` |
| **converse** (ring uncovered) | Farey neighbour `1/2` | Farey neighbour `1/3, 2/3` |

The forward direction reads `6` off the two **neighbouring denominators** `2n, 2n+1`; the
converse reads the same `6` off the two smallest **Farey denominators** `2, 3`. 

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
