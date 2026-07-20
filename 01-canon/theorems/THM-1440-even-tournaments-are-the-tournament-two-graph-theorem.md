---
id: THM-1440
title: "EVEN TOURNAMENTS ARE THE TOURNAMENT TWO-GRAPH THEOREM — the analogue kind-pasteur's THM-1415 declared nonexistent does exist; it was compared against the wrong object. A tournament IS a SKEW-Seidel matrix, so the analogue of an even GRAPH (all degrees even) is an even TOURNAMENT (all SCORES even), not an even graph. THE SCORE-PARITY LAW: switching at U changes s_v by |Uᶜ| mod 2 for v ∈ U and |U| mod 2 for v ∈ Uᶜ, so for n ODD the score-parity vector flips exactly on the EVEN-SIZED member of {U, Uᶜ}, and U ↦ that set is a BIJECTION from the 2^{n−1} switchings onto the even-weight code. Since Σs_v = C(n,2): at n ≡ 1 (mod 4) every switching class contains a UNIQUE tournament with all scores even; at n ≡ 3 (mod 4) none, but exactly n with a single odd score; at n even the parity vector is pinned only up to global complement. COROLLARY: #even tournaments up to iso = #switching classes up to iso = A049313(n) for n ≡ 1 mod 4 — verified as a BIJECTION (not merely matching counts) at n = 5, and predicting 792 at n = 9. SECOND COROLLARY, touching a stated open problem: every automorphism of a switching class fixes that unique even member, so Babai–Cameron's Remark 7.4 failure set is EMPTY at n ≡ 1 mod 4. AND THE ODD-VALUED SIDE: det S(T) = 0 for n odd, while for n even det S = Pf(S)² with Pf(S) always ODD — because S ≡ J − I (mod 2) and Pf(J−I) = (n−1)!!, a product of odd numbers. A companion to Rédei's hp(T) odd."
status: >
  PROVED (the parity law and the bijection are three short lemmas, given below) and
  VERIFIED EXHAUSTIVELY.  n = 5: all 64 labelled switching classes, each of size 16,
  each with EXACTLY ONE all-even member, and the induced map on iso classes checked to be
  a BIJECTION 2 <-> 2 rather than a coincidence of counts.  n = 7: exactly 7 single-odd
  representatives per class, 64 = 2^(n-1) distinct parity vectors, zero all-even -- as
  predicted for 3 mod 4.  n = 6: exactly 2 parity vectors per class, as predicted for n even.
  Switching-class counts recomputed independently here: 1, 2, 2, 6, 12 at n = 3..7,
  reproducing A049313 (Babai-Cameron) exactly.
  Pfaffian parity: proved in one line, verified at n = 4, 6, 8.
  HONEST: the argument is ELEMENTARY and may well be folklore.  A literature search
  (OEIS, Babai-Cameron, Seidel, Mallows-Sloane, Hage-Harju-Welzl) found NO published
  statement of it; "even tournament" returns no OEIS results.  Absence of a found source
  is not a priority claim.
source: klein-2026-07-20-S338 (owner: think skew-Seidel and chase the high leverage question; see the relation between odd valued functions and tournament adjacent ideas; they both relate also to even concepts like even graphs and even functions)
depends_on:
  - THM-474   # mac-mini: tilings ARE switching classes of tournaments (the Gauge Theorem)
corrects:
  - THM-1415  # kind-pasteur: "the graph two-graph theorem has no tournament analogue" -- wrong comparand
external:
  - "Seidel (1974): a graph of ODD order has a UNIQUE even graph in its switching class."
  - "Mallows-Sloane (SIAM J. Appl. Math. 28, 1975): two-graphs = switching classes = Euler graphs, A002854."
  - "Babai-Cameron, EJC 7 (2000) #R38: switching classes of tournaments = A049313 = 1,1,1,2,2,6,12,79,792,...; Remark 7.4 poses the failure-set enumeration and states 'We cannot do this'."
  - "Moon: the automorphism group of a tournament has ODD order."
script: 04-computation/skew_seidel_even_tournaments_klein_S338.py (+ .out)
---

# THM-1440 — even tournaments, and the odd-valued side

## 0. The frame the owner asked for

A tournament **is** a **skew**-Seidel matrix: `S_ij = +1` if `i→j`, `−1` if `j→i`, `0` on the
diagonal, so `Sᵀ = −S`. A graph is a **symmetric** Seidel matrix. Seidel switching
`S ↦ D_U S D_U` acts on both and preserves the symmetry type. So:

| | symmetric / **even** side | skew / **odd** side |
|---|---|---|
| object | graphs | **tournaments** |
| switching classes | two-graphs, `A002854` | oriented two-graphs, `A049313` |
| canonical representative | **even graph** (all degrees even), `n` odd — *Seidel 1974* | **even tournament** (all scores even), `n ≡ 1 mod 4` — **§2** |
| determinant | — | `det S = 0` (`n` odd); `Pf(S)²`, **`Pf` odd** (`n` even) — **§4** |

## 1. The score-parity law

**Lemma 1.** Switching at `U` changes the score by `Δs_v ≡ |Uᶜ| (mod 2)` for `v ∈ U`, and
`Δs_v ≡ |U| (mod 2)` for `v ∈ Uᶜ`.

*Proof.* For `v ∈ U`, exactly its `|Uᶜ|` arcs to `Uᶜ` reverse; if `a` were out-arcs, the new
count is `|Uᶜ| − a`, so `Δs_v = |Uᶜ| − 2a`. ∎

**Lemma 2 (`n` odd).** Exactly one of `|U|, |Uᶜ|` is even; call it `W(U)`. Then switching at
`U` flips the score-parity vector **exactly on `W(U)`**.

*Proof.* `n` odd ⟹ `|U|, |Uᶜ|` have opposite parities. If `|U|` is even then `|Uᶜ|` is odd, so by
Lemma 1 the vertices of `U` flip and those of `Uᶜ` do not; flipped set `= U = W(U)`. Symmetric
in the other case. ∎

**Lemma 3.** `U ↦ W(U)` is a **bijection** from switchings (subsets modulo `U ∼ Uᶜ`, of which
there are `2^{n−1}`) onto the even-sized subsets of `[n]` (also `2^{n−1}` for `n` odd).

*Proof.* `W(U) = W(Uᶜ)`, so it is well defined on the quotient, and `W(U)` determines
`{U, Uᶜ}`, so it is injective; the two sets have equal size. ∎

> **Consequence.** For `n` odd the `2^{n−1}` members of a switching class realise `2^{n−1}`
> **distinct** score-parity vectors, namely the whole coset `p₀ + 𝓔` where `𝓔` is the
> even-weight code. So the achievable parity vectors are *exactly* those of weight
> `≡ |p₀| (mod 2)`, each attained **once**.

## 2. The trichotomy

`|p₀| = #\{v : s_v \text{ odd}\} ≡ Σ_v s_v = C(n,2) \pmod 2`, and so:

| `n` | `C(n,2)` | canonical representative |
|---|---|---|
| `≡ 1 (mod 4)` | even | **a UNIQUE tournament with all scores even** |
| `≡ 3 (mod 4)` | odd | none all-even; **exactly `n`** with a single odd score (one per vertex) |
| even | — | `\|U\|,\|Uᶜ\|` share parity ⟹ either no score flips or all do: the parity vector is pinned only up to global complement. **No analogue.** |

**Verification.** `n = 5`: all 64 labelled switching classes, each of size 16, each with
**exactly one** all-even member. `n = 7`: 64 distinct parity vectors per class, **zero**
all-even, **exactly 7** single-odd. `n = 6`: exactly **2** parity vectors per class.

## 3. The count — and the correction

Since relabelling preserves "all scores even" and commutes with switching, the bijection of §2
is `S_n`-equivariant and descends:

> **Corollary.** For `n ≡ 1 (mod 4)`, `#\{`even tournaments`\}/≅ = #\{`switching classes`\}/≅
> `= A049313(n)`.

`n = 5`: **2 = 2**, and checked as an actual **bijection** on iso classes, not merely matching
counts. `n = 9`: predicts **792**. Switching-class counts recomputed here independently —
`1, 2, 2, 6, 12` at `n = 3..7` — reproduce `A049313` exactly.

**What this corrects.** kind-pasteur's THM-1415 §II compared `1,2,2,6` against
`A002854 = 2,3,7,16` and concluded *"the graph two-graph theorem has no tournament analogue by
this route"*, recording it "because the analogy is attractive and someone will try it again."
The analogy is fine; the **comparand** was wrong. The two-graph theorem says *switching classes
of graphs ↔ even graphs*; its tournament analogue must be *switching classes of tournaments ↔
even **tournaments***. Comparing tournament switching classes to even **graphs** compares two
different categories. That file also asks for an OEIS check on `1,2,2,6`: it is **A049313**
(Babai–Cameron), so the sequence is known and should not be treated as new.

**A stated open problem, partly answered.** Babai–Cameron Remark 7.4 (verbatim): *"every
automorphism of a switching class of graphs fixes some graph in that class … It would be
interesting to enumerate the switching classes of tournaments for which this fails … **We
cannot do this**."* At `n ≡ 1 (mod 4)` the answer is **zero**: any `g ∈ Aut(Σ)` sends the unique
all-even member to an all-even member of `Σ`, hence to itself. Verified at `n = 5`:
`Aut(Σ) = Aut(T₀)` in every class, every automorphism fixes `T₀`, and every order is odd (3 or
5) — consistent with the classical fact (Moon) that tournament automorphism groups have odd
order, which makes the oddness automatic rather than accidental.

## 4. The odd-valued side

The owner's "odd valued functions" and "even concepts" meet here:

- **`n` odd:** a skew-symmetric matrix of odd order is singular, so `det S(T) = 0` for
  **every** tournament on an odd number of vertices. *(Verified `n = 3,5,7`.)*
- **`n` even:** `det S = Pf(S)²`, and **`Pf(S)` is always ODD**.

  *Proof.* Every off-diagonal entry is `±1 ≡ 1 (mod 2)`, so `S ≡ J − I (mod 2)`; hence
  `Pf(S) ≡ Pf(J−I) = #\{`perfect matchings of `K_n\} = (n−1)!!`, a product of odd numbers. ∎
  *(Verified `n = 4,6,8`: `Pf ∈ {±1,±3,±5,…}`, `det ∈ {1,9,25,49,81,121,…}`.)*

So the project's parity mandate has **three** faces, all forced and all odd-valued:
`hp(T)` odd (Rédei, all `n`), `Pf(S(T))` odd (`n` even), and the vanishing `det S(T) = 0`
(`n` odd) that is the same odd/even split which makes the score-parity law of §1 work only at
odd `n`.

## 5. Scope

Elementary throughout, and possibly folklore — the literature search found no statement of §2,
but absence of a found source is not a priority claim. Nothing here is asymptotic and nothing
touches LRC. The concrete next test is `n = 9`: **792** even tournaments up to isomorphism.

*Files: `04-computation/skew_seidel_even_tournaments_klein_S338.py` (+ `.out`).*
