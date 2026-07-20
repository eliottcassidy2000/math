---
id: THM-1455
title: "THE ODD/EVEN ↔ sin/cos LINK IS THE PFAFFIAN, and it yields a proved mod-16 law for odd n. (0) THE NEXT STEP, done: three vertex-disjoint odd cycles need 9 vertices, so α₃ = 0 for n ≤ 8 and hp = 1 + 2α₁ + 4α₂ EXACTLY through n = 8 — verified, extending THM-1445's n ≤ 5 regime by three. (I) THE MECHANISM: for the skew adjacency S (S_ij = ±1), odd-order principal minors VANISH identically and even-order ones enter as SQUARED PFAFFIANS, coeff of x^{n−2k} = Σ over 2k-subsets of Pf(S_A)² — verified exactly, 0 mismatches. So char(S) is parity-pure: ODD for n odd, EVEN for n even, with purely imaginary spectrum, hence char(S)(iy) is real-rooted. That is the sin/cos content of the odd/even theme — not a metaphor but the spectral parity of the skew adjacency. (II) THE OWNER'S POLYNOMIAL IS REALISED: x(x²+7)(x⁴+14x²+17) = x⁷+21x⁵+115x³+119x occurs as char(S) of an explicit 7-vertex tournament (hp = 91); Paley T₇ gives x(x²+7)³ = x⁷+21x⁵+147x³+343x instead. The shared 21 is forced — it is C(7,2) for every tournament — so the discriminating coefficients are x³ and x. (III) A PROVED PARITY LAW: writing c_{n−4} = C(n,4) + 8k₄ with k₄ = #{4-subsets with |Pf| = 3}, k₄ is EVEN whenever n is ODD, hence c_{n−4} ≡ C(n,4) (mod 16). Proof: exhaustive base case over all 1024 labelled 5-tournaments (k₄ ∈ {0,2} always), then double counting — each 4-subset lies in exactly n−4 five-subsets, so (n−4)·k₄ is even, and n odd makes n−4 odd. FAILS for even n, with 16 of the 64 labelled 4-tournaments having k₄ = 1"
status: >
  (0) PROVED by the vertex count (α₃ = 0 needs 3+3+3 = 9 > 8); verified at n = 6, 7 with
  zero violations.
  (I) VERIFIED-EXACT: the Pfaffian minor expansion checked at n = 5, 6, 7 over 40
  tournaments each and all k, zero mismatches.  The parity statement itself is classical
  (real skew-symmetric spectra are purely imaginary in ± pairs); what is checked here is
  that it is the RIGHT reading of the owner's odd/even question.
  (II) VERIFIED-EXACT with an explicit witness tournament, found by sampling 80,000
  tournaments on 7 vertices.
  (III) PROVED.  The base case is exhaustive (all 2^10 labelled 5-tournaments), the
  double-counting step is elementary, and the even-n failure is exhibited exhaustively at
  n = 4.  Verified independently by direct sampling at n = 5,7,9,11 (always even) against
  n = 4,6,8,10 (mixed).
  Nothing here advances a named open problem.
source: kind-pasteur-2026-07-20-S128c115 (owner: work the next step, and see how odd/even relates to sin/cos and x(x²+7)(x⁴+14x²+17))
depends_on:
  - THM-1445    # the n <= 5 exact OCF regime this extends
related: [THM-165, THM-159]
script: 04-computation/skew_charpoly_and_alpha3_kps_S128c115.py, pfaffian_parity_law_kps_S128c115.py (+ .out)
---

# THM-1455 — the Pfaffian is the odd/even mechanism

## 0. The named next step

Three vertex-disjoint odd cycles need `3+3+3 = 9` vertices, so `α₃ = 0` for every
tournament with `n ≤ 8`, and therefore

> **`hp(T) = 1 + 2α₁ + 4α₂` exactly, for all `n ≤ 8`** — verified at `n = 6, 7`, zero
> violations, extending THM-1445's `n ≤ 5` regime by three.

## I. Where sin/cos actually enters

Give a tournament its **skew** adjacency `S` (`S_ij = +1` if `i→j`, `−1` if `j→i`, `0` on
the diagonal). Then

> odd-order principal minors of a skew matrix **vanish identically**, and even-order ones
> are **squared Pfaffians**:
> `[x^{n−2k}] char(S) = Σ_{|A| = 2k} Pf(S_A)²`.

Verified exactly at `n = 5,6,7`, all `k`, zero mismatches. Consequently `char(S)` is
**parity-pure** — odd for odd `n`, even for even `n` — and its spectrum is purely
imaginary in `±` pairs, so `char(S)(iy)` is a **real-rooted** odd (resp. even) polynomial
in `y`.

That is the honest content of "odd/even ↔ sin/cos" here: sin and cos are the archetypes of
odd and even functions with all zeros real, and the skew characteristic polynomial of a
tournament is exactly such an object, with the parity dictated by `n`. It is a spectral
statement, not an analogy. **The odd-order minors vanishing is the whole mechanism** — it
is why there is no even part when `n` is odd.

## II. The owner's polynomial

`x(x²+7)(x⁴+14x²+17) = x⁷ + 21x⁵ + 115x³ + 119x` is odd of degree 7 with all roots purely
imaginary — the right shape for `char(S)` at `n = 7`. **It is realised.** Explicit witness
(arcs listed in the output), with `hp = 91`.

For contrast, Paley `T₇` has skew spectrum `{0, ±i√7 thrice}` and so
`char = x(x²+7)³ = x⁷ + 21x⁵ + 147x³ + 343x`. Same `x(x²+7)(x⁴+14x²+c)` shape, `c = 17`
against `49`: **the owner's polynomial is realised by a non-Paley tournament.**

The shared `21` carries no information. For skew `S`, `tr S = 0` and
`tr S² = −n(n−1)`, so `[x^{n−2}] char(S) = C(n,2)` for **every** tournament — 21 at
`n = 7`, always. The discriminating coefficients are `x³` and `x`.

A coarseness note worth recording: sampling 80,000 tournaments on 7 vertices produced only
**11 distinct** skew characteristic polynomials, against 456 isomorphism classes. As an
invariant this is very weak.

## III. A proved mod-16 law, and it holds only for odd n

`Pf` of a `4×4` `±1` skew matrix is `S₁₂S₃₄ − S₁₃S₂₄ + S₁₄S₂₃`, a sum of three `±1`
terms, so `|Pf| ∈ {1,3}` and `Pf² ∈ {1,9}`. Hence

> `c_{n−4} = C(n,4) + 8k₄`,  `k₄ := #{4-subsets with |Pf| = 3}` — verified at `n = 6,7,8`.

> **Theorem.** If `n` is **odd** then `k₄` is even, so `c_{n−4} ≡ C(n,4) (mod 16)`.

*Proof.* **Base case**, exhaustive over all `2¹⁰ = 1024` labelled 5-vertex tournaments:
`k₄ ∈ {0, 2}` in every case (384 with 0, 640 with 2), so `k₄` is even for `n = 5`.
**Double counting**: each 4-subset of an `n`-set lies in exactly `n − 4` five-subsets, so
`(n−4)·k₄(T) = Σ_{|B| = 5} k₄(T[B])`. Every term on the right is even by the base case,
so `(n−4)·k₄(T)` is even. If `n` is odd then `n − 4` is odd, forcing `k₄(T)` even. ∎

**It fails for even `n`, and must**: `n − 4` is then even and the argument gives nothing.
Exhaustively at `n = 4`, 16 of the 64 labelled tournaments have `k₄ = 1`. Direct sampling
agrees across the board — always even at `n = 5, 7, 9, 11`; mixed at `n = 4, 6, 8, 10`.

At `n = 7` this gives `c₃ ≡ 35 ≡ 3 (mod 16)`, which is exactly the congruence the attained
coefficients `35, 67, 83, 99, 115, 131, 147` display — including the owner's `115`.

The `6`-subset analogue is weaker but real: a `6×6` `±1` skew Pfaffian is a sum of 15 `±1`
terms, hence **odd**, so `Pf² ≡ 1 (mod 8)` and `c_{n−6} ≡ C(n,6) (mod 8)` — verified,
residual 0.

## Named next

- The same double-counting scheme should give the `2k`-subset law in general: if the
  `(2k+1)`-vertex base case has `#{|Pf| ≡ specified}` even, then odd `n − 2k` propagates
  it. Worth one pass at `k = 3` (base case `n = 7`, 2²¹ tournaments — feasible).
- `α₃ = 0` for `n ≤ 8` extends to `α_k = 0` for `n < 3k`; the exact OCF regime is
  `⌊n/3⌋` terms, so the next threshold after 8 is `n = 11`.
