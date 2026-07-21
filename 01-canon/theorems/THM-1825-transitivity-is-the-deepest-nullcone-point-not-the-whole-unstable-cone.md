---
id: THM-1825
title: "TRANSITIVITY IS THE DEEPEST NULLCONE POINT OF THE CHARACTERISTIC BINARY FORM — but NOT the whole GIT-unstable cone, which leaks at n = 7 (an honest correction of my own S128c131 conjecture). char_A(x) homogenizes to a degree-n binary form in Sym^n(ℂ²); SL₂ acts. (0) ALWAYS: transitive ⟺ char_A = xⁿ = the DEEPEST point of the GIT nullcone (a single root of multiplicity n = the rational-normal-curve vertex; transitive A is strictly triangular = nilpotent, THM-895) — this direction is unconditional and matches death-star-S75, klein-S385, mac-mini-S154, boxeph-S188. (1) FOR n ≤ 6, the UNSTABLE locus (root multiplicity > n/2, Hilbert–Mumford) is EXACTLY the n! transitive tournaments (exhaustive: 6,24,120,720, no non-transitive unstable). (2) REFUTED AT n = 7: the reducible non-transitive tournament with score sequence [0,2,2,2,4,5,6] (one 3-cycle, c₃ = 1) has char_A = x⁴(x−1)(x²+x+1), so root 0 has multiplicity 4 > n/2 = 3.5 — GIT-UNSTABLE but NOT transitive. Its transitive backbone is a single nilpotent Jordan block of size 4 (geometric mult of 0 is 1, algebraic 4). So 'unstable = transitive' is a small-n coincidence that BREAKS at n = 7, another n ≥ 7 phase transition. (3) TWO LEMMAS, PROVED, that characterise the leak: LEMMA A — a root of multiplicity μ > n/2 is an INTEGER eigenvalue (its minimal polynomial f satisfies fᵘ | char_A so deg f · μ ≤ n, forcing deg f = 1, and a rational eigenvalue of an integer matrix is an integer); the counterexample's is 0. LEMMA B — since (A−λI)+(A−λI)ᵀ = J−(1+2λ)I, and rank M ≥ ½ rank(M+Mᵀ), the geometric multiplicity g(λ) = n − rank(A−λI) ≤ ⌊n/2⌋ for every λ ∉ {(n−1)/2, −½} (in particular rank A ≥ ⌈n/2⌉); so the excess multiplicity of an unstable eigenvalue is ALWAYS Jordan structure, never a big eigenspace. The counterexample confirms: geometric mult 1, algebraic 4, one Jordan block. (4) UNCHANGED: tr(Aᵏ) = power sums = SL₂-invariants (the moment ladder is the Sym^n coefficient map); the fibers of T ↦ char_A are co-spectral classes, where H splits (THM-1780); char_S is the even companion, nullcone on the A-side only (tr S² = −n(n−1); the half-dictionary ½ moves off the Weyl axis). The Vandermonde SURVIVORS (klein-S385) = transitive is a DIFFERENT, exact-for-all-n statement about the √-discriminant covariant — not the stability cone, which is why it does not leak"
status: >
  RENUMBERED THM-1805 -> THM-1810 -> THM-1825 (TWO first-pusher bumps: klein-S385 Vandermonde took 1805; opus-S434 subquestions took 1810). Originally from THM-1805 (klein-S385's Vandermonde THM-1805 first-pushed 20:32:03, mine
  20:32:22 — 19 s later).
  (0) PROVED (transitive ⟺ nilpotent ⟺ char_A = xⁿ, classical + THM-895).
  (1) VERIFIED exhaustively n = 3..6 (unstable count = n!, all acyclic).
  (2) REFUTED at n = 7 by an EXACT witness (integer char poly x⁴(x−1)(x²+x+1), c₃ = 1,
  strongly-connected = False) — this CORRECTS the "conjecture for n ≥ 7" that my S128c131
  version hedged; the conjecture is false.  The numpy sweep that first flagged it also threw
  two clustering artifacts (true mult 2), caught by exact factoring — recorded so the method
  is trusted only after the exact check.
  (3) LEMMA A and LEMMA B are PROVED (the arguments are complete and elementary).  They do
  not prove "unstable = transitive" (which is false); they characterise the unstable locus as
  {integer eigenvalue with a Jordan block of total size > n/2}, and the deepest such point is
  transitive.
  (4) Unchanged identifications, exact.
  Honest: the headline of the S128c131 version ("unstable = exactly transitive") is downgraded
  to "transitive = the deepest unstable point; unstable = transitive for n ≤ 6 only."
source: kind-pasteur-2026-07-20-S128c132 (owner: work the clean next steps of THM-1805/1810)
depends_on:
  - THM-895     # λ = 0 ⟺ transitive
  - THM-1775    # the moment-nullcone template
  - THM-1780    # H leaves the ladder at n=6 (the fiber forgets the permanent)
  - THM-1555    # the half-dictionary
related: [THM-1805, THM-1725]   # THM-1805 = klein's Vandermonde (the covariant companion)
script: 04-computation/tournament_binary_form_git_kps_S128c131.py, unstable_is_transitive_proof_kps_S128c132.py (+ .out)
---

# THM-1825 — transitivity is the deepest nullcone point, not the whole unstable cone

The clean next step from the binary-form frame was to prove "the GIT-unstable tournaments are
exactly the transitive ones" for all `n`. **It is false — it breaks at `n = 7`** — and finding
exactly where, with two proved lemmas that explain it, is the result.

## (0) Always: transitive is the deepest nullcone point

`char_A(x)` homogenizes to a degree-`n` binary form. Transitive `A` is strictly triangular,
hence nilpotent, so `char_A = xⁿ` — a single root of multiplicity `n`, the **deepest point of
the GIT nullcone** (the vertex of the rational normal curve). This direction is unconditional
and is the shared content of the day's convergence (death-star-S75 "rational-normal-curve
vertex", klein-S385, mac-mini-S154, boxeph-S188 Kempf–Ness).

## (1)–(2) The unstable cone: transitive for `n ≤ 6`, larger from `n = 7`

By Hilbert–Mumford a binary `n`-ic is **unstable** iff it has a root of multiplicity `> n/2`.

- **`n ≤ 6`:** the unstable tournaments are *exactly* the `n!` transitive ones (exhaustive
  census: `6, 24, 120, 720`, none non-transitive). So "unstable = transitive" holds — but only
  here.
- **`n = 7`: it fails.** The reducible tournament with score sequence `[0,2,2,2,4,5,6]` (a
  source, a sink, one embedded 3-cycle, `c₃ = 1`) has

  > `char_A = x⁴(x−1)(x²+x+1)`,

  so `0` is a root of **multiplicity 4 > n/2 = 3.5**: **GIT-unstable, but not transitive.** Its
  transitive backbone contributes a single nilpotent Jordan block of size 4 (`x⁴`), while the
  3-cycle contributes `(x−1)(x²+x+1)`. So the unstable cone strictly contains the transitive
  locus from `n = 7` on — one more entry in the project's `n ≥ 7` phase-transition ledger.

(The numpy sweep that surfaced this also produced two false "multiplicity ≥ 4" flags that were
really `(x−1)²(x²+x+1)²`-type clusters of true multiplicity 2; exact factoring caught them.
Method note: trust the cluster only after the exact char-poly check.)

## (3) Two proved lemmas, and they explain the leak

**Lemma A — an unstable eigenvalue is an integer.** If `λ` is a root of `char_A` of
multiplicity `μ > n/2`, its minimal polynomial `f` over `ℚ` is irreducible with `f^μ | char_A`,
so `deg(f)·μ ≤ n`, forcing `deg f = 1`; `λ ∈ ℚ`, and a rational eigenvalue of an integer matrix
is an integer. ∎ (Here `λ = 0`.)

**Lemma B — geometric multiplicity `≤ ⌊n/2⌋` off the Perron and `−½` values.** Because
`(A−λI)+(A−λI)ᵀ = A+Aᵀ − 2λI = J − (1+2λ)I`, and `rank M ≥ ½·rank(M+Mᵀ)`,

> `rank(A−λI) ≥ ½·rank(J − (1+2λ)I) = n/2`  for `λ ∉ {(n−1)/2, −½}`

(`J − (1+2λ)I` has full rank `n` unless `1+2λ ∈ {0, n}`). Hence `g(λ) = n − rank(A−λI) ≤
⌊n/2⌋`, and in particular `rank A ≥ ⌈n/2⌉`. ∎

**Together they pin the mechanism.** An unstable eigenvalue is an integer (Lemma A), so `λ ≠
−½`; off the Perron value its geometric multiplicity is `≤ ⌊n/2⌋ < μ` (Lemma B), so the excess
`μ − g` is **Jordan structure** — never a large eigenspace. The `n = 7` witness realises the
extreme: geometric multiplicity `1`, algebraic `4`, one Jordan block. So the unstable locus is

> `{` integer eigenvalue carrying a Jordan block of total algebraic size `> n/2 }`,

and its **deepest point** (a single block of size `n` at `0`) is the transitive tournament. The
lemmas turn "unstable = transitive" from a false equality into a correct *description* of a cone
whose vertex is transitivity.

## (4) What still stands

`tr(Aᵏ) = ` power sums `= ` SL₂-invariants; the moment-nullcone ladder is the `Sym^n`
coefficient map (THM-1775); `T ↦ char_A` has co-spectral fibers, where `H` splits at `n = 6`
(THM-1780); `char_S` is the even (skew) companion with the nullcone on the `A`-side only
(`tr S² = −n(n−1) ≠ 0`; the half-dictionary `½` moves off the Weyl axis, THM-1555). None of
these depend on the false converse.

**The Vandermonde companion (klein-S385, THM-1805).** `∏(xᵢ−xⱼ) = Σ_T sgn(T)·x^{score(T)}`, and
the **surviving** terms are exactly the `n!` transitive tournaments (permutation score
sequences), for *all* `n`. That is a statement about the `√`-discriminant *covariant*, not the
stability *cone* — which is precisely why it does not leak at `n = 7`: survival (a permutation
score sequence) is a sharper condition than instability (a root of multiplicity `> n/2`). The
two "transitive is special" facts are genuinely different, and only the coarser one (stability)
degenerates.

## Named next

- **Characterise the `n ≥ 7` unstable non-transitive tournaments.** The witness is reducible
  (a dominating and a dominated vertex). Conjecture: the unstable non-transitive tournaments are
  exactly those with a large transitive "core" contributing a nilpotent Jordan block of size
  `> n/2` — i.e. reducible tournaments whose condensation is near-total.
- **The strictly-semistable stratum** (max root multiplicity `= n/2`, `n` even; `960` at
  `n = 6`) is still worth identifying.
- **Turn Lemma A/B into a stability criterion**: `char_A` is semistable iff every integer
  eigenvalue's Jordan blocks total `≤ n/2` — a finite, checkable condition.
