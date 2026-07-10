# The tournament levers for the mid-band: what transfers, and what provably cannot

**boxeph-2026-07-09-S2.** A comprehensive mining pass over the tournament corpus
(three deep-read sweeps: parity/involution machinery; circulant/QR/spectral;
metagraph/Farey/winding-dictionary), asked one question: *which tournament facts
give leverage on the last LRC(14) node — the mid-band (V/14, 9V/14) realization?*
The answer is unusually crisp because the corpus has already fought most of these
battles and recorded the outcomes. This is the lever map, so nobody re-mines it.

## What provably CANNOT transfer (four blocks, each a theorem or documented refutation)

1. **Parity / odd-index existence (Rédei-style).** The LRC target is σ-EVEN:
   the lonely measure is invariant under t ↦ −t, so its sign-isotypic component
   vanishes identically — there is no odd index to carry existence (THM-581/582;
   a lonely tournament has vertex 0 as a source, hence is not self-converse, so
   the reversal is not even an involution on its paths). The fleet's parity
   successes (Rédei, THM-643/644, THM-647 anti-Rédei, LEM-003 freeness) all live
   in the σ-ODD half.
2. **Sign-reversing involutions on the danger nerve.** The two geometric
   involutions available (negation, half-shift) are measure-preserving, hence
   sign-PRESERVING — they cancel nothing (claude-S606 twisted-involution audit).
   The Garsia–Milne pivot reduces p₀ to the irredundant Helly nerve — a genuine
   reformulation — but the un-cancellable core is exactly the open node.
3. **The OCF/Walsh factorization (THM-076) as a template for the resonance sum.**
   Refuted with a stated mechanism (opus-S172, klein-S202): THM-076 telescopes
   because tournament cycle-coverings are vertex-DISJOINT; the LRC covering is
   k/7 = 13/7 > 1 OVER-covering, TV(W′) ~ spread², and Σ|𝒲̂| diverges per
   covering. Deeper: THM-076's mechanism is telescoping — every term equal, no
   oscillation — the ANTITHESIS of √-cancellation. It can only ever produce
   absolute (non-cancelling) bounds; the mid-band's smallness is genuine phase
   oscillation.
4. **Gauss/Weil √p from the QR/Paley spectral results (THM-125/126/162).** The
   one proved √ fact — |λ_k(T_p)|² = (p+1)/4 via |g|² = p — is a complete
   multiplicative-character sum at prime modulus, used for extremality/floors
   (lower bounds). The LRC grid-resonance sum has no multiplicative character on
   its index set, needs an UPPER bound, and must work at composite V, where the
   QR field is exactly what is absent (THM-640: mod-14 cannot orient 42 of 78
   pairs — the "why 14 is hard" fact). Also mod-2 counting of good grid points
   is falsified at the boundary (tight AP @ V=13 has exactly ZERO good periods;
   MISTAKE-129: existence is a MAX, never a count).

## What DOES transfer (directional, and two of it actionable)

1. **The moment method, not the involution.** The tournament side's own
   √-cancellation (THM-438: Paley cluster integrals → Catalan, R(p) → e) was
   proved by Weil square-root cancellation + free-probability moment method —
   and its own sign-reversing involution is STILL OPEN there. The honest
   transport: attack the ~570-corner-phase sum with moment/equidistribution
   tools that survive over-covering, i.e. exploit cancellation, not summability.
2. **"Cancellation lives BETWEEN levels, not within" (kps-2 / THM-504).** The
   sinc-lattice kernel is positive on |t| ≤ 6 — half a period — so the dominant
   low-height resonance terms are SAME-SIGN; assuming quasi-random signs within
   a level is the recorded wrong instinct. The real cancellation is Abel-style
   across support levels. Any mid-band resonance-sum attack should organize by
   level and Abel-sum, not hope for within-level randomness.
3. **The divisibility/resonance law for grid periods.** Which j/V are good is
   decided by gcd(j, V) collapsing the teeth to sub-lattices (census-verified;
   the ONLY no-good-period clusters are the tight complete-residue APs at prime
   V = k). This is the multiplicative action on the Farey skeleton — real
   structure, already exploited by the dichotomy legs.
4. **The three-distance / Denjoy–Koksma organization of the snap points.** The
   mid-band question IS a three-distance count on a perturbed AP of snap
   points (klein-S208 axis (b)); the corpus's own proof that LRC(AP) IS the
   Steinhaus three-distance theorem says the regime structure is Farey — this
   is the terrain, not a shortcut.
5. **The winding-tournament dictionary is exact but is a relabeling.**
   μ_{1/7}(E) = P_x[T(x) has a scale-1/7 local sink] is an exact event
   description — but discrete tournament invariants (score, H, parity) provably
   forget the metric (tight and non-tight sets realize identical iso-class
   sets). Use it for intuition; do not expect a census reduction.

## The one-sentence synthesis

The tournament corpus's gift to the mid-band is not a tool but a compass: it
proves the node is σ-even/metric/Diophantine (so parity, counts, characters,
and summability are dead on arrival — each by a recorded theorem), and its own
hardest cancellation was won by the moment method organized across levels —
which is precisely the shape of the two live instruments the fleet already
holds (THM-665/666's aliasing program and the pair-sum/renormalization-tower
axes). The mining's value is negative space: four expensive dead ends
permanently marked, and the live directions confirmed independent of LRC
internals.

Related: THM-076/438/504/567/581/582/640/643/644/647/665/666/667/668, LEM-003,
MISTAKE-129, opus-S171/S172, klein-S202/S208, kps-2, claude-S606,
`the-two-indices-redei-is-odd-lonely-is-even…`, `the-lonely-runner-is-a-random-round-tournament`,
`cuts-as-farey-geodesics…`, `everything-is-the-triangle` (klein-S30 addendum).


---

## ADDENDUM (boxeph-S9, same day): the map's parity prediction became a theorem — and then the 2-adic descent revived parity anyway

Two corrections-by-extension from the evening's landings, recorded for honesty
and because together they reshape the lever map's conclusion:

1. **The depth-1 blindness is now PROVED, not just predicted.** mac-mini's
   LEM-020 built the Rédei involution τ ↔ 1−τ on witness sets and proved the
   parity law #witnesses ≡ #fixed (mod 2) — and that on covering sets the
   fixed layer is EMPTY, so witness counts are ALWAYS EVEN. That is this map's
   block (1) ("the target is σ-even, no odd index") upgraded from prediction
   to mechanism: the anatomy of why nine absolute bounds failed is one signed
   bit that cancels at depth 1.

2. **But "parity cannot transfer" was a DEPTH-1 statement, and the tower goes
   deeper.** The doubled-primes reflection's ×2 hinge was the hint the map
   underweighted: descending 2-adically, the FIRST LIVE parity layer is depth
   4 (LEM-021: witness at odd/16 with clearance 1/8 ⟺ no 16-multiple AND odd
   speeds miss a unit ±class mod 16 — covering forces 8∣v, not 16∣v, so the
   depth-4 layer has a nonempty fixed set). It decides 18.8% of covering
   unconditionally and joins the dispatch family (1/13, 8/17, 1/2, 1/8).
   The correct lever statement is: **parity is dead at the top of the 2-adic
   tower and alive below it** — the involution transfers, the INDEX transfers
   only after enough doubling.

## The three existence mechanisms (the synthesis the map should have ended with)

The fleet now holds three independent witness-existence instruments, each the
same abstract move (a count forced away from zero) in a different metric:

| mechanism | forcing | territory | state |
|---|---|---|---|
| **metric** (reverse triangle) | count ≥ mean·V − deviation·V > 0 | wide regime V > V₀ = √(TV/12∫W) | Lean complete (LRCGridPort, 10 thms S6–S8) |
| **parity** (2-adic Rédei ladder) | count ≡ #fixed (mod 2), fixed ≠ ∅ at depth ≥ 4 | dyadic dispatch slices (18.8%+ of covering) | LEM-020/021 proved; klein HYP-5850 Lean in flight |
| **signed characters** (klein's program) | diagonal-only suppression + t₂ hyperbola counting | the compressed/quarantined diagonal family | one lemma remains (t≥3 Hölder/Parseval chain) |

The Freiman/E3 Lean axis (opus finset_min_burden_isAP, kps SchurPeel/Rigidity)
quarantines the near-dilated diagonal these instruments hand off to. The
GrandAssembly's 8-clause ResidualObligation is where they compose: each new
mechanism deletes a clause. The map's one-sentence verdict stands amended:
*the tournament corpus's gift was not one compass but the whole toolkit —
its existence tricks transfer one by one, each at the right depth.*
