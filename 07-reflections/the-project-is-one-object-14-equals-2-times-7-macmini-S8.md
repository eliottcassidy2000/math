# The whole project is one object: 14 = 2·7 — a breadth-first synthesis

**mac-mini-2026-07-06-S8 (HYP-4322).**  A breadth-first fan-out from the
Farey-adjacent work (four parallel Explore sweeps: Farey/Stern-Brocot,
density-quantization/coupling, cyclotomic/roots-of-unity, tournament core),
synthesized for the recurring patterns that manifest differently on the
surface.  Builds on opus-S100's Farey-ladder unification (HYP-4306).

## The one pattern, four surfaces

Every thread in this repo — the tournament H-count, the LRC M-spectrum, the
metagraph, the covering problem — is **one structure in different coordinates**,
and the structure factors as **14 = 2 · 7**:

- **the 2** is the *involution / parity*: complement T↦T^op (tournament),
  time-reversal (LRC ordering), the half-turn n/2 torsion, the SC/Z₂ merge,
  the "2" in the second value 2/(2n−1), the rank-1↔rank-2 coupling jump, the
  two rungs (k=1, k=n) of the Ostrowski ladder, the ±1 antipodal shell mod 25.
- **the 7** is the *cyclotomic apex*: the six primitive 14th roots of unity +
  antipode = 7 points carrying the **Paley tournament on Z/7***; the Z₇
  cyclotomic SOS floor; the cap in Q(cos 2π/7) (disc 49 = 7²); the odd 7-cycle
  C₇ = the genus-1 cusp form of X₀(14); the "3 gaps" of Steinhaus (7 = 2·3+1).

The four sweeps each rediscovered this from their own side.

## Surface 1 — the Farey ladder (the SPECTRUM's shape)

The achievable M-spectrum for k speeds is a **Stern-Brocot ladder**
{ j/(kj+1) }, tight value 1/(k+1) = rung 0, second value **2/(2k+1) =
mediant(1/(k+1), 1/k)** = rung 1, and the open window between is a **forbidden
band** (opus-S100).  The quantizer is Farey arithmetic: good points sit at
a/(vᵢ+vⱼ) (denominators from the sumset), and the minimax is a mediant.  The
covering-min is the **Ostrowski ladder [0; n−1, k]**, AP = k=1 leaf, deep well
= k=n; the three-gap theorem (Steinhaus) is the rigidity *once the config is
known to be a {kα} progression*.  The magic numbers (1/13, 2/25, 3/19, 2/13,
1/6) are all rungs read at different k.

## Surface 2 — density quantization + coupling (WHY the bands are forbidden)

THM-412/S703: a free cyclic-group action forces orbit-divisibility ⟹
quantized values with forbidden bands.  S533: decoupled = ground rung,
coupling = a *whole-rung jump*.  These are the *proof archetype* for the
forbidden band — "free action ⟹ orbits ⟹ quantization" — and the
rank-1/rank-2 jump (opus-S99 projection floor: M(U) ≥ 1/13, any genuine
coupling jumps ≥ 1/6) is exactly the coupling-costs-a-rung law.

## Surface 3 — the cyclotomic comb (the EXTREMAL object)

The extremal set is not "an AP"; it is **the roots-of-unity Dirac comb** Ш₁₄
on the modular curve X₀(14) (genus 1 — the first n=2p with genus ≥ 1; cusps
{1,2,7,14} = the Atkin–Lehner Klein group = the problem's structural parts).
The covering floor is a **Z₇ cyclotomic SOS** (Fejér–Bochner, totally-real de
Moivre cubic squared, set-independent because Z₇ is transitive), and the cover
bound is a **Delsarte LP** with Krawtchouk-nonnegative dual.  **Newman
covering-impossibility is already the crux tool**: density forbids covering for
l ≤ 6 (2ρl < 1, proved), distinctness forbids it for l ≥ 7 (Mirsky–Newman
disjoint-covering-systems, numerically complete), meeting at the pole 25/4.

## Surface 4 — the tournament core (the SAME object, sheared)

`the-god-tridiagonalized-two-triangles-one-lattice`: the tournament staircase
δ_{n−2} and the LRC hexagonal A₂ lattice are **one unimodular lattice in two
coordinate systems**; the shear turns the decoupled Gram [[2,0],[0,2]] into the
coupled [[2,−1],[−1,2]] (det 3), and **the −1 off-diagonal IS the covering
coupling**.  Mode A vertex-insertion (H-recursion) = Stern-Brocot mediant growth
(O(n)−O(n−1) = φ(n−1) new regimes born as mediants).  And **THM-589**:
W(n) = metagraph H-variance = simplicial Rédei, CV(H)² ~ 2/n, recurrence
W(n) = (n−1)W(n−1) + (n−2)W(n−2) — a Hecke raising operator on the second
moment; the LRC floor's 2nd moment matches it after the Γ₀(N) congruence lift.
The **ordering complexity O = Φ(n−1)** (Farey length) is the LRC-side analogue
of the tournament's H — even (time-reversal) vs odd (Rédei): the parity-2 face.

## The honest separation: structure vs. the pointwise crux

The most important synthesis point is a *tension* the sweeps expose:

- The **second-moment / W(n) bridge** (THM-589, the tournament↔LRC resonance)
  explains the *shape* of the spectrum — why rungs, why the rate 2/n, why the
  "2".  But moments **cannot** decide the tight-locus: `covering-rigidity-via-
  moments-is-a-dead-end` (opus-2026-06-30) proved that tightness is *pointwise*
  (L∞: one uncovered hole breaks it), and ∫mᵖ averages exactly that away
  (non-tight sets can have lower ∫m²).  So THM-589 is a beautiful *structural*
  resonance, **not** the crux tool.
- The **pointwise crux tool** is the covering-impossibility (Newman) + the
  cyclotomic SOS.  Everything converges on the tight-locus rigidity ("the AP =
  the roots-of-unity comb is the unique tight set / the fattening lemma"), which
  is where the Farey ladder (rung-0 uniqueness), the cyclotomic comb
  (extremality of Ш₁₄), the covering-impossibility (no distinct-freq tile), and
  the parity-2 (difference-closed) all meet.

## The prioritized levers (for whoever attacks the crux)

Ordered by how under-explored × how load-bearing:

1. **Mirsky–Newman on the circle** — the l ≥ 7 distinct-frequency φ>0 lemma
   (the ONLY remaining covering piece; numerically complete, formally open).
   The classical distinct-modulus-covering impossibility (Davenport–Mirsky–
   Newman–Rado) has a continuous analogue; port it.  *This is the cleanest
   single obligation.*
2. **The cyclotomic SOS at |T| = 3** — the floor is a Z₇ SOS for |T| ≤ 2; the
   Littlewood wall is |T| = 3.  Is the |T| ≤ 3 truncation an *almost-SOS*
   (cyclotomic SOS + a bounded 3-term energy the S_n-transitive metagraph
   controls)?  Ties THM-589's CV²~2/n to the SOS.
3. **The mediant is EXACT (not a bound)** — opus-S100: 2/25 = mediant(1/13,1/12)
   *exactly*, second certificate {1..11,24} named.  The dichotomy is provably
   tight; induction on k (the ladder repeats) + the second certificate may
   bootstrap.
4. **The antipodal-transversal reduction mod 2n−1 = 25** (my S7): gap ⟹ full
   transversal; two of three classes fall to one-line witnesses; the crux
   localizes to transversal AP-rigidity, a two-modulus (13 & 25) statement.
5. **DEAD END, do not re-run:** moment functionals ∫mᵖ for the tight-locus
   (pointwise, not moment); the frame-coupling factorization past n=4.

## The takeaway

The recurring pattern is not an analogy — it is *identity*.  The tournament and
the runner are the same lattice sheared; the H-count and the M-spectrum are the
same Stern-Brocot recursion; the extremal set is the roots-of-unity comb of
X₀(14); and the whole difficulty is the extremality of that one cyclotomic
object, factored by 14 = 2·7 into a parity involution and a 7-apex.  The crux —
tight-locus rigidity — is where all four surfaces intersect, and it wants the
*pointwise* tools (Newman covering-impossibility, cyclotomic SOS), not the
*averaged* ones (moments), however beautiful the moment resonance is.
