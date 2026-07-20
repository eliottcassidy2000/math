# The Involution Must Live Outside: Moon = the tournament fixed-point property, and the one dihedral law

**boxeph-2026-07-20-S150** · companions: THM-1365 (Galois reduction of JC),
`dihedral_quaternion_tournament_census_boxeph_S150.py` (HYP-8205)

The owner asked for the tournament↔dihedral relationship "thoroughly — there
will be multiple instances." The instances are indeed everywhere in canon —
Paley D_{2p} (THM-127), the staircase grid reflection = converse (THM-280,
THM-643), negation on circulants, the involutory anti-automorphism theorem
(THM-024) and anti-Rédei (THM-647), the complement-antipodal map (THM-584).
This session found the single law behind all of them, and it is the same law
that runs the new Jacobian salvage theorem THM-1365.

## 1. The law

**In a tournament, an involution can never be a symmetry — only a
reflection.** Moon: |Aut(T)| is odd (an involutive automorphism would have to
fix the arc between some swapped pair both ways). So whenever a Z₂ acts, it is
forced OUTSIDE Aut, into the anti-automorphism coset — and then the extension

  1 → Aut(T) → Aut±(T) → Z₂ → 1   (Aut± = autos ∪ anti-autos, T self-converse)

**always splits** (Cauchy: an order-2 element exists in Aut±, and it cannot
lie in odd-order Aut — THM-024's proof). There is no quaternionic option:

**Lemma (anti-order law, new small).** Every anti-automorphism σ has order
≡ 2 mod 4. (σ² ∈ Aut has odd order m, so |σ| = 2m.) The order-4 twist that
builds Q₈ out of D₄'s material is arithmetically impossible for tournaments.

**Census (this session, exact, full n ≤ 6):** classes 2/4/12/56 (= A000568),
SC classes 2/2/8/12, quaternionic 0/0/0/0; anti-order spectra all ⊆ {2, 6},
first order-6 anti-autos at n = 6 (classes with |Aut| = 3, 9 — an order-3
automorphism commuting with the reflection); #anti-autos = |Aut| in every SC
class (coset size). QR₅ shows Aut± = D₅ literally: 5 rotations + 5 involutive
reflections. Circulant verification n ≤ 13: negation x ↦ −x is an involutive
anti-automorphism of every circulant tournament (all 2^{(n−1)/2} of them) —
rotations are automorphisms, reflections are anti-automorphisms: **the
dihedral group is the unique extension pattern tournament theory permits.**

## 2. The same law in Keller-land (THM-1365)

Replace (tournament, Aut, converse) by (Keller map, deck group, Galois
closure). THM-1365(A): a polynomial deck group acts FREELY on C^n — a deck
element fixing a point has identity differential there (dF invertible) and
Cartan's lemma kills it. Freeness is the Keller analogue of Moon: **the
symmetry cannot sit inside** (no fixed points ↔ no fixed arc-pairs). In dim 2,
linearizability (Kambayashi) says a free finite group cannot exist at all, so
polynomial-Galois Keller maps are automorphisms — the extension either splits
into an automorphism or the deck group was never there. The kernel
counterexample K evades with deck group N_{S₃}(S₂)/S₂ = 1: **deck-poor**,
exactly as non-SC tournaments (the sea, 96% at n = 8) are reflection-poor.

The dictionary:

| tournament world | Keller world |
|---|---|
| Aut(T) odd (Moon) | deck acts freely (Cartan, THM-1365 A) |
| involution forced outside Aut | Z₂ forced outside Aut(C²) (linearization) |
| Aut± splits, always dihedral (THM-024) | polynomial-Galois ⟹ automorphism (THM-1365 B) |
| anti-order ≡ 2 mod 4 (this session) | no free finite groups = FP(n) bridge |
| SC classes = the spine (reflection-rich) | Galois maps = automorphism locus |
| non-SC sea = generic, reflection-free | non-Galois counterexamples, deck-trivial |
| anti-Rédei: reflected Ham-paths odd (THM-647) | odd-degree conjecture (S146 dictionary) |

The last row is the striking one: THM-647 counts σ-anti-palindromic
Hamiltonian paths and finds ODD; the odd-degree conjecture says Keller
counterexample cover-degrees are odd. Both are Rédei parity surviving a
reflection quotient. If the h-monoid ⟷ Keller-degree-monoid dictionary (S146)
extends to the reflection structure, the anti-Rédei theorem is the tournament
shadow of a provable "odd-degree theorem for Galois-stable degree data" —
worth a dedicated session.

## 3. Where the dihedral thread should go next

1. **Anti-order spectrum law:** conjecture — the spectrum of anti-auto orders
   of SC classes at size n is exactly {2·d : d odd, d | |Aut| realizable};
   first 10 at n = 11? (needs an order-5 auto commuting with a reflection —
   QR₁₁ has Z₁₁ rotations, reflections… 2·11 = 22 at n = 11 via σ·rotation.)
   Cheap census extension.
2. **The Aut± Burnside count:** count SC classes by Aut±-type (dihedral D_m
   vs Z₂ × odd) — refine the SC census the way THM-586 did Paley arithmetic.
3. **Anti-Rédei ⟷ odd-degree:** make row 7 of the dictionary precise (the
   reflection-equivariant h-monoid).
