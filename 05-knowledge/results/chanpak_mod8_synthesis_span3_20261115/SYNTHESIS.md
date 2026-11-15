# A mod-8 dictionary: Collatz and 3x−1 edges, multiplication tables mod odd numbers, Paley tournaments, and the norm form of disjoint path pairs

Session klein-2026-11-15-synthesis. **No open problem is solved here.** The aim is to say precisely where the
tournament work meets other open problems:
- §1–§5 are proved (elementary or classical) and verified by `verify_mod8.py`;
- §6 marks analogies as analogies.

## 1. Three multiplication tables and their signs (PROVED, classical)

For odd m and a prime to m, the permutation x ↦ ax of ℤ/m has sign equal to the Jacobi symbol (a/m)
(Zolotarev–Frobenius). Three multipliers come from Collatz edges:
- **doubling ×2**, the inverse of the halving edge n → n/2;
- **negation ×(−1)**, the mirror. The identity −(3(−x) + 1) = 3x − 1 conjugates the 3x+1 map to the 3x−1 map, the
  Collatz map on the negative integers;
- **their product ×(−2)**, one "back and forth" between a Collatz edge and its mirror.

The affine steps T₊(x) = (3x+1)/2 and T₋(x) = (3x−1)/2 on ℤ/m satisfy T₋(x) = −T₊(−x), and both have sign (6/m).
Verified for every odd m < 400.

## 2. The three quadratic characters mod 8 (PROVED)

(ℤ/8)* = {1, 3, 5, 7} = ⟨−1⟩ × ⟨3⟩: **the Collatz multiplier and the mirror generate the 2-adic units mod 8**, and
every odd residue mod 2^k is ±3^j. The three nontrivial characters are the three ways to split this Klein group.

| character | kernel | primes p with value +1 | sign of | tournament side | Collatz side |
|---|---|---|---|---|---|
| χ₋₁ = (−1/·) | {1, 5} | p ≡ 1 mod 4 | x ↦ −x | Paley tournaments exist exactly when χ₋₁(p) = −1 | the mirror |
| χ₂ = (2/·) | {1, 7} | p ≡ ±1 mod 8 | x ↦ 2x | ×2 is an automorphism of P_p for p ≡ 7 (8), and reverses it for p ≡ 3 (8) | halving edges |
| χ₋₂ = (−2/·) | **{1, 3} = ⟨3⟩** | p ≡ 1, 3 mod 8, i.e. **p = x² + 2y²** | x ↦ −2x | diagonal joins N = (H + P₂)² + 2H² ≡ 3 (8); ×(−2) is an automorphism of P_p iff p ≡ 3 (8) | mirror ∘ doubling |

Three facts follow.
- **The norm form of disjoint path pairs is the Collatz-multiplier character.**
  - A prime is x² + 2y² exactly when it is congruent mod 8 to a power of 3.
  - Diagonal joins (THM-5222) have H odd and P₂ even, so they land in the class 3 itself.
  - Here P₂ is the number of ordered pairs of disjoint nonempty paths, and β = (H + P₂) + H√−2.
- **For Paley primes (p ≡ 3 mod 4), three conditions are equivalent:**
  - doubling reverses P_p;
  - mirror ∘ doubling preserves P_p;
  - p = x² + 2y², the diagonal-join class p ≡ 3 (8).
- **In Collatz, mod 8 enters through 2-adic valuations.** v₂(3n+1) is 2, 1, ≥ 3, 1 for n ≡ 1, 3, 5, 7 (mod 8),
  while v₂(3n−1) is 1, ≥ 3, 1, 2.
  - The mirror n ↦ −n swaps 1 ↔ 7 and 3 ↔ 5 (mod 8). So it keeps χ₂(n) and flips χ₋₁(n) and χ₋₂(n).
  - The deep class (v₂ ≥ 3) is n ≡ −3 for 3x+1 and n ≡ 3 for 3x−1, that is ∓3⁻¹, since 3² ≡ 1 (mod 8).
  - The merges 3(4n+1)+1 = 4(3n+1) and 3(4n−1)−1 = 4(3n−1) are the "4 ± 1" laws.

So "mod 8 in Collatz" (the residue of the *argument*) and "mod 8 in the norm form" (the residue of the *prime*) are
two readings of one group (ℤ/8)* = ⟨−1, 3⟩. The diagonal-join character χ₋₂ is the unique nontrivial character
that fixes the Collatz multiplier.

## 3. Paley tournaments are the multiplication-table tournaments (PROVED and computed)

For p ≡ 3 mod 4, P_p has i → j iff j − i is a nonzero square. Multiplication by a is an automorphism iff (a/p) = 1,
and an isomorphism onto the converse iff (a/p) = −1 (verified p < 80 for a ∈ {2, −1, −2, 3, −3}).
- *The mirror is a converse map.* Negation always reverses P_p. So the Collatz mirror acts on Paley tournaments as
  the converse, whose path count is the same: H(T) = H(T^op).
- *Divisibility.* The affine maps x ↦ ax + b with (a/p) = 1 form a group of order p(p−1)/2. It acts freely on
  Hamiltonian paths, so p(p−1)/2 divides H(P_p).

Computed with `hp2.c`:

| p | H(P_p) | orbits H/(p(p−1)/2) | P₂(P_p) | diagonal join mod 8 | prime? |
|---:|---:|---:|---:|---:|---|
| 3 | 3 | 1 | 6 | 3 | no |
| 7 | 189 | 9 | 1,638 | 3 | no |
| 11 | 95,095 | 1,729 | 1,586,310 | 3 | no |
| 19 | 1,172,695,746,915 | 6,857,869,865 | 38,153,145,325,662 | 3 | no |

The orbit count 1,729 at p = 11 is a coincidence with the Hardy–Ramanujan number.

## 4. Artin's conjecture: when the doubling table is one cycle

The doubling table is a single cycle on (ℤ/p)* exactly when 2 is a primitive root mod p. Artin's conjecture, proved
under GRH by Hooley, predicts this for a proportion A = 0.3739… of primes. It is possible only when χ₂(p) = −1,
that is p ≡ 3, 5 (mod 8).

Measured for p < 2·10⁵: 2 is a primitive root for 0.7406 of the primes p ≡ 3 (mod 8) and 0.7488 of the primes
p ≡ 5 (mod 8), close to 2A = 0.7479.

For a Paley prime p ≡ 3 (mod 8) with 2 primitive, doubling is a converse map whose powers alternately preserve and
reverse P_p and run through every multiplier.

## 5. The metagraph and its norm equations (PROVED laws, assembled)

Take tournaments as vertices and host operations as edges. Each edge acts on the path-cover polynomial P_T(q), or on
β = (H + P₂) + H√−2 ∈ ℤ[√−2]:

| operation | law |
|---|---|
| ordinal sum A ⇒ B | P_(A⇒B) = P_A · P_B, so (H, F₂) multiplies coordinatewise |
| cyclic join A → B → c → A | H = Re(β_A β̄_B) (THM-5222) |
| diagonal join (B a copy of A) | H = N(β_A) = (H + P₂)² + 2H², a norm equation in ordered pairs of disjoint paths |
| double cone a → B → c → a | P_D(q) = q[P(q+1) + qP(q) − P(q−1)] ([LLT.md](../chanpak_generic_lines_llt_20261115/LLT.md)) |
| padding by TT_k | (H, F₂) ↦ (H, 2^k F₂) |

So the arithmetic of the counts lives in ℤ[√−2]. Its prime splitting is χ₋₂, the character §2 ties to the Collatz
multiplier.
- *The diagonal sub-family has a genuine reciprocity obstruction.* It only reaches N ≡ 3 (mod 8). This is the kind
  of obstruction Haag–Kertzer–Rickards–Stange found for Apollonian packings.
- *The full family has none.* Every residue class is equally dense ([LLT.md](../chanpak_generic_lines_llt_20261115/LLT.md)
  §1), so Chan–Pak is consistent with a local–global principle.

## 6. Where the tournament work meets other open problems

| open problem | meeting point here | kind |
|---|---|---|
| Chan–Pak (all odd except 7, 21) | single-kernel conjecture HYP-9772, plus a 1D local limit law | reduction (proved) |
| Collatz (3x+1) | C2 births 3H + 2σ = T∘D (clone_clocks_negative_collatz_20261113); (ℤ/8)* = ⟨−1, 3⟩ | proved identities; the rest is analogy |
| negative-integer Collatz (3x−1) | the mirror conjugacy; the mirror acts as the converse on Paley tournaments | proved |
| Artin's primitive root conjecture (base 2) | the doubling table is one cycle; it pairs with Paley primes p ≡ 3 (8) | analogy, with measured densities |
| quadratic reciprocity | Zolotarev signs of the three Collatz multipliers | proved, classical |
| primes of the form x² + 2y² (Fermat–Euler) | diagonal joins; the kernel of χ₋₂ is generated by 3 | proved |
| Apollonian local–global (Bourgain–Kontorovich; Haag–Kertzer–Rickards–Stange) | thin-family coverage; reciprocity obstructions for sub-families | analogy |
| Zaremba (Bourgain–Kontorovich density one) | the bilinear circle-method route to a local limit law | method transfer (proposed) |
| Lehmer's problem / Mahler measure | Riesz products of HYP-9782 (M(1 + z) = 1); growth constants ρ*, tribonacci 1.839, 5^(1/3) | analogy |

**The one-sentence synthesis.** The same Klein group (ℤ/8)* = ⟨−1, 3⟩ organizes all of these:
- the Collatz multiplier and its mirror generate it;
- its characters are the signs of the doubling and mirror multiplication tables;
- one of those characters, the one fixing 3, is exactly the splitting law of the norm form (H + P₂)² + 2H² that
  tournament path counts obey under diagonal joins.

The Collatz side asks about orbits of ⟨3, 2⁻¹⟩ on integers. The tournament side asks which integers a semigroup of
host operations reaches. Both are coverage questions for thin semigroups, with 2-adic and ℚ(√−2) local structure.
That is the shared shape; neither problem is reduced to the other.
