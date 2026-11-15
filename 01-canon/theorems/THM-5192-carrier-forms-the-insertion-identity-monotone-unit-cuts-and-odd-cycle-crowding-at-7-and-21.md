---
id: THM-5192
title: "Carrier forms, the insertion identity, monotone unit cuts, and odd-cycle crowding at 7 and 21. (A) For every marked block C with data (a,q,s,e,t,f,n), the strong tournament made of a marked guard, an X-tape, C, a Y-tape and an unmarked guard, with one apex, has exactly F_C(X,Y) = qXY + tX + sY + a(p(X)+p(Y)) + n + f(X mod 2) + e(Y mod 2) Hamiltonian paths, for all integers X, Y >= 1, where p(2u+1) = p(u), p(2u) = p(u) + u + (u mod 2), p(0) = p(1) = 0. (B) Insertion identity: adding a vertex that beats exactly X changes the count by H(T+x) - 2H(T) = E + R - Z, with Z the Hamiltonian paths that finish X^c before entering X, E the excess X->X^c steps along Hamiltonian paths and R the reverse-junction covers. (C) FINITE-EXACT: through order 8 every cut with a double X->X^c crossing has defect at least 3, so every unit cut is monotone; unit children of one sign have at most two types through order 9. (D) FINITE-EXACT: the carrier forms of the strong marked cores of orders 3-5 attain every odd N in [333, 10^9], those of order 6 every odd N in [615, 10^8]. (E) FINITE-EXACT with classical inputs: 7 and 21 are the values skipped by odd-cycle crowding at the orders to which Busch's bound confines them"
status: >
  PROVED (A), (B) (elementary). FINITE-EXACT (C), (D), (E).
  (A) generalizes the fixed two-tape carrier of chanpak_absolute_steering_20261115/TWO_TAPE_STEERING.md
  (credited: its (B,Q,T,S,h) = (2855,5760,4034,3914,45) form for one count-45 core) to every block, with exact
  constants. Checks:
  - the closed form against the cut transfer on 1,058,400 (core, X, Y) triples, all 216 strong marked cores of
    orders 3-5, X, Y <= 70;
  - the transfer against direct Hamiltonian-path DP on 48 random carrier tournaments of up to 18 vertices.
  (C):
  - all 1,747,520 cuts of the 6,880 order-8 classes (and all cuts of orders 4-7); a random sample of 1,000
    order-9 classes, minimum 15;
  - child types: all armed classes of orders 5-9 (307 at order 9).
  (D): exhaustive enumeration of carrier values (4.9*10^10 configurations at orders 3-5 to 10^9, 8.3*10^9 at order 6
  to 10^8).
  (E): the (c3, c5, d) census of all tournaments of orders 4-6, with Moon's pancyclicity and triangle bound and
  Busch's bound (KNOWN, cited).
  NOT CLAIMED: Chan-Pak (OPEN); the Finite Carrier Conjecture (HYP-9722); the general crossing lemma (HYP-9732).
session: mac-mini-2026-11-15-carriers
source: 05-knowledge/results/chanpak_waves_20261113/README.md#wave-27-carriers-crossings-and-crowding--a-finite-list-of-digit-quadratic-forms-covers-333-to-108
scripts:
  - 05-knowledge/results/chanpak_waves_20261113/scripts/closedform.c, verify_union.c, unioncar.c, greedy_cores.c
  - 05-knowledge/results/chanpak_waves_20261113/scripts/zer.c, zer_all.c, childtypes.c, sigpat.c, strongarm.c
  - 05-knowledge/results/chanpak_waves_20261113/scripts/cycles.py
related:
  - THM-5182 (pair expansion, passive channels), HYP-9712 (sign freezing), HYP-9722, HYP-9732,
    chanpak_absolute_steering_20261115 (two-tape steering, O(N^(3/8)) brackets),
    chanpak_viability_shape_20261115 (first moment of the fixed carrier), chanpak_reservoir_slide_20261115
    (unit slide, raw-unit components), busch_staircase_fermat_block_cuts_20261113 (Busch staircase)
---

# Carrier forms, the insertion identity, monotone unit cuts, and odd-cycle crowding

## (A) The carrier form

**Construction.** Let C be any marked block (a tournament with a mark set; the apex beats the marked vertices),
with cut-transfer data (a, q, s, e, t, f, n) as in THM-5182. For integers X, Y ≥ 1 let 𝒯_C(X, Y) be the
following ordinal sum, read left to right, plus one apex:
- a marked guard g₁;
- the X-tape: ⌊log₂X⌋ transitive vertices carrying the binary digits of X below its leading 1, most significant
  first, with mark = digit;
- the block C;
- the Y-tape: ⌊log₂Y⌋ transitive vertices, the j-th of which is unmarked iff bit j−1 of Y is 1;
- an unmarked guard g₂, playing the leading 1 of Y.

The apex beats g₁, and every spine vertex reaches g₂, which beats the apex. So 𝒯_C(X, Y) is strong for every
block C.

**Theorem (A).** H(𝒯_C(X, Y)) equals

    F_C(X, Y) = qXY + tX + sY + a(p(X) + p(Y)) + n + f·(X mod 2) + e·(Y mod 2),

where p(0) = p(1) = 0, p(2u + 1) = p(u) and p(2u) = p(u) + u + (u mod 2).

*Proof.* A singleton with mark μ acts on the cut state by

    (h, N, S, d) ↦ (h, N + (1 − μ)(S + d), 2S + μh, μh).

1. **Prefix.** After g₁ and the X-tape the state is (1, p(X) + 1, X, X mod 2). This holds at X = 1 (just g₁), and
   appending a digit b turns X into 2X + b exactly as p's recursion requires.
2. **Core.** The block C turns this into (a, a·p(X) + n + tX + f(X mod 2), qX + s, e).
3. **Suffix.** The Y-tape and g₂ send any state (h, N, S, d) to N + SY + h·p(Y) + d(Y mod 2). This holds at Y = 1
   (g₂ alone adds S + d). Prepending a vertex for the lowest bit of Y gives the two cases p(2Y′ + 1) = p(Y′) and
   p(2Y′) = p(Y′) + Y′ + (Y′ mod 2).

Combine the three steps. ∎

The shape is a hyperbola qXY + tX + sY with a binary surcharge a(p(X) + p(Y)) and two parity bits. The function p
is the 2-regular inversion surcharge: p(Y) sums ⌊Y/2^(j+1)⌋ plus one over the 0-digits j of Y lying below a 1. The
3-cycle with one marked vertex has data (3, 12, 3, 1, 6, 2, 5), so

    F(X, Y) = 12XY + 6X + 3Y + 3(p(X) + p(Y)) + 5 + 2(X mod 2) + (Y mod 2).

## (B) The insertion identity

Let x beat exactly X ⊆ V(T), with X ≠ ∅, V. For a Hamiltonian path P of T, x can be inserted:
- at the front iff first(P) ∈ X;
- at the end iff last(P) ∉ X;
- between consecutive v_i, v_(i+1) iff v_i ∉ X and v_(i+1) ∈ X.

For a binary word, [first ∈ X] + [last ∉ X] + #(X^c X) = 1 + #(X X^c). So P has exactly 1 + #XX^c(P) insertion
slots. Every other Hamiltonian path of T + x has x between some u ∉ X and some w ∈ X with w → u; call their number
R. Hence

    H(T + x) − 2H(T) = E + R − Z,

where:
- Z = #{P : no X → X^c step}, the paths that cover X^c first and then X;
- E = Σ_P max(0, #XX^c(P) − 1).

So a unit cut is one where the reverse-junction covers balance the paths of shape [X^c][X] up to ±1.

## (C) Monotone unit cuts and essentially unique unit children (FINITE-EXACT)

**Crossing minimum.** Every cut X of every tournament of order 4–8 such that some Hamiltonian path crosses from X to
X^c twice (E ≥ 1) has H(T + x) ≥ 2H(T) + 3. At order 8 that is 1,620,135 such cuts out of 1,747,520. On 1,000
random order-9 classes the minimum is 15.

Consequently every cut with defect ≤ 1, in particular every unit cut, is **monotone**: each Hamiltonian path
visits X in at most two blocks, the second at the very end (shape [X^c][X][X^c][X], parts possibly empty). For unit
cuts R = Z + ε exactly. All 2,783 unit cuts of order 8 have E = 0.

**Essential uniqueness.** Classify the unit children T + x of a given sign by a converse-symmetrized invariant (H,
sorted Π row and column sums, scores). Every armed seed of orders 5–9 has at most two child types per sign.

| order | armed seeds | +children: 1 type / 2 types | −children: 1 type / 2 types |
|---:|---:|---|---|
| 8 | 133 | 113 / 20 | 120 / 13 |
| 9 | 307 | 281 / 26 | 289 / 18 |

**Joint signatures.** Tabulate, per armed seed, the set of signatures of its +children and of its −children:
- the same-sign children of a seed almost always share one signature;
- the pair (armed, armed) never occurs (HYP-9712);
- (+only, +only) never occurs, and (−only, −only) occurs twice at order 8.

**Armed does not imply strong** once arbitrary cuts are allowed: 81 of the 133 armed order-8 classes are
reducible. So the leaf-model raw-unit law (chanpak_reservoir_slide_20261115/RAW_UNIT_COMPONENTS.md) does not
extend to arbitrary out-sets.

## (D) A finite list of carrier forms covers 333 to 10^9 (FINITE-EXACT)

Take the union, over strong marked cores C, of all values F_C(X, Y) with X, Y ≥ 1
([`unioncar.c`](../../05-knowledge/results/chanpak_waves_20261113/scripts/unioncar.c)):

| cores | marked cores (distinct forms) | Σ 2/q over distinct forms | odd misses above 21 |
|---|---:|---:|---|
| order 3 | 8 (4) | 0.667 | about 11% of all odd N (to 2·10^7) |
| orders 3–4 | 24 (20) | 1.667 | 10,853 to 10^8, density about N^(−0.58), unbounded |
| orders 3–5 | 216 (168) | 4.780 | **32 to 10^9, all ≤ 331** |
| order 6 alone | 2,240 (1,950) | 12.597 | 219 to 10^8, all ≤ 613 |

A greedy cover of [333, 10^7] uses 48 of the 216 cores, and the same 48 forms cover [333, 10^9]:
- the 3-cycle carrier above covers 84% alone;
- ten masks of the count-5 core of order 4 bring it to 99.9%;
- 36 order-5 cores take the sporadic rest.

**Multiplicity.** For orders 3–5 the number of representations averages:

| N | mean / ln N | window minimum |
|---|---:|---:|
| 10^4 | 5.03 | 4 |
| 10^5 | 4.43 | 4 |
| 10^6 | 4.48 | 8 |
| 10^7 | 4.70 | 14 |
| 10^8 | 4.89 | 18 |

This matches the first-moment prediction Σ_C (2/q_C)·ln N = 4.78·ln N.

The only local bias is mod 3: near 10^7 the classes N ≡ 0, 1, 2 average 87.9, 63.0 and 76.2 representations.
Primes average 69.1, the mean of their two classes.

## (E) Odd-cycle crowding at 7 and 21 (FINITE-EXACT, classical inputs)

Three ingredients confine and then skip the two values.
1. **Multiplicativity.** H is multiplicative over strong components. 7 is prime and 21 = 3·7, so both would need a
   strong witness.
2. **Busch's bound.** A strong tournament of order n has H ≥ b(n) = 3, 5, 9, 15, 25, … (Busch 2006). So a strong
   witness has order ≤ 4 for 7 and ≤ 6 for 21.
3. **Odd-cycle crowding.** By the odd-cycle expansion, H = 1 + 2c + 4d + 8t + … with c odd cycles, d disjoint
   pairs and t disjoint triples. Moon: a strong tournament is pancyclic and has at least n − 2 triangles.

The census of orders 4–6 ([`cycles.py`](../../05-knowledge/results/chanpak_waves_20261113/scripts/cycles.py)) gives:
- **Order 4:** c ≤ 2, so H ∈ {1, 3, 5}. The value 7 = 1 + 2·3 would need exactly three pairwise-meeting odd
  cycles. At order 5, strong tournaments already have c ≥ 4 (at least 3 triangles and a 5-cycle), so c = 3 never
  occurs.
- **Order 6:** the strong classes realize (c, d) ∈ {(6,1), (7,0), (8,0), (9,0), (9,1), (10,2), (11,0), (11,1),
  (12,0), (12,1), (12,2), (13,1), (14,1), (14,2), (14,4), (16,0), (16,1), (16,2), (19,1), (20,1)}. So c + 2d takes
  the values 7, 8, 9, 11, 12, 13, 14, 15, 16, 18, 20, 21, 22. The value c + 2d = 10, which is exactly H = 21, is
  missed. Its neighbours 9 (H = 19) and 11 (H = 23) are hit. The non-strong order-6 counts are products of smaller
  ones and never 21.

So 7 and 21 are the values that odd-cycle crowding skips at the only orders allowed to them. From order 7 on, the
Busch staircase has long runs, and from 333 on the carrier forms of (D) cover every odd value tested.
