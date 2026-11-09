---
id: THM-4930
title: "Eternal colours of odd Hamiltonian-path counts. An odd h is the path count of some self-converse tournament (eternally blue) iff h = d^2 e with d a path count (or 1) and e = 1 or the path count of a strong self-converse tournament. Otherwise (eternally black) every tournament with H = h, at every order, has a distinct converse, so its classes come in converse pairs at every size. By Busch's minimum this is decidable value by value. Exactly 29 odd values below 225 are never self-converse: 7, 21 and 19, 23, 35, 39, 43, 47, 53, 55, 59, 63, 101, 103, 107, 109, 113, 119, 131, 163, 167, 181, 183, 185, 203, 205, 207, 209, 215. Of the 44 values below 225 that are black at birth (THM-4922), 17 turn blue later, 121 = 11^2 last, at order 10 as B11 => B11^op. Through 625 the eternally black values are exactly 7, 21 and 59 attainable values (27 below 225; 32 in [225, 625), all born at order 8), except 377, which needs order 12. Only two eternally blue values below 625 have no strong self-converse realisation: 121 = 11^2 (first self-converse at order 10) and 275 = 5^2 * 11 (first at order 13, five orders after birth)."
status: >
  PROVED: statements 1-2 (the criterion; Busch decidability; parity).
  FINITE-EXACT, on complete self-converse inventories:
  - statement 3, from orders <= 9 (THM-4922's census generator), order 10 (all 2^25 anti-invariant labelled
    tournaments) and order 11 with H < 625 (all anti-invariant extensions);
  - statement 4: exact for every odd h < 375 and for every h = 3 (mod 4) below 625; h = 377 is OPEN (it needs order 12).
  FINITE: witnesses for values above the exact range.
  KNOWN, credited:
  - the palindromic condensation (THM-4866(E));
  - the self-converse generating function D_h(z) = sum_(d^2 e = h) F_d(z^2)[1_(e=1) + Q_e(z)] (Session24 PROOF eq. (2);
    Session224 section 5);
  - the zero self-converse series for h in {7, 19, 21, 23, 35, 39, 43, 47, 53, 55, 59, 63} below 73 (Session24
    low_H_formulas.json);
  - Busch's minimum f (THM-4744, CITED);
  - even-order self-converse tournaments have H = 1 (mod 4) (THM-791);
  - insertion monotonicity (THM-1370);
  - every odd value to 2^32 except 7, 21 is a path count (THM-4875, CITED);
  - birth colours (THM-4922).
  New here:
  - the eternal colour as the stable form of the owner's two kinds;
  - the criterion stated as a membership test with Busch decidability;
  - the exact lists above 73;
  - the comparison of birth colours with eternal colours;
  - the order-10 and order-11 self-converse enumerations.
  Self-audited only (PRE-AUDIT).
session: mac-mini-2026-11-08-firstrealizers (owner directive 2026-11-08: first creators of odd values, one class or a converse pair, two kinds, the splitting pattern; prior work on 33)
source: 05-knowledge/results/eternal_colours_self_converse_values_20261108.md
depends_on: [THM-4866, THM-791, THM-1370, THM-4744, THM-4875, THM-4922]
verification:
  - 04-computation/experiments/eternal_colours_20261108/sc_holes.py (+ .out)
  - 04-computation/experiments/eternal_colours_20261108/sc_even10.c (+ sc_even10.out; -DN=6, 8 validation builds)
  - 04-computation/experiments/eternal_colours_20261108/sc_odd11.c, sc_odd11_fast.c (+ sc_odd11.out; -DCHECK self-test), check_order11.py (+ .out)
  - 04-computation/experiments/eternal_colours_20261108/sc_witness.c (+ .out)
---

# THM-4930 — eternal colours: which odd path counts are ever self-converse

## Setting

- `H(T)` is the number of Hamiltonian paths of `T`, and `T^op` is the converse. Since `H(T^op) = H(T)`, the
  tournaments with a given `H` at a given order are a union of converse orbits, of size 1 (self-converse) or 2.
- `Σ` is the set of path counts. It contains every odd number up to `2^32` except 7 and 21 (THM-4875). `Σ_sc` is the
  set of path counts of self-converse tournaments, and `Σ_sc^str` that of strong self-converse tournaments.
- **Busch's minimum** (THM-4744, CITED): a strong `m`-tournament has `H ≥ f(m)`, with `f(3q) = 3·5^(q−1)`,
  `f(3q+1) = 5^q`, `f(3q+2) = 9·5^(q−1)`; so `f = 3, 5, 9, 15, 25, 45, 75, 125, 225, 375, 625` at `m = 3, …, 13`.
  Put `B(e) = max{m : f(m) ≤ e}`.
- **Colours.** THM-4922's *birth colour* of `h` is blue if some first creator (a class of least order `n(h)` with
  `H = h`) is self-converse, black otherwise. The *eternal colour* is blue if `h ∈ Σ_sc`, black otherwise.

## Statement

1. **Criterion.** `Σ_sc = {d²·e : d ∈ Σ, e ∈ Σ_sc^str ∪ {1}}`.
   - If `h ∉ Σ_sc` (eternally black), then at every order every tournament with `H = h` lies in a converse orbit of
     size 2.
   - Every self-converse realisation is generated from smaller structures as a palindrome
     `B_1 ⇒ … ⇒ B_j ⇒ Z ⇒ B_j^op ⇒ … ⇒ B_1^op`, with `Z` strong self-converse or a single vertex or empty.
2. **Decidability and parity.**
   - `e ∈ Σ_sc^str` iff a strong self-converse tournament of order at most `B(e)` has `H = e`. So the eternal colour of
     every `h` is decided by the strong self-converse inventories through order `B(h)`.
   - Every self-converse tournament of even order has `H ≡ 1 (mod 4)` (THM-791).
   - So a squarefree `h ≡ 3 (mod 4)` is eternally blue iff a strong self-converse tournament of *odd* order at most
     `B(h)` has `H = h`.
3. **Below 225** (FINITE-EXACT + PROVED). The odd `h < 225` with `h ∉ Σ_sc` are exactly

       7, 21; 19, 23, 35, 39, 43, 47, 53, 55, 59, 63, 101, 103, 107, 109, 113, 119, 131,
       163, 167, 181, 183, 185, 203, 205, 207, 209, 215.

   - Up to 73 this is Session24's (KNOWN). Above 73 it is new.
   - Since `f(11) = 225`, orders ≤ 10 decide everything below 225.
   - **Birth against eternal colour.** Of the 109 attainable odd `h` with `3 ≤ h < 225`, 65 are blue at birth and 44
     black. Of the 44, 27 are eternally black. The other 17 turn blue later:
     - 25, 27, 31 at order 7;
     - 57, 73, 89 at order 8;
     - 97, 179, 187, 191, 195, 199, 211, 213, 219, 223 at order 9;
     - 121 at order 10, as `B₁₁ ⇒ B₁₁^op`.
   - 22 of the 44 are black at birth only because they are born at an even order with `h ≡ 3 (mod 4)` (THM-791). Ten of
     those 22 are later blue.
4. **From 225 to 625** (FINITE-EXACT + PROVED; `h = 377` OPEN). The eternally black values in `[225, 625)` are exactly

       231, 235, 239, 251, 267, 271, 291, 299, 301, 303, 311, 317, 321, 327, 331, 335, 347, 359,
       415, 427, 431, 439, 451, 463, 483, 487, 503, 515, 519, 523, 527, 563,

   with 377 undecided.
   - All 32 are born at order 8. Below 625 there are 61 eternally black values (with 7 and 21); 49 are `≡ 3 (mod 4)`.
   - Of the 309 values `3 ≤ h < 625`, 135 are blue at birth and 174 black at birth. The 174 split into 59 eternally
     black, 114 later blue, and 377.
   - **Delays.** For the eternally blue values the first self-converse order exceeds the birth order by:
     - 0 (135 values), 1 (95), 2 (5) or 3 (13);
     - 5 for `275 = 5²·11`. It has no self-converse realisation through order 12, and its first is
       `B₅ ⇒ C₁₁ ⇒ B₅^op` at order 13.
   - **Composite-only blue.** 121 and 275 are the only eternally blue values below 625 that are not the path count of
     any strong self-converse tournament (exact, since `B(121) = 9` and `B(275) = 11`).

## Proofs

**1.**
- **The condensation of a self-converse tournament is a palindrome** (THM-4866(E); Session24 PROOF eq. (2)).
  - Let `φ: T → T^op` be an isomorphism. It maps strong components to strong components. It also preserves the
    condensation order of `T^op`, which is the reverse of `T`'s.
  - So `φ(C_i) = C_(k+1−i)`, and `T[C_(k+1−i)] ≅ T[C_i]^op`. For odd `k` the central component is self-converse.
- **Hence the formula.** `H` is multiplicative over the condensation and `H(C^op) = H(C)`, so
  `H(T) = (∏_(i ≤ ⌊k/2⌋) H(C_i))² · H(C_centre)`. The last factor is 1 when `k` is even or the centre is a single
  vertex.
- **Conversely**, for any tournaments `B_i` and any self-converse `Z`, the palindrome displayed in 1 is self-converse.
  The anti-isomorphism exchanges the copies of `B_i` and `B_i^op` by the identity on labels, and acts on `Z` by an
  anti-automorphism.
- **The value sets.** `∏ H(B_i)` ranges over `Σ`, because `Σ` is closed under products (`A ⇒ B`) and contains 1. A
  centre that is not strong is itself a palindrome, and can be absorbed into the `B_i`.
- **Eternally black values.** If `h ∉ Σ_sc`, no tournament with `H = h` is self-converse, so every converse orbit has
  size 2.

**2.**
- **Busch.** A strong tournament with `H = e` has `f(|T|) ≤ e`, so `|T| ≤ B(e)`.
- **Parity of the centre.** For squarefree `h`, `d = 1`, so the centre is strong with `H = h`. Its order has the parity
  of `|T|` (THM-4866(E)). At even order THM-791 forces `h ≡ 1 (mod 4)`.

**3–4. The exact lists.**
- The self-converse classes of order `≤ 9` come from the census.
- **Order 10.** Every self-converse tournament has an anti-involution: `Aut(T)` has odd order, so an odd power of an
  anti-automorphism squares to 1. At even order it is fixed-point-free, and all such involutions are conjugate. So the
  `2^25` labelled tournaments anti-invariant under one fixed involution contain every self-converse class of order 10
  (`sc_even10.c`). Run at orders 6 and 8, the same program reproduces the census's self-converse value sets and strong
  flags exactly.
- **Order 11.** Let `σ` fix `v0`. Then `T − v0` is anti-invariant of order 10, and `H(T − v0) ≤ H(T)` by insertion
  monotonicity (THM-1370). So every self-converse 11-tournament with `H < 625` is one of the 32 anti-invariant
  extensions of an order-10 anti-invariant tournament with `H < 625` (`sc_odd11.c`).
  - `sc_odd11_fast.c` evaluates the extensions by splitting each path at `v0`. It was checked against direct DP on
    every extension of the first `2^17` codes, and the full run used it in four parallel parts of the code range. The merged output passes `check_order11.py`: every self-converse value of orders 3, 5, 7 and 9 below 625 recurs at order 11, as padding by a source and a sink requires. (The direct-DP program `sc_odd11.c` is kept as the reference implementation.)
- **The decision.** `sc_holes.py` applies 1 and 2:
  - for `e < 375`, `B(e) ≤ 11`;
  - for `e ≡ 3 (mod 4)` below 625, order 12 is even (excluded by THM-791) and order 13 needs `e ≥ f(13) = 625`;
  - for 275, orders ≤ 11 are complete, and order 12 is excluded because `275 ≡ 3 (mod 4)`. ∎

## Remarks

- **The owner's two kinds.** The birth colour is not stable: black at birth can turn blue later, often just because of
  the parity at an even least order. The eternal colour is stable.
- **Generating larger structures.** Each eternally blue value is realised by a palindrome built from smaller pieces
  (statement 1). The delay from birth to a self-converse realisation is not bounded by 3. It is 5 for 275, whose only
  self-converse realisations are sandwiches around the order-5 self-converse atom with `H = 11`.
- **THM-791's parity, proved.** At even order the anti-involution has no fixed vertex, so `C ↦ reverse(σ C)` pairs up
  the directed odd cycles. Hence `c_odd` is even, and OCF gives `H ≡ 1 + 2c_odd ≡ 1 (mod 4)` (THM-791(a); Session29
  section 4).
