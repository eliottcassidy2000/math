---
id: THM-4930
title: "Eternal colours of odd Hamiltonian-path counts. An odd h is the path count of some self-converse tournament (eternally blue) iff h = d^2 e with d a path count (or 1) and e = 1 or the path count of a strong self-converse tournament. Otherwise (eternally black) every tournament with H = h, at every order, has a distinct converse, so its classes come in converse pairs at every size. By Busch's minimum this is decidable value by value. Exactly 29 odd values below 225 are never self-converse: 7, 21 and 19, 23, 35, 39, 43, 47, 53, 55, 59, 63, 101, 103, 107, 109, 113, 119, 131, 163, 167, 181, 183, 185, 203, 205, 207, 209, 215. Of the 44 values below 225 that are black at birth (THM-4922), 17 turn blue later, 121 = 11^2 last, at order 10 as B11 => B11^op. Below 1131 the eternally black values are exactly 7, 21 and 81 attainable values: 27 below 225; 33 in [225, 625), all born at order 8 (the last one decided, 377, by a complete order-12 enumeration below 625); and 21 in [625, 1131), all born at order 9 (671, 683, 707, 719, 723, 755, 785, 795, 815, 851, 869, 887, 907, 931, 949, 983, 1003, 1009, 1067, 1087, 1115; complete self-converse inventories of orders 12 and 13 below 1125). Only four eternally blue values below 1131 have no strong self-converse realisation: 121 = 11^2 (first self-converse at order 10), 275 = 5^2 * 11 (first at order 13, five orders after birth), 711 = 3^2 * 79 and 855 = 3^2 * 95 (first at order 13, four orders after birth)."
status: >
  PROVED: statements 1-2 (the criterion; Busch decidability; parity).
  FINITE-EXACT, on complete self-converse inventories:
  - statement 3, from orders <= 9 (THM-4922's census generator), order 10 (all 2^25 anti-invariant labelled
    tournaments) and order 11 with H < 625 (all anti-invariant extensions);
  - statement 4: exact for every odd h < 625. The last undecided value, h = 377, is settled at order 12: audit K1's
    aud_sc12.c (section 1.6) enumerates every self-converse 12-tournament with H < 625 (117 values, 377 not among
    them), and this session's independent sc_even_pair.c, run with H <= 377, finds exactly the audit's values up to
    377.
  - statement 5 (added after audit K1; audit K2: PASS WITH CORRECTIONS, confirmed by an independent enumeration of
    orders 12 and 13 below 1125): exact for every odd h < 1131, from complete self-converse inventories of orders 12
    and 13 with H <= 1124 (sc_odd13.c, validated against every complete inventory of orders 6-11) and the complete
    order 11 (sc_odd13.c -DM=4, identical to audit K1's).
  FINITE: witnesses for values above the exact range.
  KNOWN, credited:
  - the palindromic condensation (THM-4866(E); THM-4743 (1b));
  - the self-converse generating function J_h(z) = sum_(d^2 e = h) F_d(z^2) C^+_e(z), with C^+_1 = 1 + z (Session24
    PROOF eq. (2); Session224 section 5);
  - the membership test and its Busch decidability (Session24 PROOF section 1: strong components have
    |C| <= 1 + 3 log_5 q; its "Completeness in H" paragraph);
  - the zero self-converse series for h in {7, 19, 21, 23, 35, 39, 43, 47, 53, 55, 59, 63} below 73 (Session24
    low_H_formulas.json);
  - the enumeration of all self-converse tournaments through one fixed anti-involution (THM-4743 (6), audit S49, all
    2^30 labellings at n = 11; existence and conjugacy of anti-involutions also in Session29's half-path and
    converse-congruences reports, section 1);
  - the order-10 self-converse value set (THM-4922 audit, aud_o10_byH.txt);
  - Busch's minimum f (THM-4744, CITED);
  - even-order self-converse tournaments have H = 1 (mod 4) (THM-791(a); Session29 converse-congruences report,
    section 3);
  - insertion monotonicity (THM-1370);
  - every odd value to 2^32 except 7, 21 is a path count (THM-4875, CITED);
  - birth colours (THM-4922).
  New here:
  - the eternal colour as the stable form of the owner's two kinds (J_h != 0 in Session24's notation);
  - the exact lists above 73;
  - the comparison of birth colours with eternal colours;
  - the strong flags of the order-10 self-converse values, and the order-11 self-converse values below 625;
  - statement 5 and the order-12 and order-13 self-converse inventories below 1125.
  Audit K1 (2026-11-08): PASS WITH CORRECTIONS, no mathematical error. Its credit and typing corrections are applied
  (process lesson: MISTAKE-738). Its supplementary order-12 enumeration settled 377.
  Audit K2 (2026-11-08, statement 5 and the other post-K1 additions): PASS WITH CORRECTIONS, no mathematical error;
  its own enumeration of orders 12 and 13 below 1125 reproduces sc_odd13_1124.out exactly. Corrections applied
  (typing, credit, records; addendum to MISTAKE-738).
session: mac-mini-2026-11-08-firstrealizers (owner directive 2026-11-08: first creators of odd values, one class or a converse pair, two kinds, the splitting pattern; prior work on 33)
source: 05-knowledge/results/eternal_colours_self_converse_values_20261108.md
depends_on: [THM-4866, THM-4743, THM-791, THM-1370, THM-4744, THM-4745, THM-4875, THM-4922]
verification:
  - 04-computation/experiments/eternal_colours_20261108/sc_holes.py (+ .out)
  - 04-computation/experiments/eternal_colours_20261108/sc_even10.c (+ sc_even10.out; -DN=6, 8 validation builds)
  - 04-computation/experiments/eternal_colours_20261108/sc_odd11.c, sc_odd11_fast.c (+ sc_odd11.out; -DCHECK self-test), check_order11.py (+ .out)
  - 04-computation/experiments/eternal_colours_20261108/sc_even_pair.c (+ sc_even12_377.out; -DM=3 validation build at order 8)
  - 04-computation/experiments/eternal_colours_20261108/sc_odd13.c (+ sc_odd13_1124.out, sc_full1011_le2001.out; -DM=2, 3, 4 validation builds; -DCHECK)
  - 04-computation/experiments/eternal_colours_20261108/sc_witness.c (+ .out), sc_walk.c (+ sc_walk.out)
  - 04-computation/experiments/eternal_colours_20261108/audit_K1/ (independent audit K1: complete orders 10, 11; order 12 below 625)
audit: 05-knowledge/results/eternal_colours_audit_K1_20261108.md (K1, PASS WITH CORRECTIONS); 05-knowledge/results/eternal_colours_audit_K2_20261108.md (K2, statement 5 and post-K1 additions, PASS WITH CORRECTIONS)
---

# THM-4930 — eternal colours: which odd path counts are ever self-converse

## Setting

- `H(T)` is the number of Hamiltonian paths of `T`, and `T^op` is the converse. Since `H(T^op) = H(T)`, the
  tournaments with a given `H` at a given order are a union of converse orbits, of size 1 (self-converse) or 2.
- `Σ` is the set of path counts. It contains every odd number up to `2^32` except 7 and 21 (THM-4875). `Σ_sc` is the
  set of path counts of self-converse tournaments, and `Σ_sc^str` that of strong self-converse tournaments.
- **Busch's minimum** (THM-4744, CITED): a strong `m`-tournament has `H ≥ f(m)`, with `f(3q) = 3·5^(q−1)`,
  `f(3q+1) = 5^q`, `f(3q+2) = 9·5^(q−1)`; so `f = 3, 5, 9, 15, 25, 45, 75, 125, 225, 375, 625, 1125, 1875` at
  `m = 3, …, 15`.
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
   - Every self-converse tournament of even order has `H ≡ 1 (mod 4)`: it has an anti-involution (`|Aut|` is odd, so
     an odd power of an anti-automorphism is one), with no fixed vertex at even order, and THM-791(a) applies (also
     Session29 converse-congruences report, section 3).
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
   - 22 of the 44 are forced black at birth by parity: they are born at an even order and `h ≡ 3 (mod 4)` (THM-791).
     Ten of those 22 are later blue.
4. **From 225 to 625** (FINITE-EXACT + PROVED). The eternally black values in `[225, 625)` are exactly

       231, 235, 239, 251, 267, 271, 291, 299, 301, 303, 311, 317, 321, 327, 331, 335, 347, 359,
       377, 415, 427, 431, 439, 451, 463, 483, 487, 503, 515, 519, 523, 527, 563.

   - All 33 are born at order 8. Below 625 there are 62 eternally black values (with 7 and 21); 49 are `≡ 3 (mod 4)`
     and 13 are `≡ 1 (mod 4)`.
   - 377 is the one value here that orders ≤ 11 do not decide (`377 ≡ 1 (mod 4)` and `B(377) = 12`). It is
     squarefree, and no self-converse 12-tournament has `H = 377` (complete enumeration below 625 by audit K1;
     confirmed here up to 377).
   - Of the 309 values `3 ≤ h < 625`, 135 are blue at birth and 174 black at birth. The 174 split into 60 eternally
     black and 114 later blue.
   - **Delays.** For the eternally blue values the first self-converse order exceeds the birth order by:
     - 0 (135 values), 1 (95), 2 (5) or 3 (13);
     - 5 for `275 = 5²·11`. It has no self-converse realisation through order 12, and its first is
       `B₅ ⇒ C₁₁ ⇒ B₅^op` at order 13.
   - **Composite-only blue.** 121 and 275 are the only eternally blue values below 625 that are not the path count of
     any strong self-converse tournament (exact, since `B(121) = 9` and `B(275) = 11`).
5. **From 625 to 1131** (FINITE-EXACT + PROVED; added after audit K1, confirmed by audit K2). The eternally black values in
   `[625, 1131)` are exactly

       671, 683, 707, 719, 723, 755, 785, 795, 815, 851, 869, 887, 907, 931, 949, 983, 1003, 1009, 1067, 1087, 1115.

   - All 21 are born at order 9; 17 are `≡ 3 (mod 4)` and four (785, 869, 949, 1009) are `≡ 1 (mod 4)`. Below 1131
     there are 83 eternally black odd values (with 7 and 21).
   - **Exactness.** For `625 ≤ e < 1125`, `B(e) = 13`, so orders ≤ 13 decide `e`. Orders ≤ 11 are complete for every
     `H ≤ 2001`, and orders 12 and 13 for every `H < 1125`. The values 1125, 1127 and 1129 are strong self-converse
     values already at order 9.
   - **The strong self-converse spectra.** Below 1125, order 12 has 41 strong self-converse values, from 405, and order
     13 has 26: 625, 675, 729, 775, 819, 825, 841, 867, 873, 881, 937, 945, 947, 961, 999, 1011, 1017, 1025, 1031,
     1035, 1049, 1053, 1083, 1089, 1111, 1123.
   - **Birth against eternal colour.** All 253 values in `[625, 1131)` are born by order 9. 189 are blue at birth and
     64 black. The 64 split into 21 eternally black and 43 later blue, with delays 1 (27 values), 2 (13) and 4 (3).
     The three delay-4 values are 711, 855 and 1011, all first self-converse at order 13. `1011 = 3·337` (squarefree)
     is born at order 9 and is first self-converse at order 13, as a strong atom.
   - **Composite-only blue.** Below 1131 exactly four eternally blue values are not the path count of any strong
     self-converse tournament: 121, 275, `711 = 3²·79` and `855 = 3²·95`. The last two are first self-converse at
     order 13, as `C₃ ⇒ Z ⇒ C₃^op` around the order-7 strong self-converse atoms `Z` with `H = 79` and `H = 95`.

## Proofs

**1.**
- **The condensation of a self-converse tournament is a palindrome** (THM-4866(E); Session24 PROOF eq. (2);
  THM-4743 (1b)).
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
- **Anti-involutions.** A self-converse tournament has an anti-involution: `|Aut(T)|` is odd, so if an
  anti-automorphism `α` has order `2s`, then `s` is odd and `α^s` is an involutive anti-automorphism. It fixes at most
  one vertex (two fixed vertices would have to be joined by an arc reversed by `α^s`), so at even order it is
  fixed-point-free, and THM-791(a) gives `H ≡ 1 (mod 4)` (also Session29 converse-congruences report, section 3).
- **Busch.** A strong tournament with `H = e` has `f(|T|) ≤ e`, so `|T| ≤ B(e)`.
- **Parity of the centre.** For squarefree `h`, `d = 1`, so the centre is strong with `H = h`. Its order has the parity
  of `|T|` (THM-4866(E)). At even order THM-791 forces `h ≡ 1 (mod 4)`.

**3–5. The exact lists.**
- The self-converse classes of order `≤ 9` come from the census.
- **Order 10.** Every self-converse tournament has an anti-involution (proof of 2). At even order it is
  fixed-point-free, and all such involutions are conjugate. So the `2^25` labelled tournaments anti-invariant under
  one fixed involution contain every self-converse class of order 10 (`sc_even10.c`; the method is THM-4743 (6) and
  its audit S49, which ran all `2^30` labellings at order 11). Run at orders 6 and 8, the same program reproduces the
  census's self-converse value sets and strong flags exactly.
- **Order 11.** Let `σ` fix `v0`. Then `T − v0` is anti-invariant of order 10, and `H(T − v0) ≤ H(T)` by insertion
  monotonicity (THM-1370). So every self-converse 11-tournament with `H < 625` is one of the 32 anti-invariant
  extensions of an order-10 anti-invariant tournament with `H < 625` (`sc_odd11.c`).
  - `sc_odd11_fast.c` evaluates the extensions by splitting each path at `v0`. It was checked against direct DP on
    all 32 extensions of each of the 7,469 parents with `H(T') < 625` among the codes up to `2^17` (audit K1 extended
    the check to all extensions of the first 4,097 codes, of any H). The full run used it in four parallel parts of
    the code range. The merged output passes `check_order11.py`: every self-converse value of orders 3, 5, 7 and 9
    below 625 recurs at order 11, as padding by a source and a sink requires. (The direct-DP program `sc_odd11.c` is
    kept as the reference implementation.)
  - Audit K1 enumerated all `2^25` reduced anti-invariant codes at order 11 (within-pair arcs fixed by conjugation,
    so every class occurs), with no H bound and no pruning (`audit_K1/aud_sc.c`); its values below 625 agree, and
    audit K2 reproduced its complete order 11 count for count.
- **Order 12, `H < 625`** (only 377 needs it). Remove a `σ`-pair `{a, b}`: `T − {a, b}` is anti-invariant of order
  10, and `H(T − {a, b}) ≤ H(T − b) ≤ H(T)` by insertion monotonicity twice. So every self-converse 12-tournament with
  `H < 625` extends an order-10 anti-invariant tournament with `H < 625` by `a` (`2^10` choices, pruned at `H < 625`),
  then by `b = σ(a)` with its forced arcs and the arc `ab`. Two independent implementations agree:
  - audit K1's `aud_sc12.c` (another labelling of `σ`, split evaluation at `a`, direct DP for the order-12 value;
    validated with no bound at orders 4-10);
  - this session's `sc_even_pair.c`, which uses no reduction: all `2^25` labelled order-10 bases (validated at order 8
    with `-DM=3` against the census; 6,000 order-12 values checked by direct DP, 2,000 in each of three workers).
  Audit K1 finds 117 self-converse values below 625, all `≡ 1 (mod 4)`, and 377 is not among them.
  `sc_even_pair.c`, run with `H ≤ 377`, finds exactly the audit's 66 values up to 377, and so does `sc_odd13.c`
  below.
- **Orders 12 and 13, `H < 1125`** (statement 5; `sc_odd13.c`). Take `σ` with pairs `(a_p, b_p)`, `p = 0, …, 5`, and
  fixed vertex `v0`. The transpositions `(a_p b_p)` commute with `σ`, so conjugating by them every within-pair arc may
  be taken as `a_p → b_p` (the reduction of audit K1). Then `T − v0` (order 12) and `T − v0 − {a_5, b_5}` (order 10)
  are anti-invariant, and by insertion monotonicity neither has larger `H`. The program extends each of the `2^20`
  reduced order-10 codes with `H < 1125` by `a_5` (`2^10` arc sets, Gray code, split evaluation), then by `b_5` (forced
  arcs), then by `v0` (`2^6` arc sets, Gray code, split evaluation), pruning at `H ≤ 1124` at every level.
  - **Validation.** The same program with 2, 3 and 4 base pairs and no pruning reproduces, value for value and with
    strong flags, the census's self-converse values at orders 6, 7, 8 and 9, `sc_even10.out` at order 10, and audit
    K1's complete unpruned order 11 (41,912 values, equal SHA-256 digests). Below 625 its order-12 values are audit
    K1's 117. All comparisons are recorded in `sc_odd13_validation.out` (`validate_sc_odd13.py`).
  - **Checks.** 18,000 sampled order-11, 12 and 13 values were recomputed by direct DP, and their tournaments checked
    anti-invariant, with no failure (run logs in `sc_odd13_validation.out`).
  - **Audit K2** enumerated orders 12 and 13 below 1125 with its own program (`audit_K2/aud_sc13.c`: another
    labelling, a full DP at order 12, Landau's strongness test) and found the same 210 and 433 values with the same
    strong flags.
- **The decision.** `sc_holes.py` applies 1 and 2 with the inventories' completeness bounds:
  - for `e < 375`, `B(e) ≤ 11`;
  - for `e ≡ 3 (mod 4)` below 625, order 12 is even (excluded by THM-791) and order 13 needs `e ≥ f(13) = 625`;
  - for `e ≡ 1 (mod 4)` in `[375, 625)`, only order 12 is added, and the order-12 inventory below 625 is complete;
  - for `625 ≤ e < 1125`, `B(e) = 13`, and orders 12 and 13 are complete below 1125 (order 11 is complete to 2001);
  - for 275, 711 and 855, every order up to 12 is complete at `h` and has no self-converse realisation, so the first is
    the order-13 sandwich. ∎

## Remarks

- **The owner's two kinds.** The birth colour is not stable: black at birth can turn blue later, often just because of
  the parity at an even least order. The eternal colour is stable.
- **Generating larger structures.** Each eternally blue value is realised by a palindrome built from smaller pieces
  (statement 1). The delay from birth to a self-converse realisation is not bounded by 3. It is 5 for 275, whose only
  self-converse realisations are sandwiches around the order-5 self-converse atom with `H = 11`, and 4 for 711, 855
  and 1011. A delay can come from a sandwich (275, 711, 855) or from a strong atom that first appears late (1011).
- **The least strong self-converse value per order** (KNOWN in substance). By the Moon–Yang theorem (THM-4745 (6);
  HYP-9356, KNOWN) the strong `m`-tournaments with `H = f(m)` are the optimal interval chains, and a chain is
  self-converse iff its composition is a palindrome (HYP-9356: "the self-converse ones are the palindromes"; THM-4744 (4)
  tabulates them). A palindromic optimal composition exists iff `m ≢ 0 (mod 6)`. At `m = 3q` the optimal multiset is
  `{2, 3^(q−1)}`, and the single 2 must sit in the middle, which needs `q` odd. At `m = 3q+1` the composition
  `(3, …, 3)`, and at `m = 3q+2` the composition `(2, 3, …, 3, 2)`, are palindromes. Hence `f_sc(m) = f(m)`
  iff `m ≢ 0 (mod 6)`, consistent with THM-791 since `f(6k) ≡ 3 (mod 4)`. The inventories give `f_sc(6) = 17` and
  `f_sc(12) = 405` (FINITE-EXACT).
- **THM-791's parity, proved.** At even order the anti-involution has no fixed vertex, so `C ↦ reverse(σ C)` pairs up
  the directed odd cycles. Hence `c_odd` is even, and OCF gives `H ≡ 1 + 2c_odd ≡ 1 (mod 4)` (THM-791(a); Session29
  converse-congruences report, sections 1 and 3:
  `05-knowledge/results/session29_framed_profiles_20261022/session29-converse-congruences/REPORT.md`).
