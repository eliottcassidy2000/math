# Eternal colours of odd Hamiltonian-path counts: which values are ever realised by a self-converse tournament

2026-11-08, mac-mini-2026-11-08-firstrealizers. Owner directive: "this session merge in prior work on 33 and the idea
that for every odd number besides 7 and 21 we can keep track of the structure that first created them as we slowly
increase tournament size, and for each, they will either be associated with a single isomorphism class or exactly two
classes which are converses of each other. it seems possible to generate larger structures from the smaller ones once
deeper patterns are understood to be established, as the odd values minus the {7,21} holes are split into two kinds, and
the splitting pattern studied".

Canon:
- [THM-4930](../../01-canon/theorems/THM-4930-eternal-colours-which-odd-path-counts-are-ever-self-converse.md):
  the eternal colours; the exact lists below 1131 (377 settled at order 12; `[625, 1131)` by complete orders 12 and
  13 below 1125, added after audit K1, PRE-AUDIT).
- [HYP-9520](../hypotheses/HYP-9520-the-eternally-black-values.md): the eternally black set.
- Audit: [K1](eternal_colours_audit_K1_20261108.md), PASS WITH CORRECTIONS (credit and typing; no mathematical
  error). Its corrections are applied below.

Claims are typed PROVED, KNOWN, CITED, FINITE-EXACT, FINITE, HEURISTIC or NEGATIVE.

## 0. What was already done (credited, not re-derived)

The same prompt reached two lanes on 2026-11-07, and both were read before any claim here.
- **Session224** (`session224_first_occurrence_20261107/`):
  - the dichotomy is false at `h = 9`, which is born at order 5 with two self-converse classes;
  - the `(b, s, p)` birth ledger through order 8;
  - the merged H = 33 fiber: 5 birth classes (one self-converse class and two converse pairs), and at every order
    `n ≥ 6`, `n² − 5n − 1` classes, exactly one self-converse;
  - the strong-component grammar.
- **kind-pasteur-2026-11-07-S16**, THM-4922 / HYP-9512:
  - first creators are converse-closed and strong (strong through order 11);
  - the value colours at birth: **blue** if some first creator is self-converse, **black** otherwise;
  - the even-order black cell (THM-791) and the odd-order self-converse balance mod 4 (PROVED);
  - the genealogy (every first creator but the orphan pair of 63 grows from a smaller first creator, through order 10);
  - the blue share outside the even-order black cell, 82% in aggregate (0.825 at order 10).
- **Earlier:**
  - Session24 multiplicities (2026-10-21): the self-converse generating functions
    `J_h(z) = Σ_(d²e=h) F_d(z²) C⁺_e(z)` (with `C⁺_1 = 1 + z`) and complete low-H formulas. Its PROOF §1 already bounds
    strong components by Busch (`|C| ≤ 1 + 3 log_5 q`), and its "Completeness in H" paragraph declares the strong and
    strong self-converse inventories complete for every odd `h ≤ 73`. Its table gives self-converse generating function
    0 exactly for `h ∈ {7, 19, 21, 23, 35, 39, 43, 47, 53, 55, 59, 63}` (`h ≤ 73`).
  - THM-4743 (6) and its audit S49: all self-converse tournaments through `n = 11` via one anti-involution (all `2^30`
    labellings at `n = 11`; automorphism data, not H). This is the method of §2.
  - Session29's half-path and converse-congruences reports: existence and conjugacy of anti-involutions; even-order
    `H ≡ 1 (mod 4)`.
  - THM-4922's audit: the self-converse class count of every H at order 10 (`aud_o10_byH.txt`).
  - S50 / THM-4744: complete strong inventories below 625, and Busch's minimum `f`.
  - THM-4866(E): the palindromic condensation.
  - THM-791: even-order self-converse tournaments have `H ≡ 1 (mod 4)`.
  - THM-1370 (7 and 21 never occur; insertion monotonicity); THM-4875 (every odd value to `2^32` except 7, 21);
    klein THM-4889 (least orders, strong windows).
- **Concurrent Codex Sessions 225–228**: graft and recursive physical constructions toward Chan–Pak. No overlap with
  the self-converse question here.

## 1. The owner's two kinds, made stable

The birth colour of THM-4922 depends on the least order. A value can be black there and blue one order later, so it
does not split the odd values into two lasting kinds. The stable split is the **eternal colour**:
- **eternally blue** if some self-converse tournament, of any order, has `H = h`;
- **eternally black** otherwise. Then, at every order, every tournament with `H = h` comes with a distinct converse,
  so its classes always come in converse pairs. This is the owner's "exactly two classes which are converses of each
  other", in the only form that holds at every size.

**The criterion** (PROVED; KNOWN in generating-function form, Session24 and Session224: it says `J_h ≠ 0`). `h` is
eternally blue iff `h = d²e`, where `d` is a path count (or 1) and `e` is 1 or the path count of a *strong*
self-converse tournament.
- This is the palindrome: the strong components of a self-converse tournament read the same backwards up to converse.
- The construction "generates larger structures from smaller ones": `B ⇒ Z ⇒ B^op` is self-converse whenever `Z` is.

**Decidability** (KNOWN: Session24 PROOF §1 and its "Completeness in H" paragraph; restated here with the exact
`B(e)`). A strong tournament with `H = e` has order at most `B(e) = max{n : f(n) ≤ e}`, by Busch's
minimum `f = 3, 5, 9, 15, 25, 45, 75, 125, 225, 375, 625` at orders 3–13. So the eternal colour of each value is a
finite computation.
- For squarefree `h ≡ 3 (mod 4)`, only odd-order centres count, because even-order self-converse tournaments have
  `H ≡ 1 (mod 4)` (THM-791).

## 2. The data

- **Order ≤ 9:** all classes, with self-converse and strong flags (THM-4922's `classgen.c`).
- **Order 10:** all `2^25` labelled tournaments anti-invariant under one fixed-point-free involution (`sc_even10.c`;
  the method is THM-4743 (6) / audit S49, and the value set is in THM-4922's audit; new here are the strong flags).
  Validated: at orders 6 and 8 the same enumerator reproduces the census's self-converse values and strong flags
  exactly.
- **Order 11, `H < 625`:** the 32 anti-invariant one-vertex extensions of every order-10 anti-invariant tournament with
  `H < 625` (`sc_odd11.c`). This is complete by insertion monotonicity. Extensions are evaluated by splitting paths at
  the new vertex, checked against direct DP.
- **Order 12, `H ≤ 377`:** pair extensions of the order-10 anti-invariant tournaments (`sc_even_pair.c`). Removing a
  pair `{a, b} = {a, σ(a)}` leaves an order-10 anti-invariant tournament, and `H` does not increase (insertion
  monotonicity twice), so pruning at `H ≤ 377` is complete. Validated at order 8 (`-DM=3`) against the census.
- **Orders 12 and 13, `H < 1125`:** `sc_odd13.c`. With the audit's reduction (every within-pair arc `a_p → b_p`, by
  conjugating with the transpositions `(a_p b_p)`, which commute with `σ`), it extends each of the `2^20` reduced
  order-10 codes by `a` (`2^10` arc sets), `b = σ(a)` (forced) and the fixed vertex `v0` (`2^6` arc sets), with Gray
  codes and split evaluation, pruning at `H ≤ 1124` at every level. Run with 2, 3 and 4 base pairs and no pruning, it
  reproduces the census's self-converse values and strong flags at orders 6-9, `sc_even10.out`, and audit K1's
  complete order 11. Its complete order-10/11 output up to 2001 is kept as `sc_full1011_le2001.out`.
- **Audit K1** independently enumerated all self-converse tournaments of orders 10 and 11 (all `2^30` labellings, no
  pruning) and of order 12 below 625 (`audit_K1/`). Every value used here agrees.
- **Above the exact range:** near-transitive self-converse witnesses. The transitive tournament is self-converse under
  the mirror `v ↦ n−1−v`, and so is every tournament obtained by reversing whole mirror orbits of arcs.
  - `sc_witness.c` tries every set of at most four reversed orbits at orders ≤ 13.
  - `sc_walk.c` runs random walks over the reversed-orbit sets at orders 12–15, accepting a move iff `H` stays below a
    cap, and records every strong value met. At order 9 it recovers all 1,173 strong self-converse values of the
    census, and nothing else.

## 3. Results below 225 (exact)

The odd values below 225 that are **never** the path count of a self-converse tournament are 7 and 21 together with 27
attainable values:

    19, 23, 35, 39, 43, 47, 53, 55, 59, 63, 101, 103, 107, 109, 113, 119, 131, 163, 167, 181, 183, 185, 203, 205, 207, 209, 215

- Up to 73 this is Session24's list (KNOWN). From 101 on it is new.
- 20 of the 29 values are `≡ 3 (mod 4)`.

**Birth colour against eternal colour** (odd `h < 225`). There are 109 values; 65 are blue at birth and 44 black at
birth.

| black at birth | eternally black | blue later |
|---|---|---|
| parity-forced (even birth order, `h ≡ 3 mod 4`) | 12: 19, 23, 43, 63, 107, 119, 163, 167, 183, 203, 207, 215 | 10: 27, 31 (at 7); 179, 187, 191, 195, 199, 211, 219, 223 (at 9) |
| not forced | 15: 35, 39, 47, 53, 55, 59, 101, 103, 109, 113, 131, 181, 185, 205, 209 | 7: 25, 57, 73, 89, 97, 121, 213 |

- **How late values turn blue.** 16 of the 17 get a strong self-converse realisation one or two orders after birth.
  `121 = 11²` gets none by order 9 and is first blue at order 10, as the sandwich `B₁₁ ⇒ B₁₁^op` built from two copies
  of the order-5 first creator of 11.
- **What the owner's picture needs.** Half of the black-at-birth values are forced black by THM-791's even-order parity.
  Of those, 10 turn blue at the next odd order and 12 never do.

## 4. Results from 225 to 625 (exact)

The order-11 enumeration (all self-converse 11-tournaments with `H < 625`) completes the strong self-converse
inventories:
- for every value below 375 (`f(12) = 375`);
- for every value `≡ 3 (mod 4)` below 625, since order 12 is even and `f(13) = 625`.

Every other value below 625 is then decided by a realisation of order ≤ 11, except 377, which order 12 settles.

Its merged output passes `check_order11.py`: every self-converse value of orders 3–9 below 625 recurs at order 11, as
padding by a source and a sink requires.

**The 33 eternally black values in `[225, 625)`:**

    231, 235, 239, 251, 267, 271, 291, 299, 301, 303, 311, 317, 321, 327, 331, 335, 347, 359,
    377, 415, 427, 431, 439, 451, 463, 483, 487, 503, 515, 519, 523, 527, 563

- All 33 are born at order 8.
- **377 is eternally black.** It is `≡ 1 (mod 4)` and `B(377) = 12`, so a strong self-converse 12-tournament could
  have realised it. None does: audit K1 enumerated every self-converse 12-tournament with `H < 625` (117 values, all
  `≡ 1 mod 4`), and this session's `sc_even_pair.c`, a separate implementation run with `H ≤ 377`, finds exactly the
  audit's values up to 377. 377 = 13·29 is squarefree,
  so no sandwich realises it either.
- **Totals.** Below 625 there are 62 eternally black odd values (with 7 and 21). Of the 309 attainable odd values
  `3 ≤ h < 625`, 135 are blue at birth and 174 black at birth. The 174 split as 60 eternally black and 114 later blue.
- **Delays.** A late-blue value reaches its first self-converse realisation 1 (95 values), 2 (5) or 3 (13) orders
  after birth, except `275 = 5²·11`. That value has no strong self-converse realisation. It cannot be self-converse at
  order 12 (`≡ 3 mod 4`), so it is first blue at order 13, as `B₅ ⇒ C₁₁ ⇒ B₅^op`, five orders after birth.
- **Composite-only blue.** 121 and 275 are the only blue values below 625 without a strong self-converse realisation.

## 4b. Results from 625 to 1131 (exact; added after audit K1, PRE-AUDIT)

For `625 ≤ e < 1125`, `B(e) = 13`. Orders ≤ 11 are complete up to 2001 and orders 12 and 13 below 1125, so every
value below 1125 is decided; 1125, 1127 and 1129 are strong self-converse values at order 9.

**The 21 eternally black values in `[625, 1131)`:**

    671, 683, 707, 719, 723, 755, 785, 795, 815, 851, 869, 887, 907, 931, 949, 983, 1003, 1009, 1067, 1087, 1115

- All 21 are born at order 9: 17 are `≡ 3 (mod 4)`, and four (785, 869, 949, 1009) are `≡ 1 (mod 4)`.
- They are exactly the values of `[625, 1131)` that the witness searches had left unresolved.
- **Totals.** Below 1131 there are 83 eternally black odd values (with 7 and 21). On `[625, 1131)`, 189 of the 253
  values are blue at birth; the 64 black at birth split as 21 eternally black and 43 later blue.
- **Delays.** 1 (27 values), 2 (13) and 4 (3): `711 = 3²·79` and `855 = 3²·95`, first self-converse as
  `C₃ ⇒ Z ⇒ C₃^op` at order 13 around order-7 atoms, and 1011, a strong self-converse atom first found at order 13.
- **Composite-only blue below 1131:** 121, 275, 711, 855.
- **The strong self-converse spectra.** Order 13 has 26 strong self-converse values below 1125, from `625 = f(13)`;
  order 12 has 41, from 405. The least strong self-converse value of order `m` equals Busch's `f(m)` iff
  `m ≢ 0 (mod 6)`. This is KNOWN in substance: the Moon–Yang extremal chains (THM-4745 (6)) are self-converse iff
  palindromic. So the parity hole of THM-791 at orders 6 and 12 is visible already at the bottom of the spectrum
  (`f_sc(6) = 17`, `f_sc(12) = 405`).

**Above 1131** (FINITE). Explicit realisations clear 422 of the 435 odd values in `[1131, 2001)`. The remaining 13 are
listed in HYP-9520, which conjectures that the black set has density zero, with finiteness as an exploratory
strengthening. The share of black values falls: 24%, 16.5%, 8.3%, then at most 3.0%.

## 5. The dichotomy, refined (NEGATIVE)

THM-4922 showed that refining the first creators by the least number of 3-cycles or of Hamiltonian cycles does not
restore "one self-converse class or one converse pair". Further refinements were tested on all 1,519 values with birth
order ≤ 9 (`refine_dichotomy.py`):

| refinement of the birth fibre | values where the dichotomy then holds |
|---|---|
| none | 68 |
| largest automorphism group | 208 |
| smallest automorphism group | 74 |
| fewest 3-cycles | 476 |
| most 3-cycles | 282 |
| largest automorphism group, then fewest 3-cycles | 616 |

No natural refinement restores it. The eternal colour of §1 is the form of the owner's idea that does hold.

## 6. Next steps

- **Order 12 for 377:** done (audit K1; `sc_odd13.c`).
- **HYP-9520:** density zero (or finiteness) of the eternally black set. The 13 values in `[1131, 2001)` without a
  realisation need complete inventories above 1124 at orders 13 (and 15 for 1891), and at orders 12-14 for 1145 and
  1245.
- **The delay.** It is not bounded by 3: 5 at 275, and 4 at 711, 855 and 1011. Is it bounded at all?

## 7. Files

`04-computation/experiments/eternal_colours_20261108/`:
- `sc_holes.py`: the decision and all lists;
- `sc_even10.c`: order 10, with `-DN=6/8` validation builds;
- `sc_odd11.c`, `sc_odd11_fast.c`: order 11 below 625. The fast variant has a `-DCHECK` self-test and was used in
  four parallel parts; `check_order11.py` checks it;
- `sc_even_pair.c`: order 12 with `H ≤ 377` (pair extension; `-DM=3` validation build at order 8);
- `sc_odd13.c`: orders 12 and 13 with `H ≤ 1124` (`sc_odd13_1124.out`), and with `-DM=4` the complete orders 10
  and 11 (`sc_full1011_le2001.out`, lines with `H ≤ 2001`); `-DM=2, 3, 4` validation builds and `-DCHECK`;
- `sc_witness.c`: near-transitive witnesses (exhaustive up to four reversed orbits);
- `sc_walk.c`: random walks over reversed-orbit sets at orders 12–15;
- `refine_dichotomy.py`: §5;
- `audit_K1/`: audit K1's independent programs and outputs (orders 1-9 re-verified, orders 10 and 11 complete,
  order 12 below 625).

The census files `cls_n.txt` (`n ≤ 9`) are regenerated by THM-4922's `classgen.c` (`classgen 9 cls`, 20 s).
