# Eternal colours of odd Hamiltonian-path counts: which values are ever realised by a self-converse tournament

2026-11-08, mac-mini-2026-11-08-firstrealizers. Owner directive: "this session merge in prior work on 33 and the idea
that for every odd number besides 7 and 21 we can keep track of the structure that first created them as we slowly
increase tournament size, and for each, they will either be associated with a single isomorphism class or exactly two
classes which are converses of each other. it seems possible to generate larger structures from the smaller ones once
deeper patterns are understood to be established, as the odd values minus the {7,21} holes are split into two kinds, and
the splitting pattern studied".

Canon:
- [THM-4930](../../01-canon/theorems/THM-4930-eternal-colours-which-odd-path-counts-are-ever-self-converse.md):
  the eternal colours; the exact lists below 375, and below 625 for `h ≡ 3 (mod 4)` (377 open).
- [HYP-9520](../hypotheses/HYP-9520-the-eternally-black-values.md): the eternally black set.

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
  - the blue share near 0.83.
- **Earlier:**
  - Session24 multiplicities (2026-10-21): the self-converse generating functions
    `D_h(z) = Σ_(d²e=h) F_d(z²)[1_(e=1) + Q_e(z)]` and complete low-H formulas. Their table gives self-converse
    generating function 0 exactly for `h ∈ {7, 19, 21, 23, 35, 39, 43, 47, 53, 55, 59, 63}` (`h ≤ 73`).
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

**The criterion** (PROVED; KNOWN in generating-function form, Session24 and Session224). `h` is eternally blue iff
`h = d²e`, where `d` is a path count (or 1) and `e` is 1 or the path count of a *strong* self-converse tournament.
- This is the palindrome: the strong components of a self-converse tournament read the same backwards up to converse.
- The construction "generates larger structures from smaller ones": `B ⇒ Z ⇒ B^op` is self-converse whenever `Z` is.

**Decidability** (PROVED). A strong tournament with `H = e` has order at most `B(e) = max{n : f(n) ≤ e}`, by Busch's
minimum `f = 3, 5, 9, 15, 25, 45, 75, 125, 225, 375, 625` at orders 3–13. So the eternal colour of each value is a
finite computation.
- For squarefree `h ≡ 3 (mod 4)`, only odd-order centres count, because even-order self-converse tournaments have
  `H ≡ 1 (mod 4)` (THM-791).

## 2. The data

- **Order ≤ 9:** all classes, with self-converse and strong flags (THM-4922's `classgen.c`).
- **Order 10:** all `2^25` labelled tournaments anti-invariant under one fixed-point-free involution (`sc_even10.c`).
  Validated: at orders 6 and 8 the same enumerator reproduces the census's self-converse values and strong flags
  exactly.
- **Order 11, `H < 625`:** the 32 anti-invariant one-vertex extensions of every order-10 anti-invariant tournament with
  `H < 625` (`sc_odd11.c`). This is complete by insertion monotonicity. Extensions are evaluated by splitting paths at
  the new vertex, checked against direct DP.
- **Above the exact range:** near-transitive self-converse witnesses (`sc_witness.c`). The transitive tournament is
  self-converse under the mirror, and every mirror-closed set of at most four reversed arc orbits is tried at orders
  ≤ 13.

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
- **What the owner's picture needs.** Half of the black-at-birth values are an artefact of THM-791's even-order parity.
  Of those, 10 turn blue at the next odd order and 12 never do.

## 4. Results from 225 to 625 (exact, except 377)

The order-11 enumeration (all self-converse 11-tournaments with `H < 625`) completes the strong self-converse
inventories:
- for every value below 375 (`f(12) = 375`);
- for every value `≡ 3 (mod 4)` below 625, since order 12 is even and `f(13) = 625`.

Its merged output passes `check_order11.py`: every self-converse value of orders 3–9 below 625 recurs at order 11, as
padding by a source and a sink requires.

**The 32 eternally black values in `[225, 625)`:**

    231, 235, 239, 251, 267, 271, 291, 299, 301, 303, 311, 317, 321, 327, 331, 335, 347, 359,
    415, 427, 431, 439, 451, 463, 483, 487, 503, 515, 519, 523, 527, 563

- All 32 are born at order 8.
- **377 is open.** It is `≡ 1 (mod 4)`, so a strong self-converse 12-tournament could realise it. Order 12 is not
  enumerated here, because the strong-ear inventory of THM-4744 needs numba.
- **Totals.** Below 625 there are 61 eternally black odd values (with 7 and 21). Of the 309 attainable odd values
  `3 ≤ h < 625`, 135 are blue at birth and 174 black at birth. The 174 split as 59 eternally black, 114 later blue,
  and 377.
- **Delays.** A late-blue value reaches its first self-converse realisation 1 (95 values), 2 (5) or 3 (13) orders
  after birth, except `275 = 5²·11`. That value has no strong self-converse realisation. It cannot be self-converse at
  order 12 (`≡ 3 mod 4`), so it is first blue at order 13, as `B₅ ⇒ C₁₁ ⇒ B₅^op`, five orders after birth.
- **Composite-only blue.** 121 and 275 are the only blue values below 625 without a strong self-converse realisation.

**Above 625** (FINITE). Near-transitive witnesses clear 652 of the 688 odd values in `[625, 2001)`. The remaining 36
(30 of them `≡ 3 mod 4`) are listed in HYP-9520, which conjectures that the black set has density zero, with
finiteness as an exploratory strengthening. The share of black values falls: 24%, then 16%, then at most 5.2%.

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

- **Order 12** for 377. This needs a strong-ear inventory with self-converse flags.
- **HYP-9520:** density zero (or finiteness) of the eternally black set, and witnesses for the 36 values in
  `[625, 2001)`.
- The delay from birth to blue is not bounded by 3: it is 5 at 275.

## 7. Files

`04-computation/experiments/eternal_colours_20261108/`:
- `sc_holes.py`: the decision and all lists;
- `sc_even10.c`: order 10, with `-DN=6/8` validation builds;
- `sc_odd11.c`, `sc_odd11_fast.c`: order 11 below 625. The fast variant has a `-DCHECK` self-test and was used in
  four parallel parts; `check_order11.py` checks it;
- `sc_witness.c`: near-transitive witnesses;
- `refine_dichotomy.py`: §5.

The census files `cls_n.txt` (`n ≤ 9`) are regenerated by THM-4922's `classgen.c` (`classgen 9 cls`, 20 s).
