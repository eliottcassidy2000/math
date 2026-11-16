# Conditional theorems: what famous conjectures can and cannot do for Chan–Pak, with 13, the triples, 67/69, and Vitali

Session klein-2026-11-16-wheels2. Scripts `scripts/conjugacy_triples.py`, `scripts/collatz_tree.py`; receipts
`conjugacy_triples.out`, `collatz_tree_1e6.out`.

## 0. The number 13, once more, and the neighbouring triples

The owner's identity is exactly THM-5254 (G): `12 | 24` puts 13 in the wheel (`ord_13(2) = 12`), in the denominator
of `B_24` (`(13 − 1) | 24`, von Staudt–Clausen) and among the pure-A Niemeier Coxeter numbers (`A_12^2`, `h = 13`).
One divisor lattice, three faces; THM-5255 adds that 13 is also the second prime whose covering wheel needs `a = 5`
(period 416).

The triples, in the square-sum reading of FENCE_AND_SQUARE.md (sums `a², b², c²` on `1…b² − 1`; receipt
`conjugacy_triples.out`): `(5, 12, 13)` and `(7, 24, 25)` are the `n = 2, 3` members of the consecutive-parameter ray
`(2n + 1, 2n(n + 1), 2n(n + 1) + 1)` of NUMBER.md, and `(5, 12, 13)` is also the Fibonacci pair `(3, 2)` of GOLDEN.md. Each
primitive triple with `b` even gives one Hamiltonian path on `1…b² − 1` whose two exposed endpoints are **`a²` and
`b²/2`** (`9, 8` for `(3, 4, 5)`; `25, 72` for `(5, 12, 13)`; `49, 288` for `(7, 24, 25)`): the "8 and 9" of the square-sum
problem are the general endpoints `b²/2` and `a²` at the first triple, and Catalan's `2³, 3²` is a coincidence of that
instance (the endpoints of the next triple are `25 = 5²` and `72`, not consecutive powers). The double `(10, 24, 26)` has
`d = gcd(100, 576) = 4`, so the orbit decomposition gives three pieces — paths on 143 and 144 vertices and one
288-cycle — exactly as the even-`d` clause predicts; doubling a triple is what first breaks the single path. The
"norm-69 Pell family" of the directive has no recorded counterpart in the repository (no file pairs 69 with a norm or a
Pell equation); 69 appears as the literal count of the 67/69 reader pair (`chanpak_pair_67_69_20261114`: `72Q + 67h`
against `72Q + 69h` from one arc reversal), as the first odd number from which the Fermat-block cuts cover every target
(`busch_staircase_fermat_block_cuts_20261113`: the odd non-cuts below 69 are exactly 1, 7, 11, 13, 21, 47, 67, and every
odd `N` in `[69, 398,129]` is a cut `H(A) + 2^j L(A, v)`), and as the third output `EAC` of the nested-recontact alphabet
whose matrices have determinant 24. 67 is the last of those non-cuts and the last exception of the near-join family
(HYP-9453) and of the `hc ≤ 2` path spectrum (HYP-9472): the two numbers are the boundary of the finite exceptional
set, not a Pell pair. The owner's observation "preserving the norm can destroy the tournament interface" is recorded
as a QUESTION: multiplying by a unit keeps a norm but not a sum graph, and no interface statement about 69 exists here.

## 1. Collatz cannot be the famous conjecture (PROVED, with a refuted prediction)

**Theorem (conjugacy).** Under `h = 2x + 1` the Collatz map `x ↦ 3x + 1` (`x` odd), `x ↦ x/2` (`x` even) becomes
`T(h) = 3h` for `h ≡ 3 (mod 4)` and `T(h) = (h + 1)/2` for `h ≡ 1 (mod 4)`. The trivial cycle `1 → 4 → 2 → 1` is
`h: 3 → 9 → 5 → 3`, and the orbit of `x = 3` is `h: 7 → 21 → 11 → 33 → 17 → 9 → 5 → 3`: **the forbidden Chan–Pak
pair `{7, 21}` is `{h(3), h(10)}`, the first two counts on the orbit of 3**, and every later count on it is realised
(this is COUNTS.md's remark "the plus sheet sends 3 → 10 → 5 to 7 → 21 → 11", stated as an exact conjugacy).

A conditional proof "Collatz ⟹ Chan–Pak" would pull realisability back along forward orbits: from `T(h) ∈ R` to
`h ∈ R`, by a dyadic lift `d ↦ 2d − 1` when `T(h) = (h + 1)/2` and by a ternary pullback `3h ∈ R ⟹ h ∈ R` when
`T(h) = 3h`. Both lifts are open in general (the first is HYP-9453's and INVERSE.md's `2d − 1` problem; the second is
not even constructive), but the structural obstruction is worse and exact:

**Theorem (the pair blocks the tree).** Every Collatz orbit that enters the final cycle through `5` passes through
`10`, because `10` is the only predecessor of `5`; in count coordinates every such orbit passes through `21`, where the
pullback fails (`21 = 2·11 − 1 ∉ R`). Only the orbits entering `16` from `32` (last odd element `(4^k − 1)/3`, i.e.
`h ∈ {43, 171, 683, …}` and their trees) avoid both 7 and 21. Numerically (receipt `collatz_tree_1e6.out`): of the
500,000 odd `h ≤ 10^6`, **373,393** have orbits through 7 or 21, and 187,200 of them are not one-apex counts. The
pre-registered P-C1 (that this exceptional set is a union of a few rays covered by the one-apex family) is
**refuted**; the exceptional set is three quarters of all odd numbers. So no pullback induction along Collatz orbits can
reach them, whatever lifting lemmas one assumes: the Collatz conjecture, true or false, does not bear on Chan–Pak, and
the matching of `{7, 21}` with the orbit of 3 is the whole relation. Conversely Chan–Pak says nothing about orbits.

## 2. A famous conjecture that does carry content: Erdős–Selfridge and the self-dual wheels (PROVED, conditional)

By THM-5255 (D), a self-dual wheel at period `2m` is a covering of `Z/m` by one class modulo `m′` for each prime of
order exactly `2m′` (`m′ | m`, `m′ > 1`), and the modulus 3 never occurs (`ν(6) = 0`).

**Theorem E.** Assume the Erdős–Selfridge conjecture (there is no covering system of the integers with distinct odd
moduli, all greater than 1). Then no wheel `W(2m)` is self-dual with at most one prime per order used, i.e. with
pairwise distinct flipper odd parts. *Proof.* A self-dual cover with distinct parts is a covering of `Z/m`, hence of
`Z` with period `m`, by classes with distinct odd moduli `m′ > 1`. ∎ Conversely, **a self-dual wheel with distinct
flipper parts, if one is ever found, refutes the Erdős–Selfridge conjecture**: the search of THM-5255 (D) for
self-dual periods is a computational attack on a famous open problem, and this is the "other way around" the owner
asked for.

Unconditionally, the known theorems constrain any self-dual wheel with distinct parts: by Hough–Nielsen (every
distinct covering system has a modulus divisible by 2 or 3) some flipper part is divisible by 3, and since 3 itself is
absent that part is at least 9; by Balister–Bollobás–Morris–Sahasrabudhe–Tiba (no distinct covering system with all
moduli odd and square-free) some flipper part is divisible by the square of an odd prime. So a distinct-part self-dual
wheel needs orders `2m′` with `9 | m′` or `25 | m′` or `49 | m′` …, i.e. primes of order 18, 50, 54, 90, 98, 150, …:
`19, 251, 4051, 87211, …`. With repeated parts allowed (several primes of the same order) the conjecture does not
apply, and the question of THM-5255 (D) stays open on its own terms.

## 3. What a conditional proof of Chan–Pak would need (typed)

The repository's two unconditional reductions are the signed descent (COUNTS.md, THM-5262: Chan–Pak follows from the
three inverse lifts `2d − 1`, `(4d − 3)/3`, `(4d + 1)/3` on realised `d`) and the Shallow-Cut programme (every odd
`N ≥ 69` as a Fermat-block cut `H(A) + 2^j L(A, v)`, verified to 398,129). Each would be settled by a lifting lemma
about tournaments, not by an arithmetic conjecture; the Collatz skeleton (§1) is strictly weaker than the signed
descent, whose forward map decreases unconditionally. The arithmetic conjectures that enter the subject at all enter
through covering systems (§2) and through the distribution of orders of 2 (Artin's conjecture would give, for a
positive density of primes `p`, a full clock `2^l (mod p)` for the residue-prescription of TWO_APEX_HOST, as the
mod-9 clock is the case `p = 3²`). Verdict: no famous open conjecture is known to imply or to follow from Chan–Pak;
the honest conditional theorem is Theorem E, about the wheels.

## 4. Bold hypotheses (CONJECTURES, labelled)

- **B1.** Self-dual wheels exist, and every one of them uses either a repeated flipper part or a part divisible by an
  odd square (forced by Erdős–Selfridge if true, and by BBMST unconditionally for the square-free case).
- **B2.** `a(q) = 6` for every prime `q ≥ 17`: no odd prime beyond 13 reaches `q` primes by level 5.
- **B3.** For composite odd `m` the odd wheel covers `Z/m` iff a box-covering inequality holds coordinate-wise on
  `∏ Z/q_i^(e_i)`; the failures at `m = 45, 33, 75` with density above 1 are the first test cases.
- **B4.** The exceptional sets of every one-port family are initial segments of `{1, 7, 11, 13, 21, 47, 67}`, and 67 is
  the last exception of any family that reaches every residue class modulo 8 (the non-cuts of the staircase, the
  near-join exceptions and the `hc ≤ 2` spectrum all end at 67).
- **B5.** The address-level exact `+2` partners of the synchronized reader have density exactly `2^(−k)` with `k` fixed
  by the reader (observed 2.2–2.9 % at all widths 12–22), and their interior fraction decays like `n^(−1)`.

## 5. Determinants and Vitali (METHOD)

The nested-recontact alphabet (`chanpak_absolute_steering_20261115/NESTED_RECONTACT.md`) is used precisely because its
letter matrices have determinant 24: a word of those letters has determinant a power of 24 times the determinants of the
other letters, never zero, so projective states never collapse. Determinants multiply, so "controlling a growing
range of nonzero determinants" cannot mean hitting every integer: the set of word determinants is a multiplicative
semigroup generated by the letter determinants. What can be controlled is the *bracket* of a reader pair,
`det [[I, Z], [I′, Z′]] = (I − I′)Z = ±2Z, ±4Z`, which grows with the slope `Z = 82·2^m`; the exact-gap theorem of
SWITCH.md makes every even gap up to `V` available at anchors of scale `V²`. The Vitali step is the one PROOF.md uses:
from a family of overlapping count windows (triples, five-count banks, gap windows) a greedy selection by increasing
left endpoint keeps a disjoint subfamily whose number is at least a third of the family (one-dimensional Vitali
covering lemma, constant 3 for intervals of comparable length), which is how `Ω(X^(2/3))` *disjoint* windows are
obtained from the collision count. The same lemma applied to the gap anchors of SWITCH.md at spacing `Θ(V)` gives
`Ω(V)` disjoint exact-gap windows of width `V` in `[V², 2V²]`; applied across dyadic scales it gives, for each even gap
`g ≤ V`, disjoint anchors at every scale `2^(2j) ≥ g²`. These are restatements of proved bounds in the covering-lemma
language, not new coverage.
