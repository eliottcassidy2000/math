# The Paley heptagon and the trivial Collatz cycle are one binary word: real and 2-adic readings give the two orientations, the sheet swap is the arc reversal, the Frobenius survives and the translations are what the carry kills; the owner's map F decoded; the bounded-approximant statement attempted

**Session:** opus, `collatz-functional-uniqueness-20261001` (S15 series, nineteenth note), 2026-10-01.
**Owner's directive (verbatim core):** "try proving the bounded 3-adic approximants statement, the lonely-runner extremal object (the Paley heptagon) is maximally symmetric, while the Collatz graph has no symmetry at all. That {1,2,4} is both the Paley connection set and the trivial cycle is not numerology and this angle should be investigated deeply. Consider the following function which is isomorphic to collatz: For even numbers 2N, they are sent to 3N and for odd numbers 2N-1, they are sent to K_N where K=0,2,1,4,2,10,3,9,4,15,5... where N in 1,2,3,4... consider that 2^{0,2,1} = {1,4,2}, the K sequence appears to have 2 copies of each element of N plus one 0, and the interval between the two occurrences of 2 in K is the microcosm".

**Status:**

* **PROVED (elementary):** Propositions 1–6, and the non-isomorphism, decoding and folding statements of §5.
* **FINITE-EXACT:** every check, in script `04-computation/experiments/collatz_paley_bridge_20261001.py` (output `.out` beside it, ending `ALL CHECKS PASSED`).
* **CITED:** the 2-adic conjugacy (Bernstein–Lagarias, Terras), Singer difference sets, Hooley.
* **Collatz OPEN.** The bounded-approximant statement is not proved (§6).
* **Independent audit OWED.**
* **Correction of my own label:** in the seventeenth note I typed the coincidence `{1,2,4}` = Paley set = trivial cycle as NUMEROLOGY. The owner is right that it is structural. The structure is Proposition 2, and that label is withdrawn there.

**Inherits:**

* The seventeenth note (THM-4523 rigidity; Proposition EH; Theorem R_a) and the eighteenth note (Theorem B; the quantifier exchange).
* The fourth note's Proposition 4: the trivial cycle is exactly linear, `T_1^(2j)(1 + 4^j t) = 1 + 3^j t`.
* The eleventh note's unit-minus-one clocks (Mersenne and Lucas counts).
* The repo's Paley/Fano/octonion note `collatz_mod6_20260922_paley_fano_octonion_design.md`: `dev{0,1,3}`, `|Aut(T_7)| = 21`.
* The flip-rank thread, in which the Paley heptagon is the `argmax |Aut|` obstruction (HYP-3805, HYP-3817).

## 0. Answers in brief

* **The Paley angle (new, exact).**
  * **One word, two readings.** The trivial cycle `1 → 4 → 2 → 1` has parity word `100`. Its three rotations, read as binary fractions, are `{1, 2, 4}/7`: the quadratic residues mod 7, i.e. the Paley connection set and a Fano line. Read 2-adically, through the conjugacy that turns Collatz into the shift, they are `−{1, 2, 4}/7`, whose numerators are the non-residues `{3, 5, 6}`.
  * **The other sheet.** The minus sheet's 3-cycle `{−5, −7, −10}` (shortcut map; equivalently `3n−1`'s `{5, 7, 10}`) has the complementary word `110` and the opposite codes.
  * **One involution.** Reversing all Paley arcs (`x ↦ −x`, an anti-automorphism since `−1` is a non-residue mod 7) is at once the place swap (real ↔ 2-adic) and the sheet swap (`3n+1` ↔ `3n−1`, bit complement, the reflection `x ↦ −1 − x`).
  * **Why length 3.** A single cycle code can be a full Paley connection set only for cycle length 3, because `2^(L−1) − 1 = L` forces `L = 3`.
  * **The symmetry contrast, made precise.** Of the 21 symmetries of the Paley heptagon, the stabiliser of `{1, 2, 4}` is the Frobenius `x ↦ 2x` (order 3). It is exactly the Collatz map moving along the trivial cycle. The 7 translations move `{1, 2, 4}` to the other Fano lines, none of which is a cycle code. These are what the carry `+1` destroys, and the rigidity theorem is their global absence.
* **The lonely-runner reading (exact).** The trivial cycle `{1, 2}` of the shortcut map, taken as a speed set, is a tight lonely-runner instance (`lon = 1/3`). Its lonely times are exactly `1/3, 2/3`, its own real parity code.
* **Your F.**
  * **As written, F is not isomorphic to Collatz.** If `K` lists every positive integer twice, as you observe, multiples of 3 get three preimages (e.g. `9 ← 6, 15, 37`), and no Collatz map has more than two. Even from the eleven given values, the 2-cycle `{2, 3}` has two preimages at each vertex, while every integer Collatz 2-cycle has in-degrees `(1, 2)`.
  * **The even rule `2N ↦ 3N` is Collatz exactly.** It is the odd branch seen from the shifted label `n + 1`. Completing it with `2N − 1 ↦ N` gives the `3n+1` map; with `2N − 1 ↦ N − 1` it gives the `3n−1` map.
  * **"Two copies of each element plus one 0"** is exactly the union of these two completions. It is the two sheets folded about `−1/2`, where Collatz becomes the two-branch relation `y ↦ ⌊y/2⌋, ⌈3y/2⌉`.
  * **Your microcosm is the trivial cycle twice.** `K_1..K_3 = 0, 2, 1` are its exponents, `2^{0,2,1} = {1, 4, 2}`, and `K_2..K_5 = 2, 1, 4, 2` is the cycle itself.
  * **Open question.** The rule behind the even-position values `2, 4, 10, 9, 15` is not recoverable from the data (§5).
* **Bounded approximants.** Not proved; the statement is equivalent to Collatz (eighteenth note, Proposition 2). §6 records what the new structure adds (a Paley-code form and the place exchange at the trivial cycle) and why it does not close.

## 1. The parity code of a cycle (PROVED; CITED conjugacy)

Let `T(n) = n/2, 3n + 1` and `T_1(n) = n/2, (3n + 1)/2`. The parity map
`Φ(x) = Σ_i (f^i(x) mod 2) 2^i` (`f = T` or `T_1`) conjugates `f` on `Z_2`
to the shift `σ(y) = (y − (y mod 2))/2`. For `T_1` this is the full 2-shift
(Bernstein–Lagarias). For `T` it is the golden-mean shift: every 1 is
followed by a 0, since `3x + 1` is even.

**Proposition 1.** Let a cycle of `f` have parity word `w = (w_0, ..., w_(L−1))`,
and put `c(w) = Σ w_i 2^i` and `r(w) = Σ w_i 2^(L−1−i)` (the reversed word).

1. **2-adic reading.** The parity image of the cycle point with word `w` is
   `Φ = −c(w)/(2^L − 1)`. Over the cycle, the numerators `−c` form one coset
   of the multiplicative group `<2>` in `Z/(2^L − 1)`, and `f` acts on them as
   multiplication by `2^(−1)`.
2. **Real reading.** The binary fraction `0.(w)_2 = r(w)/(2^L − 1)`. Over the
   cycle, the numerators form the coset of `r(w)`, and `f` acts on them as
   doubling on `R/Z`.

*Proof.* `Φ = Σ_k c 2^(kL) = c/(1 − 2^L)` in `Z_2`, and
`0.(w)_2 = Σ_(k ≥ 1) r 2^(−kL) = r/(2^L − 1)` in `R`. A rotation of the word
multiplies the numerators by 2 or by `2^(−1)` modulo `2^L − 1`. ∎

The same rational `c/(2^L − 1)` thus appears with opposite signs at the two
places, `R` and `Z_2`. This is the per-cycle form of the eleventh note's
Mersenne clock. Every integer cycle on `Z` is listed in `.out`, section A:

| map | cycle | `L` | word | real numerators | 2-adic numerators |
|---|---|---|---|---|---|
| `T` | `1, 4, 2` | 3 | `100` | `{1, 2, 4}` = QR_7 | `{3, 5, 6}` = NQR_7 |
| `T` | `−1, −2` | 2 | `10` | `{1, 2}` mod 3 | `{1, 2}` mod 3 |
| `T` | `−5, −14, −7, −20, −10` | 5 | `10100` | `{5, 9, 10, 18, 20}` mod 31 | `{11, 13, 21, 22, 26}` |
| `T` | `−17, …` | 18 | `101010100101010000` | 18 elements mod `2^18 − 1` | 18 elements |
| `T_1` | `1, 2` | 2 | `10` | `{1, 2}` mod 3 | `{1, 2}` mod 3 |
| `T_1` | `−5, −7, −10` | 3 | `110` | `{3, 5, 6}` = NQR_7 | `{1, 2, 4}` = QR_7 |
| `T_1` | `−17, …` | 11 | `11110111000` | 11 elements mod 2047 | 11 elements |

## 2. Length 3 is the Paley level (PROVED)

**Proposition 2 (the bridge).**

1. **Plus sheet.** The trivial cycle `{1, 4, 2}` of `T` has real code
   `QR_7 = {1, 2, 4}`, the connection set of the Paley tournament `P_7`, and
   2-adic code `NQR_7 = {3, 5, 6}`.
2. **Minus sheet.** The 3-cycle `{−5, −7, −10}` of `T_1`, equivalently the
   cycle `{5, 7, 10}` of `n ↦ n/2, (3n − 1)/2`, has real code `NQR_7` and
   2-adic code `QR_7`.
3. **One involution.** Three operations agree on these codes:
   * the arc reversal of `P_7`, `x ↦ −x`, an anti-automorphism because
     `−1 ∈ NQR_7`;
   * the place swap, real ↔ 2-adic, which changes the sign;
   * the sheet swap: the bit complement `w ↦ 1 − w`, i.e. `y ↦ −1 − y` in
     `Z_2`, which conjugates `T_1` to its reflection `x ↦ −1 − x`.

   Each exchanges `QR_7` and `NQR_7`.
4. **Two maps, two clocks.** The word `100` is the integer cycle `{1, 4, 2}`
   for `T`, with clock `2^2 − 3 = 1`. For `T_1` it is the rational cycle
   `{1/5, 4/5, 2/5}`, with clock `2^3 − 3 = 5`. Its Mersenne clock is
   `2^3 − 1 = 7`. On `{1, 4, 2}` the map `T` acts as `x ↦ 4x (mod 7)`, i.e.
   halving, which is a Paley automorphism.

*Proof.* Proposition 1 with `c(100) = 1`, `r(100) = 4` and `c(110) = 3`,
`r(110) = 6`. Then `<2> = {1, 2, 4} = QR_7` modulo 7 and
`−QR_7 = {6, 5, 3}`. The complement satisfies `c ↦ 2^L − 1 − c ≡ −c`. ∎

**Proposition 3 (why only length 3).** Let `p = 2^L − 1` be a Mersenne prime.
A single coset of `<2>` has `L` elements and `QR_p` has `2^(L−1) − 1`, so
they coincide only if `2^(L−1) − 1 = L`, i.e. `L = 3`. For `L = 3`, indeed
`<2> = QR_7`. For `L = 2, 5, 7, 13, ...` the coset is a proper part of
`QR_p` (`.out`, section C, to `L = 31`). So the trivial cycle's length
3 (two halvings and one odd step) is exactly what makes its code a whole
Paley connection set. Among the known integer cycles on `Z` (conjecturally
all of them), the only ones whose codes are non-trivial difference sets are
the two 3-cycles of Proposition 2.

**Proposition 4 (Fano and Singer).** `{1, 2, 4}` is a (7, 3, 1) difference
set: the Singer set of `PG(2, 2)`. It is the set of exponents `i` with
`Tr(α^i) = 0` for a root `α` of `x^3 + x + 1` over `F_2`, i.e. the Frobenius
orbit `α, α^2, α^4` of `α`. Its 7 translates are the lines of the Fano plane:
`dev{0,1,3}` of the Paley/Fano note, whose block `(1, 2, 4)` this is. So the
Collatz dynamics on the trivial cycle is the Frobenius of `F_8` acting on
exponents. FINITE-EXACT: `.out`, section D.

## 3. Maximal symmetry against rigidity, made precise (PROVED)

**Proposition 5.**

1. **The symmetry group.** `Aut(P_7) = {x ↦ ax + b : a ∈ QR_7, b ∈ Z/7}` has
   order 21. It is the Frobenius group `C_7 ⋊ C_3`, the point stabiliser of
   `PSL(2, 7)` (the flip-rank thread's apex group).
2. **The part the cycle keeps.** The setwise stabiliser of the connection set
   `QR_7` in `Aut(P_7)` is the multiplier group `<×2> ≅ C_3`. Through the
   codes of Proposition 1 it acts on `QR_7` exactly as `T` acts on the
   trivial cycle.
3. **The part the carry kills.** No non-trivial translate `QR_7 + b` is a
   coset of `<2>`, so none is the code of any cycle of any affine
   parity-shift map. Translations have no Collatz counterpart.

*Proof.* The automorphisms are checked exhaustively (21, with multipliers
`{1, 2, 4}`). A translation fixing a 3-set of `Z/7` is trivial, since 3 does
not divide 7, and the translates were checked one by one (`.out`,
section E). ∎

**Reading.** The Paley heptagon's symmetry splits as translations (7) ⋊
Frobenius (3). Collatz realises the Frobenius part as its own dynamics on
the trivial cycle; that is the only non-trivial self-map of a Collatz
structure. The translation part is broken by the carry, the `+1` that also
fixes the Eckmann–Hilton point `−1/5`. THM-4523 (rigidity) is the global
statement that nothing of it survives anywhere in the graph. So the
contrast "maximally symmetric against no symmetry" is the same object seen
twice: the cyclic word `100` with all its affine symmetries mod 7, and the
cyclic word `100` embedded in a dynamics whose carry kills every symmetry but
the shift along the cycle.

## 4. The lonely-runner reading (PROVED + FINITE-EXACT)

**Proposition 6.**

1. **Speeds `{1, 2}`.** These are the elements of the trivial cycle of
   `T_1`. They form a tight lonely-runner instance, `lon = 1/3 = 1/(k+1)`.
   The lonely times are exactly `t = 1/3, 2/3`, i.e. `0.(01)_2` and
   `0.(10)_2`: the real codes of that same cycle. Its 2-adic codes are
   `−1/3 = Φ_(T_1)(1)` and `−2/3 = Φ_(T_1)(2)`.
2. **Speeds `{1, 2, 4}`.** At `t = 1/7` the runners sit on `QR_7/7`, and at
   `t = 3/7` on `NQR_7/7`. The loneliness is `1/3`, not tight for `k = 3`,
   again attained at `t = 1/3`.

*Proof.* `min(‖t‖, ‖2t‖) ≤ 1/3` for every `t`: if `‖t‖ > 1/3` then
`‖2t‖ < 1/3`. Equality holds exactly at `t = 1/3, 2/3`. The grid check is
in `.out`, section F. ∎

So the lonely-runner problem (the real circle) and the Collatz problem (the
2-adic integers) read the same periodic parity words. The tight instance
`{1, 2}` is lonely precisely at its own Collatz code: a small exact instance
of the owner's thesis that the two problems are local and global faces of
one structure. The repo's stronger phrase "the Paley heptagon is the LRC
extremal object" is the flip-rank thread's typed claim (HYP-3805,
HYP-3802 "LRC atoms"); it is not re-proved here.

## 5. The owner's map F (PROVED: decoding and non-isomorphism; the even-position rule unresolved)

The owner's map is `F(2N) = 3N`, `F(2N − 1) = K_N`, with
`K = 0, 2, 1, 4, 2, 10, 3, 9, 4, 15, 5, ...`. The odd positions show
`K_(2j+1) = j`, i.e. `F(4j + 1) = j`. The even positions show
`2, 4, 10, 9, 15`.

* **Not isomorphic to Collatz as written.** Every Collatz map (`T` or `T_1`,
  on `N`, `Z` or the rationals with odd denominator) has in-degree at most 2.
  * **The extended pattern.** With the odd-position pattern extended,
    `9 ← 6, 15, 37` and `15 ← 10, 19, 61`: in-degree 3. If `K` lists every
    positive integer twice, every multiple of 3 has three preimages.
  * **The given data alone.** The 2-cycle `{2, 3}` has `2 ← 3, 9` and
    `3 ← 2, 13`. Every integer Collatz 2-cycle (`{1, 2}` for `T_1`,
    `{−1, −2}` for `T`) has in-degrees `(1, 2)`, and `T` has no 2-cycle on
    `N`.
* **The even rule is Collatz exactly.** `F^+(2N) = 3N`, `F^+(2N − 1) = N` is
  `T_1` in the labels `n + 1` (the `3n+1` sheet). `F^−(2N) = 3N`,
  `F^−(2N − 1) = N − 1` is `n ↦ n/2, (3n − 1)/2` in the labels `n − 1` (the
  `3n−1` sheet). Checked to `N = 2000`. Among relabellings by a translation
  `x ↦ x + c` these are the only two completions: `c = +1` for `3n+1` and
  `c = −1` for `3n−1`. `K` agrees with `F^+` at `N = 2, 4` and with `F^−` at `N = 1`.
* **"Two copies of each element plus one 0".** This is exactly the multiset
  of the odd-label values of `F^+` and `F^−` together. In other words, it is
  the two sheets glued: folding `Z` about `−1/2` (`z ↦ z` for `z ≥ 0`,
  `z ↦ −1 − z` otherwise) turns `T_1` into the two-branch relation
  `y ↦ ⌊y/2⌋` and `y ↦ ⌈3y/2⌉` (checked on `|z| < 500`). That is Mahler's
  `⌈3g/2⌉` (S15b) beside halving: the owner's "two identical
  multiplications" (×2 and ×3/2 with rounding), one per sheet.
* **The microcosm.** `K_1, K_2, K_3 = 0, 2, 1` are the exponents of the
  trivial cycle in order, `2^0 → 2^2 → 2^1`. `K_2, ..., K_5 = 2, 1, 4, 2` is
  the cycle `2 → 1 → 4 → 2` itself. In exponents the cycle is the rotation
  `j ↦ j − 1` on `Z/3`, and `2^j mod 7` turns that rotation into `<2> = QR_7`
  (Proposition 2). Exponent and value thus meet exactly in the Paley set.
  The odd step `1 → 4 = 2^2` is the carry closing the rotation
  (`3·2^0 + 1 = 2^2`, clock `2^2 − 3 = 1`).
* **Unresolved.** No simple rule (affine, two-step, folded) reproduces
  `K_6 = 10`, `K_8 = 9`. The owner's rule for the even positions is needed
  to say what `F` is.

## 6. The bounded-approximant statement (attempted; not proved)

The statement is: for every `n`,
`δ_k(n) = min{m ∈ Gamma_1 : m ≡ n (mod 2·3^k)}` is bounded in `k`. It is
equivalent to Collatz (eighteenth note, Proposition 2), so a proof would be
a proof of the conjecture. What the new structure adds:

* **Paley form.** By Lagarias, `n ∈ Gamma_1` iff the parity code of `n`
  eventually equals that of the trivial cycle, i.e. iff `Φ_T(n)` eventually
  lies in `−QR_7/7`. So the conjecture says every positive integer's 2-adic
  code eventually enters the Paley set `NQR_7/7` and stays there under the
  shift. This is a reformulation, not a reduction.
* **The place exchange at the trivial cycle.** `T_1^(2j)(1 + 4^j t) = 1 + 3^j t`
  (fourth note): 2-adic proximity to the cycle (radius `4^(−j)`) is
  converted affinely, with size decreasing, into 3-adic proximity to 1
  (radius `3^(−j)`). So bounded approximants propagate along the 2-adic
  shadow of the trivial cycle, but that shadow is only the classical
  descent class `n ≡ 1 (mod 4)`.
* **Why it does not close.**
  * The Paley bridge is a statement about *cycles*: finite words, the cycle
    half. Bounded approximants concern the whole *component*: the
    divergence half.
  * The bridge sees the sheet (QR against NQR), which rigidity did not, but
    it is drift-blind: the same dictionary codes `5n+1`'s cycles (cosets
    modulo `2^7 − 1` for its trivial cycle `1010000`).
  * Theorem B of the eighteenth note therefore still applies in its drift
    half, and no argument from the codes alone can bound approximants.

## 7. Verdicts

| claim | status |
|---|---|
| Proposition 1 (cycle codes = cosets of `<2>` mod `2^L − 1`; real and 2-adic readings differ by sign and reversal) | PROVED + FINITE-EXACT |
| Proposition 2 (trivial cycle: real QR_7, 2-adic NQR_7; minus-sheet 3-cycle reversed; arc reversal = place swap = sheet swap) | PROVED + FINITE-EXACT |
| Proposition 3 (a single 2-coset = QR_p only for `L = 3`) | PROVED + FINITE-EXACT to `L = 31` |
| Proposition 4 (Singer/Fano: trace-zero exponents, Frobenius orbit) | PROVED + FINITE-EXACT |
| Proposition 5 (stabiliser of QR_7 = Frobenius = cycle dynamics; translations have no counterpart) | PROVED + FINITE-EXACT |
| Proposition 6 (`{1, 2}` tight, lonely exactly at its own real code) | PROVED + FINITE-EXACT |
| owner's F: not isomorphic as written; `F^±` exact; two copies = folded sheets; microcosm placed | PROVED + FINITE-EXACT; even-position rule unresolved |
| bounded approximants / Collatz | OPEN |
| seventeenth note's NUMEROLOGY label for `{1,2,4}` | WITHDRAWN (structural, Proposition 2) |

## 8. Directions

* **D65.** Cycle codes of the 3x+d sheets: which rational cycles have codes
  that are difference sets, or unions of cyclotomic classes forming Paley-type
  or Hadamard tournaments? Mersenne-prime lengths are the natural
  candidates.
* **D66.** The Sierpiński tournament tower (HYP-9162: doubly regular
  tournaments of every Mersenne order `2^L − 1`, `T_3` = Paley) against the
  Collatz codes modulo `2^L − 1`: is there an `L`-th analogue of
  Proposition 2 in which a cycle code is a block of `T_L`?
* **D67.** Lonely times at own codes: classify speed sets `S` that are cycles
  of some `n/2, (an+b)/2` and are lonely exactly at their real parity codes
  (`{1, 2}` is one).
* **D68.** The owner's `K`: once the even-position rule is known, test the
  isomorphism claim against the folded two-sheet relation.
