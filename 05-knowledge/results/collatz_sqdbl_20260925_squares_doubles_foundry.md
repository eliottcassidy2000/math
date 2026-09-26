# The squares/doubles transport foundry: an additive-to-multiplicative axis for Collatz reframes, its 128-cell deck, six exact lemmas, and the one two-place statement it leaves open

**Status: METHOD / TYPOLOGY (the deck and its typing are modelling
judgements), plus PROVED elementary lemmas (Lemma 1, Theorem 2, Lemmas 3-5,
Proposition 6, proofs below) and FINITE-EXACT probes with printed universes
and controls. Collatz, the 3n-1 sheet and the Periodicity Conjecture remain
OPEN. No new proof mechanism is claimed; the session's honest product is one
exact two-place reformulation of the divergence half and one new hypothesis
(HYP-9160). Session `collatz-squares-doubles-20260925` (opus), 2026-09-25.**

Scripts: `04-computation/experiments/collatz_sqdbl_20260925_foundry.py`
(the deck) and `collatz_sqdbl_20260925_probes.py` (probes P1-P6).
Outputs: `collatz_sqdbl_20260925_foundry.out`, `collatz_sqdbl_20260925_probes.out`.

## 0. The seed, and the inheritance pass

The owner's seed: *multiplication relates to the squares the way addition
relates to the doubles*. Read literally it is a statement about diagonals:
the doubles `x+x` are the diagonal of addition, the squares `x*x` the
diagonal of multiplication, and the exponential carries one pair to the
other. Applied to the shortcut map `T(n) = n/2, (3n+1)/2`, the transport
`n -> X = q^n` sends

| additive world | multiplicative world |
|---|---|
| `n` is a double (even) | `X` is a square |
| halving `n -> n/2` | square root `X -> sqrt X` |
| `3n+1` | cube and multiply, `X -> q X^3` |
| `3n-1` (the other sheet) | `X -> q^-1 X^3` |
| the clock `3^L / 2^K` | the exponent word `3/2^v` |

The session used this as a *transport axis* and asked, procedurally, what
each Collatz ingredient becomes under sixteen transports of this kind
(exponential, discrete-log, monoid extension, polarization, quarter squares,
difference of squares, quadratic characters, Jacobi symbols, local squares,
global squares, Hilbert symbols, the real value of the digit series, cubic
theta polarization, geometric-mean fairness, elliptic heights, sums of
squares).

**Inheritance pass** (AGENTS.md): the closest proved mechanisms are the
strategy cube THM-4474 and the price theorem THM-4475, the word-function
theorem (wave six: any invariant that is a function of the parity word is
sheet-blind), the Eliahou product identity on cycles (S4 of the
[counterexample portrait](collatz_mod6_20260922_counterexample_portrait.md))
and Proposition T of the [hard-class lane](collatz_procgen_20260922_hard_class.md)
(the 2-adic and real values of the same Bernstein series are not both
rational, for bounded-discrepancy words with `Dio > 2`). The canonical
hostiles are the `3n-1` cycles `{5,7,10}` and `{17,...}` (SHEET), `5n+1`
(DRIFT) and planted flips (DEFECT, THM-4470/4475). The corrected near miss
is the procgen thread's own 2026-09-22/23 corrections (choice is blind to the
drift; the `32/27` escape floor). The least-used sidecar is the *real* value
of the Bernstein series, which the atlas uses only for supercritical words.

## 1. The grammar and the deck

A cell is `TRANSPORT x INGREDIENT` with

```text
TRANSPORT   in {EXP, DLOG, MONOID, POLAR, QSQ, DIFFSQ, QCHAR, JACOBI, LOCSQ,
                GLOBSQ, HILB, REALVAL, THETA3, GMFAIR, HEIGHT, SUMSQ}
INGREDIENT  in {parity gate, halving, 3n+1, valuation v_2(3n+1), cycle
                equation/clock, parity word/Bernstein series, sign law, drift}
```

Each transport is typed by its kind (lossless, residue, sparse, sign
transfer, two-place, ...), the ambient structure it gains, and whether it
can see SHEET, DRIFT and DEFECT. Lossless transports inherit every blindness
of the statement they transport; residue transports are word functions and
sheet-blind by the word-function theorem; density-zero objects carry no
measure statement. The 128 cells classify as

| class | cells | meaning |
|---|---|---|
| EMPTY | 81 | the transport does not touch the ingredient |
| COSMETIC | 9 | exact restatement, no new predicate |
| RESIDUE LAW | 8 | exact word-function lemmas (Lemmas 3-5b) |
| NEW OBJECT | 5 | the multiplicative Collatz map (Theorem 2) |
| TWO-PLACE | 5 | real value vs 2-adic value of the digit series (Proposition 6, HYP-9160) |
| EXACT MODEL | 4 | the parity graph as a square-root graph on `F_p^*` (Lemma 1) |
| SPARSE, SPARSE/SIGN | 6 | squares and sums of squares in the orbit (Lemma 5a) |
| FAMILY VIEW, SIGN TRANSFER, DEGENERATE, AMBIENT ONLY | 7 | typed and parked |

The 47 curated cards, each with statement, status, class and cheapest
decisive test, are printed by the deck script (`.out`, section "curated
cards"). Only the TWO-PLACE class sees the sheet, the drift and the defect
at once. That is the synthesis's stated requirement for a mechanism
([section 4](collatz_procgen_20260922_synthesis.md)), so the deck's ranking
is forced: REALVAL first, then the new object, then the exact model, then
the residue laws. Everything below that is recorded so it is not regenerated.

## 2. Exact results

### Lemma 1 (the exponential transport and its Fermat-prime model; PROVED, FINITE-EXACT)

(a) Let `q` have infinite order in an abelian group written multiplicatively.
Under `n -> X = q^n`: `n` is even iff `X` is a square in `<q>`; `n/2`
corresponds to the square root of `X` inside `<q>`; `3n+1` corresponds to
`q X^3` and `3n-1` to `q^-1 X^3`. A cycle with `L` odd and `K` even steps
is `X^(3^L/2^K - 1) = q^(-c)`, i.e. `X = q^(c 2^K/(2^K - 3^L))`.

(b) Let `p = 2^k + 1` be prime and `g` a primitive root. The parity graph
`G_+` of THM-4474 (nodes `Z/2^k`, edges from `s` to the two lifts of
`T(s) mod 2^(k-1)`) is carried by `s -> g^s` to the graph on `F_p^*` with
edges `y -> +-sqrt(y)` when `y` is a square and `y -> +-sqrt(g y^3)` when it
is not. The `3n-1` graph uses `g^-1`; negation `s -> -s` is `y -> y^-1`.

*Proof.* (a) is the definition of the transport. (b) `g^(2^(k-1)) = -1`, so
the two lifts `t, t + 2^(k-1)` of `T(s) mod 2^(k-1)` are carried to
`+-g^(T(s))`. For `s` even these are the two square roots of `g^s = y`. For
`s` odd, `g^((3s+1)/2)` squares to `g^(3s+1) = g y^3`. ∎

P2 checks the edge sets for `p = 17, 257, 65537` (`g = 3`, 32, 512 and
131,072 edges): identical. The loops at `y = 1` (node `0`, density `0`) and
at `y = g^-1` (node `-1`, density `1`) are the two loops that give Collatz
its maximal window `[0,1]` in THM-4474. A non-primitive `g` (2 mod 17) fails
the construction, as it must. **Verdict.** An exact picture of THM-4474:
"provable iff every square-root cycle has non-square density `< log_3 2`".
The field addition of `F_p` is new structure that no Collatz predicate uses
(card DLOG x DRIFT); nothing is gained beyond the picture.

### Theorem 2 (the multiplicative Collatz map; PROVED, FINITE-EXACT to 10^5 on both sheets)

Define on integers `X >= 1`

```text
M(X) = sqrt(X)          if X is a perfect square,
       rad(X) * X^3     otherwise,
```

with `rad(X)` the product of the distinct primes dividing `X`. On the
exponent vector `e = (e_p)` of `X` (support invariant) this is: if all `e_p`
are even, `e -> e/2`; otherwise every nonzero `e_p -> 3 e_p + 1`. So `M`
restricted to the powers of one squarefree `m` is the Collatz map `C`
(`x/2` or `3x+1`) on the exponent, and the seed's dictionary is literal:
the square gate is the parity gate taken simultaneously in every coordinate.

**Theorem.** (i) If the exponents of `X` are not all equal, the `M`-orbit of
`X` diverges: after finitely many steps the gate never opens again and
every exponent follows `x -> 3x+1` forever. (ii) If all exponents equal `e`,
the orbit reaches `rad(X)` iff the Collatz orbit of `e` reaches `1`. Hence
Collatz is equivalent to "every perfect power of a squarefree number
returns to its radical under `M`", and `M` is provably divergent off that
diagonal. The same holds for the sheet `M_-(X) = X^3/rad(X)` (exponents
`3e - 1`), with the `3n-1` cycles appearing as `m^5 -> m^14 -> m^7 -> m^20 -> m^10 -> m^5`.

*Proof.* (ii) is immediate. (i) Take two coordinates `e_p != e_q`. If at
some step their parities differ, then the gate is closed, both receive
`3x+1`, which flips parity, so the parities differ again; the gate can open
only when all coordinates are even, hence never again, and the coordinates
grow without bound. If the parities never differ, the gated dynamics
coincides on both coordinates with `C`, and the two `C`-orbits have the same
infinite parity sequence. The parity vector determines the 2-adic integer
(Terras; Lagarias 1985, Theorem B), so `e_p = e_q`, a contradiction. ∎

P1 runs the exponent dynamics for all `X <= 10^5`: 38,949 off-diagonal `X`
diverge and 61,050 diagonal `X` return to the radical (plus sheet); on the
minus sheet 61,035 return, 15 enter the exponent cycles, and again all
38,949 off-diagonal `X` diverge. Zero violations. **Verdict.** A lossless
reframe with a proved companion theorem, and a control: the "square gate"
is a synchronisation constraint. It sees nothing the parity gate does not.

### Lemma 3 (quadratic characters read the first two valuation bits; PROVED, FINITE-EXACT n < 2^20)

For odd `n` and `s = +-1`: `v_2(3n + s) >= 2` iff `s = chi_{-4}(n)`, and
`v_2(3n + s) >= 3` iff moreover `chi_8(n) = -1`.

*Proof.* `3n + s = 0 mod 4` iff `s = -3n = n mod 4`, i.e. `s = chi_{-4}(n)`.
`3n + s = 0 mod 8` with `s = +-1` needs `3n = -+1 mod 8`, i.e. `n = 5`
(`s = +1`) or `n = 3` (`s = -1`) mod 8: exactly `chi_8(n) = -1` together
with `s = chi_{-4}(n)`. ∎

Consequences. The unique class-(i) strategy of THM-4474 at level 2 is the
character `chi_{-4}` (every odd `n` is a two-step down point), the unique
class-(ii) strategy is `-chi_{-4}` (every odd `n` is an up point with
`v = 1`, so the odd residues form a closed class of density `1`), and
Collatz and `3n-1` are the two constant characters. P3 confirms the four
behaviours on odd `n < 10^5`. Since `(Z/2^k)^*` has only four real
characters, bits `>= 3` of the valuation are not quadratic-character data.

### Lemma 4 (the Jacobi symbol along a Syracuse orbit; PROVED, FINITE-EXACT)

Let `m_(i+1) = (3 m_i + b)/2^(v_(i+1))`, `b = +-1`, all `m_i` positive odd
and prime to 3. Then on both sheets

```text
(3 / m_i) = (-1)^( v_i + [v_(i+1) = 1] ).
```

For `5n +- 1` the law is `(5 / m_i) = (-1)^(v_i)` with no coupling term.

*Proof.* `m_i = b 2^(-v_i) mod 3`, so `(m_i/3) = b (-1)^(v_i)`. Jacobi
reciprocity gives `(3/m_i) = (m_i/3)(-1)^((m_i-1)/2)`. For `b = +1`,
`m_i = 3 mod 4` iff `v_(i+1) = 1`; for `b = -1`, `m_i = 3 mod 4` iff
`v_(i+1) >= 2`, and the extra sign `b = -1` restores the same formula.
For `q = 5`, `(m_i/5) = (-1)^(v_i)` and `5 = 1 mod 4` removes the
reciprocity sign. ∎

P4: 1,100,000 checks on each of the four maps, zero failures; the
law with the coupling term moved fails about half the time (control).
**Verdict.** The reciprocity sign, which is where the real place enters
quadratic reciprocity, couples adjacent valuations; but for positive `m`
the law is a function of the word on both sheets, so it is sheet-blind, as
the word-function theorem predicts. It is an exact `mod 12` law, not a
mechanism.

### Lemma 5 (squares in the orbit; PROVED, FINITE-EXACT s <= 2*10^6)

(a) For `s = 2t+1`, the odd Syracuse successor of `s^2` on the plus sheet is
`3t(t+1) + 1 = 6 T_t + 1` (`T_t` triangular), with `v = 2` exactly. So every
odd square is a down point of `T_+` and, since `3s^2 - 1 = 2 mod 8`, an up
point of `T_-`.

(b) The successor of an odd square is again a square iff `2t+1 = y_j` with
`j` odd, where `x_j + y_j sqrt3 = (2 + sqrt3)^j`: `s = 1, 15, 209, 2911,
40545, 564719, ...` (`225 -> 169`, `43681 -> 32761`). No chain of three
squares exists except `1 -> 1`.

(c) An odd iterate `m` is a 2-adic square (`m = 1 mod 8`) iff its next
valuation is exactly `2`; it is a 3-adic square unit (`m = 1 mod 3`) iff its
previous valuation is even. A global square in the orbit forces both.

*Proof.* (a) `3s^2 + 1 = 4(3t^2 + 3t + 1)` and `3t(t+1) + 1` is odd.
(b) `3t^2 + 3t + 1 = u^2` iff `(2u)^2 - 3(2t+1)^2 = 1`; the solutions of
`x^2 - 3y^2 = 1` have `x_j` even iff `j` is odd (the parities of `(x_j, y_j)`
cycle `(0,1), (1,0)`), and then `y_j` is odd. A third square needs
`x_j = 2 y_i` for odd `i, j`, i.e. `sqrt3 (e^j + e^-j) = 2 (e^i - e^-i)` with
`e = 2 + sqrt3`. If `j >= i+1` the left side exceeds `sqrt3 e^(i+1) > 2 e^i`;
if `j <= i-1` the left side is at most `sqrt3 e^(i-1)(1 + e^-2) < 2 e^i (1 - e^-2i)`
for `i >= 1`; so `i = j`, which forces `y_i = 1`. (c) is Lemma 3 with
`s = +1`, and `m_i = 2^(-v_i) mod 3`. ∎

P5 finds exactly the six Pell values below `2*10^6` and no length-3 chain.
**Verdict.** Exact and sign-aware (squares are positive), but squares are a
density-zero set with no closed dynamics; the Pell bridge is a curiosity for
the Pell/Pythagorean threads (THM-3357 family), not a Collatz tool.

### Proposition 6 (the real value of the Bernstein series, orbitwise; PROVED, FINITE-EXACT)

Let `q >= 3` be odd, `b = +-1`, `n` a positive odd integer whose Syracuse
orbit `m_(l+1) = (q m_l + b)/2^(v_(l+1))` stays positive, `d_L = v_1 + ... + v_L`,
`d_0 = 0`, and `B_L = sum_(l<L) q^(L-1-l) 2^(d_l)`. Then for every `L`:

1. `m_L 2^(d_L) = q^L n + b B_L`;
2. `m_L 2^(d_L) / q^L = n prod_(l<L) (1 + b/(q m_l))`;
3. hence the partial real sums of the digit series satisfy
   `R_L(d) := sum_(l<L) 2^(d_l)/q^(l+1) = B_L/q^L = b n ( prod_(l<L)(1 + b/(q m_l)) - 1 )`.

Consequences, with `R(d) = lim R_L(d) in (0, infinity]` and the growth
constant `c := lim m_L 2^(d_L)/q^L`:

4. **Minus sheet (`b = -1`).** `R_L(d) <= n` for every `L`, so `R(d)` converges
   and `R(d) <= n` for **every** positive orbit; `c = n - R(d) = n prod_l (1 - 1/(q m_l)) >= 0`,
   and `c = 0` iff `sum_l 1/m_l = infinity`, in particular for every
   eventually periodic orbit.
5. **Plus sheet (`b = +1`).** `R(d) < infinity` iff `sum_l 1/m_l < infinity`;
   then `c = n + R(d) > n` and the orbit diverges with `m_L ~ c q^L/2^(d_L)`.
   Convergent orbits have `R(d) = infinity`.
6. **2-adically** the same series has value `R_2(d) = sum_l 2^(d_l) q^(-l-1) in Z_2`
   and `n = -b R_2(d)` (Bernstein's formula).

*Proof.* (1) by induction, `m_(L+1) 2^(d_(L+1)) = q m_L 2^(d_L) + b 2^(d_L)`.
(2) multiply `m_(l+1) 2^(v_(l+1)) = q m_l (1 + b/(q m_l))` over `l < L`.
(3) divide (1) by `q^L` and compare with (2). (4) the factors lie in `(0,1)`,
the product decreases, and a product of factors `1 - x_l` with `0 < x_l <= 1/3`
tends to `0` iff `sum x_l` diverges. (5) similarly with `1 + x_l`. (6) from
(1), `q^L n = -b B_L mod 2^(d_L)` and `d_L -> infinity`. ∎

P6 verifies (3) with exact rationals along every orbit of odd `n <= 2000` on
both sheets, `R_L <= n` on the minus sheet, and the closed-form value
`R(d) = n` for all 1,000 orbits that enter a `3n-1` cycle (the Eliahou
identity is the periodic case). On `5n +- 1` the presumed divergent orbits of
`n = 37, 47` (plus) and `n = 9, 11` (minus) have `c_L` stabilised to twelve
digits after 100 odd steps (`37.3743738502`, `47.8391985283`, `8.5283668520`,
`10.6604585650`) with `sum 1/m_l` about `0.05`-`0.27`: full-rate divergence.

**What is new and what is not.** Item (2) at a cycle is Eliahou's identity
(S4 of the counterexample portrait); Proposition T of the hard-class lane
compares the two places for bounded-discrepancy words. New here are the
orbitwise identity (3) for *every* word, the inequality `R(d) <= n` for
*every* minus-sheet orbit, and the dichotomy `c = 0 iff sum 1/m_l = infinity`.

**Reading.** The sign law of the divergence half is one real inequality. A
positive plus-sheet orbit is a word whose 2-adic value is `-n` while its real
value is positive and unconstrained. A positive minus-sheet orbit (equivalently
a negative plus-sheet orbit, `T_-(n) = -T_+(-n)`) is a word whose 2-adic
value is `n` and whose real value is at most `n`, with equality exactly when
the reciprocals of the odd iterates are not summable. Since `T_-(n) = -T_+(-n)`,
item 4 also says: every negative `T_+` orbit has `R(d) <= |n|`. The drift is
visible in the same formulas through `q`. This is the only cell of the deck
that is sheet-, drift- and defect-aware, which is why it is ranked first.

### HYP-9160 (no slow divergence; the real and 2-adic values never coincide off periodic words)

**Update, later the same day: PROVED as [THM-4476](../../01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md)** (thin divergence: `N(X) = O(X^(0.95+eps))` for every non-periodic orbit; see [the thin-divergence note](collatz_thin_20260925_thin_divergent_orbits.md)). The text below is the original formulation.

**Statement.** Every divergent Syracuse orbit, on either sheet and of either
sign, has `sum_l 1/m_l < infinity`. Equivalently on the minus sheet: for a
word `d` that is not eventually periodic, the real value `R(d)` and the
2-adic value `R_2(d)` of `sum_l 2^(d_l)/3^(l+1)` are never both equal to
the same positive integer. On the plus sheet: every divergent positive orbit
has a finite real Bernstein value.

**Why this is the right size.** It isolates the slow-divergence part of the
Periodicity Conjecture, the only part where the real place gives an
*equality* (`c = 0`) rather than an inequality. It is implied by "no
divergent orbit" on the relevant sheet and by nothing weaker that is known.
It holds for words of bounded critical discrepancy (they carry no integer
orbits at all; the in-house capacity-plus-ordered-carry theorem). A slow
divergence `m_l ~ l^alpha`, `alpha <= 1`, has discrepancy of order `log l`,
just beyond that theorem. Cheapest partial: extend the bounded-discrepancy
mechanism to discrepancy `<= C log l`. No decisive finite test exists (it
needs a divergent orbit); the file is
[HYP-9160](../hypotheses/HYP-9160-no-slow-divergence-real-2adic-coincidence.md).

## 3. The rest of the deck, briefly

* **POLAR.** Syracuse is the slice `m = 3` of the commutative operation
  `n o m = (nm + 1)/2^(v_2(nm+1))`; the `qn+1` maps are the other slices, and
  the diagonal `n o n = (n^2+1)/2` has `v = 1` always, so it is monotone and
  useless (no balance). Polarization `nm = ((n+m)^2 - (n-m)^2)/4` needs
  addition and halving, the two Collatz operations, which is the seed's own
  content; the drift lives in the slice index.
* **QSQ.** `(3n+1)/2 = (Q(n+3) - Q(n-3) + 1)/2` and `n/2 = sqrt(Q(n))` with
  `Q(x) = floor(x^2/4)`, the Turán number `t(x,3)`. Cosmetic.
* **DIFFSQ.** At a doubled clock `(2a, 2b)`, `2^(2a) - 3^(2b) = (2^a - 3^b)(2^a + 3^b)`,
  so a positive cycle there needs `2^a + 3^b | B`. A side constraint with no
  leverage below the Hercher bound. For `L` even a positive denominator is
  `7 mod 8`, never a sum of three squares (Legendre); for `L` odd it is `5 mod 8`.
* **HILB.** The product formula `prod_v (2^K - 3^L, -1)_v = 1` moves the sign
  of `2^K - 3^L` to the finite places: the sign is the parity of the number
  of odd-multiplicity primes `3 mod 4` of `|2^K - 3^L|`, corrected by `L mod 2`.
  It needs the factorisation of `2^K - 3^L`, so it gives nothing over Baker.
* **THETA3.** The cubic theta `Theta_3(x,y) = sum q^(k^3) x^(k^2) y^k` satisfies
  `Theta_3(x,y) = 1 + qxy Theta_3(q^3 x, q^3 x^2 y)`: the shift acts on cubes
  through squares, the seed's polarization. HYP-9127's 2026-09-23 update
  already records this as a unipotent Mahler system and rules it out.
* **GMFAIR.** No geometric-mean-fair pairing exists for `T`
  (`(a/2)(3b+1)/2 = ab` forces `b = 1`); AM-fairness (THM-4470) becomes
  GM-fairness only inside the exponential image. Degenerate.
* **HEIGHT, SUMSQ.** The parallelogram law of the canonical height is the
  seed's polarization, and on a rank-one curve nothing beyond `n` appears.
  Odd sums of two squares are `1 mod 4`, hence down points of `T_+` and never
  down points of `T_-`; positivity is algebraic globally (Lagrange) but every
  2-adic integer is a sum of four squares, so nothing local sees it.

## 4. Was a winning reframe found?

By the synthesis's own criterion (a mechanism must see `log 3 < 2 log 2` and
the sign together, orbit by orbit), the only candidate is the two-place
formulation of Proposition 6. It is exact, sheet- and drift-aware, and it
turns the divergence half into a comparison of the real and the 2-adic value
of one digit series. It does **not** supply the missing transversality
statement; it names its slow-divergence corner (HYP-9160) and shows that the
full-rate corner is invisible to the real place. So: no winning reframe in
the sense of a proof route; one clean reformulation, one hypothesis of the
right size, one exact model of THM-4474, one new dynamical object with a
proved theorem, and four residue laws that will not need regenerating.

**Non-consequences.** Nothing here bears on Collatz, the `3n-1` sheet, the
Periodicity Conjecture, E-SCC or HYP-9127. Lemmas 3-5 are word functions.
Theorem 2 is equivalent to Collatz on its diagonal. Lemma 1 is a picture.

## 5. Frontier and next probes

1. HYP-9160 on the plus sheet with the bounded-discrepancy mechanism relaxed
   to `O(log l)`: the first regime in which `sum 1/m_l` can diverge.
2. The DLOG ambient addition (card DLOG x DRIFT): a predicate on the
   square-root graph of `F_65537^*` that uses `y + 1` would be genuinely new
   structure; none was found.
3. The doubled-clock divisibility `2^a + 3^b | B` as a filter in the
   cycle-gate census of wave 11 (cheap; expected vacuous).
4. Proposition T with the orbitwise identity: the real tail of a word is
   `n prod_(l>=L)(1 +- 1/(3 m_l)) - ...`; whether this sharpens the
   `Dio > 2` threshold is untested.

## 6. Reproduction

```bash
python3 04-computation/experiments/collatz_sqdbl_20260925_foundry.py > 05-knowledge/results/collatz_sqdbl_20260925_foundry.out
python3 04-computation/experiments/collatz_sqdbl_20260925_probes.py  > 05-knowledge/results/collatz_sqdbl_20260925_probes.out
```

Runtime about one minute (P6 dominates). SHA-256 (raw LF bytes):
`collatz_sqdbl_20260925_foundry.py` `ff434d4a880c174c00d159a551de45a4f7af94b1a6c59a957b496aff141f5915`;
`collatz_sqdbl_20260925_probes.py` `9086597813e5b1796a556f0bf2d0c7dba44e06f66f2d48f4ac254a012fab038c`;
`collatz_sqdbl_20260925_foundry.out` `5aaab54cc217b440eb860081113ee3714644686c87a6520a9a79146b4d0c1c25`;
`collatz_sqdbl_20260925_probes.out` `71f70f25a088852d5f6d81d214802bc426a8fc12be7f9acdcdfd37b2d1c55706`.
