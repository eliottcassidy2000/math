# Collatz among one-out-edge systems: the power maps `n^k + 1` are bare rays (the multiplicand side), the affine family is the only one that can balance (the summand side), Collatz on `Q` is all cycles while polynomial dynamics on `Q` is all trees (`-7/4` and `c = -29/16`), the shape atlas of countable functional graphs, the shift and cellular-automaton reading, the Eckmann–Hilton defect is the carry cocycle, and the glued real line of Chamberland's extension

**Session:** opus, `collatz-posets-zeta5-20260927` (S15, fifteenth note), 2026-09-30.
**Owner's directive:** "consider instead of collatz n^3+1 for odds, and other
dynamical systems related ideas. Roughly, N^3+1 collatz feels roughly
analogous to prior work on the multiplicand graph while normal collatz
matches up with the summand graph. Consider collatz like structures on the
rationals in addition to just the natural numbers. think of possible
families of infinite directed graphs each with one outwards edge per node,
what shapes those patterns can take, maybe they relate to cellular automata.
Consider which structures emit cycles while which are a single infinite tree
and see how that relates at a deep creative level with past work regarding
-7/4 and -29/14 etc and possible structures in those types of systems.
understand that collatz is special to the fact that integer multiplication,
by eckmann-hilton, is commutative so is really made up of two identical
copies of itself, and this relates to prior ideas regarding collatz and
gluing the positive and negative integers together to make the real number
line."

**Status: PROVED (elementary) for Theorem P, Proposition S, Proposition E
and the shape lemma; FINITE-EXACT for every census (power maps to `10^5`,
rational preperiodic sets to height 60, the rational Collatz box, the
six-system census, the real fixed points); CITED for Northcott, Poonen's
conjecture with Morton and Flynn–Poonen–Schaefer (the
[Zsigmondy triad note](collatz_mod6_20260917_zsigmondy_triad.md) already
cites them), Bernstein–Lagarias, Moore–Myhill, Conway, Chamberland (from
memory), and the repo notes named by path; ANALOGY typed; Collatz OPEN.
Independent audit OWED (weekly subagent limit); self-checks: the cubic
identity against direct factorisation, the preperiodic sets by two height
bounds, the box census by two step caps.** Script
`04-computation/experiments/collatz_one_out_edge_systems_20260930.py`
(outputs `.out`, `_p3to6.out`).

**What the repo had.** The `-7/4` of the directive is the rational
three-cycle `-7/4 → 5/4 → -1/4` of `x ↦ x^2 - 29/16` in the
[arithmetic-braids geometry note](arithmetic_braids_20260917_geometry.md)
(so the second number is `-29/16`), together with the seventh-root bridge:
the doubling map `t ↦ t^2 - 2` permutes `2cos(2πk/7)` in a three-cycle whose
affine time reversal is the parabolic cubic three-cycle of `x^2 - 7/4`. The
gluing of the two integer sheets is the
[glued number lines note](glued_xor_20260921_synthesis.md) (signed
Collatz conjugacy with its valuation shell, the clock dichotomy between the
plus and minus sheets) and THM-4472 (the two diamonds are time reversal, not
the two sheets). Chamberland's real extension appears in the
[sixth note](collatz_generic_price_topology_20260927.md) §5 (not odd, so no
Borsuk–Ulam) and in the synthesis ("every continuous extension has all
periods on both sides; the sign decides stability"). The "summand" and
"multiplicand" graphs are the Farey summand/multiplicand bridge of the LRC
thread (HYP-2934/2935, LTI-079) and the user's summand/doubling rows of the
[arithmetic braids Collatz note](arithmetic_braids_20260917_collatz.md).
New here: the power-map theorem and its reading of the dichotomy, the
rational Collatz box against rational polynomial dynamics, the shape
lemma, the shift/automaton typing, the Eckmann–Hilton reading of the carry
cocycle, and the fixed-point symmetry of the real extension.

## 1. The multiplicand side: power maps are bare rays (PROVED + FINITE-EXACT)

**Theorem P.** For odd `n ≥ 3` and `k ≥ 2`, `oddpart(n^k + 1) > n`. For `k = 2`,
`v_2(n^2 + 1) = 1`. For `k = 3`, `oddpart(n^3 + 1) = oddpart(n + 1) · (n^2 - n + 1)`.
*Proof.* `n^2 + 1 ≡ 2 (mod 4)`. `n^3 + 1 = (n + 1)(n^2 - n + 1)` with the second
factor odd, and `n^2 - n + 1 > n` for `n ≥ 2`. For general `k`, `n^k + 1 ≤
2^{v} m` with `2^v ≤ n + 1` when `k` is odd (`v = v_2(n+1)`) and `v = 1` when `k`
is even, so `m ≥ n^{k-1} (n + 1)/(n + 1) > n`. ∎

So the functional graph of `n ↦ oddpart(n^3 + 1)` on odd `N` is the fixed
point `1` and a disjoint union of infinite rays: no cycle, no branching
(the images of odd `n < 10^5` are `50000` distinct numbers, in-degree `≤ 1`
throughout), and almost every odd number is an orphan (only `11` of the
`1000` odd `m ≤ 2000` have a preimage). This is the "multiplicand graph"
extreme: the multiplicative growth `n^{k-1}` is never matched by the
division, because the 2-adic valuation of `n^k + 1` is at most that of
`n + 1`.

**Why the affine family is the only balanced one (Proposition S).** For
`n ↦ oddpart(P(n))` with `P` a polynomial of degree `k`, the expected
division is `2^{E[v]}` with `E[v] = 2` (Terras's geometric valuations when
the residues are balanced), a constant, while the multiplier is `n^{k-1}`
times the leading coefficient. Balance (`log` multiplier `≈ E[v] log 2`) is
possible only for `k = 1`, and there it reads `log_2 q` against `2`: `q = 3` is
below (`log_2 3 = 1.585`, contracting on average, the Collatz case), `q = 5, 7`
above (expanding). The "summand" `d` in `qn + d` never enters the drift; it
enters only as the carry, the cocycle of the twelfth note. So Collatz is the
degree-one, below-critical member of the only family in which the sum term
can matter at all — the owner's summand/multiplicand intuition, made
precise.

## 2. Collatz on `Q` against polynomial dynamics on `Q` (FINITE-EXACT; CITED)

**Polynomial side.** For `x ↦ x^2 + c` on `Q`, Northcott's theorem gives
finitely many preperiodic rationals for each `c`; Poonen's conjecture says
rational cycles have period `≤ 3` (period 4 excluded by Morton, period 5 by
Flynn–Poonen–Schaefer, period 6 by Stoll under BSD); every other rational
wanders to infinity. Computed (height `≤ 60`): `c = -29/16` has exactly the
eight preperiodic rationals `±1/4, ±3/4, ±5/4, ±7/4` (the three-cycle, its
negatives, and `±3/4 ↦ -5/4`); `c = -1`: `{-1, 0, 1}`; `c = -2`: `{0, ±1, ±2}`;
`c = 1/4`: `{±1/2}`; `c = -3/4, -7/4`: `{±1/2, ±3/2}`; `c = -6`: `{±2, ±3}`. The
functional graph of `x^2 + c` on `Q` is finitely many unicyclic components
and infinitely many tree components, each with one ray to infinity.

**Collatz side.** On `Z_(2)` (rationals with odd denominator; the map
`x ↦ (3x+1)/2^{v}` on odd-numerator points, `-1/3 ↦ 0`), every one of the
`634` odd-numerator points `m/d` with `|m| ≤ 40`, `d ≤ 40` odd is eventually
periodic within 600 steps; the cycles reached, by denominator, are
`d = 1: -1, -7, -91, 1` (the four integer cycles with their least elements),
`5: 1/5, 19/5, 23/5`, `7: 5/7`, `11: -29/11, 1/11, 13/11`, `13: 1/13, 131/13`,
`17: -179/17, 1/17, 23/17`, `19: 5/19`, `23: 41/23, 5/23, 7/23`, `25: 17/25, 7/25`,
`29: 1/29, 11/29`, `31: 13/31`, `35: 13/35, 17/35`; no cycle has `3 | d`
(clocks are coprime to 3; preperiodic points may have `3 | d`). These are the
`3x + d` cycles of the fourteenth note's Theorem F; the periodicity conjecture
(Lagarias) says the box is typical: **Collatz on `Q` is all cycles, polynomial
dynamics on `Q` is all trees.** The opposition is a height statement: a
polynomial of degree `≥ 2` increases height, the Collatz map on `Z_(2)` is
height-recurrent on every rational (conjecturally), because its multiplier
`3/2^v` has average `log` below zero on the plus sheet and the denominator is
preserved.

**The `-7/4` bridge, placed.** The parabolic three-cycle of `x^2 - 7/4` comes
from the seventh roots of unity under doubling, the rational three-cycle of
`x^2 - 29/16` from its affine time reversal; both are the *finite* cycle part
of a map whose rational graph is otherwise trees. The Collatz analogue of
"which `c` emit cycles" is Theorem F: every denominator `d` coprime to 6
emits cycles (all `d ≤ 100` checked), i.e. every sheet `3x + d` has a cycle,
the opposite regime.

## 3. The shape atlas of countable functional graphs (PROVED, elementary)

**Shape lemma.** In a functional graph on a countable set (one out-edge per
vertex) every connected component is either *unicyclic* (one cycle with
in-trees attached) or *acyclic with exactly one forward ray* (an infinite
tree in which every vertex's forward orbit is injective and all forward
rays merge). *Proof.* A component with a cycle has exactly one (two cycles
would need a vertex with two out-edges); a component without a cycle has
injective forward orbits, so every vertex starts an infinite ray, and two
rays from the same component merge (the component is connected through
forward edges only). ∎

| system | cycles | tree (ray) components | in-degree profile (odd `m < 2·10^4`, preimages `< 8·10^4`) |
|---|---|---|---|
| `3x+1` on `N` | `{1}` (conjecturally the only one) | none (conjecturally) | orphans exactly `1/3` (the multiples of 3), in-degree 2: `40%`, 3: `14%`, 4: `3.5%` |
| `3x-1` on `N` | `1, 5, 17` | none (conjecturally) | the same profile (the sheets are conjugate) |
| `5x+1` on `N` | `1, 13, 17` | `96%` of starts escape the bound | orphans `20%`, in-degree 1: `57%`, 2: `22%` |
| `7x+1` on `N` | `1` | `99.8%` escape | orphans `57%` |
| `x^2+1`, `x^3+1` on odd `N` | `{1}` | everything else, bare rays | orphans `99%`, in-degree `≤ 1` |
| `3x+1` on `Z_(2)` | infinitely many (one family per `d`) | none (conjecturally) | — |
| `x^2 + c` on `Q` | finitely many, period `≤ 3` | infinitely many | `≤ 2` (`±√(x - c)`) |

The drift decides the regime (the tenth note's typology): below critical,
everything is unicyclic; above, rays dominate and the in-degree profile
thins toward orphans; polynomial growth is the limit with bare rays.

## 4. The shift and the automaton reading (CITED + ANALOGY)

On `Z_2` the map `T` is conjugate to the one-sided full 2-shift
(Bernstein–Lagarias; the parity-vector bijection of Terras), hence
surjective, two-to-one, and Haar-measure-preserving: the balance that
surjective cellular automata have by the Moore–Myhill (Garden-of-Eden)
theorem. It is not a cellular automaton: no sliding-block rule in the
2-adic digits computes `3n + 1` (the carry propagates; the cocycle of the
twelfth note is the non-locality). `N ⊂ Z_2` is a countable forward-invariant
set, the Garden-of-Eden states of the Syracuse graph on `N` are exactly the
odd multiples of 3 (density `1/3`, the orphans of §3), and Conway's theorem
that the generalized Collatz family is Turing-complete is the counterpart
of the undecidability of long-term behaviour for automata (Rule 110). The
two commuting "shifts" are multiplication by 3 (a 3-adic shift) and
division by 2 (a 2-adic shift); their parity-controlled alternation is
Collatz, and the pair `×2, ×3` is the setting of Furstenberg's measure
rigidity, a known heuristic bridge that no one has crossed (ANALOGY).

## 5. Eckmann–Hilton and the two sheets (PROVED, elementary)

**Proposition E.** The multiplier monoid `<2, 3> ⊂ Q^×` is `N × N`, two
commuting copies of `N`; the interchange law of Eckmann–Hilton holds for
multipliers and fails for the affine letters by exactly the carry: the
words `(1, 2)` and `(2, 1)` have carries `5` and `7`, so the two composites
`x ↦ (9x + 5)/8` and `x ↦ (9x + 7)/8` differ by the translation `-1/4`, and in
general the commutation defect is `β(u, v) = D_u D_v (x_v - x_u)` (twelfth
note). The defect is antisymmetric, so the structure is *symmetric*, which
is what Eckmann–Hilton would force if the interchange held; the obstruction
is the cocycle, and Terras's bijection says the affine letters generate a
free monoid (no relation at all). *Proof.* Twelfth note, Propositions 3–4. ∎

**The two copies.** On `Z_(2)` negation conjugates `3x + 1` to `3x - 1`, and
the fixed point of a word for `3x - 1` is minus the fixed point for `3x + 1`
(checked on the cycle words): the plus and minus sheets are two copies of
one map glued at the fixed point `0`, and the fourteenth note's Theorem F
glues all sheets `3x + d` into one tree through the mediant law. The glued
number lines note's clock dichotomy is the arithmetic content: the clocks
behave differently on the two sheets (the `-17` cycle has a prime clock, the
plus sheet none beyond the unit clock of `{1}`).

**The real line.** Chamberland's extension
`f(x) = (x/2) cos^2(πx/2) + ((3x+1)/2) sin^2(πx/2)` carries both sheets on one
line (`f(1) = 2, f(2) = 1, f(-5) = -7, f(-7) = -10, f(-10) = -5`). Its fixed
points in `[-6, 6]` are `-5.468, -4.540, -3.446, -2.577, -1.278, 0.278,
1.577, 2.446, 3.540, 4.468, 5.526`, and the fixed-point equation
`cos^2(πx/2) = (x + 1)/(2x + 1)` is invariant under `x ↦ -1 - x`, which swaps
the integer fixed points `0` and `-1`: the real extension glues the sheets
about `-1/2`, not about `0`. The only attracting real fixed point is
`0.278` (multiplier `0.386`); its mirror `-1.278` is repelling (`1.614`). The
synthesis's record stands: every continuous extension has all periods on
both sides (Sharkovskii gives nothing), and the sign decides stability.

## 6. Verdicts

| claim | status |
|---|---|
| Theorem P (power maps increase; cubic identity); rays with in-degree `≤ 1` to `10^5`; orphan census | PROVED + FINITE-EXACT |
| Proposition S (only degree one can balance) | PROVED (heuristic average, exact for the drift sign) |
| rational preperiodic sets of `x^2 + c` (eight points at `-29/16`) | FINITE-EXACT (Northcott, Poonen CITED) |
| rational Collatz box: all 634 points preperiodic; cycles by denominator | FINITE-EXACT |
| shape lemma; the six-system census | PROVED + FINITE-EXACT |
| shift/automaton typing | CITED + ANALOGY |
| Proposition E; sheet symmetry `x_w(3x-1) = -x_w(3x+1)` | PROVED |
| Chamberland fixed points; the `x ↦ -1 - x` symmetry of the fixed-point equation | FINITE-EXACT + PROVED |
| Collatz | OPEN |

## 7. Directions

* **D49.** The balanced family beyond `qn + d`: maps `n ↦ (qn + d)/p^{v_p}`
  for other primes `p` (Matthews–Watts) and mixed divisions; the critical
  surface `log q = (p/(p-1)) log p` and which integer points lie below it.
* **D50.** The rational box to height `200` with the cycle inventory per
  denominator, against Theorem F's prediction that every `d` coprime to 6
  emits a cycle; the first rational, if any, whose orbit is not eventually
  periodic would refute the periodicity conjecture.
* **D51.** The exact symmetry behind the `x ↦ -1 - x` invariance of
  Chamberland's fixed-point equation, and whether the attracting real fixed
  point `0.278` has a basin that meets `Z` (it cannot, since integers map to
  integers), i.e. how the integer cycles sit among the real basins.
* **D52.** A mixed-radix digit system in which `T` is a sliding-block code,
  or a proof that none exists (the carry's non-locality quantified by the
  cocycle's growth).
