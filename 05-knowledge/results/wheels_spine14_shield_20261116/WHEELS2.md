# Zsigmondy wheels II: the complete criterion, the least periods, self-dual wheels, and Erdős's cover

Session klein-2026-11-16-wheels2, continuing THM-5254 / HYP-9823 (`two_three_wheels_synthesis_20261116/WHEELS.md`).
Notation as there: `W(N)` = primes dividing `2^N − 1`, one class modulo `ord_p(2)` each; `ν(d)` = number of primes of order
exactly `d`; `N = 2^a·m`, `m` odd; the *odd wheel* of `N` is the multiset `{oddpart(ord_p(2)) : p ∈ W(N), ord_p(2) not a
power of two}` with one class modulo each entry. Scripts `scripts/wheel_theorems_check.py`, `nu_factor.py`,
`classical_covers.py`; receipts in `receipts/`.

## 1. The complete criterion (PROVED)

**Theorem A (characterization).** `W(N)` covers `Z/N` if and only if **either** the dyadic classes alone cover `Z/2^a` —
which happens exactly when `a ≥ 6` — **or** the odd wheel covers `Z/m`.

*Proof.* By the Chinese remainder theorem `Z/N ≅ Z/2^a × Z/m`, and a class modulo `d = 2^i·m′` (`i ≤ a`, `m′ | m`) is a
product `A × B` with `A` a class modulo `2^i` in `Z/2^a` and `B` a class modulo `m′` in `Z/m` (`B = Z/m` for a dyadic
prime). Fix any cover and any `u ∈ Z/2^a` not in a dyadic `A`. The fibre `{u} × Z/m` must be covered by the `B`'s of the
non-dyadic primes whose `A` contains `u`; those form a sub-multiset of the odd wheel, and a sub-multiset cannot cover `Z/m`
if the whole odd wheel cannot (extra classes never hurt). So either every `u` is covered dyadically or the odd wheel covers
`Z/m`. Conversely, if the odd wheel covers `Z/m`, the nested dyadic classes `0 (2), 1 (4), …, 2^(a−1) − 1 (2^a)` leave one
`u`, and each non-dyadic prime's class is placed over that `u` with the prescribed `B` (THM-5254 (B)); if the dyadic
classes cover `Z/2^a`, nothing else is needed. Finally, classes modulo `2^i` are subtrees of depth `i` in the binary tree of
`Z/2^a`, so (Kraft) they can cover iff `Σ_{dyadic} 2^(−i_p) ≥ 1`; with `ν(2^i) = 1` for `i ≤ 5` and `ν(64) = 2` (THM-5254
(A)) this sum is `1 − 2^(−a) < 1` for `a ≤ 5` and `≥ 1` for `a ≥ 6`. ∎

This settles Conjecture W1 of HYP-9823 in its corrected form: "nested is without loss" is **true for `a ≤ 5`** and the
only other way to cover is the dyadic part alone, available from `a = 6` on (the composite Fermat number `F_5`). The
exact table of THM-5254 (E) agrees with Theorem A at every decided `N ≤ 150` (receipt `wheel_theorems_check_150.out`, no
mismatch).

**Theorem B (prime-power odd part).** If `m = q^e` (`q` an odd prime) and `a ≤ 5`, then `W(N)` covers `Z/N` iff
`Σ_{p non-dyadic} q^(−j_p) ≥ 1`, where `q^(j_p) = oddpart(ord_p(2))`.

*Proof.* Classes modulo `q^j` in `Z/q^e` are the depth-`j` subtrees of the `q`-ary tree of depth `e`; a multiset of depths
covers the leaves iff `Σ q^(−j) ≥ 1` (place subtrees in nondecreasing depth into the uncovered region, which is always a
union of subtrees of the current depth or deeper; necessity is measure). Apply Theorem A. ∎ (Checked against the exact
table for every prime-power odd part `≤ 150`, no mismatch.) For `e = 1` this is THM-5254 (C). Examples: `N = 36 = 4·9`
covers because `2/3 + 4/9 = 10/9 ≥ 1`; `N = 18 = 2·9` fails (`1/3 + 2/9 < 1`); `N = 108 = 4·27` covers
(`2/3 + 4/9 + 4/27 = 34/27`); `N = 54` fails (`17/27`). For `m` with two or more prime factors the odd wheel is a genuine
box-covering problem on `∏ Z/q_i^(e_i)` and the density is not sufficient (`N = 90, 132` have density above 1 and fail).

## 2. The least periods (PROVED where stated)

Let `a(q) = min{a : W(2^a q) covers Z/(2^a q)}`. By Theorem A, **`a(q) ≤ 6` for every odd prime `q`** (the dyadic part of
`W(64·q)` covers), and for `a ≤ 5` the criterion is `c(a, q) = Σ_{i≤a} ν(2^i q) ≥ q`. The cyclotomic values `Φ_(2^i q)(2)`
were factored exactly (sympy; receipts `nu_11.out`, `nu_13.out`, `nu_17.out`, `nu_19.out`, `nu_23.out`):

| `q` | `ν(q), ν(2q), ν(4q), ν(8q), ν(16q), ν(32q)` | `c(a, q)`, `a = 0..5` | `a(q)` | least period |
|---|---|---|---|---|
| 3 | 1, 0, 1, 1, 2, — | 1, 1, 2, 3 | 3 | 24 |
| 5 | 1, 1, 1, 1, 1, — | 1, 2, 3, 4, 5 | 4 | 80 |
| 7 | 1, 1, 2, 1, 2, — | 1, 2, 4, 5, 7 | 4 | 112 |
| 11 | 2, 1, 2, 2, 3, 2 | 2, 3, 5, 7, 10, 12 | 5 | 352 |
| 13 | 1, 1, 3, 2, 1, 5 | 1, 2, 5, 7, 8, 13 | 5 | 416 |
| 17 | 1, 1, 3, 2, 3, 3 | 1, 2, 5, 7, 10, 13 | 6 | 1088 (dyadic) |
| 19 | 1, 1, 3, 3, 3, 6 | 1, 2, 5, 8, 11, 17 | 6 | 1216 (dyadic) |
| 23 | 2, 1, 4, 1, 3, 3 | 2, 3, 7, 8, 11, 14 | 6 | 1472 (dyadic) |

(`Φ_352(2)` is a 49-digit number with the prime factor 275509565477848842604777623828011666349761; `Φ_544(2)` has 78
digits and three prime factors; `Φ_416(2)` five; `Φ_608(2)` has 87 digits and six prime factors, so `c(5, 19) = 17 < 19`; `Φ_736(2)` has 106 digits and exactly three prime factors — 76392570609857, a 34-digit and a 59-digit prime, all of order 736 — taken from factordb and verified locally (`nu23_verify.py`, receipt `nu_23_factordb_verify.out`; the level-5 sympy run `nu_23.out` was stopped at level 4), so `c(5, 23) = 14 < 23`.) So the odd primes whose wheel covers before the dyadic part takes over
(`a ≤ 5`) are **3, 5, 7, 11, 13** among the primes tested, and `a(q) = 6` for 17, 19, 23. Conjecture
W2 of HYP-9823 (`a(q) ≪ √q`) was therefore trivially true and badly stated; the right statement is the dichotomy above,
with the conjecture that `a(q) = 6` for every prime `q ≥ 17` (`c(5, q) < q`; the primitive parts `Φ_(2^i q)(2)` rarely
have more than five prime factors).

**Reading of 24, 36, 64 against 24, 80, 112.** The three classical covering sets are the three cheapest wheels in a
different sense each: odd part 3 at `a = 3` (24: `c(3, 3) = 3`), odd part 9 at `a = 2` (36: Kraft `10/9`), odd part 1 at
`a = 6` (64: the composite `F_5`). The least periods per prime odd part are 24, 80, 112, 352, 416, then the dyadic
1088 = 64·17 and `64·q` for the larger primes. Sierpiński's 1960 set `W(64)` is thus not an accident of history but the
first point where the dyadic part alone suffices; Riesel's `W(24)` is the first point where any odd part suffices.

## 3. Self-dual wheels and Brier numbers (PROVED)

By THM-5254 (D) a prime serves both signs of one coefficient iff its order is even. Call an even-order prime a *flipper*
if `ord/2` is odd and a *keeper* if `ord/2` is even; the prime 3 (`ord 2`) is the flipper whose two classes are the two
parities. A wheel is *self-dual* if residues `r_p` for its even-order primes make both sheets `{r_p}` and
`{r_p + ord_p/2}` cover `Z/N`.

**Theorem C (`N = 2m`).** For `N = 2m` with `m` odd, `W(N)` is self-dual iff the flipper wheel `{m′ : some prime has
order exactly 2m′, m′ | m, m′ > 1}` covers `Z/m` by one class modulo `m′` per prime.
*Proof.* Without loss 3 covers the even residues on the minus sheet, hence the odd residues on the plus sheet. A flipper
of order `2m′` placed on an odd class `r` covers, on the minus sheet, the odd `x ≡ r (mod 2m′)`, i.e. in the coordinate
`x = 2t + 1` the class `t ≡ (r − 1)/2 (mod m′)`; on the plus sheet its antipode `r + m′` is even and covers `x = 2t` with
`t ≡ (r + m′)/2 ≡ (r − 1)/2 + (m′ + 1)/2 (mod m′)`. The shift `(m′ + 1)/2 ≡ 1/2 (mod m′)` is the restriction of the single
element `(m + 1)/2 ∈ Z/m` to every `m′`, so the plus-sheet classes are a common translate of the minus-sheet classes, and
a translate of a cover is a cover. Odd-order primes serve one sign only (THM-5254 (D)). ∎
Since no prime has order 6, the modulus 3 is absent from every flipper wheel: self-duality must be bought with the
moduli 5, 7, 9, 11, 13, 15, … (orders 10, 14, 18, 22, 26, 30, …).

**Theorem D (parity-sheet bound; Brier bound).** Let `S` be a covering set of period `N` containing 3 that serves both
signs of one coefficient (odd-order primes serving one sheet each). Then

`Σ_{keepers} 1/d_p + 2·Σ_{flippers, d_p > 2} 1/d_p + (1/2)·Σ_{odd order} 1/d_p ≥ 1`

(if `3 ∉ S` the left side must be `≥ 2`). *Proof.* With 3 placed as above the minus sheet still needs its `N/2` odd
residues and the plus sheet its `N/2` even ones. A keeper's class and antipode have the same parity, so it helps at most
one sheet, by `N/d`; a flipper on an odd class helps both by `N/d`; an odd-order class has `N/d` elements of mixed parity
and helps its one sheet by `N/(2d)`. ∎ Values: `W(24)`: `0.667`; `W(36)`: `0.722`; `W(64)`: `0.5` — none can host a Brier
number (the CLOCKS.md remark on 7 generalized); the period-720 Brier set of Filaseta–Finch–Kozek scores `1.476`
(receipt `classical_covers.out`). For the self-dual question with `4 | N` the bound reads `Σ_keep 1/d + 2Σ_flip 1/d ≥ 1`;
it is **below 1 for every `N ≤ 150` with `4 | N`**, and for `N = 2m ≤ 150` the flipper wheel never reaches density 1
(receipt `wheel_theorems_check_150.out`). Hence **no wheel of period `≤ 150` is self-dual**: pre-registered P-W2 of the
previous packet is PROVED outright, including the periods 96, 120, 140, 144 where the exact two-sheet search had timed
out. Whether any wheel is self-dual (equivalently, whether some flipper wheel covers `Z/m`) remains OPEN; it needs an
odd `m` with many small divisors other than 3 and the factorisations of `2^(2m′) − 1` for `m′ | m`.

## 4. Erdős's covering system is the Riesel wheel (verified)

Erdős's 1950 system for the odd integers not of the form `p + 2^k` — `0 (mod 2), 0 (mod 3), 1 (mod 4), 3 (mod 8),
7 (mod 12), 23 (mod 24)` with the primes 3, 7, 5, 17, 13, 241 — covers `Z/24` and is, class for class, the first of the
two covering patterns of `W(24)` recorded in CLOCKS.md (phases `(0, 1, 0, 7, 3, 23)` in the order of the moduli
`2, 4, 3, 12, 8, 24`); receipt `classical_covers.out`. The repository's "covering fear" reflection (klein-20261102) cited
this system by its moduli; it is now identified with the Riesel wheel and explained by Theorem A: Erdős needed period 24
for the same reason Riesel did — the missing level 6.
