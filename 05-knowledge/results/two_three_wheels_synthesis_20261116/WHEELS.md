# Zsigmondy wheels: which periods carry a covering set for k·2ⁿ ± 1, and what the Riesel wheel shares with 24

Session klein-2026-11-16-two3. Scripts `scripts/wheels2.py`, `wheels3.py`, `brier_sheets.py`, `exact_one.py`,
`zsig_wheels.py`, `leech_numerics.py`; receipts in `receipts/`. Everything computational is exact (bitmask exact-cover
search over residues; sympy factorisations of `2^d − 1`, `d ≤ 150`).

## 0. Definitions

Fix the base 2. For `N ≥ 1` the **wheel** `W(N)` is the set of primes dividing `2^N − 1`, i.e. the primes `p` with
`ord_p(2) | N`. For `d ≥ 2` let `ν(d)` be the number of primes with `ord_p(2) = d` (the primitive prime divisors of
`2^d − 1`). Zsigmondy–Bang: `ν(d) ≥ 1` for every `d ≥ 2` except `d = 6` (`2^6 − 1 = 63 = 3²·7`); `ν(1) = 0`.

A **cover of period `N`** is a choice of one residue class `r_p (mod ord_p(2))` for each `p ∈ W(N)` whose union is all of
`Z/N`. By the Chinese remainder theorem a cover gives a Riesel number (`k ≡ 2^(−r_p) mod p` makes `p | k·2^n − 1` for
`n ≡ r_p`) and a Sierpiński number (`k ≡ −2^(−r_p)`), and conversely every covering set of period `N` for `k·2^n ± 1` is a
subset of `W(N)` with such classes; unused primes may take any class, so "`W(N)` covers `Z/N`" is the right existence
statement. The classical sets: Riesel's `{3, 5, 7, 13, 17, 241}` for 509203 (also Sierpiński's 271129) is **all of
`W(24)`**; Selfridge's `{3, 5, 7, 13, 19, 37, 73}` for 78557 is `W(36) ∖ {109}`; Sierpiński's 1960 set
`{3, 5, 17, 257, 641, 65537, 6700417}` is **all of `W(64)`**, with `Σ 1/ord = 1/2 + 1/4 + 1/8 + 1/16 + 1/32 + 2/64 = 1`
exactly (receipt `wheels2_base2_150.out`, rows 24, 36, 64; CLOCKS.md of `logarithmic_repairs_covers_20261116` for the
two patterns at 24).

## 1. Theorems (PROVED)

**Theorem A (the dyadic levels are always present).** For every `i ≥ 1`, `ν(2^i) ≥ 1`: every prime factor of the Fermat
number `F_(i−1) = 2^(2^(i−1)) + 1` has order exactly `2^i`. Indeed `2^(2^(i−1)) ≡ −1 (mod p)` gives `2^(2^i) ≡ 1` and
`ord ∤ 2^(i−1)`. So `ν(2) = ν(4) = ν(8) = ν(16) = ν(32) = 1` (the primes 3, 5, 17, 257, 65537) and `ν(64) = 2`
(641 and 6700417, because `F_5` is composite). ∎

**Theorem B (reduction to the odd part).** Write `N = 2^a·m` with `m` odd, and call a prime of `W(N)` *dyadic* if its
order is a power of two. Let `L` be the set of residues not covered by the dyadic classes.

(i) *Necessary condition.* If `W(N)` covers `Z/N` and `L ≠ ∅` (which holds whenever `Σ_{i≤a} ν(2^i)/2^i < 1`, in
particular for every `a ≤ 5`), then

`Σ_{p ∈ W(N) non-dyadic} 1 / oddpart(ord_p(2)) ≥ 1.`

(ii) *Sufficient condition (nested-dyadic criterion).* If `Z/m` is covered by one class modulo `oddpart(ord_p(2))` for
each non-dyadic `p ∈ W(N)`, then `W(N)` covers `Z/N`.

*Proof.* Every class modulo `2^i` with `i ≤ a` is a union of cosets of `2^a·Z/N`, so `L` is such a union. A non-dyadic
prime has order `d = 2^i·m′` with `i ≤ a` and `m′ | m`, `m′ > 1`; its class meets a coset `C` of `2^a·Z/N` either not at
all (if the residues modulo `2^i` disagree) or in exactly `|C|/m′` elements (inside `C ≅ Z/m`, via `x = c + 2^a·t`, the
condition is `t ≡ s (mod m′)` because `2^a` is invertible modulo `m′`). Hence the class covers at most `|L|/m′` elements
of `L`, and a cover forces `Σ 1/m′ ≥ 1`: this is (i). For (ii) choose the dyadic classes nested — `0 (mod 2)`, `1 (mod 4)`,
`3 (mod 8)`, …, `2^(a−1) − 1 (mod 2^a)` — which leave uncovered exactly the coset `C = {x ≡ 2^a − 1 (mod 2^a)} ≅ Z/m`;
for each non-dyadic prime pick, by the Chinese remainder theorem, the class that is `≡ 2^a − 1 (mod 2^i)` and whose
restriction to `C` is the prescribed class modulo `m′`; the hypothesis says these restrictions cover `Z/m`. ∎

**Corollary C (prime odd part).** Let `N = 2^a·q` with `q` an odd prime and `a ≤ 5`. Then `W(N)` covers `Z/N` if and only
if `c(a, q) := Σ_{i=0}^{a} ν(2^i·q) ≥ q`. (Each non-dyadic class covers exactly one residue of `Z/q` in the nested picture;
(i) gives `c/q ≥ 1`, (ii) assigns `q` distinct residues.) Consequently:

| `q` | `ν(q), ν(2q), ν(4q), ν(8q), ν(16q)` | `c(a, q)` | least `a` | least period `2^a q` |
|---|---|---|---|---|
| 3 | 1, **0**, 1, 1, 2 | 1, 1, 2, 3 | 3 | **24** |
| 5 | 1, 1, 1, 1, 1 | 1, 2, 3, 4, 5 | 4 | **80** |
| 7 | 1, 1, 2, 1, 2 | 1, 2, 4, 5, 7 | 4 | **112** |

The zero in the first row is Zsigmondy's exception: with a prime of order 6 the period 12 would cover, so **the smallest
covering period 24 for `k·2^n ± 1` exists because `2^6 − 1 = 63` has no primitive prime divisor** (compare the hand
argument: classes modulo 2, 3, 4, 12 cover at most five of the six odd residues of `Z/12`). For `a ≥ 6` every `N = 2^a·q`
is covered already by `W(64) ⊆ W(N)`. Since `ν(2^i q) ≥ 1` for all `i` with `2^i q ≠ 6`, `c(a, q) ≥ a + 1 − [q = 3]`, so the
least period for `q` is at most `2^(q−1)·q` (`2^q·q` for `q = 3`): every odd prime is the odd part of some covering
period. ∎

**Lemma D (one prime, two signs).** A prime `p ∤ k` divides both `k·2^n − 1` (some `n`) and `k·2^m + 1` (some `m`) if and
only if `ord_p(2)` is even; then the two phases differ by `ord_p(2)/2`. (`k·2^n ≡ 1`, `k·2^m ≡ −1` give `2^(m−n) ≡ −1`, so
`−1 ∈ ⟨2⟩`, i.e. the order is even; conversely `2^(ord/2) ≡ −1`.) This is the general form of CLOCKS.md's remark that 7
(`ord 3`) blocks a Brier use of `W(24)`. ∎

## 2. The exact tables (FINITE-EXACT)

**Base 2, `N ≤ 150`** (`wheels3_base2_150.out`; exact search branching on the smallest uncovered residue). `W(N)` covers
`Z/N` exactly for

`N ∈ {24, 36, 48, 60, 64, 72, 80, 84, 96, 108, 112, 120, 128, 140, 144}`

and for no other `N ≤ 150`, where the five periods the first search left undecided are settled as follows: `N = 90, 126,
150` have `a = 1`, so the single dyadic class leaves exactly one coset and the exact problem *is* the nested criterion,
which fails; `N = 100 = 4·25` fails by Theorem B (i) (`Σ 1/m′ = 3/5 + 7/25 = 0.88 < 1`); `N = 132 = 4·33` does not cover: the two dyadic classes (modulo 2 and modulo 4) are either nested — then the exact problem is the odd wheel on `Z/33` with odd parts `{3: 2, 11: 5, 33: 5}`, which the exact search refutes — or overlapping, leaving a class modulo 2 of 66 residues against a capacity of 65 (a prime whose order is divisible by 4 reaches only one of the two cosets); the inconclusive long search is `exact_132.out`.
The nested-dyadic criterion of Theorem B (ii) agrees with the exact answer at **every** `N ≤ 150` (pre-registered P-W1,
confirmed; see PREREGISTERED.md). The density `Σ_{p ∈ W(N)} 1/ord_p(2)` alone is not the criterion: it exceeds 1 at
`N = 12, 18, 20, 28, 30, 40, 42, 44, 54, 56, 66, 70, 78, 88, 90, 100, 102, 104, 110, 132, 136, 150`, none of which
covers. `N = 64` and `128` cover with no odd part at all because `ν(64) = 2` (Theorem A).

**Base 3, `N ≤ 120`** (`wheels2_base3_120.out`; wheel = odd primes `p` with `ord_p(3) | N`). The coverable periods are
`N ∈ {48, 72, 90, 96, 108, 120}`. Base 3 lacks the modulus **2** (`3 + 1 = 4` is a power of two, the `n = 2` exception of
Zsigmondy) while base 2 lacks **6**: the two bases' wheels are missing complementary small moduli, and the least
covering period moves from 24 to 48. Pre-registered P-W3 (base-3 periods = multiples of 24 that are ≥ 48, plus 90)
is **refuted by one period**: 108 = 4·27 also covers.

**Self-dual wheels.** By Lemma D only even-order primes can serve both signs of one coefficient; a single wheel serves
both signs ("self-dual") iff residues `r_p` can be chosen so that both sheets `{r_p}` and `{r_p + ord_p/2}` cover `Z/N`.
Exact search (`brier_sheets.py`): no self-dual cover for any `N ≤ 150` that the search decided (24, 40, 48, 56, 60, 64,
72, 80, 84, 112, 128 all fail; 96 and 140 fail (`brier_96.out`, `brier_140.out`); 120 is undecided after `2·10^8` nodes (`brier_120.out`)). The known Brier numbers use two
different sub-wheels of one period (e.g. the period-720 set of Filaseta–Finch–Kozek contains the odd-order primes 7, 31,
73, 151, each serving one sign), which is consistent with this.

## 3. Conjectures and questions (new constructions, as asked)

- **W1 (nested is without loss).** For every `N`, `W(N)` covers `Z/N` iff the odd wheel covers `Z/m` (Theorem B (ii) is
  also necessary). True for all `N ≤ 150`. A proof would need that a non-nested choice of dyadic classes never helps;
  the coset argument of B (i) shows it never helps *in density*, not yet *in structure*.
- **W2 (least period for a prime odd part).** `a(q) := min{a : c(a, q) ≥ q}` satisfies `a(3) = 3`, `a(5) = 4`, `a(7) = 4`,
  and `a(q) ≤ q − 1` (Corollary C). Heuristically `ν(2^i q)` grows like `log(2^i q)`, so `c(a, q) ≈ a²/2 + a·log q` and
  `a(q) ≪ √q`. Question: compute `a(11), a(13), a(17)` (needs the primitive parts of `2^(2^i q) − 1` to be factored).
- **W3 (self-dual wheels).** Is there any period whose wheel is self-dual? None below 150 among the decided cases. If
  none exists, every Brier number needs two sub-wheels, which is how the known examples are built.
- **W4 (other bases).** For base `b` the missing moduli are exactly Zsigmondy's exceptions: `n = 2` when `b + 1` is a power
  of two (`b = 3, 7, 15, 31, …`) and `n = 6` for `b = 2`. Question: is the least covering period in base `b` determined by
  the first level where the base-`b` analogue of `c(a, q)` reaches `q`, as in Corollary C?
- **W5 (the Riesel pair as a two-prime phenomenon).** `W(24) = {p : (p − 1) | 24, p odd} ∪ {17, 241}`: the four small primes
  are forced into every wheel of period divisible by 24 by `(p − 1) | 24 ⇒ ord_p(2) | 24`; 17 and 241 are the two primes that
  make the covering possible (`c(3, 3) = 3` needs 7, 13, 241; the dyadic chain needs 3, 5, 17).

## 4. The Riesel wheel against {4, 2, 5, 7, 7, 13} and the Leech lattice

Facts (receipt `leech_numerics`): Riesel primes `(3, 5, 7, 13, 17, 241)`, orders `(2, 4, 3, 12, 8, 24)`, phases of
`509203·2^e − 1`: `(0, 1, 2, 7, 7, 3)`; `∏ primes = 5592405 = (2^24 − 1)/3`. Pure-A Niemeier lattices `A_n^k`
(`nk = 24`) have Coxeter numbers `h = n + 1 ∈ {2, 3, 4, 5, 7, 9, 13, 25} = {d + 1 : d | 24}`; those with `h − 1 | 12` are
`{2, 3, 4, 5, 7, 13}`, the owner's earlier set (THM-4838 typed its "twelve" reading NUMEROLOGY). `196560 = 2^4·3^3·5·7·13
= 48·(2^12 − 1)`; `|Co_0| = 2^22·3^9·5^4·7^2·11·13·23`.

**What is structural.** The primes `{2, 3, 5, 7, 13}` are exactly the primes with `(p − 1) | 24`, i.e. the von
Staudt–Clausen primes of the Bernoulli number `B_24` (denominator 2730 = 2·3·5·7·13). Three sets coincide for one
reason: (a) the prime Coxeter numbers of the pure-A Niemeier lattices are `{h = d + 1 prime : d | 24}` = these primes;
(b) the odd ones lie in the Riesel wheel because `(p − 1) | 24` forces `ord_p(2) | 24`; (c) they divide `196560` and
`|Co_0|` because the Leech lattice's numerology is built on 24. The common mechanism is the divisor lattice of 24, not
a Leech-lattice action on the covering. The covering structure — which classes, the two patterns, the role of 17 and
241, the failure of period 12 — is invisible to the lattice side and is governed by Theorem B and Zsigmondy's exception.

**What is coincidence.** `196560 = 48·4095` with 4095 the number of nonzero Golay codewords (the three Leech shell
shapes `97152 + 98304 + 1104` do not factor through 4095); the literal multiset `{4, 2, 5, 7, 7, 13}` (sum 38, the
two-fragment constant of THM-5253; a repeated 7 like the Riesel phases `7, 7`): no mechanism, NUMEROLOGY under
MISTAKE-559. Nothing here supports a Leech-lattice construction of covering sets, and nothing is needed: Theorem B
produces every covering period from the orders of 2 alone.
