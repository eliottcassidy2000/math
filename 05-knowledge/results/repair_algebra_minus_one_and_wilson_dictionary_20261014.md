# The repair algebra at the −1 anchor, the Wilson / {2,4,8} / Fermat-multiplier dictionary, and what has escaped us (THM-4640, THM-4641, HYP-9270)

Session `mac-mini-2026-10-14-wilson`, 2026-10-14. Not yet independently audited.

Owner's brief: "build the repair algebra at the −1 anchor to enumerate primitive collisions, and assemble the Wilson/{2,4,8}/Fermat-multiplier dictionary with exact checks. and mainly try to deeply understand the underlying themes that have escaped our grasps so far" — with focus on Collatz, Hadwiger–Nelson, the completeness of the Wilson primes {5, 13, 563}, and the planar Jacobian conjecture.

Scripts are in `04-computation/experiments/repair_algebra_20261014/`.

## 0. How the brief was read

* The **−1 anchor** is the 2-adic fixed point of `f_1(x) = (3x+1)/2`, the end of the Mersenne direction. **Collisions at −1** are THM-4555's `f_u(−1) = f_(u')(−1)`, which are exactly its uniform trailing-ones switches.
* A **repair** is the term that makes two different words meet: the ladder offset `(4^i − 1)/3` of THM-4601's `+1` compiler, transported to −1. The **repair algebra** is the algebra these relations generate.
* **{2, 4, 8}** is read in the way the computations bear out: the powers of 2 adjacent to a power of 3 (`2 = 3 − 1`, `4 = 3 + 1`, `8 = 9 − 1`). They turn out to index the generators of the repair algebra.
* The **Fermat multiplier** has two exact meanings in the corpus:
  * the multiplier in `a^(p−1) = 1 + p·q_p(a)` (Fermat quotients);
  * the character `χ_8(p) = (2/p)` through which THM-307's tournament count sees the Wilson and Fermat quotients.

## 1. The repair algebra at −1 (THM-4640)

**Collisions are coincidences of increasing {2,3}-representations.**
* For a reduced word `u` with head partial sums `S_j`, `f_u(−1) = K(s)/2^(A−1)`, where `K(s) = −3^m + Σ 3^(m−j) 2^(e_j)` and `e_j = S_j − 1` is increasing.
* So two words collide iff they have the same total and the same `K`.
* The first generation is the anchor row: `K((a)) = 2^a − 3 = 1/c_a`. That is `−1, 1, 5, 13, 29, 61, …`, which contains every anchor of THM-4600, including the two Wilson primes 5 and 13.

**Primitive collisions are ladders.** A collision whose words end with different letters is primitive; every collision reduces to one by cancelling the common suffix. Its heads satisfy

    f_h'(−1) = 4^i f_h(−1) + (4^i − 1)/3,

THM-4601's ladder with the anchor moved from `+1` to `−1`. Each primitive collision is an infinite ladder of word collisions, hence a uniform switch family (THM-4555 (iv)) and, through the root twin, a second family on reset-2 sources.

**Parity.** Ladder gaps `2i` are even. `K mod 8` determines whether a canonical representation has even length and whether it starts with exponent 2. So canonical length differences, and the switch distances `D`, are even; root twins give odd `D`.

**The Gersonides carries.**
* The only two-position relations between exponent sequences are:
  * the unit carry `(b+1, b) ~ (b, b+2)` (`3(2 − 1) = 2^2 − 1`);
  * the Catalan carry `(b+3, ·, b) ~ (b, ·, b+6)` (`3^2(2^3 − 1) = 2^6 − 1`).
* The proof reduces to `3^Δ(2^x − 1) = 2^y − 1`, then uses `x | y` and Bang–Zsigmondy.
* With the anchor relation `f_1(−1) = −1`, the generators are indexed by `{2, 4, 8}`.
* Both carries need a descent. So they are **repairs** of inadmissible sequences, and no primitive collision between admissible words differs in fewer than three positions.
* The carries do **not** generate the primitive collisions: none of the 133 primitive pairs up to exponent 14 is reachable by carries. The first one, `125 = 2^7 − 3`, `(7) ~ (3, 4, 5)`, rests on `27 = 3 + 8 + 16`.
* For each number of terms, the underlying vanishing `{2,3}`-unit sums come in finitely many non-degenerate types (S-unit theorem, cited).

**Census (exact).**

| max exponent `E` | 14 | 16 | 18 | 20 | 22 | 24 |
|---|---|---|---|---|---|---|
| primitive values `N(E)` | 132 | 462 | 1,666 | 6,129 | 23,037 | 87,836 |
| `N(E)·√E/2^E` | 0.0301 | 0.0282 | 0.0270 | 0.0261 | 0.0258 | 0.0257 |
| switch coverage `c(E)` | 0.0328 | 0.0375 | 0.0417 | 0.0454 | 0.0489 | — |

* Up to `E = 22`:
  * 20,393 primitive values have two canonical representations, 2,611 three and 32 four;
  * gaps `2i = 2, 4, 6` occur 24,732, 874 and 138 times.
* The smallest primitive values are `−1` (the root), `125`, `259`, `439`, `599`, `905`, `1063`, `1081`, `1225`, `1445`, `2053`, … (full list in `repair_algebra.out`).
* HYP-9270 records the growth laws, `N(E) ~ 0.0255·2^E/√E` and `c(E) ≈ 0.036 ln E`.

**What it buys.**
* The primitive collisions are exactly the uniform trailing-ones switches.
* Their Haar coverage of odd sources is about 5% at exponent 22 and grows only logarithmically. So the −1 mechanism cannot by itself descend most sources.
* The `+1` anchor and the anchor row `2^a − 3` are where the rest of the mass lives.

## 2. The Wilson / {2,4,8} / Fermat-multiplier dictionary (THM-4641)

| row | content | type |
|---|---|---|
| 1 | `{2, 4, 8}` are the powers of 2 adjacent to a power of 3 | KNOWN (Gersonides) |
| 2 | they index the generators of the repair algebra at −1; anchors `−1, 1, 1/5`; `K((a)) = 1/c_a` | PROVED |
| 2b | they also index the Collatz cycles with unit denominator `|2^A − 3^k| = 1`: `−1`, `+1`, `{−5, −7}` (and `0`); the only other known integer cycle (through `−17`) has denominator 139 | PROVED + FINITE-EXACT |
| 3 | Wilson's `−1` is the product of a group's involutions; for `⟨2⟩ mod p` the product is `−1` iff `ord_p(2)` is even | KNOWN / PROVED |
| 4 | Lerch: `w_p ≡ Σ q_p(a)`; `(q_p(2), q_p(4), q_p(8)) ≡ (1, 2, 3) q_p(2)` | KNOWN |
| 5 | `f_p ≡ χ_8(p) q_p(2)/2`, so `T(p) ≡ χ_8(p)(q_p(2)/2 − w_p) − A_p (mod p)`. Tournaments see both quotients through the same Fermat multiplier `χ_8(p)`; re-derived from exact Burnside counts for all primes `p < 48` | PROVED |
| 6 | 5 and 13 are the only Wilson primes among the primes `2^a − 3`, `a ≤ 44` | COINCIDENCE (exact to `a = 24`, cited beyond) |
| 7 | 2 is a primitive root mod 5, 13, 563 (`2^a − 3 ≡ 5 mod 8` forces part of this) | COINCIDENCE (baseline 0.052) |
| 8 | `p − 2 ∈ {3, 11, 561}` absolutely Fermat | COINCIDENCE (baseline 0.10 per prime) |
| 9 | 563, 561, 1093, 3511 are not numerators at −1 | FINITE-EXACT |
| 10 | `λ(684106401) + 1 = 27·563` | COINCIDENCE |
| 11 | Wilson primes are expected to be infinite (`log log x`); a theorem in `F_q[t]` (Thakur) | HEURISTIC + CITED |

The structural core is rows 1–5: Gersonides (as repairs and as unit-denominator cycles), Lerch, the group-theoretic `−1`, and the tournament character. Rows 6–10 are small-number coincidences; their baselines are computed in the script.

## 3. What has escaped us, problem by problem

The four problems share one shape. The repository proves statements about **almost all** objects or about **finite certificates**. The open problems quantify over **every** integer, prime or polynomial, and the remaining set is invisible to measure and too large for finite search.

* **Collatz.**
  * **What we have.** Coalescence and descent hold almost surely (THM-4581, THM-4606/4607, THM-4627), and many rule banks are proved (switches, ladders, anchors, block certificates).
  * **The gap.** Every rule bank covers 2-adic cylinders. Their union grows like `log` (this note) or tends to measure one, yet it never covers every point: the anchors `c_a = 1/(2^a − 3)` and the rational cycle points are 2-adic points whose neighbourhoods need ever deeper rules.
  * Integers are Haar-null one at a time. A proof must show that positive integers avoid the uncovered set, and that is a Diophantine statement about `2^A` versus `3^k`.
  * Today's tool for that is Baker's theory. It bounds cycle lengths (THM-4633: about `7·10^10` odd steps from the budget alone), but it says nothing about divergence.
  * The repair algebra makes the obstruction explicit: collisions at −1 are `{2,3}`-unit relations, the two-term ones are exactly Gersonides', and the rest are higher S-unit relations with no generating pattern.
  * **Anchor obstruction (PROVED, elementary).** Call a certificate *uniform* if it is a clopen 2-adic cylinder together with a descent rule valid for every positive integer in it. No uniform certificate contains the 2-adic point 1: the rule would have to make `n = 1` descend.
    * A finite bank has a clopen union. So if the bank contained every odd `n > 1`, it would contain their 2-adic closure, which includes 1.
    * Hence every finite bank misses all `n ≡ 1 (mod 2^K)` for some `K`, and certificates near the trivial cycle must have unbounded depth.
    * Near `+1` the infinite family is explicit and trivial: every `n ≡ 1 (mod 4)` with `n > 1` has `U(n) ≤ (3n+1)/4 < n`. So the cylinders `n ≡ 1 + 2^j (mod 2^(j+1))`, `j ≥ 2`, cover every such source except 1 itself.
    * So the anchor obstruction is real but harmless. The genuine gap is the generic part, the stopping-time problem: control of `3^k` against `2^A` along orbits with no anchor structure. That is the same Diophantine question that bounds cycles (THM-4633).
* **Wilson primes.**
  * Nothing structural forces `W_p ≢ 0 (mod p)`. Lerch, Glaisher and the lift trees all express `W_p` as a sum of "random-looking" Fermat quotients.
  * Completeness of `{5, 13, 563}` is expected to be **false** (heuristic `log log x`; theorem in `F_q[t]`). It is not a viable proof target.
  * The viable projects are an exact search beyond `2·10^13` (the Costa–Gerbicz–Harvey algorithm) or structural analogues, such as Wilson primes in other rings.
* **Hadwiger–Nelson.**
  * The repository works in field planes (`Q(ζ55)`, THM-4614/4628), where residue colourings are rigid, and in finite Beckman–Quarles (THM-4624).
  * The lower bound 5 rests on explicit graphs. The smallest known 5-chromatic unit-distance graph has 509 vertices (Parts; Heule's SAT pipeline reached 510), and a 2026 preprint (Haugland, arXiv:2608.04542) still cites 509. The upper bound 7 rests on one tiling.
  * Moving either bound needs a construction, not a density argument.
* **Planar Jacobian.**
  * The logic-lens program (THM-4620–4626) reformulates JC2 through Keller dual fields. HYP-9255 is equivalent to JC2.
  * The classical verified range is degree `≤ 100` (Moh). Every known partial result needs an idea specific to the Newton polygon.
  * No reformulation so far has produced a new verified case.

**What could be completed in weeks.** No proof of any of the four headline problems is realistic in that time, and the Wilson one is probably false. Notable but realistic targets, best first:

1. **Hadwiger–Nelson, constructive.** A 5-chromatic unit-distance graph with fewer than 509 vertices (the Parts record), or a smaller Moser-spindle-free one than Haugland's 2,131 vertices (arXiv:2608.04542). These are checkable by SAT (`pysat` is installed), have a clear success criterion, and would be noticed. The 509 record came from heavy SAT-guided search by Heule and Parts, so success on this machine is uncertain; the spindle-free variant is younger and likely softer.
2. **Collatz, structural.** A classification of the uncovered set of the anchor rule banks. Specifically: every 2-adic point not in the union of switch and ladder cylinders lies in the closure of the anchor orbits. That is a compactness statement; it is not a proof of Collatz.
3. **Wilson analogues.** For example, Wilson primes in `Z[i]` or `Z[ω]` with an exact search, where nothing is known numerically in this corpus.
4. **JC2.** Keep the dual-field lens, aiming at one new verified family (a Newton-polygon class) rather than the full conjecture.

## 4. Reproduction

All runs use `04-computation/experiments/repair_algebra_20261014/`, Python 3.10 with numpy and sympy.

| command | output | runtime |
|---|---|---|
| `python3 repair_algebra.py 16 24 14` | `repair_algebra.out`: word-level validation, census, carries, padded differences, integer switches; 82,115 exact checks | about 1 min |
| `python3 switch_coverage.py 22` | `switch_coverage.out`: exact Haar coverage | about 1 min |
| `python3 gaps_census.py 22` | `gaps_census.out`: gaps, length differences, representations per value | about 2 min |
| `python3 carry_derivation.py 14 30` | `carry_derivation.out`: the carries do not derive primitive collisions | about 1 min |
| `python3 wilson_248_fermat_dictionary.py` | `wilson_248_fermat_dictionary.out`: 1,827 exact checks | about 10 s |
