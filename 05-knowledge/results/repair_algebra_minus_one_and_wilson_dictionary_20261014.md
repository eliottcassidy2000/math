# The repair algebra at the −1 anchor, the Wilson / {2,4,8} / Fermat-multiplier dictionary, and what has escaped us (THM-4640, THM-4641, HYP-9270)

Session `mac-mini-2026-10-14-wilson`, 2026-10-14. Independently audited the same day (audit A, `04-computation/experiments/repair_algebra_20261014/audit_A/REPORT.md`). Its corrections are applied below and logged as MISTAKE-606.

Owner's brief: "build the repair algebra at the −1 anchor to enumerate primitive collisions, and assemble the Wilson/{2,4,8}/Fermat-multiplier dictionary with exact checks. and mainly try to deeply understand the underlying themes that have escaped our grasps so far" — with focus on Collatz, Hadwiger–Nelson, the completeness of the Wilson primes {5, 13, 563}, and the planar Jacobian conjecture.

Scripts are in `04-computation/experiments/repair_algebra_20261014/`.

## 0. How the brief was read

* The **−1 anchor** is the 2-adic fixed point of `f_1(x) = (3x+1)/2`, the end of the Mersenne direction. **Collisions at −1** are THM-4555's `f_u(−1) = f_(u')(−1)`, which are exactly its uniform trailing-ones switches.
* A **repair** is the term that makes two different words meet: the ladder offset `(4^i − 1)/3` of THM-4601's `+1` compiler, transported to −1. The **repair algebra** is the set of all such relations, closed under suffixes and padding (THM-4640 (2)).
* **{2, 4, 8}** is read in the way the computations bear out: the powers of 2 adjacent to a power of 3 (`2 = 3 − 1`, `4 = 3 + 1`, `8 = 9 − 1`). They turn out to index its two-position relations and its padding relation, and the integer cycles of unit denominator. They are not generators: the carries do not generate the collisions.
* The **Fermat multiplier** has two exact meanings in the corpus:
  * the multiplier in `a^(p−1) = 1 + p·q_p(a)` (Fermat quotients);
  * the character `χ_8(p) = (2/p)` through which THM-307's tournament count sees the Wilson and Fermat quotients.

## 1. The repair algebra at −1 (THM-4640)

**Collisions are coincidences of increasing {2,3}-representations.**
* For a reduced word `u` with head partial sums `S_j`, `f_u(−1) = K(s)/2^(A−1)`, where `K(s) = −3^m + Σ 3^(m−j) 2^(e_j)` and `e_j = S_j − 1` is increasing.
* So two words collide iff they have the same total and the same `K`.
* The first generation is the row of reciprocal anchors: `K((a)) = 2^a − 3 = 1/c_a`. That is `−1, 1, 5, 13, 29, 61, …`, the reciprocals of the anchors of THM-4600, including the two Wilson primes 5 and 13.
* The formula is classical (Böhm–Sontacchi 1978, Lagarias 1985: `K = (B_u − 3^k)/2`), and the collision criterion is already in THM-4555's Remarks. What is new here is the algebra built on it.

**Positivity and rigidity.**
* Canonical numerators are positive: `K(s) ≥ 3^(m+1) − 2^(m+2) ≥ 1`, and `2^(e_m) ≤ K(s) + 3` with equality only for one-letter heads.
* So every value has finitely many representations, all with exponents `≤ log_2(K + 3)`.
* The smallest value with two canonical representations is the anchor `125 = 2^7 − 3`. For `a ≤ 31`, the anchor `2^a − 3` has a second representation only for `a ∈ {7, 13, 17, 21, 28, 31}`. So 5 and 13 are rigid, but only because they lie below 125 (`anchor_rigidity.py`).

**Primitive collisions are ladders.** A collision whose words end with different letters is primitive; every collision reduces to one by cancelling the common suffix. Its heads satisfy

    f_h'(−1) = 4^i f_h(−1) + (4^i − 1)/3,

THM-4601's ladder with the anchor moved from `+1` to `−1`. Each primitive collision is an infinite ladder of word collisions, hence a uniform switch family (THM-4555 (iv)) and, through the root twin, a second family on reset-2 sources.

**Parity.** Ladder gaps `2i` are even. `K mod 8` determines whether a canonical representation has even length and whether it starts with exponent 2. So canonical length differences, and the switch distances `D`, are even; root twins give odd `D`.

**The Gersonides carries.**
* The only two-position relations between exponent sequences are:
  * the unit carry `(b+1, b) ~ (b, b+2)` (`3(2 − 1) = 2^2 − 1`);
  * the Catalan carry `(b+3, ·, b) ~ (b, ·, b+6)` (`3^2(2^3 − 1) = 2^6 − 1`).
* The proof reduces to `3^Δ(2^x − 1) = 2^y − 1`, then uses `x | y` and Bang–Zsigmondy.
* With the padding relation `f_1(−1) = −1`, these relations are indexed by `{2, 4, 8}`.
* Both carries need a descent. So they are **repairs** of inadmissible sequences, and no primitive collision between admissible words differs in fewer than three positions.
* The carries do **not** generate the primitive collisions: none of the 133 primitive pairs up to exponent 14 is reachable by carries. The first one, `125 = 2^7 − 3`, `(7) ~ (3, 4, 5)`, is the non-degenerate 6-term relation `2^7 + 27 = 3 + 2^5 + 3·2^4 + 9·2^3`.
* As unit sums, though, about 40% of primitive collisions are degenerate, mostly through the 3-term relation `4 − 3 − 1 = 0`: 737 of 1,749 canonical pairs with exponents `≤ 18` (audit A). The first is `259`, `(2, 8) ~ (2, 3, 5, 6)`.
* For each number of terms, the underlying vanishing `{2,3}`-unit sums come in finitely many non-degenerate types (S-unit theorem, cited).

**Census (exact).**

| max exponent `E` | 14 | 16 | 18 | 20 | 22 | 24 |
|---|---|---|---|---|---|---|
| primitive values `N(E)` | 132 | 462 | 1,666 | 6,129 | 23,037 | 87,836 |
| `N(E)·√E/2^E` | 0.0301 | 0.0282 | 0.0270 | 0.0261 | 0.0258 | 0.0256 |
| coverage without the root, `c(E)` | 0.0328 | 0.0375 | 0.0417 | 0.0454 | 0.0489 | 0.0521 |
| coverage with the root (reset switch) | 0.2718 | 0.2750 | 0.2778 | 0.2803 | 0.2826 | 0.2847 |

* Up to `E = 22`:
  * 20,393 primitive values have two canonical representations, 2,611 three and 32 four;
  * gaps `2i = 2, 4, 6` occur 24,732, 874 and 138 times.
* The smallest primitive values are `−1` (the root), `125`, `259`, `439`, `599`, `905`, `1063`, `1081`, `1225`, `1445`, `2053`, … (full list in `repair_algebra.out`).
* Audit A's independent census continues to `E = 29` (`N(29) = 2,563,623`). HYP-9270 records the growth `N(E) = 2^E E^(−β+o(1))`, with local exponent `β ≈ 0.47` at `E = 28, 29`. The coverage is a theorem (below).

**What it buys.**
* The primitive collisions are exactly the uniform trailing-ones switches.
* **The root is the reset switch.** The root collision is THM-4555's reset switch `n ↦ (n−1)/2`, which covers exactly 1/4 of odd sources.
* **Canonical and twin switches.** Canonical switches lie inside the reset domain. Their twins act on the reset-2 sources, with exactly twice the canonical mass.
* **The coverage, with and without the root.** Without the root, `c(E) ≤ 3/8`. With it, the coverage is 0.2826 at exponent 22, i.e. 56.5% of trailing-ones sources. By THM-4581 6(e) it tends to 1/2: almost every trailing-ones source has a −1 switch (THM-4640 (8)).
* **Correction.** A pre-audit draft omitted the root, reported 5% coverage growing logarithmically, and concluded that the −1 mechanism cannot descend most sources. That conclusion was false (MISTAKE-606).
* Sources with no trailing one are never switched, but they descend in one step. THM-4601 treats one residual branch of the trailing-ones sources, the one whose post-run point lies 2-adically near `+1`, with a second anchor.

## 2. The Wilson / {2,4,8} / Fermat-multiplier dictionary (THM-4641)

| row | content | type |
|---|---|---|
| 1 | `{2, 4, 8}` are the powers of 2 adjacent to a power of 3 | KNOWN (Gersonides) |
| 2 | they index the padding relation and the two two-position relations (the carries) of the repair algebra at −1. These repair descents but do not generate the collisions. `K((a)) = 2^a − 3 = 1/c_a` | PROVED |
| 2b | they also index the Collatz cycles with unit denominator `|2^A − 3^k| = 1`: `−1`, `+1`, `{−5, −7}` (and `0`); the only other known integer cycle (through `−17`) has denominator 139 | PROVED + FINITE-EXACT |
| 3 | Wilson's `−1` is the product of a group's involutions; for `⟨2⟩ mod p` the product is `−1` iff `ord_p(2)` is even | KNOWN / PROVED |
| 4 | Lerch: `w_p ≡ Σ q_p(a)`; `(q_p(2), q_p(4), q_p(8)) ≡ (1, 2, 3) q_p(2)` | KNOWN |
| 5 | `f_p ≡ χ_8(p) q_p(2)/2`, so `T(p) ≡ χ_8(p)(q_p(2)/2 − w_p) − A_p (mod p)`. Tournaments see both quotients through the same Fermat multiplier `χ_8(p)`; re-derived from exact Burnside counts for all primes `p < 48` | PROVED |
| 6 | 5 and 13 are the only Wilson primes among the primes `2^a − 3`, `a ≤ 44` | COINCIDENCE (exact to `a = 24`, cited beyond) |
| 7 | 2 is a primitive root mod 5, 13, 563 (`2^a − 3 ≡ 5 mod 8` forces part of this) | COINCIDENCE (baseline 0.052) |
| 8 | `p − 2 ∈ {3, 11, 561}` absolutely Fermat | COINCIDENCE (baseline about 0.09 with size-dependent rates) |
| 9 | 563, 561, 1093, 3511 are not numerators at −1 | FINITE-EXACT |
| 10 | `λ(684106401) + 1 = 27·563` | COINCIDENCE |
| 11 | Wilson primes are expected to be infinite (`log log x`). Infinitude in `F_q[t]` is a theorem (Thakur), by an algebraic criterion that does not test the count | HEURISTIC + CITED |

The structural core is rows 1–5: Gersonides (as repairs and as unit-denominator cycles), Lerch, the group-theoretic `−1`, and the tournament character. Rows 6–10 are small-number coincidences; their baselines are computed in the script.

## 3. What has escaped us, problem by problem

The four problems share one shape. The repository proves statements about **almost all** objects or about **finite certificates**. The open problems quantify over **every** integer, prime or polynomial, and the remaining set is invisible to measure and too large for finite search.

* **Collatz.**
  * **What we have.**
    * Coalescence holds almost surely (THM-4581), and `n`, `n + h` merge for a set of density one (THM-4627). Coalescence fails for expanding maps (THM-4606/4607).
    * Many rule banks are proved: switches, ladders, anchors, block certificates.
    * The −1 switches cover almost every trailing-ones source (THM-4640 (8)).
  * **The gap.**
    * Every rule bank covers 2-adic cylinders. Their union can tend to measure one, but at each finite depth it misses a set of positive measure, and the limit misses a null set.
    * The anchors `c_a = 1/(2^a − 3)` and the rational cycle points are 2-adic points whose neighbourhoods need ever deeper rules.
  * Integers are Haar-null one at a time. A proof must show that positive integers avoid the uncovered set, and that is a Diophantine statement about `2^A` versus `3^k`.
  * Today's tools bound cycles only. Eliahou's continued-fraction budget with Barina's verification gives `k ≥ 72,057,431,991` odd steps (THM-4633), and Baker-type linear forms give weaker bounds. Nothing comparable is known for divergence.
  * The repair algebra makes the obstruction explicit: collisions at −1 are `{2,3}`-unit relations, the two-term ones are exactly Gersonides', and the rest are higher S-unit relations (about 40% degenerate, mostly through `4 − 3 − 1 = 0`) with no known generating set.
  * Covering every integer `n > 1` by a one-step descent or a −1 switch would prove Collatz outright: a switch merges `n` with a smaller integer. So the exceptional null set of THM-4640 (8) carries the whole difficulty.
  * **Anchor obstruction (PROVED, elementary).** Call a certificate *uniform* if it is a clopen 2-adic cylinder together with a descent rule valid for every positive integer in it. No uniform certificate contains the 2-adic point 1: the rule would have to make `n = 1` descend.
    * A finite bank has a clopen union. So if the bank contained every odd `n > 1`, it would contain their 2-adic closure, which includes 1.
    * Hence every finite bank misses all `n ≡ 1 (mod 2^K)` for some `K`, and certificates near the trivial cycle must have unbounded depth.
    * Near `+1` the infinite family is explicit and trivial: every `n ≡ 1 (mod 4)` with `n > 1` has `U(n) ≤ (3n+1)/4 < n`. So the cylinders `n ≡ 1 + 2^j (mod 2^(j+1))`, `j ≥ 2`, cover every such source except 1 itself.
    * So the anchor obstruction is real but harmless. The genuine gap is the generic part, the stopping-time problem: control of `3^k` against `2^A` along orbits with no anchor structure. That is the same Diophantine question that bounds cycles (THM-4633).
* **Wilson primes.**
  * Nothing structural forces `W_p ≢ 0 (mod p)`. Lerch, Glaisher and the lift trees all express `W_p` as a sum of "random-looking" Fermat quotients.
  * Completeness of `{5, 13, 563}` is expected to be **false** (heuristic `log log x`; infinitude in `F_q[t]` is a theorem, though by an algebraic criterion). It is not a viable proof target.
  * The viable projects are an exact search beyond `2·10^13` (the Costa–Gerbicz–Harvey algorithm) or structural analogues, such as Wilson primes in other rings.
* **Hadwiger–Nelson.**
  * The repository works in field planes (`Q(ζ55)`, THM-4614/4628), where residue colourings are rigid, and in finite Beckman–Quarles (THM-4624).
  * The lower bound 5 rests on explicit graphs. The smallest known 5-chromatic unit-distance graph has 509 vertices (Parts; Heule's SAT pipeline reached 510), and a 2026 preprint (Haugland, arXiv:2608.04542) still cites 509. The upper bound 7 rests on one tiling.
  * Moving either bound needs a construction, not a density argument.
* **Planar Jacobian.**
  * The logic-lens program (THM-4620–4626) reformulates JC2 through Keller dual fields. HYP-9255 is equivalent to JC2.
  * The verified range is degree `≤ 100` (Moh 1983), extended to 104 (Nguyen, arXiv:1902.05923) and 108 (Guccione et al., arXiv:2204.14178). Every known partial result needs an idea specific to the Newton polygon.
  * No reformulation so far has produced a new verified case.

**What could be completed in weeks.** No proof of any of the four headline problems is realistic in that time, and the Wilson one is probably false. Notable but realistic targets, best first:

1. **Hadwiger–Nelson, constructive.** A 5-chromatic unit-distance graph with fewer than 509 vertices (the Parts record), or a Moser-spindle-free one with fewer than 1,441 vertices. That record is Heule's ("Odd-distance virtual edges in unit-distance graphs", Geombinatorics 31(2), 2021), a SAT reduction of the 64,513-vertex graphs of Voronov, Neopryatnaya and Dergachev (arXiv:2106.11824). It is cited in v5 of Haugland's arXiv:2608.04542 (6 October 2026), whose structured spindle-free graph has 2,131 vertices; Haugland's v5 also cites a certified reduction of that graph to 1,501 vertices (K. T. Wong, Zenodo, 2026). The 509 record is Parts, arXiv:2010.12665. These are checkable by SAT (`pysat` is installed), have a clear success criterion, and would be noticed. The 509 record came from heavy SAT-guided search by Heule and Parts, so success on this machine is uncertain; the spindle-free variant has had less optimization (one SAT reduction in 2021, and 2026 work that has not beaten it), so it is likely softer. These references were read from the arXiv abstract and HTML pages only.
2. **Collatz, structural.**
   * Make THM-4640 (8) quantitative: the uncovered reset-2 mass at depth `E`. THM-4581 (4) suggests order `E^(−1/2)` up to logarithms.
   * Describe the exceptional null set, and find which integers fall in it.
   * Covering every integer would prove Collatz, so the realistic goal is the rate and the structure, not completeness.
3. **Wilson analogues.** For example, Wilson primes in `Z[i]` or `Z[ω]` with an exact search, where nothing is known numerically in this corpus.
4. **JC2.** Keep the dual-field lens, aiming at one new verified family (a Newton-polygon class) rather than the full conjecture.

## 4. Reproduction

All runs use `04-computation/experiments/repair_algebra_20261014/`, Python 3.10 with numpy and sympy.

| command | output | runtime |
|---|---|---|
| `python3 repair_algebra.py 16 24 14` | `repair_algebra.out`: word-level validation, census, carries, padded differences, integer switches; 82,115 exact checks | about 1 min |
| `python3 switch_coverage.py 24` | `switch_coverage.out`: exact Haar coverage, split into canonical, twin and root-inclusive parts; the identities `twin = 2·canonical` and `with root = 1/4 + twin` and the bound `3/8` are asserted as rationals | about 2 min |
| `python3 gaps_census.py 22` | `gaps_census.out`: gaps, length differences, representations per value | about 2 min |
| `python3 carry_derivation.py 14 30` | `carry_derivation.out`: the carries do not derive primitive collisions | about 1 min |
| `python3 wilson_248_fermat_dictionary.py` | `wilson_248_fermat_dictionary.out`: 1,835 exact checks | about 10 s |
| `python3 anchor_rigidity.py 20 31 16` | `anchor_rigidity.out`: positivity lemma on 524,287 representations, complete anchor search, recursion = brute force below `2^15`; 540,712 exact checks | about 1 s |

Audit A's independent scripts and outputs (`a1`–`a9`) and its report are in `audit_A/`.
