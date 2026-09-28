# The graph tower against the tournament tower (D18), Emma Lehmer's quintic and the Brocard triple: the squares 25, 121, 5041 are f(0), f(−1)², f(1)² for the quintic's conductor polynomial, Brocard's equation is a divisor involution with a missing fixed point, and the fives typed

**Session:** opus, `collatz-posets-zeta5-20260927` (eighth note), 2026-09-27.
**Owner's directive:** "pursue D18 iterated graph doubling against the tower;
consider the Lehmer 5 and how they might be analogous to previous work on
the 5 Fermat primes and the Platonic solids and other deep 5 structures;
consider how {4, 5, 7} are the only factorials one less than a square, and
those squares are {5, 11, 71}, and how these pairs of triples relate
abstractly with other triples we have studied; get a handle on the
underlying structure."
**Interpretation of "the Lehmer 5" (CORRECTED in the ninth note).** The owner
meant the *Lehmer five* `276, 552, 564, 660, 966`, the five smallest numbers
whose aliquot sequences are not known to terminate or cycle (named after D. H.
Lehmer); this note's reading as Emma Lehmer's quintic was wrong, and is logged
in MISTAKES. The coincidence below stands under its own name. The reading
taken here was **Emma Lehmer's simplest quintic** (Math. Comp. 50, 1988; Schoof–Washington
1988), the degree-5 member of the simplest-fields ladder, whose conductor
polynomial `f(n) = n^4 + 5n^3 + 15n^2 + 25n + 25` takes the values `11, 25,
71` at `n = -1, 0, 1`: exactly the owner's `5, 11, 71` (with `25 = 5^2`).
This note proceeds on that reading and says where it is exact.
**Inherits (cited):** the seventh note (Theorem 4, the Seidel doubling law;
Paley rows), THM-447/455/483 (the skew-Sylvester tower and its trans values
`3, 5, 7, 11`; the zigzag law), the Kuratowski note's triple dictionary
(KT-a "twins + container", KT-b "dual pair + self-dual"), the S10 note (the
five Fermat primes and the five Platonic solids are not the same five;
Gauss–Wantzel; THM-871's Fermat rungs), HYP-3771/3772 (the `(2,3,n)`
angle-defect spine), the third note (the golden field, `Q(zeta_5)`,
Petersen and `K_6` as folds). Classical facts from memory, verified where
stated: Brocard 1876/1885, Ramanujan 1913, Overholt 1993 (finitely many
solutions under abc), Berndt–Galway 2000 and Matson 2017 (no further
solution below `10^15`), Ogg's supersingular primes, the binary polyhedral
groups, Wilson primes `5, 13, 563`, Mordell's `((p-1)/2)! = +-1 mod p`.

**Status: PROVED elementary (Proposition 1's law re-verified along the
towers; Proposition 3's divisor-involution reading) + FINITE-EXACT (the
towers to order 80, the Lehmer quintic's conductors, discriminants and
Galois groups by sympy, Brocard solutions to `n = 60`, the half-Wilson
census to 2000, random benchmarks) + CITED + NUMEROLOGY typed + DIRECTION.
Collatz is not addressed; this leg is the tournament/Ramsey side of the
session. Independent audit OWED.** Script
`04-computation/experiments/collatz_doubling_tower_brocard_20260927.py`,
output beside it.

## 0. The answer in one paragraph

The graph tower behaves like the tournament tower. Iterating the Seidel
doubling from the pentagon `P_5` gives graphs of orders `5, 10, 20, 40, 80`
with `(omega, alpha) = (2,2), (3,3), (4,4), (5,6), (7,8)`, so `max(omega,
alpha) = 3, 4, 6, 8` at orders `10, 20, 40, 80`, at or below the best of a
dozen random graphs at every order (random means `4.3, 5.7, 7.2, 9.0`), just
as THM-455 found the tournament tower `3, 5, 7, 11` below the random median
at orders `7, 15, 31, 63`; the seventh note's law `omega(D) = cs`, `alpha(D)
= max(alpha + 1, sep)` holds at every level of every tower computed, and
`cs(D(G))` is strictly less than `cs(G) + sep(G)` (Proposition 1). The
Brocard triple has one exact structural reading and one exact numerical
coincidence. Exact: `n! + 1 = m^2` says the divisor involution `x -> n!/x`
of `n!` has an orbit `{m-1, m+1}` of diameter two around the missing fixed
point `sqrt(n!)` (which is never an integer): `24 = 4·6`, `120 = 10·12`,
`5040 = 70·72`; in the Kuratowski dictionary this is KT-b, a dual pair with
the self-dual point absent (Proposition 3). Coincidence, verified exactly:
the three squares are `f(0) = 25`, `f(-1)^2 = 121`, `f(1)^2 = 5041` for the
conductor polynomial of Lehmer's simplest quintic, whose members at `n = -1,
0, 1` are the cyclic quintic fields of conductors `11, 25, 71` (irreducible,
Galois group `C_5`, discriminants `11^4`, `5^8·7^2`, `23^2·71^4`, all checked);
the three roots `5, 11, 71` are supersingular primes, `71` the largest; `4!`
and `5!` are the orders of the binary tetrahedral and binary icosahedral
groups; `5` is a Wilson prime (`4! + 1 = 5^2`) and `11` satisfies the
half-Wilson congruence at level `p^2` (`5! = -1 mod 11^2`; below 2000 only
`11` and `47` do), while `7! + 1 = 71^2` has no Wilson reading. No mechanism
connects any of these, and the note says so. The fives are inventoried
against the repo: the five Platonic solids have one exact source (the sign
of `1/2 + 1/3 + 1/n - 1`), the five Fermat primes another (Gauss–Wantzel),
the golden `Q(sqrt5) ⊂ Q(zeta_5)` a third, and Lehmer's quintic is a fourth
(the degree, not a count); they are four different fives.

## 1. D18: the graph tower (PROVED law + FINITE-EXACT)

**Proposition 1.** Along the towers `D^k(G)` for `G = K_1, K_2, P_5, P_13`
(orders up to 80) the law `omega(D(G)) = cs(G)`, `alpha(D(G)) = max(alpha(G)
+ 1, sep(G))` holds at every level (all four parameters recomputed by
brute force), and `cs(D(G)) < cs(G) + sep(G)` in every case computed
(`K_3`: `4 < 6`; `C_5`: `4 < 6`; `P_13`: `5 < 8`).

| start | order | `omega` | `alpha` | `cs` | `sep` |
|---|---|---|---|---|---|
| `P_5` | 5, 10, 20, 40, 80 | 2, 3, 4, 5, 7 | 2, 3, 4, 6, 8 | 3, 4, 5, 7, 9 | 3, 4, 6, 8, 12 |
| `P_13` | 13, 26, 52 | 3, 4, 5 | 3, 4, 6 | 4, 5, 7 | 4, 6, 8 |
| `K_1` | 1, 2, 4, 8, 16, 32 | 1, 1, 2, 3, 4, 6 | 1, 2, 3, 4, 6, 8 | 1, 2, 3, 4, 6, 8 | 1, 2, 4, 6, 8, 12 |
| `K_2` | 2, 4, 8, 16, 32 | 2, 2, 3, 4, 5 | 1, 2, 3, 4, 6 | 2, 3, 4, 5, 7 | 2, 3, 4, 6, 8 |

Against the tournament tower (THM-455: `trans(T_7, T_15, T_31, T_63) = 3, 5,
7, 11`, increments `+2, +2, +4`) the graph tower from `P_5` has `max(omega,
alpha) = 3, 4, 6, 8` at orders `10, 20, 40, 80` (increments `+1, +2, +2`) and
the random benchmark `G(n, 1/2)` has mean `max(omega, alpha) = 4.3, 5.7, 7.2,
9.0` with ranges `[3,5], [5,7], [6,8], [8,10]`: the tower sits at the bottom
of the random range at every order, as the tournament tower does, and like
it is far from the Ramsey records (`R(9,9) > 80` from `D^4(P_5)`, against
`R(9,9) >= 565`). Why `cs(D(G)) < cs(G) + sep(G)`: a complete split subgraph
of `D(G)` is a clique `K ∪ I'` (clique `K` of `G` on the plain copy,
independent set `I` on the prime copy, `K` complete to `I`) joined to an
independent set `J ∪ L'` (independent `J` on the plain copy, clique `L` on
the prime copy, no edges between `J` and `L`), with the eight cross
conditions `K–J` complete, `K–L` complete, `I–J` complete, `I–L`
anticomplete; the two split structures of `G` must therefore be
compatible in all four cross pairs, which they never are at full size on
the graphs computed.

Reading. The doubling adds one copy and one reversed copy; both towers
suppress the monochromatic substructures relative to random graphs but
cannot suppress the mixed ones (a clique joined to an independent set, a
chain shuffled with a reversed chain), and the mixed ones are what the
Ramsey question counts. The exact laws (THM-483 for tournaments, Theorem 4
of the seventh note for graphs) say the mixed structure of level `k+1` is a
compatible pair of mixed structures of level `k`; the growth `3, 4, 6, 8` is
the growth of the largest compatible pair, and no closed form is offered.

## 2. Lehmer's simplest quintic and the Brocard squares (FINITE-EXACT; NUMEROLOGY typed)

Emma Lehmer's quintic `P_n(x) = x^5 + n^2 x^4 - (2n^3 + 6n^2 + 10n + 10) x^3 +
(n^4 + 5n^3 + 11n^2 + 15n + 5) x^2 + (n^3 + 4n^2 + 10n + 10) x + 1` has
discriminant `(n^3 + 5n^2 + 10n + 7)^2 f(n)^4` with `f(n) = n^4 + 5n^3 + 15n^2
+ 25n + 25`, and when `f(n)` is squarefree its splitting field is the cyclic
quintic field of conductor `f(n)` (the roots are units and Gaussian-period
combinations). Checked with sympy: for `n = -1, 0, 1` the polynomials `x^5 +
x^4 - 4x^3 - 3x^2 + 3x + 1`, `x^5 - 10x^3 + 5x^2 + 10x + 1`, `x^5 + x^4 - 28x^3
+ 37x^2 + 25x + 1` are irreducible with Galois group `C_5` and discriminants
`11^4`, `5^8·7^2`, `23^2·71^4`, matching the formula; the values `f(n)` for
`n = -6..6` are `631, 275, 101, 31, 11, 11, 25, 71, 191, 451, 941, 1775,
3091`, every prime among them `= 1 mod 5`.

**The coincidence (exact as a numerical fact).** The Brocard solutions
`4! + 1 = 25`, `5! + 1 = 121`, `7! + 1 = 5041` are `f(0)`, `f(-1)^2`, `f(1)^2`:
the owner's `5, 11, 71` are the conductors of Lehmer's quintics at `n = 0`
(as `5^2`), `-1`, `1`, i.e. `Q(zeta_25)^+`'s quintic subfield, the quintic
subfield of `Q(zeta_11)` and that of `Q(zeta_71)`. `11` is the smallest
conductor of any cyclic quintic field, `25` the second, `71` the sixth (`11,
25, 31, 41, 61, 71`; Lehmer's family also hits `31` at `n = -3` and misses `41,
61`). No mechanism is known to this note that would make `f(0) - 1 = 4!`,
`f(-1)^2 - 1 = 5!`, `f(1)^2 - 1 = 7!` anything but a coincidence of three small
numbers, and the classification is NUMEROLOGY; it is recorded because it is
exact and because it is presumably what the directive meant by "the Lehmer
5".

## 3. Brocard's problem: the exact structure (PROVED reading + FINITE-EXACT)

**Proposition 3.** `n! + 1 = m^2` iff `n! = (m - 1)(m + 1)`, i.e. iff the
divisor involution `sigma(x) = n!/x` on the divisors of `n!` has an orbit of
diameter two. Since `n!` is not a square for `n >= 2` (Bertrand: a prime in
`(n/2, n]` divides `n!` exactly once), `sigma` has no fixed point and the
Brocard solutions are the cases where the two divisors closest to
`sqrt(n!)` are as close as they can be: `(4, 6)`, `(10, 12)`, `(70, 72)`. In
the Kuratowski dictionary this is KT-b, "dual pair + self-dual", with the
self-dual point removed: the duality is `sigma`, the pair is `{m-1, m+1}`,
and the would-be fixed point `sqrt(n!)` is irrational. The solutions found
for `n <= 60` are exactly `(4, 5), (5, 11), (7, 71)` (conjecturally all;
`10^15` by Matson; finitely many under abc by Overholt).

*Proof.* Rewrite and factor. ∎

Nearby equations (`n <= 30`): `n! + 2 = m^2` only at `n = 2`; `n! + 3 = m^2` at
`n = 1, 3`; `n! + 4` never; `n! - 1 = m^2` only at `n = 1, 2`. The "twins +
sporadic third" shape (KT-a) also fits the arguments: `4` and `5` are
consecutive solutions and `7` is isolated; but no order relation makes `7`
a container, so that reading is ANALOGY.

**Other exact facts about the same numbers (each true, none connected).**
`4! + 1 = 5^2` is the statement that `5` is a Wilson prime (`(p-1)! = -1 mod
p^2`; the Wilson primes are `5, 13, 563`). `5! + 1 = 11^2` is `((p-1)/2)! = -1
mod p^2` for `p = 11`, the half-Wilson congruence of Mordell at level `p^2`,
which below 2000 holds only for `p = 11` and `p = 47`. `7! + 1 = 71^2` has `7 =
(71-1)/10` and no Wilson reading. The roots `5, 11, 71` all lie in Ogg's list
of supersingular primes `{2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59,
71}` (the primes dividing the Monster's order; `X_0(p)^+` of genus 0), `71`
being the largest. `4! = 24 = |2T|` and `5! = 120 = |2I|` are the orders of the
binary tetrahedral and binary icosahedral groups (McKay: `E_6` and `E_8`;
the binary octahedral group has order `48 = 2·4!`), while `7! = |S_7|`.

## 4. The fives, inventoried (CITED + typed)

| five | its exact source | in the repo | relation to the others |
|---|---|---|---|
| Platonic solids `{3,3},{3,4},{4,3},{3,5},{5,3}` | `(p-2)(q-2) < 4`, i.e. the sign of `1/2 + 1/3 + 1/n - 1` for `n = 3, 4, 5` | HYP-3771/3772; S10 section 4 | two dual pairs + one self-dual (KT-b twice, and once with a fixed point) |
| Fermat primes `3, 5, 17, 257, 65537` | Gauss–Wantzel (constructible `p`-gons); THM-871's rungs; the sea kernel (S10) | S10 section 2 | "not the same five" as the solids (S10) |
| the golden five: `Q(sqrt5) ⊂ Q(zeta_5)`, `phi`, the pentagon `P_5 = C_5` | the cyclotomic field of the pentagon; `zeta_5 + conj = 1/phi` | third note | the pentagon is the seed of the graph tower above and the Ramsey witness `R(3,3) > 5` |
| Lehmer's quintic (degree 5) | cyclic quintic fields, conductors `f(n)`; subfields of `Q(zeta_p)` for `p = 1 mod 5` | this note | its first conductors `11, 25, 71` are the Brocard roots |
| Catalan `C_3 = 5` | binary trees with three nodes; the critical spine blocks of length 7 (seventh note) | seventh note | — |
| the negative cycle `-5` and the `3n-1` cycle `5` | `9/8 = 3^2/2^3`, the convergent `3/2` of `log_2 3` | S12, fifth note | `5 = F_1` and `17 = F_2` are Fermat primes: a coincidence of two small primes, typed NUMEROLOGY |

Verdict on "deep 5 structures": each five has its own exact source and the
sources are different (angle defect, constructibility, cyclotomy at 5,
cyclotomy at primes `1 mod 5`, ballot paths); the repo's earlier finding
that the Platonic and Fermat fives are not the same five (S10) extends to
all six rows. What is exact across rows is the KT-b shape (dual pair +
self-dual) in the Platonic five and in the Brocard factor pairs, and the
cyclotomic thread `Q(zeta_5) -> Q(zeta_11), Q(zeta_71)` behind the golden
and Lehmer rows.

## 5. Directions (DIRECTION; none pursued)

* **D21.** A closed form for `cs(D(G))` from the four cross-compatible
  split structures of `G`, and hence for the graph tower's `omega` and
  `alpha` sequences; compare with THM-483's zigzag numbers of the tower
  levels (`z(T_31) = 11`).
* **D22.** The Brocard coincidence with Lehmer's conductors: is there any
  reason `f(n)^2 - 1` or `f(n) - 1` should be a factorial for small `n`?
  Compute `f(n)^2 - 1` and `f(n) - 1` for `|n| <= 50` against factorials and
  near-factorials (none expected beyond the three).
* **D23.** The half-Wilson primes at level `p^2` (`11, 47` below 2000): a
  sequence worth checking against OEIS and against the Brocard/Lehmer
  numbers; `47` is also supersingular.

## 6. Reproduction

    cd 04-computation/experiments
    python3 collatz_doubling_tower_brocard_20260927.py > collatz_doubling_tower_brocard_20260927.out

Standard library only (sympy was used once, in the scratchpad, for the
Galois groups and discriminants quoted in section 2); about ten minutes.
