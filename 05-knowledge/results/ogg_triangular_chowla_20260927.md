# Ogg's fifteen primes, the nodal cubic against triangular numbers, and the aliquot map without its unit: genus-zero rigidity as the repo's exact-invariant regime, class numbers as Paley half-row imbalances, an elliptic curve behind "triangular y", and the third carry regime (parity alternation)

**Session:** opus, `collatz-poset-dag-20260927` (S17 continuation), 2026-09-27.
**Owner's directive:** the attached paper (Valamontes, *Supersingular Primes,
Genus-Zero Rigidity, and the Limits of Distributional Explanations*, v1.0.0,
4 pages, sha256 `05d26aa3afcdd114...`), "elliptic curves like `y^2 = (x+1)x^2`
vs the triangular numbers", the pasted Hurwitz/Shimura snippet (42, 168,
504, 1092; the `(2,3,7)` volume `1/42`), "the 15 Ogg primes, supersingularity,
the number `30 = 2·3·5`", and "aliquot sequences under the condition
`F = U + S`, the number itself and `1` both discounted", read here as the
map `s'(n) = sigma(n) - n - 1` (Chowla's function iterated). Asked: themes,
crossovers, surprising connections, recurring integers, fundamental
techniques abstractly.
**Inherits (cited, not re-derived):** the barrier atlas's typing rule and the
STICKY control ([atlas](collatz_procgen_20260922_barrier_atlas.md) sections
0 and 7; [STICKY note](collatz_sticky_20260927_size_coupled_persistence.md)),
the parallel session's aliquot notes ([ninth](collatz_aliquot_lehmer_five_20260927.md),
[tenth](collatz_two_carries_typology_20260927.md): the driver lock, the
typology "provability tracks an exact invariant"), its eighth note (the
roots `5, 11, 71` of the Brocard cases lie in Ogg's list), S10's Platonic
groups ([gilbreath_fermat_platonic](gilbreath_fermat_platonic_20260926.md):
`PSL(2,3), PGL(2,3), PSL(2,5)` of orders `12, 24, 60`; `PSL(2,7)` of order
`168` as the Klein quartic's group), the Paley thread (THM-640 the Paley
bridge; THM-448 the DRT doubling tower with `Aut = F_21`; HYP-3805, the Paley
heptagon as the LRC extremal object), the owner's moonshine exercise in the
synthesis (`6 = 5 + 1 -> S_5 = PGL(2,5) -> M_12 -> M_24 -> Golay -> Leech ->
196884`), HYP-2220 (even perfect numbers as `T_(2^p - 1)`). Classical facts
from memory unless marked: Ogg 1975, Deuring 1941, Dirichlet's class number
formula, Nagell–Lutz, Siegel, Cattaneo 1951 (quasiperfect numbers are odd
squares), Hagis–Lord 1977 (quasi-amicable pairs), Lucas/Watson/Ma
(cannonball), Hurwitz, Takeuchi.

**Status: PROVED (Proposition 1 the parity theorem of the Chowla map and its
corollaries; Proposition 2 the nodal cubic's two triangular conditions as a
conic and an elliptic curve; the class-number reading of Ogg's condition is
classical and re-derived) + FINITE-EXACT (Ogg's list reproduced for `p < 120`
by the genus formula and for `5 <= p <= 97` by the field of definition of the
supersingular `j`-invariants; the Paley half-row class-number formula to
`p < 200`; the Chowla map on all `n <= 10^6`; integer points of the curve to
`X <= 10^6`; Pell and triangular searches) + CITED (Ogg, Dirichlet, moonshine)
+ NUMEROLOGY typed (378 = T_27 and friends) + DIRECTION. The paper is typed
in section 1: correct on Ogg's theorem, a valid methodological point, an
overclaim in its framing. Independent audit: section 8.** Script
`04-computation/experiments/ogg_triangular_chowla_20260927.py`, output
beside it.

## 0. The answer in one paragraph

Four exact things, one abstract theme, and a pile of typed coincidences.
(1) Ogg's fifteen primes are the `p` with `g(X_0(p)^+) = 0`, i.e. `2 g(X_0(p))
+ 2 = h(-4p) + h(-p)·[p ≡ 3 (4)]` (the fixed points of the Fricke involution
are CM points, counted by class numbers), and for `p ≡ 3 (mod 4)` Dirichlet
writes `h(-p) = (sum of (a/p) over 0 < a < p/2)/(2 - (2/p))`: the class
number is the residue imbalance of the first half-row of the Paley
tournament, the repo's THM-640 object; so the nine odd-type members of
Ogg's list (`3, 7, 11, 19, 23, 31, 47, 59, 71`) are decided by a Paley
statistic (class numbers `1, 1, 1, 1, 3, 3, 5, 3, 7`), and both routes (genus
formula; supersingular `j`'s in `F_p`, computed from the Hasse polynomial on
the Legendre line) reproduce the list exactly (section 2). (2) The numbers
`30` and `42` are the two sides of the Gauss–Bonnet formula in the pasted
snippet: `1/2 + 1/3 + 1/5 - 1 = +1/30` (spherical excess: the icosahedral
group of order `2·30`, S10's `PSL(2,5)`) and `1/2 + 1/3 + 1/7 - 1 = -1/42`
(hyperbolic deficit: Hurwitz's `42(2g - 2)`, the orders `168, 504, 1092 =
42·4, 42·12, 42·26`), with the Platonic ladder `1/6, 1/12, 1/30` and the
Euclidean `(2,3,6)` in between; "the smallest product of three distinct
primes" is the last spherical case (section 3). (3) On the nodal cubic `y^2 =
x^2(x+1)`, whose points are `(t^2 - 1, t^3 - t)` and whose group is `G_m`, "x
triangular" is a conic (`u^2 - 8t^2 = -7`, infinitely many: `x = 3, 15, 120,
528, 4095, 17955, ...`) while "y triangular" is the elliptic curve `Y^2 = X^3
- 4X + 1` (`X = 2t`, conductor `916`), of rank at least two (the audit proved
`P = (0, 1)` and `Q = (2, 1)` independent), all of whose eleven integer points
found are `±(nP + mQ)`; its solutions give
`y = 6, 120, 210, 990, 185136, 258474216` (`T_3, T_15, T_20, T_44, T_608,
T_22736`), and the only point with both coordinates triangular below `t =
10^6` is `(3, 6) = (T_2, T_3)`: genus zero versus genus one decides infinite
versus finite, the paper's rigidity theme in miniature (section 4). (4) The
aliquot map with the unit divisor removed, `s'(n) = sigma(n) - n - 1`,
**alternates parity**: `s'(even)` is odd unless the odd part is a square, and
`s'(odd)` is even unless `n` is a square (Proposition 1). So the aliquot
driver lock of the STICKY note is literally the divisor `1`: keep it and
parity persists (aliquot, sticky, conjecturally divergent); add it on the
other side and parity is fresh each step (Collatz, memoryless); remove it
and parity flips (Chowla, anti-sticky: persistence `≍ N^(-1/2)`, measured
`0.025, 0.008, 0.0025, 0.0008` at `N = 10^3..10^6`). Every start `n <= 10^6`
ends at `0` (through a prime, `98.2%`) or in a 2-cycle (`1.8%`; eighteen
betrothed pairs, all of opposite parity, as the theorem forces unless a
square-type member is present); fixed points must be squares or twice
squares (the quasiperfect problem: none known; Cattaneo's theorem, "odd
squares", is strictly stronger, the parity theorem recovers only its
parity half); the odd step's mean drift is `-2.8` bits. (5)
The abstract theme is the paper's own: rigidity (a discrete invariant, a
finite list) against distribution (density, asymptotics). That is the repo's
discipline in the Collatz thread (density-one theorems never reach the
pointwise statement; the typology's "provability tracks an exact
invariant"; STICKY), and it reappears in (3) as genus. The paper is right
about Ogg's theorem and about separating the two notions of "supersingular
prime", but its "moonshine-independent complete explanation" explains the
*list* (a genus computation Ogg already did), not the coincidence with the
Monster's primes, which is what moonshine explains; and its reference for
the list is the 1974 hyperelliptic paper rather than the 1975 note
(recollection). Recurring integers, typed: `378 = T_27` is the sum of the
fifteen primes (NUMEROLOGY); `24, 120 = 5!, 210 = 7#` are the `t = 5, 6`
points of the elliptic family (`x = 24` beside Lucas's `24` and the Leech
vector; NUMEROLOGY); `189` is a Pell index (`T_189 = 134^2 - 1`; the memory's
`189`, NUMEROLOGY); `168 = 8·21` contains the Paley heptagon's `F_21` as the
Sylow-7 normalizer of the Klein quartic's group (STRUCTURAL, classical);
genus `14` of the first Hurwitz triplet against LRC(14) (NUMEROLOGY); `118 =
2·59` (NUMEROLOGY).

## 1. The paper, typed

The paper distinguishes (i) Ogg's finite list `{2, 3, 5, 7, 11, 13, 17, 19,
23, 29, 31, 41, 47, 59, 71}`, characterized by "every supersingular elliptic
curve in characteristic `p` is defined over `F_p`" iff "`X_0(p)^+ = X_0(p)/w_p`
has genus zero", from (ii) the Lang–Trotter supersingular reduction primes of
a fixed curve over `Q` (`a_p(E) = 0`; infinitely many by Elkies; conjectured
`~ c sqrt(X)/log X`), and states a "regime-separation principle": a rigidity
classification (discrete invariant, finite, complete) cannot explain a
distributional law without an explicit bridge.

* **CITED, correct:** Theorem 1 (Ogg) as stated. Reproduced here two ways
  (section 2).
* **Valid methodological point, and the repo's own:** the separation of a
  finite classification by a discrete invariant from an asymptotic
  statistical law is exactly the Collatz thread's separation of "exact
  invariant / pointwise" from "density / distributional" (the parallel
  tenth note's typology: provability tracks an exact invariant; the barrier
  atlas: Terras, Korec, Tao are distributional and cannot reach the
  pointwise statement; STICKY: what the distributional mechanisms use is
  memorylessness at every scale). The paper's Proposition 1 ("the finite
  classification does not imply, approximate or constrain the Lang–Trotter
  asymptotic") is true and, in the repo's language, a blindness statement.
* **Overclaim in framing:** "A complete explanation of the fifteen primes
  exists within arithmetic geometry itself, without appeal to moonshine."
  Ogg's genus computation explains *why the list is what it is*; it does not
  explain why it coincides with the set of primes dividing the order of the
  Monster, which is the content of Ogg's 1975 remark (the Jack Daniels
  question) and of monstrous moonshine (Conway–Norton's genus-zero
  Hauptmoduln for `Gamma_0(p)+`, Borcherds). The paper never mentions the
  Monster, so it does not contradict this; it only calls "complete" an
  explanation of the wrong question.
* **Citation imprecision (recollection):** the list and the genus-zero
  characterization are in Ogg's 1975 Séminaire Delange–Pisot–Poitou note
  "Automorphismes de courbes modulaires", not in the 1974 Bull. SMF paper
  "Hyperelliptic modular curves" cited as [2]. Not checked against the
  sources here.
* **No new mathematics** is claimed by the paper beyond the distinction.
  Two imprecisions (audit): its Lang–Trotter sentence (`sqrt X/log X`) omits
  the non-CM hypothesis (CM curves have supersingular reduction at a set of
  primes of density `1/2`), and "the finiteness of the list follows from the
  classification of genus-zero modular curves" is a loose justification
  (finiteness is immediate from `g(X_0(p)^+) >= (p/12 - O(sqrt p log p))/2`).
  Its Theorem 1 and Proposition 1 are unaffected.

## 2. Ogg's list, reproduced twice, and the Paley reading (PROVED classical + FINITE-EXACT)

**Genus route.** For prime `p`, `g(X_0(p)) = 1 + (p+1)/12 - nu_2/4 - nu_3/3 -
1` with `nu_2 = 1 + (-4/p)`, `nu_3 = 1 + (-3/p)` (Kronecker symbols; `nu_3 = 0`
at `p = 2`, `nu_2 = 0` at `p = 3`). The Fricke involution `w_p` fixes the CM
points of discriminants `-4p` and, when `p ≡ 3 (mod 4)`, `-p`, so it has `f =
h(-4p) + h(-p)·[p ≡ 3 (4)]` fixed points (for `p = 2`: `h(-8) + h(-4) = 2`),
and Riemann–Hurwitz gives `g(X_0(p)^+) = (2g + 2 - f)/4`. With class numbers
counted as reduced primitive forms, `g(X_0(p)^+) = 0` exactly for the fifteen
primes among all `p < 120` (script, part 1; the table `p : g, f, g^+` is in
the output: `11: 1, 4, 0`; `37: 2, 2, 1`; `71: 6, 14, 0`; `73: 5, 4, 2`; `97: 7,
4, 3`). **Ogg's condition in one line (`p >= 3`; at `p = 2` the count `f =
h(-8) + h(-4) = 2` is the convention used above):** `2 g(X_0(p)) + 2 = h(-4p) +
h(-p)·[p ≡ 3 (4)]`. The audit extended the check to all `p < 1000` (genus zero
exactly at the fifteen) and confirmed the fixed-point count without class
numbers: `f = 2·#{supersingular j in F_p}` for all `5 <= p < 200`.

**Supersingular route.** On the Legendre line, `E_lambda: y^2 = x(x-1)(x-
lambda)` is supersingular iff `H_p(lambda) = sum_(k <= m) C(m,k)^2 lambda^k =
0`, `m = (p-1)/2` (Deuring; all roots lie in `F_(p^2)`), and `j = 256
(lambda^2 - lambda + 1)^3 / (lambda^2 (lambda - 1)^2)`. Enumerating `lambda in
F_(p^2)` for `5 <= p <= 97`: the number of supersingular `j`-invariants is
`1, 1, 2, 1, 2, 2, 3, 3, 3, 3, 4, 4, 5, 5, 6, 5, 6, 7, 6, 7, 8, 8, 8` for `p = 5,
7, ..., 97`, and *all* of them lie in `F_p` exactly for `p in {5, 7, 11, 13, 17,
19, 23, 29, 31, 41, 47, 59, 71}` (e.g. `p = 37`: 3 supersingular `j`'s, 1 in
`F_p`; `p = 97`: 8 and 2). This is Ogg's theorem checked on both sides (the
equivalence of the two routes is verified numerically here, not derived;
the derivation is `#ss = g + 1` and `#{ss in F_p} = f/2`).

**The Paley reading (Dirichlet, re-verified).** For `p ≡ 3 (mod 4)`, `p > 3`:
`h(-p) = (1/(2 - (2/p))) sum_(0 < a < p/2) (a/p)`, verified for all such `p <
200`. The sum is the excess of quadratic residues over non-residues in the
first half-row of the Paley tournament `T_p` (arc `i -> j` iff `j - i` is a
residue: THM-640's construction), i.e. the out-degree minus in-degree of
vertex `0` into the "first half" `{1, ..., (p-1)/2}`. So for the odd-type
primes the genus-zero condition reads: `2 g(X_0(p)) + 2 - h(-4p)` equals the
Paley half-row imbalance divided by `2 - (2/p)`. Values on Ogg's list: `h(-7)
= h(-11) = h(-19) = 1`, `h(-23) = h(-31) = h(-59) = 3`, `h(-47) = 5`, `h(-71) =
7`. **Typing:** STRUCTURAL and classical (Dirichlet's class number formula;
verified by the audit for every `p ≡ 3 (mod 4)` in `(3, 2000)`, and false at
`p = 3`); the connection to the repo is that its Paley objects (THM-640,
THM-448, and HYP-3805, whose status in the repo's investigation backlog is
REFUTED-BUT-CLOSE with the corrected formula open) carry the class number
as a first-half-row statistic. It is not a mechanism for LRC or for Collatz;
it is a shared object.

**The nodal cubic as the boundary of the same family.** The owner's `y^2 =
x^2(x+1)` is the Legendre curve `y^2 = x(x-1)(x-lambda)` at `lambda = 1` after
`x -> x + 1` (exactly; the fibre at `lambda = 0` is its twist by `-1`, a
non-split node at `p ≡ 3 (mod 4)`): the singular fibre of the family whose
supersingular locus is the zero set of `H_p`; it has split multiplicative
reduction at every odd prime (`a_p = 1`, `#E_ns(F_p) = p - 1`; at `p = 2` the
node is a cusp, `E_ns = G_a`), the opposite extreme from supersingular (`a_p
= 0`). So the two curves the owner names are the two ends of the Legendre
line.

## 3. Thirty and forty-two: one formula (PROVED, classical)

The snippet's Gauss–Bonnet/Riemann–Hurwitz identity for a triangle group
`(2, 3, q)` has area `2 pi (1 - 1/2 - 1/3 - 1/q)` up to sign. The excess `1/2 +
1/3 + 1/q - 1` is `1/6, 1/12, 1/30, 0, -1/42` for `q = 3, 4, 5, 6, 7`: the three
spherical cases give the finite rotation groups of orders `2·6 = 12`, `2·12
= 24`, `2·30 = 60` (S10's `PSL(2,3) = A_4`, `PGL(2,3) = S_4`, `PSL(2,5) =
A_5`), `q = 6` is Euclidean (the hexagonal lattice, the CM curve `y^2 = x^3 +
1` with `j = 0`, supersingular exactly at `p = 3` and at `p ≡ 2 (mod 3)`), and `q = 7` is the
first hyperbolic case with the minimal deficit `1/42`, whence Hurwitz's `#Aut
<= 42(2g - 2)` and the orders `168 = 42·4` (`g = 3`, Klein), `504 = 42·12`
(`g = 7`, Fricke–Macbeath), `1092 = 42·26` (`g = 14`, the first Hurwitz
triplet). So `30 = 2·3·5` and `42 = 2·3·7` are the last spherical and the
first hyperbolic values of the same denominator: "the smallest number with
three distinct prime factors" is the excess of the icosahedron, and its
successor `42` is the deficit of the Klein quartic. Factorizations: `168 =
2^3·3·7`, `504 = 2^3·3^2·7`, `1092 = 2^2·3·7·13`, all supported on Ogg primes,
as any `|PSL(2,q)|` with `q in {7, 8, 13}` must be (`q(q^2 - 1)/2` for odd `q`).

**Klein quartic and the Paley heptagon (STRUCTURAL, classical).** `|Aut(T_7)|
= 21` (brute force over `S_7`), the Frobenius group `F_21 = 7:3`, which is the
normalizer of a Sylow 7-subgroup in `PSL(2,7) = Aut(Klein quartic) =
Aut(Fano plane)`, of index `8`. The repo's HYP-3805 (Paley heptagon as the
LRC extremal object) and THM-448 (the doubling tower with `Aut = F_21`)
therefore sit inside the Klein quartic's symmetry; `168 = 8·21`.

## 4. The nodal cubic and the triangular numbers (PROVED small + FINITE-EXACT)

**Proposition 2.** On `C: y^2 = x^2(x + 1)`:
(a) the rational points are `(t^2 - 1, t^3 - t)`, `t in Q`; the nonsingular
points form a group isomorphic to `G_m` via `(x, y) -> (y + x)/(y - x) = (t +
1)/(t - 1)` (the tangents at the node are `y = ± x`); the integer points have
`y = (t-1) t (t+1) = 6 C(t+1, 3)`, six times a tetrahedral number.
(b) `x = t^2 - 1` is triangular iff `(2m + 1)^2 - 8 t^2 = -7`: a conic, with
infinitely many integer solutions, `t = 1, 2, 4, 11, 23, 64, 134, 373, 781,
2174, 4552, 12671, ...` (two orbits under the unit `3 + sqrt 8`), giving `x =
T_m` for `m = 0, 2, 5, 15, 32, 90, 189, 527, ...`.
(c) `y = t^3 - t` is triangular iff `8t^3 - 8t + 1` is a square iff `(X, Y) =
(2t, 2k + 1)` lies on the elliptic curve `E: Y^2 = X^3 - 4X + 1` (discriminant
`2^4·229`). `P = (0, 1)` has infinite order (`3P` has `x = -7/4`, so `P` is not
torsion by Nagell–Lutz). The integer points with `-2 <= X <= 10^6` (complete
to `10^7`, audit) are exactly `(-2, 1), (-1, 2), (0, 1), (2, 1), (3, 4), (4, 7),
(10, 31), (12, 41), (20, 89), (114, 1217), (1274, 45473)`, and with `Q = (2, 1)`
and the sign of `Y` respected they are `-P - Q, -P + Q, P, Q, -2P - Q, 2P, -2P
+ Q, -2Q, 2P + 2Q, -4P - Q, -2P + 3Q` (the session script compared `|Y|` and
listed the combinations up to sign). **`P` and `Q` are independent, so `rank
E(Q) >= 2`** (audit): `E(Q)` has trivial torsion (`gcd(#E(F_3), #E(F_5)) =
gcd(7, 9) = 1`), so a dependence could be taken as `aP + bQ = O` with `a, b`
not both even; but in `E(F_37) = Z/2 x Z/22` the reductions of `P`, `Q` and `P
+ Q` all lie outside `2E(F_37)`, which rules out each parity class of `(a,
b)` (the same at `p = 53, 173, 241`; no relation with `|a|, |b| <= 30` either).
Whether the rank is exactly two is not settled here (a descent, or a table
at conductor `916 = 2^2·229`, the conductor found by the audit's Tate
algorithm: type IV at `2`, `I_1` at `229`, minimal model). The even `X >= 2`
give `t = 1, 2, 5, 6, 10, 57, 637` and `y = 0, 6, 120, 210, 990, 185136,
258474216 = T_0, T_3, T_15, T_20, T_44, T_608, T_22736`; by Siegel the list of
`t` is finite, and the search to `t = 10^6` (audit: `5·10^6`) found no others.
(d) Both coordinates triangular, `t <= 10^6`: only `(0, 0)` and `(3, 6) = (T_2,
T_3)`.

*Proof.* (a) Put `y = t x`. (b) `t^2 - 1 = m(m+1)/2` iff `(2m+1)^2 = 8t^2 - 7`.
(c) `t^3 - t = k(k+1)/2` iff `(2k+1)^2 = 8t^3 - 8t + 1`; with `X = 2t` this is
`Y^2 = X^3 - 4X + 1`. The group law computations are in the script (exact
rationals). (d) Intersection of the lists in (b) and (c). ∎

**Reading.** The same object carries a genus-zero condition (infinitely
many, a Pell family) and a genus-one condition (finitely many integer
points): the paper's rigidity/distribution dichotomy at the level of one
Diophantine question, decided by genus. The integers that appear are the
owner's: `x = 24` at `t = 5` (with `y = 120 = 5! = T_15 = |2I|`, the binary
icosahedral order of the parallel eighth note), `x = 35`, `y = 210 = 7# =
T_20` at `t = 6`, and `24` is also Lucas's cannonball number (`1^2 + ... + n^2`
is a square only for `n = 1, 24`, re-checked to `10^6`; `70^2` is the Leech
vector `(0, 1, ..., 24; 70)` in `II_(25,1)`, the owner's moonshine chain).
These are recurrences of small integers with distinct mechanisms
(NUMEROLOGY), recorded because the owner asked for them; the structural
content of the section is (a)–(d).

**Perfect numbers.** `T_(2^p - 1) = 2^(p-1)(2^p - 1)` (HYP-2220) is the
triangular reading of Euclid–Euler; on the nodal cubic `x = T_m` with `m =
2^p - 1` needs `2^(p-1)(2^p - 1) + 1` to be a square (`p = 2`: `7`, no; `p =
3`: `29`, no; none for prime `p <= 100`, audit): the perfect numbers are not
on the Pell family, although the Euclid-shaped non-prime case `k = 4`, `T_15
= 120 = 2^3(2^4 - 1)`, is (`121 = 11^2`, the point `t = 11`). No connection.

## 5. The aliquot map without its unit (PROVED + FINITE-EXACT)

Read "`F = U + S`, the number itself and `1` discounted" as `s'(n) = sigma(n) -
n - 1`, the sum of the divisors strictly between `1` and `n` (Chowla's
function iterated; its 2-cycles are the quasi-amicable or betrothed pairs,
its fixed points the quasiperfect numbers, none known).

**Proposition 1 (parity alternation).** Let `n >= 2`.
(a) If `n` is even, `n = 2^a m` with `m` odd, then `s'(n)` is odd iff `m` is
not a perfect square.
(b) If `n` is odd, then `s'(n)` is even iff `n` is not a perfect square.
(c) Consequently: a fixed point of `s'` is a square or twice a square (this
is only the parity half of Cattaneo's theorem that quasiperfect numbers are
odd squares, which also excludes `2^a m^2` with `a` odd, since `sigma(n) = 2n
+ 1` would force `m^2 ≡ -1 (mod 2^(a+1) - 1)`, impossible modulo `4`); a cycle
of odd length contains an element whose odd part
is a square; a 2-cycle whose two members have the same parity contains
such an element (so betrothed pairs have opposite parity unless one member
is of square type); and the proportion of `n <= N` with `s'(n) ≡ n (mod 2)`
is `O(N^(-1/2))`.

*Proof.* `sigma(n)` is odd iff `n` is a square or twice a square (each odd
prime power `p^e` contributes `sigma(p^e) ≡ e + 1 (mod 2)`; the power of two
contributes an odd number). For even `n`, `s'(n) ≡ sigma(n) + 1 (mod 2)`, and
`n = 2^a m` is a square or twice a square iff `m` is a square. For odd `n`,
`s'(n) ≡ sigma(n) (mod 2)`. (c): a fixed point or an odd cycle cannot
alternate parity all the way round; the count of square-type `n <= N` is
`O(sqrt N)`. ∎

*FINITE-EXACT (script, part 4).* Parity theorem checked on all `n <=
3·10^5`. Dynamics of all `n <= 10^6` (divisor sieve to `4·10^6`, three values
beyond it factored): `982099` sequences end at `0` (the last nonzero term a
prime), `17900` enter a 2-cycle, none enters a fixed point or a longer
cycle; the longest sequence has `48` steps (`n = 948375`), the largest peak
is `6.3 n` (`n = 980100`), the mean length `9.3`. The eighteen 2-cycles
reached: `(48, 75), (140, 195), (1050, 1925), (1575, 1648), (2024, 2295),
(5775, 6128), (8892, 16587), (9504, 20735), (62744, 75495), (186615, 206504),
(196664, 219975), (199760, 309135), (266000, 507759), (312620, 549219),
(526575, 544784), (573560, 817479), (587460, 1057595), (1139144, 1159095)`,
every one of opposite parity, none containing a square-type member (the
audit found nine more pairs with smaller member in `(10^6, 2·10^6]`, none
reached from starts below `10^6`, all of opposite parity). Mean `log_2(s'(n)
/n)`: `-0.048` on even `n` (the aliquot value: `s' = s - 1`), `-2.82` on odd
composite `n` (against `-5.26` for `s` on odd `n <= 10^6`, audit, where the
primes give `s = 1`; the STICKY note's `-5.34` was on `n <= 2·10^6`); `210` odd
`n <= 10^5` have `s'(n) > n` (the odd abundant numbers, `945` first). Parity
persistence in `[M, 2M)`: `0.0254, 0.0079, 0.0025, 0.0008` for `M = 10^3, ...,
10^6`; by the theorem the persisting `n` are exactly the square-type numbers
`2^a k^2` (`k` odd), `22, 71, 224, 707` of them in those ranges, so the share
is `(2M)^(-1/2)(1 + o(1))` among all `n` and the measured proportions are
forced, not an independent statistic.

**The three regimes of one carry (the crossover with STICKY).** The
aliquot map's driver lock (STICKY note, Theorem 1) is the statement that
`sigma(2^a m)` has the parity of `sigma(m)`, even for non-squares, so `s(n) =
sigma(n) - n` keeps the parity of `n`; that parity-preservation is exactly
the contribution of the divisor `1` to `sigma`. The three maps built on the
same multiplicative engine differ only in what they do with a unit:

| map | carry | parity of the next term | memory of the driving quantity | conjectured fate |
|---|---|---|---|---|
| aliquot `s(n) = sigma(n) - n` | keeps the divisor `1` | persists (unless square type) | sticky, size-coupled (`N^(-1/2)` loss at `a = 1`) | many diverge (Guy–Selfridge) |
| Collatz `(3n + 1)/2^v` | adds `+1` | fresh: `v` geometric at every scale | memoryless (exactly `1/2`) | all terminate |
| Chowla `s'(n) = sigma(n) - n - 1` | removes the divisor `1` | flips (unless square type) | anti-sticky: persistence `O(N^(-1/2))` | all terminate or 2-cycle (all `n <= 10^6`) |

So the unit divisor is to the aliquot map what the `+1` is to Collatz: the
carry that sets the parity regime, with opposite effects (persist / fresh /
alternate). Removing it turns the sticky, conjecturally divergent map into
an alternating, empirically convergent one whose only growth steps are the
even steps of drift `-0.05` (rarely positive) and the odd abundant numbers
(density about `0.002`): every second step is a collapse by `-2.8` bits on
average. **Typing for the atlas:** the Chowla map is not a barrier (its
conjectured fate agrees with Collatz's), it is the *complement control* of
STICKY: a map on which the aliquot mechanism has been switched off by
deleting one divisor, and on which the residue-averaging conclusions
(Terras-type descent on a set of density 1) hold visibly. No theorem of
termination is claimed; odd abundant numbers and even growth steps can in
principle alternate for ever, and the parity theorem alone does not exclude
it.

## 6. Recurring integers, typed

| integer | where it appears here | where it appears in the repo | type |
|---|---|---|---|
| `378 = T_27 = 2·3^3·7` | sum of Ogg's fifteen primes | `27` the Collatz record; `189 = 378/2` in the memory's integer list | NUMEROLOGY (the sum of a rigidity list is not an invariant of anything) |
| `637 = 7^2·13` | sum of the Monster's prime factors with multiplicity | — | NUMEROLOGY |
| `30`, `42` | excess of `(2,3,5)`, deficit of `(2,3,7)` | S10's Platonic groups; the Hurwitz snippet | STRUCTURAL (one formula) |
| `168 = 8·21` | `|PSL(2,7)|`, Klein quartic | `F_21 = Aut(Paley heptagon)`, HYP-3805, THM-448; the memory's `168` | STRUCTURAL (Sylow-7 normalizer) |
| `13`, `1092` | `PSL(2,13)`, the Hurwitz triplet | `13` in the memory's list; the `-17` cycle's `(7, 11)` reading | NUMEROLOGY |
| `14` | genus of the Hurwitz triplet | LRC(14) | NUMEROLOGY |
| `118 = 2·59` | genus of the `PSL(2,27)` Hurwitz curve (`9828 = 84·117`), the snippet's "next arithmetic" one; genus `17` (group of order `1344`) lies between, and every Hurwitz curve is a quotient of the arithmetic `(2,3,7)` group, the snippet's "non-arithmetic" meaning "not a congruence quotient" | `59` is an Ogg prime | NUMEROLOGY |
| `24, 120, 210` | `(x, y)` at `t = 5, 6` on the nodal cubic; cannonball `24`; `5!`, `7#` | the parallel eighth note's `|2I| = 120`, `4!`, `5!`; the Leech chain | NUMEROLOGY (distinct mechanisms) |
| `189` | Pell index: `T_189 = 134^2 - 1` | memory's `189 = 3^3·7` (a small-number coincidence there too) | NUMEROLOGY |
| `5, 11, 71` | Ogg primes | the parallel eighth note's Brocard roots | NUMEROLOGY (recorded there) |
| `3, 7, 31` | Mersenne primes in Ogg's list | THM-448's tower orders `7, 31` | NUMEROLOGY (`127` is not an Ogg prime) |

## 7. Directions (DIRECTION; none pursued beyond the probes)

* **Class numbers as tournament statistics.** The repo's Paley theorems
  never mention `h(-p)`; Dirichlet's half-row formula makes it a first-row
  statistic of `T_p`. A cheap check with content: whether the doubly regular
  tournaments of THM-448's tower (non-Paley beyond order 7) have a
  half-row imbalance with an arithmetic meaning; and whether LRC's Paley
  extremality (HYP-3805) sees `h(-7) = 1`. No claim.
* **The Chowla map as a solvable cousin.** With parity alternating, a
  termination theorem would need to control the odd abundant steps; the
  density of odd abundant numbers (`~ 0.002`) against the even step's drift
  is a Markov model like STICKY's Proposition B, with no size-coupling
  (the square-type exceptions decay like `N^(-1/2)`). If the model's
  conclusion is "all terminate", the Chowla map is the first member of the
  aliquot family in the typology's convergent corner, and a proof attempt
  would face the same pointwise wall as Collatz. Cheapest test: the
  longest sequences below `10^8`.
* **The rank of `E: Y^2 = X^3 - 4X + 1`.** Rank at least two is proved
  (section 4); a descent or a table lookup at conductor `916` would settle
  whether it is exactly two, and Siegel's finite list of `t` with `t^3 - t`
  triangular (the search to `5·10^6` gives seven).
* **What does not transfer.** Genus-zero rigidity of `X_0(p)^+` and the
  Collatz problem share only the abstract regime distinction; there is no map
  from modular curves to the Syracuse map here, and none is claimed.

## 8. Independent audit (2026-09-27)

Auditor subagent, independent re-derivation of every statement before
comparison, then an own exact-arithmetic script (930 lines, 19 s):
`04-computation/experiments/ogg_triangular_chowla_20260927_audit.py`,
output `ogg_triangular_chowla_20260927_audit.out`, report
`ogg_triangular_chowla_20260927_audit.md`. **Verdict: SOUND WITH
CORRECTIONS, all applied above.** CONFIRMED: Proposition 1 (exhaustive to
`2·10^6`; the persisting set is exactly the square-type numbers); the genus
formula, the fixed-point count and `g(X_0(p)^+) = 0` exactly at the fifteen
for all `p < 1000` (class numbers cross-checked against Dirichlet's analytic
formula on all 305 fundamental discriminants in `(-1000, 0)`); the
supersingular counts by three methods (point counting, Hasse invariant,
Legendre roots) with the Eichler–Deuring mass formula; Dirichlet's half-row
formula for all `p ≡ 3 (mod 4)` in `(3, 2000)`; Proposition 2 with own
chord-tangent code, the Pell orbits, the integer points complete to `X =
10^7`, the `t`-list to `5·10^6`, `(3, 6)` alone; the Chowla dynamics with
identical endings, cycles, basin counts, longest sequence, peak and drifts;
`|Aut(T_7)| = 21` and the Sylow-7 normalizer of an explicitly constructed
`PSL(2,7)`; every integer of sections 2, 3, 6; the reading of the paper.
CORRECTED: (C1) the eleven integer points are `±(nP + mQ)`, the session
script compared `|Y|`; (C2) the rank is at least two, proved by the audit
through `E(F_37) ≅ Z/2 × Z/22` and trivial torsion, and relations were not
tested by the session script; (C3) the nodal cubic is the Legendre fibre at
`lambda = 1` after `x -> x + 1`, the `lambda = 0` fibre being its `-1` twist;
(C4) the node is a cusp at `p = 2`; (C5) Proposition 1(c) recovers only the
parity half of Cattaneo's theorem; (C6) the `118` row: genus `17` lies
between and every Hurwitz curve is arithmetic; (C7) two imprecisions of the
paper noted; (C8) the one-line Ogg condition needs `p >= 3`; (C9) "even `X >=
2`"; (C10) the aliquot odd drift compared on one range (`-5.26` at `10^6`) and
no longer a hard-coded string in the script; (C11) the persisting set is the
square-type numbers with share `(2M)^(-1/2)`, so the measured proportions
are forced; (C12) `j = 0` is also supersingular at `p = 3`; (C13) `T_15 =
120` is a Euclid-shaped Pell point; (C14) HYP-3805's status carried. ADDED
by the audit: conductor `916`; `f = 2·#{supersingular j in F_p}` for `5 <= p <
200`; nine further betrothed pairs in `(10^6, 2·10^6]`. UNVERIFIED
(recollections, consistent between author and auditor): Ogg 1975 versus
1974 as the home of the observation; Cattaneo 1951; Hagis–Lord 1977; the
genus-17 Hurwitz group of order `1344`. OPEN: the exact rank of `E`; the
completeness of its integer points beyond `10^7`; termination of the Chowla
map. The mechanisms of C1, C3 and C10 are logged in MISTAKES.

## 9. Reproduction

```text
python 04-computation/experiments/ogg_triangular_chowla_20260927.py > 05-knowledge/results/ogg_triangular_chowla_20260927.out
```

Needs numpy and sympy; about six minutes (the Chowla dynamics of `10^6`
starts and the supersingular `lambda` search over `F_(p^2)` dominate).
