# Problems that live on the interaction of 2 and 3: the compiled lists, and fifty more

Session klein-2026-11-16-two3. The owner asked other sessions to compile this question; the compiled lists in the
repository are the **20 lanes** of `two_three_problems_census_20261113/README.md`, the **23 open problems** of
`session335_three_slot_forests_pythagorean_addresses_20261114/inventory/TWO_THREE_OPEN_PROBLEMS.md`, the three-role
synthesis `two_three_axis_20261019_numbers_fixed_points_relations.md` (THM-4720, HYP-9339) and the archive
`chanpak_two_three_archive_20261115` (which itself says "not a new list of 43 independent conjectures": 20 + 23 = 43).
Part A repeats the 43 in one line each (cited, not re-derived). Part B adds fifty items not on those lists, each tagged
`[C]` classical theorem, `[O]` classical open problem, `[R]` repository result not yet on the lists, `[Q]` a question
posed here, with the mechanism it exhibits (Θ1–Θ6 of SYNTHESIS.md). Small-number coincidences are excluded
(MISTAKE-559) unless the mechanism is named.

## A. The compiled 43 (census lanes 1–20, inventory 1–23)

Census: 1 Collatz; 2 Chan–Pak (3 = H(C₃), 2^α3^β frames); 3 lonely runner (HYP-9174); 4 Mersenne (7 = 2³ − 1);
5 Pillai/Catalan; 6 Berggren (THM-4807/4808); 7 Zeckendorf; 8 Mahler 3/2; 9 {(3/2)^k}; 10 Wieferich (THM-4644);
11 Hadwiger–Nelson (7 = Φ₃(2)); 12 3-smooth representations; 13 ideal Waring g(k); 14 S-units {2, 3}; 15 Furstenberg
×2×3; 16 abc for 2 and 3; 17 Gelfond digit sums; 18 Erdős ternary digits of 2^n; 19 Fermat primes (THM-4847);
20 ord_p(2) odd, density 7/24 (Hasse).
Inventory: 1–2 Collatz plus/minus sheets; 3 Mahler Z-numbers; 4 Erdős (2⁰, 2², 2⁸); 5 Lagarias 3-adic exceptional set;
6 Furstenberg measure rigidity; 7 Artin for 2 and 3; 8 ×2×3 lonely spectrum; 9 Collatz Fourier gate; 10 relaxed Collatz
reachability; 11 Chan–Pak; 12 two-hub Chan–Pak; 13 Goldbach/Lemoine same pair; 14 Erdős–Straus; 15 aliquot 3-cycles;
16 perfect cuboid; 17 Beal; 18 lonely-runner tightness; 19 AMM 12592 (HYP-9129); 20 Markov unicity; 21 Fermat–Catalan;
22 twin-prime Goldbach (mod-6 wheel); 23 prime hypotenuses on the (t + 1, t) branch.

## B. Fifty more, by mechanism

### Θ1 — multiplicative independence of 2 and 3 (`log 2 / log 3` irrational)

1. `[C]` **Cobham's theorem** (1969): a sequence automatic in base 2 and in base 3 is eventually periodic. The
   independence of the two bases is the whole content.
2. `[C]` **Furstenberg's orbit theorem** (1967): for irrational `x` the `×2, ×3` orbit on the circle is dense (proved),
   against the open measure-rigidity conjecture already listed; Rudolph (1990) under positive entropy.
3. `[C]` **Normality is base-dependent** (Cassels 1959, Schmidt 1960): numbers normal in base 2 but not in base 3
   exist; Schmidt's criterion is exactly `log 2 / log 3 ∉ Q`.
4. `[C]` **Collatz cycle lengths from the convergents of `log₂ 3`** (Eliahou 1993): any nontrivial cycle has length
   `301994a + 17087915b + 85137581c`; the convergents `3/2, 8/5, 19/12, 65/41, 84/53, 485/306, 1054/665, …` are the
   Pythagorean comma (`3^12/2^19`), the 41- and 53-tone temperaments (`2^84 ≈ 3^53`). `[Q]` the repository's 27→53→105 ray
   and the boundary `{13, 27, 41}` of the forbidden pair contain the denominators 41 and 53: coincidence or carry?
5. `[C]` **Baker–Tijdeman**: `|2^a − 3^b| > 2^a / a^C` effectively; this is what makes Pillai's equation decidable for
   the pair `(2, 3)` while the general Pillai conjecture stays open.
6. `[C]` **Størmer for `{2, 3}`**: the consecutive 3-smooth pairs are `(1,2), (2,3), (3,4), (8,9)` only; Catalan's `8, 9`
   is the last. (Known to HYP-9339 via the balanced cycles; listed here as the general form.)
7. `[C]` **Cantor versus Sierpiński**: the Cantor set (base-3 digits in `{0, 2}`) has dimension `log 2 / log 3`; the
   Sierpiński gasket (Pascal's triangle mod 2, three-peg Hanoi) has dimension `log 3 / log 2`. The two canonical fractal
   dimensions of the 2–3 world are reciprocal.
8. `[C]` **Behrend-type 3-AP-free sets**: integers whose base-3 digits lie in `{0, 1}` contain no three-term progression
   (no carries), giving `N^(log 2 / log 3)` elements; the same digit set is Kummer's criterion for `3 ∤ C(2n, n)` (item 14).
9. `[C]` **Josephus with skip 3** (Odlyzko–Wilf 1991): skip 2 has the binary-rotation closed form `J(n) = 2l + 1`; skip 3
   is governed by the iteration `a ↦ ⌈3a/2⌉` and a constant `K ≈ 1.6222` of unknown nature — Mahler's `3/2` again.
10. `[C]` **Radix economy**: `r / ln r` is minimised at `r = e`, so base 3 (`2.731`) beats base 2 (`2.885`) and is the most
    economical integer base; balanced ternary is the integer base closest to `e`.
11. `[C]` **The `2^n ≡ 3 (mod n)` problem**: the least solution is `n = 4700063497` (Lehmer's question); the first base-2
    pseudoprime is `341 = 11·31`, the first base-3 pseudoprime `91 = 7·13`.
12. `[Q]` **Three-peg Hanoi as a 2–3 automaton**: `2^n − 1` moves through `3^n` states; the state graph is the Sierpiński
    gasket of item 7. Is the Hanoi graph's Hamiltonian-path count a span-bounded ribbon in the sense of THM-4880?

### Θ2 — dual primitive roots, lifting the exponent, and Zsigmondy

13. `[C]` **2 generates `(Z/3^k)^×` and 3 generates `(Z/2^k)^×/{±1}`** for every `k`; equivalently `ord_{3^k}(2) = 2·3^(k−1)`
    and `ord_{2^k}(3) = 2^(k−2)`. The repository's `(Z/8)^× = ⟨−1, 3⟩` (mod-8 synthesis) is the `k = 3` case.
14. `[C]` **Lifting the exponent**: `v₃(2^n − 1) = v₃(n) + 1` for even `n`, `v₂(3^n − 1) = v₂(n) + 2` for even `n` and `= 1`
    for odd `n`. These are the carries of LANDING.md's Collatz construction (`41·2^e ≡ −1 (mod 3^m)`) stated as a law.
15. `[C]` **Kummer's theorem, asymmetric**: `C(2n, n)` is even for all `n ≥ 1`, but `3 ∤ C(2n, n)` exactly when every
    base-3 digit of `n` is 0 or 1 — the fractal exception set of item 8. `[O]` Graham's \$1000 problem: are there
    infinitely many `n` with `gcd(C(2n, n), 105) = 1`? (Known for 15, Erdős–Graham–Ruzsa–Straus.)
16. `[C]` **Pisano periods**: `π(2^k) = 3·2^(k−1)`, `π(3^k) = 8·3^(k−1)`: the Fibonacci period modulo a power of 2 carries a
    3 and modulo a power of 3 carries `8 = 2³`; `π(6) = 24` and `π(24) = 24` (THM-486).
17. `[C]` **Jacobsthal numbers** `J_n = (2^n − (−1)^n)/3`: `0, 1, 1, 3, 5, 11, 21, 43, 85, 171, 341`; the forbidden
    Chan–Pak value `21 = J_6 = 63/3` and Sarrus's pseudoprime `341 = J_10` are on this line. `[Q]` is `J_6` the only
    Jacobsthal number that is not a Hamiltonian-path count?
18. `[C]` **3 divides exactly one of `2^n − 1`, `2^n + 1`**; in general a prime serves both signs of `k·2^n ± 1` iff
    `ord_p(2)` is even (WHEELS.md Lemma D). The Mersenne numbers `2^p − 1` and Fermat numbers `2^(2^n) + 1` are never
    divisible by 3; Pépin's test uses base 3 for that reason (THM-4800 context).
19. `[R]` **Zsigmondy wheels** (WHEELS.md): the covering periods for `k·2^n ± 1` are governed by the odd part of `N`; the
    least period 24 exists because `2^6 − 1 = 63` is Zsigmondy's exception; the classical covering sets are the full
    wheels `W(24)`, `W(36) ∖ {109}`, `W(64)`; base 3 lacks modulus 2 instead of 6.
20. `[C]` **Perfect numbers**: every even perfect number other than 6 is `≡ 1 (mod 9)`; the Catalan–Mersenne chain
    `3 → 7 → 127 → 2^127 − 1` (next term unknown); Ramanujan–Nagell `2^n − 7 = x²` has solutions `n = 3, 4, 5, 7, 15`.
21. `[O]` **The New Mersenne conjecture** (Bateman–Selfridge–Wagstaff): for `p = 2^k ± 1` or `4^k ± 3`, `2^p − 1` prime
    iff `(2^p + 1)/3` prime, with the Wagstaff numbers `(2^p + 1)/3` dividing by 3 by item 18.
22. `[Q]` **Base-3 wheels and the Chan–Pak ternary lanes**: THM-4950/4952 show base-3 ribbon words miss exactly 10 and 31.
    Is there a covering-type explanation (a base-3 wheel on the ternary lane) for the missed residues?

### Θ3 — the modular group `Z/2 * Z/3` and the number 24

23. `[C]` **`PSL(2, Z) ≅ Z/2 * Z/3`** and `χ(SL(2, Z)) = ζ(−1) = −1/12`: the `q^(1/24)` of Dedekind's `η`, `Δ = η^24`,
    the weight-12 cusp form and the rank 24 of the Niemeier lattices all descend from the free product of the orders 2 and 3
    (the repository's Bass–Serre 2–3 tree, THM-3145, is this group acting on its tree).
24. `[C]` **`Q(ζ_24) = Q(i, √2, √3)`** with Galois group `(Z/2)³`: 24 is the largest `n` with `x² ≡ 1 (mod n)` for all units
    (THM-484), equivalently the largest `n` whose cyclotomic field is multiquadratic; its seven quadratic subfields
    `Q(√−1), Q(√±2), Q(√±3), Q(√±6)` are the points of a Fano plane.
25. `[C]` **von Staudt–Clausen for `B_24`**: the primes with `(p − 1) | 24` are `2, 3, 5, 7, 13`; they are the prime Coxeter
    numbers of the pure-A Niemeier lattices and the small primes of every period-24 covering wheel (WHEELS.md §4).
26. `[C]` **The binary tetrahedral group `SL(2, F_3)` has order 24** = the 24 Hurwitz units = the vertices of the 24-cell;
    `GL(2, F_3)` has order 48. `[Q]` the coverable wheel periods 24 and 48 are these orders: numerology, or is there a
    group acting on the two covering patterns of `W(24)`?
27. `[C]` **2 and 3 are the only primes with `SL(2, F_p)` solvable** (`PSL(2, F_2) ≅ S_3`, `PSL(2, F_3) ≅ A_4`); this is
    why Wiles could start at the prime 3 (Langlands–Tunnell) and switch to 5, and why 2 and 3 are the hard residual
    primes in modularity lifting.
28. `[C]` **The two Golay codes glue the two extreme Niemeier lattices**: the binary Golay code `[24, 12, 8]` glues
    `A_1^24` (Coxeter 2) and the ternary Golay code `[12, 6, 6]` glues `A_2^12` (Coxeter 3); `M_24 ⊃ M_12`.
29. `[C]` **`j = 1728 = 2^6·3^3` and `Δ = −16(4a³ + 27b²)`**: the constants of elliptic curves are 2–3 numbers; 2-descent is
    routine, 3-descent hard; 2-torsion of class groups is Gauss's genus theory, 3-torsion is Scholz reflection and
    Davenport–Heilbronn.
30. `[C]` **Monstrous moonshine's first coefficient** `196884 = 196560 + 324` with `324 = 2²·3⁴ = 300 + 24`: the Leech
    minimal-vector count plus the symmetric square of the 24 oscillators (`dim Sym²(R^24) = 300`) plus the 24 themselves.
31. `[C]` **Partition congruences fail exactly at 2 and 3** (Ahlgren–Boylan 2003: Ramanujan-type congruences exist only
    for 5, 7, 11); `[O]` the parity of `p(n)` (Parkin–Dyson) and its distribution mod 3 are open.
32. `[C]` **Fermat–Catalan at exponent 3 versus 2**: `x² + y² = z²` has infinitely many primitive solutions, `x³ + y³ = z³`
    none (Euler); the modular proof of Fermat starts at the 3-adic representation (item 27).

### Θ4 — "two is easy, three is hard" (complexity, dynamics, geometry)

33. `[C]` 2-SAT, 2-colouring, bipartite matching and 2-dimensional matching are polynomial; 3-SAT, 3-colouring and
    3-dimensional matching are NP-complete (Karp 1972).
34. `[C]` **Hilbert's third problem**: equidecomposability of polygons holds (Wallace–Bolyai–Gerwien); for polyhedra it
    fails (Dehn 1900), the Dehn invariant being the obstruction.
35. `[C]` **Two bodies integrable, three chaotic** (Poincaré, Bruns); `[C]` **Sharkovskii**: period 3 forces every
    period, powers of 2 sit at the bottom of the order (Li–Yorke "period three implies chaos"); the logistic map's
    period-3 window opens at `r = 1 + √8`.
36. `[C]` **Frobenius–Hurwitz**: real division algebras exist in dimensions 1, 2, 4, 8 only; there is no
    three-dimensional one.
37. `[O]` **Turán's `(3, 4)`-problem** (3-uniform hypergraphs; conjectured density 5/9) against the solved graph case;
    `[O]` **Tutte's 3-flow conjecture** (6-edge-connected suffices, Lovász–Thomassen–Wu–Zhang) against the trivial
    2-flows (Eulerian); `[O]` **Berge–Fulkerson**: six perfect matchings covering every edge of a bridgeless cubic graph
    exactly twice.
38. `[C]` **Ternary problems yield to the circle method, binary ones do not**: three primes (Vinogradov, Helfgott 2013)
    against binary Goldbach; three squares (Legendre) and three cubes (`33`, `42`) against two; the minor arcs need two
    `L²` copies and a third variable to save.
39. `[C]` **Steiner triple systems** exist iff `v ≡ 1, 3 (mod 6)` (Kirkman 1847); the Fano plane `STS(7)` and the
    Hamming codes `[7, 4, 3]` (binary) and `[13, 10, 3]` (ternary) have lengths `(q³ − 1)/(q − 1)` for `q = 2, 3`.
40. `[C]` **Lieb's square ice**: the residual entropy of 3-colourings of the square lattice is `(4/3)^(3/2)` per vertex
    (the six-vertex model, 6 = 2·3); the 2-state Ising and 3-state Potts models have different critical universality.
41. `[C]` **Jumping champions**: 6 is the most frequent prime gap up to about `1.7·10^35`, then 30 (Odlyzko–Rubinstein–
    Wolf); the repository's base-six prime-gap clocks (six_tutte_collatz) are the local form.
42. `[C]` **Erdős discrepancy** (Tao 2015): the extremal low-discrepancy `±1` sequences are built from the Dirichlet
    character modulo 3 with free choices at the powers of 2; `[Q]` is the Chan–Pak signed-descent union (COUNTS.md's
    3-edge-colouring without a residue colouring) a discrepancy phenomenon of the same type?

### Θ5 — additive 2 against multiplicative 3 in the tournament counts

43. `[R]` **The triangle sandwich** `1 + 2τ ≤ H ≤ 3^τ` (triangular_backbone_chanpak_20261116/CORE.md): each directed
    triangle adds at least 2 and multiplies by at most 3; the first positive induced repair costs `2⌈log₃ H⌉`.
44. `[R]` **The sharp per-letter growth constants** for span-3 ribbons, `1, 9, 17, 19, 33, 33, 39, 31` (SHARP_EIGHT.md,
    THM-5253): for a single beaten offset `a` the constant is `2^(a+2) + 1` (`9, 17, 33`), from the deletion count
    `1 + 2^(a+d)` with window out-degree `d = 2` and SHARP_EIGHT's witnesses; the maximum `39 = 5·2³ − 1` is the gadget
    `TT_6` plus a vertex beating ranks 1 and 3 (THM-5243 (2′)).
45. `[R]` **Doubling against tripling inverses**: `2d − 1`, `(4d − 3)/3`, `(4d + 1)/3` (six_tutte_collatz COUNTS.md,
    THM-5262, THM-5272: the commutator of `3h` and `2h − 1` is 2 and the rays meet only at `9 = 2³ + 1 = 3²`); the forbidden
    21 sits between its realised ternary preimage 27 and image 63 (THM-5253 T2).
46. `[Q]` **The sealed-gadget constant in general span**: for span 3 the maximal one-vertex growth factor is
    `39 = 5·2³ − 1`, with the exponent 3 the span-3 rank bound (THM-5243 (2′)); is the maximal factor of span-`w` ribbons
    always of the form `c_w·2^w − 1`? (SHARP_GROWTH gives only the crude `1 + w·2^(2w−1)`.)
47. `[Q]` **Which Jacobsthal, Mersenne and `2^a3^b ± 1` values are Hamiltonian-path counts?** Session 333 proved every
    `1 + 2^A3^B` with `A ≥ 3` is a strong count; THM-1370 excludes `7 = 2³ − 1` and `21 = J_6`; is `2^a − 1` a count for
    every `a ≥ 4` (15, 31, 63 are)?

### Θ6 — the dual ramification of 2 and 3 (`Z[i]` against `Z[ω]`)

48. `[C]` In `Z[i]` the prime 2 ramifies (`2 = −i(1 + i)²`) and 3 is inert; in `Z[ω]` the prime 3 ramifies
    (`3 = −ω²(1 − ω)²`) and 2 is inert. Sums of two squares are governed by `p mod 4`, Loeschian numbers `x² + xy + y²`
    by `p mod 3`; the square and hexagonal lattices are the two sides. The repository's Eisenstein grids (THM-4825/4826,
    the split primes 13 and 37) and `Z[√−2]` joins (THM-5242, MSG-3385) are instances.
49. `[C]` **`Φ₆(x) = x² − x + 1 = Φ₃(−x)`**: the sixth cyclotomic polynomial couples the sign (2) with the cube root of unity
    (3); the LRC(14) "Φ₆ resonance" and the Eisenstein hexagonal lattice are the same object.
50. `[C]` **D₄ triality**: the only Dynkin diagram with an `S_3` symmetry, among diagrams whose symmetries are otherwise
    `Z/2`; the trefoil `T(2, 3)` is detected by Fox 3-colourings, and its group `⟨a, b | a² = b³⟩` is the braid group `B_3`
    modulo its centre, again `Z/2 * Z/3` up to centre (item 23).

Count: 50 items in Part B against 43 in Part A. Items marked `[Q]` (4, 12, 17, 22, 26, 42, 46, 47) are questions posed
here; none is claimed as a theorem.
