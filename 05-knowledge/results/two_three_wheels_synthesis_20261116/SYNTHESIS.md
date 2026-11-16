# Six mechanisms behind the 2–3 problems, and what they prove here

Session klein-2026-11-16-two3. LIST.md sorts ninety-three problems by these six mechanisms; this note states each
mechanism as a theory with a proof obligation, says which of the owner's claims it settles, and types the rest.

## Θ1. Incommensurability: `log 2 / log 3 ∉ Q` (indeed transcendental)

*Theory.* Any statement comparing powers of 2 with powers of 3 reduces to the Diophantine approximation of
`log₂ 3 = [1; 1, 1, 2, 2, 3, 1, 5, 2, 23, …]`. Effective versions (Baker) decide fixed-base Pillai; ineffective density
versions (Furstenberg, Cobham, Schmidt) say the two bases do not see each other. The convergent denominators
`1, 2, 5, 12, 41, 53, 306, 665, …` are where `2^a` and `3^b` nearly meet.
*What it proves here.* Nothing new; it classifies: the Collatz cycle bounds (THM-4633, Eliahou), the Mahler `3/2`
problems and the Josephus-3 constant are one family. The repository's 41 and 53 (the boundary `{13, 27, 41}` of
`{7, 21}`; the ray 27 → 53 → 105) coincide with two convergent denominators; no carry between the two has been found, so
this is typed COINCIDENCE pending a mechanism.

## Θ2. Dual primitive roots, lifting the exponent, Zsigmondy

*Theory.* 2 generates the units modulo every `3^k`; 3 generates the units modulo every `2^k` up to sign;
`v₃(2^n − 1) = v₃(n) + 1` and `v₂(3^n − 1) = v₂(n) + 2` for even `n`; `2^d − 1` has a primitive prime divisor for every
`d ≥ 2` except `d = 6`. Everything about "which prime divides `k·2^n ± 1` when" is a corollary.
*What it proves here.* WHEELS.md Theorems A–D: the dyadic levels of a covering wheel are supplied by the Fermat numbers;
covering a period `N = 2^a m` reduces to covering the odd part `m` by the odd parts of the orders (sufficient, and
necessary in density); for a prime odd part `q` the period `2^a q` covers iff at least `q` primes have order with odd part
`q`; hence the least covering period is 24 for `q = 3` **because `ν(6) = 0`**, 80 for `q = 5`, 112 for `q = 7`; the classical
covering sets are the full wheels `W(24)`, `W(36) ∖ {109}`, `W(64)`; a prime serves both signs iff its order is even.
The base-3 wheels lack the modulus 2 instead of 6, and their least period is 48. The exact table to 150 confirms the
nested criterion at every period.

## Θ3. The modular group `Z/2 * Z/3`

*Theory.* `PSL(2, Z)` is the free product of a 2 and a 3; its Euler characteristic `−1/12 = ζ(−1)` is the `1/24` of
Dedekind's `η`, hence the weight 12, the 24 of `Δ = η^24`, the rank of the Niemeier lattices and the von Staudt–Clausen
primes `{2, 3, 5, 7, 13}` of `B_24`. Separately, `(Z/24)^×` has exponent 2, so `Q(ζ_24) = Q(i, √2, √3)` is the largest
multiquadratic cyclotomic field (THM-484), and `SL(2, F_p)` is solvable exactly for `p = 2, 3`.
*What it proves here.* The owner's question "how the Riesel wheel relates to `{4, 2, 5, 7, 7, 13}` and the Leech lattice"
has a precise answer: the primes `{3, 5, 7, 13}` common to the Riesel wheel, the prime Coxeter numbers of the pure-A
Niemeier lattices and the Leech shell count `196560 = 2^4·3^3·5·7·13` are the odd primes with `(p − 1) | 24`, forced into
every period-24 wheel because `(p − 1) | 24 ⇒ ord_p(2) | 24`. That is the whole relation; the covering structure itself
(the two patterns, the roles of 17 and 241, the failure of 12) comes from Θ2 and is invisible to the lattice. The literal
multiset `{4, 2, 5, 7, 7, 13}` carries no mechanism (NUMEROLOGY, MISTAKE-559); the earlier `{4, 2, 5, 3, 7, 13}` is
`{d + 1 : d | 12}` (THM-4838). Two different 24s meet here — `lcm(ord_p 2)` of a wheel and `2/|χ(PSL₂Z)|` — through the
divisor lattice of 24 alone. New constructions from the lattice side are not expected and none is proposed.

## Θ4. "Two is easy, three is hard"

*Theory.* Across complexity (2-SAT/3-SAT, 2-/3-colouring, 2D/3D matching), geometry (Dehn), dynamics (two-body/three-
body, Sharkovskii), algebra (no three-dimensional division algebra) and additive number theory (ternary problems yield
to the circle method, binary ones do not), the passage from two to three objects is where the first genuine
obstruction appears. In the repository the same line separates the realisable two-port readers from the open
three-way selection (Chan–Pak, two-hub Chan–Pak, binary Goldbach-type "same pair" problems).
*What it proves here.* It is a classifier, not a theorem. It predicts which of the owner's claims are provable by direct
construction (two-port, binary-word programmable readers: yes) and which need a third ingredient (an exact target from
an arbitrary odd count: no construction yet). The variable-reader packet of 22:32 (below) is exactly a two-object
result: two synchronized readers, one binary address.

## Θ5. Additive 2, multiplicative 3 (the tournament sandwich)

*Theory.* In Hamiltonian-path counts each directed triangle adds at least 2 and multiplies by at most 3:
`1 + 2τ ≤ H ≤ 3^τ` (CORE.md); the child rule `3H + 2σ`, the inverses `2d − 1` and `(4d ± c)/3`, the forbidden
`{7, 21} = {2³ − 1, (2⁶ − 1)/3}` between the realised 27 and 63, and the one-vertex growth ceiling `39 = 5·2³ − 1` are its
faces. The per-letter sharp constants `9, 17, 33 = 2^(a+2) + 1` for a single beaten offset `a` follow from the deletion
count `1 + 2^(a+d)` with `d = 2` and SHARP_EIGHT's witnesses (PROVED by citation of both halves).
*What it proves here.* The owner's second sentence — "the leading coefficient is programmable by an arbitrary binary
word, while the same small reader swaps preserve consecutive exact counts; at least `cX^(2/3)` disjoint repair windows
below `X`, and a bracket of any large target within `O(√X)`" — is the content of
`variable_reader_coverage_20261116/PROOF.md` (written 22:32 by another session): three strong rooted readers share a
continuation state `(Z, I)`; the binary address `b` of `n` bits gives `Z = 82D + 9b` (`D = 2^n`), a complete arithmetic
progression of slopes; the exact count is `H(b) = I + Z(2D − R_n(b)) + g(2b − J(b))`; a collision count over prefixes
(at most `D·⌊18L/164⌋` eligible pairs, each colliding for at most one prefix) leaves `Ω(X^(2/3))` disjoint three-count
windows; the approximate product law `H/D = (82 + 9x)(2 − y) + O(n/D)` with the two half-addresses as coordinates
brackets any large target within `O(√N)`. We replayed `verify.py` ("All seed profiles, literal count triples,
collision lower bound, and target bracketing pass"; receipt `replay_variable_reader_verify.out`) and read the collision
argument; it is sound at the level read. LIMITS.md of the same packet caps the two banks at under `0.173·X` counts, so
this is a density and a bracketing result, not coverage. Verdict: PROVED (cited, replayed), scope as stated by its
authors; Chan–Pak remains open. Mechanism: Θ5 supplies the `+2` steps (reader swaps), Θ1/Θ2 the binary address.

## Θ6. Dual ramification: `Z[i]` against `Z[ω]`

*Theory.* 2 ramifies in `Z[i]` and is inert in `Z[ω]`; 3 ramifies in `Z[ω]` and is inert in `Z[i]`; `Φ₆(x) = Φ₃(−x)`
couples the sign with the cube root of unity. Sums of two squares (`mod 4`), Loeschian numbers (`mod 3`), the square and
hexagonal lattices, the Eisenstein grids and `Z[√−2]` joins of the repository, and the LRC(14) `Φ₆` resonance are
instances.
*What it proves here.* Nothing new; it is the lens under which the repository's two grid families (THM-4825/4826,
THM-5242) are one phenomenon seen from the two quadratic fields of discriminant −3 and −4.

## What was proved in this session, in one place

1. WHEELS.md Theorem A (Fermat levels), Theorem B (reduction to the odd part; necessary density, sufficient nested
   criterion), Corollary C (prime odd part; least periods 24, 80, 112; Zsigmondy's exception as the reason for 24),
   Lemma D (even order ⟺ both signs). Elementary, complete.
2. The identification of the three classical covering sets as full Zsigmondy wheels (`W(24)`, `W(36)∖{109}`, `W(64)`
   with density exactly 1).
3. The von Staudt–Clausen identity behind the Riesel/Niemeier/Leech prime overlap (Θ3).
4. The closed form `2^(a+2) + 1` for the single-offset sharp growth constants (Θ5), by citation of SHARP_GROWTH's count
   and SHARP_EIGHT's witnesses.

What was refuted or typed: the naive guess "coverable iff `12 | N`" (112 = 16·7 covers); pre-registered P-W3 (base-3
periods) by the extra period 108; the Leech-lattice reading of the Riesel wheel (NUMEROLOGY beyond the divisor lattice
of 24); the 41/53 coincidence (COINCIDENCE pending a mechanism).

What remains open: Conjecture W1 (nested is without loss, true to 150), the least periods for `q ≥ 11`, the existence of
a self-dual wheel, the base-`b` analogue of Corollary C, and every `[O]` and `[Q]` item of LIST.md.
