# The fifth root, the golden slope, fences and squares, {2, 3, 11}, and the count modulo 9

Session klein-2026-11-16-wheels2. Typed readings of the remaining items of the directive, each anchored in a repository
packet; verdict vocabulary as in MISTAKE-559.

## 1. Where `Ω(N^(1/5))` comes from, and what it shares with fences and square packing

The live exterior reader (`live_exterior_reader_20261116/EXTERIORS.md`) splits an `n`-bit address into two exteriors
`x, y` of `p = ⌈n/2⌉` and `q = ⌊n/2⌋` bits, reserves `k = ⌊n/10⌋` low phase bits on each (`S = 2^k`), and shows that for
each of the `S²` phase pairs a tangent walk along the product level `A_t B_t` crosses the target with error
`O(S²√R + M)`, `R = P/S`, `M = O(R + T)`. The exponent `1/10` is the balance of the two error terms `S²√R ≍ S^(3/2)D^(1/4)`
and `R ≍ D^(1/2)/S`: equality at `S ≍ D^(1/10)`, giving error `D^(2/5)` and `S² ≍ D^(1/5) ≍ N^(1/5)` free phases per side
(LIVE.md reaches the same exponent by a rational-slope crossing with `L ≍ B^(1/5)`). It is an optimisation exponent, not an
arithmetic one.

The fence and packing bounds the directive names are balances of the same shape. Friedman's axial fences
(`inverse_fences_landing_20261116/FENCE_AND_SQUARE.md`, THM-4505, THM-4509): `A ≤ K(n) = max{k : 2k + ⌈2√k⌉ ≤ n}`, the
area term against the perimeter term `√k`, sharp on quasi-square ominoes. Erdős–Graham's packing of unit squares into a
square of side `s` wastes `O(s^(7/11))` area: the exponent `7/11` is the balance between the waste of a strip of tilted
squares and the waste of the untilted remainder (Roth–Vaughan's lower bound `s^(1/2)` is the other side). All three are
"two-term balance" exponents; none of `1/5`, `7/11`, `1/2` carries number-theoretic content, and the appearance of 11 in
`7/11` is not related to the 11 of `{2, 3, 11}` below. Verdict: ANALOGY of method (balance two errors), no transfer.

## 2. The golden slope

Two golden items exist in the repository and are different objects. `coalescent_27_63_20261116/GOLDEN.md`: consecutive
Fibonacci parameters `(u, v)` give the Euclid triples `(uv ± 1, 2uv, u² + v²)` whose square-sum rotation on `1…b² − 1` is a
two-digit carry, one full low cycle advancing the high digit by one, with `s/m → 1/4` (not the golden ratio); THM-4717:
the balanced rational Collatz cycles converge to circles on which the Collatz step is a rotation, with the golden
conjugacy checked numerically. The live-exterior LIVE.md route crosses fractional phases with a *rational* slope so that
integer samples revisit all phases in a short segment. A golden slope would do the opposite job: by the three-distance
theorem, rotation by `φ` (or by a Fibonacci convergent `F_(k−1)/F_k`) places `L` samples with at most two gap lengths in
ratio `φ`, the best possible equidistribution; this would make the phase-crossing step of LIVE.md uniform in the phase
rather than "any prescribed phase, here 0". Verdict: PROPOSAL (ANALOGY); not implemented, and it does not touch the
exactness gap, which is integer-valued.

## 3. `{2, 3, 11}`

The repository's one exact statement on this set: `procgen_numerology_20260923_odd_square_brackets.md` proves that the
primes `p` for which `2p` lies in the same odd-square bracket `[(2k − 1)², (2k + 1)²)` as `p` are exactly `2, 3, 11`. In
this session's wheel theory 11 appears twice: the period `132 = 4·3·11` has density above 1 and still fails to cover
(the odd wheel on `Z/33` with odd parts `{3: 2, 11: 5, 33: 5}` cannot cover, Theorem A of WHEELS2.md), and `q = 11` is the
first odd prime whose least covering period needs `a = 5` (`352 = 32·11`, since `c(4, 11) = 10`). `2024 = 2³·11·23` and
the Erdős–Graham `7/11` are unrelated appearances of 11. Verdict: the bracket theorem is PROVED (cited); the wheel facts
are PROVED here; any link between them is NUMEROLOGY.

## 4. Prescribing the count modulo 9, and 8 against 9

Two packets prescribe residues modulo 9: `chanpak_cone_collision_placement_20261115/TWO_APEX_HOST.md` (the cofactor
`Q = H/7` of a fixed-arc host hits every residue mod 9 with eight repair bits; the key coordinate "is exactly `2^l` modulo 9
at left length `l`, detecting the parity of `l`", with a period-6 accepting row) and EXTERIORS.md (both strict witnesses
can be required to have `H ≡ r (mod 9)` in every phase cell, because only the term `−82·rev(u_t)` varies modulo 9). Our
census (SHIELD.md §2) confirms that every phase cell of the synchronized reader reaches all nine residues for `n ≤ 22`.
The mechanism is Θ2 of the previous packet: 2 is a primitive root modulo 9, so `2^l (mod 9)` has period 6 and reads the
length parity and residue of a run; the two lowest ternary digits of a count are a function of run lengths modulo 6.
The dual statement modulo 8: `H` is odd (Rédei) and `H ≡ 1 + 2·(#odd cycles) + 4·(#disjoint pairs) (mod 8)` by the
odd-cycle packing identity, so the two lowest binary digits above the forced 1 count single odd cycles and disjoint
pairs; the seed swaps `+2, +4` move exactly these digits. Prescribing the count modulo `72 = 8·9` would therefore fix the
dyadic digit pair by a packing count and the ternary digit pair by a run-length clock — two different kinds of data, as
in every 2–3 interaction of the previous packet.

The "8 and 9" of the square-sum problem (`triangular_backbone_chanpak_20261116/NUMBER.md`: the path on `1…15` with sums
9, 16, 25 has exposed endpoints 8 and 9; GOLDEN.md: the `8 → 9` path is the `(u, v) = (2, 1)` case; THM-4505: the square-sum
chain first exists at 15 with this unique path) and Catalan's `8, 9 = 2³, 3²` (THM-4516's "`2³ − 3² = −1`" clock; the
`K = 3` member of a Pythagorean family) are the same two integers but not the same mechanism as the moduli 8 and 9 of the
residue prescriptions: the endpoints 8 and 9 are degree-one vertices of a sum graph (8 + 1 = 9 and 9 + 7 = 16 are their
only square sums below 16), and Catalan's pair is the unique consecutive perfect powers (Mihăilescu). Verdict: the
mod-9 and mod-8 prescriptions are STRUCTURAL (Θ2, Rédei, odd-cycle packing); their link to the square-sum endpoints and
to Catalan is NUMEROLOGY, as the repository already types `2^3 − 3^2 = −1` outside its Pythagorean-family statement.
