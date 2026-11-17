# Plus or minus one: a catalogue of the ±1 structures in the programme, typed

Owner: *"identify as many as you can that could be viewed as ±1, and plus or minus anything is just plus or minus 1 with a
multiplier. See how these sets of digits relate to number theory and tournament structures and other concepts beyond."*
Each item: the structure, where it is proved or recorded, its type.

1. **The tournament itself is a ±1 matrix.** `S = A − Aᵀ` is skew with entries ±1 off the diagonal; skew-Hadamard
   tournaments (`S Sᵀ = (n+1)I − J`, orders `≡ 3 (mod 4)`) are the nuclear-norm extremals (THM-4814). The contact variables
   of THM-5266 (A) are spins `s_x = 2·[x ∈ X] − 1`, and the extension count is the Ising-type energy
   `h + Σ (B_x − e_x)(1 + s_x)/2 − Σ M(y, x)(1 + s_y)(1 + s_x)/4`; a seam is one spin flip with energy change `±2`
   ([SEAMS.md](SEAMS.md)). PROVED (the form) / ANALOGY (the physics).
2. **Fermat and Mersenne.** Near-joins over transitive bases give `2^m + 1` (3, 5, 9, 17, 33, …); the holes are
   `7 = 2^3 − 1` and `21 = 3(2^3 − 1)`, the Mersenne shadows of 9 and 27 (THM-5266 (C)); `63 = 2^6 − 1` is Zsigmondy's
   exception and the period-7 Mandelbrot count, divisible by both holes exactly for periods `p ≡ 1 (mod 6)` (THM-5265 (E)).
   PROVED.
3. **Busch's floors are products of Fermat-type numbers.** Moon's interval chains have `H(IC(p_1, …, p_k)) =
   Π (2^{p_i − 1} + 1)` (checked on all 33 compositions of `N ≤ 8`, [receipts/heavy_ties_ic.out](receipts/heavy_ties_ic.out)),
   so the strong floors `3, 5, 9, 15, 25, 45, 75, 125` are products of `3 = 2 + 1`, `5 = 4 + 1`, `9 = 8 + 1`: "+1 with a
   multiplier" twice over. KNOWN (Moon 1972, THM-4744) / arithmetic.
4. **The tower criteria `2^k B_t = h ± 1`.** `+1` is the digit-1 condition that forces `3h + 2` from a tower (THM-5260 (C));
   `−1` is the charge-0 condition (`a_t = B_t − e_t − h + 1`, wave 45 / THM-5260). The same quantity, one unit either side.
   PROVED.
5. **The rooted append count `2h − 1 + a_t`** (wave 45; `D(h) = 2h − 1` in THM-5272): "−1 with multiplier 2". PROVED.
6. **Digit and response.** `B_t = h + σ_t`, `σ_t ≥ 0` (THM-5263 (I1)); digit 0 is out-degree ≤ 1 (THM-5264); the stability
   `e = g` of the concurrent `center_126` packet is the same `e_t` against the deletion count. PROVED.
7. **Counts are `≡ ±1 (mod 4)`.** Rédei: `H` odd; THM-4922 (4): at odd orders the self-converse classes split evenly
   between `H ≡ 1` and `H ≡ 3 (mod 4)`; THM-791: values `3 (mod 4)` with even least order are black. PROVED.
8. **Collatz `3n ± 1`.** The window digits `3q + 2 = 3(q + 1) − 1` and `3q + 4 = 3(q + 1) + 1` (THM-5256): the +2 window is
   a seam with centre `3(q + 1)` and multiplier 3 — `30021/30027 = 3·(10008 ± 1)`. PROVED (the identity) / ANALOGY (Collatz).
9. **Twin primes are seams of the integers.** `11/13` (centre 12, both near-join exceptions), `41/43` (the restart row),
   `10007/10009`, `10^9 + 7/10^9 + 9`, `30011/30013`; one-sided: `7 = 8 − 1`, `47 = 48 − 1`, `67 = 68 − 1`, `21 = 20 + 1`.
   The contact census places every even centre in `[264, 10000]`; 12 and 42 are placed, 20 and 22 cannot be. NUMEROLOGY
   (primes) / FINITE-EXACT (centres).
10. **Paley and Legendre.** Tournaments from quadratic residues need `−1` to be a non-residue, `p ≡ −1 (mod 4)`; the
    Legendre tournament on primes `≡ 3 (mod 4)` (HYP-9375). KNOWN.
11. **Riesel and Sierpiński numbers `k·2^n ± 1`** and the Zsigmondy wheels covering them (THM-5254, THM-5255). PROVED there.
12. **Signed binary digits.** Every integer is `Σ ±2^k` (non-adjacent form); a near-join pencil `h + 2^k B` is one signed
    digit with multiplier `B`; the `+2` of a seam is the digit `2^1`. KNOWN / the owner's phrasing made literal.
13. **The alternating identity** `Σ_S (−1)^{|S|} H(D[S]) H(D[V∖S]) = 0` (wave 56) and the odd-cycle packing weights
    `2^{#cycles}` with inclusion–exclusion signs (THM-002, THM-4114). PROVED.
14. **Deletion monotonicity `H(T − v) ≤ H(T)`**: every Hamiltonian path of `T − v` admits an insertion of `v` (at the front
    if `v → p_1`, at the end if `p_{n−1} → v`, otherwise before the first `p_i` with `v → p_i`); this is the "+ at least 0"
    behind the heaviest-parent tree ([EXT10.md](EXT10.md)). PROVED (classical).
15. **Contact and flip digits.** The contact digits `(B_x − e_x) − Σ (M + Mᵀ)` take every value in `[0, 3638]`, the flip
    digits `M(b, a) − M(a, b)` every value in `[−153, 153]` at orders ≤ 8: "plus or minus anything" is a one-contact move
    with that multiplier. FINITE-EXACT.
16. **`n(h) = m(h)` or `m(h) + 1`** (THM-4889 / HYP-9481: least order against the cyclic-triangle level). FINITE-EXACT.
17. **Locker parity**: `H(D_11) ≡ 3 (mod 4)`, the τ-parity involution (THM-865/866). PROVED.
18. **The seam itself**: a prescribed even centre `c`, counts `c ± 1`, realised by one port contact (this packet).

What these share: in every tournament item the `±1` is a single arc — the arc between the new vertex and one old vertex
(contacts), one arc of `T` (flips), the one back arc of an interval-chain block (`2^{p−1} + 1`), the one reversed arc of
the near-join (`h + 2^k B`) — and the multiplier is a path count through that arc (`M`, `B`, `2^k`). The number-theoretic
items (`2^m ± 1`, twin primes, `3n ± 1`, signed digits) share only the arithmetic shape; no mechanism joins them to the
tournament side beyond the Fermat/Mersenne reading of the holes (item 2). Typed accordingly.
