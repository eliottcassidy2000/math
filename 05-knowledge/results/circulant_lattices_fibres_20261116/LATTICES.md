# The packing-polynomial lattices: delayed choices at every level, and what the diagonal cannot see

Owner: *"Also consider related to the packing-polynomial lattices."* THM-5269 found that two tournaments of order `n` with
equal `(h, P_2)` can differ in `P_3` only by elements of `120ℤ` (orders 7, 9) or `360ℤ` (order 8). Here the general law,
its proof, and the census of what the diagonal `P(q) = G(2q, q)` forgets.

## 1. The delayed-choice lattice at every level (PROVED)

At order `n` the diagonal polynomial is `P(q) = q^n + Σ_{m ∈ M} a_m q^m` with `M = {m ≡ n (mod 2), 1 ≤ m ≤ n − 2}` and
`a_m = Σ_{c+s=m} 2^c N_{c,s}` (THM-5269). Two tournaments with equal `P_1, …, P_k` have a coefficient difference `Δa` in
the integer kernel of the `k` rows `(c^m)_{m ∈ M}`, `c = 1, …, k`. Writing `m = m_0 + 2j`, row `c` is `c^{m_0}` times the
geometric row `(u_c^j)` with `u_c = c²`, so the kernel is that of the Vandermonde rows in `u_1, …, u_k`.

**Theorem.** The integer kernel has the ℤ-basis of *band vectors*: the coefficients of `Π_{c=1}^{k} (X − c²)` placed on
`k + 1` consecutive elements of `M`, for each starting position. Consequently the lattice of `ΔP_{k+1}` over pairs of
order `n` with equal `P_1, …, P_k` is

    L_n^{(k+1)} = (k+1)^{m_0} · Π_{c=1}^{k} ((k+1)² − c²) · ℤ,     m_0 = 1 (n odd) or 2 (n even),

whenever `|M| ≥ k + 1`, i.e. `n − 2 − m_0 ≥ 2k`.

*Proof.* A band vector `v` at start `j_0` has `Σ_i v_i u^{j_0+i} = u^{j_0} Π_c (u − c²) = 0` for `u = c²`, so it is in the
kernel; the band vectors are triangular with leading entry 1 (the coefficient of `X^k`), so any integer kernel vector is
reduced by them, from its highest position downward, to a vector supported on the lowest `k` positions, which must be
zero (a square Vandermonde system). Hence they are a ℤ-basis. The image of the band at start `j_0` under
`Δa ↦ Σ_m Δa_m (k+1)^m` is `(k+1)^{m_0} (k+1)^{2j_0} Π_c ((k+1)² − c²)`, and the gcd over `j_0 ≥ 0` is the `j_0 = 0`
value. ∎

| level `k+1` | band vector (increasing `m`) | `Π ((k+1)² − c²)` | lattice, `n` odd | lattice, `n` even | first `n` |
|---|---|---|---|---|---|
| 3 | `(4, −5, 1)` | 40 | `120ℤ` | `360ℤ` | 7 / 8 |
| 4 | `(−36, 49, −14, 1)` | 1260 | `5040ℤ` | `20160ℤ` | 9 / 10 |
| 5 | `(576, −820, 273, −30, 1)` | 72576 | `362880ℤ` | `1814400ℤ` | 11 / 12 |

([receipts/lattices2.out](receipts/lattices2.out): band basis against the formula for `n ≤ 14`; the first-batch
`lattices.out` found `120/360` for level 3 by an independent enumeration.) The constants are
`Π_{c=1}^{k} ((k+1)−c)((k+1)+c) = k!·(2k+1)!/(k+1)!`: `40 = 2!·5!/3!`, `1260 = 3!·7!/4!`, `72576 = 4!·9!/5!`. The
concurrent pair of order 8 differs by `720 = 2·360`; at order 9 the 338 split groups of THM-5269 differ by multiples of
120 (240, 2160, 2400, 2640).

## 2. What the diagonal forgets ([receipts/fibres2.out](receipts/fibres2.out))

The concurrent `packing_hadamard_flatness_20261117` proves that a difference of bivariate packing polynomials divisible
by `(x − 2y)` is invisible to every uniform Hamiltonian read, and exhibits two order-6 tournaments with identical
`P(q) = q^6 + 8q^4 + 8q^2` and `G_A − G_B = x(x − 2y)`. Counting that phenomenon:

| order | classes | distinct diagonals `P(q)` | distinct bivariate `G` | classes whose diagonal is shared with a different `G` |
|---|---|---|---|---|
| 5 | 12 | 8 | 8 | 0 |
| 6 | 56 | 24 | 29 | 16 |
| 7 | 456 | 156 | 203 | 196 |
| 8 | 6880 | 1728 | 3007 | 4503 |

At order 8 two thirds of all classes carry a diagonal that some class of different cycle structure also carries; the
example `(95, 2552, 26127)` is shared by four packing polynomials differing in `N_{1,3}` and `N_{2,2}` by multiples of
`(−4, +2)`, i.e. by `2 x y^3 (x − 2y)`-type terms. H-3 and H-5 were refuted: the diagonal is far coarser than
predicted (1728 of 6880), and the hidden share far larger (4503). So the "packing-polynomial lattice" has two parts:
along the diagonal, the levels of §1 (what a next evaluation can still reveal); transverse to it, multiples of
`(x − 2y)` (what no uniform read reveals, recovered by vertex deletions as the concurrent packet shows).

## 3. The circulant lattice

For a circulant of prime order `p`, THM-5270 (B) gives `G ≡ y^p + ((p−1)/2) x (mod p)` coefficientwise, so the whole
coefficient vector `(a_m)` lies in the coset `(…, 0, …, 0, −1) + pℤ^{M}` and every `P_q` in `pℤ`; two circulants of order
`p` with equal `P_1, P_2` differ in `P_3` by an element of `L_p^{(3)} ∩ pℤ = 120p ℤ` (`p ≥ 7`): at order 7 the pair
`175/189` has `ΔP_3 = 42` only because their `P_2` differ (1988 vs 2016).
