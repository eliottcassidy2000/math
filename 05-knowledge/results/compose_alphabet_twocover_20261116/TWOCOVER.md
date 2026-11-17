# The two-cover count is the packing polynomial at (4, 2): selecting the smaller block, frozen fibres, delayed choices

Owner (the "hadamard" line): *"hold the two-cover count fixed while changing the smaller Hamiltonian count, and hold
both counts fixed while storing choices that only show up after another wrap."* The concurrent
`two_cover_feedback_20261117` proves this constructively over a bank of dyadic tiles (`u, R = 5, 2^5; 19, 2^8; 223, 2^12;
599, 2^13`), reads `15u + 62R` through its two-run wrapper, and stores delayed choices with two order-8 blocks of equal
`(P(1), P(2)) = (345, 5640)` and different `P(3)` (45225, 44505). Here the same objects are explained by one identity and
censused over all classes.

## 1. The identity (PROVED)

For a tournament `T` let `P_q(T)` be the number of ordered `q`-tuples of disjoint, possibly empty vertex sets covering
`V(T)` weighted by the product of their Hamiltonian counts (`P_1 = h`, `P_2` the two-cover count, `P_3` the three-cover
count), and let `G_T(x, y) = Σ_π x^{c(π)} y^{s(π)}` over the spanning packings `π` of `T` by directed odd cycles of length
≥ 3 (`c` of them) and singletons (`s` of them). THM-002 is `h = G_T(2, 1)`. Then

    P_q(T) = G_T(2q, q)      for every q ≥ 1.

*Proof.* Expand each slot's count by THM-002 and interchange the sums: a `q`-slot cover together with one packing per
slot is a spanning packing of `T` together with an assignment of each of its components to a slot; a cycle carries the
weight 2 and `q` choices, a singleton the weight 1 and `q` choices. ∎ Replay: all 530 classes of orders 3–7, packings
enumerated directly, `q = 1, 2, 3`: 0 failures ([receipts/packing.out](receipts/packing.out)). Examples:
`G = y^3 + x` for the 3-cycle (`3, 12, 33`), `y^4 + 2xy` for the strong order-4 tournament (`5, 32, 117`),
`y^5 + 5xy^2 + 2x` for `RT_5` (`15, 120, 525`); the Paley tournament `P_7` has `(189, 2016, 12159)`.

**Consequences.** (i) `P_q` is a polynomial in `q` of degree `n` whose coefficients `a_m = Σ_{c+s=m} 2^c N_{c,s}` vanish
unless `m ≡ n (mod 2)` (an odd cycle of length `ℓ` lowers `m` by the even number `ℓ − 1`); `a_n = 1`. (ii) `P_2` is even
and `P_q ≡ 0 (mod q)`. (iii) The signed two-cover sum `Σ_W (−1)^{|W|} h(W) h(W^c)` vanishes for every nonempty `T`: each
component contributes `(−1)^{|component|} + 1 = 0` (wave 56's alternating identity; 0 failures on 530 classes). (iv) `P_q`
is multiplicative over ordinal sums (no cycle crosses the cut), which is the concurrent packet's bank arithmetic.
(v) **Frozen fibres and delayed choices are rank statements.** Two tournaments of order `n` with equal `h` and `P_2`
have coefficient differences `Δa` orthogonal to `(1, 1, …)` and `(2^m)`; the lattice of possible `ΔP_3 = Σ Δa_m 3^m` is
`120ℤ` at order 7, `360ℤ` at order 8 (`Δa = t(1, −5, 4)` on `(a_6, a_4, a_2)`) and `120ℤ` at order 9
(`ΔP_3 = 1680Δa_7 + 120Δa_5`). So a block that keeps `(h, P_2)` and stores a choice for the next wrap must change its
packing numbers along that lattice — the concurrent pair differs by `720 = 2·360`.

## 2. Census ([receipts/twocover.out](receipts/twocover.out) orders ≤ 8, [receipts/collisions9.out](receipts/collisions9.out) order 9)

- Orders ≤ 8: 7,412 classes, 978 distinct two-cover counts (13 %); 602 of them (62 %) are shared by classes of
  different count — the frozen-cover fibres, including the dyadic ones `2^5: {1, 5}`, `2^8: {1, 5, 19}`, `2^9: {5, 19}`,
  `2^{10}: {19, 25}`, `2^{12}: {217, 223}`, `2^{13}: {599}` (the concurrent bank plus 25 and 217). Of the 1,816 `(h, P_2)`
  collision groups exactly **one** is split by `P_3`: `(345, 5640)`, the concurrent pair (C-3's "50 %" refuted).
- Order 9: 191,536 classes, 50,173 distinct `(h, P_2)`, 48,870 collision groups, **338** split by `P_3` (C-8's "≤ 100"
  refuted), with differences 240, 2160, 2400, 2640 — all in `120ℤ` as (v) requires.
- Selecting the smaller block with the right two-cover count: the least-order block for each `P_2 ≤ 150` is
  `2: M_1, 4: M_2, 8: M_3, 12: C_3, 16: M_4, 24: (order 4, h = 3), 32: strong-4 or M_5, 48: (5, 3), 64: (5, 5),
  84: (5, 9), 96: (6, 3), 100: (5, 11), 104: (5, 13), 108, 120: (5, 15), 128: (6, 5), 144: (6, 9)`; every `P_2 ≥ 4` is
  divisible by 4.

## 3. What this adds to the concurrent packet

Its frozen-fibre selection theorem and the non-dyadic bridge `(11, 100)` stand; the identity above says what a block's
two-cover count *is* (its packing data at `x = 4, y = 2`), why ordinal juxtaposition multiplies it, why the alternating
sum vanishes, and that delayed choices live in a lattice of dimension one at order 8 and two at order 9. Nothing here
is absolute coverage. Chan–Pak OPEN.
