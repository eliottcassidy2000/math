---
id: THM-432
name: Average degree is additive under the Cartesian/Minkowski product; the
      unit-distance product family ties 3N at 27=3^3 and 30, first beats 3N at
      32, so the crossover N* in [25,28] is necessarily NON-PRODUCT
status: VERIFIED (elementary proof + exhaustive factorization over the PROVEN
        small-n optima u(n), n<=21, AMP arXiv:2412.11914)
date: 2026-06-07
session: monad-explorer-2026-06-07-S1
depends_on:
  - THM-431   # u(21)=57; N* in [25,28]; extremal n=21 graph = K3 [] W7
  - THM-421   # combinatorial common-neighbour floor N* >= 17
external:
  - "Erdos, 'On sets of distances of n points', AMM 1946 (the product/Minkowski
     unit-distance construction)."
  - "Alexeev, Mixon, Parshall, arXiv:2412.11914 (2024): exact u(n) for n<=21."
---

# THM-432: avgdeg is additive under [] ; N* is non-product; 27 = K3^[]3

Let `u(N)` be the Erdos unit-distance maximum (OEIS A186705) and `N*` the
smallest `N` with `u(N) > 3N` (average degree `> kappa = 6`). THM-431 placed
`N* in [25,28]`. This theorem isolates WHAT KIND of construction can sit at the
crossover, and gives the exact structure of the clean tie at `27 = 3^3` that the
S710 reflection flagged as "too-clean."

## (A) Average degree is additive under the Cartesian product  [PROVED, elementary]

For finite graphs `G, H`, the Cartesian product `G [] H` (vertex set
`V(G) x V(H)`, with `(u1,u2) ~ (v1,v2)` iff `u1=v1 & u2~v2` or `u2=v2 & u1~v1`)
has
```
   n(G [] H) = n(G) * n(H),
   e(G [] H) = e(G) * n(H) + n(G) * e(H),
```
hence
```
   avgdeg(G [] H) = 2 e/n = 2 e(G)/n(G) + 2 e(H)/n(H) = avgdeg(G) + avgdeg(H).
```
**Average degree is additive under [].** (Equivalently `avgdeg = 2u/n` is a
homomorphism from the multiplicative monoid of graphs-under-[] to (R, +).)

**Unit-distance realizability (Erdos 1946).** If `G, H` are unit-distance graphs,
the generic-angle Minkowski sum `{ g + R_theta h }` realizes `G [] H` as a
unit-distance graph: at a generic relative rotation `theta` no two of the
`n(G)n(H)` points coincide and the only forced unit distances are the product
edges (any extra coincidences only RAISE the count). Therefore
```
   u( n(G) n(H) )  >=  e(G) n(H) + n(G) e(H).            (Erdos product bound)
```
The n=21 extremal graph `K3 [] W7` (THM-431) is exactly this with
`avgdeg = 2 + 24/7 = 38/7 = 5.43 < 6`.

## (B) The product family vs 3N  [PROVED — exhaustive over proven optima]

Define the **product cap** `P(N)` = the maximum edge count of any iterated
Cartesian product of unit-distance graphs on `N` points whose atomic factors are
the PROVEN small optima. Because every nontrivial factorization `N = prod n_i`
(`n_i >= 2`) of any `N <= 42` has all factors `<= N/2 <= 21`, and `u(n)` is
proven EXACT for `n <= 21`, the cap is itself an exact, proven quantity computed
by the recurrence
```
   P(1) = 0,
   P(N) = max_{ d | N, 2 <= d <= 21 }  u(d) * (N/d) + P(N/d) * d.
```
(`unit_distance_product_cap_s1.py`, all integer arithmetic.) The result:

```
   P(N) <= 3N  for EVERY N <= 31.
   P(N) = 3N   (a tie)  iff  N in {27, 30}   (for N <= 31).
   First strict beat  P(N) > 3N  at  N = 32:  P(32) = 98 > 96  via  W16 [] K2.
```

So **no Cartesian-product unit-distance graph beats `3N` for any `N <= 31`**; the
product family first overtakes the penny/kissing rate only at `N = 32`.

## (C) Consequence: N* is NON-PRODUCT  [PROVED]

`N* in [25,28]` (THM-431) lies strictly below the product-beat threshold 32. For
each `N <= 31`, `P(N) <= 3N`, so any unit-distance graph on `N` points with more
than `3N` edges has more edges than ANY product on `N` points — it cannot be a
Cartesian product of two smaller unit-distance graphs. In particular the
realizable `u(28) >= 85 > 84 = 3*28` (Engel/Schade Moser lattice; product cap
`P(28) = 83 < 84`) is an **irreducible** unit-distance graph.

```
   The crossover N* is a NON-PRODUCT (Moser-lattice / irreducible) phenomenon.
```

This **corrects the implicit suggestion** in the S710 handoff that "a product
construction giving 82 at n=27" could lower the ceiling: by (A), every product on
27 points has avgdeg = a sum of atomic average degrees, and for `27 = 3^3` the
maximum such sum is exactly `6` (see (D)) — **no product can give 82 at n=27.**
Only a non-product config can decide `u(27) > 81`.

## (D) The tie at 27 = 3^3 is the Cartesian CUBE of the triangle  [PROVED]

The triangle `K3` is the unique 3-point unit-distance graph; `avgdeg(K3) = 2 =
kappa/3` exactly. By additivity,
```
   avgdeg(K3 [] K3 [] K3) = 2 + 2 + 2 = 6 = kappa,    e = 81 = 3 * 27,
```
so `K3^[]3` is an explicit 27-point, 81-edge unit-distance graph: `u(27) >= 81`.
The other tie factorization `G9 [] K3` (where `G9` is the proven-optimal 9-point
graph, `u(9)=18`, `avgdeg(G9) = 4` exactly) gives the same `4 + 2 = 6`. The two
`avgdeg`-exactly-4 atoms are `G9` and `G10` (`u(9)=18, u(10)=20`); `G9 [] K3` and
`G10 [] K3` are precisely the ties at `27` and `30`. The "3" of "`3N`" is
`kappa/2 = 3`; the triangle carries `kappa/3` of average degree, and **three
triangle-equivalents stacked multiplicatively land exactly on the threshold** —
that is why the tie is at `3^3` and why the deficit sequence touches 0 there.

## (E) Bonus: a better lower bound at N = 32

`W16 [] K2` (the proven-optimal 16-point graph times a unit segment, a "prism")
has `e = 41*2 + 1*16 = 98`, so `u(32) >= 98`. This improves the repo's previous
best `u(32) >= 97` (THM-421(B), a sqrt(7) Eisenstein subset) by one, via a
cleaner construction.

## Relation to canon

Strictly sharpens THM-431 / HYP-2298 (no contradiction; no court case). The
nested ladder of S710 gains an interior rung with a clean meaning:
```
   combinatorial floor 17  <  N* in [25,28]  <  product-beat 32  <  lattice-disk 43.
                              (PROVED non-product)
```
The gap `[27..28] -> 32` between the product TIE and the product BEAT is exactly
the regime the Moser lattice owns. The transferable shape (companion to S710's
"first crossing is boundary-dominated"): **the first crossing of a product-stable
threshold is also IRREDUCIBLE — it is not a product of smaller near-crossings,
because the threshold quantity is a product-homomorphism (here avgdeg) and the
crossing happens strictly inside a product gap.**
