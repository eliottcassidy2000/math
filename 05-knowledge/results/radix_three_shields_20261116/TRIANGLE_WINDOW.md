# Four exact states through a live triangle shield; 17, 68, and the inverse guard

## Positive factor-three test

Write \(q_j(U)\) for the number of covers of \(U\) by \(j\) ordered
directed path slots, including empty slots, and \(C_2=q_2-2H\) for
the two-slot count with both slots nonempty. Let \(W(U)\) add two
vertices \(x,y\), with \(U\to x\to y\to U\). This is a strongly
connected output for any nonempty \(U\), not an ordinal sum with a
disjoint triangle. Any successful input switch leaves the shield and
all its contacts literally unchanged.

The archived order-seven atlas gives a five-state palette. Exact subset
DP verifies these aligned outgoing row masks independently:

| H(U) | Rows (vertices 0 through 6) | H(W(U)) |
|--:|:--|--:|
| 91 | (62,60,88,112,96,68,3) | 1359 |
| 93 | (62,28,88,112,32,70,19) | 1365 |
| 95 | (62,28,72,80,36,78,19) | 1371 |
| 97 | (62,76,72,80,38,78,17) | 1377 |
| 99 | (62,76,104,16,38,74,25) | 1383 |

Successive states differ by exactly two internal arc flips, with all
shield arcs retained. Their profiles are

\[
q_j(U_i)=j^7+20j^5+60j^3+(10+2i)j,\quad 0\le i\le4.
\]

Thus \(C_2=1086\) **constant**. The order-nine shielded outputs form five
exact states at step six (and so include a width-four window). These are
not asserted to lie in a single independent two-arc square. The checker
also gives a separate order-nine input palette, at H=367,369,371,373 and
C₂=10122, with literal shielded counts 11223,11229,11235,11241.

A separate complete check of all 33,867 labeled tournaments of orders
one through six finds no constant-\(C_2\) width-four palette, even
allowing mixed orders in the bank. The only constant-\(C_2\) consecutive
three-term run is 11,13,15 at \(C_2=78\). Reproduce with `minpal.cpp`:
the only line is `78: 11,13,15 and next=0`. Hence order seven is the
minimum allowed input order for an entire width-four window in this
**constant-\(C_2\), live triangle wrapper** architecture.

In fact the complete readout is
\[
q_j(W(U))=j\{q_{j+1}(U)+j q_j(U)-q_{j-1}(U)\},\qquad j\ge1.
\]
If \(x,y\) are in different slots, there are \(j(j-1)\) choices and
the inside factor is \(q_j(U)\). If in the same slot, there are \(j\)
choices. Write \(S\) for the vertices of \(U\) using that slot:
for nonempty \(S\), a path either uses \(x\to y\) consecutively,
leaving two optional \(S\)-runs, or goes \(y\to S\to x\), a single
nonempty \(S\)-run. Thus the contribution is \(q_{j+1}(U)+q_j(U)\)
before subtracting \(q_{j-1}(U)\) for the double count when \(S\) is empty.
In particular \(H'=q_2+H=C_2+3H\) and \(C_2'=2q_3+2C_2\).

## Why the first factor is not an odometer return

Define \(Mf(z)=z(f(z+1)+zf(z)-f(z-1))\). The five profile states
have adjacent differences exactly \(2z\); after \(k\) identical
triangle wraps they still form five equally spaced counts, but the
steps are \(2(M^kz)(1)=2,6,30,246,\ldots\). In particular the first
wrap leaves \(C_2'=19566,19578,19590,19602,19614\), so the
interface that certified multiplication by three was not returned.
The next round's exact step is thirty rather than eighteen.

Indeed **no** two fixed nonempty input tournaments with different
path counts can yield exact multiplication by three at every depth
under this common wrapper. Let \(f\) be the difference of their
slot polynomials. These are polynomials by grouping nonempty slots;
they vanish at zero for a nonempty tournament. If the difference at
height \(k\) were \(3^k f(1)\) for all \(k\), put \(g=(M-3)f\).
Then \((M^kg)(1)=0\) for all \(k\), and \(g(0)=0\). The highest
evaluation in \((M^kg)(1)\) is \(k!g(k+1)\), from always shifting up;
all other terms use \(g(0),\ldots,g(k)\) (coefficients at zero stop
because of the prefactor \(z\)). Induction gives \(g\equiv0\).
But \(M\) raises the degree of a nonzero polynomial by exactly two,
a contradiction. This excludes this *common wrapper*, not a sorted
or state-dependent factor-three repair.

This explains why the independent-shield requirement matters. In
particular two separate single-arc toggles at the *same* apex never
have zero mixed difference: the one-apex odd-cycle cut formula has
mixed coefficient \(-2(K_{ab}+K_{ba})\) where \(K_{ab}\) is the
positive-weighted even-order path kernel from one contact to the
other in the deleted tournament. At least the oriented two-vertex
path contributes a positive term. Multi-arc swaps or independent
shields are not excluded.

## The exact 17/68 bridge, and the different things named 68

- Reverse the long endpoint arc of a transitive six-vertex
  tournament. Its directed odd cycles are the reversed arc plus odd
  subsets of four intervening vertices. No two coexist in a pack.
  Its slot polynomial is \(j^6+8j^4+8j^2\), hence \(H=17\),
  \(q_2=224\) and \(C_2=190\). The **same live triangle shield**
  transforms this literal 17 into \(17+224=241\), strong at order
  eight. In the archived period-24 Riesel wheel, 17 and 241 were
  actual labels, of orders eight and twenty-four, linked by the
  arithmetic identity \(17\cdot241=2^{12}+1\). Indeed
  \(2^4=-1\pmod{17}\) and \(2^{12}=-1\), \(2^8=15\pmod{241}\).
  The count identity is exact; it does not make clock order a
  tournament invariant.
- \(68=4\cdot17\) is an **even midpoint**, not a Rédei count.
  The [67/69 seam](../chanpak_pair_67_69_20261114/README.md)
  is a physical move straddling 68, not a realization of it.
  For the earlier [signed Collatz cycle](../collatz_three_mirrors_20260929.md),
  \(-68\to-34\to-17\) is the two-edge return along the same
  odd class. Its dyadic valuation is a coordinate forgotten by the
  odd-inducing projection. The independent
  [union-graph minor study](../procgen_kuratowski_20260925_tait_kempe_triples.md)
  also finds a \(K_5\) threshold at 68; there is no map transferring
  that graph threshold into a path-count shield.
- \((8,15,17)\) is primitive; quadrupling gives \((32,60,68)\),
  on the Pythagorean cone but outside the primitive tree. Another
  occurrence is \((51,68,85)=17\cdot(3,4,5)\), also imprimitive. The
  [archived Berggren chart](../ternary_berggren_20260925.md)
  has positive inverse branches selected by disjoint ranges on a
  marked ratio, with a decreasing parent. In that chart 17's
  primitive triple corresponds to roots \((5,3)\), whose parent is
  \((3,1)\) via the branch \((2s-t,s)\).
  [Price's other ternary tree](https://arxiv.org/abs/0809.4324)
  likewise distinguishes an equation from a selected primitive child.
  The shield analogue is to keep the marked continuation data,
  not just the scalar gap, in any claimed return. Exact gap-six
  alone is no more a proof of repeatable shielding than an
  imprimitive Pythagorean point is a tree node.

## Multiplicative reduction (odd primes only)

Using the known exclusion of 7 and 21 and closure under ordinal
products, Chan–Pak is equivalent to realizing all **odd** primes
other than 7, all \(7p\) for odd primes \(p\ne3,7\), and the seeds
49, 63, 343. Indeed for \(N=7^v M\), \(7\nmid M\), use odd-prime
products when \(v=0\); when \(v\ge2\), decompose \(v\) into twos
and threes and use 49 and 343. For \(v=1\), use \(7p\) for a
prime factor \(p\ne3\) of \(M\), or 63 if \(M\) is a power of
three with exponent at least two. The other cases are exactly
7 and 21. The seeds are arithmetically indispensable to **this
product deduction**: each proper factorization of them contains
7 or 21. Their three-way role is not an identification with
the geometric root triple \((3,4,5)\).
