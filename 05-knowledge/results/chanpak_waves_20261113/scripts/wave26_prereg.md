# Wave 26 pre-registration (mac-mini-2026-11-15-slides)

Registered before running any wave-26 computation. Timestamp below is the write time.

Context: owner prompt of 2026-11-15 (stable plus/minus fan pairs; the growing window as an
inner repair with a different left reservoir; slides that preserve what the surrounding
tournament reads; number-theoretic invariants). Prior art read first:
chanpak_leaf_fans_20261115 (fan law, cross-port law, movable-port census, opposite fans,
parent compiler) and chanpak_window_descent_20261115 (matched slides 13/29, half-steps
145/39, CRT amplifier, exposure theorem, primitive lanes).

The leaf-fan notes treat only single-vertex ports (newborn beats exactly one old vertex, or
loses to exactly one). Here a port is an arbitrary out-set X: the newborn beats X and loses
to the rest. Odd-cycle form: H(T + x_X) = H + 2 W(X), W(X) = sum over u in X, v not in X of
Pi(u,v), Pi(u,v) = sum over even-order paths u -> v of H(complement).

P77 (arming is common at small order). At order 6, at least half of the isomorphism
classes with H not in {3, 11} admit both unit extensions 2H+1 and 2H-1 by one new vertex
with some out-set.

P78 (adaptive depth). With the out-set re-chosen freely at every node, some seed of
order <= 6 carries a full binary tree of unit extensions to depth 5 (orders up to 11).

P79 (decay). No seed of order <= 6 carries an adaptive full binary tree to depth 9.

P80 (fixed rules fail). With fixed port sets X+, X- in the seed (newborns transitive in
birth order, older beats younger), no seed of order <= 7 carries a full binary tree to
depth 4.
registered 2026-11-14 20:51:52 MST

## Addendum 1 (written while the order-7 slide search and the (8,7) coverage run were running)

Exploratory, NOT pre-registered (already run): the unit-extension census and adaptive trees
(P77-P80 outcomes recorded in the README), the base-4 two-vertex trees from orders 3-5, the
passive-palette census (orders 1-8), the mod-3 law t = (-1)^m s + [m odd] a (found in data,
then proved), the slide census for orders 3-6, and the marked two-block coverage for block
orders (6,6) and (7,7): every odd N in [23, 67967) resp. [23, 1542343) except none.

P81. The marked two-block family with block orders (8, <=7) covers every odd N in
[23, X) with 3*10^7 <= X <= 2*10^8.

P82. The order-7 matched-slide search (tile width <= 28) finds no unit pump (|Delta| = 2).

P83. Every nonzero slide step found at orders <= 7 is even, and no core of order 6 has a
nonzero matched slide at any width (extrapolating the width-30 search).
addendum 1 registered 2026-11-14 21:25:31 MST
