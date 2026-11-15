# How balanced can the two addresses be?

This applies to the physical six-site three-choice ripple bank in
[INVERSE.md](INVERSE.md), not to unrestricted tournaments.

Fix any constant \(\kappa>1\) and insist that
\(\kappa^{-1}\le x/y\le\kappa\). Up to a cap X, the bank can hit
at most

    (log κ)/72 + o(1)

of the odd targets, with natural logarithm. The same upper bound on the
fraction holds within any fixed fractional band [aX,bX], 0<a<b.

Indeed, each of the three counts equals 432xy+O(x+y+1), uniformly
because 0≤p(z)≤z. The comparable-address constraint makes x+y=O_κ(sqrt X)
for any presented count ≤bX, since the positive term 432xy dominates.
Every contributing pair therefore lies in the product annulus

    aX−O_κ(sqrt X) ≤ 432xy ≤ bX+O_κ(sqrt X).

Use product coordinate u=xy and ratio coordinate v=log(x/y).
The positive-quadrant area element is (1/2)du dv, and the v-interval
has length 2 log κ. Thus the wedge has area
(b−a)X log κ/432+O_κ(sqrt X). Lattice count agrees up to O_κ(sqrt X):
slice at integer x; each permitted y interval incurs at most two
rounding errors, and the endpoints are piecewise monotone with
total variation O_κ(sqrt X). The same reasoning applies to a=0.
At most three targets can come from each pair. Dividing by the
number of odd targets proves the bound.

So even if all possible collisions vanished, a decoder confined to any
fixed aspect ratio κ<exp(72) would fail to hit density one for this
bank. This conclusion involves capacity, not a congruence obstruction;
choosing the center by residue class cannot rescue the shortfall.
Practical searches should either admit substantially unbalanced
addresses or change the palette of center packets rather than search
only square-like baselines. The exact inverse already admits all
imbalances; its worst-case bound does not use this result.
