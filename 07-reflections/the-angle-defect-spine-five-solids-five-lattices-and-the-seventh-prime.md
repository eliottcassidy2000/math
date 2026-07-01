# The angle-defect spine: five solids, five lattices, and why the seventh prime is the frontier

*mac-mini-2026-06-30-S65. Prompted by "the 5 plane tilings correspond to the 5 Platonic solids — where does
their count come from, and what does it connect to in the proofs?"*

Two famous fives. Five Platonic solids. Five two-dimensional lattices (the crystallographic orders
`{1,2,3,4,6}`). It is tempting to read the shared "5" as coincidence or mysticism. It is neither. Both are the
same inequality read at two signs.

Take the `(2,3,n)` triangle and its angle defect `1/2 + 1/3 + 1/n - 1`. Positive defect (`n = 3,4,5`) is a
**sphere**, and its regular figures are the five Platonic solids. Zero defect (`n = 6`) is the **plane**, and its
regular figure is the hexagonal tiling. Negative defect (`n ≥ 7`) is the **hyperbolic plane**, and `n = 7` is
the Klein quartic, `PSL(2,7)`, `168`, the Fano plane. The count of five solids and the count of five planar
lattices are the *same* Gauss–Bonnet dichotomy — curvature positive versus curvature zero — sorted by whether
`2cos(2π/n)` is an integer (`φ(n) ≤ 2`). Nothing is coincidental; the fives are two faces of one trichotomy.

Now the thing this project has been circling for a year drops into place. LRC(2p) is governed by its apex prime
`p`, and `p` walks the very same `(2,3,p)` spine. `p = 3`: order three, crystallographic, Euclidean — LRC6 is
trivial, `X_0(6)` has genus 0. `p = 5`: the icosahedron — still spherical, still a Platonic solid, `X_0(10)`
still genus 0, still tractable. `p = 7`: the first prime that is **neither** — it tiles no sphere and no plane,
it is purely hyperbolic, `X_0(14)` jumps to genus 1, and the cusp form `f_14` is born. The modular genus
`0, 0, 1, 2, 2` for `p = 3,5,7,11,13` is not a separate arithmetic fact bolted onto the geometry; it *is* the
geometry — the genus jump at 7 is the spherical-to-hyperbolic transition.

So LRC14 is hard for a reason you can point at on a triangle. The five Platonic solids use up the spherical
orders `{3,4,5}`. Order 6 is the lone Euclidean one — and it is exactly the covering-min: the hexagonal `A2`
lattice, whose Dedekind anomaly `s(n, Φ_6) ≠ 0` *is* the safety margin (last session). Order 7 is the first
escape into the hyperbolic, and the apex floor gap `4cos²(3π/7) = 2 + λ_min(C_7) > 0` is nothing but the
crystallographic *impossibility* of a seven-fold lattice — the heptagon never closes, so its cycle has a
spectral gap, and that gap is the floor. The reason there is no seven-fold crystal and the reason the lonely
runner has room at the seventh prime are one reason.

The whole picture is a seam. On one side, order six, hexagonal, Euclidean, the covering-min, the Dedekind
anomaly, `B2` its anomaly-free square shadow. On the other side, order seven, Klein quartic, hyperbolic, the
apex, the Fano plane, the genus-1 cusp form. LRC14 lives exactly on the seam between the last plane tiling and
the first hyperbolic surface — which is precisely why it is the first case the sphere and the plane can no
longer help with. The five solids and the five lattices were the two shores; `7` is where the water starts.

*See [[HYP-3771]], [[HYP-3768]] (the order-6 Dedekind margin), [[HYP-3586]] (X₀(14) genus / apex cusp),
[[HYP-3606]] (the 4cos²(3π/7) apex gap), and [[everything-is-the-triangle]].*
