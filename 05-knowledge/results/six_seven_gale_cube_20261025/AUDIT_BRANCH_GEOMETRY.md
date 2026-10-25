# Independent audit: branch cubes and five elliptic factors

Reviewer: arithmetic_shapes. Date: 2026-10-25.
**Verdict: PASS in the stated characteristic-zero and no-fourfold scope.**
No mathematical blocker or required correction found.

Reviewed frozen drafts in `/tmp/six_seven_drafts/geometry`:

- `BRANCH_CUBE_AND_DOUBLE_COVER.md`, SHA256
  `e8cc2b8b3d4c29b0d5713ae0fdb1f2bea9fd664fc83b5725e5ea002e0410be5e`.
- `GENUS_FIVE_ELLIPTIC_FACTORS.md`, SHA256
  `00946c30fb3422db14dec33778fa6aa066664aba0a614ecfe9deabcbd6e6187f`.
- `check_geometry.py`, independently rerun with **10,576 passing gates**.

## General B1/B2 mechanism

The repeated pairs on the new line are disjoint: an old line has only one
intersection with the new line, and an old triple there would violate the
no-fourfold hypothesis. Thus there are exactly d=n-r distinct restricted
zeros. After selecting the r proportional-root signs, valuations at these
zeros give independent square classes and component degree 2^(d-1).
This proves connectedness and the Riemann–Hurwitz genus without assuming
the numerical component count. The total degree also checks:
2^r times 2^(d-1)=2^(n-1).

Locally at a repeated pair the equation alpha*u^2+beta*v^2 factors into two
transverse lines over the geometric ground field. The two labels differ
in precisely that pair coordinate. Each component has 2^(d-2) points over
the corresponding place. Hence the cube multigraph and the absence of
self-nodes follow from local branches, not just node totals. Its edge
count is r*2^(r-1)*2^(d-2)=r*2^(n-3), agreeing with the independent fiber
count. The divisor B is 2H and misses the old exceptional curves.

Blowing up exactly those branch crossings gives B'=2(sigma*H-sum E).
The local model w^2=uv resolves by the stated normalized double cover;
the exceptional curve has self-intersection -2. Old nodes split etale
locally, and their exceptional P1 curves lift twice over the algebraically
closed field. The canonical blowup terms cancel exactly. The resulting
canonical divisor is nef for n>=6, excluding exceptional (-1)-curves.
Thus the finite map is to the blown-up base, not incorrectly to Y_n.

## Actual hourglass and field sidecars

For outer-root deletion the new line has one repeated pair, c and h,
with h^2=2c^2. This yields two genus-five components over Q(sqrt(2)).
Their eight common points c=h=0 instead lie over Q(i), as the displayed
independent sign coordinates show. The argument excluding Q-points of
the branch forces every coordinate to vanish, while real points away
from c=0 do exist. These two component/intersection fields are not conflated.

For center deletion the branch conics satisfy the three signed linear
relations and h^2=a^2+d^2. Differing in one sign forces one of a,d,h to
zero and gives two transverse intersection points; differing in two signs
forces all three zero and gives none. There are no triple conic crossings.
All 24 crossings are Q(i)-rational, and conjugation changes all signs.

The three Picard relations are supported by actual hyperplane cuts:
b-ia=0, e-id=0, k-ih=0, together with their opposite signs. They are not
deduced merely from a singular intersection matrix. The matrix
-2I+2A(Q3) therefore has the stated geometric rank-five span and true
three-dimensional relation space. No full Picard-rank claim is made.

The dual conic is nonsingular, its six marked points and A2 vectors are
correct, and the displayed rational parametrization is exact. Under the
Mobius matrix [[2,-1],[1,1]], the sextic multiplier is -27. Multiplying
the matrix by a rational scalar changes this multiplier by a sixth power,
so cannot remove its nonsquare class. A projective branch-set symmetry
therefore does not silently become a Q-automorphism of the hyperelliptic
cover. The Kummer step is explicitly OPEN, as it should be.

## Elliptic appendix and cohomological direction

The four cubic character products and one quartic product have four odd
branch places each; all other nontrivial characters have two. Each genus-one
quotient has a rational point. The five displayed changes of variables
preserve actual equations and twists, including x=-t for the second E0
and (x,Y)=(2/t,2y/t^2) for E2. Equality of j alone is not used.

The five pullbacks of regular one-forms have distinct deck characters and
are nonzero by separability, so they span H^0(C,omega_C). This is exactly
the cotangent-space isomorphism of the product norm map Jac(C)->product E.
That map is an isogeny, and duality gives the displayed pullback isogeny
in the other direction. All maps are over Q. The proof does not claim an
integral polarization decomposition or a specified kernel.

The computation covers all 252 six-old/one-new transitions among the nine
magic labels, plus symbolic identities and explicit hostile boundaries.
This audit certifies the supplied elementary proofs and their numerical
applications; it does not turn the finite census into a proof for untested
line arrangements. Those are handled by B1/B2 themselves.

## Integration amendment

The Kummer OPEN sentence above records the reviewed draft, not the current
frontier. Concurrent Session60 subsequently closed the central identification
over C and determined Picard ranks 20/19/58. The integrated proof routes to
those exact results and retains both rational twist factors with their maps.
Status labels and reproduction paths were updated after this audit; B1/B2
and the elliptic proof are unchanged. See [AUDIT.md](AUDIT.md) for final bytes.
