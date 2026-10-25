# Independent review of the support/Cremona packet

2026-10-25. **Verdict: PASS**, in the stated pointed-monoid convention and
complex surface scope. Reviewed `support_reflections/RESULT.md` and
`verify.py`, including the auxiliary family. No substantive mathematical
discrepancy found. The conclusions below rest on the coordinate proofs;
the finite controls supplement them.

## 1. The geometric seven and the degree-six surface

The two equations `x0*y0=x1*y1=x2*y2` have no extra surface component.
Away from the three coordinate points the second projective point is
uniquely determined; over each coordinate point its opposite line is
the exceptional P1. For example near `(e0,e1)`, in charts x0=y1=1,
the equations eliminate `y0=x1=x2*y2`, giving a smooth affine plane.
This verifies the blowup graph and its smoothness, not only its support
count. The monoid charts are correspondingly toric.

The 13 support pairs are precisely the disjoint nonempty pairs and the
full/full pair. Deleting the six singleton/singleton pairs leaves seven
and makes the swap exactly modified complement. The rank and order
change are necessary: the inverse support bijection cannot be continuous.
The finite-field counts, including 16 versus 13 at q=3, correctly expose
the failure of a scheme isomorphism. The original U carrier is not
silently replaced.

Adjunction as a `(1,1),(1,1)` complete intersection gives `-K=O(1,1)`
and degree six. Its seven anticanonical monomials and the 19 products in
degree two give nine quadratic relations. These differ from the hourglass
canonical ring as stated. The six coordinate flags have the adjacent-swap
hexagon and tournament-converse action; boundary components are a distinct
six-element object and are not interchanged without incidence data.

The final fan addition also checks: the ordered primitive vectors
`(1,0),(1,1),(0,1),(-1,0),(-1,-1),(0,-1)` have consecutive determinant
one and are obtained by inserting adjacent sums into the P2 fan. Their
opposition is precisely sign reversal of the normalized hourglass slopes.
The cyclic map `(a,b)->(a-b,a)` and a coordinate exchange generate the
hexagon's dihedral lattice action, with the cube of the cyclic map equal
to `-I`. The seventh label is the zero cone. This is an action on labels
and the lattice, not preservation of branch-line intersection incidence;
the opposite-line versus opposite-boundary hostile is correctly retained.

## 2. Actual hourglass K3 involutions

The matrix B has determinant five. Substituting it into the six specified
line forms gives exactly the displayed transformed list. Multiplying each
of the three paired products under `(X,Y,Z)->(YZ,XZ,-XY)` yields
respectively `-X^2,-Y^2,+Z^2` times the same pair product. Their product
is a square, establishing the lift identity. The squared base map is
projectively the identity.

Blowing up the three selected double crossings subtracts twice each
exceptional divisor from the sextic, producing `-2K` on the degree-six
del Pezzo surface. Local normalization `w=xv`, `v^2=t` resolves the
three corresponding A1 nodes. The twelve remaining crossings are all
in the base torus and form six exchanged pairs. Thus their exceptional
curves cannot contribute fixed points after resolution.

The four fixed base points have `U=±1,V=±i` and are unbranched. The
eta-preserving lift therefore fixes both sheets above each, giving eight;
the eta-negating lift fixes neither. Opposite boundary components of the
hexagon are disjoint, so no boundary fixed point was omitted. Inversion
negates each logarithmic one-form, preserving their wedge; consequently
the two-form signs are exactly as claimed. The free anti-symplectic
involution on the smooth K3 gives an Enriques quotient. The explicit
local/fixed-locus argument handles this nodal specialization rather than
appealing to an unverified smooth-branch substitution.

The weighted-plane signs are also consistent: the product of the three
new coordinates is `-(XYZ)^2`, so preserving eta sends w to `-XYZ*w`.

## 3. Obstruction to lifting back to individual roots

All radicand ratios defining the six-root multiquadratic field have zero
valuation at X=0. Their squareclass span therefore has even valuation
there. The pullback ratio
`X*(Y-2Z)/(Y*(X+2Z))` has valuation one. A base-field element becomes a
square in an elementary multiquadratic extension exactly when its
squareclass belongs to the generating span, by its quadratic subfields.
Thus no field automorphism above this Cremona map can exist. Adding the
seventh root does not repair the witness because E=X-Y-Z has valuation
zero on the same axis. This is a specific nonlift, not a claim about all
K3 automorphisms.

## 4. Auxiliary lambda family, including the boundary

On the torus the square-coordinate map to a,b is etale. The sum of the
logarithmic critical equations factors as `(a+b)(ab-1)/(ab)`. On ab=1,
the remaining equation is `(a-1)(a+1)(a^2+a+1)/a^2`; on a+b=0 it is
`a-a^-1`. Evaluating g gives precisely 6,-2,-3, hence the three excluded
lambda values -6,2,3. This split exhausts the torus critical locus.

On each open coordinate line, the remaining equation has the simple
roots of `1+t^2`. At an exceptional divisor, the quadratic tangent cone
has the same two distinct roots. The strict transform avoids each
corner of the boundary hexagon. Consequently there are no additional
boundary singularities, and the stated iff smoothness condition holds.
The lambda=2 factorization checks the singular boundary positively.

The smooth double cover of the del Pezzo has trivial canonical bundle
and q=0, so is K3. Pullback of -K has square 12 and eight sections via
the cover algebra, not seven; the manuscript keeps this distinction.
The resulting genus-seven polarization is an auxiliary construction,
not the canonical polarization of the seven-root general-type surface.

## 5. Independent replay

The producer's final 126 explicit checks passed. The separate
`check_character_support_independent.py` uses integer polynomial
dictionaries, not SymPy or producer imports, to verify the sextic identity;
exact Fraction cross-products enumerate all 15 crossings and their six
remaining paired orbits. It also counts the four unbranched fixed fibres
over F17, where both sheets split, obtaining 8 and 0 for the two lift types.
This finite-field check is a control, not the complex fixed-locus proof.
Normal and optimized independent runs must agree.

External-source priority is not audited here. The standard classical
identifications are used only in their conventional scopes, and the
specific geometry of the displayed hourglass follows from the proof above.
