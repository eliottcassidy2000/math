# Move the centre, touch the connector, control the digit

**Session269, 11 November 2026. Research complete.**

The session develops three native mathematical advances, rather than auditing
old claims or treating a shared formula as a solution. The strongest new planar
Jacobian step is an explicit certificate operation: moving the Gaussian centre
turns two pure-direction first hits into an entire mixed vanishing tail. The
strongest tournament step is exact target-directed dyadic control inside a
fixed-score, rigid, strong family, with a matching order exponent. The geometry
step supplies both a general interacting Euclidean construction and a positive
connector-sensitive equality family.

**JC(2), Hadwiger–Nelson and Chan–Pak remain open.** JC(2) here is the constant
nonzero Jacobian problem over C in dimension two. Chan–Pak's target is the
positive odd Hamilton-path spectrum with the known omissions 7 and21, not a
claim that those two values should be attained. No literature-priority claim
or promotion of inherited pre-audit results is made.

## 1. Planar Jacobian: a missing implication becomes a certificate compiler

For a Keller pair F=(f,g), let D,E be its commuting dual derivations and put

    U=f-a, V=g-b,
    M_(a,b)=(D-V) C[x,y] + (E-U) C[x,y].

This is a vector-space quotient, not a quotient algebra. Multiplying a zero
class by an arbitrary polynomial would be invalid.

The new elementary theorem says: suppose that for every monomial h and every
Gaussian centre in some nonempty open set, at least ONE positive D-derivative
and ONE positive E-derivative of h lie in M_(a,b). The orders and witness
degrees may vary without any assumed bound. Then there are rational
certificates producing

    [D^m E^n h]=0 whenever m>=R_h or n>=S_h

away from one algebraic denominator locus for that h. At very general centres,
every polynomial therefore has only finitely many surviving mixed derivatives.

The mechanism has two genuinely different parts. Finite witness membership is
a constructible rank condition. A countable-cover argument extracts one
rational certificate per direction. Then the flat centre connections

    nabla_a = partial_a+V,       nabla_b = partial_b+U

commute with both Stein operators. If the parameter-independent target k has
certificate k=(D-V)p+(E-U)q, two explicit local rules give certificates for
Dk and Ek. Differentiation raises powers of existing denominators but creates
no new pole locations. Thus the same permitted centres support the entire
infinite compiled tail.

There is a sharper single-centre version. Let k be the countable field generated
by the coefficients of f,g. At ONE centre algebraically independent over k,
pure first hits for every monomial already imply mixed finite chaos at that
same centre. Nonzero minors remain nonzero at such a centre, so the rational
certificate descends directly by linear algebra. This is not detection at an
arbitrary numerical centre. If a finite Weyl generating set is supplied, it
is enough to establish the first hits for those generators.

**Conditional consequence.** The repository's THM-5138 finite-chaos criterion
is explicitly PRE-AUDIT and uses an inherited comparison theorem. Assuming it,
the new propagation theorem shows that universal K-GMC, cone vanishing, and
one-sided Stein nilpotence are equivalent to JC(2), with their full quantifiers.
One no longer needs a separate arbitrary-single-centre detection theorem to
close that UNIVERSAL implication. The first-hit existence assertion itself is
not proved for all Keller pairs.

[Full proof, finite-generator sharpening, generic-centre theorem and controls](jacobian/proof.md).

## 2. Hadwiger–Nelson: interacting geometry and the precise support boundary

Session254 realized independent four-edge paths at two ports; their six-colour
channels were 105^r and104^r. This session admits arbitrary interactions in
the middle of those paths.

If an induced unit-distance core H lies in an open disk of radius2, there is
an exact induced Euclidean realization with two generic ports A,B and ears

    A - x_w - w - y_w - B

for every core vertex w, retaining every edge of H and adding no unintended
unit edges or coincident points. The graph has3|H|+2 vertices. The proof
eliminates a finite union of lower-dimensional incidence conditions; the
all-pairs experiment uses exact radicals rather than distance tolerances.

For fixed port colours a,b, its response is exactly

    Z_H(a,b) = sum_(proper c of H) product_w
               (q-2+1[c(w)=a])(q-2+1[c(w)=b]).

At every q>=3 all factors are positive. Consequently **every port assignment
extends if and only if the core is q-colourable**. Even arbitrarily complicated
middle-core interactions cannot delete the unequal-port state while preserving
the equal state. The next interaction must reach the connector layer or change
the attachment itself. This is a proved location-of-interaction result, not a
claim that every geometric approach fails.

There is further retained structure: Z is a Gram matrix, with
`f_a(c)=product_w(q-2+1[c(w)=a])`. Its equal-channel count D is at least
its unequal-channel count O, and

    D-O = (1/2) sum_c (f_a(c)-f_b(c))^2.

Equality occurs exactly when every proper q-colouring of the core has all
q colour classes of equal size (vacuously also for an uncolourable core).
A genuine unit-triangle core at q=3 has D=O=24. The interacting square core
at q=6 has D=57807360 and O=56451584. The boundary statistic is the full
colour-class-size pattern, not merely the number of proper colourings.

**Positive forcing calibration.** The connector edge of a unit diamond really
does force equality at three colours. Rotating and gluing two diamonds gives
seven-vertex induced equality gadgets at every port separation in (0,2sqrt(3))
except finitely many values, including arbitrarily small positive distances.
An exact instance at squared separation1/18 has response(4,0) at three colours
and(48,32) at four. The latter is an explicit guard against promoting a
three-colour gadget to the five-colour obstruction needed for a six-colour
lower bound. No improved chromatic bound for the plane is claimed.

[Geometric construction, Gram theorem, equality family and exact tests](hadwiger/PROOF.md).

## 3. Chan–Pak: control every low binary digit without changing scores

For the inherited Moon interval chains, a composition of N into parts p_i>=2
has Hamilton count product_i(1+2^(p_i-1)). The new ingredient is a bank of
same-size switches whose count differences have exact dyadic valuations:

    j=1: (5) versus(2,3),             H17 versus15;
    j=2: (7) versus(3,4),             H65 versus45;
    j>=3: (2j) versus(j,j),           valuation j.

Choose one switch at each level1,...,k-1. The resulting2^(k-1) tournaments
are strong and rigid, all on exactly k^2-k+7 vertices, all with the same score
sequence and the minimum n-2 cyclic triangles among strong n-tournaments.
Their H values represent **every odd residue modulo2^k exactly once**.

A prescribed residue is constructed by k-1 exact switch decisions. The first
changed switch determines the first changed binary digit; later switches
cannot disturb it. This is an adjacency-producing operation, not signed
interpolation or a count-only existence table.

The quadratic order exponent is optimal within ALL Moon interval chains.
At fixed total N their distinct counts are bounded by the partition number
p(N)<=exp(pi sqrt(2N/3)). Covering2^(k-1) residues therefore forces

    N >= (3(log2)^2/(2pi^2)) (k-1)^2.

Infinite switch words even give an exact2-adic isometry onto the odd2-adic
integers, carrying fair bits to uniform Haar measure. This is a limit of
counts, not one finite tournament with a noninteger path count.

**The integral boundary is visible.** Every chain count is coprime to7.
The family approximates63 to arbitrary dyadic precision but never realizes63,
even though the inherited repository has actual H63 tournaments. At precision8
the decoder returns a63-vertex chain with H12644053665875775, congruent to63
modulo256. Thus complete local control, strongness and fixed scores still do
not supply exact integer landing. This family is a useful arithmetic testing
ground, not a proof of the spectrum conjecture.

[Proof, lower bound, isometry and physical target constructor](chanpak/PROOF.md).

## 4. A concrete bridge: keep state-specific budgets before scalar elimination

Session267 gave a reduced-rational criterion for the source-row cone
`deg F_m<=floor((m+5)/2)`. A new forward certificate makes that criterion local.
For a polynomial state system y=b+A(t,x)y, with positive t-degree in A, assign
weight w(t^r x^j)=2j-r. State budgets h_v satisfying

    w(b_v)<=h_v,
    w(A_vs)+h_s<=h_v,
    w(c_v)+h_v<=5

prove the desired bound for every row of c^T y. Different budgets are essential:
a positive-weight transition can be paid by a negative return transition.

For nonnegative coefficients this is an exact criterion on the trimmed path
system: no accepting positive-weight cycle, and no overweight simple accepted
path. For signed coefficients it is sufficient but not necessary, because
unstable paths can cancel in the scalar output. This is a finite certificate
for a given carrier; deriving the actual Jacobian source carrier remains open.

The bridge is deliberately typed. In the three native settings, the missing
sidecar is respectively a rational centre certificate, connector incidence,
and integer height beyond residue. There is no identification of a tournament
with a Keller map or an arbitrary colouring graph with a unit-distance graph.
The geometric Gram positivity is actual sum-of-squares positivity; complex
Gaussian/Stein vanishing does not inherit it merely from the word “Gaussian.”

[State-budget theorem, exact positivity boundary and controls](bridges/FILTERED_TRANSFER.md).

## 5. Best next attacks

1. **Jacobian:** construct one pure first-hit certificate in each direction
   from forward Keller data for a supplied Weyl generator. The new compiler
   then supplies the mixed tail; deriving the initial certificates is the
   substantive remaining task.
2. **Hadwiger–Nelson:** constrain a connector using five distinct forced
   neighbour colours in an exact geometric network. Middle-core changes
   alone have now been ruled out for the two-ear construction. A unit circle's
   unit-chord graph is bipartite, so a local clique cannot furnish that palette.
3. **Chan–Pak:** add a positive non-Moon seam that permits factors of7 while
   retaining a useful digit inverse and controlling integer size. A dyadic
   residue bank can diagnose such a move; it cannot replace exact landing.
4. **Cross-lane:** derive a finite source recurrence with a checkable
   state-budget potential, or show which genuinely necessary memory makes
   a finite recurrence impossible. A finite-prefix rational fit is not enough.

## 6. Evidence and publication

All new infinite assertions have written proofs. Computations check only their
explicit operations and hostile controls, not inherited conjectural claims:
32,764 dyadic residue targets,58 literal tournament controls,80 high-precision
arithmetic targets,112 exact geometric point-pair checks, and exact symbolic
Stein compiler/commutator identities. Normal and optimized Python runs execute
the same explicit checks. [Reproduce](REPRODUCE.md) · [Results](FINAL_CHECKS.json)
· [Pinned sources](SOURCE_PINS.json) · [Publication scope](PUBLICATION.md).
