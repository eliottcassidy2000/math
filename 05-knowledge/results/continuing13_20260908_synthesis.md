# Continuing synthesis: fixed units, iterated fibres and local repair cores

**PREVIOUS CHECKPOINT, 2026-09-08. Status: PROVED scoped results,
FINITE-EXACT controls and independent audits.** A complete univariate
family lowers the fixed-unit hidden-arm construction from order six to
order three. A polynomial mutation produces a genuine global quintic
with moving affine-fibre moduli and a complete target-collision law.
Two more LRC clocks close. Original-cell separators localize exact grid
repair, including connected instances with unbounded global concurrency.
A fixed critical orbit then keeps three unit supports and every source
parameter unchanged while the hidden component multiplicities grow.
General LRC(14), planar JC and extremal no-three-in-line remain **OPEN**.
No external priority claim is made.

The [current synthesis](continuing14_20260908_synthesis.md) classifies the complete zero-budget clock slice, reconstructs geometric generic fibres, excludes source-linear logarithmic witnesses and separates repair blocks by line-capacity states.

The [previous synthesis](continuing12_20260908_synthesis.md) supplies the
complete quadratic rational-mate classification, the order-six hidden-arm
family, native middle-depth obstruction and exact hyperedge branching.
This session starts from bcfbfa23d. No incoming commit was present at
startup or the natural research checkpoint. The
[manifest](continuing13_20260908_manifest.json) pins all proofs and replays;
the original response ring and the finite computation universes remain
explicit in their suppliers.

## 1. Portfolio, inheritance and changed coordinates

The anchor was the next two declared LRC clocks. The niche was original
source degree three after the completed quadratic table. The wildcard
was the incidence topology omitted by global exceptional-cell branching.
The six live concepts were complete physical graphs, original atoms,
polynomial primitives, target-value collisions, complete unit modules,
and the state passed through a separator.

Closest mechanisms were the full 126-profile forest consumer, coefficient
spans in the component quotient, the two actual surface charts, and exact
capacity-two matching. Hostiles were a disconnected zero graph with no
isolated vertex, an affine submersion that fails on the added divisor,
mistaking a rational mate for a rational parametrization, and inseparable
cell choices whose benefit appears only jointly. Underused sidecars were
an arbitrary univariate polynomial in the second global coordinate, its
actual primitive, target multiplicities, and articulation CELL states.

The first useful cubic did not have unit order one. Instead its coordinate
choice generalized to every degree and exposed a much simpler fixed-unit
family. Applying a polynomial repair to that family then changed the
fibration itself. This goes beyond taking more derivatives of its old
unit, which the previous session proved insufficient.

## 2. The same order-three pointed module has arbitrarily many ambient arms

On fixed W2 put

    u=x-h, z=u-2h+u^3t, A=u f(z).

For this ENTIRE ansatz, globality is exactly f(0)=0. Global submersion is
then exactly squarefreeness of f together with f(-2h)!=0. For q=deg f>=2,
let Q'=f^2 and Q(0)=0. The
[complete family](continuing13_20260908_univariate_fixed_unit.md) and
[independent combined audit](continuing13_20260908_quintic_mutation_audit.md)
prove

    J(A,Q)=A^3, G=Q/A^3, J(A,G)=1,
    C(x,t)=C(A)(z),
    [C(x,t):C(A,G)]=2q+1.

Every nonzero source fibre is exactly C[z,f(z)^-1]. The zero fibre has
q+1 reduced disjoint components: u=0 and one for each root of f. Hence
the complete original-source torsion has q arms. The degree-three case
already exceeds the two-arm ceiling of the completed quadratic class.

The complete principal parts in the SAME target parameter g=A are

    Eu: Q(-2h)/g^3 + f(-2h)/g^2,
    E_rho: Q(rho)/g^3.

There is no simple term. The root rho=0 supplies a regular component and
must remain before quotienting by the common diagonal. If all Q(rho)
vanished, f^3 would divide Q, contradicting 3q>2q+1. Thus the two actual
coefficient directions are always independent. After a constant basis
choice the unit is always

    theta=C3/g^3+C2/g^2.

Its exact pointed Weyl presentation is

    D / ( Dg^3 + D[1+partial*g+(partial^2/2)*g^2] ).

The arbitrary-operator converse uses the determinant-minus-one columns
g*theta and g^2*theta as a basis over C[partial]. This is the SAME pointed
module for every admissible f,h and every q>=2. The unit generates exactly
two arms, while q ambient arms and rational pair degree 2q+1 grow without
bound. Every rational constant-Jacobian mate gives the same pair field.
At q=1 the coefficient directions coincide; this is an explicit boundary.

The theorem exhausts this univariate ansatz, not arbitrary global cubics.
Unit order one on W2 remains open outside the previously closed classes.

## 3. A polynomial mutation reaches an actual quintic with moving punctures

The [mutation theorem](continuing13_20260908_quintic_mutation.md) uses the
actual polynomial Q, whose information the old unit quotient discarded.
For lambda!=0 and lambda+f(-2h)!=0, define

    T=lambda*A+Q(z), H=1/(2A^2).

Then J(T,H)=1 and T is a global submersion on W2 of original t-degree
2q+1. Away from A=0 this follows from J(A,Q)=A^3. On Eu the actual
source derivative is f(-2h)[lambda+f(-2h)]; on the entire added divisor
the tangential derivative is -lambda*f'(0). The second forbidden lambda
is a genuine critical-line hostile despite the rational Jacobian identity.

Let

    S={Q(-2h)} union {Q(rho):f(rho)=0}.

For c outside S the WHOLE source fibre is

    C[z,f(z)^-1,(c-Q(z))^-1].

It is P1 minus 3q+2 points. The roots of f and infinity are fixed; the
2q+1 roots of Q(z)=c move. These punctured curves are non-isotrivial even
without labels: an isomorphism is determined by the images of three fixed
punctures, giving only finitely many possibilities for any fixed reference
curve. Distinct c give different remaining root sets. The projective
completions all have genus ZERO; completed genus is not what varies here.

Every special fibre has one irreducible regular component in addition to
the A=0 components assigned to its target value. Their complete primitive
principal parts are pure double poles:

    E_rho: lambda^2/[2(T-Q(rho))^2],
    Eu: [lambda+f(-2h)]^2/[2(T-Q(-2h))^2],
    regular component: 0.

There is no simple pole, including on Eu. Consequently the full torsion
has q+1 arms, but the unit generates exactly |S| arms, one per distinct
target value. Its exact scalar annihilator is

    product_(c in S)(T-c)^2.

The degree is 2|S|; only the LOCAL pole order is two. The pointed unit
module remembers S and loses its component multiplicities. Unit cyclicity
is equivalent to all q+1 values being distinct. The intrinsic rational
pair degree is 2(2q+1), and no change of rational mate lowers it.

For an explicit quintic take h=lambda=1 and f=z(z-1):

    u=x-1, z=u-2+u^3t,
    T=u z(z-1)+z^5/5-z^4/2+z^3/3,
    H=1/[2u^2z^2(z-1)^2].

It has three full and generated arms, exact annihilator

    [T(T-1/30)(T+256/15)]^2,

pair degree ten and non-isotrivial eight-punctured ordinary source fibres.
Choosing -2h on 6z^2-15z+10=0 gives a genuine admissible collision: three
ambient arms remain, but the unit generates two. This first function is
an actual object at the live quintic frontier; its nonzero unit explicitly
rules out a polynomial mate. It does not constitute a Keller entry.

## 4. Complete zero cuts close the two declared clocks

The [full probe](continuing13_20260908_lrc_zero_clock_probe.md),
[native cut theorem](continuing13_20260908_lrc_zero_cuts.md), and
[independent unpruned audit](continuing13_20260908_lrc_independent_audit.md)
close both predeclared clocks 9240 and 11088. Their full-profile word
counts are 22,469 and 18,780, with minimum margins 108 and 156. Both have
zero overlap budget, so hypothetical unsafety requires every actual
strict edge to have zero overlap.

Every possible zero-edge graph is disconnected. Almost every word has
an isolated vertex, but exactly FOUR do not. Their full nontrivial
component cuts are the required hostile controls against reducing graph
disconnection to isolation. No native-depth continuation is needed for
these two clocks; the previous 7560 depth obstruction remains valid.

The independent audit starts from 1,545,830 unpruned seven-multisets,
reconstructs all 41,249 proper-profile words and every one of 101,836
compatible native ratios, and reproduces all 4,084 zero ratios and
positional cuts. Every proper profile and the actual connected-complement
scope remain in force.

The necessary array now has **7,618 clocks, maximum 11,935**, SHA256
`89a3e7545cc77467c5f85fbe0bdaab1f71071cff65184bbc1212eb86c9f07177`.
Only 9240 and 11088 were removed here. General and disconnected-complement
LRC entry remain open.

## 5. Separators localize exact original-cell repair

The [cell-separator theorem](continuing13_20260908_no3_cell_separators.md)
and [independent audit](continuing13_20260908_no3_separators_audit.md)
replace a single global branch set by exact local continuation states.
Start with the bipartite graph of ORIGINAL overfull lines and original
active cells. Merge biconnected blocks through articulation LINE vertices,
so each bag owns every incidence of its line constraints. Articulation
CELL vertices connect the bags in a forest and retain two states.

If h_B counts cells meeting at least three owned lines inside bag B, the
exact computation uses at most

    2 sum_B 2^(h_B)

weighted integral b-matching calls. The messages count each original cell
once. Their weights can be negative, already on the smallest connected
example; replacing them by cardinalities destroys the continuation cost.

For three directions, exceptional articulation cells never contribute to
h_B. An explicit connected integer-geometric family has k such cells,
2k+1 overfull triples, no local branching, and exact repair k. Its occupied
rows and columns have one point; it is not a saturated two-regular board.
In four or more directions an articulation cell may still have three
incidences inside one bag and must remain in the local branch set.

Two geometric hostiles identify the boundary. An inseparable two-cell
theta has deletion-state costs 3,4,4,2: neither unilateral improvement
reveals the optimum. Another board has three genuine articulation cells
on a cyclic block; removing all articulation cells first and reconnecting
components incorrectly recreates a cycle. The actual block-cut construction
keeps that bag whole. The audit's cycle-excess corollary makes cactus
incidence graphs sufficient for zero local branching, without an iff claim.

Every one of the 2,137 complete two-regular boards of sizes two through
five has h_B=0 for the stated three directions. A dimension-six board
with rows (01,01,23,45,35,24) has one local exceptional cell and repair four.
The first fourteen lexicographic dimension-six boards suffice to discover
it; no full dimension-six census is claimed. No random-board frequency,
mean improvement or extremal bound follows from these finite controls.

## 6. A finite critical orbit fixes the unit while multiplicities grow

The [iteration theorem](continuing13_20260908_fixed_support_iteration.md)
and its [independent referee](continuing13_20260908_fixed_support_iteration_audit.md)
turn target collisions into an unbounded family with FIXED source parameters
and the SAME complete pointed unit module. Choose

    zeta^2+zeta+1=0, alpha^2=zeta-1, P(z)=z^3+alpha.

The critical orbit is 0 -> alpha -> zeta*alpha -> zeta*alpha. Zero never
returns to zero. Let eta be either root of

    z^2+zeta*alpha*z-zeta^2=0,

so eta is a fixed point distinct from 0, alpha and zeta*alpha. Fix
h=-eta/2, lambda=1 and sigma^2=3 once. For every integer k>=2 put

    Q_k=P^k-zeta*alpha,
    f_k=sigma^k product_(j=0)^(k-1) P^j,
    q_k=(3^k-1)/2.

Here powers of P mean composition. Disjoint simple preimage levels of
zero prove that f_k is squarefree, while the chain rule gives Q'_k=f_k^2.
At the fixed source point, f_k(eta)=(sigma*eta)^k is nonzero and never -1:
eta is an algebraic integer, whereas a hypothetical equality would force
eta^(2k)=1/3^k, a rational noninteger algebraic integer. Thus the global
mutation hypotheses hold for ALL k with the same h and lambda.

Apply Section 3 with T_k=A_k+Q_k, A_k=u f_k(z). Its complete support is
the same set {0,a,b}, where a=(1-zeta)*alpha and b=eta-zeta*alpha. The
ambient arm counts at these three supports are exactly

    ((3^(k-1)-1)/2, 3^(k-1), 1).

Each special fibre also has one regular component. The unit sees just
one pure double-pole arm at each support, with exact scalar annihilator
[tau(tau-a)(tau-b)]^2. Rescaling its three nonzero coefficient vectors
identifies the FULL pointed unit-generated Weyl modules, with no target
coordinate change. Their ambient embeddings and multiplicities differ.

| k | Original t-degree | Ambient arms | Intrinsic rational-pair degree | Ordinary-fibre punctures |
|---|---:|---:|---:|---:|
| 2 | 9 | 5 | 18 | 14 |
| 3 | 27 | 14 | 54 | 41 |
| 4 | 81 | 41 | 162 | 122 |
| all k>=2 | 3^k | (3^k+1)/2 | 2*3^k | (3^(k+1)+1)/2 |

For each fixed k, the punctured ordinary fibres are non-isotrivial as
the target varies; their projective completions remain rational. Every
rational mate has the same embedded pair field. The finite controls at
k=2,3,4 check 115 always-active exact gates; the all-k claim has a separate
analytic proof and independent audit.

This supplies a precise dynamics-to-response connection: compositional
iteration preserves the critical-value set and multiplies its preimage
levels; square roots of derivatives feed the global mutation theorem.
The missing sidecar is component multiplicity over each retained value.
The cheap hostile P=z^3 has the same small-support appearance, but its
second derivative square-root product is z^4, which is not squarefree.
No-return of the critical point, not finite support alone, pays the map.

## 7. What the new connections require next

The univariate family gives a stronger hostile to recovering ambient
geometry from the unit. The polynomial mutation restores target values,
yet still loses their multiplicities. In grid repair, global concurrency
overcounts the hard part, whereas local core choices genuinely interact.
In LRC, an isolated vertex is a useful sufficient witness but the whole
cut is the complete state needed by these clocks. Each connection has an
actual preserved predicate and an explicit lost coordinate.

Next tests: retain critical-value multiplicities and branch configuration
when comparing rational pairs; search outside the closed univariate and
quadratic classes for unit order one; use entire native zero cuts before
lifting a clock's missing absolute depths; and measure exceptional cells
inside inseparable incidence bags before making probabilistic claims.

The completed packet contains 37 frozen artifacts, nine exact engines
and eight certificates. Normal and optimized runs reproduce identical
raw output across 613,410 always-active gates per mode. All 958 inherited
artifact pins and all 32 artifacts of the first continuing13 checkpoint
remain unchanged. The five additional iteration artifacts match their
independently accepted pins. These finite validations supplement the
scoped analytic proofs; they do not close the general conjectures.
