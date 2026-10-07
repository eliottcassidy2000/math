# Local defects, endpoint charges, and finite grounded receipt networks

2026-10-07. **PROVED:** the elementary finite-network equivalences, exact
residual identity, receipt export and refinement statements below.
**FINITE-EXACT:** the declared small controls. The six supplied paper
headlines are **CITED / ACCEPTED PREMISES**, not re-audited. No universal
Collatz completion, independent all-source positive floor, or new source
coverage is claimed by this package.

[Program](../../04-computation/experiments/paper_local_global_grounding_20261007.py)
and [exact output](paper_local_global_grounding_20261007.out).

## 1. Inheritance and the paper mechanisms

The closest proved mechanism is the undirected ROOT component of an
authenticated actual-word graph in
[fair frontier extension, section 2](fair_frontier_extension_20261004.md).
The [source-sensitive Poisson dual](collatz_poisson_source_dual_20261005.md)
already distinguishes a genuine global subsolution from a finite witness
that contains its own ROOT path. The
[Poisson patch theorem](collatz_poisson_patch_20261005.md) already retains
transported forcing during Schur elimination. We do not claim these
mechanisms as new. Here the finite graph is symmetrized only because a
checked common future gives a bidirectional **certificate** implication;
we do not symmetrize the physical Collatz transition kernel.

Canonical hostile: a component of compatible receipts disconnected from
ROOT. Corrected near miss: small coordinatewise defects do not imply a
small global dual-energy error. Least-used sidecar: the original endpoint
charge and the resource convention attached to a compressed edge. The
board is **source / authenticated edge / terminal charge / local residual /
global energy / retained path / actual deadline**.

The supplied files were read at the following specific interfaces:

| Trusted source | Mechanism retained here | Transfer boundary |
|---|---|---|
| [Geometric Erdős similarity](C:/Users/Eliott/Downloads/geometric-erdos-similarity.pdf), sections1.2,5, PDFpp2--3,10--12 (repair p12) | Finite local routing controls normalized scales; an open neighborhood repairs the closed set of missed centers because translated geometric points approach their center | Collatz terminal equality has no such open-neighborhood property |
| [First lattice correction to Bloch's law](C:/Users/Eliott/Downloads/first-lattice-correction-bloch-law.pdf), section3, Lemma3.2, PDFpp7--8 | A trial potential needs a scalar residual and a global dual-energy bound; both residual orientations occur for a nonselfadjoint transition | Our graph has a deliberately symmetric Laplacian; its simpler identity cannot replace both estimates in the paper's nonsymmetric model |
| [Planar XY logarithmic correction](C:/Users/Eliott/Downloads/paper-9.pdf), introduction and section5, PDFpp52,65--68 | A marginal coordinate has reciprocal drift and size about1/j, with summable higher-order errors | Decay to zero is not a finite ROOT hit |
| [Complete-graph crossings](C:/Users/Eliott/Downloads/paper-10.pdf), section2, PDFpp3--5 | The cycle space is the kernel of the signed endpoint boundary; local endpoint terms accompany interior intersections | A cycle has zero boundary, not the required source-to-ROOT charge |
| [Complete-bipartite crossings](C:/Users/Eliott/Downloads/paper-11.pdf), section3, PDFpp5--8 (Figure1 p7) | Endpoint order terms survive local gluing; Figure1 was inspected visually | We retain incidence charges, without transferring a crossing-number bound |
| [Half-plane double dimers](C:/Users/Eliott/Downloads/paper-12.pdf), introduction and sections2,4, PDFpp3--11,18--24 | Puncture/topological data alone omit excursions and repeated traversal; traversal control is additional data | A quotient or a near-return is not a finite authenticated route with clocks |

The useful positive transfer is a complete finite grounding checker that
retains these boundary data. The capacity and flow arguments below are
standard elementary network mathematics, proved here for the exact typed
object; no literature-priority claim is made.

## 2. The authenticated graph and its signed boundary

Let U be the accelerated Collatz map on positive odd integers. A stored
edge between x and y carries two exact finite valuation words a,b with

    U_a(x)=U_b(y),

and their actual positive odd endpoints. The verifier checks the words,
integer labels and matching endpoint. The words in this package stop before
any padding after ROOT; an empty word is allowed. It follows that

    x reaches ROOT iff y reaches ROOT.

Reversing an edge reverses a certificate implication, not time in an
unverified integer orbit. The graph has a distinguished vertex ROOT1.
Arbitrary integer labels on a synthetic test graph are not authenticated
Collatz edges.

Choose an orientation for each edge e=(u,v). Its incidence column is

    B_e = e_v-e_u.

For a fixed source s different from ROOT r, the following statements are
equivalent on any finite graph:

1. There is an undirected path from s to r.
2. There is a rational edge vector f with Bf=e_r-e_s.
3. No vertex vector z satisfies z(s)=1,z(r)=0 and B^T z=0.

**Proof.** A path gives a signed unit flow. If no path exists, the indicator
of s's component is such a z; pairing it with Bf would give both0 and-1,
contradiction. This proves all directions and supplies either a path-flow
witness or an exact separating cut. A cycle vector lies in ker B. Adding
cycles to any receipt flow does not change its endpoint obligation and
cannot create e_r-e_s from zero.

The negative conclusion is **no path in this stored bank**, not that s
diverges or lacks another proof. New authenticated edges or a genuinely
new certified terminal may change the component.

### Export retains the two clocks

An oriented receipt contributes odd-clock pair (a,b) with
U^a(x)=U^b(y). Across the same intermediate integer, the inherited
[clock-composition rule](collatz_clock_holonomy_20261007.md) is

    (a,b)*(c,d)=(a+max(c-b,0), d+max(b-c,0)).          (1)

A path to r=1 gives U^A(s)=U^B(1)=1. Replaying at most A actual odd edges
therefore exports a strict first-hit word, cutting off any formal root
padding. This is a bounded certificate conversion, not an unbounded ROOT
search. The entry clocks cannot be replaced by their difference alone.

For a simple path with at most V-1 edges, let M bound both word lengths
on every stored edge. Equation(1) gives A<=(V-1)M. This explicit finite-bank
deadline comes from the authenticated path labels, not from capacity alone.

The small checked example is

    7 -- 3 -- 5 -- 1,

where 7 and3 meet at5 via words1123 and1. The accumulated clocks are(5,0)
and the exported ROOT word is11234. Reversing the first edge and using the
complete checked route from7 gives3--7--1 with clocks(2,0), exporting14.
The public `export_from_receipts` returns None only when no stored path
exists, never as a divergence verdict.

## 3. A positive capacity is exactly a finite grounding certificate

Assign each edge a positive rational conductance c_e. This is an explicit
proof-network weight, not an imported Collatz invariant. Write

    E(h)=sum_e c_e (h(v)-h(u))^2,
    L=B diag(c) B^T.

For distinct terminals s,r define

    C(s,r)=min { E(h): h(s)=1, h(r)=0 }.               (2)

There is no assumption that the graph is connected. Components meeting
neither terminal have a harmless additive gauge. The minimum exists in
finite dimension and is unchanged by those gauges.

**PROVED iff.** C(s,r)>0 exactly when s and r are connected. A disconnected
source-component indicator has zero energy. On a connected finite path,
zero energy forces equal values at its two ends, so it cannot meet the
boundary values in(2). Equivalently Cauchy--Schwarz on that path gives an
explicit positive bound.

More generally every signed unit flow from section2 supplies

    C(s,r) >= 1 / sum_e f_e^2/c_e.                    (3)

Indeed <Bf,h>=-1; weighted Cauchy--Schwarz gives
1<=E(h) sum f_e^2/c_e. Conversely the harmonic potential's electrical
current, divided by C, is a unit flow attaining equality. Thus a flow,
capacity margin, and ROOT-connected component are three exact finite
carriers for the same fact. Computing one does not independently prove
that an unknown source has a path.

### Exact residual energy and the condition a local estimate must meet

Let h be any rational trial vector with the two pinned boundary values.
Let I be the other vertices and rho=(Lh)_I. Solve

    L_II y_I = rho,   y(s)=y(r)=0.                   (4)

This system is always compatible. Its kernel consists of constants on
components disjoint from the two terminals; rho pairs to zero with those
constants. Different solutions differ only by such a gauge. Then

    h_*=h-y,
    C(s,r)=E(h)-rho^T y,
    rho^T y=E(y)>=0.                                (5)

**Proof.** Equation(4) says h_* is harmonic at every interior vertex.
For any vector z vanishing at both terminals, <z,Lh_*>=0. Expanding the
quadratic energy therefore shows that h_* minimizes(2) and that the gap
between the trial and minimum is E(y). Pairing(4) with y gives the last
identity. Floating components do not affect any term.

Define the residual's dual-energy norm on those zero-boundary vectors by

    ||rho||_* = sup |<rho,z>|/sqrt(E(z)).

The ratio is evaluated on positive-energy vectors, and rho annihilates
the zero-energy kernel. If no positive-energy vector exists, this norm
is defined to be0. Equation(4) and Cauchy--Schwarz give
||rho||_*^2=E(y). An independently checked **global** estimate
||rho||_*<=epsilon consequently proves

    C(s,r)>=E(h)-epsilon^2.                           (6)

A strictly positive right side certifies a stored route. A scalar or
coordinatewise estimate is not the hypothesis of(6). The symmetric
identity here is the exact finite counterpart of the source paper's
dual-energy requirement; it does not erase the adjoint condition in its
nonsymmetric setting.

**Small local residual hostile.** Take an L-edge unit-conductance chain
from s to an ungrounded leaf, with ROOT isolated. Put h(s)=1 and let h
decrease linearly to0 at that leaf. The sole nonzero interior residual is
-1/L at the leaf. Nevertheless

    E(h)=1/L,   ||rho||_infinity=1/L,
    ||rho||_*^2=1/L,   C(s,r)=0.

Replacing the dual bill by the square of the coordinatewise residual
would falsely produce the positive floor1/L-1/L^2 for L>1. The missing
factor is the graph's long-range conditioning. This example also shows
why omitting a terminal defect can make a locally correct calculation
look grounded.

## 4. Interior elimination is lossless only with boundary and path data

For an interior vertex v with conductances c_vi and total d=sum_i c_vi>0,
minimizing its star energy sets

    h(v)=sum_i c_vi h(i)/d.

The minimized energy replaces the star by edge conductances

    c_ij <- c_ij + c_vi*c_vj/d,   i<j.               (7)

An isolated interior vertex can be discarded. Substitution proves(7)
exactly. Keeping s,r as protected ports makes capacity invariant under
successive elimination, regardless of order. This is the elementary
star-mesh/Kron version of the inherited Schur-patch mechanism.

Each new positive edge also retains a path witness through v, or the
original authenticated graph is kept so that a route can be recovered
there. An unlabelled reduced matrix is insufficient for a Collatz ROOT
word: it has forgotten the actual source seams and clock pairs. No terminal
can be silently eliminated, and no residual forcing can be declared zero
merely because its coordinate was hidden.

**Resource convention under refinement.** One conductance1 macro edge has
capacity1. Replacing it by L conductance1 edges gives capacity1/L, although
grounding is unchanged. To preserve its numerical capacity, retain the
total resistance: use L edges of resistance1/L, i.e. conductance L.
Likewise duplicating a proof record as a parallel unit-conductance edge
changes the numerical observable without adding a proof. Provenance and a
declared weight-allocation rule are needed when comparing capacities.

There is no uniform path-length deadline from a positive capacity floor
alone. L parallel, internally disjoint, unit-conductance paths of length L
have capacity1, while every terminal path has length at least L. A bound
on graph size and the actual edge clocks supplies the deadline in section2;
the scalar1 by itself does not.

## 5. Why geometric neighborhood repair is a different operation

The geometric-similarity source repairs its closed set R of missed centers
by an open neighborhood V. For each x in R and normalized t, the points
x+t q^n eventually belong to V because they converge to x. These later
indices are allowed by the stated hitting problem. That convergence and
open membership are the premises that upgrade an estimate to every center.

A fixed exact Collatz word of length r and valuation cost A has endpoint

    U_w(n)=(3^r n+B_w)/2^A.

Its native source cylinder is infinite, but equality of the endpoint to1
selects at most **one** integer. On n=n0+2^(A+1)t the endpoints change by
2*3^r t. For example the exact letter4 on5mod32 sends

    5+32t -> 1+6t.

The t=0 ROOT proof is not a ROOT certificate for the other members of its
open binary cylinder. The profinite sequence1+2^j also approaches1 while
never equalling1; it is a sequence of observations, not asserted to be one
Collatz orbit. Thus an open neighborhood cannot simply be declared a
terminal. A separately proved rooted family, with an authenticated parameter
decoder and a termination proof, would supply the missing permission.
The point/family distinction is inherited from
[finite-seed kernels](finite_seed_kernel_20261005.md).

The XY source similarly tracks a marginal coordinate of order1/j rather
than asserting a finite zero. The elementary recurrence
u_(j+1)=u_j/(1+u_j), u_0=1, has u_j=1/(j+1)>0 forever. A decreasing real
quantity needs an additional well-founded rank, a positive minimum decrement,
or a proved terminal rule before it gives finite exhaustion. Finally the
dimer source explains why a finite topological record does not record all
traversals; our clock pairs and expanded witness path are the corresponding
extra data, without asserting a map from dimer loops to integer orbits.

## 6. Relation to the restored partial-edge bank and the remaining task

The concurrent [terminal-basis package](collatz_terminal_basis_20261007.md)
owns the actual large-graph reconstruction and finite census. Its gain
comes from retaining authenticated partial observations, not from inventing
new edges with positive weights. In the reported frozen universe it grounds
16,258 of16,384 requests before the additional terminal certificates, and
identifies37 unresolved sinks. Those are results of that package, not a
second independently recounted census here.

Our small concrete unresolved control retains only its authenticated prefix
4591 ->2717873 of26 odd edges, plus isolated ROOT. Its capacity to ROOT is
zero. This says exactly that the prefix alone has an unpaid endpoint; it
does not question the validity of the prefix or classify4591's true fate.
An independently checked suffix of that sink legitimately adds a path.

For a finite receipt bank, the next constructive action is therefore to
retain its missing boundary sources, add authenticated routes or terminal
proofs, and recompute the grounded components. The capacity dual records
why a proposed local repair has or has not grounded a source. To make an
infinite compressed graph yield a theorem for every input, one must still
prove the complete edge/guard specification, boundary conditions and a
sourcewise positive certificate or finite extraction bound. The six trusted
headlines do not supply those arithmetic obligations.

## 7. Reproduction and declared universe

```text
python -B -X utf8 04-computation/experiments/paper_local_global_grounding_20261007.py
python -B -O -X utf8 04-computation/experiments/paper_local_global_grounding_20261007.py
```

The complete graph universe is every simple graph on2..5 labelled vertices
(1,098 graphs), with fixed distinct terminals, two rational trial potentials,
unit-flow/cut certificates and both forward/reverse interior elimination
orders. Additional exact controls cover disconnected/connected chains of
length1..32, resistance-preserving refinement, L parallel L-edge paths for
L=1..7, the actual7/3 ROOT splices, the stated26-edge unpaid prefix, and
14 malformed or ungrounded receipt cases. The synthetic cycle is explicitly
not asserted to be a positive Collatz cycle. All checks remain active under
`-O`. The reproduction does not use an unbounded orbit search, and its
capacity is not identified with the earlier source measure or Green weight.
