# Odd-prime unramifiedness for the graph Kummer carrier

2026-10-26. **PROVED deduction; independent audit PASS.** The scope is the
actual compact-support cohomology carrier of the raw graph character sums.
No pure decomposition, rank, Euler factor, or conductor exponent at 2 is
identified. The integral construction below closes the support gate; it does
not close those other gates. The [independent audit](../certificates/AUDIT_CONDUCTOR_SUPPORT.md)
records the reviewed proof hashes; this promotion changes status and adds
the audit link without changing the mathematical proof.

## Statement and inherited carrier

The six graph sums in [EXTENSION_TWIST](../motives/EXTENSION_TWIST.md) have
coordinates x_0=0, x_1=1, and five variable coordinates. Four base coordinates
must be distinct; the three cone coordinates may coincide with nonneighbors.
Every summand is the quadratic character of an oriented product of coordinate
differences. The orientation contributes only a constant sign.

More generally, take finitely many labelled affine coordinates, two fixed at
0 and 1, a set of required pairwise inequalities, and

    f = eta product_(i,j) (x_j-x_i)^(e_ij),
    eta in {1,-1}, e_ij nonnegative integers.

Let U be the open locus satisfying the required inequalities and f != 0,
over S=Spec Z[1/2]. Let V -> U be the finite etale double cover y^2=f.
For every prime ell, every H_c^i(V_Qbar,Q_ell), and hence its minus eigenspace
for y -> -y, is unramified at every odd prime p != ell. The same assertion
holds for H_c^i(U_Qbar,Q_ell).

The weighted finite-field sum is the alternating Frobenius trace of the
minus eigenspaces: at each U-point the cover has 1+chi(f) points. Assigning
chi(0)=0 makes the omitted graph-edge collisions contribute zero. Thus this
is a carrier for the raw sums, with their actual allowed collisions retained.
The alternating trace does not itself specify an individual pure constituent.
Any actual Galois-stable subquotient of this cohomology inherits the same
unramifiedness. If a proposed pure compatible constituent is subsequently
identified there, its conductor has no odd-prime support. Its exponent at 2
and the existence of that constituent remain separate questions.

## Integral construction

**1. Separate collisions.** Stratify U by the equality partition of all
coordinate labels. A partition that identifies 0 with 1, a required unequal
pair, or an edge with positive exponent is inadmissible. On any remaining
stratum there are r variable blocks, all distinct from one another and from
0,1. The stratum is the ordered configuration space M_(0,r+3), with the
three fixed markings 0,1,infinity. Merged edge factors acquire integer
multiplicities; reordering their differences changes eta only by a sign.
This description is integral and does not depend on the residue characteristic.

**2. Compactify integrally.** Put X=bar M_(0,r+3) over Z. Knudsen's Theorem
2.7 gives smoothness, properness, and relative normal crossings of the
boundary. In genus zero this is a scheme: bar M_(0,3)=Spec Z and the universal
stable n-pointed curve is bar M_(0,n+1), inductively a scheme. The boundary
has globally labelled components, indexed by stable partitions; their
intersections are smooth proper products of smaller genus-zero moduli
spaces (or disjoint unions of such products). Consequently the boundary is
relative *simple* normal crossings, as needed here. These statements also
cover r=0, when the boundary is empty. See [Knudsen, printed p179, Theorem
2.7](https://www.dam.brown.edu/people/mumford/alg_geom/papers/1976b--DetDiv%2BKnudsenSequel-GT.pdf).

The radicand is a unit on the entire integral open stratum. Its divisor on
X is therefore supported on boundary components, with no vertical component.
The generic point of every vertical fibre lies in the open stratum, so its
valuation is zero. Write

    div(f)=sum_D m_D D.

The integers m_D may be negative. Every boundary component at infinity is
included. Odd pole orders, just like odd zero orders, are branch components.

**3. Separate odd branch components.** Over S, repeatedly blow up the
intersection of any two distinct boundary components whose m_D are odd.
Such a center is smooth over S and has relative codimension two; if
disconnected its components may be blown up together. In local SNC coordinates
the center is (x_1,x_2), and the generic exceptional multiplicity is m_1+m_2,
which is even. The two selected strict transforms become disjoint. No new
intersection between previously disjoint odd components is introduced, and
no exceptional component joins the odd list. Hence at most one blowup for
each original odd pair suffices. Blowups of these SNC strata preserve the
smooth proper model and its relative SNC total boundary.

Call the result X'. Its odd boundary components are now pairwise disjoint.
For the pulled-back divisor write div(f)=2A+B, where B is their reduced sum.
This is valid also for negative m_D: use m_D=2 floor(m_D/2)+(m_D mod 2).
The principal-divisor relation gives L=O(-A) with L^2 isomorphic to O(B).
The associated double cover

    Y=Spec_(X')(O_(X') direct-sum L^(-1))

is finite over X' and proper over S. Its local equations are z^2=u or
z^2=u x_1, with u a unit. Since 2 is invertible and B is a disjoint smooth
divisor, these charts are smooth over S. Etale-locally one may take a square
root of u; the second chart then replaces x_1 by z^2. This also proves that
the reduced preimage E of the *entire* boundary is relative SNC. The preimage
of each labelled boundary component, and every intersection of these
preimages, is smooth proper over S, possibly disconnected. Y minus E is
the original double cover of this collision stratum. A square radicand may
give a disconnected cover, which causes no difficulty.

**4. Pass to compact supports over the integral base.** Work over
S_ell=Spec Z[1/(2 ell)] with coefficients Lambda=Z/ell^a. For j:Y-E -> Y,
the augmented boundary complex

    j_! Lambda -> Lambda_Y -> direct-sum Lambda_(E_i)
                  -> direct-sum Lambda_(E_i intersect E_j) -> ...

is exact in the usual augmented sense: at a geometric boundary point it is
the augmented simplex complex on the boundary components through that point;
off the boundary it is Lambda -> Lambda. It thus computes j_!Lambda from a
finite complex of constant sheaves on smooth proper S_ell-schemes.
[Proper-smooth base change, Stacks Proposition63.15.2](https://stacks.math.columbia.edu/tag/0GLG)
makes their higher direct images finite locally constant. Kernels, cokernels
between such sheaves, and extensions remain finite locally constant by
[Stacks Lemma18.43.5](https://stacks.math.columbia.edu/tag/093U).
Thus R^i(Y->S_ell)_* j_!Lambda is finite locally constant.
[Proper base change, Theorem59.91.11 and Lemma59.91.12](https://stacks.math.columbia.edu/tag/095S)
identifies its geometric fibres with the required compact-support cohomology.
Passing through all a and then tensoring with Q_ell proves unramifiedness for
this stratum.

Finally order the finite collision stratification by the number of blocks,
so the successive closed/open localization triangles glue its compact-support
complexes over S_ell. The same locally-constant closure argument proves the
claim for V itself. Equivalently apply a common proper compactification to
these localization triangles; proper base change identifies their fibres.
The identical argument without the double cover proves the U assertion.
The deck minus eigenspace is a direct summand over Q_ell, including ell=2.

The relative-base argument is essential: merely observing unramified graded
pieces in a spectral sequence of generic Galois representations would not
exclude a ramified extension. Here all terms and extensions exist as locally
constant sheaves over S_ell before restricting to its generic point.

## Boundaries, hostile controls, and remaining gates

The carrier here has freely varying labelled coordinates subject only to
equalities and inequalities between them. It is not the norm-seven hourglass
line arrangement. The inherited [reduction result, R2](../../normseven_arithmetic_f1_20261026/reduction/RESULT.md)
has genuine new triples 136 and 245 in characteristic 3, and the four marked
quotient kernels cease to act freely there. Restricting coordinates by linear
relations or specializing extra markers can introduce such odd-prime bad
fibres. Transferring the present support conclusion to that different carrier
requires an actual cohomological map and a good integral model; a shared
seven-label count does not supply either.

- Infinity parity is retained by div(f) on the proper model. For f=x on
  P1 both 0 and infinity are odd branch points; forgetting infinity fails.
- At z^2=x_1 x_2 x_3, blowing up (x_1,x_2) gives exceptional multiplicity 2,
  not 3: x_3 does not contain the generic center. Blowing up a triple center
  instead can create an odd exceptional component, so that is not the
  terminating algorithm used here.
- The sign restriction matters. Replacing eta by 3 already gives the
  zero-dimensional cover Spec Q(sqrt(3)), ramified at 3. Arbitrary rational
  constant twists require enlarging the allowed prime set.
- The fixed markers matter. Prescribing an extra marker at 3, rather than
  letting it vary, can introduce collision with 0 in characteristic 3.
  The theorem concerns coordinate differences with only fixed 0 and 1.
- Collision partitions preserve the actual graph product, including repeated
  factors after mergers. They do not assert that a raw cone sum is a fully
  injective seven-vertex count. Exact merge corrections are still required.
- This is an Artin-ramification support statement for cohomological Galois
  representations. It is unrelated to the conductor ideal measuring the
  difference between a singular ring and its normalization.

**Transport ledger:** source = literal graph Kummer cover; target = integral
smooth proper covers with SNC boundary, one for each collision stratum;
map = stratification, boundary blowups, and double-cover normalization;
preserved predicate = the open cover and its compact-support cohomology;
lost information = none about the carrier, but no canonical pure splitting
is produced; needed sidecar = collision multiplicities, orientation sign,
all boundary valuations, and cohomological degree. The cheapest falsifier
is an odd-prime degeneration of one declared model, not a fitted conductor.

The imported inputs above were read at the cited theorem statements on
2026-10-26. The blowup construction and its application to these graph sums
are deductions given here, not assertions quoted from those sources.
