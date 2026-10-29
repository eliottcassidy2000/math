# Schur boundary carriers for recursive tournament circulations

## Status and inheritance

**PROVED in the stated operation menu**, with exact finite checks. This is a new
repository synthesis of classical Cauchy/Littlewood–Richardson identities with
the older circulation-cone lane, not a priority claim about symmetric functions.
The source snapshot is `d3741fe6ed02f6072dbe6bfda2f3d04193d75ce0`.

The least-used inherited mechanism is
[`REGULAR_FLOW_RECIPROCITY.md`](../../regular_tournament_flow_reciprocity_20261021/session23b-reciprocity/REGULAR_FLOW_RECIPROCITY.md):
a regular tournament's nonnegative integral circulation series has inversion
shift equal to its number of arcs. Its companion
[`COROLLARIES.md`](../../regular_tournament_flow_reciprocity_20261021/session23b-reciprocity/COROLLARIES.md)
keeps fixed arc orbits and their actual integer weights. Those are raw flows
and orbit flows, not Hamilton paths. We recover the raw-flow object, retaining
its *charged boundary responses* so that tournament substitution becomes closed.

The earlier [flow-decoration packet](../../flow_decorated_tournament_series_20261021/session23b-flows/PROOF.md) already proves rational fixed cones and exact factorization over ordered strong components. Its tetrahedral 4+3 Fano support records lose oriented integral flow multiplicities, and its strong order-four cone is free on cycles of weights3,3,4. We do not reprove those results or confuse mod-two support with the present nonnegative integer carrier. The new step concerns *cyclic interchild attachment*, where cross-child flow need not vanish and SCC scalar factorization no longer applies. As in the inherited packet, forgetting a flow decoration means q=0, not the generally divergent q=1.

The other inherited model is the recursive ternary tournament family from
[Session94](../../session94_recursive_ports_crt_clocks_20261029/README.md).
The current operation is complete directed interchild attachment. It is not an
arbitrary named-port observation menu. Named individual probes require the full
Laurent boundary dictionary, discussed below.

Classical inputs: the Schur Cauchy identity and product/coproduct rules.
Primary references include [Molev, *Comultiplication rules for the double Schur
functions and Cauchy identities*](https://arxiv.org/abs/0807.2127), specializing
to ordinary Schur functions, and [Zinn-Justin, *Littlewood–Richardson
coefficients and integrable tilings*](https://arxiv.org/abs/0809.2392).
The identities are explicitly stated below; the verifier materializes both
sides independently in small finite alphabets.

## 1. The carrier

Let T be any finite directed loopless graph, with n named vertices, and let

    F_T(q; x) = product_(u→v in T) (1 − q x_u/x_v)^(-1).

The coefficient of q^d x^b counts nonnegative integral assignments to arcs
of total use d and signed divergence b. The ordinary circulation series is
R_T(q)=CT_x F_T(q;x). Here CT is a coefficient operation in the formal
q-series; every q coefficient is a finite Laurent polynomial.

For partitions λ,μ, define

    M_T(λ,μ;q) = CT_x F_T(q;x) s_λ(x) s_μ(x^(-1)).             (1)

This is zero unless |λ|=|μ|, by total homogeneity. A Schur polynomial with
more than n rows is zero. Let M_T[t] be the finite matrix indexed by partitions
of t with at most n rows. Entries lie in Z_≥0[[q]]. The empty-partition entry
is R_T. No division by an automorphism order occurs.

**What is kept:** incoming and outgoing boundary partitions are paired against
one and the same child flow. Replacing M by separate marginal responses would
forget their joint compatibility. The partition coordinates are sufficient for
complete bipartite attachments because of the exact Cauchy identity

    product_(i,j)(1 − z x_i/y_j)^(-1)
      = sum_α z^|α| s_α(x) s_α(y^(-1)).                      (2)

This can also be read through the classical matrix/paired-tableau
correspondence: both tableaux have the same shape. Keeping only one marginal
shape or replacing the two tableaux by unrelated sources is not the identity.

## 2. A trace recurrence for the cyclic three-child composition

Let C3[T1,T2,T3] orient every vertex of child i toward every vertex of child
i+1 (indices modulo3). Expand each cross kernel by (2). The local constant
term in child i is M_Ti(α_i,α_(i−1)). Homogeneity forces all three α_i to
have the same size t. Consequently

    R_(C3[T1,T2,T3])(q)
       = sum_(t≥0) q^(3t) tr(M_T1[t] M_T3[t] M_T2[t]).       (3)

The matrices are extended by zero to one common partition bank if child sizes
differ. The displayed ordering follows directly from the local factors;
for three identical children it simplifies to

    R_(C3[T,T,T])(q) = sum_(t≥0) q^(3t) tr(M_T[t]^3).        (4)

A coefficient through q^D uses only t≤floor(D/3) and finitely many child
coefficients. This is an all-order exact recurrence for specified tournaments,
not a new formula for unrestricted A000568.

Internal arc weights may be retained as independent formal variables in F_T.
Then the same formulas hold coefficientwise. Giving the three complete cross
blocks their own variables z1,z2,z3 replaces q^(3t) by (z1 z2 z3)^t. Arbitrary
per-cross-arc weights do not in general have the uniform Cauchy factorization;
then the full named Laurent kernel remains the exact fallback. Rank-one
separable cross weights admit the corresponding weighted alphabets, but should
not silently be replaced by the unweighted M.

## 3. The matrix carrier itself has a closed LR recurrence

Writing only a scalar recurrence is insufficient if its matrix inputs cannot
be generated. The following closes the carrier.

Let c^λ_(a1,a2,a3) be the iterated Littlewood–Richardson coefficient, so

    s_λ(X1+X2+X3) = sum_(a1,a2,a3) c^λ_(a1,a2,a3)
                                      product_i s_ai(Xi).

Let c^ρ_(a,α) be the ordinary Schur product coefficient. For equal children T,

    M_(C3[T,T,T])(λ,μ)
      = sum_(a,b,α) c^λ_(a1,a2,a3) c^μ_(b1,b2,b3)
          q^(|α1|+|α2|+|α3|)
          product_(i=1)^3 [ sum_(ρ,σ)
              c^ρ_(ai,αi) c^σ_(bi,α_(i−1)) M_T(ρ,σ) ].     (5)

All indices are partitions. Terms with too many rows vanish in the relevant
child alphabet. Equation (5) follows by splitting the two external Schur
polynomials, applying (2), multiplying the two Schur factors at each end,
and taking the child constant terms. The same proof permits distinct T_i.

**Formal finiteness:** to compute a coefficient through q^D, only
Σ|α_i|≤D contributes. At external degree r=|λ|=|μ|, every local partition
has size at most r+D. Every LR splitting is finite. Induction over a finite
substitution tree therefore terminates. This does not assert bounded state
size independent of degree, a fixed-order recurrence in total flow, or a
measured complexity advantage over all existing flow enumerators.

For any fixed quotient tournament Q, replace the three α_i by one partition
α_e per quotient arc. At child i multiply all outgoing Schur factors into the
positive boundary and all incoming ones into the negative boundary. Repeated
LR multiplication and (1) again close the recurrence. Thus a tournament
substitution grammar can carry an exact flow inventory in addition to its
object counts. This is a decoration of a supplied grammar, not a solution to
finding all prime tournament representatives.

## 4. A genuine root-symmetry recurrence, with its scope

For C3[T,T,T], identify corresponding child vertices and let r rotate the
three children. An r-fixed flow has the same internal child flow in all three
copies and the same n×n cross-flow matrix on each of the three cross blocks.
Every arc orbit therefore costs three units. Collapsing the balance equations
leaves one child plus a complete ordered-pair cross matrix, including its
n diagonal entries. Cauchy gives

    R_fixed(r)(q) = sum_(t≥0) q^(3t) tr(M_T[t](q^3)).         (6)

The quotient by this particular root C3 action has series

    (R_(C3[T,T,T])(q) + 2 R_fixed(r)(q))/3.                 (7)

This is not the quotient by the entire recursive automorphism group. The
latter needs all its fixed edge-orbit constraints. With formal source weights,
(7) is a weighted orbit formula only when those weights are invariant under
the root action; (6) then cubes repeated internal weights as well as q.
Noninvariant source mass must not be silently replaced by a Burnside average.

## 5. Exact low-degree consequences in the recursive ternary family

Put T0=K1 and T_(k+1)=C3[T_k,T_k,T_k]. Let N=3^k, and let c_j(k) count
directed simple j-cycles without a chosen starting point. For j=3,4,5 these
are exactly the q^j coefficients of R_Tk: a nonzero tournament circulation
has weight at least3, so two nonempty circulations cannot occur below degree6.

For the child n=3^(k−1), put m=n(n−1)/2 and a2=n((n−1)/2)^2. The one-box
boundary entry begins

    M_T[(1),(1)] = n + m q + a2 q² + O(q³).

Indeed the q coefficient chooses one arc; the q² coefficient chooses a
length-two directed path, counted by Σ_v indeg(v)outdeg(v). T_k is regular.
Using (4), with no t≥2 contribution below degree6, gives

    c3(k)=3c3(k−1)+n³,
    c4(k)=3c4(k−1)+3n²m,
    c5(k)=3c5(k−1)+3n²a2+3nm².

Starting from zero, these solve to

    c3(k) = N(N²−1)/24,
    c4(k) = (N⁴−N)/52 − (N³−N)/16,
    c5(k) = (N⁵−N)/160 − (N⁴−N)/26 + (N³−N)/16.           (8)

At N=9 the counts are30,81,162; at N=27 they are819,8991,70470.
These exact family formulas are not universal formulas for arbitrary regular
tournaments. At degree6, cycle multisets and compatibility begin to matter;
the matrix recurrence handles them without identifying a flow with an
arbitrarily chosen cycle decomposition.

## 6. Hostile controls and the boundary of reciprocity

Take T=T2, the directed two-vertex chain. Its scalar R_T=1. Nevertheless

    R_(C3[T2,T2,T2]) = 1+8q³+12q⁴+6q⁵+41q⁶+O(q⁷),
    R_fixed(r)       = 1+2q³+5q⁶+O(q⁷).

The root-C3 quotient has four triangle-flow orbits, not one. The eight raw
triangles arise from independently choosing one of two vertices in each
child. Treating the children as scalar empty-flow objects and using
1/(1−q³) misses those ports. This disproves that particular scalar-only
composition; it does not prove that order plus scalar R can never determine
an attachment in any restricted family.

A stronger **FINITE-EXACT** witness keeps even the entire scalar series and
child order unchanged. Let T=C3⊕T2, on vertices0,...,4 with arcs

    (1,0),(2,0),(3,0),(4,0),(2,1),(3,1),(4,1),(3,2),(2,4),(4,3).

Both T and its opposite have order5 and exact scalar series1/(1−q³): all
circulation stays in the unique cyclic strong component. Yet

    [q³] M_T((3),(2,1)) = 878,
    [q³] M_Top((3),(2,1)) = 874.

Reversing all arcs transposes M. The coefficient comparison can be verified
by enumerating the220 weak compositions of3 across the10 arcs and the explicit
Schur monomials. Consequently (order, entire scalar circulation series) does
not determine even this symmetric boundary carrier.

The actual cyclic attachments C3[T,K1,T2] and C3[Top,K1,T2] have identical
raw coefficients through degree11 but degree12 coefficients9470 and9462.
Their common prefix through degree11 is

    1,0,0,11,25,34,96,270,573,1127,2344,4843.

Equation (3) computes these from child boundaries, and an independent direct
cross-matrix enumeration checks the assembled counts. Thus this witness is
also an attachment failure, not merely a distinction in unused matrix data.
The exact coefficient decomposition locates the first difference at cross
size t=3 and partition(2,1): the internal response difference starts at4q³,
and the T2 pairing has constant term2. Nine cross-edge uses plus three
internal uses therefore produce8q¹². Every other retained term contributes
zero difference through degree12. This explains the observed11/12 boundary
without asserting a universal attachment threshold.

The old regular-flow inversion law also needs a guard. On the directed C3,

    R_C3(q)=1/(1−q³),
    M_C3((1),(1);q)=3/(1−q).

For the second identity, any boundary flow consists of a directed segment of
length0,1,2 plus a nonnegative number of full cycles. There are three choices
for every total length. Thus the raw scalar has inversion shift3, whereas
this boundary response has shift1:

    M_C3((1),(1);q^−1)=−q M_C3((1),(1);q),

not −q³ times itself. A raw cone reciprocity theorem must not be copied to
all charged responses. The boundary sum changes the affine chambers and
inversion expansion. This is an explicit corrected near miss, not a
retraction of the inherited zero-divergence theorem.

Finally, (1) deliberately tests symmetric boundary polynomials. Arbitrary
named-port operations can distinguish more than these matrices. The exact
uncompressed carrier is [q^d x^b]F_T, with divergence vector b still named.
The Schur compression is certified only for the complete-block/LR menu.

## Verification

Run `python check_schur_flows.py` and `python -O check_schur_flows.py`.
The standard-library verifier builds arc-flow kernels by independent geometric
factor multiplication, Schur polynomials by semistandard tableaux, and LR
products by triangular Schur decomposition. It checks Cauchy in one/two
variables through degree4; trace composition for children K1,T2,C3;
full carrier closure for external partitions through degree2 over K1 and
degree1 over T2; literal root-rotation fixed flow equations through degree12;
Burnside integrality; and the low-cycle recurrences. These are finite checks
of the proofs, not an exhaustive tournament classification.
