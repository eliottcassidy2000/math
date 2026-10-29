# Phase2: distinct K3 surfaces and fifteen disposable secant bits

Status: **FINITE-EXACT + mathematical deductions from standard good-reduction
and elliptic-surface facts**. Dependencies: Python3 and SymPy1.14.0; no Sage.
All paths resolve relative to the script. Parent first-packet files are frozen.

## 1. Three named fibrations are not one Q-isomorphism class of K3 surfaces

The formulas attached to ICARM273,302,4239 define Weierstrass fibrations with
coefficient weights2,4,6,8,12. `check_fibration_geometry.py` computes
c4,c6,Δ of degrees8,12,24, proves c4³−c6²=1728Δ and that Δ is squarefree and
coprime to c4. Thus all24 finite singular geometric fibres are nodal I1;
there is no singular fibre at infinity. The minimal elliptic surfaces have
χ(O)=2 and baseP1, hence are K3 surfaces. This does not establish their Picard
rank or a17-dimensional Mordell–Weil group.

`count_k3_reductions.py` admits a prime only if the degree24 discriminant
remains degree24 and squarefree modp. For p>3 this gives a smooth proper
elliptic-surface model with the same fibre geometry. Counting every fibre,
including infinity and its identity point, gives:

|p|273 surface|302 surface|4239 surface|
|---:|---:|---:|---:|
|233|58509|58404|58374|
|239|61593|61620|61326|
|251|67512|67620|67521|
|263|73920|73884|74025|
|269|77493|77196|77166|

`verify_k3_counts.py` independently counts y-values directly using finite-field
frequency tables, without quadratic characters; all15 counts agree.

**Deduction:** these three K3 surfaces are pairwise non-isomorphic overQ.
A Q-isomorphism would identify their ℓ-adic cohomology, so smooth-proper good
reduction would force equal finite-field point counts. Already p233 excludes
each pair. This is stronger than showing different elliptic fibrations, but
we do **not** claim non-isomorphism overQbar: different rational forms of one
geometric surface can have different Frobenius traces. In particular, the
273 and302 formulas are not alternative Q-fibrations of the same Q-surface
as4239. Identification with named geometric lattices1092/1152/etc remains open
for273/302 here.

Sources of the equations:
[273](https://elliptic-rank.icarm.cloud/curve/273),
[302](https://elliptic-rank.icarm.cloud/curve/302),
[4239](https://elliptic-rank.icarm.cloud/curve/4239).

## 2. Three generic named sections proved independent, not seventeen assumed

For the common equation

    y²+(Lx+B)y=x(x+D)(x+E), B=pE+qD+pq(L−p−q),

the sections

    P=(0,0), Q=(−D,0), R=(−pq,p(pq−E))

are exact polynomial solutions. The additional obvious section(−E,0) lies
on the same horizontal line asP,Q and therefore is their negative sum.
ForR, y=p(pq−E) and y+Lx+B=q(D−pq), so their product equals
−pq(D−pq)(E−pq), the required cubic.

`check_named_sections.py` specializes these three sections at each record's
parameter, transports them into the certified minimal curve, and independently
proves their three character columns independent. A generic integer relation
would specialize to a relation, so all three triples are generically
independent. This establishes rank≥3 for each named fibration. It does not
upgrade the subset to generic rank17 merely because the specialized curve
has rank≥30.

## 3. Reproduce the classic rank17 surface, while keeping it separate

The primary [Elkies2026 paper](https://arxiv.org/html/2608.25406v1) supplies
17 polynomial sections and their height Gram matrix for the classic
rank17 fibration of determinant948. It states273 and302 are not in that
original fibration. That statement alone does not identify their surfaces.

`elkies_948_mathdata.json` retains the mathematical coefficient arrays, source,
and published sign choices. `check_elkies_948.py`:

1. reconstructs all17 y-polynomials from exact squares;
2. proves every section identity;
3. checks the squarefree degree24 discriminant;
4. recomputes all height pairings from homogeneous section intersections;
5. recovers the published17×17 Gram matrix, determinant948, and positive
   leading principal minors;
6. verifies all15 quadratic secant slopes relative toP1.

For this no-reducible-fibre K3 model, integral sections have height4 and
h(P,Q)=2−deg gcd(xP−xQ,yP−yQ), counting the homogeneous intersection at
infinity. We use this standard height formula as proved in the source;
we independently replay all polynomial and matrix arithmetic. The exact
nonzero determinant proves these17 generic sections independent. The known
upper bound17 for a Q-rational ellipticK3 fibration is cited, not reproved.
Nothing here transfers the determinant948 certificate to theX1152 surface.

## 4. A precise fifteen-bit quotient theorem

Let P1,…,Pn be a basis in a free abelian group. For any subset of indices
2,…,n−1, replace Pi with −P1−Pi, fixing the other vectors. The transition
matrix has determinant±1, so every choice is again a basis of the same group.
For a pairing G, its matrix changes by integral congruence UᵀGU; determinant
and positive definiteness survive, but named inner products need not.

On the classicK3, each of15 quadratic secant slopes throughP1 determines
Pi only up to exactly that involution. The line's third intersection is
−P1−Pi. If its slope is m, its other coordinates are

    x' = m²−x1−xi,   y' = y1+m(x'−x1).

`check_secant_bit_quotient.py` verifies all15 alternate section identities,
that no alternative equals its original, and that the named secant is the
same. It replays full height matrices for19 masks (all singleton changes,
none/all, and two alternating masks). The unimodular argument proves the
statement for all32768 independent choices, not a numerical extrapolation.

**Information conclusion:** fifteen branch bits are unnecessary to certify
rank and lattice determinant, but necessary to identify this ordered tuple
of named sections from the slope data. This is a concrete arithmetic version
of the distinction between a short decision certificate and a faithful
carrier. It does not establish a relation with the15 supersingular primes
or15 symplectic Lagrangians solely from cardinality.

## 5. Hostile audit across66 high-rank point lists

`highrank_inputs.json` is a compact public snapshot of66 point lists with
claimed rank≥30 (ICARM attribution retained). The saved local-character
certificates are checked by `verify_highrank_snapshot.py` in132458 exact
checks. For64 lists, the matrix has full column rank. For3567 and3568 only
rank17 was found among good primes below20000, although each lists30 points.
This certifies an independent17-dimensional image, **not a rank upper bound**.
A nonsaturated even-index generator list can lose rank modulo2; halving or
saturation may repair it. Absence of a full direct-character certificate does
not refute the database's claimed independence. These two lists are reserved
as hostile controls for any claim that the mod2 certificate is always complete.

## Reproduction

Run all check/verify scripts in this folder. Every completed normal run and
`python -O` run agrees. The downloadable `elkies_source.html` is a local reading
cache, not a publication artifact. The compact mathematical data and exact
certificates are the intended research inputs.

Open obligations: certify the generic rank of the273/302/4239 fibrations
beyond the three named sections; classify geometric rather than Q-isomorphism;
recover the missing mod2 information for3567/3568 through exact halving;
prove any broader family-height statistical statement rather than extrapolate
from record selection.

## 6. Closed the two hostile cases by changing the certificate prime

**The outstanding certificate gap in§5 is now closed without halving.**
For3567 and3568, use observation maps modulo3 instead of modulo2.
At a good prime p≠3 with v3(#E(Fp))=1, multiplication by#E(Fp)/3 maps the
reduction group to its cyclic order3 subgroup. Naming a nonzero generator
identifies this subgroup withF3 and gives a homomorphism E(Q)→F3.

Thirty such observations give full rank30 for each point list. The largest
primes used are563 and821 respectively. Separate good reductions of orders26
atp17, and28 atp19, exclude rational3-torsion. The same coefficient-divisibility
argument as for the binary certificate now uses repeated divisibility by3.
Thus each list is independent modulo rational torsion, with no BSD or GRH.

`check_odd_primary_certificates.py` constructs the certificates with the native
Weierstrass group law. `verify_odd_primary_certificates.py` independently:

* transforms to short Weierstrass coordinates;
* counts finite groups using direct square multiplicities;
* uses a different scalar-multiplication direction and short-form group law;
* performs dense F3 elimination;
* verifies the torsion-excluding reductions;
* rejects duplicate-point, wrong-group-order, and wrong-scalar hostiles.

It performs7876 checks, agrees under `python -O`, and rejects all6 hostiles.
Consequently **all66 provided high-rank point lists have independent exact
rank certificates**:64 by F2 observations, two by F3 observations. We do not
assert that the two original point subgroups actually have even index, since
that would require a saturation argument; it was only a possible explanation
for insufficient F2 data. Nor is any rank upper bound established.

This is a concrete source-preserving re-routing theorem: change the observing
prime while retaining the named curve, points, finite-reduction group orders,
subgroup generators, and torsion guards. A modulus-dependent information loss
is not an intrinsic rank obstruction. All2005 claimed independent point
directions across the66 lists are now certified. The next research question
is how to choose a cheap observation prime adaptively, not to force all
certificates through2-adic saturation.
