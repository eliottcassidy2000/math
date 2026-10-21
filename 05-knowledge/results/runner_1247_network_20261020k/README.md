# The 1247 clock, triadic colouring rigidity, and an exact contact triangle

2026-10-20. **PROVED:** the exact runner/colouring comparison below, its
all-parameter triadic extension, the marked voltage lift, and the displayed
electrical certificate. **FINITE-EXACT:** the declared controls in
[check.py](check.py) and [receipt.json](receipt.json). General LRC, LRC(14),
and plane six-versus-seven colouring remain **OPEN**. No literature-priority
claim is made. The general ellipsoid construction is inherited, not new.

## 1. Inheritance and the retained objects

The closest proved mechanisms are the repaired
[rotation-colouring boundary](../runner_area_clock_20261020j/boundary_and_colouring.md),
the [projected-cube star certificate](../projected_cube_ellipsoid_20261021/README.md),
its [one-pair triangle repair](../triangle_ellipsoid_repair_20261021/README.md),
and the [electrical contact certificate](../electrical_ellipsoid_certificate_20261021/README.md).
The canonical hostile is the now-refuted identification
`chi_f(G(Z,D))=1/kappa(D)`: the source `{1,3,4,7}` has respective values
4 and5. The current correction is recorded at the front of
[MISTAKES](../../../01-canon/MISTAKES.md) and in
[HYP-9350](../../hypotheses/HYP-9350-polygonal-order-kepler-flow-and-the-hadwiger-nelson-bridge-of-the-lonely-runner-conjecture.md).

An older correction also distinguishes the connection set `{1,2,4}` of the
seven-vertex Paley tournament from a tight runner set: its clearance is at
least1/3, not the four-runner threshold1/4. Here the additional speed7 is
retained, including the loop it creates under reduction modulo7.

Live board: **literal integer speeds / rotation equivariance / independent
sets / winding phase / weighted contact network / arithmetic time**.
The anchor is the exact four-speed comparison. The niche is a complete
rigidity and quantitative obstruction theorem for sources containing1,2.
The wildcard is a sharp three-state restoration of a tournament quotient.

Connection contract: an integer source D gives both a distance graph and
an ellipsoid carrier. The rotation map `n -> {nt}` preserves every forbidden
distance only when its single common clock satisfies the retained inequalities.
Unrestricted colourings forget that group-action condition. Electrical
potentials preserve the labelled pair capacities and covariance, but forget
the integer lattice and absolute time. The sidecars are respectively the
rotation clock, winding coordinate, and full source/coordinate gauge.

## 2. The two requested sets have different exact answers

For a finite nonempty set of positive integers D, put

    kappa(D)=max_(t mod1) min_(d in D) ||dt||,
    G_D: vertices Z, edges {n,n+d} for d in D.

The ordinary circular and fractional chromatic numbers use unrestricted
graph colourings. The inherited elementary inequalities are

    chi_f(G_D) <= chi_c(G_D) <= 1/kappa(D),
    chi_c(G_D) <= chi(G_D).                            (1)

The restricted rotation-colouring optimum is exactly `1/kappa(D)`.

For `{1,2,4,7}`, the pair1,2 bounds the clearance by1/3. Indeed both
`||t||>=1/3` and `||2t||>=1/3` force `t=1/3` or `2/3` modulo1.
Every speed in the set is a nonzero residue modulo3, so both clocks attain
1/3. The graph has the triangle `{0,1,2}`, and `c(n)=n mod3` is proper.
Thus all of the following values are exact:

| Source D | kappa | chi_f | chi_c | chi | restricted rotation optimum |
|---|---:|---:|---:|---:|---:|
| 1,2,4,7 | 1/3 | 3 | 3 | 3 | 3 |
| 1,3,4,7 | 1/5 | 4 | 4 | 4 | 5 |

The second row is the inherited corrected example: `{0,3,4,7}` is a K4,
and the period-eight colouring `2 floor(r/4)+(r mod2)`, `r=n mod8`, is
proper. Its exact five-interval rotation obstruction and witnesses at
`t=1/5,2/5,3/5,4/5` are independently replayed here. Nothing permits
replacing the first source's speed2 by3 while retaining its clock.

## 3. Complete triadic rigidity and a sharp quantitative family

**Rigidity theorem.** Let D be any finite positive integer source containing
1 and2. Every proper three-colouring of G_D, if one exists, is `n mod3`
up to a single palette permutation. Consequently

    chi(G_D)=3 iff 3 does not divide any member of D
               iff kappa(D)=1/3.                     (2)

Proof: every three consecutive integer vertices form a triangle. In the
next triangle the colour at n+3 must be the one at n. This forces the
entire bi-infinite colouring from one triangle. The resulting colouring
works exactly when no forbidden distance is divisible3. The assertion
about kappa follows from the two possible maximizing clocks of the pair1,2.
Compactness of the real clock makes the maximum attained, so failure of
both clocks means strict inequality, not just failure of a selected test.

There is a stronger fractional obstruction. If `3m in D`, m>=1, then

    chi_f(G_D) >= 3+1/m,
    kappa(D) <= m/(3m+1).                             (3)

Proof: in the induced window `{0,...,3m}`, an independent set has gaps at
least3. Having m+1 members would force exactly the points0,3,...,3m;
the first and last are adjacent by distance3m. Thus its independence
number is at mostm. Sum all vertex-cover constraints of a fractional
colouring to obtain the first inequality. For an independent direct
rotation argument, move an open arc of length kappa around the circle
containing the3m+1 points `{nt}`. Its contained labels are independent,
so it contains at mostm. Integrating the count gives
`(3m+1) kappa <= m`, including equality boundaries because arcs are open.

Both bounds are sharp for the whole family

    D_m={1,2,3m},   q=3m+1,   t=m/q:
    kappa=m/q,  chi_f=chi_c=q/m,  chi=4.              (4)

At this clock the three cyclic increments are m,2m,-m moduloq, all at
cyclic separation at leastm. Thus `n -> (mn mod q)/m` is a circular
q/m-colouring. The m-subset `{0,3,...,3(m-1)}` is independent moduloq;
its q translates, each with weight1/m, give an explicit fractional
colouring. Every vertex occurs exactlym times. The finite-window lower
bound gives equality. The ordinary period-q colouring

    c(n)=floor((mn mod q)/m)

uses four colours and is proper because each colour class occupies at most
m consecutive residues, spanning fewer thanm unit gaps. Three colours are forbidden by(2).

For a larger source D containing1,2,3m, (3) remains a lower/upper bound;
its chromatic number need not be4. The packet never promotes(4) to such
supersets. This is a family of exact source certificates, not a general
fractional-colouring identity or a new runner-count theorem.

## 4. The honest tournament map retains one winding coordinate

On the seven residues, orient `r -> r+d` for `d in {1,2,4}`. Exactly one
of each nonzero residue and its negative occurs, so this is an intrinsic
tournament with no arbitrary ties. The directed triangle

    0 ->1 ->3 ->0      has increments1,2,4.

Its integer lift is `0 ->1 ->3 ->7`, not a closed integer path. The
quotient has discarded one winding. The extra forbidden distance7 in
the requested source also becomes a loop at every residue. Deleting those
loops is a source change, so the loopless tournament cannot directly colour
G_1247. Nor does its undirected K7 certify seven colours for the integer graph.

Write the actual integer as `n=r+7z`, 0<=r<7. Retain z modulo3. The edge
of length d transports the marked state by

    (r,z) -> ((r+d) mod7, z+floor((r+d)/7) mod3).      (5)

This is exactly the quotient modulo21. The colour `(r+z) mod3` is the
original `n mod3`, because7=1 modulo3. It properly colours all distances
1,2,4,7. In particular, the displayed triangle ends at `(0,1)` rather
than returning to `(0,0)`.

The extra phase size3 is sharp **among the periodic lifts modulo7m**:
by(2), any three-colouring pulled back from such a lift must have period3,
hence3 divides7m, or3 dividesm. Conversely every such m permits it.
This is a restriction on this explicitly marked quotient family, not a
claim about every possible encoding or every tournament construction.

## 5. The same source lies one unit outside the electrical star cone

Use the inherited projected-cube carrier `P[-1,1]^4`, with P orthogonal
projection onto `(1,2,4,7)^perp`. The dominant-star condition would be

    7^2 >= 7(4+2)+4*2.

It fails by exactly1:49<50. Only the leaf pair2,4 violates it; the next
pair4,1 passes because49>=39. The inherited one-pair repair therefore
has normalized correlation

    c=1/[(7+4)(7+2)]=1/99.                            (6)

In grounded electrical coordinates
`y_i=x_i/v_i-x_7/7`, for the leaves1,2,4, its covariance is

    Sigma = [[64/49,     0,      0],
             [    0,81/196,  1/392],
             [    0, 1/392,121/784]].                 (7)

The exact positive-conductance graph is a triangle on2,4,7 with the leaf1
attached to7:

    lambda_(1,7)=49/64,  lambda_(2,4)=1/25,
    lambda_(2,7)=119/50, lambda_(4,7)=161/25.          (8)

All other conductances are zero. Its grounded Laplacian L satisfies
`L^-1=Sigma`. Every pair resistance obeys
`R_ij<=(1/v_i+1/v_j)^2`, and every positive edge has equality. Hence the
inherited electrical certificate proves global optimality, including
translated ellipsoids: every feasible covariance T has `tr(LT)<=3`, so
`det(T)<=1/det(L)` by arithmetic-geometric mean, with equality only at
Sigma. Tight positive contacts and connectivity force the centre to zero.

The exact outputs are

    det Sigma = 200/2401,
    sum_(spanning trees) product lambda_e = det L = 2401/200,
    Euclidean projected determinant = 1280/343.       (9)

The first determinant uses the grounded potential coordinate gauge; the
last uses an orthonormal basis of the projected space. They are not the
same volume normalization. Removing the2--4 contact leaves a connected
tree but violates a pair-resistance constraint. Thus connectivity or an
unweighted triangle count is not the certificate.

This optimum cannot come from any real ambient diagonal covariance.
In the normalized star coordinates, such a covariance has off-diagonal
entries `v_i v_j q_7/[(7+v_i)(7+v_j)]`. They are either all zero or all
nonzero. The optimum has only the2--4 entry nonzero. This is a concrete
source in the inherited theorem's strict higher-dimensional failure of
the diagonal class, despite its elementary three-colouring and clock.

The electrical cycle has no preferred orientation and is not a cyclic
runner/tournament obstruction. Its role is to repair a covariance
constraint. The integer clock witness is supplied separately by Section2.

## 6. Exact geometry still loses the absolute time gauge

Multiplying every speed by3 leaves v-perp, P, and the literal projected
cube unchanged. It also leaves kappa unchanged after reparametrizing time.
But the original optimal clocks1/3,2/3 become the six preimages modulo1
under `t ->3t`. The integer clock grid is not retained by the shape alone.

There is an even sharper finite denominator readout. At the five-runner
threshold1/5,

    N_min({1,2,4,7})=3,
    N_min({3,6,12,21})=5.                             (10)

The first source fails denominator2 and succeeds at1/3. For the scaled
source, denominators2,3,4 collide through6,3,12 respectively;1/5 succeeds.
Both denominators are finite by the displayed witnesses, not an appeal
to general LRC. The scaled source is intentionally nonprimitive: the
lost gcd/time normalization is precisely the sidecar being tested.

Thus an optimal ellipsoid and even its same unlabelled contact pattern
cannot be read as a named rational runner time without restoring the
original source and clock. No assertion is made that convex geometry is
unable to contribute after those data and an arithmetic hitting argument
are supplied.

## 7. Exact validation

The checker uses two independent clock readers: all critical rational
breakpoints/intersections of tent functions, and closed-interval intersection.
It retains singleton maximizing times. The complete bounded triadic source
universe consists of all128 subsets of `{1,...,9}` containing1,2. Independent
three-colouring enumeration covers every word of lengths3 through8. The
sharp family is checked for m=1,...,12; complete finite-window independence
enumeration is limited to m<=5 and is stated as such.

The electrical receipt validates exact types, source speeds, every resistance,
positive contacts and the inverse matrix; the spanning-tree expansion is an
independent determinant path. Wrong-source, Boolean and floating aliases are
rejected. The voltage checks retain each edge's actual source lift. None of
these finite controls substitutes for the quantified proofs above.

```text
python -B 05-knowledge/results/runner_1247_network_20261020k/check.py --write
python -B -O 05-knowledge/results/runner_1247_network_20261020k/check.py
```

The default reads and compares the immutable receipt; only `--write` saves it.
No orbit discovery, numerical optimizer, or floating-point comparison is
used. Normal and optimized replay both pass **1,765 checks**, with saved
LF SHA256
`e621094a5a324397c843fcbf3c8590091deb95d4d7dc3e0d233451347186dda9`.

Independent audit: root reviewed the quantified proofs and source boundaries;
tiling_fold independently audited the full proof and production APIs and
replayed the optimized checker against the saved receipt. Both report PASS.
The review includes the sharp all-m family, the winding sidecar and its scoped
minimal phase, the exact electrical contact certificate, and the nonprimitive
absolute-clock hostile. All three owned files are frozen after this review.
