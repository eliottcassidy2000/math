# Sharp local deletion counts for every parameter type of the minimal AND

**Status: PROVED local sharpness, with exact finite circuit checks.** These
families show that AP1's three deletions and GP's two deletions at one child
cannot be increased unconditionally for the specified legal choices. Together
with [the AP2 four-deletion theorem and example](AP2_FOUR_DELETIONS.md), all
three minimal-AND parameter cases have sharp local gate counts. This is not a
claim that every elimination strategy must incur these counts, or that the
worst rows can accumulate without compensating losses.

The model and legal-state invariants are those of the pinned provisional
elimination manuscript linked in the AP2 proof. Write

    IP_M(a,b)=sum_(j=1..M) a_j b_j,  M>=1.

It has 2M essential inputs and a circuit of M AND gates and M-1 XOR gates,
hence exactly 2M-1 full-B2 gates. All sums below are binary.

## 1. AP1: three deletions, at minimum-size endpoints

Let n=4+2M and

    f(p,w,r,z,a,b)=pw+rz+IP_M(a,b).

Mark the pair (p,r), excluding corner (0,0); all other parameters are unpaired.
The affine map is the identity and there are no affine equations. The polar
form of f is nondegenerate on n=2(M+2) coordinates. Therefore every constant
affine flat has dimension at most d=M+2: its direction is totally isotropic
in this nondegenerate alternating space. Setting one coordinate of each
quadratic pair to zero attains dimension d. Thus the global affine-disperser
hypothesis is satisfied for this d.

The state has q=1 and safe-core dimension a_state=n-1=2M+3. Since M>=1,

    a_state >= d+2.

Use gates G=pw, H=rz, K=G+H, a usual circuit for IP_M, and output K+IP_M.
There are 3+(2M-1)+1=2M+3=n-1 gates. All are nondegenerate and live; no
preparatory paired or affine row applies. The paired parent p of G and its
partner r have distinct unique children G,H. The other parent w is unpaired.
Thus G is a legal AP1 choice.

Fix r=1 safely and then w=0 to kill G. The gates G,H,K disappear, and exactly
three original gates are deleted. The output becomes z+IP_M, with the usual
2M-gate circuit remaining. The old paired p becomes unused but is retained
as a parameter. The exact row data are

    (Delta s, Delta I, Delta q, Delta a_state)=(-3,-1,-1,-1).

The initial circuit is minimal even among circuits agreeing only on the
marked domain: all n coordinates are essential on legal same-fiber pairs.
For p or r flips, set the other marked coordinate to one and the quadratic
mate w or z to one. For w or z, set p or r to one and make the marked pair
legal. For an IP coordinate, set its mate to one and p=r=1. Thus at least n-1
gates are required. The residual z+IP_M has 2M+1 essential coordinates and
requires at least 2M gates, attained above. Neither endpoint contains an
avoidable gate.

## 2. GP with one child: two deletions, at minimum-size endpoints

Let n=2M+3 and

    f(x,y,z,a,b)=xy+z+IP_M(a,b).

Initially q=0, the affine map is the identity, and every parameter is used.
The polar form has rank 2M+2 and one-dimensional radical spanned by z.
Its quadratic function has value one on that radical vector.

A constant affine-flat direction U is totally isotropic. In dimension
2M+3 such a direction has dimension at most M+2. Any isotropic subspace of
that maximum dimension contains the radical: otherwise adjoining the radical
would produce an isotropic space of dimension M+3. But a direction containing
the z vector cannot support a constant affine restriction, since translation
by z always flips f. Hence every constant affine flat has dimension at most

    d=M+1.

This dimension is attained by fixing y=z=0 and one coordinate of each IP pair
to zero, while leaving x and the other M IP coordinates free. The global
affine-disperser hypothesis therefore holds. The loop guard is satisfied:
a_state=n=2M+3 >= d+2.

Use G=xy, K=G+z, a usual IP_M circuit, and output K+IP_M. There are
2+(2M-1)+1=n-1 gates, which is minimal because all n inputs are essential.
G has exactly one child K. All gates are live and nondegenerate, and no
preparatory row applies.

GP marks (x,y), excluding corner (1,1), and makes G constant zero on its new
three-point domain. Exactly G and K disappear. The residual circuit computes
z+IP_M in the minimum 2M gates. The exact row data are

    (Delta s, Delta I, Delta q, Delta a_state)=(-2,-2,+1,-1).

There is an immediate **global-accounting caveat**: x and y become unused
paired parameters. Therefore P0 applies next, removes the pair, and supplies
additional potential loss while preserving a_state. This family proves
sharpness of the individual GP row only. It is explicitly not a sequence
with zero compensating loss, and cannot refute a stronger amortized theorem.

## 3. Complete local gate-count conclusion

For a legal chosen minimal AND, the proved guaranteed deletion counts are

- AP1, one paired and one unpaired parent: three, attained above;
- AP2, two parents from distinct pairs: four, attained in the companion proof;
- GP, two unpaired parents with k children: at least 1+k, with k=1 attaining two.

The corresponding parameter changes matter. Gate counts alone do not measure
the total potential charge, and the examples do not certify arbitrary
transition sequences or eliminate reachability and strategy questions.

`check_all_and_rows.py` evaluates literal circuits and legal restrictions for
M=1,...,4, finds same-fiber essentiality witnesses, checks gate counts, and
verifies the polar ranks and radical value. Its finite checks support the
all-M proofs; they are not a proof-assistant build. Normal and optimized
Python runs must match `all_rows_receipt.json`.
