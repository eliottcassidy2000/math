# Two residue measures give a stronger fourth-coefficient floor

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
The [full independent audit](continuing3_20260906_residue_floor_audit.md)
accepts the proof and exact companions. Root separately reconstructed the
formal residue series, determinant identities, complete floor chain,
bounded-support hostile and substituted tail envelope. The manifest records
the producer's pre-audit report identity.

    e4 > 161875/888583 > 9/50.                            (1)

Combining this floor with the unchanged original-root response identity
proves Q(-s)<0 at every original positive phase s>=2500. Neither the floor
nor the cutoff is claimed optimal. The remaining finite-phase question and
general actual Laurent noncancellation remain OPEN.

## 1. Inheritance, source, and the weak-root scope

The incoming proved
[residue-envelope theorem](long_frontier_sep06_anchor.md)
supplies the C/B moment representation and its first coefficient bounds.
The earlier
[effective-anchor theorem](continuing2_20260906_effective_anchor.md)
proves e4>1/100; the
[finite-phase theorem](continuing3_20260906_laurent_finite_phase.md)
closes the smallest phase and gives a tail at 75000. The new floor below
does not depend on the earlier floor. It uses one additional C moment,
the separate D measure, and the common bound on their support.

The canonical hostile is the finite-phase theorem's coefficient surrogate:
it has both anchors, all displayed Newton inequalities, and four positive
first phases, but positive full response at one phase. The corrected near
miss is to identify truncated Hankel positivity with the actual bounded
beta-root predicate. The least-used sidecar is the upper endpoint of the
residue supports, before discarding their common beta-root origin.

Let 0<=a_1<=...<=a_5, with sum 13 and square sum 59. Write x=e3, y=e4,
z=e5, so e1=13 and e2=55. Set

    B(v)=v^5-13v^4+55v^3-xv^2+yv-z,
    C(v)=v^4-12v^3+45v^2-(2/3)xv+(3/7)y,
    D(v)=v^4-11v^3+36v^2-(5/12)xv+(1/7)y.

Assume that both C and D weakly interlace B. Zero roots, repeated roots,
and common factors are allowed throughout this proof.

For either monic quartic interlacer A, cancel common factors in A/B.
At a B root of multiplicity m, interlacing forces multiplicity at least
m-1 in A, so no pole of order greater than one survives. The alternating
root order makes each surviving residue positive; canceled nodes can be
assigned zero weight. The monic degree difference gives residue sum one.
Thus C/B and D/B are transforms of two positive probability measures,
each supported on a subset of the nonnegative beta nodes:

    C(v)/B(v)=sum_i w_i/(v-a_i),
    D(v)/B(v)=sum_i h_i/(v-a_i),
    w_i,h_i>=0, sum w_i=sum h_i=1.

The measures need not have the same weights. No mixed C/D Gram matrix is
used. Each ordinary moment matrix is positive semidefinite, as is each
shifted moment matrix because its nodes are nonnegative. This direct
cancellation argument includes the weak repeated-root boundary; no strict
generic approximation is required as a hidden hypothesis.

The source-to-target map retains these separate measures and their common
support bound, then extracts scalar inequalities. It loses the full beta
root locations, so the final scalar packet is used only as a necessary
condition, not as an equivalence with interlacing. The cheapest decisive
test is the old exact surrogate in Section 5.

## 2. Full division moments and the two determinant identities

Write mu_j for the C moments and nu_j for the D moments. Division at
infinity gives the following complete initial sequences:

| j | mu_j | nu_j |
|---|---|---|
| 0 | 1 | 1 |
| 1 | 1 | 2 |
| 2 | 3 | 7 |
| 3 | x/3-16 | 7x/12-19 |
| 4 | 16x/3-373-4y/7 | 115x/12-632-6y/7 |
| 5 | 54x-59y/7+z-3969 | 199x/2-92y/7+z-7171 |

The source reconstructs both rows by formal inversion of
1-13u+55u^2-xu^3+yu^4-zu^5 and checks the entire product through order
five against each numerator. These are polynomial identities in x,y,z.

Set delta=x-75. The C shifted two-by-two determinant is delta/3, so
delta>=0. The next C shifted Gram determinant is

    0 <= det(mu_(i+j+1))_(i,j=0)^2
       = -x delta^2/27 +(15/7)delta y
         -(16/49)y^2 +delta z/3.                       (2)

The ordinary D three-by-three determinant is

    0 <= det(nu_(i+j))_(i,j=0)^2
       = (-333+2334delta-49delta^2)/144 -(18/7)y.        (3)

Both determinants are computed in their own positive measure. If delta=0,
(2) forces y=0, whereas (3) equals -37/16. Consequently delta>0.

## 3. The common support endpoint creates the strict floor

Put M=max a_i. Cauchy on the other four roots yields

    5M^2-26M-67<=0.

Since M>=13/5 and the left side is increasing thereafter, its value 9/20
at M=71/10 proves M<71/10. The coefficient identity

    M y-5z=sum_i (M-a_i) product_(j!=i) a_j>=0            (4)

is valid also when some roots vanish. In particular

    z <= (71/50)y.                                     (5)

Use this weak inequality before establishing y>0; there is no circular
strict-support assumption at the zero boundary.

Divide (2) by delta>0 and use (5). Dropping the nonpositive square term
gives

    x delta/27 <= (15/7)y+z/3
               <= (2747/1050)y.

Its left side is positive, so y>0. Since x=75+delta>75, this implies the
strict slope bound

    y > (8750/8241)delta.                              (6)

Meanwhile (3) gives

    2334delta >= 333+49delta^2+(2592/7)y.

Substitute delta<(8241/8750)y from (6). The coefficient left after collecting
y is positive, and therefore

    y > 333/[2334*(8241/8750)-2592/7]
      =161875/888583.

Its difference from 9/50 is exactly 96503/44429150>0. This proves (1),
including strictness at all zero/common/repeated-root boundaries. The bound
comes from positive moments and an elementary support inequality; it is
independent of a search for near-minimizing shapes.

## 4. The coupled residue bounds start the same tail at 2500

Keep the complete original rows and variable t:

    G_B=(1,13,55,x,y,z),
    G_C=(1,12,45,(2/3)x,(3/7)y),
    G_D=(1,11,36,(5/12)x,(1/7)y),
    beta=t^-1 G_B, C_raw=t^-1 G_C, D_raw=t^-1 G_D,
    O=sum_j binom(14,2j+1)t^j, E=sum_j binom(14,2j)t^j,
    P=O star beta,
    Q=(O^2+t^-1 E^2) star(beta^2+2t C_raw D_raw).

Star is coefficientwise multiplication. In particular the exponent -1
coefficient of Q is 28; the relative factor 2t and both Laurent shifts are
unchanged. At an original positive phase s, with P(-s)=0,

    P(-s)=182-20020s+2002x s^2-3432y s^3+2002z s^4,
    z=(12/7)y u-xu^2+10u^3-u^4/11,  u=1/s.             (7)

The exact response after (7) is the same polynomial R=Q(-s)/s^7 as in the
proved finite-phase theorem. Keeping its negative constant and every
positive monomial, and dropping only nonpositive monomials, gives

    R <= -(26075790/7)y^2+A(u)xy+B(u)x+C_1(u)y+D_1(u),
    A(u)=(153780300/7)u+(647843760/7)u^2,
    B(u)=1122025905u^3+3467704710u^4+6175260u^6,
    C_1(u)=(5521932000/11)u^4,
    D_1(u)=3654364350u^6+565082u^8.                     (8)

For completeness x=(sum a_i^3-52)/3<1223/10<123 follows from the same
M<71/10 and the square-sum anchor. Set h=50/9. By the newly proved y>9/50,
division of (8) by y^2 bounds R/y^2 above by

    F(u)=-26075790/7+123h A(u)+123h^2 B(u)
                     +h C_1(u)+h^2 D_1(u).            (9)

Every nonconstant coefficient is positive. Exact rational evaluation gives

    F(1/4100)
      =-5961609014706655321683472522194343
        /99603957308055354000000000000 < -59000 <0.

As s increases, u decreases, so (9) already proves a uniform tail from
4100. Retain this independent-cap calculation as a control. The proved
slope (6) gives a stronger bound without any new root or interlacer input:

    x/y=75/y+delta/y < 75*(50/9)+8241/8750
                       =10962223/26250 =r_1,
    x/y^2 < r_1*(50/9)=10962223/4725 =r_2.              (10)

Use these coupled ratios directly in the same positive monomials (8):

    F_c(u)=-26075790/7+r_1 A(u)+r_2 B(u)
                         +(50/9)C_1(u)+(2500/81)D_1(u).

Every nonconstant coefficient remains positive, and exact evaluation gives

    F_c(1/2500)
      =-653528391359305169367452997041401
        /13323669433593750000000000000 < -49000 <0.

Consequently **Q(-s)<0 at every original positive phase s>=2500**.
The gain retains the dependence of x on y supplied by the residue slope;
it does not optimize over a new model. Multiple first roots are included:
no derivative or division by a first-row discriminant was used. The
original-root condition and both interlacer hypotheses remain attached.

Combining this tail with the independently proved smallest-phase interval
theorem leaves only original phases in (1/80,2500) unresolved in the model,
at most three per shape counted with multiplicity. This is a reduction of
the remaining interval, not a proof of its full negativity.

## 5. Exact hostiles and nonvacuity

The previously frozen surrogate

    (x,y,z)=(104,50,37435088/3898125), s=15/2

satisfies the original first equation and has
Q(-15/2)=78541969368658673/18480>0. Its first polynomial has four simple
positive roots, and it passes the four previously displayed Newton
inequalities. It is not a member of the model: its beta polynomial has
only one real root.

It also passes the lower unbounded-support Hankel checks for both measures.
The ordinary three-by-three and shifted two-by-two determinants are,
respectively,

    C: 2693/63 and 29/3,
    D: 3338/63 and 103/3.

Their leading principal minors are positive as well. Thus these truncated
Hankel matrices are positive definite. But any residue support in [0,M]
must satisfy mu_4<=M mu_3 and nu_4<=M nu_3, by integrating
v^3(M-v)>=0. The anchor-derived upper bound M<71/10 already rejects this
surrogate at order four:

    mu_4-(71/10)mu_3=2159/105>0,
    nu_4-(71/10)nu_3=1091/42>0.

The C shifted three-by-three determinant also rejects it, with value
-70052935886/81860625. The restored coordinate is the actual support
endpoint, or stronger moment data retaining its effect; plain lower
Hankel positivity is insufficient. These necessary conditions are not
promoted to a complete characterization of beta interlacing.

The strict rational model with roots (1,3,9,22,30)/5 supplies nonvacuity.
The source computes the literal positive residues of both quotients and
recovers every moment through order five independently of series division.
The C-only repeated boundary (0,0,3,5,5) has measure weights 2/3 at zero
and 1/3 at three; it saturates (2) while failing (3). Direct evaluation
also gives D(5)=-25/4. This control verifies the weak-root cancellation
scope and the role of the second interlacer.

Finally the strict rational model has Q(-4100)>0 away from P(-s)=0.
The source checks that P(-4100) is nonzero, preserving the original-root
condition as a genuine requirement of the tail theorem.

## 6. Reproduction and audited freeze

The [standalone source](../../04-computation/continuing3_20260906_residue_floor.py)
uses only Python's standard library. It reconstructs both full division
moment rows, the determinant identities, the six-variable support identity,
the complete carried Laurent response, and every positive monomial of the
unchanged tail envelope. Its controls are the strict rational model, the
weak C-only boundary, and the exact surrogate. These identity and inequality
checks supplement the analytic proof; they are not a census of shapes.

The [JSON certificate](../../04-computation/continuing3_20260906_residue_floor_certificate.json)
stores the moments, determinants, full original-root response, both envelopes,
and the two coupled ratio bounds.
The optional path is relative to the working directory, so the source remains
portable when filed under `04-computation`.

```text
python -B 04-computation/continuing3_20260906_residue_floor.py --certificate 04-computation/continuing3_20260906_residue_floor_certificate.json
python -B -O 04-computation/continuing3_20260906_residue_floor.py
```

Both runs pass **94 always-active exact gates**, with byte-identical LF
[output](continuing3_20260906_residue_floor.out).
No producer module, numerical optimizer, or external algebra library is
imported. The output preserves its pre-audit candidate wording. The completed
independent audit and the filed promotion above establish current status;
the manifest records the unchanged source/output and prior report identity.

```text
source      8974668b0c36e129910142b03a66f95626b18861a009d037325f73f22f00f205
output      fa8b9936221b9f6cbea053fdc91c84f0eafe4faa25f19ffa8af5afe2d63033f0
certificate c2a20f53670198790529fa980b7d33913b9e5a1f64c3b206f9fc34aa42e9104e
```

The initial 91-gate candidate used the independent-cap cutoff 4100. The
94-gate revision retains that exact control and adds the two coupled bounds
in (10), giving 2500. This immediate improvement was proposed and independently
checked by the referee; the floor proof and all model hypotheses are unchanged.
No broader threshold optimization was performed.

The outside producer changed no shared repository state. Root subsequently
filed the audited package and its current navigation without a new theorem ID.
