# Clock duality, rational functions, and the cost of repeated blocks

2026-10-15. **PROVED:** the elementary clock-form dualities, local inverse,
repeat-cost identity, and a doubling-only mean-current obstruction below.
**CITED:** geometric Langlands/S-duality, Gijswijt's block/glue hierarchy,
and the genus-zero moonshine-prime classification. **OPEN:** any transfer
of those theories that proves the actual Collatz coarse-mean bound.

The anchor and new positive result are
[the source-clock mesh bounds](collatz_clock_mesh_20261015a.md).
The niche result is [all-depth unit lifting](unit_lift_clock_20261015a.md).
This note records which suggested connections give exact structures and
which supply only a search direction.

## 1. A duality that carries the source and its scale

For a nonconstant finite current block, write its positive definite form as

    G=[[L,U],[U,V]], D=V-U^2/L>0,
    tau=(U+i sqrt(det G))/L=U/L+i sqrt(D/L).

Then Q_G(j)=Lj^2+2Uj+V=L|j+tau|^2. The two matrices

    R_d=[[1,d],[0,1]], S=[[0,-1],[1,0]]

act by congruence G -> A^T G A and send tau to tau+d and -1/tau,
respectively. Under S, G'=[[V,-U],[-U,L]] and

    Q_G(j)=j^2 Q_G'(-1/j)       (j != 0).                     (1)

At j=0 use the homogeneous source vector (j,1): its transformed vector is
S^-1(j,1)=(1,-j), including the point at infinity. There is no singularity
of the original quadratic cost. Constant blocks have D=0 and lie on the
boundary of this upper-half-plane parametrization.

This is a literal elementary SL_2 action. However, after S the leading
coefficient is V, an energy, not the original physical length L. The
source vector, normalization factor and integer length must be retained.
Discarding j^2 in (1), or treating V as elapsed time, manufactures a false
energy decrease. The projective coordinate is useful for shapes of forms;
it is not a bound on their size.

The connection to the user's S-duality suggestion is deliberately precise
and limited. [Kapustin--Witten](https://arxiv.org/abs/hep-th/0604151)
relate geometric Langlands to a twisted N=4 gauge theory compactified on a
Riemann surface, using electric-magnetic duality and additional geometric
structures. Equation (1) is an elementary model for exchanging coordinates
while transporting observables. No equivalence between that gauge theory
and Collatz has been constructed here.

The Q two-generation matrix [[3,4],[2,3]] preserves X^2-2Y^2, while on
its expanding eigenvector (sqrt(2),1) it multiplies both coordinates by
3+2sqrt(2). Thus its projective state is fixed while its Euclidean energy
grows by (3+2sqrt(2))^2 per iteration. This positive control explains
why the exact Pell conjugacy in
[THM-4630's inherited model](golden_nested_trees_20261013c.md#6-incoming-q-covariance-work-a-genuine-two-step-pell-mechanism)
does not bound the Collatz current's magnitude. Golden normalization has
the same scale issue: preserving value and recording the rewrite history
are different tasks; see [the golden-memory note](conway_golden_memory_20261013d.md).

## 2. The cumulative kernel has a local inverse

Let P be the lower unit-triangular cumulative-sum matrix on a finite
ordered event list, and let D=diag(L_i) with all integer L_i>0. The
clock-energy kernel is K=P^T D P. If B=P^-1 is the first-difference
matrix, then exactly

    K^-1=B D^-1 B^T.                                         (2)

The right side is tridiagonal. Thus the nonlocal cumulative energy has a
local precision operator, provided the gap lengths have not been erased.
For any fixed source charge vector v the resolvent

    F_v(z)=v^T(K^-1+zI)^-1 v
          =sum_lambda ||projection_lambda v||^2/(lambda+z)    (3)

is a rational function on the Riemann sphere. Its poles are negative
eigenvalues with nonzero source weight. At z=0 it gives v^T K v.
An unmarked spectrum loses the weights; a boundary resolvent usually
loses an arbitrary bulk source v. The marked source is essential.

For unit gaps and N events, the smallest eigenvalue is
`2-2 cos(pi/(2N+1))`, of order N^-2. Hence localization of the inverse
does not remove the low-frequency/mean mode. The new mesh theorem pays
fine fluctuations and leaves precisely such large-scale source weight.
Restricting the variational identity
`v^T K v=max_f(2v^T f-f^T K^-1 f)` supplies lower, not upper, bounds.

This is different from the existing
[rational approximants to the Collatz transfer law](collatz_rational_approximants_riemann_sphere_20261013.md)
and [HYP-9260](../hypotheses/HYP-9260-collatz-transfer-law-on-the-riemann-sphere-rational-approximants-exact-circle-poles-unfolding-moebius-words-q-covariance-map.md).
Those already handle the rational function 1/(1-x-y), irrational winding,
and the Q Möbius map. Here the rational object is the finite source-clock
Green function, with its positive spectral weights retained.

## 3. Gijswijt suggests an authenticated block/glue grammar

The [original Gijswijt paper, sections2--3](https://cs.uwaterloo.ca/journals/JIS/VOL10/Sloane/sloane55.pdf)
proves a hierarchy in which level-m blocks repeat m+1 times and acquire
glue; concatenated glue supplies the next level. It proves unboundedness
and finiteness of the individual blocks/glues. The first4 occurs at220.
Our literal curling-number implementation independently reproduces it.
These facts suggest retaining interfaces when compressing long repeated
words; they do not identify Collatz labels with Gijswijt's sequence.

For any finite current-increment block with summary (L,d,U,V), repeating
it r times gives exactly

    L_r=rL,  d_r=rd,
    U_r=rU+Ld*r(r-1)/2,
    V_r=rV+dU*r(r-1)+Ld^2*r(r-1)(2r-1)/6.                  (4)

Proof: copy j starts at current jd, so its cost is
V+2jdU+L j^2 d^2; sum j=0,...,r-1. The same derivation gives U_r.
If d!=0, the leading energy is Ld^2 r^3/3. A late glue block can
cancel the terminal drift but cannot erase cost already incurred.
For actual Collatz use, repetition also requires verified phase, label,
and source-port compatibility; an algebraically legal word is not
automatically a native trajectory or adjacent source block.

The useful proposal is a grammar of **signed terminal-weight rows** with
authenticated source ports. Its primitives may have large or small net
drift, but every concatenation records its full quadratic cost. The
[terminal-mesh certificate](collatz_terminal_mesh_20261015a.md) implements
the first step without needing a universal grammar theorem.

## 4. A hostile showing why the odd Collatz relation is indispensable

Impose only |t|<=1, t(1)=1 and t(2n)=-t(n). On [X,2X), let M=sum J.
Summation by parts gives

    2M=sum_(1<=n<2X) a_X(n)t(n),
    a_X(n)=-w(n)*#{m in [X,2X):m>=n}+W(n)[X<=n<2X],
    w(n)=3[n=2 mod3]-1.                                    (5)

Group each n=2^e q by its odd part q. Write
`A_X(q)=sum_e (-1)^e a_X(2^e q)`. The exact largest possible
absolute value of 2M is

    |A_X(1)|+sum_(odd q>1)|A_X(q)|.                          (6)

This follows by choosing each free odd label's sign; the labels then
extend consistently to all integers by the doubling rule.
For X divisible by4, the X/4 odd q in [X,3X/2) have no further
doubled representative below2X and satisfy |a_X(q)|>=X/2-1.
Thus the admissible maximum |M| is at least X^2/16-X/8, so
the normalized mean energy M^2/X^3 is unbounded in this class.
For X=16 an explicit extremizer has 2|M|=410 and mean energy
10.260009765625, already far above1/2.

This is not a counterexample to Collatz. The first failed implication is
that even-branch cancellation plus bounded labels would force a small
source mean. The strongest survivor is the exact regrouping (5)--(6).
The missing constraint is the actual odd-edge relation
t(n)=-t((3n+1)/2), with coalescing terminals identified. This is what
the next certificate uses rather than assuming a Q covariance closure.

## 5. Where the requested prime families actually connect

| Family | Exact structure worth retaining | Loss or stopping boundary |
|---|---|---|
| 15 moonshine/supersingular primes | Genus0 of X_0(p)^+ permits a rational parameter; the primes are 2,3,5,7,11,13,17,19,23,29,31,41,47,59,71 | This is the Fricke quotient, not generally X_0(p); the covering and branch data are extra information |
| Non-generous primes | The least primitive root has a delayed order lift modp^2; the valuation h controls all later levels | The residue order alone loses h; choosing a different lift can change it |
| Wagstaff primes | ord_(W_q)(2)=2q for q>3; an exact odd-order orbit count follows | The order-q symmetry is the full Paley symmetry at11 but not at43 |
| Cullen primes n2^n+1 | p-1=n2^n constrains the split multiplicative group | The order of a selected unit still needs proof |
| Woodall primes n2^n-1 | p+1=n2^n constrains the nonsplit norm-one torus | Inertness in the particular quadratic field is a separate test |
| Stern primes | Absence of a representation p'+2b^2, b>0, is an additive witness problem | A finite list or a multiplicative order does not settle the absence for all primes |
| Leyland primes x^y+y^x, x,y>1 | Primality forces coprime x,y of opposite parity | These necessary conditions neither certify primality nor imply a Collatz estimate |

For the first row, inherit the repository's already stronger
[Ogg audit](ogg_triangular_chowla_20260927_audit.md), with primary source
[Ogg, Automorphismes de courbes modulaires](https://www.numdam.org/item/SDPP_1974-1975__16_1_A4_0/).
In particular X_0(11) has genus1 even though X_0(11)^+ has genus0.
Definitions of the additive/exponential families are recorded in
[Stern A042978](https://oeis.org/A042978),
[Cullen A002064](https://oeis.org/A002064),
[Woodall A003261](https://oeis.org/A003261), and
[Leyland-prime A094133](https://oeis.org/A094133).
No new completeness or infinitude claim is made about these families.
The non-generous/Wagstaff rows are pursued to proofs and exact checks in
the unit-lift companion; the other rows remain typed connections.

## Audit and research direction

[The exact control script](../../04-computation/experiments/clock_duality_glue_20261015a.py)
checks 4,368 repeat identities, 5,460 homogeneous duality identities,
five precision matrices, the doubling support function for k=4..16,
and Gijswijt's first4. [Output](clock_duality_glue_20261015a.json).

The useful common move is to transport the predicate together with the
coordinate that its quotient forgets: source weight for a resolvent,
scale for a projective action, lift valuation for a residue unit, and
glue cost for a repeated word. The Collatz target is now a small bound
on authenticated coarse terminal-weight rows. Calling two coordinate
systems dual does not itself make that bound true.
