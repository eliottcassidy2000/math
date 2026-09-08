# All graph degrees: rationally integrable source-linear submersions and their unit-response orders

**Status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.**
This classifies a specified layer on the explicit surfaces W_m. It is not
an entry theorem for arbitrary Keller maps or a new case of JC(2).

## 1. Inheritance, theorem and the retained source

Work over C. For every integer m>=1 use the actual two-chart surface

    W_m=(P1_x x P1_Z) minus {Z=x^m},
    U0=A2_(x,t), Uinf=A2_(r,b),
    x=1/r, t=-r^m-r^(2m)b, b=-x^m-x^(2m)t,
    D=W_m minus U0={r=0},
    omega=dx wedge dt=r^(2m-2) dr wedge db.

The charts and complete section filtration are inherited from
[planar_jc48_sep08_dg_genus.md](planar_jc48_sep08_dg_genus.md), section 2.
The closest m=2 mechanism is the complete layer, source-degree descent and
fourfold-root hostile in
[continuing10_20260907_dg_linear_carrier.md](continuing10_20260907_dg_linear_carrier.md).
The fixed-source response and its intrinsic connection are the proved
[component-jet theorem](planar_jc48_sep06_torsion.md), sections 2--4,
with antecedents **THM-3770**,
`01-canon/theorems/THM-3770-vertical-principal-part-equalizer-and-log-canonical-dressing-gate.md`,
and **THM-3412**,
`01-canon/theorems/THM-3412-hamiltonian-principal-part-differential-and-prufer-torsion-arms.md`.
No external-priority claim is made for these elementary rational-integration
or relative-exactness mechanisms.

Consider nonconstant F in O(W_m) of degree at most one in the original
source variable t. A rational mate means G in C(x,t) with
`J(F,G)=F_xG_t-F_tG_x=1`. A global submersion means dF is nonzero at every
point of W_m, including D. These are distinct predicates.

**Classification.** Write `F=A(x)t+B(x)`. Apart from the constant-A
exception below, every global submersion with a rational mate has

    A=a(x-h)^n, a!=0,
    B=B0+a sum_(j=m+1)^n binom(n,j)(-h)^(n-j) x^(j-m).       (1)

The complete allowed parameters are:

| Graph degree | Allowed A and parameters |
|---|---|
| m=1 | A=a!=0, B=B0; or (1) with n=2 and arbitrary h |
| m>=2 | (1) with h!=0 and either m+1<=n<=2m-2 or n=2m |

An interval whose lower endpoint exceeds its upper endpoint is empty.
In particular m=2 has only n=4. The case n=2m-1 always has one critical
point on D and is excluded, even when its entire affine source is smooth.
All displayed parameters supply actual functions, not just necessary
root patterns. The rational mates are

    G0=-x/a                              if A=a,
    G0=1/[a(n-1)(x-h)^(n-1)]             in (1).              (2)

For every admitted pure-power case, put c=B(h), g=F-c, and use the fixed
affine response module

    C_F=C[x,t]/D_F(C[x,t]), D_F=F_x partial_t-F_t partial_x.

The unit theta=[1] has **exact c-primary order n-1**. Its j-th canonical
connection derivative has exact order n-1+j for every j>=0. The constant-A
m=1 case has theta=0 in this affine source module.

The five live concepts are the two-chart section space, the degree of a
rational primitive, the normal derivative at D, the actual labelled special
fibre, and the ring in which a response is required to be regular. The
canonical hostile is the missing n=2m-1 case; the corrected near miss is
replacing affine submersion by global submersion. The least-used sidecar is
the boundary component omitted by restriction to the original source chart.

## 2. Complete global layer and rational-integrability lemma

For an arbitrary global source-linear function, write
`A=sum_(j=0)^(2m) a_j x^j`, allowing zero high coefficients. The complete
two-chart intersection is

    deg A<=2m,
    B=B0+sum_(j=m+1)^(2m) a_j x^(j-m),                     (3)
    F_inf=B0-sum_(j=0)^m a_j r^(m-j)
                -b sum_(j=0)^(2m) a_j r^(2m-j).            (4)

Indeed the b-coefficient is `-r^(2m)A(1/r)`, so globality forces the
degree bound. The remaining negative powers cancel precisely as in (3).
Conversely (4) is polynomial. A=0 therefore forces F constant and is
outside the theorem. Equivalently B is the positive-degree part of
`A/x^m` plus a free constant; its constant term is not counted twice.

For every A!=0, `C(x,t)=C(F)(x)` via `t=(F-B)/A`, and

    D_F=-A partial_x |_F.                                (5)

Its constants are exactly C(F). A rational mate exists if and only if
`-dx/A(x)` has a rational primitive over C(F). In this polynomial-A
setting, this is equivalent to A being constant or a single pure power
`a(x-h)^n`, n>=2.

Here is a complete proof of the last assertion, including coefficient-field
extensions. If a nonconstant A has a simple finite root, 1/A has a nonzero
residue there and cannot be a rational derivative. Otherwise let its degree
be n and let it have r distinct roots, of multiplicities m_i>=2. A putative
rational primitive R over C(F) has a pole of order exactly m_i-1 at each
root and has no other finite poles: differentiation creates a pole at any
finite pole of R, and no such extra pole occurs in 1/A. At infinity,
`R=R_infinity+gamma x^(1-n)+...`, gamma!=0, so R has no infinity pole
and R-R_infinity has a zero of order n-1 there. Its total pole degree is
`sum_i(m_i-1)=n-r`. Equality of the degrees of zeros and poles of a rational
function gives `n-1<=n-r`, hence r<=1. Since A is nonconstant, r=1.
This argument is unchanged over the larger constant field C(F); no new
finite poles or extra independent component constants can appear.
The primitives (2) prove sufficiency. No bounded degree search for a mate
or restriction on the number of roots was used.

The inherited m=2 hostile `F=(x^2-1)^2t+x^2` is globally submersive:
at x=1,-1 its x-derivative is respectively 2,-2, and its boundary
b-derivative is -1. But `1/(x^2-1)^2` has residues -1/4 and 1/4 there.
Thus smoothness alone does not put a function in the classified list.

## 3. Exhaustion of all source and boundary critical points

For a pure power with n>=2, the only possible affine critical locus is
x=h. It is critical exactly when B'(h)=0, since `F_t=a(x-h)^n` and
`F_x|x=h=B'(h)`. If n<=m then B is constant and this derivative vanishes.
If k=n-m>=1, formula (1) gives

    beta=B'(h)
        =a(-1)^(k-1) binom(n-2,m-1) h^(k-1).             (6)

For completeness, the binomial identity is the coefficient of T^(k-1)
in `(1-T)^n(1-T)^(-2)=(1-T)^(n-2)`: the coefficient on the product side
is `sum_(i=0)^(k-1)(k-i)(-1)^i binom(n,i)`. Thus (6) holds in all degrees.
It implies that affine smoothness requires n>m and, unless n=m+1,
also h!=0. For n=m+1 the derivative is the nonzero constant a.

The whole boundary test is already in (4):

    F|D=B0-a_m-a_(2m)b,
    F_r|D=-a_(m-1)-a_(2m-1)b,
    F_b|D=-a_(2m).                                      (7)

For pure powers this gives three exhaustive cases.

* If n<=2m-2, then F_b=0 and
  `F_r|D=-a binom(n,m-1)(-h)^(n-m+1)`. In the source-smooth range n>m,
  it is nonzero exactly when h!=0. In particular, the h=0 source-smooth
  n=m+1 branch still fails on the entire boundary.
* If n=2m-1, F_b=0 while F_r is affine in b with nonzero slope -a.
  It vanishes at exactly one point of D for every h. Thus no parameter
  in this degree is globally submersive.
* If n=2m, F_b=-a is nonzero throughout D. Combining this with (6) permits
  h!=0 for m>=2 and every h for m=1.

For constant A=a, one has F=a t+B0 and
`F_inf=B0-ar^m-ar^(2m)b`. Its boundary normal derivative is -a when m=1
and zero when m>=2; its tangential derivative is always zero on D.
This proves the exceptional constant-A line in the table and completes
the classification.

The missing exponent 2m-1 is an actual geometric obstruction, not a
consequence of an omitted coefficient box. It is the precise transition
where the b term first appears in the boundary normal derivative while
the boundary restriction remains constant.

## 4. Exact unit order with every affine component retained

For any admitted pure-power F, put z=x-h and c=B(h). Then

    g=F-c=zK,
    K=a z^(n-1)t+C(z),
    C(z)=(B(z+h)-B(h))/z, C(0)=beta!=0.                 (8)

The fibre g=0 has exactly the two reduced components E1={z=0} and
E2={K=0}. They are disjoint because C(0)!=0. K is irreducible as a
primitive linear polynomial in t with coprime coefficient and constant
term; its coordinate ring is C[z,z^-1]. Every other affine fibre is
irreducible because the coefficients of `a z^n t+B(z+h)-d` are coprime
when d!=c. Thus no further special component is being discarded.

The primitive q=G0 from (2) has its only affine pole on E1, of order n-1;
it is regular on E2 and every other affine divisor. The scalar leading
principal coefficient in the actual uniformizer g is

    alpha=beta^(n-1)/[a(n-1)] != 0.                    (9)

For the upper order, the actual polynomial

    P_(n-1)=g^(n-1)q=K^(n-1)/[a(n-1)]

satisfies `D_F P_(n-1)=g^(n-1)`.
For the lower order, `g^(n-2)q` has a genuine simple pole at E1 and is
regular on E2. Every rational primitive of g^(n-2) differs from it by
H(F), by (5). Cancelling the E1 pole forces H to have a simple pole at c;
the same pole appears at E2, where g has order one, and cannot cancel.
A higher-order H pole already fails at E1. This proves
`g^(n-2)theta!=0` with no degree bound on a possible polynomial primitive.
It also covers n=2, where this is the nonvanishing of theta itself.

The gradient is a unit ideal by the proved affine smoothness, and the
constant-field hypothesis of the component-jet theorem is (5). Its
canonical map therefore identifies the unit's torsion with a scalar
principal part supported on [e_E1] modulo the common diagonal, with
top term `alpha g^(-(n-1))`. The entire torsion submodule has one arm at c.
The intrinsic connection differentiates all its scalar coefficients;
characteristic zero and alpha!=0 imply exact order n-1+j after j
derivatives. No formula for a bounded packet is being extrapolated.

**Order spectrum in the classified layer.** For fixed m>=2, the positive
unit orders are exactly

    {m,m+1,...,2m-3} union {2m-1}.                     (10)

Over all m>=1 they are exactly `{1} union {3,4,5,...}`. Order two is
absent. For every e>=3, choose m=e, n=e+1 and h!=0 to realize e; m=1,
n=2 realizes one. This is a spectrum only for the classified source-linear
global submersions with rational mates. The incoming order-two control
`v=x+x^3t` is m=2,n=3,h=0 and has its missing boundary critical point.
Its higher-source-degree thickenings are outside this classification.

## 5. Generic fibres, polynomial mates and the omitted boundary component

For a pure power, a generic fibre on U0 has x!=h and is isomorphic to G_m.
If n<2m, (7) makes F constant on D, so a generic fibre of W_m has no
boundary point and remains G_m. If n=2m, F|D is an affine bijection and
adds exactly the x=infinity point; the generic W_m fibre is P1 minus
{h}, hence A1. For constant A in the m=1 exception, the source generic
fibre is already A1 and no generic boundary point is added. These are
geometric generic-fibre descriptions; they do not change the response ring.
For n<2m the fibre at the boundary value contains D as an additional
component. The assertion in section 4 that c is the only reducible fibre
is explicitly about the original affine source; it is not an assertion
about all fibres of W_m.

There is no global regular mate for any nonconstant source-linear global
F on any W_m, even without the submersion or rational-integrability
conditions. Indeed the inherited all-degree source descent says

    J(F,P) in C[x]  iff  P=H(F)+Q(x),  H,Q polynomial.

To recall its mechanism, a highest t-coefficient P_N in a proposed
polynomial P satisfies `N A' P_N-A P_N'=0`; hence P_N is a scalar times
A^N. Subtracting that multiple of F^N lowers the t-degree and terminates.
Thus a nonzero constant Jacobian requires A constant. By globality B is
then constant, and every polynomial mate is `-x/a+H(F)`. The first term
has a nonzero pole on D while H(F) is global, so no such mate is global.

The constant m=1 exception gives an explicit demonstration that changing
the ring restores a lost component. Here omega=dr wedge db has no zero
at D, so D_F is a global derivation of O(W_1). For F=a t+B0 define the
global response quotient `C_F^W=O(W_1)/D_F(O(W_1))`. Then

    theta_source=0,                    theta_W!=0,
    (F-B0) theta_W=0.                                  (11)

The source primitive is q=-x/a. The global polynomial witness is
`(F-B0)q=-xt`, whose second-chart expression is `1+rb`. It satisfies
`D_F(-xt)=F-B0`. To prove theta_W nonzero, retain both components of
the special fibre t=0: D={r=0} and the closure T of {t=0} in U0.
They are reduced and disjoint, since in Uinf their factors are r and
1+rb. The only pole of q on W_1 is a simple pole on D. A correction
H(F) cancelling it introduces a pole on T, where q is regular. The paid
rational constant field excludes any different correction. This proves
exact global order one directly, without transferring the affine-plane
torsion theorem to a different ring.

For verification, in the second chart

    D_F=a r^2 partial_r-a(1+2rb) partial_b,

which is polynomial and compatible with the source derivation. Restriction
to U0 therefore induces a map of response quotients that kills the
nonzero class theta_W. For m>=2 omega vanishes on D and the corresponding
source Hamiltonian derivation need not preserve O(W_m); no global response
module or global connection is asserted there.

The source-target transfer in this classification is exact: F maps to
the rational primitive of -dx/A, then to its labelled fibre poles. It
preserves rational integration and the fixed affine derivation. It loses
global critical points and boundary poles, restored respectively by (7)
and by the actual second chart. Order two and the m=1 constant example
are cheap hostile tests for those two losses.

## 6. Exact controls and remaining scope

The standalone source reconstructs the entire linear coefficient basis in
both charts for m=1..10, all pure-power degrees n=2..2m, and all h in
{0,1,-2}: 300 exact cases, 93 admitted, plus ten constant-A controls.
It checks the binomial derivative and both boundary derivatives without
floating arithmetic. All 39 multiplicity words from {2,3,4} on one,
two or three fixed distinct roots check the rational-primitive obstruction;
a simple-root residue is retained separately. The special-fibre factor
identity checks the arbitrary-power upper witness without expanding large
powers. The m=1 ring-change witness is checked in both charts.

These finite checks support the mechanisms; the complete all-m and
all-degree statements are the proofs above. Run

    python -B 04-computation/continuing10_20260908_dg_linear_all_m.py
    python -B -O 04-computation/continuing10_20260908_dg_linear_all_m.py

Before installation the source runs standalone from any working directory.
It imports no other producer, reads no mutable repository data, and writes
no files. Both modes pass 1,844 always-active exact gates with identical
raw LF stdout. Independent audit is complete; see the linked referee below.

The next question is outside this closed layer: whether higher-source-degree
global submersions with rational mates can realize the missing affine unit
order two, or whether a new invariant explains its absence beyond the
linear layer. The present theorem proves no such higher-degree exclusion.

The [independent complete audit](continuing10_20260908_linear_classification_audit.md) accepts the classification, unit-order spectrum, and source/global response contrast. Filing changes only status and routing; frozen source and output bytes are unchanged.
