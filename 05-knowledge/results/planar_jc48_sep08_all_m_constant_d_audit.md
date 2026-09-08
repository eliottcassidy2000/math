# Independent audit of the all-m constant-boundary gate

**Status: INDEPENDENT FULL ANALYTIC / SOURCE / REPLAY AUDIT PASS.**
This sidecar audits the
[primary proof](planar_jc48_sep08_all_m_constant_d.md). The accepted
analytic mechanism is uniform for every integer m>=2. It gives a
necessary condition for rational mates of the specified global
square-prefix family, not a sufficient condition or a polynomial-mate
classification on these surfaces. Parent owns promotion and Git.

## 1. Exact scope, maps and inherited data

The actual surface, source chart and boundary chart are

    W_m=(P1_x x P1_z) minus {z=x^m},
    t=1/(z-x^m), x=1/r, t=-r^m-r^(2m)b,
    D={r=0}, omega=dx wedge dt=r^(2m-2) dr wedge db.

Let `H in L2`, `L in L1`, `F=H^2+L`, and let the leading
t-coefficient N of H have exact degree `4m-2`. The theorem states
that any rational constant-Jacobian mate forces both H and L to
have constant restrictions to D. In this class that is equivalent
to constancy of F on D. Every coefficient and every finite root
of N is allowed, including multiple roots and square N.

The complete global filtration is the proved
[all-m DG theorem](planar_jc48_sep08_dg_genus.md). Its numerator
boxes are the actual global section spaces, not a selected collection
of native carriers. The other mechanism is the faithful formal
inverse and coefficientwise primitive equation in
[leading-field exactness](planar_jc48_sep08_leading_exactness.md).
The proof here rederives the two needed coefficients for every m.
Its m=2 specialization recovers the T1,T5 obstruction in the separate
degree-six note; it does not need that note as a provisional dependency.

The proof map goes from an actual rational Jacobian mate to an exact
differential on one chosen point of its coefficient field. It retains
all coefficient-field denominators and the primitive equation. It loses
finite-place residues and rational algebraization. The positive and
negative controls below exhibit that loss. In particular, constancy
on D does not mean that L itself is constant on the source.

I directly checked the volume exponent `2m-2`. The derivative
`dx/dr=-r^-2` times `dt/db=-r^(2m)` gives the displayed positive
factor. The proof then works in the original x-coordinate, so it
does not silently change the exactness equation to an unweighted
equation on the compact chart.

## 2. Full numerator boxes and all degree estimates

Write the complete numerators as

    H=(A+Bz+Qz^2)/(z-x^m)^2,
    L=(S+Rz)/(z-x^m),
    deg A,deg B,deg Q<=2m, deg S,deg R<=m.

Thus in the original source chart

    N=A+x^m B+x^(2m)Q,
    P=B+2x^m Q, M=S+x^m R.

If the leading coefficient of N is n!=0, the normalization
`(H,L,F,G)->(H/n,L/n^2,F/n^2,n^2G)` makes N monic while
preserving the exact mate equation, globality, degree and constancy
on D. There is no parameter-dependent rational source translation.

The top two coefficients of Q must vanish. Indeed their contributions
to N would have degrees at least `4m-1`, exceeding both `3m` and
`2m` for m>=2. Consequently `deg Q<=2m-2`. The m=2 endpoint
needs a separate observation: `3m=4m-2`, so the remaining top Q
coefficient can combine with B's top coefficient to give monicity.
It need not itself be one. This does not permit either of the two
already excluded higher coefficients of Q. For m>=3 monicity
sets the remaining top Q coefficient to one. The source retains
both alternatives explicitly.

Let `b0=B_(2m)`, `r0=R_m`, and choose `kappa^2=N`. Completing
the original quadratic and linear expressions gives

    D0=Q-P^2/(4N)=(4AQ-B^2)/(4N),
    E0=R-MP/(2N)
      =(2AR+x^m BR-BS-2x^m QS)/(2N),
    C0=M/kappa.

These are literal identities; I re-expanded both numerators
independently. On the chosen positive-leading infinity branch,
`kappa=x^(2m-1)(1+O(x^-1))`. The full boxes give

    D0=-b0^2 x^2/4+O(x), E0=O(x^2),
    C0=r0 x+O(1).

In the first identity, `deg AQ<=4m-2`, so it cannot cancel the
degree-4m term `-b0^2 x^(4m)` in the numerator. The numerator
of E0 has degree at most 4m, and M has degree at most 2m.
No finite-root or squarefreeness condition is used in these estimates.

When b0=0, both terms of the numerator of D0 have degree at most
`4m-2`, and the highest possible E0 term has degree `4m-1`.
The other E0 terms have degree at most `4m-2`; in particular
`deg AR<=3m<=4m-2`, including equality for m=2. Hence

    D0=O(1), E0=O(x), C0=r0 x+O(1).

When r0 also vanishes, E0 becomes bounded and C0 becomes bounded.
These are exact uniform degree implications of the original boxes,
not formal generic-coefficient assumptions.

The restrictions on D have signs determined by the actual denominator:

    H|D=A_(2m)+b0 b,
    L|D=-S_m-r0 b.

For H the square denominator gives a positive leading factor; for
L the denominator contributes the minus sign. There is no quadratic
b-term for H because Q_(2m)=0. Thus the quadratic coefficient of
F|D is `b0^2`; after it vanishes, F's linear coefficient is `-r0`.
Over C this proves the claimed equivalence between the three boundary
constancy conditions. A sign error in a preliminary transfer message
has no surviving occurrence in the inspected primary or source.

The restriction m>=2 is essential to this degree proof: for m=1,
the B and Q degrees have additional leading ties, and the first chosen
positive inverse index would not be available. No m=1 extension is
asserted here.

## 3. Faithful inverse field and exact coefficient extraction

The coefficient field is `K=C(x)(kappa)`. If N is a square, use
the field C(x) with one chosen square root, not a disconnected
double-cover presentation. Otherwise K is a quadratic field. In
either case the characteristic-zero derivation extends uniquely.

Formal inversion at v=0 gives

    F(x,T)=v^-4,
    T=kappa^-1 v^-1+T0+sum_(j>=1) T_j v^j.

The highest nonzero t-degree of any polynomial has uniquely least
v-valuation after substitution. Thus the substitution is injective
on the original rational field and pays every denominator of a
proposed rational mate. The chain rule, with v held fixed, is

    partial_x G(x,T)=(1/4)v^5 T_v.

Therefore

    ([v^(j+4)]G(x,T))'=(j/4)T_j

for the occurring nonzero indices, in particular every positive j
used here. This makes `T_j dx` exact in K. T0 is not constrained
by this equation and is not used.

For the exact coefficient formula, center with
`y=t+P/(2N)`, and put

    Phi(z)=1+2D0 z^2+C0 z^3+(D0^2+E0)z^4.

Then `q=1/(kappa y)` satisfies `v=q Phi(q)^(-1/4)`.
Changing variables in the formal residue for `[v^j]q^-1`
gives

    [z^(j+1)]Phi^(j/4)
     -(1/4)[z^j]Phi' Phi^(j/4-1).

Differentiating the formal power makes the second coefficient
`(j+1)/j` times the first one. Their difference is therefore
`-1/j` times it. Hence, for every j>=1,

    T_j=-(1/(j kappa))[z^(j+1)]Phi(z)^(j/4).

The centering contributes only T0. This derivation is formal and
uniform in j; it does not assume convergence or infer a general
formula from finitely many inverse coefficients.

## 4. First moving index: only pure D can have the top degree

Take `j=2m-3` and put `k=m-1`, so j+1=2k.
Before any forced vanishing, each D0 contributes at most two
x-degrees using two z-degrees, each C0 at most one using three,
and each E0 at most two using four. The D0 squared term remains
pure D0. A monomial `D0^a C0^b E0^c` contributing at the required
z-degree has

    2a+3b+4c=2m-2,
    x-degree <=2a+b+2c=(2m-2)-2b-2c.

Thus equality is possible only for b=c=0. No other term, or
lower coefficient inside D0, can contribute to that top degree.
Setting C0=E0=0 gives the exact pure-D inverse

    kappa v y=sqrt(1-D0 v^2).

The x^-1 coefficient of T_j is consequently

    binom(1/2,m-1)(b0^2/4)^(m-1).

I independently obtained the same coefficient from Lagrange
extraction as

    -(1/(2m-3)) binom((2m-3)/2,m-1)
                       (-b0^2/4)^(m-1).

These expressions agree exactly. Their generalized binomial
coefficient is nonzero for every m>=2. The highest numerator
degree is precisely bounded by `2m-2`, one less than the leading
degree of kappa. Therefore lower terms of `1/kappa` cannot mix
with a higher numerator term and cancel the x^-1 coefficient.
Exactness at infinity forces b0=0.

## 5. Second moving index: only pure C can have the top degree

After b0 vanishes, D0 is bounded, C0 has degree at most one,
and E0 has degree at most one. Take `j=6m-7`; its required
z-degree is `6m-6`. A contributing monomial with b C0-factors,
c E0-factors, and a D0-factors satisfies

    2a+3b+4c=6m-6,
    x-degree<=b+c.

If c is positive, `3(b+c)<=6m-6-c`, so the integer x-degree
is at most `2m-3`. If c=0 but a is positive, equality at
`2m-2` is again impossible. The sole top-degree monomial is
therefore `C0^(2m-2)`. Its coefficient in T_j, and hence its
x^-1 coefficient, are respectively

    -binom((6m-7)/4,2m-2) C0^(2m-2)/((6m-7)kappa),
    -binom((6m-7)/4,2m-2) r0^(2m-2)/(6m-7).

The upper binomial argument has odd numerator and denominator four;
none of the falling factors is zero. Again no numerator degree
above `2m-2` is possible, so corrections from kappa cannot cancel
the residue coefficient. Exactness forces r0=0. This completes
the actual boundary-constancy theorem.

Both residues are taken at an actual point of the field K.
For nonsquare N, its even degree gives two unramified infinity
points with local parameter r=1/x and kappa leading signs plus
and minus. Choose the plus branch. For square N, choose the
polynomial square root with that leading sign in K=C(x). In
either case the residue of an x^-1 term times dx is its negative.
The proof does not take an unweighted trace: the two nonsquare
residues can be opposite, but exactness in K requires each one
individually to vanish. Repeated finite roots do not change this
normalization argument at infinity.

There is also a precise stopping boundary for this operation.
After both slopes vanish, D0,E0,C0 are all bounded. Every
positive-index inverse coefficient is then
`O(x^-(2m-1))`, and T_{-1}=1/kappa has the same estimate.
Their infinity residues automatically vanish because m>=2.
This says nothing about T0 or finite-place exactness and is not
a sufficiency claim. It explains why extending this same
infinity-residue table cannot pay a further obstruction by itself.

## 6. Actual positive and negative controls

Put `y=x^(m-1)(1+x^m t)`. The proved basis gives y in L1,
so H=y squared is in L2. Direct compact substitution gives
`y=-r b`, and thus H and `L=lambda y` have zero restrictions
to D. H's leading N is exactly `x^(4m-2)`.

Independently in the original source coordinates,

    J(y,1/((2m-2)x^(2m-2)))=1.

Applying the chain rule to `F=y^4+lambda y` gives the stated
rational mate with denominator `4y^3+lambda`. It is a nonzero
rational function for every lambda. In particular, lambda nonzero
gives a nonconstant L with an actual rational mate. The theorem
cannot be strengthened to global constancy of L or to rational
nonexistence after the boundary slopes vanish. No global or
polynomial regularity of the mate is asserted.

For the opposite control, `H=y^2-t^2`, `L=t` are both global
and have constant restrictions to D. Their leading coefficient
is `N=x^(4m-2)-1`. At its simple root x=1 the rational inverse
coefficient `T2=-1/(4N)` has residue
`-1/(4(4m-2))`, which is nonzero. Normalized trace of a
primitive in the coefficient field would make this rational
coefficient exact in C(x), or equivalently its ramified residues
would have to vanish. Hence no rational mate exists. This control
uses a finite place that the two infinity tests deliberately
discard. The positive is a square-field case; this negative has
a nonsquare N with simple roots.

A further literal control confirms why the indices must move.
On W3 put `y=x^2(1+x^3t)`, `h=x^3(1+x^3t)`,
`H=y^2+h`, and `L=0`. Both y and h are global L1 sections.
Its complete numerator is

    A=0, B=-x^6, Q=x^4+x^3, N=x^10,
    D0=-x^2/4, C0=E0=0, kappa=x^5, H|D=-b.

The exact pure-D square-root inverse therefore gives

    T1=1/(8x^3), T5=x/1024, T3=-1/(128x).

The old T1 and T5 have zero infinity residue on this actual
nonconstant-boundary carrier, while the correct moving index
`j1=2m-3=3` detects it. This is an additional independently
checked audit control supplied during the final read; it makes
no minimality claim and requires no producer modification.

## 7. Independent computational path and frozen manifest

I have read the standalone producer. Its universe consists of the
universal numerator identities; inverse recursion versus Lagrange
extraction through its named indices; all admissible top-weight
monomials for m=2..9; all symbolic global coefficient boxes for
m=2..5; and the literal positive and negative source families for
m=2..9. Its special m=2 monicity relation is retained explicitly.
Finite controls support the uniform proof rather than replacing
it with an all-m census.

A separate temporary calculation passed **1,189 alternative controls**.
It groups the power as

    ((1+D z^2)^2+C z^3+E z^4)^a

and first chooses the C and E factors. For b such factors of C,
c of E, and a remaining D-exponent k, the exact coefficient is

    falling(a,b+c)/(b! c!) * binom(2a-2b-2c,k),
    2k+3b+4c=j+1.

Here `falling(a,n)=a(a-1)...(a-n+1)`.

This independently reconstructs all **570 admissible grouped
monomials** at both indices for m=2..12, verifies their exact
nonzero coefficients and the unique top-degree terms, and checks
the equality of the two first-coefficient formulas. It differs
from the producer's expansion into the four individual nonconstant
terms of Phi. Further direct source checks for m=2..5 verify the
rational positive Jacobian, its exact compact chart, and the
negative finite residue. The formulas above make the independent
path explicit; its temporary script is not a new proof dependency.

The complete frozen standalone source uses always-active exceptions
and imports no inherited mathematical program. I inspected the final
primary, including its corrected restriction sign, chosen infinity
branch, and explicit exclusion of T0 from the stopping statement.
Fresh independent normal and optimized runs both exited zero and
passed **439 gates**. Both outputs are byte-identical to the
**329-byte** frozen output. The verified pins are:

| Artifact | Bytes | SHA256 |
| --- | ---: | --- |
| Primary before status promotion | 11302 | `1a3506e3f25521122d5f1b61ae2714b5df07156ae27771405bec23788a4e3b33` |
| Source | 7034 | `7b2956263b1c03298c348e180b782cb3761f9e93e03cfcd86bbd4b818fa3b00d` |
| Frozen, independent normal, independent optimized outputs | 329 each | `28ddd67d3d7aac3ddb7abd508aea1be1bb857b3f7f5e0953d1d2e406570fb8f9` |

Semantic record SHA256:
`1552692247545b143e22e9391d5ad8d0ec21f67938462a72c36542bba652bd0f`.

Reproduction commands, from the worktree root:

```sh
python3 -B 04-computation/planar_jc48_sep08_all_m_constant_d.py
python3 -B -O 04-computation/planar_jc48_sep08_all_m_constant_d.py
```

The fresh output files were temporarily retained at
`/tmp/planar_jc48_sep08_all_m_constant_d_audit_normal.out` and
the matching `_optimized.out`. The separate grouped-weight
calculation is additional corroboration, not a new proof dependency.

**Final acceptance:** exact all-m scope, complete numerator boxes,
monic normalization, m=2 degree tie, both restriction signs, faithful
formal field substitution, all-index Lagrange formula, unique top
weights and nonzero coefficients, both radical-field cases, stopping
boundary, sharp controls, source and both replay modes PASS. No
correction remains. Producer files have not been edited by this
auditor; parent owns promotion and checkpointing.
