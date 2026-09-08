# Independent referee: three whole global quadratic pencils

**Status: PASS — independent analytic audit and FINITE-EXACT controls.**
The primary packet `continuing12_20260908_quadratic_pencils.md` is accepted
without mathematical repair. The first two specified whole pencils contain
no globally submersive function. The seventh-order pencil has exactly the
displayed global submersion locus, exact original-source unit order three,
two full torsion arms, and the stated two complex rank-one branches. Its
full Weyl annihilator is the displayed relation-generated left ideal.

This referee concerns the fixed W2 surface and the three named entire
two-dimensional discriminant pencils. It does not classify the other
pencils, change the response ring to O(W2), or assert a regular global
mate or a general Jacobian-conjecture result.

## 1. Target, dependencies, and independent path

Audited primary SHA256:

    report 02cf3a3322ad1fcd2b51cb21e1f2b9fcb79275810a4e0ac04908d2bdb9b53815
    source 755ca42c481de6c6467d46ba881e2a2be309249d61e5e09d4fecb5d58d76c25c
    output 058e75ed9c7f4de9836ae0c50d2bc7c990b121defc94ea598a052b1aadd4f0c7
    cert   e083d60b782384a91d03538fa25065664f1bbd491aa1cc9ef12767fb8dbce0af

The actual supplier `planar_jc48_sep06_torsion.md`, Sections2–3, was read.
It requires a unit polynomial gradient and rational constants C(F), and
identifies source response torsion with all labelled component principal
parts modulo the common diagonal. Its canonical connection acts as
ordinary differentiation of those scalar principal parts. Both
hypotheses and the entire component lists are paid below. The two-chart
quadratic coefficient criterion was checked against
`planar_jc48_sep08_dg_quadratic.md`. The six possible independent exact
pencils are inherited from `planar_jc48_sep08_exact_pencils.md`; this audit
does not reclassify all exact radical differentials.

The companion imports no mathematical producer. It recovers the whole
global coefficient kernel from negative Laurent coefficients, derives
the three matching problems, integrates the field derivation directly,
extracts local scalar coefficients by triangular Taylor jets, and tests
the Weyl compiler through a separate finite Laurent-action implementation.
The finite controls support the analytic arguments; finite parameter
samples are not substituted for the universal assertions.

## 2. Complete global boxes and the two excluded pencils

Write F=N(x)t²+P(x)t+Q(x), N nonzero. Under x=1/r and
t=-r²-r⁴b, the b² coefficient gives deg N<=8. The b coefficient then
gives deg P<=6, and the constant coefficient gives deg Q<=4. Thus no
higher-degree terms have been omitted before solving the finite box.
The complete negative Laurent kernel is

    P6=2N8, P5=2N7,
    Q4=N8, Q3=N7, Q2=P4-N6, Q1=P3-N5.

It has dimension15. The referee recovers these six independent conditions
from literal substitution, rather than assuming the primary's matching
formula. A target translation remains free; below s denotes Q at u=0,
where u=x-h. Changing the original x-constant to this value is a target
translation, not an additional coefficient constraint.

For the third pencil, N=u³(au+d), so deg P<=4. Since u³ divides both
N and P²-4NQ, it divides P²; hence P=u²(pu²+qu+v). Globality determines
Q from this P up to s. The discriminant's u⁸ coefficient is p², and
after p=0 its u⁶ coefficient is q². Both must vanish in the specified
pencil. Consequently

    F=s+u³(au+d)t²+v u²t.

The determinant of the two pencil coefficient vectors is dv², so the
entire pencil condition is exactly dv!=0. Both actual source derivatives
vanish along the entire line u=0. This excludes every member under
consideration before global submersion is tested.

For the fourth pencil, N=u⁴(au²+d), and the same argument gives
P=u²(pu²+qu+v). Globality gives

    Q=(p-a)u²+[q+(4a-2p)h]u+s.

The discriminant's u⁸ coefficient is (p-2a)², so p=2a. The only
remaining off-pencil coefficient is 2q(v-2d)u⁵. Independence is
a(v-2d)²-dq²!=0. This splits exhaustively and disjointly into:

* q=0, k=v-2d, with ak!=0;
* v=2d, with dq!=0.

The intersection of the two branches has zero determinant and is not a
whole-pencil member. With w=1+u²t, their respective functions are

    s-d-k+(au²+d)w²+kw,
    s-d+(au²+d)w²+quw.

The first has an entire source critical line u=0. The second is source
submersive: on u!=0 the two (u,w) derivatives are
w(2auw+q) and 2(au²+d)w+qu. If w=0, the second is qu!=0. If a!=0
and w=-q/(2au), the second is -dq/(au)!=0. If a=0, the first vanishes
only at w=0. On the omitted source line the actual derivative is q!=0.

Nevertheless literal second-chart substitution gives

    F_b|D=0,
    F_r|D=-(4ah+q)(3h²+b).

Thus r=0,b=-3h² is critical for every member of the second branch. If
4ah+q=0, the entire added divisor is critical. This is a complete
two-chart obstruction, and includes arbitrary h and a=0.

## 3. Seventh-order matching and all-point submersion

For N=u⁷(au+d), u⁷|P² and deg P<=6 force P=u⁴(c6 u²+c5 u+c4).
The two global top equations uniquely give c6=2a and c5=2d-4ah.
Put c4=k-4dh. The four remaining global equations determine Q up to s.
Their unique result is

    z=u-2h+u³t,
    F=s+u z[(au+d)z+k].

Conversely direct substitution pays every global coefficient condition,
and the full discriminant identity is

    P²-4N(Q-s-g)=u⁷[(k²+4ag)u+4dg].

The coefficient determinant is dk². Thus the whole pencil has rank two
exactly when dk!=0. No condition a!=0 is introduced; the degree-drop
case a=0 is part of the classification.

On u!=0 set y=uz. This is an actual localized coordinate with
J(u,y)=u⁴, and

    F=s+(a+d/u)y²+ky.

Its u derivative at fixed y is -dy²/u²; vanishing forces y=0, where
the y derivative is k!=0. On the omitted line u=0, literal source
differentiation gives F_t=0 and F_u=-2h(k-2dh). Therefore source
submersion in the whole-pencil locus is equivalent to h(k-2dh)!=0.
At either excluded factor, the entire u=0 line is critical.

For the full second chart put

    R=h²(3-hr)+b(1-hr)³.

Here z=-rR and

    F=s+(1-hr)R[(a(1-hr)+dr)R-k].

At r=0 let R0=3h²+b. Then F_b=2aR0-k. For a=0 it is -k!=0
everywhere. For a!=0, at its unique zero R0=k/(2a), the actual normal
derivative is dR0²!=0. Thus every added-chart point is paid, and the
global submersion locus is exactly

    d*k*h*(k-2dh)!=0, with a and s arbitrary.

All affine gradient claims are literal polynomial gradient claims; over
C, absence of a common zero gives the unit gradient ideal required by
the torsion supplier.

## 4. Primitive, constants, component completeness and order

Write g=F-s and L=(au+d)z+k. On the rational field,

    u=d y²/(g-ay²-ky),
    D_F y=-d u²y²=-d³y⁶/(g-ay²-ky)².

The nonzero Jacobian shows that g,y are algebraically independent, and
the displayed inverse recovers C(x,t)=C(g)(y). The derivation is a
nonzero multiple of partial_y over C(g). In characteristic zero its
rational constants are exactly C(g), as required.

The referee integrates - (g-ay²-ky)²/(d³y⁶) with respect to y. This gives

    G=g²/(5d³y⁵)-gk/(2d³y⁴)+(k²-2ag)/(3d³y³)
      +ak/(d³y²)+a²/(d³y).

Direct substitution then verifies D_F G=1 and that g³G is polynomial in
the original source. Equivalently,

    G=V/(30d³u³z³),
    V=(16a²u²-8adu+6d²)z²+k(7au-3d)z+k²,
    g³G=L³V/(30d³).

The special source fibre consists of exactly u=0, z=0, and L=0.
The two linear polynomials in t are primitive: the z coefficient is
u³ and its constant at u=0 is -2h!=0; the L coefficient is
u³(au+d), while its constant at u=0 is k-2dh!=0 and at u=-d/a,
when this root exists, is k!=0. They are therefore irreducible by
Gauss's lemma. Their pairwise separation values are -2h,k-2dh,k,
so they are disjoint. The fibre is reduced, also following from source
submersion.

For g=c!=0 the quadratic discriminant has odd u-valuation seven, hence
is not a square in C(u). The polynomial in t is primitive: the only
possible common root of N and P is u=0, and its constant term on that
fibre is -c!=0. If a!=0, P(-d/a)=k d⁴/a⁴!=0 pays the other root
of N. Thus every other source fibre is irreducible. Only the support of
gcd(N,P) is fixed here; its exponent can jump. The referee retains the
admitted controls (a,d,h,k)=(0,1,1,4) and (1,2,1,8), with gcd u⁵ and
u⁶ respectively. Neither invalidates the primitivity argument.

The denominator of G is supported only on u=0 and z=0; G is regular
at L=0. At z=0,

    residue_z(z*g²G)=k⁴/(30d³u)!=0.

Thus g²G has a genuine pole on z=0 and is regular on L=0. A common
rational target correction H(g) cannot remove that non-diagonal pole
tuple. Since the rational constants have been proved to be C(g), no
different rational primitive repairs it. Combined with the polynomial
g³ witness, this proves the distinguished source unit annihilator is
exactly (g³). The actual component theorem now gives exactly two full
torsion arms, both at target value s. No component or fibre of W2 is
silently substituted into this original-source statement.

## 5. Full jets and the coefficient-rank wall

Let H=g³G. In an actual local parameter e, expand
g=f1 e+f2 e²+O(e³) and H=H0+H1 e+H2 e²+O(e³). The independent
extraction is triangular:

    C3=H0, C2=H1/f1, C1=(H2-C2*f2)/f1².

At u=0 this uses the actual source jet z=-2h+u+u³t before expansion.
At z=0 it uses (u,z), where u is a unit. The resulting component
coefficient vectors, in the gauge setting the regular L component to
zero, are exactly

    C3=((k-2dh)³(k²+6dhk+24d²h²)/(30d³), k⁵/(30d³)),
    C2=((k-2dh)(4ad²h²+2adhk+ak²-6d³h)/(3d³), ak³/(3d³)),
    C1=(a²k/d³,a²k/d³).

All extracted coefficients are scalar on the entire component; t cancels.
The referee also verifies the complete remainders modulo the local cube.
The omitted-jet hostile is decisive: at a=d=h=k=1, incorrectly fixing
z=-2h during the u expansion changes C2_u from -1/3 to -7/3.

For a=0, C2_u=-2h(k-2dh)!=0 and C2_z=0 while C3_z!=0. Hence the
coefficient rank is two. If a!=0, C1 is a nonzero diagonal vector, so
rank one is equivalent to both of the following differences vanishing:

    C2_u-C2_z=-2h(4ah²-6dh+3k)/3,
    C3_u-C3_z=-8h³(12d²h²-15dhk+5k²)/15.

With q=k/(dh), these give 5q²-15q+12=0 and
a=3(2dh-k)/(4h²). The quadratic has two distinct complex roots
(15±sqrt(-15))/10, neither 0 nor 2. Therefore both branches satisfy
the full smoothness assumptions and a!=0. The two independent equations
define the asserted codimension-two locus in the smooth parameter space
(with s free). For real nonzero d,h, the second equation has negative
discriminant -15d²h², so there is no real alignment locus.

The independent quotient-ring controls verify both complex branches
simultaneously, including invertibility of every required denominator.
They also retain the different complex locus k²+6dhk+24d²h²=0: here
the Eu cubic coefficient can vanish while the Ez cubic coefficient
remains nonzero. At the tested a=d=h=1 slice Eu has order two, Ez has
order three, the source is smooth, and this is not the rank-one wall.
The primary never assumes a cubic pole on every component.

## 6. Exact Weyl annihilator for arbitrary operators

The supplier's canonical connection is ordinary differentiation on the
entire principal-part module. Thus theta=C1/g+C2/g²+C3/g³ in
C² tensor J, J=g^-1 C[g^-1]. The inherited Euler interpolation argument
gives D theta=span(C1,C2,C3) tensor J. It generates both complete arms
off the alignment locus and one on it. The highest coefficient C3 is
always nonzero in the component quotient, so every canonical derivative
has exact scalar order 3+j, including on the rank-one locus.

Let D=C<g,nabla>/([nabla,g]-1). PBW with derivatives on the left gives
a unique remainder

    P=p0(nabla)+p1(nabla)g+p2(nabla)g² modulo Dg³.

Each arm is a free C[nabla]-module on 1/g, since
g^-j=(-1)^(j-1)nabla^(j-1)(1/g)/(j-1)!. Hence this remainder acts as

    [C1*p0+C2*(p1-nabla*p0)
      +C3*(p2-nabla*p1+nabla²*p0/2)](1/g).

The triangular coefficient substitution is invertible over C[nabla].
The kernel of the constant matrix with columns C1,C2,C3 over C[nabla]
is its constant relation space tensored with C[nabla]. For a constant
relation lambda1*C1+lambda2*C2+lambda3*C3=0, the inverse substitution
therefore produces

    R_lambda=lambda1+(lambda2+lambda1*nabla)g
      +(lambda3+lambda2*nabla+lambda1*nabla²/2)g².

Every remainder annihilating theta is a sum of C[nabla] multiples of
these R_lambda. Conversely every R_lambda and g³ annihilates theta.
This proves exactly, for arbitrary finite Weyl operators,

    Ann_D(theta)=Dg³ + sum_lambda D R_lambda,

where lambda ranges over any basis of the constant relation space.
There is one such relation off the wall and two on it. Each nonzero
relation produces a nonzero PBW remainder of g-degree at most two, so
Dg³ alone is strictly too small throughout this family. No minimality
claim about the chosen left-generator lists is needed.

The independent implementation tests actual finite Laurent actions for
all21 PBW basis terms g-degree0..2 and derivative degree0..6, checks the
triangular inverse coefficient by coefficient, and verifies the universal
compiler action. It also applies both wall generators and their first
three left derivatives in the exact quadratic parameter quotient, plus
five independent smooth real controls and derivative-order tests. These
finite checks are separate from the all-operator proof above.

## 7. Frozen reproduction and stopping boundary

The companion passes **316 always-active exact gates**. Normal and -O
executions produce identical raw-LF stdout and certificate bytes.

    source 1bc2bba491c7621b619e99ecd7465846c2965f42257dfb27446dd93ea4b71436
    output 77fa7db6a780f1495fdff5fc5cb105af3bbfe9e6fa19c97183f0e31bb0ac2b47
    cert   53f02db7f7956118c565b482882906ae062d80ee17fd32361e8e8efccbd3a1e5

    python -B 04-computation/continuing12_20260908_quadratic_pencils_audit.py
    python -B -O 04-computation/continuing12_20260908_quadratic_pencils_audit.py

Outside repository installation, the source writes its certificate beside
itself; under 04-computation it writes to 05-knowledge/results. No
mathematical producer is imported in either location.

The complete family classification and exact module presentation are
accepted. The response statements remain in C[x,t]. In particular the
symplectic form on the added chart is r² dr wedge db, so the original
Hamiltonian derivation is not automatically a derivation of O(W2).
Neither global regularity of a rational primitive nor unrestricted
quadratic or higher-degree exhaustion is asserted.
