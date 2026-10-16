# Cubic unit-path rings, mixed primes, and a marked Serre interface

2026-10-16, clock_duality_oct15. **PROVED, scoped:** the new norm-one unit,
full unit-path component formula in both cubic fields, the mixed-prime
example at 19, a short explicit connector, the full cubic metric, and a valuation lower bound for
a specified unit-word alphabet. **PROVED algebra + FINITE-EXACT:** the joint
circuit, induced graph controls and arithmetic computations. **CITED:**
classical algebraic-number-theory inputs and Serre modularity.
**OPEN:** four-chromatic cubic contexts, six-colour forcing, Collatz uniform
coverage deficit and planar JC. Self-audited; no independent referee or
historical priority claim.

## 1. Inheritance and the actual reroute

The closest mechanism is the complete unit-path localization theorem in
[the cyclotomic component note](cyclotomic_unit_component_20261011.md).
Its proof requires the path ring to contain the full ring of integers.
That premise cannot simply be copied to a noncyclotomic cubic field.
The hostile is the five-colour residue map on the smaller cubic base,
C1 of [cubic packets](cubic_crossedge_packets_20261016d.md). The corrected
near miss is reading a residue graph's lower chromatic bound as a lower
bound for its source field. The least-used sidecar is the **individual
prime ideal above a rational prime**, together with its conjugate.

The incoming [THM-4661, odd-degree carriers](../../01-canon/theorems/THM-4661-odd-degree-carriers-a-failed-six-colour-exterior-has-at-least-twenty-vertices-and-equality-testers-are-seven-colour-certificates.md)
was read in full, including its update acknowledging C1 and the earlier
64-term circuit. Its list-criticality and nineteen-vertex exclusion make
small six-colour searches less plausible. Its proofs of the high/low count
and the nineteen-vertex contradiction were checked directly here; this is
not an independent referee of its cited graph-theoretic and fish inputs.
In particular its threshold concerns failed exteriors, not all five-cores.

Anchor / niche / wildcard and the shared six-concept board are recorded in
[the interval-atlas companion](collatz_interval_atlas_20261016e.md).
The operative reroute is from blind word balls to a unit with a prescribed
pole-zero pair. It then supplies the missing integral coordinate needed
for a genuine universal path-coverage theorem.

Write

\[
 j=\sqrt{-11},\quad K=\mathbb Q(\zeta_5,j),\quad
 L=\mathbb Q(\zeta_{55}),\quad \xi=\zeta_{11},
\]
\[
 f(w)=w^3+2w^2-4w-4=0,\quad 1<w<2,\quad
 F=\mathbb Q(w),\quad K'=KF,\quad L'=LF.
\]

The selected w is twice the earlier theta. The degrees are
`[K':Q]=24`, `[L':Q]=120`, `[L':K']=5`, `[L':K]=15`.
All fields use their actual complex embeddings and conjugation c.

## 2. A direction that distinguishes the primes above nineteen

**B1 (PROVED).** The following is a unit direction in K':

\[
 u=\frac{w+3-j}{w+3+j}
   =\frac{11w-6+(w^2-w-12)j}{38},\qquad u\bar u=1.        \tag{1}
\]

The exact identities explaining its denominator and trace are

\[
 (w^2+6w+20)(4-w)=76,\qquad
 f(-3+j)=f(-3-j)=76,                                    \tag{2}
\]
\[
 u+\bar u=\frac{11w-6}{19},\qquad
 w=\frac{19(u+\bar u)+6}{11}.                           \tag{3}
\]

They follow by reducing polynomials modulo f and `j^2+11`. In particular
u leaves L: otherwise (3) would put w in L, contradicting the cubic degree.
The numeral 76 here is an exact polynomial evaluation and denominator
invoice; no inference is made from its other appearances in the repo.

At 19,

\[
 f(x)\equiv(x-4)(x^2+6x+1),\qquad j^2\equiv8\pmod{19}.
\]

The quadratic is irreducible over F_19 and splits over F_(19^2), with roots
`-3+j,-3-j`. The old fields K and L each have four primes above 19, of
residue degrees two and ten respectively. They are conjugation-fixed;
c acts on their residues by an odd power of Frobenius, hence sends j to -j.
The three roots of f in each residue field are simple. Over each old prime,
the branch w=4 is fixed by c and the other two branches are exchanged.
Thus in **each** of K' and L' the twelve primes above 19 consist of four
fixed primes and eight nonfixed primes in four conjugate pairs.

On the branch `w=-3-j`, `w+3+j` has valuation exactly one and its conjugate
numerator is a unit. Indeed (2) has 19-adic valuation one and the derivative
at that root is a unit; Hensel lifting preserves this exact valuation.
So u has valuation -1 there, +1 on the conjugate branch, and zero on the
fixed branch. The checker separately lifts all three roots modulo 19^2.
A rational-prime label cannot distinguish these three behaviours.

## 3. Universal unit-path coverage in the two cubic fields

For M=K' or L', let

\[
 U_M=\{z\in M:z\bar z=1\},\qquad
 \Gamma_M=\sum_{\text{finite}}\mathbb Z U_M.
\]

This is exactly the connected component of zero in the full unit-distance
graph. It is a ring, since products of unit directions are unit directions.

**B2 (PROVED).** With A the full ring of integers of M,

\[
 \boxed{\Gamma_M=\bigcap_{c(\mathfrak p)=\mathfrak p} A_{\mathfrak p}.}
                                                               \tag{4}
\]

Every component is an additive coset of this ring. Equivalently, the
allowed points are the diagonal elements of M in the finite adelic region
that is integral at every conjugation-fixed prime and unrestricted at the
other primes. This is a statement about actual finite unit walks, not only
about local approximations.

**Proof of the new premise, A subset Gamma_M.** The old connector
`eta=(-3+sqrt(-55))/8` has `eta+bar eta=-3/4`, so
`1+eta+bar eta=1/4` is in the path ring; hence 1/2 is too.
The old 82-step connector to 1/11 lies in Q(zeta5): if

\[
 \beta_r=\frac{\alpha^r-3}{\alpha^{-r}-3},\quad
 \alpha=\zeta_5,\quad 1\le r\le4,
\]

then `|beta_r|=1` and `11 sum_r beta_r-38=1/11`.
Equation (3) now puts w in Gamma_M. Also `(5+j)/6` is a unit direction,
so `j=6((5+j)/6)-5` belongs to the smaller path ring. Consequently
`Z[zeta5,(1+j)/2]` is contained in it; this is O_K because the discriminants
125 and -11 are coprime. For L, `O_L=Z[zeta55]` is already contained.

For completeness, `r=w^2/2+w-1` obeys

\[
 g(r)=r^3-r^2-3r+1=0,\quad w=r^2-3,\quad
 \operatorname{disc}(g)=148.
\]

The element w^2/2 is integral, satisfying `x^3-6x^2+8x-2=0`.
Dedekind's index criterion excludes an extra index at 2:
`g=(x+1)^3 mod 2` and `(g-(x+1)^3)/2=-2x^2-3x` is coprime to x+1 modulo 2.
The only other discriminant prime is 37, to exponent one. Thus
`O_F=Z[r]=Z[1,w,w^2/2]` and `disc(F)=148`.

The discriminants of K and L have support only {5,11}, disjoint from
{2,37}. Hence `O_(BF)=O_B O_F` for B=K,L. One direct proof of this familiar
coprime-discriminant fact uses the trace-dual integral bases: multiplying
any integral element of BF by disc(F) puts it in `O_B O_F`; exchanging B
and F does the same with disc(B). Bezout and coprimality put the original
element in that product order. All the required generators now lie in
Gamma_M, proving A subset Gamma_M.

**Proof of (4).** If c(p)=p then every unit direction has valuation zero,
since `2 v_p(z)=v_p(z bar z)=0`. Hence Gamma_M is contained in A_p.
If c(p) differs from p, choose `a in c(p) minus p`. Then `a/bar a` is a
unit direction with a pole at p. Thus Gamma_M is not contained in A_p.
The Dedekind overring lemma proved in section 2 of
[the old component note](cyclotomic_unit_component_20261011.md) identifies
an overring A subset R subset M with the intersection of exactly the A_p
that contain it. Apply it to R=Gamma_M. QED.

The full-integer-ring step is essential. Pole testing by itself does not
justify applying a Dedekind overring theorem to an arbitrary subring.
Here the new unit's **trace identity**, rather than a field-degree slogan,
pays that step.

There is no uniform short-path conclusion. Enumerating finite unit words
terminates for each point certified by (4), since arithmetic and norm tests
in these fixed number fields are exact. A bound for the length of the first
witness is a separate question.

### A 328-step connector, and an unreachable rational point

**B3 (PROVED).**

\[
 \frac{w-4}{19}=\frac{u+\bar u-2}{11}\in\Gamma_{K'},
 \qquad \frac1{19}\notin\Gamma_{L'}.                    \tag{5}
\]

Multiplying the 82 old unit steps for 1/11 separately by u, bar u, -1,
and -1 gives the first point in exactly 328 unit steps. Every multiplication
preserves unit length; repeated vertices are allowed. In contrast 1/19 has
a pole at the fixed branch w=4, so (4) excludes it even from the larger
field's zero component. The smaller field is excluded as well.

Thus Gamma_M is **not** obtained from O_M by inverting a set of rational
primes: u has a pole above 19, yet 1/19 is not in the ring. Individual prime
ideals, not rational-prime support alone, are the right denominator data.
The earlier five-colouring of K' at the fixed branch remains valid.

## 4. A valuation invoice for a restricted unit-word alphabet

Put A_0=Z[xi]. The prime 19 is inert in Q(xi), since its order modulo 11 is
ten. Consider an exact relation

\[
 \sum_{k=l}^{h} a_k u^k=0,\quad a_k\in A_0,\quad
 a_l a_h\ne0,\quad l<h.                                \tag{6}
\]

**B4 (PROVED).** Both a_l and a_h are divisible by 19 in A_0. If each
coefficient is recorded as an integer sum of signed xi powers, the top and
bottom layers each require at least nineteen steps; genuine coupling in
this alphabet therefore requires at least 38 steps.

At the pole branch of u, all coefficient valuations are nonnegative. If
a_h were a unit, its term would have the unique least valuation -h,
impossible in a zero sum. The conjugate branch gives the claim for a_l.
For the word count, write a coefficient as `sum_(j=0)^10 c_j xi^j`.
Divisibility by 19 means all c_j are equal modulo 19, since Phi_11 is
irreducible modulo 19. A nonzero coefficient has two unequal c_j, whose
difference has absolute value at least 19; hence `sum |c_j|>=19`.
This survives adding a multiple of Phi_11 to the presentation.

This is a lower bound for nonzero extreme layers of (6), not for every
unit-distance edge or every graph circuit. Old cyclotomic relations with
all a_k=0 and directions outside the specified alphabet are excluded from
its scope. It is not a replacement for the incoming exact rank checks.

## 5. The exact three-coordinate metric

Put `a=(20+11j)/19`. The unit u satisfies

\[
 u^3+a u^2-\bar a u-1=0.
\]

**B5 (PROVED).** For arbitrary A,B,C in L, the direction
`d=A+Bu+Cu^2` has norm one if and only if, with
`D_1=B bar A+C bar B` and `D_2=C bar A`,

\[
 D_1-aD_2+\overline{D_2}=0,
 \quad |A|^2+|B|^2+|C|^2+\bar aD_2+a\overline{D_2}=1.  \tag{7}
\]

Indeed `u^2=u^(-1)-a u+bar a` and
`u^(-2)=u-bar a u^(-1)+a`. Expand d bar d and compare its coefficients
in the L-independent basis `1,u,u^(-1)` to obtain (7).

A unit direction has either **one or three** nonzero coordinates in this
basis; it cannot have exactly two. If A or C vanishes the first equation
forces a second coordinate to vanish. If B=0 and AC differs from zero,
that equation would imply `|a|=1`, whereas `|a|^2=1731/361`.
Consequently the graph on any of the coordinate planes `L+Lu`, `L+Lu^2`,
or `Lu+Lu^2` is a Cartesian product of old-field unit graphs. This recovers
the low-width mechanism of G3 in
[elliptic cross-edges](elliptic_crossedges_20261016c.md) in the new basis;
(7) retains the full metric for leaving those planes.

This basis is explicit:

\[
 w=\frac{-14+11j+(39+11j)u+19u^2}{11}.
\]

The checker converts the previous generator v to this basis, verifies (7),
and finds all three coordinates nonzero. Thus the triple-support alternative
is nonempty. It also compares (7) with the independent 120-coordinate norm
on all 125 integer coefficient triples in [-2,2]^3.

An edge changing both a port address and a seed address is only a useful
first screen: another linear coordinate system can still reveal a Cartesian
product. The full metric (7), the retained ports and the actual boundary
relation remain the stronger tests. Three-coordinate support alone does
not force more colours.

## 6. Joint circuits and the cross-edge test

Keep the prior cubic generator

\[
 v=\frac{w(1-\xi)+2\xi^2}{w(1-\xi^{-1})+2\xi^{-2}}.
\]

Since `w=-3+j(1+u)/(1-u)`, define `A(u)=-3(1-u)+j(1+u)`. Then

\[
 v[A(u)(1-\xi^{-1})+2(1-u)\xi^{-2}]
 =A(u)(1-\xi)+2(1-u)\xi^2.                             \tag{8}
\]

This is an exact bilinear relation in u and v, unlike a collection of
independent first-power directions. Expand j by its quadratic Gauss sum,
multiply by xi squared, and subtract the median multiple of Phi_11 in each
of the four u/v coefficient blocks. The checked integer relation has
positive and negative mass 32 each. Thus it gives two new 32-step unit
walks with a shared endpoint. Replacing v by xi squared while leaving u
fixed fails the relation, an exact hostile to losing the common cubic
parameter. The receipt stores every step and the reduction coefficients.

The full induced graph controls, with the original eight ports literally
first in every coordinate list, are:

| Declared universe | Vertices | Unit edges | Chromatic number |
|---|---:|---:|---:|
| radius-two bank over 1,u,v,uv times xi powers, plus fifth roots | 4,708 | 9,250 | 3 |
| eight port-translates of the joint circuit packet | 1,884 | 4,052 | 3 |
| joint circuit, earlier circuit, and its u-rotation | 696 | 722 | 3 |
| eight port-translates of two orders of the 328-step connector | 5,074 | 10,807 | 3 |

The third assembly has no edges between points with disjoint copy
memberships. The connector assembly has six edges changing both port and
seed addresses with no alternative common address. **All six have
old-field displacement in L**; none changes the cubic coordinate. Merely
counting its 5,685 cross-copy edges with a cubic endpoint would be misleading:
most are ordinary port edges repeated over a new seed point. The checker
retains both addresses and tests the displacement itself.

Every graph has a verified three-colouring and contains the old pentagon.
All 570 old six-colour port patterns therefore persist, by
[THM-4614, the double-port extension theorem](../../01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md):
its nine-double-port census is unchanged in L'/K of odd degree fifteen.
This conclusion uses that census and the exterior colouring; it is not
inferred from a global colour count alone.

A separate cheap contact probe tested 224 distinct real Mobius transforms
`r=(a+bw)/(c+dw)`, coefficients in [-2,2], nonzero determinant, at rational
separations `1/4,1/2,3/4,1`. For
`u_r=(r+j)/(r-j)`, it tested squareness of `D(4-D)/11` and `/55` in F,
where `D=t^2+1-2t(r^2-11)/(r^2+11)`. All 224 hits were the trivial t=1
rhombus control; none of the other 1,568 tests hit. This is a bounded probe,
not an exclusion theorem for contacts. Its reproducible exact GP script
and receipt are included with the other artifacts.

## 7. The cubic has a genuine elliptic and Serre interface

Take the specific elliptic curve

\[
 E:\quad y^2=x^3+2x^2-4x-4.                              \tag{9}
\]

**PROVED algebra.** The chosen point `(w,0)` is a nonzero 2-torsion point,
and F is its x-coordinate field. Irreducibility modulo 3 and the nonsquare
discriminant 592 make the splitting group S3. The full 2-torsion field is
`F(sqrt(37))`, of degree six over Q; its quadratic resolvent is Q(sqrt37).
The other two roots are not elements of L': adjoining them would introduce
sqrt37, which cannot enter an odd-degree extension of L. Thus the three
real embeddings used in a trace identity are not three automorphisms of
our chosen cubic carrier.

Direct Weierstrass invariants are

\[
 c_4=256,\quad c_6=640,\quad \Delta=9472=256\cdot37,
 \qquad j(E)=65536/37.
\]

PARI independently computes field discriminant 148 and elliptic conductor
148. Exact direct point counts at the twelve good primes from 3 through 43
agree with PARI. In particular `#E(F_13)=14`, so `a_13=0`.

**CITED theorem applied to an explicit object.** The representation
`rho_(E,5):G_Q -> GL_2(F_5)` is continuous. The Weil pairing makes its
determinant cyclotomic, hence it is odd. Its Frobenius polynomial at 13 is
`T^2+13`, whose discriminant is 3 modulo 5, a nonsquare. It is therefore
irreducible over F_5. It is absolutely irreducible as well: a common
invariant line over an algebraic closure would be an eigenline of complex
conjugation, whose distinct eigenvalues +1,-1 already have F_5-rational
eigenlines. That would contradict irreducibility over F_5.

Serre modularity consequently applies to **this residual representation**.
The theorem concerns continuous, odd, absolutely irreducible two-dimensional
representations of G_Q, and is proved using compatible lifts and changes
of residual characteristic. This is a real arithmetic connection, not a
claim that our terminal-transfer matrix is a Galois representation.
[Khare--Wintenberger, *Serre's modularity conjecture I*, introduction and section 9](https://www.math.ucla.edu/~shekhar/papers/results.pdf).
For the torsion construction and Frobenius polynomial, see
[Sutherland, *Elliptic Curves*, Lecture 26, sections 26.3--26.4](https://math.mit.edu/classes/18.783/2019/LectureNotes26.pdf).

No minimal level or weight is identified for rho_(E,5) here; the conductor
148 above is that of E, not an automatically identical residual conductor.
Modularity also does not select `(w,0)`, its real embedding, or the fixed
prime over 19. Those marks are needed by the geometric construction.
Frobenius traces of the unmarked curve cannot decide a six-colour boundary
assignment or a signed Collatz source invoice. This is precisely the
missing-map hazard of MISTAKE-233, avoided by keeping the target predicates
separate.

## 8. Connection contract and the next useful questions

| Source -> target, map | Preserved | Lost / required sidecar | Decisive test |
|---|---|---|---|
| unit words -> prime-ideal integrality, diagonal embedding (4) | exact finite path existence | word length and colour relation | (5): cubic 19-denominator is reachable, rational 1/19 is not |
| cubic root -> marked 2-torsion of (9), w -> (w,0) | defining algebraic coordinate and Galois action | unmarked traces lose root, embedding and valuation branch | split versus nonsplit branches at 19; E[5] irreducibility at 13 |
| independent native controllers -> physical intervals, A1--A3 of companion | original atoms, full clocks and scalar current | no automatic small anchor residual or common Gram matrix | exact residue hostile and all-scale sampled bound still owed |
| geometric projection -> old carrier, normalized trace | old ports | unit length loses the variance identified in C3 of cubic packets | keep the full cubic displacement when counting cross-edges |

The promising next geometric target is a coupled edge whose displacement
uses the cubic coordinate and whose address changes both factors. The
current connector gives universal path coverage but still fails this
specific test. A second is to use the prime-ideal valuation vectors to
select compatible **joint** directions, rather than reject a rational prime
wholesale. A third is to seek a four-chromatic exterior before attempting
the much stronger six-colour equality tester. None is supplied by (4).

On the Collatz side, the corresponding useful move is to permit different
controllers while retaining the exact common source functional. Adelic
approximation alone cannot identify two integer terminals; the companion
constructs native hostile pairs at every finite residue precision. The
planar-JC analogue remains the existing requirement to retain the marked
differential in [the Lyness bridge](jacobian_elliptic_bridge_20261016c.md);
no new Keller map or Jacobian theorem is claimed.

## 9. Reproduction and checks

Programs:
[cubic direction and graphs](../../04-computation/experiments/cubic_adelic_20261016e.py),
[arithmetic audit](../../04-computation/experiments/cubic_adelic_arithmetic_20261016e.py).
Exact coefficient arithmetic follows the inherited 120-dimensional basis;
three finite-field screens only reject pairs, and every survivor receives
an exact norm test. All old port pairs are also checked without the screen.
Every SAT colour witness is replayed on every edge. No unverified UNSAT
claim is used. Normal and optimized executions agree byte for byte for all four graph
modes and the arithmetic audit. The contact probe's GP output is parsed as
JSON and its counts are checked explicitly; no GP database is needed.

Receipts: [direction bank](cubic_adelic_20261016e_bank.json),
[joint circuits](cubic_adelic_20261016e_circuit.json),
[two-circuit union](cubic_adelic_20261016e_two_circuits.json),
[connector assembly](cubic_adelic_20261016e_connector.json),
[arithmetic and metric](cubic_adelic_arithmetic_20261016e.json),
[contact probe](cubic_contact_probe_20261016e.json),
[hash manifest](interval_cubic_audit_20261016e.json).

```text
python3 -B 04-computation/experiments/cubic_adelic_20261016e.py --mode bank
python3 -B 04-computation/experiments/cubic_adelic_20261016e.py --mode circuit
python3 -B 04-computation/experiments/cubic_adelic_20261016e.py --mode two_circuits
python3 -B 04-computation/experiments/cubic_adelic_20261016e.py --mode connector
python3 -B 04-computation/experiments/cubic_adelic_arithmetic_20261016e.py
gp -fq 04-computation/experiments/cubic_contact_probe_20261016e.gp
```

Add `-O` after `-B` for the optimized Python checks. The direct point counts
and invariant identities are independent of the PARI calculations; PARI
supplies the exact conductor computation. B2's integral-order argument and
B5's metric equivalence are written out above, not inferred from a finite
coordinate census.

Classical background was checked against the repo's CORE-PAPERS and primary
sources. The adelic interpretation uses the usual restricted-product and
valuation conventions, not a new approximation theorem; see
[Milne, *Algebraic Number Theory*, sections 3 and 7](https://www.jmilne.org/math/CourseNotes/ANTc.pdf).
The proofs of the specific ring generation, mixed-prime behaviour, connector
and native transport are given here. No local-solubility or modularity claim
is used as an unproved global colouring or convergence principle.
