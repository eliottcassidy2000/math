# Every polynomial-coefficient y-linear carrier has generic genus at least six

**Status: PROVED + FINITE-EXACT; independent analytic/source audit PASS.**
The all-parameter genus formula below extends the incoming constant-slope
theorem to arbitrary polynomial coefficients, including common repeated
factors and vanishing of the y coefficient at p=0. The flow consequence
inherits the audited actual completion and fixed-input comparison. This is
not an exclusion of all universal carriers or of compositions of different
flows, and it is not a solution of JC(2).

## 1. Inheritance, domain, and the precise new statement

Work over any characteristic-zero field K. Put

```
D=p³-y²,  I=p²D(A(p)+B(p)y),  A,B in K[p],  B != 0.
```

The actual source coordinates are p=t(1+x²t), y=xtp, u=x²t;
K(p,y)=K(x,t), with x=yp/D and t=D/p². The bracket is the original
J_(x,t), so {p,y}=-D/p. The source carrier is K+p²D K[p,y].

The closest proved mechanism is the discriminant and invariant-curve method
in [the incoming genus theorem](continuing9_20260907_flow_carrier_genus.md)
and [its arbitrary-A, constant-B extension](continuing9_20260907_ylinear_genus.md),
with their independent audits. The canonical hostile to inferring scalar
nonrationality from non-local-nilpotence is the rational flow of u=x²t;
it is outside this carrier. The corrected near miss is A=0,B=p²: its genus
is eight, whereas the constant-B A=0 formula gives six. The underused
sidecar is the critical-*value* lemma in the incoming y-linear proof:
critical curves do not invalidate exclusion of a transcendental value.

The live concepts are: retained invariant rather than outer Hamiltonian;
fixed versus moving branch factors; Newton polygons at common coefficient
zeros; the index-two infinity branch as a primitivity certificate; and
the actual completed source flow. The source-to-target map is the invariant
I and the p-projection of its smooth proper generic curve. It preserves
the rational function field and exact scalar iteration. A displayed
discriminant alone loses the ramification at poles of its leading
coefficient; the local table below restores exactly that information.

**Theorem.** The generic curve of I is geometrically integral. Let
M=deg B, and put deg A=-infinity when A=0. Set

```
n=max(4M+13, 4 deg A+7).
```

Take all geometric roots q of B, and also q=0 whether or not B(0)=0.
For each such distinct q, let a=ord_q A, b=ord_q B, with a=infinity
when A=0. Let eta=2 for q=0 and eta=0 otherwise. Define

```
nu_q=min(3a+eta, 2b),

rho_q = 2 if 3a+eta>2b and 3 does not divide b+eta;
        0 if 3a+eta>2b and 3 divides b+eta;
        a mod 2 if 3a+eta<2b;
        0 if 3a+eta=2b.
```

All sums count geometric roots individually. Equivalently, an irreducible
factor over K contributes its degree times the displayed local value.
Then the exact geometric genus is

```
g=(n-sum_q nu_q+sum_q rho_q-3)/2 >= 6.                 (1)
```

The lower bound is sharp: A=0,B=1 and A=0,B=p both give genus six.
For coprime A,B with B(0) nonzero, (1) reduces to
g=max(2M+6,2 deg A+3). If gcd(A,B)=1 and B(0)=0, it instead gives
g=max(2M+4,2 deg A+1). No common-factor removal is assumed in (1).

## 2. Generic integrality from the infinity branch

Let c be transcendental over K, E=K(c), and L=K(p,y) with c=I.
The polynomial I-c is irreducible in K[p,y,c] since it is linear and monic
in c. Localization gives an irreducible cubic over E(p), so [L:E(p)]=3.

Here is the complete infinity calculation, which does not require
B(0) nonzero. Dividing the equation by its leading coefficient gives

```
(y²-p³)(y+A/B)+c/(p²B)=0.                              (2)
```

Write m=deg A-deg B, taking m=-infinity when A=0. In the chart
p=r^-2, y=r^-3Y, if m<=1, multiply (2) by r^9. Its initial equation is
Y(Y²-1)=0: the c term has positive order 13+2M, and r³ A/B has positive
order or is zero. The two simple branches near +1 and -1 are exchanged
by (r,Y)->(-r,-Y), so they descend to one point of index two.
The simple branch near zero is fixed by that deck action and descends
to an unramified point.

If m>=2, multiply the same transformed equation additionally by
r^(2m-3). Its initial equation is alpha(Y²-1)=0, where alpha is the
nonzero leading coefficient of A/B. These are again the exchanged pair.
The third branch uses p=q^-1,y=q^-m V: after multiplication by q^(3m),
the initial equation is V²(V+alpha)=0. The simple root -alpha is an
ordinary q-series, hence unramified. In either case the indices 2 and 1
exhaust the cubic degree. Infinity contributes exactly one to ramification.

This also supplies a geometric-integrality argument when the old p=0
valuation test fails. The relative algebraic closure E' of E in L has
degree dividing three, since E'(p) is intermediate between E(p) and L.
If [E':E]=3, then L=E'(p). After algebraic closure of E this would split
into three degree-one p-covers, all unramified, contradicting the index-two
point just constructed. Hence E'=E. In characteristic zero this regular
function field has geometrically integral generic curve.

## 3. The complete discriminant and every finite correction

Literal cubic discriminant expansion gives

```
N=disc_y(I-c)/p⁴
 =4p⁷(B²p³-A²)²+4Ac p²(9B²p³-A²)-27B²c².              (3)
```

Its p-degree is n above. The first term dominates the other two at
infinity. The degrees 2 deg B+3 and 2 deg A cannot be equal, so its
leading coefficient cannot cancel. If A=0 the first term is 4p^13 B^4.

At a nonzero root q of B the exact fixed multiplicity of N is
min(3a,2b). At p=0 it is min(3a+2,2b). To see exactness, compare the
three terms in (3); terms of the smallest order have different powers of
the transcendental c, so cannot cancel. At a root of B with a=0 this
multiplicity is zero. These statements include A=0.

We now count the missing local ramification. At q not zero the dominant
pole equation has the Newton points (0,0),(2,a),(3,b); at zero it has
(0,0),(2,a+2),(3,b+2). Omit the middle point if A=0. The other actual
cubic coefficient lies above the lower hull. If the three-term edge is
dominated by the cubic, all three roots have valuation -(b+eta)/3;
there is one index-three point when 3 does not divide b+eta, and three
unramified points otherwise. If the quadratic point is below that edge,
two roots have valuation -(a+eta)/2 and the remaining root has integral
valuation a-b. The pair contributes a mod2. In the balanced case all
root valuations are integral and all three points are unramified.

These Newton conclusions include the residue equations: on a single
segment the nonzero roots of the initial binomial are simple; in the
balanced case the initial cubic is beta Y³+alpha Y²+c times nonzero
constants. A repeated root of that cubic would force c to be algebraic
over the original coefficient field. Hensel lifting therefore gives the
listed branches. When a=0 at a nonzero root of B, the two finite roots
are also distinct at transcendental c. Thus rho_q is the full local
ramification contribution, not an estimate.

Every remaining zero of N is simple and supplies one index-two point.
Here are the details that prevent an unproved generic-coefficient premise.
A nonconstant polynomial J over a characteristic-zero field has finitely
many geometric critical values. Its critical locus has finitely many
components; on a curve component dJ=0 forces J constant (otherwise extend
the derivative on K(J) to the finite separable component function field).
Point components are finite. Hence the generic value c is not critical.

Away from pB=0 the polynomial is a cubic with invertible leading
coefficient. A triple root would satisfy
y=-A/(3B) and A²+3B²p³=0. The latter is not the zero polynomial, by the
parity of its degrees/valuations; its solutions, and the resulting y, are
algebraic over the coefficient field, incompatible with I=c. At a double
root translate to y=0 and write the monic cubic as y²(y+q), q nonzero.
Its discriminant has first variation -4q³ V(0) under perturbation by V.
Thus a repeated discriminant zero would force J_p=J_y=0 at J=c, already
excluded. Invertible discriminant gives unramified points elsewhere.

The total ramification is therefore
n-sum nu_q+sum rho_q+1. Applying the characteristic-zero
[Riemann--Hurwitz formula](https://stacks.math.columbia.edu/tag/0C1B)
to the degree-three p-map proves (1).

Finally sum nu_q<=2 sum ord_q B=2M. If M>=1, (1) gives
g>=(4M+13-2M-3)/2=M+5>=6. If M=0, nu_0=0 and rho_0=2,
so g=(n-1)/2>=6. This proves the sharp universal lower bound without
requiring any restriction on repeated common factors.

## 4. Actual scalar-time consequence and its boundary

For every nonconstant f in K[z] and lambda in K*, the inherited completed
Hamiltonian flow for S=f(I) cannot have both images of p and y rational
in K(p,y). In particular it cannot give a polynomial source
automorphism at that nonzero time.

Indeed f(I)-f(0) is divisible by p²D, so the actual completed source
carrier theorem applies. Retain the invariant I itself, rather than
f(I), whose generic fibres can split. In the logarithmic coordinates
p=s²+tau, y=sp the actual source bracket is
tau(F_s G_tau-F_tau G_s), and

```
I=tau W(s)+O(tau²),
W(s)=s⁸(A(s²)+s³B(s²)) != 0.
```

The two summands in parentheses have opposite parity. If f begins with
nonconstant term f_k z^k, then

```
exp(lambda delta)(p)-p
 =2lambda k f_k s W(s)^k tau^k+O(tau^(k+1)).              (4)
```

The derivation raises tau order by at least k, so later iterates cannot
cancel this coefficient. The same nonzero displacement holds at n lambda
for every positive integer n. The inherited fixed-input comparison
identifies this computation with the completed actual source flow; it
does not identify arbitrary completed elements with rational Laurent
expressions. A rational specialization would give a nonconstant
[selfmap of the smooth proper generic curve](https://stacks.math.columbia.edu/tag/0BY1).
Since g>=6, Riemann--Hurwitz makes it an automorphism, and its geometric
[automorphism group is finite](https://stacks.math.columbia.edu/tag/0DST).
This contradicts (4) and the scalar-time group law.

The exact primary dependencies for this last paragraph are
[the completed carrier](planar_jc_long_20260906_hamiltonian.md),
[the fixed-input rationality comparison](planar_jc_long_20260906_nonrational.md),
and their audits. The Stacks curve inputs above were read directly.
Constant Hamiltonians, time zero, and cancelling compositions are exempt.
This proves neither that each separate image is transcendental nor that
every universal carrier has high genus. Multipliers of y-degree at least
two and compositions of flows retaining different invariants remain open
to this argument. The current finite-row source response search is also
unaffected: a nonrational all-order flow can have the required finite jet.

## 5. Reproducible exact controls

[The standalone producer](../../04-computation/planar_jc48_sep07_carrier_frontier.py)
imports no research producer. It reconstructs (3) from the literal5x5
Sylvester determinant, the triple-root and discriminant variation
identities, and 364 independent lower Newton hulls: both locations,
0<=a,b<=12, and a=infinity. Twenty-four named pairs A,B include repeated
common factors, irreducible quadratic factors, both balanced thresholds,
vanishing B(0), A=0, and arbitrary-A constant-B inheritance controls.
For each it checks exact fixed discriminant orders, a good squarefree
specialization of the residual divisor, infinity ramification, parity,
the genus bound, and actual first source jets with outer orders1,2,3.
No numerical root calculation or sampled inference proves the theorem.

```
python3 -B 04-computation/planar_jc48_sep07_carrier_frontier.py
python3 -B -O 04-computation/planar_jc48_sep07_carrier_frontier.py
```

Normal and optimized runs pass1018 always-active gates with byte-identical
383-byte output. The frozen pins are:

```
source 04f6008ed12dc33dbd2ef2ce7429fab8a97bc28f77a5cb5c5d232f1a6277db68
output 5b27c5e31ea9bfdd6a4fe3976a00b37effd780f33fa8a4faa0c72e88eac6cef2
semantic f65798e59de1a2cce5a34b8fe387ee5b2087216146e231b094b27d7fbda5fb56
```

The [independent analytic/source audit](planar_jc48_sep07_carrier_frontier_audit.md) passes, including separate normal and optimized replays. The source and output are frozen.
