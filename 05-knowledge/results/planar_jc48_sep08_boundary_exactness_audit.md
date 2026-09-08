# Independent audit: the complete degree-eight radical differential classification

**Audit status: PASS — full analytic/source review, independent partition
reconstruction and uniform identity checks, and normal/optimized/frozen
replay agreement.** September 8, 2026. Root owns primary promotion.

Primary: [boundary exactness](planar_jc48_sep08_boundary_exactness.md).
Source: [standalone exact controls](../../04-computation/planar_jc48_sep08_boundary_exactness.py).
Output: [frozen replay](planar_jc48_sep08_boundary_exactness.out).

The accepted theorem classifies exactness of **the specified differential**
dx/sqrt(N), for every nonzero N in C[x] with degree at most eight, in the
actual field C(x)(sqrt(N)). It retains distinct finite root positions and
the chosen affine coordinate x. It does not assert a Jacobian mate,
classify two-variable polynomials, or use the separate proposed formal
leading-coefficient bridge as a proved dependency. No producer file was
modified in this audit.

## 1. Field types and the complete local/global reduction

If N is square, its chosen square root already belongs to C(x), and the
field is C(x). The two components of the affine equation y^2=N are not
treated as one connected quadratic curve. A different sign for the square
root changes the sign of a primitive and not exactness. Since the scalar
field is C, every nonzero scalar coefficient has a square root.

For nonsquare N I independently checked the local calculations on the
smooth projective normalization of the connected quadratic field. At a
root of odd multiplicity m, x-p=z^2 times an invertible local coordinate
change, and eta has order 1-m. More explicitly, one can choose z with
x-p=z^2 and write eta as a nonzero multiple of
z^(1-m) times an invertible power series in z^2, times dz. At an even root
of multiplicity m, there are two unramified points and eta has order -m/2
at each. Therefore the total possible primitive pole order above a finite
root is max(m-2,0), in both parity cases. Every leading local coefficient
used in this calculation is nonzero.

A nonconstant rational function with pole order e>=1 has derivative of
pole order exactly e+1 in characteristic zero. Thus a multiplicity-two
root gives an unavoidable simple pole, immediately excluding exactness.
At an odd root every exponent of eta in the parameter z is even; no
z^-1 dz term occurs. This does not eliminate the separate global
exactness obstruction in positive genus.

At infinity, for degree n, the orders are

    n odd:  n-3 at the unique ramified point;
    n even: n/2-2 at each of the two unramified points.

For n>=3 these orders are nonnegative. A hypothetical primitive would
therefore have no pole at infinity and no poles away from the finite
roots of N. If

    P=sum_{finite roots with m>=3}(m-2),

its degree as a nonconstant map to P1 is at most P. At each infinity
point, integrating the nonzero leading local term of eta gives local
mapping degree one more than the displayed differential order. Hence
P>=n-2 for odd n, and P>=n/2-1 for even n. The two even-degree local
degrees must **not** be added: the primitive can have different values
at the two points. The primary uses only the valid individual bounds.

This is the same rational-map degree mechanism recovered in
[THM-2071, quadratic-fiber-square-parity-gate, Section 5](../../01-canon/theorems/THM-2071-quadratic-fiber-square-parity-gate.md)
and
[THM-2723, split exact-square-prefix rational-primitive pole capacity](../../01-canon/theorems/THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity.md).
I checked the cited portions: the antecedent compares finite primitive
pole degree with the local degree at infinity. The new proof supplies
its own connected radical-curve orders and does not claim an external
priority result.

For odd n>=3, any partition with at least two roots has P<=n-3, unless
all roots are simple, when P=0 anyway. Thus only a pure odd power survives.
Degree one is treated separately: eta has a double pole at infinity and
does have a primitive. For degree two, the nonsquare curve has two simple
poles at infinity, while the square curve has a finite simple pole.
Both cases are nonexact. Constants are exact in C(x).

In the square case the finite multiplicity-two obstruction leaves only
the pure powers 4,6,8 and the pattern 4+4. At a root p in the latter,

    eta=dx/[sqrt(kappa)(x-p)^2(x-q)^2],
    Res_p(eta)=-2/[sqrt(kappa)(p-q)^3],

which is nonzero for distinct roots. The opposite residue at q does not
cancel this local obstruction to a primitive. This direct calculation
is sufficient; the connected quadratic argument is not misapplied to
the square case.

## 2. Exhaustion of the entire partition universe

I independently generated multiplicity-count vectors (m_1,...,m_8),
where m_j is the number of distinct roots of multiplicity j and
sum j*m_j<=8. This differs from the producer's recursive partition
enumeration. For each vector I generated the actual local places and
orders, computed the total primitive pole capacity, and compared it to
each actual infinity local degree. Square cases and low-degree cases
were separated before applying the connected-curve test.

The degree counts are 1,1,2,3,5,7,11,15,22, giving exactly **67**
partitions. The independent decisions agree with the producer:

| Decision | Number of partitions |
|---|---:|
| Nonzero constant | 1 |
| Pure-power exact type | 7 |
| Two-odd-root exact type | 5 |
| Surviving type requiring a position/coefficient condition | 4 |
| Finite simple-pole rejection | 30 |
| Infinity simple-pole rejection | 1 |
| Primitive-degree rejection | 18 |
| Square-field residue rejection | 1 |

Thus exactly 17 partition types can contain exact differentials, with
four conditional types. Before their position conditions are imposed,
the complete list by degree is:

| Degree | Candidate exact types |
|---|---|
| 0 | constant |
| 1 | 1 |
| 2 | none |
| 3 | 3 |
| 4 | 4; 3+1 |
| 5 | 5 |
| 6 | 6; 5+1; 4+1+1; 3+3 |
| 7 | 7 |
| 8 | 8; 7+1; 6+1+1; 5+3; 5+1+1+1; 4+3+1 |

This exhausts all partitions, rather than testing only the displayed
survivors. The rejection mechanisms depend only on multiplicities and
therefore cover every position of the distinct complex roots. Positions
are restored by the exact formulas below. Coalescing two roots leaves
the stated partition and is handled by another row of the exhaustive
table; no conditional row claims validity after its distinctness
hypotheses fail.

## 3. Pure powers and every two-odd-root orientation

For N=kappa(x-p)^m, direct differentiation gives the primitive

    2/[sqrt(kappa)(2-m)] * (x-p)^(1-m/2)

for m=1 or 3<=m<=8. It lies in the specified radical field; for odd m it
need not lie in C(x). Constants have primitive x/sqrt(kappa), and m=2
is precisely the missing logarithmic case.

For two odd multiplicities 2a+1 and 2b+1, put k=a+b and d=p-q. The
surviving cases have k=1,2,3. The change

    z^2=(x-q)/(x-p),  x-p=d/(z^2-1)

really parametrizes the full quadratic field: the other root difference
is x-q=z^2(x-p), and the original chosen square root is, up to its allowed
sign,

    sqrt(kappa)(x-p)^(k+1) z^(2b+1).

Substitution yields

    eta=-2/[sqrt(kappa)d^k] * (z^2-1)^(k-1) z^(-2b) dz.

Every Laurent monomial exponent is even. Termwise integration has only
the denominators 2j-2b+1, all odd and nonzero. This proves exactness of
3+1,5+1,3+3,7+1,5+3 in every labelled orientation, including cases with
a or b zero. It does not invoke a false restriction to rational functions
of x or omit the possible pole at z=0.

## 4. The three exceptional genus-zero conditions are iff statements

Let u=x-p, d1=p-q, d2=p-r and

    D(u)=(u+d1)(u+d2)=u^2+s*u+t,
    s=d1+d2, t=d1*d2, v^2=D(u).

The hypotheses d1*d2*(d1-d2)!=0 pay all denominators and distinctness.
Each exceptional curve has exactly two branch points and genus zero.
The two points above the even root have opposite residues, while the odd
roots have no residues and infinity has no poles. Nonzero individual
residue proves necessity. Explicit primitives, not a genus-zero
residue heuristic alone, prove sufficiency.

For **4+1+1**, eta=du/[sqrt(kappa)u^2v]. The normalized residue is
-s/(2t). It vanishes iff s=0, equivalently 2p=q+r. Direct differentiation
then gives the primitive -v/[sqrt(kappa)t*u].

For **6+1+1**, eta=du/[sqrt(kappa)u^3v]. The normalized residue is
(3s^2-4t)/(8t^2). It vanishes iff 3s^2=4t. With

    R=-1/(2t*u^2)+3s/(4t^2*u),

the exact identity is

    D R'+D'R/2-1/u^3=(4t-3s^2)/(8t^2*u).

Thus vR/sqrt(kappa) is a primitive precisely on the stated locus. The
condition is over C; it is not restricted to real simple-root positions.

For **4+3+1**, q is the triple root and
eta=du/[sqrt(kappa)u^2(u+d1)v]. The logarithmic derivative at zero of
the residual unit is

    -1/d1-(1/d1+1/d2)/2=-(d1+3d2)/(2d1*d2).

The necessary and sufficient condition is d1+3d2=0, equivalently
3/(p-q)+1/(p-r)=0. The weight three is attached to the triple root in
the reciprocal condition; this labelling has been checked explicitly.
The displayed primitive is

    v*(A*u+B)/[sqrt(kappa)u(u+d1)],
    A=-2/[d1*d2*(d1-d2)], B=-1/(d1*d2).

I independently checked all three uniform primitive identities with
sparse Laurent-polynomial arithmetic over `fractions.Fraction`, without
SymPy or producer imports. For the last identity, set d2=-d1/3 and clear
Q=u(u+d1). Then A=9/(2d1^3), B=3/d1^2, and the identity becomes the
literal polynomial equality

    D*(f'Q-fQ')+(D'/2)*fQ-(u+d1)=0,
    f=A*u+B.

For the sixfold case the same independent computation recovered the
full obstruction (4t-3s^2)/(8t^2*u), including its sign. An independent
binomial-series expansion recovered the opposite normalized residue
coefficient. The producer's symbolic positive controls and named
nonzero-residue controls agree with these computations.

## 5. The exceptional elliptic primitive space is complete

For the remaining partition write

    N=u^5(a*u^3+b*u^2+c*u+d),
    a*d!=0,

with the cubic squarefree. Its three roots are nonzero, distinct and
simple. The exact inversion X=1/u, Y=y/u^4 gives

    Y^2=a+bX+cX^2+dX^3,   eta=-X^2 dX/Y.

The reversed cubic is also squarefree: its finite roots are reciprocals
of the original nonzero roots. Its leading and constant coefficients
are nonzero. Its smooth completion has three finite branch points and
one branch point at infinity, so genus one. This is a connected elliptic
field, not a parametrized genus-zero curve.

At every finite point, dX/Y is regular; at Y=0 one uses Y itself as a
local parameter since the corresponding cubic root is simple. At the
unique point at infinity, ord(X)=-2, ord(Y)=-3 and ord(eta)=-4. Thus a
primitive would have no finite poles and a pole of order exactly three
at infinity. Allowing order at most three gives the same complete space
needed for the argument.

The affine cubic is smooth and normal. A rational function with no finite
poles belongs to its coordinate ring and has a unique expression
A(X)+B(X)Y. The leading pole orders of the two summands are respectively
even and odd, so cannot cancel. A bound of three forces deg A<=1 and
deg B<=0. Therefore

    L(3 infinity)=span_C{1,X,Y}.

This proves completeness of the primitive space, rather than checking
only a convenient ansatz. There is no assumption of irreducibility of
a two-variable generic fibre elsewhere in the repository.

The differential is anti-invariant under Y->-Y. For a primitive
k0+k1X+k2Y, comparing invariant parts of its differential gives k1=0.
The remaining coefficient equation is

    (k2/2)(b+2cX+3dX^2)=-X^2.

Since d!=0, the quadratic coefficient forces k2=-2/(3d); the linear
and constant coefficients then force c=b=0. Conversely these values
give the primitive

    -2Y/(3d)=-2sqrt(N)/(3d*u^4).

Consequently **b=c=0 is necessary and sufficient**. Its locus is
nonempty: a*u^3+d is squarefree whenever a*d!=0. The example
u^5(u^3+1) is exact on an elliptic curve. The two examples with an added
u^2 or u term remain elliptic and have zero residues everywhere, but
fail the coefficient equation. This explicitly pays the missing global
primitive condition beyond residues.

## 6. Scope and hostile controls

All scalar factors are preserved. In the elliptic statement a,b,c,d
include the original scalar of N; the condition b=c=0 is unchanged by
multiplication by a nonzero scalar. In the genus-zero formulas the
chosen sqrt(kappa) is retained. Translations u=x-p preserve dx; the
primary does not claim that an arbitrary projective reparametrization
preserves eta or defines an automorphism of the DG surface.

I checked the literal rational pair F=x^3t^2, G=1/(x^2t): its Jacobian
is one, while the leading radical differential has primitive
-2/sqrt(x). It is a valid hostile to insisting on a primitive in C(x).
Conversely F=t^2+x^3 has exact constant-leading differential but no
rational mate. Its generic fibre is a smooth elliptic curve and the
nonzero relative differential -dx/(2t) is holomorphic on its compact
normalization, including infinity. A rational primitive of a nonzero
holomorphic differential cannot exist. Thus passing the leading gate
is not sufficient for a full rational or polynomial Jacobian mate.

The classification is confined to the specified one-variable differential.
The separate formal leading-coefficient theorem must be independently
proved before it is used to infer restrictions on Jacobian pairs. No
claim of such a bridge is imported into this audit, and no conclusion
about a surviving partition's actual Jacobian realization is made.

## 7. Source review, reproduction and frozen pins

I read the complete standalone source. Its mathematical implementation
uses SymPy for exact symbolic arithmetic; it imports no other producer.
All checks use explicit exceptions and remain active under `python3 -O`.
The proof is analytic. The finite partition bank validates the complete
declared degree universe, and the symbolic identities cover all allowed
parameters rather than a sampling of root positions or possible mates.

The reviewed controls include:

* all 67 multiplicity partitions and the complete survivor list in every
  degree, with branch-count genus checks for every nonsquare survivor;
* local and infinity orders, the degree-one/two exceptions, square-field
  residues, and all pure-power primitives;
* all nine labelled orientations of the two-odd-root survivors, including
  the exact parameter substitution and the absence of logarithmic terms;
* all three genus-zero residues and uniform primitives, with positive and
  negative controls satisfying the required distinct-root hypotheses;
* elliptic inversion, the coefficient obstruction and complete primitive
  space argument, literal elliptic positive/negative controls, and direct
  differentiation of the positive primitive in the original coordinate;
* actual rational-pair controls and the explicit separation between a
  passed leading condition and full Jacobian integrability.

The bounded list of monomial orders in the source is only a control of
the elliptic pole-space calculation; the primary's parity argument
proves the space complete in all polynomial degrees. Likewise the
partition enumeration is exhaustive because the external degree bound
is eight, not because of a guessed bound on primitive complexity.

Independent commands from the worktree root:

```bash
python3 04-computation/planar_jc48_sep08_boundary_exactness.py
python3 -O 04-computation/planar_jc48_sep08_boundary_exactness.py
```

Both completed successfully. The independent normal and optimized
outputs are byte-identical to the frozen **195-gate**, **632-byte**
output. The semantic control manifest hash is
`91100476475b4851378fa4b6a21cf69ecd002bdabc77cc6938a750397213564e`.

| Accepted artifact | Bytes | SHA256 |
|---|---:|---|
| Source | 10,124 | `bdc72806752cdfe31e2250b48cc6e26b19ead0bcadfc7fdaea887377ee579af0` |
| Frozen output and each independent replay | 632 | `c8ea8699d0d7f13eafa3961c017fea1bfe39d1c26a26cf231b609dc9d5105d66` |
| Primary before status promotion | 16,443 | `cecc6a9877e41285211b4bdd2b296d1e074a8ed8b0a841a1e925c7f16ac7f5e1` |

No mathematical correction remains. This audit is frozen for root's
promotion/checkpoint with the differential-only scope preserved.
